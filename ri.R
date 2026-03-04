ri <- function(
    df,           # data.frame or data.table with outcome and initial treatment cols
    model,        # function(df)->fit, OR a quoted call referencing `df` (e.g. quote(feols(..., data=df)))
    tx,           # character: treatment variable name(s); must match names(T0)
    T0,           # named list: T0[[txvar]] has an ID column + R permutation columns
    id,           # character: ID column name(s) — REQUIRED, no default.
                  # Single string: same column used for every element of T0.
                  # Character vector of length(tx): id[i] is the ID column for
                  # T0[[tx[i]]], allowing e.g. unit-level id for an individually-
                  # assigned treatment and cluster-level id for a cluster-assigned
                  # treatment.  May be unnamed (positional) or named (names = tx).
    stat    = 't',  # 't' (default; Chung–Romano studentised RI), 'b' (coefficient; faster),
                    #     or function(fit)->numeric; see extract_stat below
    R       = NULL, # permutations to use; default = all columns in T0 after dropping id
    n.cores = NULL, # NULL=auto-detect (SLURM_CPUS_PER_TASK, then detectCores()-1); 1L=serial
    lean    = TRUE  # inject lean=TRUE into unqualified feols() calls inside the
                    # model function, suppressing large internal objects at construction
                    # time.  Set FALSE to disable (e.g. custom stat= needing model-matrix
                    # data).  Has no effect on fixest::feols() qualified calls or
                    # non-fixest model functions.
) {

    #--------------------------------------------------------------------------#
    #  1. Core detection                                                     ----
    #--------------------------------------------------------------------------#

    if (is.null(n.cores)) {
        slurm_env <- Sys.getenv("SLURM_CPUS_PER_TASK")
        n.cores <- if (nzchar(slurm_env)) {
            as.integer(slurm_env)
        } else {
            max(1L, parallel::detectCores(logical = FALSE) - 1L)
        }
    }
    n.cores <- as.integer(n.cores)

    #--------------------------------------------------------------------------#
    #  2. Align T0 rows to df's ID order; convert to plain integer matrices  ----
    #--------------------------------------------------------------------------#

    #  Convenience: accept a bare data.frame / data.table for T0 when there is
    #  only one treatment variable.  Wrap it into the expected named list so the
    #  rest of the function is agnostic to how the caller supplied it.
    if (is.data.frame(T0)) {
        if (!is.character(tx) || length(tx) != 1L) stop(
            "ri: T0 may be a bare data.frame only when tx is a single string"
        )
        T0 <- setNames(list(T0), tx)
    }

    #  Normalise id to a named character vector keyed by tx.
    #  Accepts two forms:
    #    (a) single string → the same ID column is used for every T0 element.
    #    (b) character vector of length(tx) → id[i] is the ID column for the
    #        i-th treatment, allowing different granularities (e.g. unit-level
    #        for an individually-assigned treatment, cluster-level for a
    #        cluster-assigned treatment).  May be unnamed (positional) or named
    #        (names must match tx).  df must contain all referenced ID columns.
    if (length(id) == 1L) {
        id <- setNames(rep(id, length(tx)), tx)
    } else {
        if (!is.character(id) || length(id) != length(tx)) stop(
            "ri: id must be a single string or a character vector of the same length as tx"
        )
        if (is.null(names(id))) {
            names(id) <- tx
        } else {
            if (!setequal(names(id), tx)) stop(
                "ri: names(id) must match tx when id is a named vector"
            )
            id <- id[tx]   # reorder to align with tx
        }
    }

    stopifnot(
        is.character(id), length(id) == length(tx), setequal(names(id), tx),
        is.character(tx), length(tx) >= 1L,
        is.list(T0),
        setequal(names(T0), tx)
    )

    T0_mats <- lapply(names(T0), function(txvar) {
        tbl    <- T0[[txvar]]
        id_col <- id[[txvar]]
        idx    <- match(df[[id_col]], tbl[[id_col]])
        if (anyNA(idx)) stop(sprintf(
            "ri: %d ID(s) in df[['%s']] not found in T0[['%s']]",
            sum(is.na(idx)), id_col, txvar
        ))
        cols <- setdiff(names(tbl), id_col)
        mat  <- as.matrix(as.data.frame(tbl)[idx, cols, drop = FALSE])
        storage.mode(mat) <- "integer"
        mat   # N x R_avail
    })
    names(T0_mats) <- names(T0)

    #--------------------------------------------------------------------------#
    #  3. Number of permutations                                             ----
    #--------------------------------------------------------------------------#

    R_avail <- ncol(T0_mats[[1L]])
    if (is.null(R)) R <- R_avail
    if (R > R_avail) stop(sprintf(
        "ri: R=%d requested but T0 only has %d permutation columns", R, R_avail
    ))

    #--------------------------------------------------------------------------#
    #  4. Model and stat interfaces                                          ----
    #--------------------------------------------------------------------------#

    #  call_model: accepts a data frame, returns a fitted model object.
    #  Supports two calling conventions:
    #    (a) model is a function:      model(dat)
    #    (b) model is a quoted call:   eval(model) with df=dat in local scope.
    #        The quoted call must reference `df` as its data argument, e.g.
    #        quote(feols(y ~ x | fe, cluster = ~cl, data = df)).
    #        Note for socket cluster users (Windows): packages required by the
    #        model must be loaded inside the model function or via clusterEvalQ.
    call_model <- function(dat) {
        if (is.function(model)) {
            model(dat)
        } else {
            df <- dat
            eval(model)
        }
    }

    #  extract_stat: pulls the scalar (or vector, if length(tx)>1) test statistic.
    #  When only.coef = TRUE was injected (stat = 'b' path), feols() returns a
    #  named numeric vector rather than a full fit object; the is.numeric() branch
    #  handles that case directly.
    extract_stat <- function(fit) {
        if (is.function(stat)) return(stat(fit))
        #  only.coef = TRUE: feols returns a named numeric vector directly
        if (is.numeric(fit)) return(fit[tx])
        b <- tryCatch(
            coef(fit)[tx],
            error = function(e) rep(NA_real_, length(tx))
        )
        if (identical(stat, 'b')) return(b)
        #  stat == 't': Chung–Romano studentised RI.  SE is re-estimated on every
        #  permutation so the statistic is approximately pivotal under the null.
        #  Prefers fit$se (fixest; reflects any cluster= supplied in the model
        #  function) then falls back to sqrt(diag(vcov(fit))).  Use cluster= in
        #  the model when treatment is assigned at the cluster level.
        se <- tryCatch(
            fit$se[tx],                          # fixest: $se slot
            error = function(e) tryCatch(
                sqrt(diag(vcov(fit)))[tx],       # base lm / sandwich vcov
                error = function(e2) rep(NA_real_, length(tx))
            )
        )
        b / se
    }

    #--------------------------------------------------------------------------#
    #  4b. Lean injection and post-fit stripping (fixest, built-in stats)   ----
    #
    #  Two complementary memory optimisations, both active when stat is 'b'
    #  or 't' (disabled for custom stat functions, which may need the full fit):
    #
    #  (1) Pre-estimation injection (lean = TRUE, default): ri() rewrites
    #      call_model to shadow any unqualified feols() call with a wrapper
    #      that injects lean = TRUE (suppresses score matrix, fitted values,
    #      residuals at construction time — the deeper saving).
    #      For stat = 'b', also injects only.coef = TRUE, which skips SE
    #      computation entirely and returns a named numeric vector.
    #      Does NOT intercept fixest::feols() qualified calls.
    #      Environment setup is done ONCE at ri() call time, not per permutation.
    #
    #  (2) Post-fit stripping: after each fit is returned, ri() also nulls
    #      the same large fields as a safety net, covering cases where
    #      injection was disabled or the call used fixest::feols() directly.
    #      (When only.coef = TRUE the fit is a numeric vector; strip_fit is
    #      a no-op in that case since it checks inherits(fit, "fixest").)
    #--------------------------------------------------------------------------#

    do_strip     <- !is.function(stat)
    do_only_coef <- identical(stat, 'b')

    strip_fit <- function(fit) {
        if (!inherits(fit, "fixest")) return(fit)
        fit[["scores"]]        <- NULL   # N x k score matrix (used for vcov)
        fit[["fitted.values"]] <- NULL   # N-vector
        fit[["residuals"]]     <- NULL   # N-vector
        fit
    }

    if (lean && do_strip && isNamespaceLoaded("fixest")) {
        .feols_lean <- local({
            base <- fixest::feols
            oc   <- do_only_coef
            function(...) {
                args <- list(...)
                if (is.null(args[["lean"]]))      args[["lean"]]      <- TRUE
                if (oc && is.null(args[["only.coef"]])) args[["only.coef"]] <- TRUE
                do.call(base, args)
            }
        })
        if (is.function(model)) {
            #  Create the injected function ONCE at setup time.
            #  Shadowing feols in the model's lexical scope via a fresh
            #  parent environment means new.env() runs only here, not per
            #  permutation.
            .env_inject       <- new.env(parent = environment(model))
            .env_inject$feols <- .feols_lean
            .fn_injected      <- model
            environment(.fn_injected) <- .env_inject
            call_model <- function(dat) .fn_injected(dat)
        } else {
            call_model <- function(dat) {
                df    <- dat
                feols <- .feols_lean
                eval(model)
            }
        }
    }

    #--------------------------------------------------------------------------#
    #  5. Observed test statistic                                            ----
    #--------------------------------------------------------------------------#

    fit_obs <- tryCatch(call_model(df), error = function(e) NULL)
    if (is.null(fit_obs)) stop("ri: model failed on observed data")
    if (do_strip) fit_obs <- strip_fit(fit_obs)
    teststat <- extract_stat(fit_obs)

    #--------------------------------------------------------------------------#
    #  6. Per-permutation worker function                                    ----
    #--------------------------------------------------------------------------#

    run_perm <- function(r) {
        dat_r <- df
        for (txvar in names(T0_mats)) dat_r[[txvar]] <- T0_mats[[txvar]][, r]
        fit_r <- tryCatch(call_model(dat_r), error = function(e) NULL)
        if (is.null(fit_r)) return(rep(NA_real_, length(tx)))
        if (do_strip) fit_r <- strip_fit(fit_r)
        extract_stat(fit_r)
    }

    #--------------------------------------------------------------------------#
    #  7. Parallel dispatch                                                   ----
    #
    #  n.cores = 1L  → serial lapply (no overhead)
    #  n.cores > 1L, Unix/Linux/macOS → mclapply (fork; shared memory via CoW;
    #      works correctly within a SLURM job allocated with --cpus-per-task N)
    #  n.cores > 1L, Windows → PSOCK socket cluster (explicit object export)
    #--------------------------------------------------------------------------#

    on_unix <- .Platform$OS.type == "unix"

    if (n.cores > 1L && on_unix) {
        perm_list <- parallel::mclapply(
            seq_len(R), run_perm,
            mc.cores = n.cores
        )
    } else if (n.cores > 1L) {
        cl <- parallel::makeCluster(n.cores, type = "PSOCK")
        parallel::clusterExport(
            cl,
            varlist = c("df", "T0_mats", "tx", "model", "stat",
                        "call_model", "extract_stat", "strip_fit", "do_strip",
                        "run_perm"),
            envir   = environment()
        )
        on.exit(parallel::stopCluster(cl), add = TRUE)
        perm_list <- parallel::parLapply(cl, seq_len(R), run_perm)
    } else {
        perm_list <- lapply(seq_len(R), run_perm)
    }

    testdist <- do.call(rbind, perm_list)   # R x length(tx)

    #--------------------------------------------------------------------------#
    #  8. Two-sided p-value                                                  ----
    #--------------------------------------------------------------------------#

    p_left  <- colMeans(testdist < teststat, na.rm = TRUE)
    p_right <- colMeans(testdist > teststat, na.rm = TRUE)
    p       <- pmin(2 * pmin(p_left, p_right), 1)

    list(
        teststat = teststat,
        testdist = testdist,
        p.left   = p_left,
        p.right  = p_right,
        p        = p
    )
}
