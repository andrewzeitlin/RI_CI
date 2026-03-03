#  example.R — Illustrative use of ri()
#
#  All examples use synthetic data and run self-contained.
#  Requires: fixest, ggplot2, tictoc (install.packages(c("fixest", "ggplot2", "tictoc")))

#  Before running this script, source ri.R manually:
#    source("path/to/ri.R")

library(fixest)
library(ggplot2)
library(tictoc)

set.seed(20260303)


#------------------------------------------------------------------------------#
#  1. Synthetic data                                                         ----
#------------------------------------------------------------------------------#

#  Cross-sectional design: N individuals in K clusters.
#  Treatment is clustered — entire clusters are assigned en bloc,
#  then half the units within each cluster are selected.

N   <- 300   # individuals
K   <- 10    # clusters (e.g. villages, schools, blocks)
R   <- 400   # RI permutations for this example

cluster_id <- rep(paste0("cl_", formatC(seq_len(K), width = 2, flag = "0")), each = N / K)
unit_id    <- paste0("unit_", formatC(seq_len(N), width = 3, flag = "0"))

#  Within-cluster randomisation: ~half treated per cluster
treated <- as.integer(ave(runif(N), cluster_id, FUN = function(x) rank(x) > length(x) / 2))

alpha_true <- 0.3   # true direct effect

df <- data.frame(
    unit_id    = unit_id,
    cluster_id = cluster_id,
    treated    = treated,
    y          = rnorm(N) + alpha_true * treated   # DGP: Y = N(0,1) + alpha * W
)


#------------------------------------------------------------------------------#
#  2. Build permutation matrix T0                                            ----
#
#  Each column is an independent re-randomisation respecting the original
#  design: treatment is permuted within each cluster block.
#------------------------------------------------------------------------------#

perm_within_block <- function(treat_vec, block_vec) {
    out <- treat_vec
    for (b in unique(block_vec)) {
        idx      <- which(block_vec == b)
        out[idx] <- sample(treat_vec[idx])
    }
    out
}

#  T0[["treated"]]: data.frame with unit_id column + R permutation columns.
T0_treated           <- data.frame(unit_id = unit_id)
T0_treated[paste0("w_", seq_len(R))] <- replicate(R, perm_within_block(treated, cluster_id))

T0 <- list(treated = T0_treated)


#------------------------------------------------------------------------------#
#  Test statistic choice and model performance                               ----
#
#  stat = 't' (default, Chung–Romano studentised RI): the SE is re-estimated on
#  every permutation draw, making the test statistic approximately pivotal under
#  the null and yielding better size control — especially in cluster-randomised
#  designs where the OLS SE can substantially understate the true variance.
#  Cluster the SE at the randomisation unit (cluster = ~cluster_id) so the
#  studentisation reflects the design.
#
#  stat = 'b' (raw coefficient): valid for RI and identical in p-value to a
#  fixed-denominator t (dividing every draw by the observed SE merely scales the
#  null distribution without changing ranks).  Use only when speed is critical
#  and the SE is expected to be stable across permutations.
#
#  Memory optimisations (fixest only, active when stat = 'b' or 't'):
#  ri() automatically injects lean = TRUE (suppresses score matrix, fitted
#  values, and residuals at construction time) into any unqualified feols()
#  call in the model.  The same fields are also
#  nulled post-fit as a safety net.  Set lean = FALSE in the ri() call to
#  disable injection (e.g. when a custom stat needs model-matrix data).
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
#  3. Function interface (recommended)                                       ----
#------------------------------------------------------------------------------#

#  Pass a function that accepts a data frame and returns a fitted model.
#  The function receives the permuted df on each RI draw.

tic("Example 1")
res1 <- ri(
    df    = df,
    model = function(df) feols(y ~ treated | cluster_id, cluster = ~cluster_id, data = df,
                               warn = FALSE, notes = FALSE),
    tx    = "treated",
    T0    = T0,
    id    = "unit_id",
    n.cores = 1L     # serial for a quick illustration
)
toc()

cat("\n── Example 1: function interface, stat = 't' (default) ──\n")
cat(sprintf("  Observed t-stat : %.3f\n",   res1$teststat))
cat(sprintf("  RI p-value      : %.3f\n\n", res1$p))

#  Null distribution plot
null_dist <- as.numeric(res1$testdist)
t_obs     <- as.numeric(res1$teststat)
p1 <- ggplot(data.frame(t = null_dist), aes(x = t)) +
    geom_histogram(binwidth = 0.2, fill = "steelblue", colour = "white", alpha = 0.8) +
    geom_vline(xintercept =  t_obs, colour = "firebrick", linewidth = 1) +
    geom_vline(xintercept = -t_obs, colour = "firebrick", linewidth = 1, linetype = "dashed") +
    annotate("text",
             x = t_obs, y = Inf,
             label   = sprintf("t[obs] == %.2f", t_obs),
             parse   = TRUE,
             hjust   = -0.15, vjust = 1.8,
             colour  = "firebrick", size = 3.5) +
    annotate("text",
             x = mean(range(null_dist)), y = Inf,
             label   = sprintf("two-sided p = %.3f  (R = %d)", res1$p, length(null_dist)),
             hjust   = 0.5, vjust = 3.5,
             size    = 3.5) +
    labs(
        title = "RI null distribution of t-statistic (Example 1)",
        x     = "t-statistic under the null",
        y     = "Count"
    ) +
    theme_bw()

print(p1)


#------------------------------------------------------------------------------#
#  4. Quoted call interface                                                  ----
#
#  Equivalent to Example 1. The quoted call must reference `df` as the data
#  argument; ri() binds the (possibly permuted) frame to that name before eval.
#------------------------------------------------------------------------------#

tic("Example 2")
res2 <- ri(
    df    = df,
    model = quote(feols(y ~ treated | cluster_id, cluster = ~cluster_id, data = df,
                        warn = FALSE, notes = FALSE)),
    tx    = "treated",
    T0    = T0,
    id    = "unit_id",
    n.cores = 1L
)
toc()

cat("── Example 2: quoted call interface ──\n")
cat(sprintf("  RI p-value: %.3f  (same as Example 1: %s)\n\n",
            res2$p, isTRUE(all.equal(res1$p, res2$p))))


#------------------------------------------------------------------------------#
#  5. stat = 'b' (raw coefficient instead of t-statistic)                   ----
#------------------------------------------------------------------------------#

tic("Example 3")
res3 <- ri(
    df    = df,
    model = function(df) feols(y ~ treated | cluster_id, cluster = ~cluster_id, data = df,
                               warn = FALSE, notes = FALSE),
    tx    = "treated",
    T0    = T0_treated,   # bare data.frame shorthand: equivalent to list(treated = T0_treated)
    id    = "unit_id",
    stat  = 'b',
    n.cores = parallel::detectCores() - 1L   # use all but one core for faster processing
)
toc()

cat("── Example 3: stat = 'b' (coefficient) ──\n")
cat(sprintf("  Observed coef : %.3f\n",   res3$teststat))
cat(sprintf("  RI p-value    : %.3f\n\n", res3$p))


#------------------------------------------------------------------------------#
#  6. Custom stat function                                                   ----
#
#  Pass any function(fit) -> numeric to extract an arbitrary test statistic.
#  Here: use the t-value from a richer ANCOVA specification.
#------------------------------------------------------------------------------#

#  Add a pre-treatment baseline covariate to the data
df$y_base <- rnorm(N)

tic("Example 4")
res4 <- ri(
    df    = df,
    model = function(df) feols(y ~ treated + y_base | cluster_id, cluster = ~cluster_id,
                               data = df, warn = FALSE, notes = FALSE),
    tx    = "treated",
    T0    = T0,
    id    = "unit_id",
    stat  = function(fit) coeftable(fit)["treated", "t value"],
    n.cores = parallel::detectCores() - 1L   # use all but one core for faster processing
)
toc()

cat("── Example 4: custom stat function (t-value from ANCOVA) ──\n")
cat(sprintf("  Observed t-stat : %.3f\n",   res4$teststat))
cat(sprintf("  RI p-value      : %.3f\n\n", res4$p))


#------------------------------------------------------------------------------#
#  7. Subset of available permutations                                       ----
#
#  T0 has R=400 columns. Use only the first 100 for a quick check run.
#------------------------------------------------------------------------------#

tic("Example 5")
res5 <- ri(
    df    = df,
    model = function(df) feols(y ~ treated | cluster_id, cluster = ~cluster_id, data = df,
                               warn = FALSE, notes = FALSE),
    tx    = "treated",
    T0    = T0,
    id    = "unit_id",
    R     = 100L,
    n.cores = parallel::detectCores() - 1L   # use all but one core for faster processing
)
toc()

cat("── Example 5: use only R=100 of 400 available permutations ──\n")
cat(sprintf("  RI p-value (R=100): %.3f\n\n", res5$p))


#------------------------------------------------------------------------------#
#  8. Null rejection rate check                                              ----
#
#  Under the null (alpha=0), the two-sided RI p-value should be approximately
#  Uniform(0,1), so the rejection rate at 5% should be ~0.05.
#  Run M outer draws, each with its own Y0; check rejection rate.
#------------------------------------------------------------------------------#

cat("── Example 6: null rejection rate (alpha=0, M=200 outer draws) ──\n")

df_null    <- df
M          <- 200
rejections <- logical(M)
pvals_null <- numeric(M)

tic("Example 6")
for (m in seq_len(M)) {
    df_null$y  <- rnorm(N)   # pure noise; no treatment effect
    res_null   <- ri(
        df    = df_null,
        model = function(df) feols(y ~ treated | cluster_id, cluster = ~cluster_id,
                                   data = df, warn = FALSE, notes = FALSE),
        tx    = "treated",
        T0    = T0,
        id    = "unit_id",
        n.cores = parallel::detectCores() - 1L   # use all but one core for faster processing
    )
    pvals_null[m] <- res_null$p
    rejections[m] <- !is.na(res_null$p) && res_null$p < 0.05
}
toc()

cat(sprintf("  Rejection rate at alpha=0.05: %.3f  (target ~0.05; M=%d, R=%d)\n\n",
            mean(rejections, na.rm = TRUE), M, R))

#  P-value distribution plot: should be approximately Uniform(0,1) under the null
p6 <- ggplot(data.frame(p = pvals_null), aes(x = p)) +
    geom_histogram(breaks = seq(0, 1, by = 0.05),
                   fill = "steelblue", colour = "white", alpha = 0.8) +
    geom_hline(yintercept = M * 0.05, colour = "firebrick",
               linewidth = 0.8, linetype = "dashed") +
    annotate("text",
             x = 0.98, y = M * 0.05,
             label = "expected\n(uniform)", hjust = 1, vjust = -0.3,
             colour = "firebrick", size = 3.2) +
    labs(
        title    = "Distribution of RI p-values under the null (Example 6)",
        subtitle = sprintf("M = %d outer draws, R = %d permutations each; rejection rate = %.3f",
                           M, R, mean(rejections, na.rm = TRUE)),
        x        = "RI p-value",
        y        = "Count"
    ) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
    theme_bw()

print(p6)


#------------------------------------------------------------------------------#
#  7. Joint F-test of two independently randomised treatments                ----
#
#  Two binary treatments (t1, t2) are independently randomised within clusters.
#  T0 must be a named list when tx has more than one element.  On each RI draw r,
#  ri() assigns column r from T0[["t1"]] and column r from T0[["t2"]]
#  simultaneously — the permutations for the two treatments are generated
#  independently but applied jointly.
#
#  Test statistic: the cluster-robust F-stat for the joint null H0: β_t1 = β_t2 = 0.
#  A custom stat = function(fit) is required; ri() has no built-in joint test.
#  Because stat is a custom function, lean injection is automatically disabled
#  (the full fit object is needed for fitstat()), which is the correct behaviour.
#------------------------------------------------------------------------------#

set.seed(20260303 + 7L)
t1 <- as.integer(ave(runif(N), cluster_id, FUN = function(x) rank(x) > length(x) / 2))
t2 <- as.integer(ave(runif(N), cluster_id, FUN = function(x) rank(x) > length(x) / 2))

df7 <- data.frame(
    unit_id    = unit_id,
    cluster_id = cluster_id,
    t1         = t1,
    t2         = t2,
    y          = rnorm(N) + 0.3 * t1   # true effect on t1 only; t2 has no effect
)

#  Independent permutation tables — one data.frame per treatment
T0_t1 <- data.frame(unit_id = unit_id)
T0_t1[paste0("w_", seq_len(R))] <- replicate(R, perm_within_block(t1, cluster_id))

T0_t2 <- data.frame(unit_id = unit_id)
T0_t2[paste0("w_", seq_len(R))] <- replicate(R, perm_within_block(t2, cluster_id))

tic("Example 7")
res7 <- ri(
    df    = df7,
    model = function(df) feols(y ~ t1 + t2 | cluster_id, cluster = ~cluster_id,
                               data = df, warn = FALSE, notes = FALSE),
    tx    = c("t1", "t2"),
    T0    = list(t1 = T0_t1, t2 = T0_t2),   # named list required when length(tx) > 1
    id    = "unit_id",
    stat  = function(fit) fitstat(fit, "f")$f$stat,
    n.cores = 1L
)
toc()

cat("── Example 7: joint F-test of two independent treatments ──\n")
cat(sprintf("  Observed F-stat : %.3f\n",   as.numeric(res7$teststat)))
cat(sprintf("  RI p-value      : %.3f\n\n", as.numeric(res7$p)))


#------------------------------------------------------------------------------#
#  8. Interaction with a fixed covariate — only t is permuted               ----
#
#  Model: y ~ t * x | cluster_id   (R expands to y ~ t + x + t:x | cluster_id)
#
#  Only the treatment t appears in tx and T0.  On each RI draw, ri() replaces
#  df[["t"]] with the permuted column; x is untouched.  The model formula then
#  builds t:x on-the-fly from the permuted t and the original x — no separate
#  permutation table for the interaction is required.
#
#  Covariate x is not in tx, so it is never permuted.  Any derived term (t:x,
#  t^2, log(t), …) that depends only on tx variables is automatically permuted
#  via the formula; any term involving only non-tx variables stays fixed.
#
#  Test statistic: cluster-robust Wald chi-squared for joint exclusion of t
#  and t:x — i.e., H0: main effect = 0 AND heterogeneous treatment effect = 0.
#------------------------------------------------------------------------------#

set.seed(20260303 + 8L)
x8 <- rnorm(N)   # unit-level continuous covariate (not permuted)

df8 <- data.frame(
    unit_id    = unit_id,
    cluster_id = cluster_id,
    t          = treated,
    x          = x8,
    y          = rnorm(N) + 0.3 * treated + 0.5 * x8 + 0.4 * treated * x8
    # DGP: main effect alpha=0.3, covariate beta=0.5, HTE gamma=0.4
)

tic("Example 8")
res8 <- ri(
    df    = df8,
    model = function(df) feols(y ~ t * x | cluster_id, cluster = ~cluster_id,
                               data = df, warn = FALSE, notes = FALSE),
    tx    = "t",           # only the treatment — formula builds t:x automatically
    T0    = T0_treated,    # same permutation table used for the original 'treated'
    id    = "unit_id",
    stat  = function(fit) {
        #  Cluster-robust Wald chi-squared for joint exclusion of t and t:x.
        #  vcov(fit) returns the cluster-robust vcov because cluster= was
        #  specified in the model.
        coefs <- c("t", "t:x")
        b     <- coef(fit)[coefs]
        V     <- vcov(fit)[coefs, coefs]
        as.numeric(t(b) %*% solve(V) %*% b)
    },
    n.cores = 1L
)
toc()

cat("── Example 8: interaction t*x; only t permuted ──\n")
cat(sprintf("  Observed Wald stat : %.3f  (alpha=0.3, gamma=0.4)\n",
            as.numeric(res8$teststat)))
cat(sprintf("  RI p-value         : %.3f\n\n", as.numeric(res8$p)))


#------------------------------------------------------------------------------#
#  Notes on non-individual randomisation units                               ----
#
#  `id` identifies the *unit of randomisation*, which need not be the unit of
#  observation.  T0 should have one row per randomisation unit.  ri() expands
#  T0 up front — once, before the permutation loop — by evaluating
#  match(df[[id]], T0[[id]]) and using the resulting index to materialise a
#  full N_obs × R integer matrix.  No merge occurs inside the loop; each
#  permutation step is a single column extraction from this pre-expanded matrix.
#
#  Cluster-randomised trials:  pass id = "cluster_id".  T0 has J rows (one per
#  cluster) and ri() maps each observation to its cluster's column, so every
#  observation receives its cluster's permuted treatment.  The up-front
#  expansion handles the J → N fanout automatically.
#
#  Panel data (N individuals × T rounds):  two patterns work.
#
#  Pattern A — pass df at observation level (N×T rows) and T0 at individual
#  level (N rows).  ri() expands via the match, replicating each individual's
#  permuted treatment across all of their rounds.  This is fully safe regardless
#  of row ordering; no recycling occurs.
#
#  Mixed assignment granularity: when treatments are assigned at different levels
#  (e.g. t1 individually, t2 at cluster level), supply id as a character vector
#  of the same length as tx.  ri() uses each element as the ID column for the
#  corresponding T0 table.  Named or unnamed (positional) vectors both work:
#
#    ri(
#        df    = df,
#        model = function(df) feols(y ~ t1 + t2 | fe, cluster = ~cl, data = df),
#        tx    = c("t1", "t2"),
#        T0    = list(t1 = T0_individual,   # one row per unit
#                     t2 = T0_cluster),     # one row per cluster
#        id    = c(t1 = "unit_id", t2 = "cluster_id")
#    )
#
#  df must contain both unit_id and cluster_id.  The two T0 tables are expanded
#  independently to N_obs rows via their respective ID columns.
#
#  Pattern B — pass df at individual level; let the model function do the panel
#  merge internally.  More explicit, and useful when the panel frame is too
#  large to modify cheaply on every permutation:
#
#    df_ind <- df[!duplicated(df$unit_id), ]   # one row per individual
#
#    ri(
#        df    = df_ind,
#        model = function(df_ind) {
#            panel <- merge(panel_data, df_ind[, c("unit_id", "treated")],
#                           by = "unit_id")
#            feols(y ~ treated | unit_id + t, cluster = ~cluster_id,
#                  data = panel, warn = FALSE, notes = FALSE)
#        },
#        tx  = "treated",
#        T0  = T0,   # T0 has one row per individual
#        id  = "unit_id"
#    )
#------------------------------------------------------------------------------#
