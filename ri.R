ri <- function(df,model,tx,T0,stat='b',R=NULL, parallelize=TRUE, ...){
  results <- list() # Container for results. 
  
  #  Packages
  require(fixest)
  require(foreach)
  require(doParallel)
  
  #  Unpack optional arguments, if supplied
  args <- list(...)
  if (length(args) > 0){
    for(i in 1:length(args)) {
      assign(x = names(args)[i], args[[i]])
    }
  }
  
  #  Set up parallel back-end
  if (parallelize){
    if (getDoParWorkers() > 1) {
      already.set.up <- TRUE
    }else{
      already.set.up <- FALSE
      all.cores <- detectCores()
      cl <- makeCluster(all.cores - 1)
      registerDoParallel(cl) # Leaving 1 core for headroom.
    }
  }else{
    registerDoSEQ()
  }
  
  #  If user does not supply a number of repetitions, take this from the width of the first permutation matrix
  if (is.null(R)){
    R <- dim(T0[[1]])[2] -1 # Number of replications to consider, minus one column for the identifer.
  }
  
  #  Expand feasible randomizations and remove identifier (assumes this is variable `schoolid`).
  to_permute <- names(T0) 
  for (txvar in to_permute){
    T0[[txvar]] <- df %>%
      select(schoolid) %>%
      left_join(
        .,
        T0[[txvar]],
        by='schoolid'
      ) %>%
      select(-schoolid)
  }
  
  #  If using fixest, set dataset to the argument "df" 
  if (any(grepl('feols', model, fixed=TRUE))){
    setFixest_estimation(data=df)
  }
  
  #  Point estimate 
  lm1 <- eval(model) 
  if (stat=='b') {
    results$teststat <- lm1$coefficients[tx]
  }else if (stat=='t'){
    results$teststat <- lm1$coefficients[tx] / lm1$se[tx]
  }
  
  
  ################################################################
  #  Distribution of test statistic under the null -------------------------
  ################################################################
  
  results$testdistribution <- matrix(data=NA,nrow=R,ncol=length(tx))
  
  #  Begin loop (parallel)
  results$testdistribution <- foreach (r=1:R,.combine=rbind, .export=names(args)) %dopar% { # for (r in 1:R){
    
    library(dplyr)
    library(haven)
    library(fixest) 
    
    for (txvar in to_permute){
      df[,txvar] <- T0[[txvar]][,r]
    }
  
    #  If using fixest, set dataset to the argument "df" 
    if (any(grepl('feols', model, fixed=TRUE))){
      setFixest_estimation(data=df)
    }

    #  Estimate the model    
    lm0 <- eval(model) 
    
    #  Extract test stat
    if (stat=='b') {
      results$testdistribution[r,] <- lm0$coefficients[tx]
    }else if (stat=='t'){
      results$testdistribution[r,] <- lm0$coefficients[tx] / lm0$se[tx]
    }else{
      results$testdistribution[r,] <- eval(stat) # Allows passing arbitrary expressions
    }
  }
  
  #  Calculate p-value: two-sided test. 
  results$p.left <- mean(results$testdistribution < results$teststat)
  results$p.right <- mean(results$testdistribution > results$teststat) 
  results$p <- min(2*min(results$p.left, results$p.right),1)
    
  #  Release fixest data frame, if set
  if (any(grepl('feols', model, fixed=TRUE))) setFixest_estimation(reset=TRUE)
  
  #  Release parallel works, if set up under this script
  if (parallelize & !already.set.up) stopImplicitCluster() # stopCluster(cl)
  
  return(results)
}