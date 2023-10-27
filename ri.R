ri <- function(df,model,tx,T0,stat='b',clusters=NULL, parallelize = FALSE, parallel.cores = NULL){
  results <- list() # Container for results. 
  
  #  Generically required packages
  require(estimatr)
  
  #  Initialize parallel processing
  if (parallelize){
    
    #  Require packages used for parallel implementation
    require(foreach) 
    require(doParallel)
    
    #  If number of cores to use not specified, default to *all*
    if (is.null(parallel.cores)){
      parallel.cores <- detectCores()
    }
    cl <- makeCluster(parallel.cores)
    registerDoParallel(cl)
  } else{
    registerDoSEQ()
  }
  
  
  #  Extract number of repetitions. Depends on whether T0 is a list of matrices/dataframes or a single dataframe
  if (is.null(dim(T0))){
    K <- length(T0)
    R <- dim(T0[[1]])[2] # Use the first object in the list to determine number of reps.
  }else{
    K <- length(tx)
    R <- dim(T0)[2] # Number of replications to consider.
  }
  
  #  Extract number of models and associated dependent variable(s)
  if (model[[1]]=='~'){
    M <- 1
    # Extract dependent variable name.
    outcome <- as.character(model[[2]])  
  }else{
    M <- length(model)
    
    #  Dependent variable(s)
    for (m in 1:M){
      outcome[m] <-model[m][[2]]
    }
  }
  
  
  #  Point estimate 
  if (M==1){
    lm1 <- lm_robust(model,df,clusters=clusters)
  }else{
    models1 <- vector(mode='list',length=M)
    for (m in 1:M){
      models1[[m]] <- lm_robust(model[m],df,clusters=clusters)
    }
  }
  
  
  #  Test statistic under the realized treatment
  if (M==1){
    if (K==1){
      if (!is.function(stat)){
        if (stat=='b') {
          results$teststat <- lm1$coefficients[tx]
        }else if (stat=='t'){
          results$teststat <- lm1$coefficients[tx] / lm1$std.error[tx]
        }else{
          #  Pass lm results to stat() function. 
          #  If two arguments required, second is assumed to be the data frame 
          #  (this is used in residualized KS stat).
          if (length(formals(stat))==1){
            results$teststat <- stat(lm1) 
          }else{
            results$teststat <- stat(lm1,df) 
          }
        }
      }
    }else{
      if (length(formals(stat))==1){
        results$teststat <- stat(lm1) 
      }else{
        results$teststat <- stat(lm1,df)
      }
    }
  }else{
    if (length(formals(stat))==1){
      results$teststat <- stat(models1) 
    }else{
      results$teststat <- stat(models1,df) 
    }
  }
  
  #  Distribution of test statistic under the null 
  results$testdistribution <- foreach(
    r = 1:R,
    .packages = c('estimatr') 
    ) %dopar% {
    
    #  Replace treatment variables as needed
    if (K==1){
      df[,tx] <- T0[,r]
    }else{
      for (k in 1:K){
        df[,tx[k]] <- T0[[k]][,r]
      }
    }
        
    #  Estimate model
    if (M==1){
      lm0 <- lm_robust(model,df,clusters=clusters)
    }else{
      models0 <- vector(mode='list',length=M)
      for (m in 1:M){
        models0[[m]] <- lm_robust(model[m],df,clusters=clusters)
      }
    }
    
    #  Extract test statistic for this permutation
    if (M==1){
      if (K==1){
        if (!is.function(stat)){
          if (stat=='b') {
            # results$testdistribution[r,] 
            this_result <- lm0$coefficients[tx]
          }else if (stat=='t'){
            # results$testdistribution[r,] 
            this_result <- lm0$coefficients[tx] / lm0$std.error[tx]
          }
        }else{
          if (length(formals(stat))==1){
            #results$testdistribution[r,] 
            this_result <- stat(lm0)
          }else{
            #results$testdistribution[r,] 
            this_result <- stat(lm0,df)
          }
        }
      }else{
        if (length(formals(stat))==1){
          # results$testdistribution[r,] 
          this_result <- stat(lm0)
        }else{
          # results$testdistribution[r,] 
          this_result <- stat(lm0,df)    
        }
      }
    }else{
      #restults$teststat <- 
      this_result <- stat(models0) 
    }
      
    this_result
  }
  
  
  #  Close parallel session
  if (parallelize) stopCluster(cl)
  
  
  #  Calculate p-value: two-sided test. 
  results$p.left <- mean(results$testdistribution < results$teststat)
  results$p.right <- mean(results$testdistribution > results$teststat) 
  results$p <- min(2*min(results$p.left, results$p.right),1)
  
  
  return(results)
}