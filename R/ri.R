ri <- function(df,model,tx,T0,stat='b',clusters=NULL){
  results <- list() # Container for results. 
  
  #  Extract number of repetitions. Depends on whether T0 is a list of matrices/dataframes or a single dataframe
  if (is.null(dim(T0))){
    K <- length(T0)
    R <- dim(T0[[1]])[2] # Use the first object in the list to determine number of reps.
  }else{
    R <- dim(T0)[2] # Number of replications to consider.
  }
  
  #  Point estimate 
  lm1 <- lm_robust(model,df,clusters=clusters)
  if (K==1){
    if (!is.function(stat)){
      if (stat=='b') {
        results$teststat <- lm1$coefficients[tx]
      }else if (stat=='t'){
        results$teststat <- lm1$coefficients[tx] / lm1$std.error[tx]
      }else{
        results$teststat <- stat(lm1) 
      }
    }
  }else{
    results$teststat <- stat(lm1) 
  }
  
  #  Distribution of test statistic under the null 
  results$testdistribution <- matrix(data=NA,nrow=R,ncol=length(tx))
  for (r in 1:R){
    
    #  Replace treatment variables as needed
    if (K==1){
      df[,tx] <- T0[,r]
    }else{
      for (k in 1:K){
        df[,tx[k]] <- T0[[k]][,r]
      }
    }
        
    #  Estimate model
    lm0 <- lm_robust(model,df,clusters=clusters)
    
    #  Extract test statistic for this permutation
    if (K==1){
      if (!is.function(stat)){
        if (stat=='b') {
          results$testdistribution[r,] <- lm0$coefficients[tx]
        }else if (stat=='t'){
          results$testdistribution[r,] <- lm0$coefficients[tx] / lm0$std.error[tx]
        }
      }else{
        results$testdistribution[r,] <- stat(lm0)
      }
    }else{
      results$testdistribution[r,] <- stat(lm0)
    }
  }
  
  #  Calculate p-value: two-sided test. 
  results$p.left <- mean(results$testdistribution < results$teststat)
  results$p.right <- mean(results$testdistribution > results$teststat) 
  results$p <- min(2*min(results$p.left, results$p.right),1)
    
  return(results)
}