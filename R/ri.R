ri <- function(df,model,tx,T0,stat='b',clusters=NULL){
  results <- list() # Container for results. 
  
  #  Extract number of repetitions. Depends on whether T0 is a list of matrices/dataframes or a single dataframe
  if (is.null(dim(T0))){
    K <- length(T0)
    R <- dim(T0[[1]])[2] # Use the first object in the list to determine number of reps.
  }else{
    R <- dim(T0)[2] # Number of replications to consider.
  }
  
  #  Extract number of models.  
  if (model[[1]]=='~'){
    M <- 1
  }else{
    M <- length(model)
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
          results$teststat <- stat(lm1) 
        }
      }
    }else{
      results$teststat <- stat(lm1) 
    }
  }else{
    restults$teststat <- stat(models1) 
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
    }else{
      restults$teststat <- stat(models0) 
    }
  }
  
  #  Calculate p-value: two-sided test. 
  results$p.left <- mean(results$testdistribution < results$teststat)
  results$p.right <- mean(results$testdistribution > results$teststat) 
  results$p <- min(2*min(results$p.left, results$p.right),1)
    
  return(results)
}