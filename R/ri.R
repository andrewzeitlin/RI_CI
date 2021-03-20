ri <- function(df,model,tx,T0,stat='b',clusters=NULL){
  results <- list() # Container for results. 
  
  R <- dim(T0)[2] # Number of replications to consider.
  
  
  #  Point estimate 
  lm1 <- lm_robust(model,df,clusters=clusters)
  if (!is.function(stat)){
    if (stat=='b') {
      results$teststat <- lm1$coefficients[tx]
    }else if (stat=='t'){
      results$teststat <- lm1$coefficients[tx] / lm1$std.error[tx]
    }
  }else{
    results$teststat <- stat(lm1) 
  }
  
  #  Distribution of test statistic under the null 
  results$testdistribution <- matrix(data=NA,nrow=R,ncol=length(tx))
  for (r in 1:R){
    df[,tx] <- T0[,r]
    lm0 <- lm_robust(model,df,clusters=clusters)
    if (!is.function(stat)){
      if (stat=='b') {
        results$testdistribution[r,] <- lm0$coefficients[tx]
      }else if (stat=='t'){
        results$testdistribution[r,] <- lm0$coefficients[tx] / lm0$std.error[tx]
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