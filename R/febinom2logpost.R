#' Calculate log-posterior in R
#' @description Computes the log-posterior for a fixed effect model with binomial likelihood using the cloglog link and default priors
#' @param testPars A named list of parameters to use
#' @param stanObj A StanNetRun object to supply the data
#' 
#' @return The value of the log-posterior distribution at the test parameters and data
#' @export
febinom2logpost <-function(testPars, stanObj) {
  
  max.delta <- as.numeric(stanObj$model$max.delta)
  
  r <- stanObj$model$data$r # Matrix of observed counts; each row is a trial and each column is an arm
  n <- stanObj$model$data$n # Matrix of sample sizes; each row is a trial and each column is an arm
  time <- stanObj$model$data$time # Matrix of follow-up times; each row is a trial and each column is an arm
  t <- stanObj$model$data$t # Matrix of treatments; each row is a trial and each column is an arm
  ns <- stanObj$model$data$ns # Number of studies
  
  d2 <- testPars$d2 # Vector of relative difference parameters without the reference; length is number of treatments - 1
  mu <- testPars$mu # Trial-specific baselines; length equal to the number of studies
  
  ds <- c(0, d2) # adding true difference for treatment 1 vs itself (0)
  
  # Initialize values
  mupart <-0
  dpart<-0
  binompart<-0
  
  for(i in 1:length(mu)) {
    
    mupart <- mupart + dnorm(mu[i], 0, max.delta*15, log = T) # contribution in priors from mu
    
  }
  
  for (i in 1:length(d2)) {
    
    dpart <- dpart + dnorm(d2[i], 0, max.delta*15, log = T) # contribution in priors from d's (except the first d)
    
  }
  
  for(i in 1:ns){
    
    for(j in 1:2) {
      
      q <-log(time[i,j]) + mu[i]+ds[t[i,j]] - ds[t[i,1]]
      p <- 1-exp(-exp(q))
      
      binompart <- binompart + dbinom(r[i,j], n[i,j], p, log = T) # Add contribution from likelihood, study i arm j
      
    }
    
  }
  
  sum(c(mupart, dpart, binompart)) # Combine results
  
}
