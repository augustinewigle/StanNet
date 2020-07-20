#' Calculate log-posterior in R
#' @description Computes the log-posterior for a non-centred random effect model with normal likelihood using the identity link and default priors
#' @param testPars A named list of parameters to use
#' @param stanObj A StanNetRun object to supply the data
#' 
#' @return The value of the log-posterior distribution at the test parameters and data
#' @export
rprenormlogpost <-function(testPars, stanObj) {
  
  max.delta <- as.numeric(stanObj$model$max.delta)
  
  y <- stanObj$model$data$y # Matrix of response values; each row is a trial and each column is an arm
  t <- stanObj$model$data$t # Matrix of treatments; each row is a trial and each column is an arm
  se <- stanObj$model$data$se # Matrix of standard errors for each response; each row is a trial and each column is an arm
  ns <- stanObj$model$data$ns # Number of studies
  
  mu <- testPars$mu # Trial-specific baselines; length equal to the number of studies
  d2 <- testPars$d2 # Vector of relative difference parameters without the reference; length is number of treatments - 1
  sigma <- testPars$sigma # Extract the between-trial standard deviations (scalar)
  delta2raw <- testPars$delta2raw # Non-centered delta variable (vector of length equal to the number of studies)
  
  ds <- c(0, d2) # adding true difference for treatment 1 vs itself (0)
  
  # Initialize values
  mupart <-0
  dpart<-0
  deltaspart<-0
  normpart<-0
  sigmapart<- dunif(sigma,0,max.delta, log = T) # make sure this and all priors match the Stan file
  
  deltas<- matrix(0, nrow = ns, ncol = 2) # Store the trial-specific relative effects as a matrix
  
  for(i in 1:ns){ #set up delta2s
    
    for (j in 1:2) {
      
      if (j>1) {
        
        deltas[i,j] = delta2raw[i]*sigma + ds[t[i,j]]-ds[t[i,1]] # centring deltas variable
        
        deltaspart <-deltaspart + dnorm(delta2raw[i], log = T)
        
      }
      
    }
    
  }
  
  for(i in 1:length(mu)) {
    
    mupart <- mupart + dnorm(mu[i], 0, max.delta*15, log = T) # contribution in priors from mu
    
  }
  
  for (i in 1:length(d2)) {
    
    dpart <- dpart + dnorm(d2[i], 0, max.delta*15, log = T) # contribution in priors from d's (except the first d)
    
  }
  
  for(i in 1:ns){
    
    for(j in 1:2) {
      
      normpart <-normpart + dnorm(y[i,j],mu[i]+deltas[i,j], se[i,j], log = T) # Contribution from likelihood for study i, arm j
      
    }
    
  }
  
  sum(c(mupart, sigmapart,dpart, deltaspart, normpart)) # Sum up contributions
  
}
