#' Compare the log-posterior of a model in R vs. Stan
#' @description Takes in a StanNetRun object and list of parameters and tests if the log-posteriors are equal (up to a constant) between the R and the Stan code
#' 
#' @param stanObj A StanNetRun object created by running `runStan` with `run = TRUE`
#' @param testPars A list of named lists of parameters to use for testing
#' @param detail Logical, TRUE if a numerical value indicating result of the test should be returned (default), or FALSE if the values of the log-posteriors and their differences should be returned
#' 
#' @return 3 vectors containing value of the log posterior in R and Stan and their differences, OR value which indicates the outcome of the comparison: 
#'         0 if the log posteriors were equal up to the same constant (may include instances where both R and Stan evaluated the log posterior to be infinite)
#'         -1 if the log posteriors were not equal up to the same constant
#'         -2 if the parameter settings resulted in too many infinite values and the differences could not be compared to check for equality
#' @export
compare_LP <-function(stanObj, testPars, detail = F) {
  
  # make sure that testPars is a list of lists
  if (!is.list(testPars[[1]])) {
    
    stop("testPars must be a list of lists to test if log-posteriors are equal up to a constant")
    
  } else {
    
    nTest <- length(testPars) # the number of sets of parameters to test
    lpR <- numeric(nTest)
    lpStan <- numeric(nTest)
    
  }
  
  family <- stanObj$family
  link <- stanObj$link
  effects <- stanObj$effects
  
  # Set up which function to use for the log-posterior in R
  if(family == "normal" && link == "identity") {
    
    if(effects == "random") {
      
      lpRfunc <- rprenormlogpost
      
    } else {
      
      lpRfunc <- fenormlogpost
      
    }
    
  } else if(family == "binomial" && link == "logit") {
    
    if(effects == "random") {
      
      lpRfunc <- rprebinom1logpost
      
    } else {
      
      lpRfunc <- febinom1logpost
      
    }
    
  } else if (family == "binomial" && link == "cloglog") {
    
    if(effects == "random") {
      
      lpRfunc <- rprebinom2logpost
      
    } else {
      
      lpRfunc <- febinom2logpost
      
    }
    
  } else if (family == "poisson" && link == "log") {
    
    if(effects == "random") {
      
      lpRfunc <- rprepoislogpost
      
    } else {
      
      lpRfunc <- fepoislogpost
      
    }
    
  }
  
  # calculate the log posterior in R and Stan
  
  for (i in 1:nTest) {
    
    upars <-unconstrain_pars(stanObj$fit ,testPars[[i]])
    lpStan[i] <-log_prob(object = stanObj$fit, upars = upars, adjust_transform = F)
    lpR[i] <- lpRfunc(testPars = testPars[[i]], stanObj = stanObj)
    
  }
  
  # Calculate vector of differences between R and Stan log-posterior
  
  diff <- lpR-lpStan # should be constant vector, with some NaNs if we have -Infs
 
  if (detail) { # For debugging - return the vector of differences and the vectors of the log posterior in R and Stan
    
    return(list(lpR, lpStan, diff))
    
  } else {
    
    tol <- .Machine$double.eps ^ 0.5 # tolerance
    
    numeric_diffs <-numeric(0) # put together the differences that are not NaN so we can check if they are all equal eventually
    
    for (j in 1:length(diff)) { # step through and compare vectors
      
      if (!is.nan(diff[j])) { # difference is numeric, so we should add it to the numeric vector
        
        numeric_diffs <- c(numeric_diffs, diff[j])
        
      } else if (is.infinite(lpR[j]) && is.infinite(lpStan[j])) { # 
        
        if (!(lpR[j]<0 && lpStan[j]<0)&&!(lpR[j]>0 && lpStan[j]>0)) {
          
          return(-1) # -1 means that the log posteriors were not equal at one or more values (in this case, one is Inf and one is -Inf)
          
        } # else they are both -Inf or Inf and we just do nothing
        
      }
      
    }
    
    if (length(numeric_diffs) <2) {
      
      return(-2) # There weren't enough non-Inf values to test if the differences are constant
      
    } else if(max(abs(diff(numeric_diffs)))<tol) {
      
      return(0) # max difference is below tolerance indicating it is equal for all parameter settings; return 0
      
    } else {
      
      return(-1) # max difference was above tolerance indicating LPs were not equal for 1 or more values; return -1
      
    }
    
  }

}
