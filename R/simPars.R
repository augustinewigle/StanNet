#' Simulate parameters for unit testing the log posteriors
#' @param ns Tne number of studies
#' @param nt The number of unique treatments
#' @param nsim The number parameter lists to simulate
#' @param effects A character string either "fixed" or "random" indicating what model to simulate parameters for
#' @param max.delta The parameter for the standard deviation of the distributions used to simulate parameters mu and d2
#' 
#' @return A list of named lists of parameters

simPars <- function(ns, nt, nsim = 10, effects, max.delta) {
  
  Pars <- list()
  
  for (i in 1:nsim) {
    
    if (effects == "fixed") {
      
      Pars[[i]] <- list(mu = rnorm(ns), d2 = rnorm(nt-1))
      
    } else if (effects == "random") {
      
      Pars[[i]] <- list(mu = rnorm(ns, 0, max.delta), d2 = rnorm(nt-1, 0, max.delta), 
                        sigma = runif(1,0,max.delta), delta2raw = rnorm(ns))
      
    } else {stop("effects must be either fixed or random")}
    
  }
  
  return(Pars)
  
}
