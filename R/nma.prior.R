#' Determine the appropriate value to be used for default priors 
#' @description Finds the largest maximum likelihood estimator from individual trials as in van Valkenhoef, G et al. 2012. This function was originally created in BUGSnet (https://github.com/audrey-b/BUGSnet) without documentation or comments. This version has small changes to avoid warnings during the CRAN check.
#' @param data.nma A data object created by running data.prep()
#' @param outcome A string for the column name containing the outcome
#' @param scale A string indicating the scale on which the relative treatment effects are being measures
#' @param N A string for the column name containing the sample sizes
#' @param sd A string for the column name containing the standard errors (only required for continuous outcomes)
#' @param time A string for the column names containing the follow-up times or person-time at risk for rate or rate2 data
#' 
#' @return The largest MLE to be used in determining default priors for mu, d and sigma
nma.prior <- function(data.nma, outcome, scale, N, sd=NULL, time = NULL){
  # Get rid of CRAN check note
  outcome.e <- N.e <- outcome.c <- N.c <- adj_r.e <- adj_r.c <- theta.e <- theta.c<-delta<-time.e <-time.c <- NULL
  
  if (scale =="Odds Ratio" ){
    type.outcome = "binomial"
  } else if (scale =="Risk Ratio"){
    type.outcome = "binomial"
  } else if (scale =="Mean Difference"){
    type.outcome = "continuous"
  } else if (scale =="Rate Ratio"){
    type.outcome = "rate"
  } else if (scale =="Hazard Ratio"){
    type.outcome = "rate2"
  }
  
  table <- suppressWarnings(by.comparison(data.nma, outcome, type.outcome = type.outcome, N, sd=sd, time = time))
  names(table)[names(table) == paste0(outcome,".e")] <- "outcome.e"
  names(table)[names(table) == paste0(outcome,".c")] <- "outcome.c"
  
  if (scale == "Odds Ratio"){
    names(table)[names(table) == paste0(N,".e")] <- "N.e"
    names(table)[names(table) == paste0(N,".c")] <- "N.c"
    
    deltas <- table %>% 
      mutate(adj_r.e = (outcome.e+0.5)/(N.e + 1), # add 0.5 to ensure ratio is non-zero
             adj_r.c = (outcome.c+0.5)/(N.c + 1)) %>%
      mutate(theta.e = log(adj_r.e/(1-adj_r.e)),
             theta.c = log(adj_r.c/(1-adj_r.c))) %>%
      mutate(delta = theta.e-theta.c) %>%
      select(delta)
    
  } else if (scale =="Risk Ratio"){
    names(table)[names(table) == paste0(N,".e")] <- "N.e"
    names(table)[names(table) == paste0(N,".c")] <- "N.c"
    
    deltas <- table %>% 
      mutate(adj_r.e = (outcome.e+0.5)/(N.e + 1),
             adj_r.c = (outcome.c+0.5)/(N.c + 1)) %>%
      mutate(theta.e = log(adj_r.e),
             theta.c = log(adj_r.c)) %>%
      mutate(delta = theta.e-theta.c) %>%
      select(delta)
    
  } else if (scale == "Mean Difference"){
    
    deltas <- table %>% 
      mutate(delta = as.numeric(outcome.e)-as.numeric(outcome.c)) %>%
      select(delta)
    
  } else if (scale == "Hazard Ratio"){
    names(table)[names(table) == paste0(N,".e")] <- "N.e"
    names(table)[names(table) == paste0(N,".c")] <- "N.c"
    names(table)[names(table) == paste0(time,".e")] <- "time.e"
    names(table)[names(table) == paste0(time,".c")] <- "time.c"
    
    deltas <- table %>% 
      mutate(adj_r.e = (outcome.e+0.5)/((N.e+1)*time.e),
             adj_r.c = (outcome.c+0.5)/((N.c+1)*time.c)) %>%
      mutate(theta.e = log(-log(adj_r.e)),
             theta.c = log(-log(adj_r.c))) %>%
      mutate(delta = theta.e-theta.c) %>%
      select(delta)
    
  } else if (scale =="Rate Ratio"){
    names(table)[names(table) == paste0(time,".e")] <- "time.e"
    names(table)[names(table) == paste0(time,".c")] <- "time.c"
    
    deltas <- table %>% 
      mutate(adj_r.e = (outcome.e+0.5)/(time.e+1),
             adj_r.c = (outcome.c+0.5)/(time.c+1)) %>%
      mutate(theta.e = log(adj_r.e),
             theta.c = log(adj_r.c)) %>%
      mutate(delta = theta.e-theta.c) %>%
      select(delta)
    
  }
  
  return(max(abs(deltas)))
}

