#' Create a Stan model and data
#' @description Takes in data and builds appropriate Stan code representation of model
#'
#' @param data a BUGSnetData object produced from data.prep()
#' @param outcome A string indicating the column name of the outcome variable
#' @param N A string indciating the column name of the variable containing number of participants in each arm
#' @param sd A string indicating the column name of the standard deviation for each arm (NOT the standard error - convert by multiplying by square root of N) of the outcome (only required for continuous outcomes).
#' @param reference A string for the name of the treatment to be used as the "referent" comparator, and to be labelled as treatment 1 in the Stan code
#' @param type If `type`="inconsistency", an inconsistency model will be built. By default, type="consistency" and a consistency model is built. Inconsistency models are not yet supported in StanNet.
#' @param time A string (only required for binomial-cloglog or poisson-log models) indicating the name of variable indicating person-time followup (e.g person years) or study followup time.
#' @param family string indicating family of the distribution of the outcome. Options are: "binomial", "normal", "poisson".
#' @param link The link function for the nma model. Options are "logit" (binomial family), "cloglog" (binomial family), "log" (poisson family), "identity" (normal family).
#' @param effects A string indicating the type of treatment effect relative to the baseline, either "fixed" or "random"
#' @param prior.mu A string of Stan code which defines priors on the baseline treatment effects, or if "default", independent normal priors are used with mean 0 and standard deviation 15u, where u is the largest maximum likelihood estimator in single trials (see van Valkenhoef, G et al. 2012)
#' @param prior.d A string of Stan code that defines priors on the relative treatment effects. If "default", independent normal priors with mean 0 and standard deviation 15u where u is the largest maximum likelihood estimator in single trials (see van Valkenhoef, G et al. 2012) (BUGSnet)
#' @param prior.sigma A string of Stan code that defines the prior on the standard deviation of the relative treatment effects. By default, uniform distribution on (0, u) where u is largest MLE from single trials (see van Valkenhoef, G et al. 2012) (BUGSnet)
#'
#' @return An object of class StanNetModel which includes a string of Stan code which represent the model, and the model data, to be passed to runStan for compilation and MCMC sampling
#'
#' @export
makeStanmodel = function(data, outcome, N, sd = NULL, reference, type = "consistency", 
                         time = NULL, family, link, effects, prior.mu = "DEFAULT", 
                         prior.d = "DEFAULT", prior.sigma = "DEFAULT") {
  # To get rid of CRAN check note
  trt.ini <- trt <- trial <- trt.jags <- arm <-value <-variable <- NULL
  
  # First check that inputs are as expected, and set scale based on family and link

  if (class(data) != "BUGSnetData") 
    stop("'data' must be a valid BUGSnetData object created using the data.prep function.")
  if (family == "normal" & is.null(sd)) 
    stop("sd must be specified for continuous outcomes")
  if (family == "normal" & link != "identity") 
    stop("This combination of family and link is currently not supported in StanNet.")
  if (family == "poisson" & link != "log") 
    stop("This combination of family and link is currently not supported in StanNet.")
  if (family == "binomial" & !(link %in% c("logit", "cloglog"))) 
    stop("This combination of family and link is currently not supported in StanNet.")
  if (link == "logit" & family %in% c("binomial", 
                                      "binary", "bin", "binom")) {
    scale <- "Odds Ratio"
  }
  else if (link == "log" & family %in% c("binomial", 
                                         "binary", "bin", "binom")) {
    scale <- "Risk Ratio"
  }
  else if (link == "identity" & family == "normal") {
    scale <- "Mean Difference"
  }
  else if (link == "cloglog" & family %in% c("binomial", 
                                             "binary", "bin", "binom")) {
    if (is.null(time)) 
      stop("time must be specified when using a binomial family with the cloglog link")
    scale <- "Hazard Ratio"
  }
  else if (link == "log" & family == "poisson") {
    if (is.null(time)) 
      stop("time must be specified when using a poisson family with the log link")
    scale <- "Rate Ratio"
  }
  
  # Convert data into appropriate format
  
  data1 <- data$arm.data
  varlist <- c(trt = data$varname.t, trial = data$varname.s, 
               r1 = outcome, N = N, sd = sd, timevar = time)
  data1 <- data$arm.data[, varlist]
  names(data1) <- names(varlist)
  trt.key <- data1$trt %>% unique %>% sort %>% tibble(trt.ini = .) %>% 
    filter(trt.ini != reference) %>% add_row(trt.ini = reference, 
                                             .before = 1) %>% mutate(trt.jags = 1:dim(.)[1])
  data1 %<>% mutate(trt.jags = mapvalues(trt, from = trt.key$trt.ini, 
                                         to = trt.key$trt.jags) %>% as.integer)
  bugstemp <- data1 %>% arrange(trial, trt.jags) %>% group_by(trial) %>% 
    mutate(arm = row_number(), .name_repair = NULL) %>% ungroup() %>% select(-trt) %>% 
    gather("variable", "value", -trial, -arm) %>% 
    spread(arm, value)
  bugsdata2 <- list()
  
  # Add appropriate variable names
  
  for (v in unique(bugstemp$variable)) bugsdata2[[v]] <- as.matrix(bugstemp %>% 
                                                                     filter(variable == v) %>% select(-trial, -variable))
  names(bugsdata2)[names(bugsdata2) == "trt.jags"] <- "t"
  names(bugsdata2)[names(bugsdata2) == "N"] <- "n"
  names(bugsdata2)[names(bugsdata2) == "covariate"] <- "x"
  if (family == "binomial" && link %in% c("log", 
                                          "logit")) {
    names(bugsdata2)[names(bugsdata2) == "r1"] <- "r"
    bugsdata2 <- bugsdata2[names(bugsdata2) %in% c("ns", 
                                                   "nt", "na", "r", "n", "t", 
                                                   "x")]
  }
  else if (family == "normal" && link == "identity") {
    names(bugsdata2)[names(bugsdata2) == "r1"] <- "y"
    bugsdata2$se <- bugsdata2$sd/sqrt(bugsdata2$n)
    bugsdata2 <- bugsdata2[names(bugsdata2) %in% c("ns", 
                                                   "nt", "na", "y", "se", "t", 
                                                   "x")]
  }
  else if (family == "poisson" && link == "log") {
    names(bugsdata2)[names(bugsdata2) == "r1"] <- "r"
    names(bugsdata2)[names(bugsdata2) == "timevar"] <- "E"
    bugsdata2 <- bugsdata2[names(bugsdata2) %in% c("ns", 
                                                   "nt", "na", "r", "E", "t", 
                                                   "x")]
  }
  else if (family == "binomial" && link == "cloglog") {
    names(bugsdata2)[names(bugsdata2) == "r1"] <- "r"
    names(bugsdata2)[names(bugsdata2) == "timevar"] <- "time"
    bugsdata2 <- bugsdata2[names(bugsdata2) %in% c("ns", 
                                                   "nt", "na", "r", "n", "t", 
                                                   "x", "time")]
  }
  bugsdata2$nt <- data$treatments %>% nrow()
  bugsdata2$ns <- data$studies %>% nrow()
  # bugsdata2$na <- data$n.arms %>% select(n.arms) %>% t() %>% 
  #   as.vector
  bugsdata2 <- bugsdata2[names(bugsdata2) != "x"] # x is for metaregression - not implemented in this package yet
  
  # Create the treatment key to be printed after the model
  add.to.model <- trt.key %>% transmute(Treatments = paste0("// ", 
                                                            trt.jags, ": ", trt.ini, "\n")) %>% t() %>% 
    paste0() %>% paste(collapse = "") %>% paste0("\n\n// Treatment key\n", 
                                                 .)
  max.delta <- paste0(nma.prior(data, outcome = outcome, scale = scale, 
                                N = N, sd = sd, time = time)) # max.delta is the sd for prior distributions based on van Valkenhoef, G et al. 2012
  
  # Create strings for prior distributions
  
  # Prior for mu
  if (prior.mu == "DEFAULT") {
    prior.mu.str <- sprintf("mu ~ normal(0,%s*15)", 
                            max.delta)
  } else {
    
    prior.mu.str <- sprintf("mu ~ %s", prior.mu)
    
  }
  
  # Prior for true relative treatment differences
  if (prior.d == "DEFAULT") {
    if (type == "consistency") {
      prior.d.str <- sprintf("d2 ~ normal(0,%s*15)", 
                             max.delta)
    }
    else if (type == "inconsistency") {
      prior.d.str <- sprintf("for (c in 1:(nt-1)) {\n                           for (k in (c+1):nt)  { \n                           d[c,k] ~ dnorm(0,(%s*15)^(-2))\n                           } \n    }", 
                             max.delta)
      stop("Inconsistency models have not been inplemented in StanNet yet.")
    }
  } else {
    
    if (type == "consistency") {
      prior.d.str <- sprintf("d2 ~ %s", 
                             prior.d)
    }
    else if (type == "inconsistency") {
      prior.d.str <- sprintf("for (c in 1:(nt-1)) {\n                           for (k in (c+1):nt)  { \n                           d[c,k] ~ %s\n                           } \n  }", 
                             prior.d)
      stop("Inconsistency models have not been implemented in StanNet yet.")
    }
    
  }
  
  # Prior for sigma 
  if (prior.sigma == "DEFAULT") {
    prior.sigma.str <- sprintf("sigma ~ uniform(0,%s)", 
                                max.delta)
  }
  else {
    prior.sigma.str <- sprintf("sigma ~ %s", 
                                prior.sigma)
  }
  
  # Create the Stan code representation of the model and add the treatment key
  model <- makeStancode(family = family, link = link, effects = effects, 
                        inconsistency = (type == "inconsistency"), prior.mu.str, 
                        prior.d.str, prior.sigma.str) %>% 
    paste0(add.to.model)
  
  smodel <- structure(list(stan = model, data = bugsdata2, 
                           scale = scale, trt.key = trt.key, family = family, link = link, 
                           type = type, effects = effects, 
                           prior.mu = prior.mu, prior.d = prior.d, prior.sigma = prior.sigma, 
                          reference = reference, time = time,max.delta = max.delta, 
                           outcome = outcome, N = N, sd = sd), 
                      class = "StanNetModel")
  
  return(smodel)

}
