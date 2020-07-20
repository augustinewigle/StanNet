#' Compile and run Stan model
#' @description Takes data and Stan code created by makeStanmodel and compiles and/or runs the model using Stan
#' 
#' @param model A StanNetModel created by running makeStanmodel
#' @param pars A vector of the character strings giving names of parameters of interest. Default is relative treatment effects, and sigma if a random effects model is used
#' @param DIC Logical indicating if quantities required for calculating DIC will be monitored in MCMC sampling. Default is `TRUE`.
#' @param n.iter Gives the total number of iterations to run the MCMC chains. Default is 5000.
#' @param n.warmup Gives the number of warmup/burn-in iterations for the MCMC chains. Default is half of `n.iter`.
#' @param thin The period for saving samples. Default is 1 (samples will be saved for each MCMC sample).
#' @param n.chains The number of different MCMC chains to run. Default is 4.
#' @param inits Named lists of initial values for the MCMC chains. Default "random" which uses random values
#' @param compilation Logical indicating if compilation is needed, default is `TRUE`. It may be set to `FALSE` if a compiled Stan model from a previous call to runStan is supplied as `comp.model`.
#' @param run Logical indicating if the model should be run or not. Default is`TRUE`, may be set to `FALSE` if the model should only be compiled but not run through Stan.
#' @param verbose Logical indicating if Stan should print intermediate output while running the model. Default is TRUE and warnings and output will be printed
#' @param comp.model Optional argument supplying a previously compiled Stan model to be run if `compilation` is set to `FALSE`
#' @param testing Logical if the sampling should be done using the "Fixed_param" algorithm and single iterations in Stan, default is FALSE. Can be set to true for the purposes of testing the logposterior
#' 
#' @return A StanNetRun object, containing the following elements:
#'        fit - The stanfit object resulting from MCMC sampling
#'        model - The StanNetModel object
#'        comp.model - The compiled Stan model
#'        scale, family, link, effects - The scale, family, link, and type of effects from the StanNetModel
#'        trt.key - A key indicating the reference treatment and indexing the other treatments
#' 
#' @export
runStan <- function(model, pars = "DEFAULT", DIC = TRUE, n.iter = 5000, 
                    n.warmup = floor(n.iter/2), thin = 1, n.chains = 4, 
                    inits = 'random', compilation = TRUE, run = TRUE, 
                    verbose = TRUE, comp.model = NULL, testing = FALSE) {
  
  # Check model class and set inits if not specified
  if (class(model) != "StanNetModel") 
    stop("'model' must be a valid StanNetModel object created using the makeStanModel function.")
  
  if (compilation) {
    
    # Compile the model
    
    st_model <- stan_model(model_code = model$stan)
    
  } else if (!is.null(comp.model)){
    
    st_model <- comp.model # comp.model is a previously compiled model from a call to runStan
    
  } else {
    
    stop("If 'compilation' is FALSE, a compiled Stan model must be provided.")
    
  }
  
  if(run) {
    
    # set up initial values
    
    # if (!is.null(inits) && inits == "DEFAULT") {
    #   seeds <- sample(.Machine$integer.max, n.chains, replace = FALSE)
    #   inits <- list()
    #   for (i in 1:n.chains) inits[[i]] <- list(.RNG.seed = seeds[i], 
    #                                            .RNG.name = "base::Mersenne-Twister")
    # }
    
    # Set up parameters to monitor based on the type of model
    
    if (pars == "DEFAULT") {
      
      makepars <- "d2"
      
      if (model$effects == "random") {
        
        makepars <-c(makepars, "sigma")
        
      }
      
    } else {
      
      makepars <- pars
      
    }
    
    if (DIC) { # Include generated quantities needed to calculate pD, Dres, DIC
      
      if (model$family == "binomial") {
        DIC.pars <- c("dev", "totresdev", "rhat")
      }
      else if (model$family == "poisson") {
        DIC.pars <- c("dev", "totresdev", "theta")
      }
      else if (model$family == "normal") {
        DIC.pars <- c("theta", "totresdev", "dev")
      }
      
      makepars <- c(makepars, DIC.pars)
      
    }
    
    # Carry out sampling
    
    if(testing) { # Use the Fixed_param algorithm and only run a single chain with 1 iteration for calculated the log posterior
      
      samples <- sampling(object = st_model, data = model$data, pars = makepars, chains =1, 
                          iter = 1, thin = thin, init = inits, verbose = verbose,
                          algorithm = "Fixed_param")
      
    } else { # Do full sampling routine with specified iterations, chains, etc.
      
      samples <- sampling(object = st_model, data = model$data, pars = makepars, chains = n.chains, 
                          iter = n.iter, warmup = n.warmup, thin = thin, init = inits, verbose = verbose)
      
    }
    
  } else {
    
    samples <- NULL
    
  }
  
  run <- structure(list(fit = samples, model = model, comp.model = st_model,
                        scale = model$scale, family = model$family, link = model$link, effects = model$effects,
                        trt.key = as.character(t(model$trt.key[1]))), class = "StanNetRun")
  
  return(run)
  
}
