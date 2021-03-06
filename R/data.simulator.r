#' Generate simulated meta-analysis data for testing purposes.
#' @description Creates a raw data frame with specified number of trials, arms per trial, and features dependent on the type of response variable.
#' This raw data can then be passed through to the data.prep and `net.tab()` functions for modeling purposes.
#' The number of patients per trial are generated uniformly between the minimum and maximum numbers specified.
#' Note that the resultant network is not guaranteed to be connected, and this criteria should be verified before model fitting.
#'
#' @param nstudies The number of studies in the simulated data set. Must be specified.
#' @param response The type of response variable in each study; one of ("rate", "binary", "continuous", "rate2"), similar to `net.tab()`
#' "rate" data includes number of events, sample size, and person-time at risk. Number of events for treatment i are generated from a Poisson distribution with mean from the ith element of `count_mean`, and person-time at risk from an exponential distribution with mean `pt.rate`.
#' "binary" data includes number of events and sample size. The number of events are generated from a binomial distribution with treatment-dependent probability and sample size equal to the number of patients in each arm/trial.
#' "rate2" data includes number of events, sample size, and median follow-up time. The number of events are generated from a binomial distribution with treatment-dependent probability as in the binomial case, and follow-up time generated from a Uniform(min_time, max_time) distribution for all treatments.
#' "continuous" data includes mean response metric, sample size, and response metric standard error. Values for treatment i are simulated from Normal(cont_mean, cont_var) with corresponding standard errors generated uniformly from 0 to cont_error.
#' @param ntreatments The total number of possible unique treatments in the network, which are randomly sampled and assigned to treatment arms. Must be specified.
#' @param seed The seed value of the function for reproducibility (optional parameter).
#' @param min.patients The minimum number of patients per trial. Must be specified.
#' @param max.patients The maximum number of patients per trial. Must be specified.
#' @param event_p The probability vector of an event occurring per individual, used for binary and rate2 data. Default is 0.5, and elements must be greater than 0 and less than or equal to 1.
#' The vector length must be equal to the number of treatments, in which the ith element is the probability of outcome for treatment i.
#' @param exp_rate The mean of the exponential distribution for follow-up times for rate2 data, parameterized as 1/lambda. Default is 0.01.
#' Specified as a scalar; the follow-up times are assumed to be equal for both arms of the same trial.
#' @param cont_mean The mean of the response variable for 'continuous' response data. Default is 0.
#' Specified as a vector with length equal to the number of treatments, in which the ith element is the mean for treatment i.
#' @param cont_var The standard error of the response variable for 'continuous' response data. Default is 1.
#' Specified as a vector with length equal to the number of treatments, in which the ith element is the variance for treatment i.
#' @param cont_error The maximum standard error of the response metric for each trial. Default is 1.
#' Specified as a vector with length equal to the number of treatments, in which the ith element is the error for treatment i.
#' @param pt_rate The mean of the exponential distribution for time-at-risk per person for 'rate' data, parametrized as 1/lambda. Default is 0.01.
#' Specified as a scalar; the follow-up times per patient are generated to be the same for each arm within the same trial, and then get scaled by the number of patients in each arm.
#' @param count_mean The mean of the Poisson distribution used to generate count data for 'rate' response data. Counts are generated by multiplying the rate
#' parameter by the time-at-risk per person. Specified as a vector equal to the length of unique treatments.
#' @param min_time The minimum of the uniform distribution which will be used to generate follow-up times for 'rate2' data.
#' @param max_time The maximum of the uniform distribution which will be used to generate follow-up times for 'rate2' data.
#'
#' @return Data frame with studies, arms, sample sizes, and the parameters dependent on the selected response.
#' @export
#'
#' @examples
#' \dontrun{
#' data.simulator(nstudies=30, response="binary", ntreatments=5, seed=20,
#' min.patients=100, max.patients=500, 
#' event_p=c(0.05, 0.10, 0.15, 0.07, 0.09))
#'
#' data.simulator(nstudies=12, response="continuous", 
#' ntreatments=6, seed=13, min.patients=30, max.patients=200)
#'
#' data.simulator(nstudies=40, response="rate", ntreatments=5, seed=30,
#' min.patients=1000, max.patients=1200, 
#' pt_rate = 1000*0.5, count_mean = c(45, 35, 47, 40, 42))
#'
#' data.simulator(nstudies=20, response="rate2", ntreatments=6, seed=13,
#' min.patients=40, max.patients=500, min_time=3, max_time=5,
#' event_p = c(0.21, 0.5, 0.2, 0.3, 0.5, 0.55))
#' }


data.simulator <- function(nstudies, response, ntreatments, seed, min.patients, max.patients,
                           event_p, exp_rate, cont_mean, cont_var, cont_error, pt_rate, count_mean,
                           min_time, max_time){

  # Error messages and condition checks
  if(missing(nstudies))
    stop("The number of studies, nstudies, must be specified")
  if(missing(response))
    stop("response must be specified")
  if(response %in% c("binary", "rate", "rate2", "continuous") == F)
    stop("response selection must be binary, rate, rate2, or continuous")
  if(missing(ntreatments))
    stop("The number of treatments, ntreatments, must be specified")
  if(missing(min.patients))
    stop("The minimum number of patients per trial, min.patients, must be specified")
  if(missing(max.patients))
    stop("The maximum number of patients per trial, max.patients, must be specified")
  if(missing(seed)){
    print("Warning: No seed has been set in the function call")}
  if(min.patients > max.patients){
    stop("The minimum number of patients cannot exceed the maximum")}
  if(ntreatments > 2*nstudies){
    stop("The number of unique treatments cannot exceed the total number of study arms")
  }

  # Set the seed for consistency if specified
  if(!missing(seed)){set.seed(seed)}

  # Create the basic data frame
  study <- sort(rep(1:nstudies, 2)) # Create list of studies
  n <- length(study)
  treatment.labels <- LETTERS[1:ntreatments] # Create list of treatments
  treatment <- as.vector(replicate(nstudies, sample(treatment.labels, 2, replace=F))) # Assign treatments to each arm
  sampleSize <- ceiling(runif(n, min.patients, max.patients)) # Generate number of patients
  
  # Create the initial data frame
  df <- data.frame(study = as.factor(study), treatment = as.factor(treatment),
                          sampleSize = as.numeric(sampleSize))

  # Add data dependent on the response type specification

  # Binomial response data
  if(response == "binary"){

    # Add an events column (initialized as empty)
    df <- cbind(df, events = vector(length=n))

    # Check if parameters missing, if yes insert default values for probabilities
    if(missing(event_p)){event_p <- as.vector(rep(0.5, length(levels(df$treatment))))}
    if(FALSE %in% between(event_p,0,1))
      stop("The event probabilities must be valid probabilities")
    if(length(event_p) != ntreatments)
      stop("Event probability vector must be the same length as the number of treatments")

    # Need to generate the number of events associated with each unique treatment
    # These are specified with a binomial distribution with user-specified probabilities for each treatment
    for(i in 1:length(levels(df$treatment))){

      # Get the number of arms which used the ith unique treatment
      num.trt <- nrow(df[df$treatment == levels(df$treatment)[i],])
      # Get the corresponding sample sizes from each arm
      sizes <- df[df$treatment == levels(df$treatment)[i],]$sampleSize

      # Update the events for treatment i with a binomial distribution with size from each row of data frame
      # Success probability is the ith element of the event probability vector
      df[df$treatment == levels(df$treatment)[i],]$events <- as.vector(sapply(sizes,function(size) rbinom(1,size,event_p[i])))

    }
    if(0 %in% df$events)
      stop("Trials with 0 events have been simulated; change your simulation parameters")
  }

  # Rate2 response type includes number of events and median follow-up time
  else if(response == "rate2"){

    # Add an events and follow-up time column column
    df <- cbind(df, events = vector(length=n), follow.up = vector(length = n))

    # Check if parameters missing, if yes insert default
    if(missing(event_p)){event_p <- as.vector(rep(0.5, length(levels(df$treatment))))}
    if(missing(exp_rate)){exp_rate <- 0.01}
    # Stopping conditions
    if(length(exp_rate) != 1)
      stop("Exponential rate vector must be a scalar")
    if(length(event_p) != ntreatments)
      stop("Event probability vector must be the same length as the number of treatments")

    # Need to generate the follow-up times for each study
    df$follow.up <- as.vector(replicate(nstudies, rep(runif(1, min_time, max_time),2)))

    # Need to generate the number of events associated with each unique treatment
    for(i in 1:length(levels(df$treatment))){

      # Get the number of arms which used the ith unique treatment
      num.trt <- nrow(df[df$treatment == levels(df$treatment)[i],])
      # Get the corresponding sample sizes
      sizes <- df[df$treatment == levels(df$treatment)[i],]$sampleSize

      # Update the events for treatment i with a binomial distribution
      # Success probability is the ith element of the event probability vector
      df[df$treatment == levels(df$treatment)[i],]$events <- as.vector(sapply(sizes,function(size) rbinom(1,size,event_p[i])))

    }
    if(0 %in% df$events)
      stop("Trials with 0 events have been simulated; change your simulation parameters")
  }

  # Continuous data requires mean of response variable (i.e. treatment effect) and standard error of this estimate
  else if(response == "continuous"){

    # Add an events and followup time column column
    df <- cbind(df, response = vector(length=n), st.error = vector(length = n))

    # Check if parameters missing, if yes insert default
    if(missing(cont_mean)){cont_mean <- as.vector(rep(0, length(levels(df$treatment))))}
    if(missing(cont_var)){cont_var <- as.vector(rep(1, length(levels(df$treatment))))}
    if(missing(cont_error)){cont_error <- as.vector(rep(1, length(levels(df$treatment))))}
    # Check conditions
    if(length(cont_mean) != ntreatments)
      stop("Continuous response mean vector must be the same length as the number of treatments")
    if(length(cont_var) != ntreatments)
      stop("Continuous response variance vector must be the same length as the number of treatments")
    if(length(cont_error) != ntreatments)
      stop("Continuous response standard error vector must be the same length as the number of treatments")

    # Need to generate the number of events associated with each unique treatment
    for(i in 1:length(levels(df$treatment))){

      # Get the number of arms which used the ith unique treatment
      num.trt <- nrow(df[df$treatment == levels(df$treatment)[i],])

      # Update the continous response for treatment i with a normal distribution
      # The parameters are the ith value of the cont_mean and cont_var vectors, respectively
      df[df$treatment == levels(df$treatment)[i],]$response <- rnorm(num.trt, cont_mean[i], cont_var[i])

      # Update the standard errors associated with each of the response measurements
      # Parameters specified using the given vector
      df[df$treatment == levels(df$treatment)[i],]$st.error <- runif(num.trt, min=0, max=cont_error[i])
    }
  }

  # Rate data with number of events and person-time-at-risk
    else if(response == "rate"){

      # Add an events and person-time-at-risk column
      df <- cbind(df, events = vector(length=n), pt.risk = vector(length = n))

      # Check if parameters missing, if yes insert default
      if(missing(count_mean)){count_mean <- as.vector(rep(50, length(levels(df$treatment))))}
      if(missing(pt_rate)){pt_rate <- 100}

      # Check conditions
      if(length(count_mean) != ntreatments)
        stop("Count mean vector must be the same length as the number of treatments")
      if(!missing(pt_rate) & length(pt_rate) != 1)
        stop("Person-time-at-risk rate vector must be of length 1")

      # Need to generate the follow-up times for each study
      df$pt.risk <- as.vector(replicate(nstudies, rep(rexp(1, pt_rate),2))) * sampleSize

      # Need to generate the number of events associated with each unique treatment
      for(i in 1:length(levels(df$treatment))){

        # Get the number of arms which used the ith unique treatment
        num.trt <- nrow(df[df$treatment == levels(df$treatment)[i],])
        # Get the corresponding person-time-at-risk
        pt.risk <- df[df$treatment == levels(df$treatment)[i],]$pt.risk

        # Update the events for treatment i with a binomial distribution
        # Success probability is the ith element of the event probability vector
        df[df$treatment == levels(df$treatment)[i],]$events <- as.vector(sapply(pt.risk,function(pt.risk) rpois(1,pt.risk*count_mean[i])))
      }
      if(TRUE %in% c(df$events > df$sampleSize))
        stop("Check your simulation parameters; the observed events cannot exceed sample size")
      if(0 %in% df$events)
        stop("Trials with 0 events have been simulated; change your simulation parameters")
    }
  df # Print the data frame of results
}

