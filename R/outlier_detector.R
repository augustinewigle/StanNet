#' Detect outliers in a rainbow plot for arbitrary constant c
#' @description Takes the result of an \code{nma.fit} call and produces a list of trials and arms which are outliers.
#' Outliers are defined by lying outside the curve x^2+y=c for constant c.
#' By default, c is chosen to be 3 as in Dias et. al 2011, but can be specified to any other constant.
#'
#' @param obj An \code{nma.fit} object called upon an \code{nma.model}.
#' @param c The constant c used to define outliers in the curve x^2+y=c.
#' @return A data frame with the trials and arms that are poorly fit, along with the deviance residual and leverage.
#' If no outliers exist for the specified c, a printed message appears.
#' @export
#'
#' @examples
#' \dontrun{
#' data("thrombolytic")
#' data_thrombo <- data.prep(arm.data = thrombolytic, varname.t = "treatment", varname.s = "study")
#' fe_model <- nma.model(data = data_thrombo, outcome = "events", N = "sampleSize",
#' reference = "tPA", family = "binomial", link = "logit",
#' effects = "fixed")
#' fe_results <- nma.run(model = fe_model, n.adapt = 1000, n.burnin = 1000, n.iter = 10000)
#' diabetes.fit <- nma.fit(fe_results, main = "Consistency model")
#' outlier_detector(diabetes.fit)
#' }


outlier_detector <- function(obj, c){

  # Condition checks
  if(is.null(obj$leverage) | is.null(obj$w))
    stop("Call 'Stan.nma.fit' function on your model object first")
  if(missing(c)){c <- 3} # Use default from technical support manual (Dias et. al 2011)

  # Step 0: Do pre-processing on strings to extract corresponding trial and arm labels
  split_results <- strsplit(names(obj$pmdev), "[.]")
  trials <- c(sapply(split_results, "[", 2))
  arms <- c(sapply(split_results, "[", 3))

  # Step 1: Create a data frame with summary information
  x <- as.vector(unname(obj$w)) # x values on plot are w_ik = sign(resid) * residual_deviance
  y <- as.vector(unlist(obj$leverage)) # y values on plot are the leverage values
  temp_data <- data.frame(trials, arms, x, y)

  # Step 2: Check whether or not the data points lie outside x^2 + y = c for user-defined c
  p_vals <- mapply(function(x,y){if((x - 0)^2 + (y-c)>0) "Outside" else "Inside"},
                   x=temp_data$x, y=temp_data$y)

  # Summarize and clean up results table; output is the outlying points and values
  results <- subset(temp_data, p_vals == "Outside")
  colnames(results) <- c("Trial", "Arm", "w_ik", "leverage_ik")
  # Give printed message if no outliers detected
  if(nrow(results)>0){results}
  else print("No outliers detected for given c")
}
