#' Visualize individual treatment results.
#' @description Plots three graphs using data from a \code{net.tab()} object.
#' In the first graph, the average outcome or average events per person (depending on type of response data) is plotted against treatment.
#' In the second graph, the number of studies which examined each treatment is plotted.
#' In the third graph, the number of events relative to the total number of patients for each treatment across all studies.
#' Note that the third graph does not appear for continuous response data.
#' The number of comparisons made can be controlled with the cutoff and metric parameters.
#' Default value for these arguments are all ptreatments, and the number of patients (descending), respectively.
#' The type.outcome parameter must be specified, which is the same as given in the \code{net.tab} documentation.
#'
#' @param data The output of the \code{net.tab()} function applied to a dataset of interest. Must be specified.
#' @param cutoff The number of pairwise treatments to be plotted. Default is all.
#' If a different value is specified, the cutoff in combination with the specified metric is used to select treatments.
#' @param metric A string which specifies the metric used to determine which pairwise treatments to be included based on cutoff.
#' The default metric used is the maximum number of patients.
#' Any numeric column name of the net.tab object may be used for this parameter.
#' @param type.outcome A string which denotes the type of response variable examined in each study, similar to the net.tab specification.
#' Options are 'binary' (binary data with no follow-up), 'rate' (count data with person-time-at-risk), 'rate2' (binary data with follow-up time) and 'continuous'
#' @return Two or three ggplot objects.
#' @export
#'
#' @examples
#' \dontrun{
#' data("diabetes.sim")
#' data_diabetes <- data.prep(arm.data = diabetes.sim,
#' varname.t = "Treatment",
#' varname.s = "Study")
#' summary_diabetes <- net.tab(data = data_diabetes,
#'                           outcome = "diabetes", N = "n", type.outcome = "rate2", time = "followup")
#' intervention_viz(summary_diabetes, type.outcome = "rate2")
#'}
#'

intervention_viz <- function(data, cutoff, type.outcome, metric){
  # To deal with CRAN Check note
  treatment <- n.studies <- Outcomes <- Value <- Treatments <- NULL
  # Early stopping condition
  if(is.null(data$intervention))
    stop("Call the net.tab function on your data first.")

  # Extract only the comparison data frame
  data <- data$intervention
  # Convert to lower case column names 
  colnames(data) <- tolower(colnames(data))

  # Check if optional argument cutoff entered
  if(missing(cutoff)){cutoff <- nrow(data)}
  if(missing(metric)){metric <- 'n.patients'}
  if(metric %in% colnames(data) == F | is.numeric(data[[metric]]) == F)
    stop("Please select a metric from the numeric column names of the data")

  # Get the appropriate number of comparisons; update dataframe
  data <- data[order(data[[metric]], decreasing=TRUE),][1:cutoff,]

  # Determine the type of data
  if(type.outcome %in% c("bin","binom","binomial","binary", "cont", "continuous")){outcome <- 'av.outcome'}
  else if (type.outcome %in% c("rate2","rate")){outcome <- 'events.per.person'}

  # Create meaningful graph labels based on the response data
  if(type.outcome %in% c("cont", "continuous")){y_lab <- "Mean Treatment Response"}
  else if(type.outcome %in% c("bin","binom","binomial","binary")){y_lab <- "Proportion of Events"}
  else if(outcome == 'events.per.person'){y_lab <- "Average Events per Person"}

  # Average response for each treatment where response is dependent on distribution
  p1_data <- data.frame(data)
  p1 <- ggplot(p1_data, aes(y=p1_data %>% pull(outcome), x=treatment)) +
    geom_bar(stat="identity", color = "black", size=0.4, fill="pink") +
    coord_flip() + xlab("Treatment") + ylab(y_lab)

  # Summarize the number of studies examining each treatment
  p3 <- ggplot(data, aes(treatment, n.studies)) +
    geom_bar(stat="identity", color = "black", size=0.4, fill="darkolivegreen") +
    coord_flip() + xlab("Treatments") + ylab("Number of Studies")

  if(type.outcome %in% c("cont", "continuous") == FALSE){
    # Only create the last graph for data that is not continuous

    # Create a temporary data frame for compatibility with ggplot2
    results <- append(data$n.patients, data$n.events)
    labels <- append(rep("Total patients", cutoff),
                     rep("Observed outcomes", cutoff))
    temp_data <- data.frame(rep(data$treatment,2), labels, results, rep(data$n.studies,2))
    colnames(temp_data) <- c("Treatments", "Outcomes", "Value", "Studies")

  # Stacked bar-plot by pairwise treatments with observed outcomes and total patients; pooled data
  p2 <- ggplot(temp_data, aes(fill=forcats::fct_rev(Outcomes), y=Value, x=reorder(Treatments, -Value))) +
    geom_bar(stat="identity", color="black", size=0.40) + coord_flip() + xlab("Treatment") +
    ylab("Number of Events") +
    theme(legend.title = element_blank(), legend.position = "top")

  grid.arrange(p1, p3, p2, nrow=2)}

  else if(type.outcome %in% c("cont", "continuous") == TRUE){grid.arrange(p1, p3, nrow=2)}

}

