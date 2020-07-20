#' Visualize pairwise treatment comparisons
#' @description Plots two graphs from the output of a \code{net.tab()} function call.
#' In the first graph, the number of total patients and observed events are plotted for all treatment pairs in the network.
#' This data is pooled from all studies, and organized by each pairwise comparison among treatments of interest.
#' This graph is omitted for continous response data since the number of events is not meaningful for this response type.
#' In the second graph, the number of studies which examined each treatment pair (summed across the network) is plotted as a bar graph.
#' The number of comparisons made can be controlled with the cutoff and metric parameters, as in \code{intervention.viz()}
#' Default value for these arguments are all pairwise treatments, and the number of events (descending), respectively.
#'
#' @param data The output of the \code{net.tab()}function applied to a dataset of interest.
#' @param cutoff The number of pairwise treatments to be plotted. Default is all.
#' @param metric Metric used to determine which pairwise treatments to be included based on cutoff. The default 'cutoff' number of
#' treatments with the highest number of patients. Any numeric column name of the \code{net.tab()} object may be used, specified as a string.
#' @return One or two ggplot objects.
#' @export
#'
#' @examples
#' \dontrun{
#' data("diabetes.sim")
#' data_diabetes <- data.prep(arm.data = diabetes.sim, varname.t = "Treatment", varname.s = "Study")
#' summary_diabetes <- net.tab(data = data_diabetes,
#'                           outcome = "diabetes", N = "n", type.outcome = "rate2", time = "followup")
#' comparison.viz(summary_diabetes)
#' }
#'

comparison.viz <- function(data, cutoff, metric){
  #adding this line to get rid of CRAN check notes
  comparison <- n.studies <- Outcomes <- Value <- Treatments <- NULL
  
  # Early stopping condition
  if(is.null(data$comparison))
    stop("Call the net.tab function on your data first.")

  # Extract only the comparison data frame
  data <- data$comparison
  # Convert to lower case
  colnames(data) <- tolower(colnames(data))

  # Check if optional argument cutoff entered
  if(missing(cutoff)){cutoff <- nrow(data)}
  if(missing(metric)){metric <- 'n.patients'}
  if(metric %in% colnames(data) == F | is.numeric(data[[metric]]) == F)
    stop("Please select a metric from the numeric column names of the data")

  # Get the appropriate number of comparisons based on the cut-off and metric specified; update dataframe
  data <- data[order(data[[metric]], decreasing=TRUE),][1:cutoff,]

  # Summarize the number of studies comparing each treatment
  p2 <- ggplot(data, aes(comparison,n.studies)) +
    geom_bar(stat="identity", color = "black", size=0.4, fill="darkolivegreen") +
    coord_flip() + xlab("Pairwise Treatments") + ylab("Number of Studies")

  # Stacked bar-plot by pairwise treatments with observed outcomes and total patients; pooled data
  # This graph is not created for continuous response data
  if('n.outcomes' %in% colnames(data)){
    # Create a temporary data frame for compatibility with ggplot2
    results <- append(data$n.patients, data$n.outcomes)
    labels <- append(rep("Total patients", cutoff),
                     rep("Observed outcomes", cutoff))
    temp_data <- data.frame(rep(data$comparison,2), labels, results, rep(data$n.studies,2))
    colnames(temp_data) <- c("Treatments", "Outcomes", "Value", "Studies")

    # Create a stacked bar chart with the number of observed events relative to the number of patients
    # The values are aggregated accross all studies which compared each pair of treatments along the x-axis
    p1 <- ggplot(temp_data, aes(fill=forcats::fct_rev(Outcomes), y=Value, x=reorder(Treatments, -Value))) +
      geom_bar(stat="identity", color="black", size=0.40) + coord_flip() + xlab("Pairwise Treatments") +
      ylab("Number of Events") +
      theme(legend.title = element_blank(), legend.position = "top")
    grid.arrange(p1, p2, nrow=2)
    }

  else if('n.outcomes' %in% colnames(data) == F){p2}

}


