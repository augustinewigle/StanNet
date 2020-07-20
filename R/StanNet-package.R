#' @name StanNet
#' @title StanNet: Bayesian Network Meta-Analysis Using Stan
#'
#' @description StanNet is an R package to conduct Bayesian network meta-analyses in compliance with best practice and reporting guidelines.
#' Bayesian analyses are conducted with Stan code generated automatically based on the user's inputs. Outputs are highly customizable and include plots for data visualization and exploration and extensive visualizations for NMA results including rainbow plots, SUCRA plots, rankograms, forest plots, and league tables and heatplots.
#'
#' To use this package, you must have Stan installed on your computer.
#'
#' @import meta
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import forcats
#' @importFrom BUGSnet data.prep net.tab nma.fit nma.model nma.run
#' @importFrom purrr map
#' @importFrom stats density median offset update dbinom dnorm dpois dunif qnorm runif rbinom reorder rexp rnorm rpois
#' @importFrom scales pretty_breaks rescale
#' @importFrom rlang quo
#' @importFrom tibble tibble
#' @importFrom plyr mapvalues
#' @importFrom utils combn globalVariables data
#' @importFrom graphics lines par plot points text
#' @importFrom magrittr %>% %<>%
#' @importFrom igraph incident is.connected clusters V E graph_from_data_frame ecount vcount layout_as_bipartite layout_as_star layout_as_tree layout_in_circle layout_nicely layout_on_grid layout_on_sphere layout_randomly layout_with_dh layout_with_fr layout_with_gem layout_with_graphopt layout_with_kk layout_with_lgl layout_with_mds layout_with_sugiyama
#' @importFrom gridExtra grid.arrange
#' @importFrom Rdpack reprompt
#' @importFrom stringr str_detect
#' @importFrom grDevices colorRampPalette dev.off
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#' @importFrom rstan stan_model sampling unconstrain_pars log_prob As.mcmc.list

if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))
