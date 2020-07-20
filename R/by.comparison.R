#' Reorder data.frame for use in nma.prior
#' @description This function was originally an internal function to BUGSnet (https://github.com/audrey-b/BUGSnet). This version was altered to avoid warnings during the CRAN check but is otherwise exactly as presented in BUGSnet. 
#' @param data.nma A data object
#' @param outcome String indicating column name of outcome variable
#' @param type.outcome String indicating the type of outcome
#' @param N String indicating column name of sample size
#' @param sd A string for the column name containing the standard errors (only required for continuous outcomes)
#' @param time A string for the column names containing the follow-up times or person-time at risk for rate or rate2 data
#' 
#' @return The rearranged data
by.comparison <- function(data.nma, outcome, type.outcome="binomial", N, sd=NULL, time = NULL){
  #adding this line to get rid of CRAN check notes
  trial <- trt <- treatments <- comparisons <- V1 <- V2 <- trt.e <- trt.c <- NULL
  
  data <- data.nma$arm.data %>% select(outcome, data.nma$varname.t, data.nma$varname.s, N, sd, time)
  names(data)[names(data) == data.nma$varname.t] <- "trt"
  names(data)[names(data) == data.nma$varname.s] <- "trial"
  
  if (type.outcome=="continuous"){
    
    names(data)[names(data) == data.nma$sd] <- "sd"
    
    data %<>% select(trial, trt, outcome, N, sd)
    data.st <- select(data, trial, trt)
    data.st %<>% nest(treatments=c(trt))
    data.st %<>%
      mutate(comparisons = map(treatments, function(x) combn(x[[1]], 2) %>% t() %>% as_tibble)) %>%
      select(-treatments) %>%
      unnest(cols = c(comparisons)) %>%
      rename(trt.e=V1, trt.c=V2) %>%
      left_join(data %>% select(trial, outcome, N, trt, sd), by = c("trial", "trt.e" = "trt")) %>%
      left_join(data %>% select(trial, outcome, N, trt, sd), by = c("trial", "trt.c" = "trt"), suffix=c(".e",".c"))%>%
      mutate(comparison = ifelse(trt.e < trt.c, 
                                 paste(trt.e, trt.c, sep = " vs. "),
                                 paste(trt.c, trt.e, sep = " vs. ")))
    
  } else if (type.outcome=="binomial"){
    
    data %<>% select(trial, trt, outcome, N)
    data.st <- select(data, trial, trt) 
    data.st %<>% nest(treatments=c(trt))
    data.st %<>% 
      mutate(comparisons = map(treatments, function(x) combn(x[[1]], 2) %>% t() %>% as_tibble)) %>%
      select(-treatments) %>% 
      unnest(cols = c(comparisons)) %>%
      rename(trt.e=V1, trt.c=V2) %>%
      left_join(data %>% select(trial, outcome, N, trt), by = c("trial", "trt.e" = "trt")) %>%
      left_join(data %>% select(trial, outcome, N, trt), by = c("trial", "trt.c" = "trt"), suffix=c(".e",".c"))%>%
      mutate(comparison = ifelse(trt.e < trt.c, 
                                 paste(trt.e, trt.c, sep = " vs. "),
                                 paste(trt.c, trt.e, sep = " vs. "))) 
    
  } else if (type.outcome %in% c("rate", "rate2")){
    
    names(data)[names(data) == data.nma$time] <- "time"
    
    data %<>% select(trial, trt, outcome, N, time)
    data.st <- select(data, trial, trt)
    data.st %<>% nest(treatments=c(trt))
    data.st %<>%
      mutate(comparisons = map(treatments, function(x) combn(x[[1]], 2) %>% t() %>% as_tibble)) %>%
      select(-treatments) %>%
      unnest(cols = c(comparisons)) %>%
      rename(trt.e=V1, trt.c=V2) %>%
      left_join(data %>% select(trial, outcome, N, trt, time), by = c("trial", "trt.e" = "trt")) %>%
      left_join(data %>% select(trial, outcome, N, trt, time), by = c("trial", "trt.c" = "trt"), suffix=c(".e",".c"))%>%
      mutate(comparison = ifelse(trt.e < trt.c, 
                                 paste(trt.e, trt.c, sep = " vs. "),
                                 paste(trt.c, trt.e, sep = " vs. "))) 
  }
  return(data.st)
}
