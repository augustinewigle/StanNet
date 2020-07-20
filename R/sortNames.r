#' Sort column names by study number, then by arm number
#' @description Creates a rearranged vector of column names to be used to reorder columns for compatibility between Stan and JAGS sampling output
#' @param string A character string of the prefix for the desired column names
#' @param ns The number of studies
#' This function is used in convertToBUGSnetRun to ensure the ordering of the data is the same between Stan and JAGS output
#' @export
sortNames <- function(string, ns) {
  
  # Initialize vector and counter
  sortednames <- vector(length = 2*ns)
  count <-1
  for (i in 1:2) { # Iterate through trial arms
    
    for (j in 1:ns) { # Iterate through studies
      
      sortednames[count] <- paste0(sprintf("%s.", string), j, ".", i, ".") # Reorganize labels to iterate through studies before arms 
      count <- count +1 # Update counter
      
    }
    
  }
  
  return(sortednames) # Return a vector with labels
  
}
