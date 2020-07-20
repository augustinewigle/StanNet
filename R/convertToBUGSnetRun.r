#' Convert a StanNetRun object to a BUGSnetRun object
#' @description Takes a StanNetRun object and converts the samples object to something of the same form as the samples object from nma.run,
#'              to be compatible with nma.fit and other plotting functions from BUGSnet
#' @param StanNetobj A StanNetRun object to be converted to a BUGSnetRun object
#' 
#' @return A BUGSnetRun object
#' 
#' @export
convertToBUGSnetRun <- function(StanNetobj) {
  
  sample_obj <- As.mcmc.list(StanNetobj$fit) # convert to mcmc.list format
  nChains <- length(sample_obj)
  nRows <- nrow(sample_obj[[1]])
  
  nt <- StanNetobj$model$data$nt # number of treatments
  ns <- StanNetobj$model$data$ns # number of studies
  
  # Add d[1] column of zeros to every part of the list, and change the column names
  for(i in 1:nChains) {
    
    sample_obj[[i]] <- cbind(rep(0, nRows), sample_obj[[i]])
    
    colnames(sample_obj[[i]])[1:nt] <- paste0("d[", seq(1,nt), "]")
    
    nCols <-ncol(sample_obj[[i]])
    
    # Now need to add the data for calculating residual deviances, depending on what family is
    # Order works like n[1,1]...n[ns,1]n[1,2]...n[ns,2] since only 2 treatments, where ns is the number of studies
    # We need to iterate through the studies followed by iterating through the arms for compatibility with subsequent diagnostic functions
    
    
    # The data depends on the response type - categorize into the 4 types and pull the corresponding data associated with each response variable
    if (StanNetobj$family == "normal") {
      
      # then add y (response variable) and precision (variation associated with response measurement)
      for(k in 1:2) {
        
        for (j in 1:ns) {
          
          sample_obj[[i]] <- cbind(sample_obj[[i]], rep(StanNetobj$model$data$y[j,k], nRows))
          colnames(sample_obj[[i]])[(nCols+1)] <- paste0("y[", j, ",", k, "]")
          nCols <- nCols+1
          
        }
        
      }
      
      # Add the precision to the data 
      for(k in 1:2) {
        
        for (j in 1:ns) {
          
          sample_obj[[i]] <- cbind(sample_obj[[i]], rep(StanNetobj$model$data$se[j,k]^(-2), nRows)) # Precision
          colnames(sample_obj[[i]])[(nCols+1)] <- paste0("prec[", j, ",", k, "]")
          nCols <- nCols+1
          
        }
        
      }
      
    } else if (StanNetobj$family == "poisson") {
      
      # Add the corresponding counts to the data frame
      for(k in 1:2) {
        
        for (j in 1:ns) {
          
          sample_obj[[i]] <- cbind(sample_obj[[i]], rep(StanNetobj$model$data$r[j,k], nRows))
          colnames(sample_obj[[i]])[(nCols+1)] <- paste0("r[", j, ",", k, "]")
          nCols <- nCols+1
          
        }
        
      }
      
    } else if (StanNetobj$family == "binomial") {
      
      # Add the counts and and sample size 
      # This order is required (despite being less efficient)
      for(k in 1:2) {
        
        for (j in 1:ns) {
          
          sample_obj[[i]] <- cbind(sample_obj[[i]], rep(StanNetobj$model$data$r[j,k], nRows))
          colnames(sample_obj[[i]])[(nCols+1)] <- paste0("r[", j, ",", k, "]")
          nCols <- nCols+1
          
        }
        
      }
      
      for(k in 1:2) {
        
        for (j in 1:ns) {
          
          sample_obj[[i]] <- cbind(sample_obj[[i]], rep(StanNetobj$model$data$n[j,k], nRows)) # Sample sizes
          colnames(sample_obj[[i]])[(nCols+1)] <- paste0("n[", j, ",", k, "]")
          nCols <- nCols+1
          
        }
        
      }
      
    }
    
  }
  
  # Combine all the results as required, and convert to BUGSnetRun object
  bugsobj <- structure(list(samples = sample_obj, model = StanNetobj$model, family = StanNetobj$family, 
                  link = StanNetobj$link, scale = StanNetobj$scale, trt.key =StanNetobj$trt.key), class = "BUGSnetRun")
  
  return(bugsobj)
}
