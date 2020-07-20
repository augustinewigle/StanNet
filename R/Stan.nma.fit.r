#' Assess NMA model fit
#' @description Computes the Deviance Information Criteria, effective number of parameters, and posterior mean of the residual deviance for a given model, and produces a rainbow plot (called a leverage plot in the NICE Technical Support Document 2). Points outside of the purple curve can be said to contribute to the model's poor fit. 
#' @param nma A BUGSnetRun object produced by running `nma.run()` or `convertToBUGSnetRun()`
#' @param plot.pD Whether to include the effective number of parameters on the plot, default is `TRUE`
#' @param plot.DIC Whether to include the DIC on the plot, default is `TRUE`
#' @param plot.Dres Whether to include the posterior mean of the residual deviance on the plot, default is `TRUE`
#' @param ... Graphical arguments such as `main=` to be passed to `plot()`
#' @return `DIC` - The Deviance Information Criteria
#' `leverage` - A vector of leverages calculated for each data point (one per arm per study)
#' `w` - A vector of each data point's contribution to the posterior mean deviance
#' `pmdev` - A vector of the posterior mean residual deviance for each data point
#' `Dres` - The posterior mean of the residual deviance
#' `pD` - The effective number of parameters, calculated as the sum of `leverage`
#' 
#' @export
Stan.nma.fit <-function (nma, plot.pD = TRUE, plot.DIC = TRUE, plot.Dres = TRUE,
                    ...) {
  if (class(nma) != "BUGSnetRun")
    stop("'nma' must be a valid BUGSnetRun object created using the nma.run function.")
  jagssamples <- nma$samples
  if (class(jagssamples) != "mcmc.list") {
    stop("Object jagssamples must be of class mcmc.list")
  }
  ns <- nma$model$data$ns
  samples = do.call(rbind, jagssamples) %>% data.frame()
  dev <- samples %>% select(., starts_with("dev"))
  dev <-dev[sortNames("dev", ns)]
  totresdev <- samples$totresdev %>% mean()
  pmdev <- colMeans(dev)

  if (nma$family == "binomial") {
    rhat <- samples %>% select(., starts_with("rhat"))
    rhat <- rhat[sortNames("rhat", ns)] # reorder to sort by trials then arms

    r <- samples %>% select(., starts_with("r."))
    n <- samples %>% select(., starts_with("n"))
    rtilde <- rhat %>% colMeans() %>% matrix(nrow = 1, ncol = ncol(rhat)) %>%
      data.frame() %>% slice(rep(1:n(), each = nrow(rhat)))
    pmdev_fitted <- 2 * (r * log(r/rtilde) + (n - r) * log((n -
                                                              r)/(n - rtilde)))[1, ]
    if (TRUE %in% c(nma$model$data$r == 0))
      warning("Leverage cannot be calculated for zero cells.")
  }
  else if (nma$family == "poisson") {
    rhat <- samples %>% select(., starts_with("theta"))
    rhat <- rhat[sortNames("theta", ns)]
    r <- samples %>% select(., starts_with("r."))
    n <- samples %>% select(., starts_with("n"))
    thetatilde <- rhat %>% colMeans() %>% matrix(nrow = 1,
                                                 ncol = ncol(rhat)) %>% data.frame() %>% slice(rep(1:n(),
                                                                                                   each = nrow(rhat)))
    pmdev_fitted <- 2 * ((thetatilde - r) + r * log(r/thetatilde))[1, ]
    if (TRUE %in% c(nma$model$data$r == 0))
      warning("Leverage cannot be calculated for zero cells.")
  }
  else if (nma$family == "normal") {
    rhat <- samples %>% select(., starts_with("theta"))
    rhat <- rhat[sortNames("theta", ns)]
    r <- samples %>% select(., starts_with("y"))
    prec <- samples %>% select(., starts_with("prec"))
    prec <- prec[sortNames("prec", ns)]
    ytilde <- rhat %>% colMeans() %>% matrix(nrow = 1, ncol = ncol(rhat)) %>%
      data.frame() %>% slice(rep(1:n(), each = nrow(rhat)))
    pmdev_fitted <- ((r - ytilde) * (r - ytilde) * prec)[1,
                                                         ]
  }
  leverage = pmdev - pmdev_fitted
  DIC = sum(leverage) + totresdev
  sign = sign(colMeans(r) - colMeans(rhat))
  w = sign * sqrt(as.numeric(pmdev))
  pD = sum(leverage)
  eq = function(x, c) {
    c - x^2
  }
  x = seq(-3, 3, 0.001)
  c1 = eq(x, c = rep(1, 6001))
  plot(w, leverage, xlab = expression("w"[ik]), ylab = expression("leverage"[ik]),
       ylim = c(0, max(c1 + 3, na.rm = TRUE) * 1.15), xlim = c(-3,
                                                               3), ...)
  points(x, ifelse(c1 < 0, NA, c1), lty = 1, col = "firebrick3",
         type = "l")
  points(x, ifelse(c1 < -1, NA, c1) + 1, lty = 2, col = "chartreuse4",
         type = "l")
  points(x, ifelse(c1 < -2, NA, c1) + 2, lty = 3, col = "mediumpurple3",
         type = "l")
  points(x, ifelse(c1 < -3, NA, c1) + 3, lty = 4, col = "deepskyblue3",
         type = "l")
  if (plot.pD == TRUE) {
    text(2, 4.3, paste("pD=", round(pD, 2)), cex = 0.8)
  }
  if (plot.Dres == TRUE) {
    text(2, 3.9, paste("Dres=", round(totresdev, 2)),
         cex = 0.8)
  }
  if (plot.DIC == TRUE) {
    text(2, 3.5, paste("DIC=", round(DIC, 2)), cex = 0.8)
  }
  return(list(DIC = DIC, leverage = leverage, w = w, pmdev = pmdev,
              Dres = totresdev, pD = pD))
}
