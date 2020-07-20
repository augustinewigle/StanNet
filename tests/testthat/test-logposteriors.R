context("Log posteriors in Stan vs R")
library(StanNet)

# Set up binary data
dat <- data.simulator(nstudies=20, response="binary", ntreatments=6, seed=14,
                      min.patients=40, max.patients=500, 
                      event_p=c(0.5, 0.25, 0.2, 0.5, 0.55, 0.7)) # binary data for binomial family with logit link

dat <- suppressWarnings(data.prep(dat, varname.t = "treatment", varname.s = "study"))

# Set up Fixed effects model and parameters
fe_model <-makeStanmodel(data = dat, outcome = "events", N = "sampleSize", 
                         reference = "A", family = "binomial", link = "logit", effects = "fixed")
fe_fit <- runStan(fe_model, testing = TRUE)
fe_Pars <- simPars(ns = 20, nt=6, effects = "fixed", max.delta = as.numeric(fe_model$max.delta))

# Set up random effects model and parameters
re_model <-makeStanmodel(data = dat, outcome = "events", N = "sampleSize", 
                         reference = "A", family = "binomial", link = "logit", effects = "random")
re_fit <- runStan(re_model, testing = TRUE)
re_Pars <- simPars(ns = 20, nt=6, effects = "random", max.delta = as.numeric(re_model$max.delta))

# Test the LPs
test_that("Log posterior in Stan and R are equal up to a constant for binomial data (without followup times)", {
  
  expect_equal(compare_LP(fe_fit, fe_Pars), 0)
  expect_equal(compare_LP(re_fit, re_Pars), 0)
  
})

# Set up continuous data
dat <- data.simulator(nstudies=25, response="continuous", ntreatments=7, seed=13, cont_mean = seq(1,1.4, length.out = 7),
                      min.patients=40, max.patients=500) # continuous data for normal family with id link
dat <- data.prep(dat, varname.t = "treatment", varname.s = "study")

# Set up fixed effects model and parameters
fe_model <-makeStanmodel(data = dat, outcome = "response", N = "sampleSize", sd = "st.error",
                         reference = "A", family = "normal", link = "identity", effects = "fixed")
fe_fit <- runStan(fe_model,testing = TRUE)
fe_Pars <- simPars(ns = 25, nt=7, effects = "fixed", max.delta = as.numeric(fe_model$max.delta))

# Set up random effects model and parameters
re_model <-makeStanmodel(data = dat, outcome = "response", N = "sampleSize", sd = "st.error",
                         reference = "A", family = "normal", link = "identity", effects = "random")
re_fit <- runStan(re_model, testing = TRUE)
re_Pars <- simPars(ns = 25, nt=7, effects = "random", max.delta = as.numeric(re_model$max.delta))

# Test LPs
test_that("Log posterior in Stan and R are equal up to a constant for continuous data", {

  expect_equal(compare_LP(fe_fit, fe_Pars), 0)
  expect_equal(compare_LP(re_fit, re_Pars), 0)

})

# Set up poisson data
dat <- data.simulator(nstudies=20, response="rate",
                      ntreatments=5, seed=4, min.patients=100, max.patients = 500,
                      pt_rate = 1, count_mean = seq(5,10, length.out = 5)/40)
dat <- data.prep(dat, varname.t = "treatment", varname.s = "study")

# Set up fixed effects model and parameters
fe_model <-makeStanmodel(data = dat, outcome = "events", N = "sampleSize", time = "pt.risk",
                         reference = "A", family = "poisson", link = "log", effects = "fixed")
fe_fit <- runStan(fe_model, testing = TRUE)
fe_Pars <- simPars(ns = 20, nt=5, effects = "fixed", max.delta = as.numeric(fe_model$max.delta))

# Set up random effects model and parameters
re_model <-makeStanmodel(data = dat, outcome = "events", N = "sampleSize", time = "pt.risk",
                         reference = "A", family = "poisson", link = "log", effects = "random")
re_fit <- runStan(re_model, testing = TRUE)
re_Pars <- simPars(ns = 20, nt=5, effects = "random", max.delta = as.numeric(re_model$max.delta))

# Test LP
test_that("Log posterior in Stan and R are equal up to a constant for rate (poisson) data", {

  expect_equal(compare_LP(fe_fit, fe_Pars), 0)
  expect_equal(compare_LP(re_fit, re_Pars), 0)

})

# Set up binary data with follow-up times
dat <- data.simulator(nstudies=30, response="rate2",
                      ntreatments=6, seed=4, min.patients=100, max.patients = 500,
                      event_p = seq(0.2,0.3, length.out = 6), min_time = 1, max_time = 4)
dat <- data.prep(dat, varname.t = "treatment", varname.s = "study")

# Set up fixed effects model and parameters
fe_model <-makeStanmodel(data = dat, outcome = "events", N = "sampleSize", time = "follow.up",
                         reference = "A", family = "binomial", link = "cloglog", effects = "fixed")
fe_fit <- runStan(fe_model, testing = TRUE)
fe_Pars <- simPars(ns = 30, nt=6, effects = "fixed", max.delta = as.numeric(fe_model$max.delta))

# Set up random effects model and parameters
re_model <-makeStanmodel(data = dat, outcome = "events", N = "sampleSize", time = "follow.up",
                         reference = "A", family = "binomial", link = "cloglog", effects = "random")
re_fit <- runStan(re_model, testing = TRUE)
re_Pars <- simPars(ns = 30, nt=6, effects = "random", max.delta = as.numeric(re_model$max.delta))
# Test LP
test_that("Log posterior in Stan and R are equal up to a constant for rate2 data (binomial data with median followup times)", {


  expect_equal(compare_LP(fe_fit, fe_Pars), 0)
  expect_equal(compare_LP(re_fit, re_Pars), 0)

})
