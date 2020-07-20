context("Checking for meaningful error messages")
library(StanNet)

test_that("makeStanmodel displays meaningful error messages", {
  
  df<-data.frame(rep(0,10), nrow = 5, ncol = 2)
  
  dat1 <- data.simulator(nstudies=25, response="continuous", ntreatments=7, seed=13, cont_mean = seq(1,1.4, length.out = 7),
                 min.patients=40, max.patients=500) # continuous data for normal family with id link
  dat1 <- data.prep(dat1, varname.t = "treatment", varname.s = "study")
  
  dat2 <- data.simulator(nstudies=20, response="rate",
                        ntreatments=5, seed=4, min.patients=100, max.patients = 500,
                        pt_rate = 1, count_mean = seq(5,10, length.out = 5)/40)
  dat2 <- data.prep(dat2, varname.t = "treatment", varname.s = "study")
  
  dat3 <- data.simulator(nstudies=30, response="rate2",
                        ntreatments=8, seed=4, min.patients=100, max.patients = 500,
                        event_p = seq(0.2,0.3, length.out = 8), min_time = 1, max_time = 4)
  dat3 <- data.prep(dat3, varname.t = "treatment", varname.s = "study")
  
  
  expect_error(makeStanmodel(data = df), "'data' must be a valid BUGSnetData object created using the data.prep function.")
  expect_error(makeStanmodel(data = dat1, family = "normal"), "sd must be specified for continuous outcomes")
  expect_error(makeStanmodel(data = dat1, outcome = "response", N = "sampleSize", sd = "st.error",
                             reference = "A", family = "normal", link = "identity", effects = "fixed", type = "inconsistency"),
               "Inconsistency models have not been inplemented in StanNet yet.")
  expect_error(makeStanmodel(data = dat2, outcome = "events", N = "sampleSize",
                             reference = "A", family = "poisson", link = "log", effects = "fixed"),
               "time must be specified when using a poisson family with the log link")
  expect_error(makeStanmodel(data = dat3, outcome = "events", N = "sampleSize",
                             reference = "A", family = "binomial", link = "cloglog", effects = "random"),
               "time must be specified when using a binomial family with the cloglog link")
  expect_error(makeStanmodel(data = dat1, outcome = "response", N = "sampleSize", sd = "st.error",
                             reference = "A", family = "normal", link = "logit", effects = "fixed"),
               "This combination of family and link is currently not supported in StanNet.")
  expect_error(makeStanmodel(data = dat2, outcome = "events", N = "sampleSize", time = "pt.risk",
                             reference = "A", family = "poisson", link = "cloglog", effects = "fixed"), 
               "This combination of family and link is currently not supported in StanNet.")
  expect_error(makeStanmodel(data = dat3, outcome = "events", N = "sampleSize",
                             reference = "A", family = "binomial", link = "identity", effects = "random"), 
               "This combination of family and link is currently not supported in StanNet.")
  
})

# Generating data and results for testing functions
sim1 <- data.simulator(nstudies=30, response="binary", ntreatments=5, seed=20,
                                min.patients=100, max.patients=500, 
                                event_p=c(0.05, 0.10, 0.15, 0.07, 0.09))
sim1.prep <- data.prep(sim1, varname.t = "treatment", varname.s = "study")
sim1.nettab <- net.tab(sim1.prep, outcome = "events", 
                       N = "sampleSize", type.outcome = "binomial")
sim1.fe_model <- nma.model(data = sim1.prep, outcome = "events", N = "sampleSize", 
                                    reference = "D", family = "binomial", link = "logit", 
                                    effects = "fixed", type="consistency")
sim1.fe_results <- nma.run(model = sim1.fe_model, 
                                    n.adapt = 1000, n.burnin = 2500, n.iter = 2500,
                                    n.chains = 2)
rainbow.fe.sim1 <- nma.fit(sim1.fe_results)

test_that("Visualization function gives appropriate error messages", {
  
  expect_error(comparison.viz(sim1), "Call the net.tab function on your data first.")
  expect_error(comparison.viz(sim1.nettab, metric = "test-string"), "Please select a metric from the numeric column names of the data")
  expect_error(comparison.viz(sim1.nettab, metric = "treatment"), "Please select a metric from the numeric column names of the data")
  
  # Testing intervention.viz
  expect_error(intervention_viz(sim1), "Call the net.tab function on your data first.")
  expect_error(intervention_viz(sim1.nettab, metric = "test-string"), "Please select a metric from the numeric column names of the data")
  expect_error(intervention_viz(sim1.nettab, metric = "treatment"), "Please select a metric from the numeric column names of the data")
  
})

test_that("Checking invalid entries into data.simulator",{
  
  expect_error(data.simulator(response="binary", ntreatments=5, seed=20,
                              min.patients=100, max.patients=500, 
                              event_p=c(0.05, 0.10, 0.15, 0.07, 0.09)), "The number of studies, nstudies, must be specified")
  expect_error(data.simulator(nstudies=6, ntreatments=5, seed=20,
                              min.patients=100, max.patients=500, 
                              event_p=c(0.05, 0.10, 0.15, 0.07, 0.09)), "response must be specified")
  expect_error(data.simulator(nstudies=6, response="random", ntreatments=5, seed=20,
                              min.patients=100, max.patients=500, 
                              event_p=c(0.05, 0.10, 0.15, 0.07, 0.09)), "response selection must be binary, rate, rate2, or continuous")
  expect_error(data.simulator(nstudies=6, response="binary", seed=20,
                              min.patients=100, max.patients=500, 
                              event_p=c(0.05, 0.10, 0.15, 0.07, 0.09)), "The number of treatments, ntreatments, must be specified")
  expect_error(data.simulator(nstudies=6, response="binary", ntreatments=5, seed=20, max.patients=500, 
                              event_p=c(0.05, 0.10, 0.15, 0.07, 0.09)), "The minimum number of patients per trial, min.patients, must be specified")
  expect_error(data.simulator(nstudies=6, response="binary", ntreatments=5, seed=20, min.patients=500, 
                              event_p=c(0.05, 0.10, 0.15, 0.07, 0.09)), "The maximum number of patients per trial, max.patients, must be specified")
  expect_error(data.simulator(nstudies=6, response="binary", ntreatments=5, seed=20, min.patients=500, max.patients = 6,
                              event_p=c(0.05, 0.10, 0.15, 0.07, 0.09)), "The minimum number of patients cannot exceed the maximum")
  expect_error(data.simulator(nstudies=2, response="binary", ntreatments=10, seed=20, min.patients=500, max.patients = 600,
                              event_p=c(0.05, 0.10, 0.15, 0.07, 0.09)), "The number of unique treatments cannot exceed the total number of study arms")
  expect_error(data.simulator(nstudies=10, response="binary", ntreatments=2, seed=20, min.patients=500, max.patients = 600,
                              event_p=c(0.05, 0.10, 0.15, 0.07, 400)), "The event probabilities must be valid probabilities")
  expect_error(data.simulator(nstudies=10, response="binary", ntreatments=2, seed=20, min.patients=500, max.patients = 600,
                              event_p=c(0.05, 0.10, 0.15, 0.07)), "Event probability vector must be the same length as the number of treatments")
  expect_error(data.simulator(nstudies=10, response="rate2", ntreatments=2, seed=20, min.patients=500, max.patients = 600,
                              event_p=c(0.05, 0.10, 0.15, 0.07)), "Event probability vector must be the same length as the number of treatments")
  expect_error(data.simulator(nstudies=10, response="rate2", ntreatments=2, seed=20, min.patients=500, max.patients = 600,
                              event_p=c(0.05, 0.10, 0.15, 0.07, 0.5), exp_rate = c(1,2)), "Exponential rate vector must be a scalar")
  expect_error(data.simulator(nstudies=12, response="continuous", ntreatments=6, seed=13, min.patients=30, max.patients=200, cont_mean = c(1,2)), 
               "Continuous response mean vector must be the same length as the number of treatments")
  expect_error(data.simulator(nstudies=12, response="continuous", ntreatments=6, seed=13, min.patients=30, max.patients=200, cont_var= c(1,2)), 
               "Continuous response variance vector must be the same length as the number of treatments")
  expect_error(data.simulator(nstudies=12, response="continuous", ntreatments=6, seed=13, min.patients=30, max.patients=200, cont_error= c(1,2)), 
               "Continuous response standard error vector must be the same length as the number of treatments")
  expect_error(data.simulator(nstudies=40, response="rate", ntreatments=5, seed=30,
                              min.patients=1000, max.patients=1200, 
                              pt_rate = 1000*0.5,
                              count_mean = c(45)), "Count mean vector must be the same length as the number of treatments")
  expect_error(data.simulator(nstudies=40, response="rate", ntreatments=5, seed=30,
                              min.patients=1000, max.patients=1200, 
                              pt_rate = c(1000*0.5, 2),
                              count_mean = c(45, 35, 47, 40, 42)), "Person-time-at-risk rate vector must be of length 1")
  
})

test_that("Check output from outlier.detector", {
  
  expect_match(outlier_detector(rainbow.fe.sim1, c=100), "No outliers detected for given")
  expect_error(outlier_detector(sim1.fe_results, c=0), "Call 'Stan.nma.fit' function on your model object first")
  
})
