#' @description Creates a string which represents the full Stan code of the model
#' @title Write Stan code for a model
#' @description Create a string containing the Stan code which describes the specified model
#' @param family string indicating family of the distribution of the outcome. Options are: "binomial", "normal", "poisson", ADD?
#' @param link The link function for the nma model. Options are "logit" (binomial family), "cloglog" (binomial family), "log" (poisson family), and "identity" (normal family).
#' @param effects A string indicating the type of treatment effect relative to the baseline, either "fixed" or "random"
#' @param inconsistency Logical indicating if the model is an inconsistency model, if TRUE then an inconsistency model is built, default is FALSE
#' @param prior.mu.str A string of Stan code which defines priors on the baseline treatment effects, or if "default", independent normal priors are used with mean 0 and standard deviation 15u, where u is the largest maximum likelihood estimator in single trials (see van Valkenhoef, G et al. 2012)
#' @param prior.d.str A string of Stan code that defines priors on the relative treatment effects. By default, independent normal priors with mean 0 and standard deciation 15u where u is the largest maximum likelihood estimator in single trials (see van Valkenhoef, G et al. 2012) (BUGSnet)
#' @param prior.sigma.str A string of Stan code that defines the prior on the standard deviation of the relative treatment effects. By default, uniform distribution on (0, u) where u is largest MLE from single trials(see van Valkenhoef, G et al. 2012) (BUGSnet)
#' 
#' @return A string of Stan code representing the full model
#' @export
makeStancode <- function (family, link, effects, inconsistency, prior.mu.str, prior.d.str, prior.sigma.str) {
  
  if (inconsistency) {
    
    stop("Inconsistency models not supported in Stan")
    
  }
  
  # Set up strings for data block, transformed parameters, generated parameters based on family and link
  
  if (family == "normal" && link == "identity") {
    
    outcome.str <- "real y[ns,2]; //outcomes in ith study, jth arm"
    other.str <- "real<lower=0> se[ns,2]; //standard errors"
    
    RHS.str <- "real theta[ns, 2]; // RHS of linear predictor in GLM"
    
    mod.str <- "y[i,j] ~ normal(theta[i,j], se[i,j]); // Model the response"
      
    dev1.str <- ""
    dev2.str<- "dev[i,j] = (y[i,j]-theta[i,j])*(y[i,j]-theta[i,j])*pow(se[i,j],-2); //Deviance contribution"
    
  } else if (family == "binomial" && link == "logit") {
    
    outcome.str <- "int<lower = 0> r[ns,2]; //outcomes in ith study, jth arm"
    other.str <- "int<lower = 0> n[ns,2]; //total subjects in ith study, jth arm"
    
    RHS.str <- "real logitp[ns, 2]; // RHS of linear predictor in GLM"
    
    mod.str <- "r[i,j] ~ binomial_logit(n[i,j], logitp[i,j]); // Model the response"
    
    dev1.str <- "real rhat[ns,2];"
    dev2.str <- "rhat[i,j] = inv_logit(logitp[i,j]) * n[i,j];
      dev[i,j] = 2 * (r[i,j] * (log(r[i,j])-log(rhat[i,j]))
                       + (n[i,j]-r[i,j]) * (log(n[i,j]-r[i,j]) - log(n[i,j]-rhat[i,j]))); //Deviance contribution"
    
  } else if (family == "binomial" && link == "cloglog") {
    
    outcome.str <- "int<lower = 0> r[ns,2]; //outcomes in ith study, jth arm"
    other.str <- "int<lower = 0> n[ns,2]; //total subjects in ith study, jth arm
    real time[ns,2]; // study followup times"
    
    RHS.str <- "real cloglogp[ns,2]; // RHS of linear predictor in GLM"
    
    mod.str <- "r[i,j] ~ binomial(n[i,j], inv_cloglog(cloglogp[i,j])); // Modelling the response"
    
    dev1.str <- "real rhat[ns, 2];"
    dev2.str <- "rhat[i,j] = inv_cloglog(cloglogp[i,j]) * n[i,j];
      dev[i,j] = 2 * (r[i,j] * (log(r[i,j])-log(rhat[i,j]))
                       + (n[i,j]-r[i,j]) * (log(n[i,j]-r[i,j]) - log(n[i,j]-rhat[i,j]))); //Deviance contribution"
    
  } else if (family == "poisson" && link == "log") {
    
    outcome.str <- "int<lower = 0> r[ns,2]; //outcomes in ith study, jth arm"
    other.str <- "real E[ns,2]; //total exposure for ith study, jth arm"
    
    RHS.str <- "real loglambda[ns, 2]; // RHS of linear predictor in GLM"
    
    mod.str <- "r[i,j] ~ poisson_log(loglambda[i,j]+log(E[i,j])); // Modelling the response - log(lambda * E) = loglambda + log E"
    
    dev1.str <- "real theta[ns, 2];"
    dev2.str <- "theta[i, j] = exp(loglambda[i,j]) * E[i,j];
      dev[i,j] = 2 * ((theta[i,j]-r[i,j]) + r[i,j]*log(r[i,j]/theta[i,j])); //Deviance contribution"
    
  } else { 
    
    stop("Family and link combination not supported")
    
  }
  
  # Make model code
  
  if (effects == "random") {
    
    if (family == "normal" && link == "identity") {
      
      GLM.str <- "theta[i,j] = mu[i] + deltas[i,j]; // model for linear predictor"
      
    } else if (family == "binomial" && link == "logit") {
      
      GLM.str <- "logitp[i,j] = mu[i] + deltas[i,j]; // model for linear predictor"
      
    } else if (family == "binomial" && link == "cloglog") {
      
      GLM.str <- "cloglogp[i,j] = log(time[i,j]) + mu[i] + deltas[i,j]; // model for linear predictor"
      
    } else if (family == "poisson" && link == "log") {
      
      GLM.str <- "loglambda[i,j] = mu[i] + deltas[i,j]; // model for linear predictor"
      
    } else { 
      
      stop("Family and link combination not supported")
      
    }
    
    # Put together data block
    
    data.str <- sprintf(
      
      "data {
      
        int ns;
        int nt;
        int t[ns,2];
        %s
        %s
      
      }", outcome.str, other.str)
    
    # Put together parameters block
    
    parameters.str <- sprintf(
      
      "parameters {
      
        real d2[nt-1];
        vector [ns] mu;
        real delta2raw[ns];
        real <lower = 0> sigma;
      
      }")
    
    # Put together transformed parameters block
    
    trans.parameters.str <- sprintf(
      
      "transformed parameters {
      
        %s
        real d[nt];
        real deltas[ns,2];
        
        d[1] = 0;
        d[2:nt] = d2;
        
        deltas[,1] = rep_array(0.0, ns);
        
        for (i in 1:ns) {
        
          for (j in 1:2) {
          
          if (j>1) {
          
            deltas[i,j] = delta2raw[i]*sigma+d[t[i,j]]-d[t[i,1]];
          
          }
          
            %s
          
          }
        
        }
      
      }", RHS.str, GLM.str)
    
    # Put together model block
    
    model.str <- sprintf(
      
      "model {
      
        %s;
        %s;
        %s;
        
        delta2raw ~ std_normal();
        
        for (i in 1:ns) {
        
          for (j in 1:2) {
          
            %s
          
          }
        
        }
      
      }", prior.sigma.str, prior.d.str, prior.mu.str, mod.str)
    
    # Put together generated quantities block
    
    gen.q.str <- sprintf(
      
      "generated quantities {
      
         real dev[ns, 2];
         real resdev[2];
         real totresdev;
         %s
         
         for (j in 1:2) {
         
           for (i in 1:ns) {
           
             %s
           
           }
           
           resdev[j] = sum(dev[,j]);
         
         }
         
         totresdev = sum(resdev[]);
      
      }", dev1.str, dev2.str)
    
    
  } else if (effects == "fixed") {
    
    if (family == "normal" && link == "identity") {
      
      GLM.str <- "theta[i,j] = mu[i] + d[t[i,j]] - d[t[i,1]]; // model for linear predictor"
      
    } else if (family == "binomial" && link == "logit") {
      
      GLM.str <- "logitp[i,j] = mu[i] + d[t[i,j]] - d[t[i,1]]; // model for linear predictor"
      
    } else if (family == "binomial" && link == "cloglog") {
      
      GLM.str <- "cloglogp[i,j] = log(time[i,j]) + mu[i] + d[t[i,j]] - d[t[i,1]]; // model for linear predictor"
      
    } else if (family == "poisson" && link == "log") {
      
      GLM.str <- "loglambda[i,j] = mu[i] + d[t[i,j]] - d[t[i,1]]; // model for linear predictor"
      
    } else { 
     
      stop("Family and link combination not supported")
       
    }
    
    # Put together data block
    
    data.str <- sprintf(
      
      "data {
      
        int ns;
        int nt;
        int t[ns,2];
        %s
        %s
      
      }", outcome.str, other.str)
    
    # Put together parameters block
    
    parameters.str <- sprintf(
      
      "parameters {
      
        real d2[nt-1];
        vector [ns] mu;
      
      }")
    
    # Put together transformed parameters block
    
    trans.parameters.str <- sprintf(
      
      "transformed parameters {
      
        %s
        real d[nt];
        
        d[1] = 0;
        d[2:nt] = d2;
        
        for (i in 1:ns) {
        
          for (j in 1:2) {
          
            %s
          
          }
        
        }
      
      }", RHS.str, GLM.str)
    
    # Put together model block
    
    model.str <- sprintf(
      
      "model {
      
        %s;
        %s;
        
        for (i in 1:ns) {
        
          for (j in 1:2) {
          
            %s
          
          }
        
        }
      
      }", prior.d.str, prior.mu.str, mod.str)
    
    # Put together generated quantities block
    
    gen.q.str <- sprintf(
      
      "generated quantities {
      
         real dev[ns, 2];
         real resdev[2];
         real totresdev;
         %s
         
         for (j in 1:2) {
         
           for (i in 1:ns) {
           
             %s
           
           }
           
           resdev[j] = sum(dev[,j]);
         
         }
         
         totresdev = sum(resdev[]);
      
      }", dev1.str, dev2.str)
    
  }
  
  # Put together blocks to form the full model
  
  code.str <- sprintf("%s\n\n%s\n\n%s\n\n%s\n\n%s", data.str, parameters.str, trans.parameters.str, model.str, gen.q.str)
  
  return (code.str)
  
}
