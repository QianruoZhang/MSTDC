# ====================================================================================
# ==== Poisson simulation study ====
# ====================================================================================

## This simulation study is designed to assess the inferential performance of the three-stage approximate model fitting algorithm.

## Competitor models include:
# Sep - The separable model proposed by Knorr-Held and Besag (1998)
# KN4 - The non-separable model with spatio-temporally autocorrelated interactions proposed by Knorr-Held (2000)
# Rush - The spatially correlated temporal AR(1) process model proposed by Rushworth et al. (2014)
# ORDER-1 - The 1-order neighbourhood bigDM algorithm proposed by Orozco-Acosta et al. (2023)
# MSTDC-MB - The proposed three-stage model fitting algorithm with a model-based smoother
# MSTDC-WA - The proposed three-stage model fitting algorithm with simple average weights

## The simulated dataset is generated under Scenario ST-2, in which the spatial structure of the disease prevalence changes slowly over time. 
## It represents a disease prevalence surface with a common disease pattern and a Leroux CAR spatial autocorrelation structure.



# ====================================================================================
# ==== Read in libraries ====
# ====================================================================================

library(tidyverse)
library(spdep)
library(MASS)
library(sf)
library(INLA)
library(igraph) # Define graph distance
library(bigDM)



# ====================================================================================
# ==== Read in files ====
# ====================================================================================

## Load the simulated neighbourhood matrix containing 3,487 LSOAs
load("Simulated_W.Rdata")

## Load the simulated dataset covering 3,487 LSOAs over 10 time periods
load("Simulated_Poi data.Rdata")

## Load the simulated LSOA object of class sf, which contains LSOAs specified with potential subgroup augmentation
load("Simulated_LSOAsf.Rdata")



# ====================================================================================
# ==== Generate quantities ====
# ====================================================================================

n.S <- length(unique(dat$LSOA)) # Number of spatial units
n.T <- length(unique(dat$year)) # Number of time periods
n.all <- n.S * n.T # Number of spatial-temporal units


theta0.est <- 0.1 # Modify initial estimates of unobserved risk

# Define graph distance matrix for MSTDC-WA
graphs <- graph_from_adjacency_matrix(W, mode = "max", diag = FALSE) 
D.graph <- distances(graphs) 

epsilon <- 0 # Weight's threshold for MSTDC-WA



# ====================================================================================
# ==== Generate the remaining quantities needed by INLA ====
# ====================================================================================

## Define identifiers

# Temporal identifiers
dat$time1 <- rep(1:n.T, each = n.S)
dat$time2 <- rep(1:n.T, each = n.S)
dat$time3 <- rep(1:n.T, each = n.S)

# Spatial identifier
dat$space1 <- rep(1:n.S, time = n.T)
dat$space2 <- rep(1:n.S, time = n.T)

# Spatio-temporal identifiers 
dat$space.time <- rep(1:n.all, time = 1)


## Define the INLA graph object
W.g <- inla.write.graph(W, "LSOAgraph-WM.png")


## Specify the half normal prior function for the standard deviation
HN.prior = "expression:
  tau0 = 0.001;
  sigma = exp(-theta/2);
  log_dens = log(2) - 0.5 * log(2 * pi) + 0.5 * log(tau0);
  log_dens = log_dens - 0.5 * tau0 * sigma^2;
  log_dens = log_dens - log(2) - theta / 2;
  return(log_dens);  
"



# ====================================================================================
# ==== Simulation study  ====
# ====================================================================================

### ==== Separable Knorr-Held model ====

time.sep <- system.time({
  # Generate formula
  form.sep <- Y ~ 
    offset(log(popsize)) +
    f(time1, model = 'rw1', constr = TRUE, hyper = list(theta = list(prior = HN.prior))) + 
    f(time2, model = 'iid', constr = TRUE, hyper = list(theta = list(prior = HN.prior))) + 
    f(space1, model = 'bym', graph = W.g, constr = TRUE, hyper = list(theta1 = list(prior = HN.prior), theta2 = list(prior = HN.prior)))
  
  # Model fitting
  mod.sep <- inla(
    formula = form.sep,
    family = 'poisson', 
    data = dat,
    control.inla = list(strategy = 'gaussian'), 
    control.fixed = list(mean = 0, mean.intercept = 0, prec = 0.00001, prec.intercept = 0.00001),
    control.compute = list(dic = TRUE, mlik = TRUE, waic = TRUE),
    control.predictor = list(link = 1, compute = TRUE))
  
  pred.sep <- mod.sep$summary.fitted.values$mean / dat$popsize
  ci.sep <- mod.sep$summary.fitted.values[ , c(3, 5)] / dat$popsize
})

# Remove mod.sep from the environment
rm(mod.sep)
  

  
### ==== Knorr-Held Type IV model ====

time.kn4 <- system.time({
  # Generate formula
  form.kn4 <- Y ~ 
    offset(log(popsize)) + 
    f(time1, model = 'rw1', constr = TRUE, hyper = list(theta = list(prior = HN.prior))) + 
    f(time2, model = 'iid', constr = TRUE, hyper = list(theta = list(prior = HN.prior))) +
    f(space1, model = 'bym', graph = W.g, constr = TRUE, hyper = list(theta1 = list(prior = HN.prior), theta2 = list(prior = HN.prior))) + 
    f(space2, model = 'besag', graph = W.g, constr = TRUE, group = time3, control.group = list(model = 'rw1'), 
      hyper = list(theta = list(prior = HN.prior)))
  
  # Model fitting
  mod.kn4 <- inla(
    formula = form.kn4, 
    family = 'poisson', 
    data = dat,
    control.inla = list(strategy = 'gaussian'), 
    control.fixed = list(mean = 0, mean.intercept = 0, prec = 0.00001, prec.intercept = 0.00001),
    control.compute = list(dic = TRUE, mlik=TRUE, waic=TRUE),
    control.predictor = list(link = 1, compute = TRUE))
  
  pred.kn4 <- mod.kn4$summary.fitted.values$mean / dat$popsize
  ci.kn4 <- mod.kn4$summary.fitted.values[ , c(3, 5)] / dat$popsize
})

# Remove mod.kn4 from the environment
rm(mod.kn4)  

  

### ==== Rushworth model ====

time.rush <- system.time({
  # Generate formula
  form.rush <- Y ~ 
    offset(log(popsize)) + 
    f(space1, model = 'besagproper2', graph = W.g, constr = FALSE, group = time1, 
      control.group = list(model = 'ar1', hyper = list(rho = list(param = c(0, 0.25)))), 
      hyper = list(theta1 = list(prior = HN.prior), theta2 = list(param = c(0, 0.25))))
  
  # Model fitting
  mod.rush <- inla(
    formula = form.rush, 
    family = 'poisson', 
    data = dat,
    control.inla = list(strategy = 'gaussian'), 
    control.fixed = list(mean = 0, mean.intercept = 0, prec = 0.00001, prec.intercept = 0.00001),
    control.compute = list(dic = TRUE, mlik = TRUE, waic = TRUE),
    control.predictor = list(link = 1, compute = TRUE))
  
  pred.rush <- mod.rush$summary.fitted.values$mean / dat$popsize
  ci.rush <- mod.rush$summary.fitted.values[ , c(3, 5)] / dat$popsize
})

# Remove mod.rush from the environment
rm(mod.rush)



### ==== k-order neighbourhood bigDM model ====

time.order1 <- system.time({
  mod.order1 <- STCAR_INLA(
    carto = lsoa.sp, 
    data = dat,
    ID.area = "LSOA", 
    ID.year = "year", 
    O = "Y", 
    E = "popsize", 
    ID.group = "LAD11CD",
    spatial = "BYM2", 
    temporal = "rw1", 
    interaction = "TypeIV",
    model = "partition", 
    k = 1, 
    strategy = "gaussian",
    merge.strategy = 'original',
    plan = "cluster",
    workers = rep("localhost", 4))
  
  pred.order1 <- exp(mod.order1$summary.linear.predictor$mean) 
  ci.order1 <- exp(mod.order1$summary.linear.predictor[ , c(3,5)])
})

# Remove mod.order1 from the environment
rm(mod.order1)



### ==== MSTDC algorithm with a model-based smoother ====

time.mb <- system.time({
  
  ### Stage 1: Time-averaged spatial structure
  
  # The initial estimates of true unobserved prevalence
  theta0.mat <- matrix(dat$Y / dat$popsize, nrow = n.S, ncol = n.T)
  theta0.mat[theta0.mat == 0] <- theta0.est
  lr0.mat <- log(theta0.mat)
  
  # Time-average spatial pattern
  ave.phi <- apply(lr0.mat, 1, mean, na.rm = TRUE) 
  
  
  ### Stage 2: Time-period specific residual spatial structure#
  
  res.phi.mb <- matrix(NA, nrow = n.S, ncol = n.T)
  for(t in 1:n.T){
    # Subset data for time period t
    dat.mb <- dat %>%
      filter(year == t) %>%
      mutate(ave.phi = ave.phi)
    
    # Generate formula for a Poisson logistic model
    dat.mb$offsetnew <- log(dat.mb$popsize) + dat.mb$ave.phi
    form.mb <- Y ~ 
      offset(offsetnew) +
      f(space1, model = 'besagproper2', graph = W.g, constr = FALSE, hyper = list(theta1 = list(prior = HN.prior), theta2 = list(param = c(0, 0.25))))
    
    # Fitting a Poisson log model with offsets
    mod.mb <- inla(
      formula = form.mb, 
      family = "poisson", 
      data = dat.mb,
      control.inla = list(strategy = 'gaussian'), 
      control.fixed = list(mean = 0, mean.intercept = 0, prec = 0.00001, prec.intercept = 0.00001),
      control.compute = list(dic = TRUE, mlik = TRUE, waic = TRUE),
      control.predictor = list(link = 1, compute = TRUE))
    
    # Estimated residual spatial random effects
    res.phi.mb[ , t] <- mod.mb$summary.random$space1$mean
    
    rm(mod.mb)
    rm(dat.mb)
  }
  
  
  ### Stage 3: Final spatio-temporal model
  
  res.phi.mb <- as.numeric(res.phi.mb)
  dat$offsetnew <- log(dat$popsize) + as.numeric(res.phi.mb)
  # Generate formula for a Poisson model 
  form.final.mb <- Y ~ 
    offset(offsetnew) +
    f(time1, model = 'ar1', constr = FALSE, hyper = list(theta1 = list(prior = HN.prior), theta2 = list(param = c(0, 0.25)))) + 
    f(space1, model = 'besagproper2', graph = W.g, constr = FALSE, hyper = list(theta1 = list(prior = HN.prior), theta2 = list(param = c(0, 0.25))))
  
  # Model fitting
  mod.final.mb <- inla(
    formula = form.final.mb,
    family = 'poisson',
    data = dat,
    control.inla = list(strategy = 'gaussian'), 
    control.fixed = list(mean = 0, mean.intercept = 0, prec = 0.00001, prec.intercept = 0.00001),
    control.compute = list(dic = TRUE, mlik = TRUE, waic = TRUE),
    control.predictor = list(link = 1, compute = TRUE))
  
  # Predict estimates
  pred.final.mb <- mod.final.mb$summary.fitted.values$mean / dat$popsize
  ci.final.mb <- mod.final.mb$summary.fitted.values[ , c(3, 5)] / dat$popsize
})

# Remove mod.final.mb from the environment
rm(mod.final.mb)



### ==== MSTDC algorithm with average weights ====

time.wa <- system.time({
  
  ### Stage 1: Time-averaged spatial structure
  
  # The initial estimates of true unobserved prevalence
  theta0.mat <- matrix(dat$Y / dat$popsize, nrow = n.S, ncol = n.T)
  theta0.mat[theta0.mat == 0] <- theta0.est
  lr0.mat <- log(theta0.mat)
  
  # Time-average spatial pattern
  ave.phi <- apply(lr0.mat, 1, mean, na.rm = TRUE) 
  
  
  ### Stage 2: Time-period specific residual spatial structure
  
  # Initial estimates of remaining spatial structure
  res.phi0 <- lr0.mat - ave.phi
  
  # Define spatially decaying weights matrix
  W.star <- exp(-D.graph)
  W.star[W.star <= epsilon] <- 0
  
  num <- W.star %*% replace(res.phi0, is.na(res.phi0), 0)
  den <- rowSums(W.star)
  res.phi.wa <-  num/den
  
  
  ### Stage 3: Final spatio-temporal model
  
  res.phi.wa <- as.numeric(res.phi.wa)
  dat$offsetnew <- log(dat$popsize) + as.numeric(res.phi.wa)
  # Generate formula
  form.final.wa <- Y ~ 
    offset(offsetnew) +
    f(time1, model = 'ar1', constr = FALSE, hyper = list(theta1 = list(prior = HN.prior), theta2 = list(param = c(0, 0.25)))) + 
    f(space1, model = 'besagproper2', graph = W.g, constr = FALSE, hyper = list(theta1 = list(prior = HN.prior), theta2 = list(param = c(0, 0.25))))
  
  # Model fitting
  mod.final.wa <- inla(
    formula = form.final.wa,
    family = 'poisson',
    data = dat,
    control.inla = list(strategy = 'gaussian'), 
    control.fixed = list(mean = 0, mean.intercept = 0, prec = 0.00001, prec.intercept = 0.00001),
    control.compute = list(dic = TRUE, mlik = TRUE, waic = TRUE),
    control.predictor = list(link = 1, compute = TRUE))
  
  # Predict estimates
  pred.final.wa <- mod.final.wa$summary.fitted.values$mean / dat$popsize
  ci.final.wa <- mod.final.wa$summary.fitted.values[ , c(3, 5)] / dat$popsize
})

# Remove mod.final.wa from the environment
rm(mod.final.wa)