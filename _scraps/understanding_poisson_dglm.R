source("_packages.R")
source("_helpers/_helpers.R")

T <- 200
n <- 2

Y <- matrix(rpois(T*n, lambda = 3), T, n)

y <- rpois(T, lambda = 10)

update_model <- function(pars, model) {
  model["Q"] <- pars[1]
  model
}
#check that variances are non-negative
check_model <- function(model) {
  (model["Q"] > 0)
}

pois_mod <- SSModel(y[1:25] ~ -1 + SSMtrend(degree = 1, Q = matrix(NA)), distribution = "poisson")
pois_fit <- fitSSM(pois_mod, inits = c(1), method = "BFGS",
                   updatefn = update_model, checkfn = check_model)
pois_fc <- importSim(pois_fit$model, nsim = 1250, n.ahead = 10)

#Negative Binomial 
nb_mod <- SSModel(y[1:i] ~ -1 + SSMtrend(degree = 1, Q = matrix(NA)), 
                  distribution = "negative binomial")
nb_fit <- fitSSM(nb_mod, inits = c(0.001), method = "BFGS",
                 updatefn = update_model, checkfn = check_model)
nb_fc <- importSim(nb_fit$model, nsim= 1250, n.ahead=1)
dglm_results <- list(as.vector(pois_fc), as.vector(nb_fc))
