source("_packages.R")
source("_helpers/_helpers.R")

T <- 50
n <- 2
n.ahead = 5
nsim = 100

set.seed(1)

Y <- matrix(rpois(T*n, lambda = 3), T, n)

pois_mod <- SSModel(Y ~ -1 + SSMtrend(degree = 1, Q = matrix(NA), type = "distinct"), 
                    distribution = "poisson")
object <- fitSSM(pois_mod, inits = c(diag(n)))$model

fcast_draws = importMvSim(object, nsim / 4, n.ahead)
