# ==============================================================================
# My multivariate stuff
# ==============================================================================

importMvSim <- function(object, nsim, n.ahead){
  nvar <- ncol(object$y)
  #Update model with appropriate NAs
  timespan <- attr(object, "n") + 1:n.ahead
  n <- attr(object, "n") <- attr(object, "n") + as.integer(n.ahead)
  endtime<-end(object$y) + c(0, n.ahead)
  object$y <- window(object$y, end = endtime, extend = TRUE)
  object$u <- rbind(object$u, matrix(object$u[1, ], nrow = n.ahead,
                                     ncol = ncol(object$u), byrow = TRUE))
  
  #Draw samples from signal (alpha_t)
  imp <- KFAS::importanceSSM(object, "signal", nsim = nsim, antithetics = TRUE)
  nsim <- as.integer(4 * nsim)
  w <- imp$weights/sum(imp$weights)
  imp$samples <- exp(imp$samples[timespan, , ])
  n <- as.integer(length(timespan))
  fcast_draws = array(0, c(n.ahead, nvar, nsim))
  for(m in 1:nsim){
    alpha = imp$samples[, , sample(1:nsim, size = 1, prob = w)]
    if(object$distribution[1] == "poisson"){
      result = rpois(n.ahead * nvar, lambda = c(alpha))
    }else{
      result = rnbinom(n.ahead * nvar, size = c(object$u[timespan, ]), mu = c(alpha))
    }
    fcast_draws[, , m] = matrix(result, n.ahead, nvar)
  }
  return(fcast_draws)
}

get_poisson_mvforecasts <- function(Y, H, nsim){
  pois_mod <- SSModel(Y ~ -1 + SSMtrend(degree = 1, Q = matrix(NA), type = "distinct"), 
                      distribution = "poisson")
  object <- fitSSM(pois_mod, inits = c(diag(ncol(Y))))$model
  fcast_draws = importMvSim(object, nsim / 4, H)
  return(fcast_draws)
}

get_nbinom_mvforecasts <- function(Y, H, nsim){
  pois_mod <- SSModel(Y ~ -1 + SSMtrend(degree = 1, Q = matrix(NA), type = "distinct"), 
                      distribution = "negative binomial")
  object <- fitSSM(pois_mod, inits = c(diag(ncol(Y))))$model
  fcast_draws = importMvSim(object, nsim / 4, H)
  return(fcast_draws)
}

# ==============================================================================
# Brian's univariate stuff
# ==============================================================================

importSim <- function(object, nsim, n.ahead){
  #This means there is only one state
  j=1
  #Update model with appropriate NAs
  timespan <- attr(object, "n") + 1:n.ahead
  n <- attr(object, "n") <- attr(object, "n") + as.integer(n.ahead)
  endtime<-end(object$y) + c(0, n.ahead)
  object$y <- window(object$y, end = endtime, extend = TRUE)
  object$u <- rbind(object$u, matrix(object$u[1, ], nrow = n.ahead,
                                     ncol = ncol(object$u), byrow = TRUE))
  
  #Draw samples from signal (alpha_t)
  imp <- KFAS::importanceSSM(object, "signal", nsim = nsim, antithetics = TRUE)
  nsim <- as.integer(4 * nsim)
  w <- imp$weights/sum(imp$weights)
  imp$samples <- imp$samples[timespan, , , drop = F]
  imp$samples[,j, ] <- exp(imp$samples[, j, ])
  n <- as.integer(length(timespan))
  preds <- sapply(1:n, function(i) {
    sample_mu <- sample(imp$samples[i, j, ], size = nsim, replace = TRUE,
                        prob = w)
    if(object$distribution == "poisson")
      result= rpois(n = nsim, lambda = sample_mu)
    else
      result= rnbinom(n = nsim, size = object$u[timespan[i],j], mu = sample_mu)
    return(result)
  })
  return(preds)
}

update_model <- function(pars, model) {
  model["Q"] <- pars[1]
  model
}

#check that variances are non-negative
check_model <- function(model) {
  (model["Q"] > 0)
}

get_poisson_forecasts <- function(Y, H, nsim){
  n <- ncol(Y)
  fcast_draws = array(0, dim = c(H, n, nsim))
  for(i in 1:n){
    pois_mod <- SSModel(Y[, i] ~ -1 + SSMtrend(degree = 1, Q = matrix(NA)), 
                        distribution = "poisson")
    pois_fit <- fitSSM(pois_mod, inits = 1, method = "Brent",
                       updatefn = update_model, checkfn = check_model,
                       lower = 0, upper = 1000)
    fcast_draws[, i, ] <- t(importSim(pois_fit$model, nsim = nsim / 4, n.ahead = H))
  }
  return(fcast_draws)
}

get_nbinom_forecasts <- function(Y, H, nsim){
  n <- ncol(Y)
  fcast_draws = array(0, dim = c(H, n, nsim))
  for(i in 1:n){
    nb_mod <- SSModel(Y[, i] ~ -1 + SSMtrend(degree = 1, Q = matrix(NA)), 
                      distribution = "negative binomial")
    nb_fit <- fitSSM(nb_mod, inits = c(0.001), method = "Brent",
                     updatefn = update_model, checkfn = check_model,
                     lower = 0, upper = 1000)
    fcast_draws[, i, ] <- t(importSim(nb_fit$model, nsim = nsim / 4, n.ahead = H))
  }
  return(fcast_draws)
}

