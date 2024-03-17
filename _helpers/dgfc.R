simulate_dgfc <- function(T, Finv, G, S, L, M){
  # Simulate the dynamic Gaussian factor copula:
  #
  # e[t] = G1 * e[t-1] + ... + Gp * e[t-p] + u[t],  u[t] ~ N(0, S)
  # z[t] = L * e[t] + v[t],                         v[t] ~ N(0, M)
  # y[t, i] = Finv(Phi(z[t, i] / sd(z[t, i])))
  #
  #    T: integer
  # Finv: n-list of function handles
  #    G: k x k x p array
  #    S: k x k matrix
  #    L: n x k matrix
  #    M: n x n matrix
  
  n = length(Finv)
  k = nrow(S)
  VARparams = var_stationary_params(numeric(k), G, S)
  P = VARparams$VARcov
  E = simulate_stationary_var(T, numeric(k), G, S)
  Z = E %*% t(L) + mvrnorm(T, numeric(n), M)
  W = L %*% P %*% t(L) + M
  Y = matrix(0, T, n)
  for(i in 1:n){
    Y[, i] = Finv[[i]](pnorm(Z[, i], mean = 0, sd = sqrt(W[i, i])))
  }
  return(Y)
}

DGFC.mcmc <- function(Y, prior = NULL, init = NULL, ndraw = 1000, burn = 0, thin = 1){
  # ----------------------------------------------------------------------------
  # dimensions
  # ----------------------------------------------------------------------------
  
  # JZ: need to decide on this n vs. p business
  
  T = nrow(Y)
  n = ncol(Y)
  p = n
  M = burn + thin * ndraw
  
  # ----------------------------------------------------------------------------
  # create or unpack prior hyperparameters
  # ----------------------------------------------------------------------------
  
  if(is.null(prior)){
    
    k.star = ceiling(.7 * n)
    
    # MGP
    # local shrinkage
    nu = 3
    # global shrinkage for MGP
    a1 = 2
    a2 = 3
    
    # error Variance
    a.sigma = 1
    b.sigma = 0.3
    
    # VAR params
    var_lag = 1
    inclconst = FALSE
    x_length = k.star * var_lag + inclconst
    v0 = k.star + 1
    P0 = diag(k.star)
    B0 = matrix(0, x_length, k.star)
    invO0 = diag(x_length)
    GS_prior = list(v = v0, P = P0, B = B0, invO = invO0)
    
  } else {
    
    # unpack whatever was passed in
    
  }
  
  # ----------------------------------------------------------------------------
  # create or unpack MCMC initialization
  # ----------------------------------------------------------------------------
  
  if(is.null(init)){
    
    Z<-NULL
    for(j in 1:n) { Z<-cbind(Z, scale(Y[,j])) }
    
    # Use SVD for simplicity
    svd0 = svd(Z);
    
    # Factor loadings (n x k.star):
    Lambda = Lambda.s = as.matrix(svd0$v[,1:k.star])
    
    ##MGP initialization
    #Local and global precision parameters:
    phi.jh = 1 / Lambda^2
    tau.h = delta.h = rep(1, k.star)
    
    # Factors (T x k.star):
    eta = svd0$u[, 1:k.star] %*% diag(svd0$d[1:k.star], k.star) #y%*%Lambda
    
    # Residuals (T x n):
    eps = Z - (tcrossprod(eta, Lambda))
    
    # Diagonal error variance (n-dimensional vector):
    Sigma.diag = apply(eps, 2, var)
    
  }else{
    
    # unpack whatever was passed in
    
  }
  
  # ----------------------------------------------------------------------------
  # JZ: Why do we need this?
  # ----------------------------------------------------------------------------
  
  # location of each y-val in the sorted list of unique values
  R <- NULL
  for(j in 1:p) { R<-cbind(R, match(Y[,j], sort(unique(Y[,j])))) }
  
  # number of unique values per variable
  Rlevels <- apply(R, 2, max, na.rm=TRUE)
  
  # ----------------------------------------------------------------------------
  # preallocate storage
  # ----------------------------------------------------------------------------
  
  Lambda_draws = array(0, c(n, k.star, ndraw))
  G_draws = array(0, c(k.star, x_length, ndraw))
  S_draws = array(0, c(k.star, k.star, ndraw))
  Z_draws = array(0, c(T, n, ndraw))
  eta_draws = array(0, c(T, k.star, ndraw))
  sig_draws = array(0, c(n, ndraw))
  ma <- vector("list", n)
  for(j in 1:n){
    ma[[j]] <- array(0, dim = c(length(unique(Y[, j])), 2, ndraw))
  }

  # ----------------------------------------------------------------------------
  # off to war
  # ----------------------------------------------------------------------------
  
  draw = 0
  
  for (m in 1:M) {
    
    # --------------------------------------------------------------------------
    # draw from p(G, S | ...)
    # --------------------------------------------------------------------------
    
    GW = sample_conjugate_posterior_varp(eta, var_lag, inclconst, GS_prior, TRUE)
    G = t(GW$B) # G needs to oriented as k.star x k.star x var_lag
    S = GW$S
    
    # JZ: might help to store these, since I recompute them for the MA
    Gam <- stationary_var1_covariance(G, S)
    
    # --------------------------------------------------------------------------
    # draw from p(factors | ...)
    # --------------------------------------------------------------------------
    
    # JZ: This assumes a VAR(1).
    
    model = SSModel(Z ~ -1 + SSMcustom(Z = Lambda,
                                       T = G,
                                       R = diag(k.star),
                                       Q = S,
                                       a1 = numeric(k.star),
                                       P1 = Gam,
                                       P1inf = matrix(0, k.star, k.star)),
                    H = diag(Sigma.diag))
    eta = simulateSSM(model, type = "states", nsim = 1)[, , 1]
    
    # --------------------------------------------------------------------------
    # draw from p(Lambda, etc | ...)
    # --------------------------------------------------------------------------
    
    # JZ: this is all copied from Joe. No clue what it's doing.
    
    # 6) Sample lambda_jh
    cp.eta = crossprod(eta)
    for(j in 1:p){
      chQj = chol(diag(phi.jh[j,]*tau.h, k.star) + cp.eta/Sigma.diag[j])
      lj = crossprod(eta, Z[,j])/Sigma.diag[j]
      Lambda[j,] = backsolve(chQj,forwardsolve(t(chQj), lj) + rnorm(k.star))
    }
    
    # 7) Sample the error variances
    eps =  Z - (tcrossprod(eta, Lambda))
    
    Sigma.diag = apply(eps, 2, function(x) 1/rgamma(n = 1, shape = a.sigma + T/2,
                                                    rate = b.sigma + 1/2*sum(x^2)))
    
    # 8) sample phi.jh
    phi.jh = matrix(rgamma(n = p*k.star, shape = (nu + 1)/2,
                           rate = (nu + Lambda^2*matrix(rep(tau.h, each = p), nr = p))/2), nr = p) #for(h in 1:k.star){for(j in 1:p) phi.jh[j,h] = rgamma(n = 1, shape = (nu + 1)/2, rate = (nu + Lambda[j,h]^2*tau.h[h])/2)
    
    # 9) sample tau.h via delta.h
    delta.h = sampleMGP(theta.jh = sqrt(phi.jh)*Lambda, delta.h = delta.h, a1 = a1, a2 = a2)
    tau.h = cumprod(delta.h)
    
    # --------------------------------------------------------------------------
    # draw from p(Z | ...)
    # --------------------------------------------------------------------------
    
    for(j in 1:p){
      
      muj = Lambda[j,]%*%t(eta)
      sdj = sqrt(Sigma.diag[j])
      
      for(r in 1:Rlevels[j]){
        # JZ: Need to understand better what's going on here. Why do we need this R stuff?
        #     why are we suppressing any warnings?
        ir <- (1:T)[R[,j] == r & !is.na(R[,j])]
        lb <- suppressWarnings(max(Z[R[,j] == r-1,j], na.rm = TRUE))
        ub <- suppressWarnings(min(Z[R[,j] == r+1,j], na.rm = TRUE))
        Z[ir,j] = qnorm(runif(length(ir), pnorm(lb, muj[ir],sdj),pnorm(ub,muj[ir],sdj)),muj[ir],sdj)
        
      }
      
    }

    if ( retain_draw(m, burn, thin) ) {
      draw = draw + 1
      
      Lambda_draws[, , draw] = Lambda
      G_draws[, , draw] = G
      S_draws[, , draw] = S
      Z_draws[, , draw] = Z
      eta_draws[, , draw] = eta
      sig_draws[, draw] = Sigma.diag
      
      # Marginal covariance of each z[t]
      Zcov = Lambda %*% Gam %*% t(Lambda) + diag(Sigma.diag)
      
      for(j in 1:n){
        yj = sort(unique(Y[, j]))
        u = pnorm(Z[, j], mean = 0, sd = sqrt(Zcov[j, j]))
        ma[[j]][, 1, draw] <- yj
        ma[[j]][, 2, draw] <- sapply(yj, margin.adjustment, Y[, j], u)
      }
      
    }
  }
  
  draws = list(Lambda = Lambda_draws,
               G = G_draws,
               S = S_draws,
               Z = Z_draws,
               eta = eta_draws,
               sig = sig_draws,
               ma = ma)
  
  return(draws)
  
}

DGFC.forecast <- function(H, draws){
  # ----------------------------------------------------
  # dims
  # ----------------------------------------------------
  
  T <- dim(draws$Z)[1]
  n <- dim(draws$Z)[2]
  ndraw <- dim(draws$Z)[3]
  k.star = dim(draws$S)[1]
  
  # ----------------------------------------------------
  # preallocate storage
  # ----------------------------------------------------
  
  Ypred <- array(0, dim = c(H, n, ndraw))
  
  # ----------------------------------------------------
  # Go!
  # ----------------------------------------------------
  
  for(m in 1:ndraw){
    
    # ------------------------------------
    # unpack draws
    # ------------------------------------
    
    # JZ: VAR(1) heavily hard-coded here
    
    G <- draws$G[, , m]
    Garr <- array(G, dim = c(k.star, k.star, 1))
    S <- draws$S[, , m]
    Gam <- stationary_var1_covariance(G, S)
    L <- draws$Lambda[, , m]
    M <- diag(draws$sig[, m])
    eta <- draws$eta[, , m]
    Zcov = L %*% Gam %*% t(L) + M
    
    new_eta = simulate_stationary_var(H, numeric(k.star), Garr, S, Y0 = matrix(eta[T, ], 1, k.star))
    new_Z = new_eta %*% t(L) + mvrnorm(n = H, mu = numeric(n), Sigma = M)
    
    for(j in 1:n){
      u = pnorm(new_Z[, j], mean = 0, sd = sqrt(Zcov[j, j]))
      xj = draws$ma[[j]][, 1, m]
      Fxj = draws$ma[[j]][, 2, m]
      Fxj[length(xj)] = 1
      Ypred[, j, m] = sapply(u, inverse_cdf, xj, Fxj)
    }
    
  }
  
  return(Ypred)
  
}

sampleMGP <- function(theta.jh, delta.h, a1 = 3, a2 = 5){
  # Description: Samples MGP parameters for factor model
  #
  # Inputs:
  # theta.jh, delta.h: local and global shrikage parameters
  #
  # Outputs
  # new sample of delta.h
  
  # Store the dimensions locally
  p = nrow(theta.jh); K = ncol(theta.jh)
  
  # Sum over the (squared) replicates:
  sum.theta.l = colSums(theta.jh^2)
  
  # h = 1 case is separate:
  tau.not.1 = cumprod(delta.h)/delta.h[1]
  delta.h[1] = rgamma(n = 1, shape = a1 + p*K/2,
                      rate = 1 + 1/2*sum(tau.not.1*sum.theta.l))
  # h > 1:
  if(K > 1){for(h in 2:K){
    tau.not.h = cumprod(delta.h)/delta.h[h]
    delta.h[h] = rgamma(n = 1, shape = a2 + p/2*(K - h + 1),
                        rate = 1 + 1/2*sum(tau.not.h[h:K]*sum.theta.l[h:K]))
  }}
  delta.h #list(tau.h = cumprod(delta.h), delta.h = delta.h)
}