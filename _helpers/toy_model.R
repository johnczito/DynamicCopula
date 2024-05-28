DGC.mcmc.0 <- function(Y, 
                     v0 = ncol(Y) + 1, V0 = diag(ncol(Y)), 
                     G0 = matrix(0, ncol(Y), ncol(Y)), invO0 = diag(ncol(Y)), 
                     ndraw = 1000, burn = 0, thin = 1){
  
  # ----------------------------------------------------------------------------
  # dimensions of things
  # ----------------------------------------------------------------------------
  
  n = nrow(Y)
  p = ncol(Y)
  M = burn + thin * ndraw
  
  # ----------------------------------------------------------------------------
  # preallocate storage
  # ----------------------------------------------------------------------------
  
  G_draws = array(0, c(p, p, ndraw))
  S_draws = array(0, c(p, p, ndraw))
  Z_draws = array(0, c(n, p, ndraw))
  
  ma <- vector("list", p)
  for(j in 1:p){
    ma[[j]] <- array(0, dim = c(length(unique(Y[, j])), 2, ndraw))
    for(k in 1:ndraw){
      ma[[j]][, 1, k] = sort(unique(Y[, j]))
    }
  }
  
  # ----------------------------------------------------------------------------
  # set prior
  # ----------------------------------------------------------------------------
  
  GS_prior = list(v = v0, P = V0, B = G0, invO = invO0)
  
  # ----------------------------------------------------------------------------
  # initialize
  # Note: you have to initialize with a Z that satisfies the rank constraint
  # ----------------------------------------------------------------------------
  
  Z = matrix(0, n, p)
  for(j in 1:p){
    Z[, j] = scale(Y[, j])
  }
  
  # ----------------------------------------------------------------------------
  # off to war
  # ----------------------------------------------------------------------------
  
  draw = 0
  
  for (m in 1:M) {
    
    # --------------------------------------------------------------------------
    # draw from p(G, S | ...)
    # --------------------------------------------------------------------------
    
    GW = sample_conjugate_posterior_varp(Z, 1, FALSE, GS_prior, TRUE)
    G = t(GW$B) # G needs to oriented as k.star x k.star x var_lag
    S = GW$S
    
    O = full_var1_covariance(n, G, S)
    
    GnSn = normalize_var1_params(G, S)
    
    Gn = GnSn$G
    Sn = GnSn$S
  
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # draw from p(Z | C, Y)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # JZ: this should be equivalent to Hoff in the all-continuous case at least.
    # Numerical problems with heavy-tailed data. Truncated normal sampling gets 
    # all fucked up
    
    for(j in 1:p){
      for(i in 1:n){
        l = (i - 1) * p + j
        zj = Z[, j]
        yj = Y[, j]
        y = Y[i, j]
        zL = suppressWarnings(max(zj[yj < y]))
        zU = suppressWarnings(min(zj[y < yj]))
        O_nl_inv = solve(O[-l, -l])
        s = c(sqrt( O[l, l] - O[l, -l] %*% O_nl_inv %*% O[-l, l] ))
        mu = c(c(t(Z))[-l] %*% t(O[l, -l] %*% O_nl_inv))
        bL = pnorm((zL - mu) / s)
        bU = pnorm((zU - mu) / s)
        u = runif(1, min = bL, max = bU)
        Z[i, j] = mu + s * qnorm(u)
      }
    }
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # store?
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    if ( retain_draw(m, burn, thin) ) {
      
      draw = draw + 1
      
      G_draws[, , draw] = Gn
      S_draws[, , draw] = Sn
      Z_draws[, , draw] = Z
      
      for(j in 1:p){
        #yj = sort(unique(Y[, j]))
        u = pnorm(Z[, j], mean = 0, sd = sqrt(O[j, j])) #ecdf(Z[, j])(Z[, j])
        #ma[[j]][, 1, draw] <- yj
        #ma[[j]][, 2, draw] <- sapply(yj, joe.margin.adjustment, Y[, j], Z[, j])
        ma[[j]][, 2, draw] <- sapply(ma[[j]][, 1, draw], margin.adjustment, Y[, j], u)
      }
      
    }
  }
  
  draws = list(G = G_draws, S = S_draws, Z = Z_draws, ma = ma)
  
  return(draws)
}

DGC.mcmc <- function(Y, prior = NULL, init = NULL, ndraw = 1000, burn = 0, thin = 1){
  # ----------------------------------------------------------------------------
  # dimensions
  # ----------------------------------------------------------------------------
  
  # JZ: need to decide on this n vs. p business
  
  T = nrow(Y)
  n = ncol(Y)
  p = n
  k = n
  k.star = n
  M = burn + thin * ndraw
  
  # ----------------------------------------------------------------------------
  # create or unpack prior hyperparameters
  # ----------------------------------------------------------------------------
  
  if(is.null(prior)){
  
    # error Variance
    a.sigma = 1
    b.sigma = 0.3
    
    # VAR params
    var_lag = 1
    inclconst = FALSE
    v0 = n + 1
    P0 = diag(n)
    B0 = matrix(0, n, n)
    invO0 = diag(n)
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
    
    Lambda = diag(n)
    
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
  
  G_draws = array(0, c(n, n, ndraw))
  S_draws = array(0, c(n, n, ndraw))
  Z_draws = array(0, c(T, n, ndraw))
  eta_draws = array(0, c(T, n, ndraw))
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
    
    # 7) Sample the error variances
    eps =  Z - (tcrossprod(eta, Lambda))
    
    Sigma.diag = apply(eps, 2, function(x) 1/rgamma(n = 1, shape = a.sigma + T/2,
                                                    rate = b.sigma + 1/2*sum(x^2)))
    
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
      
      O = full_var1_covariance(2, G, S) + diag(2) %x% diag(Sigma.diag)
      C = cov2cor(O)
      C0 = C[1:n, 1:n]
      C1 = C[(n+1):(2*n), 1:n]
      Gn = C1 %*% solve(C0)
      Sn = matrix((diag(n*n) - Gn %x% Gn) %*% c(C0), n, n)
      
      G_draws[, , draw] = Gn
      S_draws[, , draw] = Sn
      #Z_draws[, , draw] = Z
      #eta_draws[, , draw] = eta
      #sig_draws[, draw] = Sigma.diag
      
      # Marginal covariance of each z[t]
      #Zcov = Lambda %*% Gam %*% t(Lambda) + diag(Sigma.diag)
      
      for(j in 1:n){
        yj = sort(unique(Y[, j]))
        u = pnorm(Z[, j], mean = 0, sd = sqrt(O[j, j]))
        ma[[j]][, 1, draw] <- yj
        ma[[j]][, 2, draw] <- sapply(yj, margin.adjustment, Y[, j], u)
      }
      
    }
  }
  
  draws = list(G = G_draws,
               S = S_draws,
               ma = ma)
  
  return(draws)
  
}