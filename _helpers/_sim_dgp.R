# ==============================================================================
# DGP factory for the margin-comparison simulation (new Section 6.2)
#
# Place in _helpers/ and source after _helpers.R
#
# make_dgp(dgp) returns ONE draw of parameters as:
#   $sim(T), function giving a T x n matrix, burn-in already discarded
#   $F, list of n true marginal CDFs
#   $discrete, logical vector, length n
#   $rho, spectral radius of the AR companion (NA for varch)
#   $n, number of variables
#
# Deliberate differences from the setup blocks in
# simulation_wrongmodel_consistency.R and simulation_wrongmodel_forecasting.R:
#
#  1. Parameters are redrawn on every call rather than fixed by set.seed(1),
#     so results describe the DGP family rather than one draw from it
#
#  2. A burn-in is discarded. simulate_stationary_varma draws its p initial
#     lags iid from the stationary marginal, which gets the marginals right but
#     not the autocorrelation between them, and the transient decays at rho^t. 
#     At T = 25 that transient is most of the sample, which is exactly the cell 
#     this section is about. VARCH does not strictly need it, since y0 is drawn 
#     from the stationary distribution and the kernel preserves it, but it is kept uniform
#
#  3. The CDF list is built with lapply. simulation_wrongmodel_consistency.R
#     gives all n closures the same i, which equals n after the loop. Harmless there
#     because the paper panels set i = 2; not harmless here
# ==============================================================================

make_dgp <- function(dgp, burnin = 200) {

  if (dgp == "varma") {

    n <- 2; p <- 3; q <- 6
    m       <- rnorm(n)
    AR      <- simulate_nonexplosive_var_params(n, p, numeric(n*n*p),
                                                0.1 * diag(n*n*p))
    ARarray <- array(c(AR), dim = c(n, n, p))
    MAarray <- array(rnorm(n*n*q), dim = c(n, n, q))
    S       <- riwish(n + 1, diag(n))
    par     <- varma_stationary_params(m, ARarray, MAarray, S)

    sim <- function(T) {
      Y <- simulate_stationary_varma(T + burnin, m, ARarray, MAarray, S)
      Y[(burnin + 1):(burnin + T), , drop = FALSE]
    }

    Flist <- lapply(seq_len(n), function(i) {
      mu <- par$VARMAmean[i] # forced in this frame, so safe
      sg <- sqrt(par$VARMAcov[i, i])
      function(x) pnorm(x, mean = mu, sd = sg)
    })

    rho      <- max(abs(eigen(companion(AR))$values))
    discrete <- rep(FALSE, n)

  } else if (dgp == "varch") {

    n <- 2
    v <- n + 2 + rexp(1)
    A <- riwish(n + 1, diag(n))
    par <- varch_stationary_params(v, A)

    sim <- function(T) {
      Y <- simulate_varch(T + burnin, v, A)
      Y[(burnin + 1):(burnin + T), , drop = FALSE]
    }

    Flist <- lapply(seq_len(n), function(i) {
      sc <- sqrt(par$VARCHscale[i, i])
      df <- par$VARCHdf
      function(x) pt(x / sc, df = df)
    })

    rho      <- NA_real_ # no VAR companion in this DGP
    discrete <- rep(FALSE, n)

  } else if (dgp == "varma_copula") {

    n <- 2; p <- 1; q <- 1

    Finv  <- list(function(u) qpois(u, lambda = 5),
                  function(u) qst(u, alpha = 2, nu = 3))
    Flist <- list(function(x) ppois(x, lambda = 5),
                  function(x) pst(x,  alpha = 2, nu = 3))

    AR      <- simulate_nonexplosive_var_params(n, p, numeric(n*n*p),
                                                0.1 * diag(n*n*p))
    ARarray <- array(c(AR), dim = c(n, n, p))
    MAarray <- array(rnorm(n*n*q), dim = c(n, n, q))
    S       <- riwish(n + 1, diag(n))

    sim <- function(T) {
      Y <- simulate_varma_copula(T + burnin, Finv, ARarray, MAarray, S)
      Y[(burnin + 1):(burnin + T), , drop = FALSE]
    }

    rho      <- max(abs(eigen(companion(AR))$values))
    discrete <- c(TRUE, FALSE) # Poisson(5), skew-t(3, 0, 1, 2)

  } else if (dgp == "dgfc") {

    # Correctly-specified case. Only used if John wants a fourth DGP; the
    # parameter choices follow dgfc_in_wrongmodel_sim.R

    Finv  <- list(function(u) qst(u, alpha = 2, nu = 3),
                  function(u) qpois(u, lambda = 5))
    Flist <- list(function(x) pst(x,  alpha = 2, nu = 3),
                  function(x) ppois(x, lambda = 5))

    n <- 2; k <- n
    Mv <- diag(n)
    L  <- rbind(diag(k), matrix(0, n - k, k))
    AR <- 0.5 * diag(k)
    G  <- array(c(AR), dim = c(k, k, 1))
    S  <- diag(k)

    sim <- function(T) {
      Y <- simulate_dgfc(T + burnin, Finv, G, S, L, Mv)
      Y[(burnin + 1):(burnin + T), , drop = FALSE]
    }

    rho      <- max(abs(eigen(companion(AR))$values))
    discrete <- c(FALSE, TRUE) # skew-t, Poisson

  } else {
    stop("unknown dgp: ", dgp)
  }

  list(sim = sim, F = Flist, discrete = discrete, rho = rho,
       n = n, dgp = dgp)
}
