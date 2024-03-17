retain_draw <- function(stage, burn, thin){
  postburn = stage > burn
  unthinned = ((stage - burn) %% thin) == 0
  retained = postburn & unthinned
  return( retained )
}

rinvgamma <- function(n, a, b){
  return( 1 / rgamma(n, a, rate = b) )
}

margin.adjustment <- function(x, y, u){
  # compute the margin adjustment estimate of F(x)
  # for a grid of x-vals, try: sapply(x_grid, margin.adjustment, y, u)
  #
  #  x: scalar 
  #  y: T-vector of input values
  #  u: T-vector of output values u = F(y), where F is unknown
  max(c(u[y <= x], u[y == min(y)]))
}

inverse_cdf <- function(u, x, Fx){
  min(x[u <= Fx])
}