library(MASS)
library(MCMCpack)
set.seed(8675309)

n <- 2000
x <- 0.75
dependent = TRUE

if(dependent == TRUE){
  df = n + 2
  #C = cov2cor(riwish(df, (df - n - 1) * diag(n)))
  C <- 0.5 * matrix(1, n, n) + (1 - 0.5) * diag(n)
}else{
  C <- diag(n)
}

z <- mvrnorm(1, numeric(n), C)
y <- pnorm(z)

hist(z, breaks = "Scott", freq = FALSE)
lines(seq(-4, 4, length.out = 1000), dnorm(seq(-4, 4, length.out = 1000)), col = "red")

hist(y, breaks = "Scott", freq = FALSE)

v <- y * (y <= x)
m <- cummax(v)
plot(1:n, m, type = "l")
abline(h = x)