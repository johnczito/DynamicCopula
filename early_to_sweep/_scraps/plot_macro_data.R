GDP = read.csv("datasets/GDPC1.csv", header = TRUE)
GDPPI = read.csv("datasets/GDPCTPI.csv", header = TRUE)
UNEMP = read.csv("datasets/UNRATE.csv", header = TRUE)

T <- nrow(GDP)

gdp_growth = 400*(log(GDP[2:T, 2]) - log(GDP[1:(T - 1), 2]))
gdp_inflation = 400 * (log(GDPPI[2:T, 2]) - log(GDPPI[1:(T - 1), 2]))
unemp = UNEMP[2:T, 2]

Y = cbind(gdp_growth, unemp, gdp_inflation)

T = nrow(Y)
n = ncol(Y)

qq = TRUE

par(mfcol = c(2, n), mar = c(2, 2, 2, 2))
labels = c("GDP growth", "Unemployment", "Inflation")
bottom_main = c("QQ plot vs. best fitting normal", "", "")
dates = seq(1948.25, 2023.75, by = 0.25)

set.seed(8675309)

for( i in 1:3){
  #plot(seq(1964, 2013.75, by = 0.25), Y[match(1964, dates):match(2013.75, dates), i], type = "l", main = labels[i])
  plot(dates, Y[, i], type = "l", main = labels[i], cex.main = 2)
  m = mean(Y[, i])
  s = sd(Y[, i])
  if(qq == TRUE){
    samp = rnorm(7500, m, s)
    qqplot(Y[, i], samp, pch = 19, cex = 0.15, col = "red",
           main = bottom_main[i])
    abline(a = 0, b = 1)
  }else{
    L = min(Y[, i])
    U = max(Y[, i])
    xvals = seq(L, U, length.out = 1000)
    plot(density(Y[, i]))
    lines(xvals, dnorm(xvals, m, s), col = "red")
  }
}