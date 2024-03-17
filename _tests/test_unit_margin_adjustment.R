# ==============================================================================
# load resources
# ==============================================================================

source("_helpers/_helpers.R")
source("_original_code/Helpers.R")

# ==============================================================================
# fake data
# ==============================================================================

n = 10
Y = cbind(rgamma(n, shape = 1, rate = 1), rpois(n, lambda = 1))
p = ncol(Y)
k.star = p

# ==============================================================================
# latent N(0, 1) data with the same ranks as the Y
# ==============================================================================

Z = matrix(0, n, p)
for(i in 1:p){
  Z[, i] = sort(rnorm(n))[rank(Y[, i], ties.method = "random")]
}

# ==============================================================================
# location of each y-val in the sorted list of unique values
# ==============================================================================

R<- NULL
for(j in 1:p) { R<-cbind(R, match(Y[,j],sort(unique(Y[,j])))) }

# ==============================================================================
# number of unique values per variable
# ==============================================================================

Rlevels<-apply(R,2,max,na.rm=TRUE)

# ==============================================================================
# Joe's implementation of the margin adjustment, which takes into account all
# the mixture stuff
# ==============================================================================

x = 1
bounds <- get_bound(R[, x], Rlevels[x], Z[, x])

z = rep(1, n)
pi_h = 1
Delta = list(diag(k.star))
mu = list(numeric(k.star))
Lambda = diag(p)
Sigma.diag = numeric(p)
alpha = numeric(p)

wts <- compute_mixprobs(bounds, Y[,x], pi_h, z, Delta, mu, alpha = alpha, Lambda, Sigma.diag, x)

# ==============================================================================
# calculate the margin adjustment directly from Joe's definition
# ==============================================================================

y <- Y[, x]
z <- Z[, x]
u <- pnorm(z)
Sy <- sort(unique(y))

my_ma <- sapply(Sy, margin.adjustment, y, u)

# ==============================================================================
# compare
# ==============================================================================

VS <- cbind(wts[, 3], my_ma)
colnames(VS) <- c("JF", "JZ1")
rownames(VS) <- NULL
VS

