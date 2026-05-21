# The Dynamic Gaussian Factor Copula

These are the replication files for [Zito & Kowal (2025 arXiv)](https://arxiv.org/abs/2502.16874).
We propose a dynamic Gaussian factor copula (DGFC) model for jointly modeling a 
collection of time series that include mix of data types (discrete versus continuous)
and/or non-Gaussian distributional features (skew, multi-modality, etc).

The two main commands in the source code are

```
DGFC.mcmc(Y, k.star, prior, init, ndraw, burn, thin)
DGFC.forecast(H, draws, use_spline)
```

`DGFC.mcmc` runs the Gibbs sampler (Algorithm 3) for approximating the pseudo-posterior
in the DGFC, and given the posterior draws, `DGFC.forecast` simulates the model forward
to forecast `H` steps into the future (Algorithm 2).

