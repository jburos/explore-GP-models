library(rstan)
library(tidyverse)
library(ggplot2)
library(lazyeval)
source('sim_data.function.R')

gp2 <- "
data {
  int<lower=1> N;
  vector[N] x;
  vector[N] y;
}
transformed data {
  vector[N] mu;
  for (i in 1:N) mu[i] = 0;
}
parameters {
  real<lower=0> eta_sq;
  real<lower=0> inv_rho_sq;
  real<lower=0> sigma_sq;
}
transformed parameters {
  real<lower=0> rho_sq;
  rho_sq = inv(inv_rho_sq);
}
model {
  matrix[N, N] Sigma;
  // off-diagonal elements
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      Sigma[i, j] = eta_sq * exp(-rho_sq * pow(x[i] - x[j],2));
      Sigma[j, i] = Sigma[i, j];
    }
  }
  // diagonal elements
  for (k in 1:N)
    Sigma[k, k] = eta_sq + sigma_sq; // + jitter
  eta_sq ~ cauchy(0, 5);
  inv_rho_sq ~ cauchy(0, 5);
  sigma_sq ~ cauchy(0, 5);
  y ~ multi_normal(mu, Sigma);
}
"

df <- sim_data(n=1000, sigma_sq=0.1, eta_sq=10, rho_sq=0.1)
ggplot(df, aes(x=x, y=y)) + geom_line()

standata <- list(N=nrow(df), x=df$x, y=df$y)
fit <- rstan::stan(model_code=gp2, data=standata)
print(fit)


