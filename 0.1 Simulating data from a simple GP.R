
library(rstan)
library(tidyverse)
library(ggplot2)
library(lazyeval)

gp1 <- "
data {
  int<lower=1> N;
  real x[N];
}
transformed data {
  matrix[N,N] L;
  vector[N] mu;
  cov_matrix[N] Sigma;
  real<lower=0> eta_sq;
  real<lower=0> rho_sq;
  real<lower=0> sigma_sq;
  eta_sq = 1;
  rho_sq = 1; 
  sigma_sq = 0.1;
  for (i in 1:N) 
    mu[i] = 0;
  for (i in 1:N) 
    for (j in 1:N)
      Sigma[i, j] = eta_sq * exp(-rho_sq * pow(x[i] - x[j], 2)) + (i==j ? sigma_sq : 0.0);
  L = cholesky_decompose(Sigma);
}
parameters {
  vector[N] z;
}
model {
  z ~ normal(0, 1);
}
generated quantities {
  vector[N] y;
  y = mu + L * z;
}
"

#' sample random groups from a df
#' @param group_by str columns to group by
#' @param n_groups int how many groups to select
sample_group_ <- function(df, group_by, n_groups) {
  df %>% semi_join(
    df %>% distinct_(group_by) %>% dplyr::sample_n(n_groups) %>% dplyr::select_(group_by)
  )
}
sample_group <- function(df, group_by, n_groups) {
  sample_group_(df, lazyeval::f_capture(group_by), n_groups)
}

N <- 20
sd1 <- list(N=N, x=seq_len(N)+rnorm(mean=0, sd=0.1, n=N))
sim1 <- rstan::stan(model_code=gp1, data=sd1)
y_sim <- rstan::extract(sim1, 'y')$y
y_df <- tbl_df(y_sim) %>%
  dplyr::mutate(iter = row_number()) %>%
  tidyr::gather(key='x_name', value='y', -iter) %>%
  dplyr::mutate(x = as.integer(stringr::str_replace(x_name, pattern='V', replacement='')))

## plot random draw from this gaussian process 
ggplot(y_df %>% sample_group(group_by=iter, n_groups=1), aes(x=x, y=y, colour=iter, group=iter)) + geom_line()

