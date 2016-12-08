
library(rstan)
library(ggplot2)
library(tidyverse)
source('sim_data.function.R')

eta_sq <- seq(from=1, to=6, by=2)
rho_sq <- seq(from=1, to=6, by=2)
sigma_sq <- seq(from=0.1, to=1.1, by=0.2)

params <- purrr::cross_d(list(n=1000, eta_sq=eta_sq, rho_sq=rho_sq, sigma_sq=sigma_sq))

sims <- params %>% 
  purrr::invoke_rows(.f = sim_data, .collate='rows', .labels = TRUE)
saveRDS(sims, 'simulated_data.Rds')

p <- ggplot(sims, aes(x=x, y=y, group=.row, color=factor(sigma_sq))) + geom_line() + facet_grid(eta_sq~rho_sq)
ggsave(p, filename='simulated_data.png', width=20, height=20, units='in')
