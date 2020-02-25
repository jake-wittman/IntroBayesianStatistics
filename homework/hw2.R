
# Load libraries ----------------------------------------------------------

pacman::p_load(tidyverse,
               rstan, 
               brms,
               bayesplot)
conflict_prefer("Position", "base")
conflict_prefer("extract", "rstan")
conflict_prefer("filter", "stats")
conflict_prefer("lag", "stats")

# data --------------------------------------------------------------------
intersections <- c(1:18)
streets <- c(rep(1, 10), rep(0, 8))
bikes <- c(16, 9, 10, 13, 19, 20, 18, 17, 35, 55, 12, 1, 2, 4, 9, 7, 9, 8)


bike <- data.frame(
  intersections = intersections,
  streets = streets,
  bikes = bikes
)

bike_street <- bike %>% dplyr::filter(streets == 1)
bike_no_streets <- bike %>% dplyr::filter(streets == 0)


# Stan model --------------------------------------------------------------


# 3a ----------------------------------------------------------------------


bike_model <- rstan::stan_model(file = "homework/hw2.stan")

out_model1_streets <- sampling(
  object = bike_model,
  data = list(y_s = bike_street$bikes,
              n_s = nrow(bike_street),
              y = bike_no_streets$bikes,
              n = nrow(bike_no_streets)),
  warmup = 0,
  iter = 4000, 
  chains = 4,
  seed = 1,
  cores = 4,
  control = list(adapt_delta = 0.9)
)

traceplot(out_model1_streets, par = c("alpha", "beta", "alpha_s", "beta_s"))


# 3b ----------------------------------------------------------------------

out_model2_streets <- sampling(
  object = bike_model,
  data = list(y_s = bike_street$bikes,
              n_s = nrow(bike_street),
              y = bike_no_streets$bikes,
              n = nrow(bike_no_streets)),
  warmup = 4000,
  iter = 6000, 
  chains = 4,
  seed = 1,
  cores = 4,
  control = list(adapt_delta = 0.9)
)

traceplot(out_model2_streets, par = c("alpha", "beta", "alpha_s", "beta_s"))

stan_hist(out_model2_streets, pars = c("posterior_mean_difference", "pred_y"))

mcmc_areas(out_model2_streets,
           pars = c("posterior_mean_difference"),
           prob = 0.95) +
  labs(title = "Posterior distribution for difference in population level means",
       subtitle = "with median and 95% credible intervals",
       x = "Difference in population level means",
       y = "Posterior density") +
  scale_x_continuous(limits = c(-50, 100),
                     labels = seq(-50, 100, 25),
                     breaks = seq(-50, 100, 25))
  
mcmc_areas(out_model2_streets,
           pars = c("pred_y"),
           prob = 0.95) +
  labs(title = "Predictive distribution for the number of bikes on a new bike route intersection",
       subtitle = "with median and 95% credible intervals",
       x = "Predicted number of bikes at a new bike route intersection",
       y = "Posterior density") +
  scale_x_continuous(limits = c(0, 100),
                     labels = seq(0, 100, 25),
                     breaks = seq(0, 100, 25))
