
# Libraries ---------------------------------------------------------------

pacman::p_load(rstan,
               tidyverse,
               rstanarm,
               loo)
conflict_prefer("extract", "rstan")
conflict_prefer("filter", "stats")
conflict_prefer("lag", "stats")
conflict_prefer("loo", "loo")
# Data --------------------------------------------------------------------

Y <- matrix(
  c(
    0.814,
    -0.203,
    -0.133,
    NA,
    -0.715,
    0.739,
    0.118,
    NA,
    NA,
    0.271,
    NA,
    -0.0023,
    -0.076,
    0.651,
    -0.249,
    0.0026,
    NA,
    1.217,
    NA,
    NA,
    NA,
    NA,
    -0.24236,
    0.00928,
    0.8073,
    -0.51149,
    1.93893,
    1.07922,
    NA,
    0.29996,
    1.41267,-0.46985,
    0.09798,
    0.29206,
    0.19483,
    0.16531,
    -0.40556,
    NA,
    0.21807,
    NA,
    -0.54369,
    NA,
    -0.04707,
    0.23272,
    0.21767,
    -0.27662,
    0.79159,
    -0.10268,
    0.6576,
    0.0604,
    -0.27151,
    0.7048,
    0.6054,
    0.38503,
    0.29848,
    NA,
    -2.20587,
    NA,
    -0.73148,
    NA,
    0.9134,
    0.13073,-0.06594,
    -0.23161,
    1.26396,
    -0.43129,
    -0.02205,
    0.42073,
    -0.16309,
    0.60758,
    0.18718,
    0.17248,
    0.2597,
    NA,
    0.35022,
    NA,
    0.60031,
    NA,
    -0.09084,
    NA,
    NA,
    0.75204,
    -0.35662,
    0.83652,-0.16441,
    -0.11157,
    0.85996,
    -0.22899,
    NA,
    0.16011,
    NA,
    NA,
    0.14491,
    NA,
    0.04111,
    0.22188,
    0.09871,
    0.01708,
    0.35535,
    0.20278,
    0.8073,
    0.37308,
    -0.64,
    -0.01021,
    0.0813,
    1.04416,
    -0.20111,
    0.20344
  ),
  18,
  6
)

P <- matrix(
  c(
    1.55472,
    0.63796,
    0.66422,
    0.0001,
    2.22103,
    3.02457,
    1.38408,
    0.0001,
    0.0001,
    2.93207,
    0.0001,
    1.33034,
    1.19442,
    2.86302,
    2.39627,
    2.20785,
    0.0001,
    1.49815,
    0.0001,
    0.0001,
    0.0001,
    0.0001,
    1.7054,
    4.22496,
    2.52965,
    3.21071,
    1.55576,
    2.06012,
    0.0001,
    4.43682,
    0.74642,
    3.74025,
    4.20441,
    3.2145,
    5.44972,
    5.87876,
    4.61706,
    0.0001,
    1.91103,
    0.0001,
    4.26425,
    0.0001,
    2.70169,
    5.39252,
    1.96031,
    6.63744,
    1.6591,
    10.12627,
    5.15778,
    12.20163,
    9.12072,
    8.56466,
    8.09651,
    6.76019,
    3.42682,
    0.0001,
    0.73295,
    0.0001,
    4.53736,
    0.0001,
    1.97782,
    5.93378,
    2.38223,
    6.65423,
    1.27929,
    10.4889,
    5.42775,
    10.75121,
    9.06588,
    9.24314,
    9.41177,
    8.3589,
    6.15393,
    0.0001,
    2.2114,
    0.0001,
    0.66407,
    0.0001,
    2.47915,
    0.0001,
    0.0001,
    0.48528,
    3.72857,
    1.14836,
    1.1933,
    0.49845,
    3.32753,
    6.14448,
    0.0001,
    2.43433,
    0.0001,
    0.0001,
    0.64786,
    0.0001,
    3.31796,
    3.30026,
    3.31,
    12.10082,
    5.82853,
    21.92399,
    4.72651,
    6.29366,
    0.74707,
    3.86906,
    4.74962,
    4.17327,
    2.72863,
    2.50757
  ),
  18,
  6
)

J <- 6 # Number of studies
K <- 18 # Number of units

# Stan can't work with NA values
dat <- data.frame(unit = rep(rep(1:K), J),
                  study = rep(1:J, each = K),
                  Y = as.vector(Y),
                  P = as.vector(P)) %>% 
  remove_missing()



model <- stan_model(file = "lectures/code/hierarchical_model_checking.stan")

n_cores <- parallel::detectCores() - 1

fit <- rstan::sampling(
  object = model,
  data = list(
    N = nrow(dat),
    K = K,
    Y = dat$Y,
    P = dat$P,
    J = J,
    studies = dat$study,
    units = dat$unit
  ),
  chains = 3,
  warmup = 5000,
  iter = 20000,
  cores = n_cores,
  control = list(adapt_delta = 0.99)
)

traceplot(fit)
summary(fit)
pairs(fit, pars = "study_a")
plot(fit, pars = "theta_out")

# Extract log likelihoods - must be named in generated quantities
log_lik1 <- extract_log_lik(fit, parameter_name = "log_lik1")
log_lik2 <- extract_log_lik(fit, parameter_name = "log_lik2")
log_lik3 <- extract_log_lik(fit, parameter_name = "log_lik3")
log_lik4 <- extract_log_lik(fit, parameter_name = "log_lik4")

# Get LOO and Waic estimates
loo(log_lik2)
waic(log_lik2)

# Compare estimates
compare(loo(log_lik1),
        loo(log_lik2),
        loo(log_lik3),
        loo(log_lik4))

compare(waic(log_lik1),
        waic(log_lik2),
        waic(log_lik3),
        waic(log_lik4))

# BRMS attempt ------------------------------------------------------------


