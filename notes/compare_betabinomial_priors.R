# Libraries
pacman::p_load(tidyverse)

# Consider three priors
# Beta(110, 7), Beta(110 * 0.5, 7 * 0.5) and Beta(110 * 0.1, 7 * 0.1)
#Assuming we observe X = 45 out of N = 50 trials
# Plot priors and corresponding posteriors
# Calculate 95% credible intervals
# make conclusions

# simulate distributions
x <- 500000
dat <- data.frame(
  prior1 = rbeta(x, 110, 7),
  prior2 = rbeta(x, 110 * 0.5, 7 * 0.5),
  prior3 = rbeta(x, 110 * 0.1, 7 * 0.1),
  posterior1 = rbeta(x, 45 + 110, 50 - 45 + 7),
  posterior2 = rbeta(x, 45 + (110 * 0.5), 50 - 45 + (7 * 0.5)),
  posterior3 = rbeta(x, 45 + (110 * 0.1), 50 - 45 + (7 * 0.1))
)

dat_long <- pivot_longer(dat, cols = everything())
dat_ci <- data.frame(
  ci1 = c(qbeta(0.025, 45 + 110, 50 - 45 + 7),
           qbeta(0.975, 45 + 110, 50 - 45 + 7)),
  ci2 = c(qbeta(0.025, 45 + (110 * 0.5), 50 - 45 + (7 * 0.5)),
           qbeta(0.975, 45 + (110 * 0.5), 50 - 45 + (7 * 0.5))),
  ci3 = c(qbeta(0.025, 45 + (110 * 0.1), 50 - 45 + (7 * 0.1)),
           qbeta(0.025, 45 + (110 * 0.1), 50 - 45 + (7 * 0.1)))
) #%>% 
  pivot_longer(cols = everything())

dat_long %>% 
ggplot(aes(value, group = name, color = name)) +
  geom_density() +
  xlim(0.7, 1) +
  geom_vline

# The graph could be better if I had organized my simulated data a bit better
# but the credible intervals show we can concluded theta > 0.85 (which is our
# hypothesis we are testing) given a very infomrative and a mildly informative
# prior. The non-informative prior ends up making it so we can't draw that claim
# at 95% CI.

# I should figure out how to use dbeta to make these graphs in ggplot rather than
# simulations with rbeta.


