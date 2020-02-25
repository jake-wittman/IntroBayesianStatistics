
# Libraries ---------------------------------------------------------------

pacman::p_load(tidyverse,
               invgamma)


# Numerical approximation of a posterior dist. ----------------------------

# Pseudocode
# Draw a sample of sigma
# draw a sample of mu given sigma and sample mean 
# Draw a sample of y* from N(drawn mu, drawn sigma)
set.seed(1)

M <- 10000
y <- 15.4
s <- 7.6
s2 <- s^2
n <- 109
dat <- data.frame(
  sigma_sq = 1 / rgamma(M, shape = (n - 1) / 2, rate = ((n - 1) * s^2) / 2)
)
dat <- dat %>% 
  mutate(mu = rnorm(M, mean = y, sd = sqrt(sigma_sq / n))) %>% 
  mutate(y_star = rnorm(M, mu, sqrt(sigma_sq)))

# Marginal posterior of mu
ggplot(dat, aes(mu)) +
  geom_density()

# Marginal posterior of y star
ggplot(dat, aes(y_star)) +
  geom_density()

# See if loop gives similar answers



sigma2_post <- rep(NA, M)
mu_post <- rep(NA, M)
ystar_post <- rep(NA, M)
for (i in 1:M) {
  sigma2_post[i] <- 1 / rgamma(1, (n -1) / 2, rate = (n-1) * s2/2)
  mu_post[i] <- rnorm(1, y, sqrt(sigma2_post[i]/ n))
  ystar_post[i] <- rnorm(1, mu_post[i], sqrt(sigma2_post[i]))
}

dat <- data.frame(
  sigma2 = sigma2_post,
  mu = mu_post,
  ystar = ystar_post
)

# Marginal posterior of mu
ggplot(dat, aes(mu)) +
  geom_density()

# Marginal posterior of y star
ggplot(dat, aes(ystar)) +
  geom_density()
