data {
  int<lower=1> n; // number of observations
  vector[n] y; // observations
  real<lower=0> sigma2;
}
parameters {
  real mu; //group-level mean
  real<lower=0> tau2; // group-level variance
  vector[n] theta; // participant effects
}
model {
  target += uniform_lpdf(tau2 | 0, 1);
  target += uniform_lpdf(mu | 0, sqrt(tau2));
  target += normal_lpdf(y | theta, sqrt(sigma2));
}
