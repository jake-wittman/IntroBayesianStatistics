

data {
  int<lower=0> n;
  int y[n];
  int<lower = 0> k;
}

parameters {
  real theta;
  real lambda;
  real b1;
  real b2;
}

model {
  for (i in 1:k) {
    b1 ~ inv_gamma(1, 1);
    theta ~ gamma(0.5, b1);
    y[i] ~ poisson(theta);
  }
  
  for(j in (k +1): n){
    b2 ~ inv_gamma(1, 1);
    lambda ~ gamma(0.5, b2);
    y[j] ~ poisson(lambda);
  }
}

generated quantities{
  real R;
  R = theta / lambda;
}

