data {
    int i; 
    int<lower=0> n[i];          // number of schools 
    real y[i]; 
    real w[i]; // estimated treatment effects
}
parameters {
    real mu;
    real<lower=0> sigma;
    real<lower = 0> m;
    real a;
    real b;
    real c;
    real d;
    real g;
    real f;
    
}
transformed parameters {
    real<lower = 0> x[i];
    x[i] = (w[i] - mu) / sigma;
    
}
model {
    m ~ gamma(a, b);
    mu ~ normal(c, d);
    sigma ~ inv_gamma(g, f);
    y[i] ~ bernoulli(x[i]);
}
//model {
//    target += normal_lpdf(eta | 0, 1);
//    target += normal_lpdf(y | theta, sigma);
//}