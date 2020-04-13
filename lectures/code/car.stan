data {
    int<lower=0> N;
    int<lower=0> N_edges;
    int<lower=1, upper=N> node1[N_edges];  // node1[i] adjacent to node2[i]
    int<lower=1, upper=N> node2[N_edges];  // and node1[i] < node2[i]
    
    int<lower=0> y[N];      // count outcomes
    vector[N] x;            // predictor
    vector<lower=0>[N] E;   // exposure
}
transformed data {
    vector[N] log_E = log(E);
}
parameters {
    real beta0;         // intercept
    real beta1;         // slope
    
    real<lower=0> tau_theta;    // precision of heterogeity effects
    real<lower=0> tau_phi;      // precision of spatial effects
    
    vector[N] theta;        // heterogeneity effects
    vector[N] phi;          // spatial effects
}
transformed parameters {
    real<lower=0> sigma_theta =  inv(sqrt(tau_theta));  // SD of the prior normal on heterogeneity effects 
    real<lower=0> sigma_phi = inv(sqrt(tau_phi));       // SD of the IAR prior on spatial effects
    
    vector[N] eta = theta * sigma_theta + phi * sigma_phi;    // marginal errors
    
    real<lower=0> sd_h = sd(theta) * sigma_theta;     // marginal SD of heterogeneity effects
    real<lower=0> sd_s = sd(phi) * sigma_phi;       // marginal SD of spatial effects
    real<lower=0> alpha = sd_s/(sd_h+sd_s);
}
model {
    y ~ poisson_log(log_E + beta0 + beta1 * x + phi * sigma_phi + theta * sigma_theta);
    
    // the following computes the prior on phi on the unite scale
    target += -0.5 * dot_self(phi[node1] - phi[node2]);
    // soft sum-to-zero constraint on phi
    sum(phi) ~ normal(0, 0.001 * N); 
    
    beta0 ~ normal(0, 100);
    beta1 ~ normal(0, 100);
    theta ~ normal(0, 1);
    tau_theta ~ gamma(3.2761, 1.81);
    tau_phi ~ gamma(1, 1);
}
generated quantities{
    vector[N] mu = exp(log_E + beta0 + beta1 * x + phi * sigma_phi + theta * sigma_theta);
}

