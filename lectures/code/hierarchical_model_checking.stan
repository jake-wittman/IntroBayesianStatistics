
data {
  int<lower=1> N; //number of obs.
  int<lower=1> J; //Number of studies
  int<lower=1> K; //Number of clinical units
  int<lower=0> studies[N];
  int<lower=0> units[N]; 
  real Y[N];
  real P[N];
  
}


parameters {
  //real theta[N];
  real study_a1[J];
  real unit_b1[K];
  real b_sigma1;
  real s_sigma;
  real s[N];
  
  real study_a2[J];
  real unit_b2[K];
  real b_sigma2;
  
  real study_a3[J];
  
  real unit_b4[K];
  real b_sigma4;
}

transformed parameters {
  
}


model {
  
      // Prior
    for(j in 1:J) {
      study_a1[j] ~ normal(0, 100);
      study_a2[j] ~ normal(0, 100);
      study_a3[j] ~ normal(0, 100);

    }
    
    for(k in 1:K) {
      unit_b1[k] ~ normal(0, b_sigma1 * b_sigma1);
      unit_b2[k] ~ normal(0, b_sigma2 * b_sigma2);
      unit_b4[k] ~ normal(0, b_sigma4 * b_sigma4);
    }
    
    b_sigma1 ~ uniform(0.01, 100);
    b_sigma2 ~ uniform(0.01, 100);
    b_sigma4 ~ uniform(0.01, 100);
    s_sigma ~ uniform(0.01, 100);

    
  
  // Full model
  for (n in 1:N) {
 
    Y[n] ~ normal(study_a1[studies[n]] + unit_b1[units[n]] + s[n], P[n]);
    s[n] ~ normal(0, s_sigma * s_sigma);
    
    Y[n] ~ normal(study_a2[studies[n]] + unit_b2[units[n]], P[n]);
    Y[n] ~ normal(study_a3[studies[n]], P[n]);
    Y[n] ~ normal(unit_b4[units[n]], P[n]);
  }

}

generated quantities{
  real theta_out1[N];
  real theta_out2[N];
  real theta_out3[N];
  real theta_out4[N];
  real log_lik1[N];
  real log_lik2[N];
  real log_lik3[N];
  real log_lik4[N];
  
  for (n in 1:N) {
    theta_out1[n] = study_a1[studies[n]] + unit_b1[units[n]] + s[n];
    log_lik1[n] = normal_lpdf(Y[n] | theta_out1[n], P[n]);
    
    theta_out2[n] = study_a2[studies[n]] + unit_b2[units[n]];
    log_lik2[n] = normal_lpdf(Y[n] | theta_out2[n], P[n]);
    
    theta_out3[n] = study_a3[studies[n]];
    log_lik3[n] = normal_lpdf(Y[n] | theta_out3[n], P[n]);
    
    theta_out4[n] = unit_b4[units[n]];
    log_lik4[n] = normal_lpdf(Y[n] | theta_out4[n], P[n]);
  }
  
  
  
}

