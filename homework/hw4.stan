
data {
  // Define data in this block
  
  //data for bike route intersections
  int<lower=0> n_bike;// sample size for bike route intersections
  int<lower=0> y_bike1[n_bike]; // bike route data model 1
  int<lower=0> y_bike2[n_bike]; //bike route data model 2
  int<lower=0> y_bike3[n_bike]; //bike route data model 3
  
    //data for non-bike route intersections
  int<lower=0> n_nonbike;// sample size for streets intersections
  int<lower=0> y_nonbike1[n_nonbike]; // nonbike routes data1
  int<lower=0> y_nonbike2[n_nonbike]; // nonbike routes data2
  int<lower=0> y_nonbike3[n_nonbike]; // nonbike routes data3
  
}


parameters {
  // Define parameters in this block
  
  // M1 parameters
  // parameters for bike route intersection
  real<lower = 0> lambda_bike1[n_bike]; //lambda for streets
  real<lower = 0> alpha_bike1; // alpha for streets
  real<lower = 0> beta_bike1; // beta for streets
  
    // parameters for non-bike route intersection
  real<lower = 0> lambda_nonbike1[n_nonbike]; //lambda for streets
  real<lower = 0> alpha_nonbike1; // alpha for streets
  real<lower = 0> beta_nonbike1; // beta for streets
  
  // M2 parameters
  real<lower = 0> lambda_bike2;
  real<lower = 0> lambda_nonbike2;
  
  // M3 parameters
  real<lower = 0> lambda_bike3[n_bike];
  real<lower = 0> lambda_nonbike3[n_nonbike];
}

// The model
model {
  // Define model in this block
  // M1
  alpha_bike1 ~ gamma(0.01, 0.01);
  beta_bike1 ~ gamma(0.01, 0.01);
  alpha_nonbike1 ~ gamma(0.01, 0.01);
  beta_nonbike1 ~ gamma(0.01, 0.01);
  
  for (i in 1:n_bike) {
    lambda_bike1[i] ~ gamma(alpha_bike1, beta_bike1);
    y_bike1[i] ~ poisson(lambda_bike1[i]);
  }
  
  for (j in 1:n_nonbike) {
    lambda_nonbike1[j] ~ gamma(alpha_nonbike1, beta_nonbike1);
    y_nonbike1[j] ~ poisson(lambda_nonbike1[j]);
  }
  
  // M2
  lambda_bike2 ~ gamma(0.01, 0.01);
  for (i in 1:n_bike) {
    y_bike2[i] ~ poisson(lambda_bike2);
  }
  
  lambda_nonbike2 ~ gamma(0.01, 0.01);
  for (j in 1:n_nonbike) {
    y_nonbike2[j] ~ poisson(lambda_nonbike2);
  }
  
  // M3 
  for (i in 1:n_bike) {
    lambda_bike3[i] ~ gamma(0.01, 0.01);
    y_bike3[i] ~ poisson(lambda_bike3[i]);
  }
  
  for (j in 1:n_nonbike) {
    lambda_nonbike3[j] ~ gamma(0.01, 0.01);
    y_nonbike3[j] ~ poisson(lambda_nonbike3[j]);
  }
}

generated quantities{
  // Define other generated quantities in this block
  
  // Difference in lambdas
 real lambda_rep_bike[n_bike];
 real y_rep_bike[n_bike];
 real lambda_rep_nonbike[n_nonbike];
 real y_rep_nonbike[n_nonbike];
 
 // log likelihoods for WAIC
 vector[n_bike] log_lik1_bikes;
 vector[n_nonbike] log_lik1_nonbikes;
 vector[n_bike] log_lik2_bikes;
 vector[n_nonbike] log_lik2_nonbikes;
 vector[n_bike] log_lik3_bikes;
 vector[n_nonbike] log_lik3_nonbikes;
 vector[n_bike + n_nonbike] log_lik1;
 vector[n_bike + n_nonbike] log_lik2;
 vector[n_bike + n_nonbike] log_lik3;
 
 for (i in 1:n_bike) {
   lambda_rep_bike[i] = gamma_rng(alpha_bike1, beta_bike1);
   y_rep_bike[i] = poisson_rng(lambda_rep_bike[i]);
   log_lik1_bikes[i] = poisson_lpmf(y_bike1[i] | lambda_bike1[i]);
   log_lik2_bikes[i] = poisson_lpmf(y_bike2[i] | lambda_bike2);
   log_lik3_bikes[i] = poisson_lpmf(y_bike3[i] | lambda_bike3[i]);
 }

 for (j in 1:n_nonbike) {
   lambda_rep_nonbike[j] = gamma_rng(alpha_nonbike1, beta_nonbike1);
   y_rep_nonbike[j] = poisson_rng(lambda_rep_nonbike[j]);
   log_lik1_nonbikes[j] = poisson_lpmf(y_nonbike1[j] | lambda_nonbike1[j]);
   log_lik2_nonbikes[j] = poisson_lpmf(y_nonbike2[j] | lambda_nonbike2);
   log_lik3_nonbikes[j] = poisson_lpmf(y_nonbike3[j] | lambda_nonbike3[j]);
   
 }

// Append the log likelihood vectors together
log_lik1 = append_row(log_lik1_bikes, log_lik1_nonbikes);
log_lik2 = append_row(log_lik2_bikes, log_lik2_nonbikes);
log_lik3 = append_row(log_lik3_bikes, log_lik3_nonbikes);


 
} // end generated quantities


