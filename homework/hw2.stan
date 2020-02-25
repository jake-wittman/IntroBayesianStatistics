
data {
  // Define data in this block
  
  //data for bike route intersections
  int<lower=0> n_s;// sample size for bike route intersections
  int<lower=0> y_s[n_s]; // bike route data
  
    //data for non-bike route intersections
  int<lower=0> n;// sample size for streets intersections
  int<lower=0> y[n]; // streets data
  
}


parameters {
  // Define parameters in this block
  
  // parameters for bike route intersection
  real<lower = 0> lambda_s[n_s]; //lambda for streets
  real<lower = 0> alpha_s; // alpha for streets
  real<lower = 0> beta_s; // beta for streets
  
    // parameters for non-bike route intersection
  real<lower = 0> lambda[n]; //lambda for streets
  real<lower = 0> alpha; // alpha for streets
  real<lower = 0> beta; // beta for streets
}

// The model
model {
  // Define model in this block
  
  // Model for bike route intersections
  target += gamma_lpdf(alpha_s | 0.01, 0.01); // gamma hyperprior on alpha 
  target += gamma_lpdf(beta_s | 0.01, 0.01); // gamma hyperprior on beta
  target += gamma_lpdf(lambda_s | alpha_s, beta_s); // gamma prior on lambda
  target += poisson_lpmf(y_s | lambda_s); // distribution of response given lambda
  
    // Model for non-bike route intersections
  target += gamma_lpdf(alpha | 0.01, 0.01); // gamma hyperprior on alpha
  target += gamma_lpdf(beta | 0.01, 0.01); // gamma hyperprior on beta
  target += gamma_lpdf(lambda | alpha, beta); // gamma prior on lambda
  target += poisson_lpmf(y | lambda); // distribution of response given lambda
  
}

generated quantities{
  // Define other generated quantities in this block
  
  // Difference in lambdas
 real posterior_mean_difference;
 real pred_lambda;
 real pred_y;
 
  
 posterior_mean_difference = gamma_rng(alpha_s, beta_s) - gamma_rng(alpha, beta);
 
 pred_lambda = gamma_rng(alpha_s, beta_s);
 pred_y = poisson_rng(pred_lambda);

  
}


