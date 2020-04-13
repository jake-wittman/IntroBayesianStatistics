data {
    int<lower=1> N1; //number of data points
    vector[N1] Y1; //number of outcomes

    int<lower=1> N0; // number of new points
    matrix[N1+N0,N1+N0] dist; //distances between points
    
    int<lower=1> p;  // number of predictors
    matrix[N1+N0,p] X;      // matrix of predictors
}
transformed data{
    int<lower=1> N = N1 + N0;
}
parameters{
    vector[N0] Y0;
    vector[p] beta;
    real<lower=0> sigma_sq;
    real<lower=0> tau_sq;
    real<lower=0> phi;
}
transformed parameters{
    vector[N1+N0] mu;
    for(i in 1:N) mu[i] = dot_product(X[i,],beta);
}
model{
    vector[N] Y;
    matrix[N,N] Sigma;
    matrix[N,N] L;
    for(i in 1:N1) Y[i] = Y1[i];
    for(i in 1:N0) Y[N1+i] = Y0[i];
    for(i in 1:(N-1)){
        for(j in (i+1):N){
            Sigma[i,j] = sigma_sq*exp((-1)*phi*dist[i,j]);
            Sigma[j,i] = Sigma[i,j];
        }
    }
    for(i in 1:N) Sigma[i,i] = sigma_sq+tau_sq;
    L = cholesky_decompose(Sigma);
    Y ~ multi_normal_cholesky(mu,L);
    
    beta ~ normal(0,100);
    sigma_sq ~ inv_gamma(0.6, 0.6);
    tau_sq ~ inv_gamma(0.1, 0.1);
    phi ~ gamma(0.1, 0.1);
    
}


