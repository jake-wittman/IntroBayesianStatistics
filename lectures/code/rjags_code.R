#### Use rjags for Bayesian Analysis

rm(list=ls())
library(rjags) # rjags automatically load the 'coda' package for MCMC diagnosis


load('D:/temp/Lin_data.RData')
ls()

# First create, initialize, and adapt the model:
jagsModel <- jags.model(file=modelfile, 
                        data=list(x=x,y=y,N=N,beta_mean=beta_mean,beta_prec=beta_prec), inits=initList, 
                        n.chains = 2, n.adapt = 1000
)

# Use the function 'update' for burnin iterations
update(jagsModel, n.iter=1000)

# Now run MCMC and collect posterior samples in coda format
codaSamples <- coda.samples(model=jagsModel, variable.names = c("beta","sigma"), 
                            n.iter = 10000, thin = 1)


# DIC calculation
dic.samples(model=jagsModel,n.iter=10000,thin=1,type='pD')


# unload package
detach('package:rjags')

