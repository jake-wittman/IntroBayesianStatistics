#### Use R2WinBUGS for Bayesian Analysis

rm(list=ls())
library(R2WinBUGS)

load('D:/temp/Lin_data.RData')
ls()

# conduct Bayesian analysis in one step, # set 'codaPkg=TRUE' to produe coda object
bugs.res <- bugs(data=list(x=x,y=y,N=N,beta_mean=beta_mean,beta_prec=beta_prec), inits=initList, 
                parameters.to.save = c("beta","sigma"), model.file=modelfile,
                n.chains=2, n.iter=10000, n.burnin=1000, n.thin=1,
                DIC = TRUE, # codaPkg=TRUE,
                bugs.directory='C:/Program Files/WinBUGS14/'
                )

print(bugs.res)
plot(bugs.res)


# unload package
detach('package:R2WinBUGS')
