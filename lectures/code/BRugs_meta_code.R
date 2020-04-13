
# Step0: save model file, data, and initial value to .RData file to repetitive uses
rm(list=ls())

# path to model file
dir <- 'D:/temp/'
modelfile <- file.path(dir,"linear_model.txt")

# specify data
x <- c( 0.0, 0.41, 0.41, 0.41, 0.92, 1.39, 1.61, 1.61, 1.95,
       2.08, 2.14, 2.20, 2.25, 2.25, 2.30, 2.48, 2.48, 2.56,
       2.56, 2.67, 2.74, 2.74, 2.80, 2.83, 3.11, 3.37, 3.45)
y <- c(1.80, 1.85, 1.87, 1.77, 2.02, 2.27, 2.15, 2.26, 2.47,
      2.19, 2.26, 2.40, 2.39, 2.41, 2.50, 2.32, 2.32, 2.43,
      2.47, 2.56, 2.65, 2.47, 2.64, 2.56, 2.70, 2.72, 2.57)
N <- length(x)
beta_mean = rep(0,2)
beta_prec = rep(0.01,2)

# specify initial values
initList <- list(list(beta = c(1,1), 
                      sigma = 1),
                 list(beta = c(0,0),
                      sigma = 0.5)
)
save.image(file='D:/temp/Lin_data.RData')

#############################################################################

#### Use BRugs meta function 'BRugsFit'
rm(list=ls())
library(BRugs)

load('D:/temp/Lin_data.RData')
ls()

# conduct Bayesian analysis in one step
BRugs.res <- BRugsFit(modelFile=modelfile, 
                    data=list(x=x,y=y,N=N,beta_mean=beta_mean,beta_prec=beta_prec), inits=initList, 
                    numChains=2, parametersToSave = c("beta","sigma"),
                    nBurnin=1000, nIter=10000, nThin=1,
                    DIC=TRUE   
                    )
samplesStats("*")  
samplesHistory("*")


#############################################################################

#### Set the argument 'coda=TRUE', the function will return posterior samples in coda format
codaSamples <- BRugsFit(modelFile=modelfile, 
                        data=list(x=x,y=y,N=N,beta_mean=beta_mean,beta_prec=beta_prec), inits=initList, 
                        numChains=2, parametersToSave = c("beta","sigma"),
                        nBurnin=1000, nIter=10000, nThin=1,
                        coda=TRUE
                        )

# the resulting 'codaSamples' has these indices: 
# codaSamples[[chainIdx]][iteratIdx, paramIdx]
codaSamples[[1]][1:10,'beta[1]']


# Extract posterior samples of certain parameter for future use
post.samples <- as.matrix(codaSamples)
beta0.samples <- post.samples[,"beta[1]"]
beta1.samples <- post.samples[,"beta[2]"]
samplefile <- file.path(dir,'post_samples.txt')
write.table( cbind(beta0.samples ,beta1.samples), file=samplefile, row.names=FALSE)


################################################################################

#### Easy for mcmc chain diagnosis using 'CODA' package for the returned coda object

library('coda')

plot(codaSamples, trace=TRUE, density=TRUE)
gelman.plot(codaSamples)
autocorr.plot(codaSamples)
crosscorr.plot(codaSamples)


# unload package
detach('package:BRugs')








