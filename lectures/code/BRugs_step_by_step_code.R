# function to use BRugs package to do Bayesian analysis in R
# step by step way as done in WinBUGS

library(BRugs)

#### STEP 1: load model

# provide the model file name and its path
dir <- 'D:/temp/'
modelfile <- file.path(dir,"linear_model.txt")

# command to load model
modelCheck(modelfile)


### STEP 2: load data
dataList <- list(x = c( 0.0, 0.41, 0.41, 0.41, 0.92, 1.39, 1.61, 1.61, 1.95,
                    2.08, 2.14, 2.20, 2.25, 2.25, 2.30, 2.48, 2.48, 2.56,
                    2.56, 2.67, 2.74, 2.74, 2.80, 2.83, 3.11, 3.37, 3.45),
                y = c(1.80, 1.85, 1.87, 1.77, 2.02, 2.27, 2.15, 2.26, 2.47,
                    2.19, 2.26, 2.40, 2.39, 2.41, 2.50, 2.32, 2.32, 2.43,
                    2.47, 2.56, 2.65, 2.47, 2.64, 2.56, 2.70, 2.72, 2.57), 
                N = 27, 
                beta_mean = rep(0,2), beta_prec = rep(0.01,2)
                )

# use the Brugs function 'bugsData' to format data and write into a txt file
datafile <- file.path(dir,'data.txt')
bugsData(dataList,fileName=datafile)

# command to load data 
modelData(datafile)

# Otherwise, you can directly load data without saving to a datafile
modelData(bugsData(dataList))


#### STEP 3: Compile model with 2 chians
modelCompile(numChain=2)


#### STEP 4: load initial values
initList <- list(list(beta = c(1,1), 
                      sigma = 1),
                 list(beta = c(0,0),
                      sigma = 0.5)
                )

# use the Brugs function 'bugsInits' to format initial lists and write into a txt file
initfile1 <- file.path(dir,'init1.txt')
initfile2 <- file.path(dir,'init2.txt')
bugsInits(initList,numChains=2,fileName=c(initfile1,initfile2))

# command to load initial sets one by one
modelInits(initfile1)
modelInits(initfile2)

# option 2: load two sets of initials at one time
modelInits(c(initfile1, initfile2)) 
# option 3: let the program generate initials 
modelGenInits()



#### STEP 5: burn in step
modelUpdate(1000)


#### STEP 6: set parameters of which posterior samples will be collected
samplesSet(c("beta","sigma"))


#### STEP 7: more iterations to collect samples
modelUpdate(10000)


#### STEP 8: check results
samplesStats("*")                   # summary statistics
samplesHistory("*",mfrow=c(3,1))    # traceplot

samplesDensity("*",mfrow=c(3,1))    # plot the densities
samplesBgr("beta",mfrow=c(2,1))     # plot bgr convergence diagonostic statistics
samplesAutoC("beta",1,mfrow=c(1,2)) # plot autocorrelations of 1st chain


#### DIC calculation
dicSet()                # start DIC monitor
modelUpdate(10000)      # run more iterations for DIC calculation
dicStats()              # display DIC statistics


#### Use 'samplesSample' function to extract posterior samples, save for future analysis
beta0.samples <- samplesSample("beta[1]")
beta1.samples <- samplesSample("beta[2]")
samplefile <- file.path(dir,'post_samples.txt')
write.table( cbind(beta0.samples ,beta1.samples), file=samplefile, row.names=FALSE)
