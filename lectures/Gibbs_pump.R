# This file is to use Gibbs sampling method to obtain posterior samples
# for the pump data 

rm(list=ls())
set.seed(12345)

library(invgamma) # this contains the function to generate inv-gamma random variable

### First, read data from txt file
data <- read.table(here::here('lectures/pump_data.txt'),header=T)

k <- dim(data)[1]
y <- data$Y
t <- data$t


### Specify priors
c <- 0.1
d <- 1
alpha <- 1


### Specify MCMC runs, burnin numbers
runs <- 10000; 
burn <- 1000; 


### Specify initial values 
theta.init <- rep(1,k)
beta.init <- 1


### Create arrays to store MCMC samples
theta.save <- array(NA,c(runs,10))
beta.save <- array(NA,runs)


### Now conduct Gibbs sampling
theta <- theta.init
beta <- beta.init

for(iter in 1:(runs+burn)){ 
    
    # at each iteration t
    
    # Step 1: generate each theta_i from Gamma full conditional
    for(i in 1:k) {
        theta[i] <- rgamma(1,y[i]+alpha,rate=t[i]+1/beta)
    }
    
    # Step 2: generate beta from Inverse-Gamma full condtional
    beta <- rinvgamma(1,k*alpha+c,scale=sum(theta)+1/d)
    
    # save the current value of theta after burnin interations
    if(iter>burn) {
        theta.save[iter-burn,] <- theta
        beta.save[iter-burn] <- beta
    }
}

### Check convergence
par(mfrow=c(3,1),mar=c(4,4.5,1,0.5))
plot(theta.save[,1],type='l',xlab='iteration',ylab=expression(theta[1]))
plot(theta.save[,1],type='l',xlab='iteration',ylab=expression(theta[2]))
plot(beta.save,type='l',xlab='iteration',ylab=expression(beta))

summary(theta.save)
summary(beta.save)



