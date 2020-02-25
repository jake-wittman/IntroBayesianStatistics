# This file is to use Metropolis-Hastings method to obtain posterior samples
# for the beetle data in Example 3.7

set.seed(12345)

### First, read data from txt file
data <- read.table('C:/temp/beetle_data.txt',header=F)
colnames(data) <- c('w','y','n')
w <- data$w
y <- data$y
n <- data$n


### Specify prior parameters
a0 <- .25
b0 <- 4
c0 <- 2
d0 <- 10
e0 <- 2
f0 <- 1000


### Define functions to calculate log(g(w)) and log(unormalized posterior) respectively
log.g <- function(theta){
    x <- (w-theta[1])/exp(theta[2])
    return(exp(theta[3])*(x-log(1+exp(x))))
}

log.posterior <- function(theta) {
    sum(y*log.g(theta)+(n-y)*log(1-exp(log.g(theta)))) +
        a0*theta[3]-2*e0*theta[2] - 0.5*(theta[1]-c0)^2/d0^2 - exp(theta[3])/b0 - exp(-2*theta[2])/f0
}


### Specify covariance matrix for proposal density q()
Sigma <- diag(c(0.00012,0.033,0.10))


### Specify MCMC runs, burnin numbers, # of chains 
runs <- 10000; 
burn <- 0; 
chains<-3;

### Specify initial values for each chain
theta.init <- matrix(c(2.1,1.8,1.5, -2.5,-3.5,-4, 0 ,-1,-2), chains,3)


### Create [runs X 3 X  chains]-dimensional arrays to store MCMC samples
theta.save <- array(NA,c(runs,3,chains))

### This is the matrix to record the number of iterations at which
### the proposed new values are accepted for each chain
accept.flag <- rep(0,chains)


### Now conduct MH sampling for each chain sequentially
for(ch in 1:chains){
    
    # give the intial values as the current values of theta for this chain
    theta <- theta.init[ch,]
    
    # collect samples for this chain
    for(t in 1:(runs+burn)){ 
        
        # at each iteration t
        
        # Step 1: generate a candidate value using the proposal distribution centered at current theta
        theta.prop <- rnorm(3)%*%chol(Sigma)+theta;
        
        # Step 2: calculate log of ratio
        log.ratio <- log.posterior(theta.prop) - log.posterior(theta)
            
        # Step 3: accept the proposed valued with probability min(ratio,1)
        if(log(runif(1)) < log.ratio) {
            theta <- theta.prop
            # add one to the count of acceptance
            accept.flag[ch] <- accept.flag[ch]+1
        }
        
        # save the current value of theta after burnin interations
        if(t>burn) theta.save[t-burn,,ch] <- theta
    }
}



### Convergence diagnosis

# calculate the acceptance rate
sum(accept.flag)/runs/3  #0.1384

# Trace plots for the three parameters
par(mfrow=c(3,1),mar=c(4,4.5,1,0.5))
matplot(theta.save[,1,],type='l',ylim=c(1.5,2.1),xlab='iteration',ylab=expression(mu))
matplot(theta.save[,2,],type='l',ylim=c(-4.5,-1.5),xlab='iteration',ylab=expression(log(sigma)))
matplot(theta.save[,3,],type='l',ylim=c(-2,3.5),xlab='iteration',ylab=expression(log(m[1])))


# Autocorrelation plots for chain 1 samples after discarding the first 1000 samples
acf(theta.save[1001:10000,1,1],lag.max=9000)
acf(theta.save[1001:10000,2,1],lag.max=9000)
acf(theta.save[1001:10000,3,1],lag.max=9000)


# Bivariate plots to check cross-correlation
par(mfrow=c(2,2),mar=c(4,4.5,1,0.5))
plot(theta.save[1001:10000,1,1], theta.save[1001:10000,2,1], xlab=expression(theta[1]),ylab=expression(theta[2]))
plot(theta.save[1001:10000,2,1], theta.save[1001:10000,3,1], xlab=expression(theta[2]),ylab=expression(theta[3]))
plot(theta.save[1001:10000,1,1], theta.save[1001:10000,3,1], xlab=expression(theta[1]),ylab=expression(theta[3]))


# Calculate the sample cross-correlations using all chians
theta1 <- c(theta.save[1001:10000,1,])
theta2 <- c(theta.save[1001:10000,2,])
theta3 <- c(theta.save[1001:10000,3,])
cor(cbind(theta1,theta2,theta3))
cov(cbind(theta1,theta2,theta3))

# Using the samples to calculate mean for mu and m1
mean(theta1) #1.81
mean(exp(theta3)) #0.39

#######################################################
### Now it's your work to try the new proposal      ###
### N(0,Sigma.new), where                           ###
### Sigma.new = 2*cov(cbind(theta1,theta2,theta3))  ###
#######################################################

