###########################################################################
#####  Now assume alpha is unknown, and assign a hyperprior for alpha #####
#####  h(alpha) = Exp(mu), with mu=1                                  #####
#####  Use MH-Gibbs hybrid algorithm to conduct posterior sampling    #####
#####  Use N(alpha^t-1,tau^2) as the proposal density for the MH step #####
###########################################################################


rm(list=ls())
# set.seed(12345)

library(MCMCpack) # this contains the function to generate inv-gamma random variable

### First, read data from txt file
data <- read.table(here::here('lectures/pump_data.txt'),header=T)
#data <- read.table('D:/4. Teaching/UM_SPH7440/Notes_Sp2016/L5.Sampling_in_R/data/pump_data.txt',header=T)

k <- dim(data)[1]
y <- data$Y
t <- data$t


### Specify priors and MH proposal stepsize
c <- 0.1
d <- 1
mu <- 1
tau <- 0.01


### define the function to calculate log (unnormalized) full conditional for alpha
log.full.a <- function(a,theta,beta){
    (exp(a)-1)*sum(log(theta)) - k*log(gamma(exp(a))) - k*exp(a)*log(beta) - exp(a)/mu + a
}


### Specify MCMC runs, burnin numbers
runs <- 10000; 
burn <- 1000; 


### Specify initial values 
theta.init <- rep(1,k)
beta.init <- 1
alpha.init <- 1


### Create arrays to store MCMC samples
theta.save <- array(NA,c(runs,10))
beta.save <- array(NA,runs)
alpha.save <- array(NA,runs)

accept.flag <- array(0)  # record movements of alpha chains


### Now conduct Gibbs sampling
theta <- theta.init
beta <- beta.init
alpha <- alpha.init
a <- log(alpha)

for(iter in 1:(runs+burn)){ 
    
    # at each iteration t
    
    # Step 1: generate each theta_i from Gamma full conditional
    for(i in 1:k) {
        theta[i] <- rgamma(1,y[i]+alpha,rate=t[i]+1/beta)
    }
    
    # Step 2: generate beta from Inverse-Gamma full condtional
    beta <- rinvgamma(1,k*alpha+c,scale=sum(theta)+1/d)
    
    # Step 3: use MH algorithm to sample a from its full conditional
    # step 3.1: generate a proposal value of a conditional on its current value
    a.prop <- rnorm(1,a,tau)
    
    # step 3.2: calculate log ratio for the unnormalized full conditional distribution
    log.ratio <- log.full.a(a.prop,theta,beta) - log.full.a(a,theta,beta)
    
    # step 3.3: accept/reject the candidate value of a
    if(log(runif(1)) < log.ratio) {
        a <- a.prop
        accept.flag <- accept.flag+1
        alpha <- exp(a)
    }
    
    
    # save the current value of theta after burnin interations
    if(iter>burn) {
        theta.save[iter-burn,] <- theta
        beta.save[iter-burn] <- beta
        alpha.save[iter-burn] <- alpha
    }
}


### Check convergence
# acceptance rate
accept.flag/(runs+burn)

# trace plot
par(mfrow=c(3,1),mar=c(4,4.5,1,0.5))
plot(theta.save[,1],type='l',xlab='iteration',ylab=expression(theta[1]))
plot(beta.save,type='l',xlab='iteration',ylab=expression(beta))
plot(alpha.save,type='l',xlab='iteration',ylab=expression(alpha))

summary(theta.save)
summary(beta.save)
summary(alpha.save)



