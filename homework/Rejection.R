# This function is to use rejection sampling algorithm to 
# sample from a unnormalized density p*(theta) = exp(-theta^2/2) (this is an unnormalized gaussian density, pretend we don't know that)
# The envelope function to use is g(theta) = Cauchy(0,gamma)
# Choose cauchy because it is heavy tailed - should cover a normal distribution



set.seed(12345)

## step 0: first plot p*(theta) and M*g(theta) for different M values
##         choose an appropriate value of M

# function to calculate log of unnormalized posterior
log.posterior <- function(theta_j){
    -theta_j^2/2
}

# calcuate p*(theta|x) and M*g(theta) for a grid of theta values
theta <- seq(-5,5,by=0.01)
p.star <- exp(log.posterior(theta))
plot(theta,p.star,type='l',ylim=c(0,1.5))

# consider different M values
M.cand <- c(1,6,6.3,7)
for(i in 1:length(M.cand)){ 
    lines(theta,dcauchy(theta, location=0, scale=2)*M.cand[i],col=i+1,lty=2)
}
legend('topright',legend=paste0('M=',as.character(M.cand)),col=2:5,lty=2)

# Based on the plot, we choose M=6.3, and proceed to rejection sampling



## Now define function to perform rejection sampling
RS <- function(M,N) {
    
    # specify half cauchy parameters
    gamma <- 2    
    
    theta.samples <- NULL
    
    # use while loop to collect samples interatively until N samples are collected
    while(length(theta.samples)<N) {
        
        theta.star <- rcauchy(1,0,gamma)
        u <- runif(1)
        log.ratio <- log.posterior(theta.star)-log(M)-dcauchy(theta.star,0,gamma,log=TRUE)
    
        if(log(u)<log.ratio) theta.samples <- c(theta.samples,theta.star)
    }
    
    return(theta.samples)
}


## Now collect samples using the defined function
M <- 6.3
theta.samples <- RS(M, 10000)


## check the empirical distribution of the samples and compare to the true density N(0,1)

hist(theta.samples, freq=FALSE, ylim=c(0,0.5), breaks=100, main="N(0,1) by rejection sampling",xlab=expression(theta))

lines(theta, dnorm(theta),col=2)


