# This function is to use importance sampling algorithm to approximate E(theta), 
# from the Bayesian model
# Likelihood:  X1,...,Xn ~ N(0,theta)
# Prior: theta ~ Gamma(3,rate=0.5)
# The importance function g(theta) = Gamma(a=c*s^2,b=c), where s^2 = Var(X)


set.seed(12345)


# first define the function to conduct importance sampling estimation
# given the data X, the constant value c, and the number of samples to generate N


# prep1: define the function to calculate log(p*(theta|x)) for each value of theta
log.posterior <- function(theta_j, X){ 
    n <- length(X)
    ( (4-n)/2 )*log(theta_j) - (1/(2*theta_j)) * (theta_j^2 + sum(X^2))
}
 
# prep2: define the function to calculate log(weights)=log(p*(theta|x))-log(g(theta))
log.weight <- function(theta_j, X, a, b){
    log.posterior(theta_j, X) - dgamma(theta_j, a,rate=b,log=TRUE)
}


IS <- function(X, c, N) {
 
    # first calculate all the parameters 
    n <- length(X) 
    s2 <- var(X)
    a <- c*s2
    b <- c
    
    # step 1: generate N samples from the importance function g(theta)
    theta <- rgamma(N,a,rate=b)
    
    
    # step 2: use the defined functions to calculate log of weights
    # We work in the log scale for the sake of computational stability
    lw <- log.weight(theta, X, a, b)
    
    
    # The following step is to avoid inifite values when taking exponential of lw
    # Using the step we have weights:  w = exp(lw) = exp(w.scaled)*exp(max(lw))
    # Since w are in both the nominator and the denominator, exp(max(lw)) will be canceled out
    lw.scaled <- lw - max(lw)
 

    ## Finally, calculate the importance sampling estimate
    theta.hat <- mean(theta * exp(lw.scaled)) / mean(exp(lw.scaled))

        
    ## Now, evaluate the importance function using the coefficient of variance
    c.v <- sd(exp(lw.scaled))/mean(exp(lw.scaled))
    
    
    return(c(theta.hat,c.v))
}



#### Apply: Generate 100 random variables from N(0,2^2)

X <- rnorm(100, mean = 0, sd = 2)

## consider a grid of c values
c.candidate <- seq(0.05, 20, length=500)

## specify vectors to store theta.hat and c.v for each value of c.candidate
theta.hats <- rep(NA, 500)
cvs <- rep(NA, 500)

## calculate estimate theta.hat and coef of var for each value of c.candidate

for(i in 1:500) {

    OUT <- IS(X = X, c = c.candidate[i], 1000) 
    theta.hats[i] <- OUT[1]
    cvs[i] <- OUT[2]
} 


## Examine where coef. of variance is low

plot(c.candidate, cvs, ylab="coefficient of variance", xlab="c", col=2, type="l")


## target mean estimator

plot(c.candidate, theta.hats, xlab="c", ylab="E(X) estimate", col=4, type="l")


## final estimate

IS(X, c.candidate[which.min(cvs)], 1000)[1]

