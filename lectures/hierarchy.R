# code for Lecture 4 hierarchical modeling
library(rstan)

set.seed(12345)

# data
n = 18
y = c(.4,.378,.356,.333,.311,.311,.289,.267,.244,.244,.222,.222,.222,.222,.222,.200,.178,.156)
ybar = mean(y)
s2 = var(y)
sigma2 = 0.06^2

# true value
theta.true <- c(.346,.298,.276,.222,.273,.270,.263,.210,.269,.230,.264,.256,.303,.264,.226,.285,.316,.200)

# define function to draw theta_i from conditional posterior
theta.sample <- function(yi, mu, tau2, sigma2){
    mean.theta <- (sigma2*mu+tau2*yi)/(sigma2+tau2)
    sd.theta <- sqrt(sigma2*tau2/(sigma2+tau2))
    return(rnorm(1,mean.theta,sd.theta))
}

# draw posterior samples
M = 10000
theta.post <- matrix(NA,n,M)
mu.post <- tau2.post <- rep(NA,M)

for(t in 1:M) {
    
    # sample tau2 from a truncated distribution
    tau2 <- 0
    while(tau2 <= 0) tau2 <- 1/rgamma(1,(n-3)/2, rate=(n-1)*s2/2) - sigma2 
    
    # sample mu and theta
    mu <- rnorm(1,ybar,sqrt((sigma2+tau2)/n))
    theta <- sapply(y, function(x) theta.sample(x, mu, tau2, sigma2))
    
    tau2.post[t] <- tau2
    mu.post[t] <- mu
    theta.post[,t] <- theta
}

est <- cbind(theta.true,y,apply(theta.post,1,mean),t(apply(theta.post,1,quantile,probs=c(0.025,0.5,0.975))))
colnames(est) <- c('True','yi','Mean','2.5%','Median','97.5%')
est <- data.frame(est)

round(est,3)


# compare MSE
# individually
table( (est$Median - est$True)^2 < (est$yi - est$True)^2 )
# overall
sum((est$Median - est$True)^2)
sum((est$yi - est$True)^2)



# Stan code ---------------------------------------------------------------
model <- stan_model(
    file = "lectures/hierachy.stan"
    
)


fit <- rstan::sampling(
    object = model,
    data = list(
        n = n,
        # one obs. per group (I think?)
        y = y,
        # number of groups
        K = 2,
        # number of columns
        sigma2 = sigma2
    ),
    chains = 1,
    warmup = 1000,
    iter = 20000,
    thin = 1,
    cores = 4,
    control = list(adapt_delta = 0.95)
)

stan_out <- summary(fit)$summary[1:n + 2, ]
stan_out <- as.data.frame(stan_out)

compare <- data.frame(
    stan_yi = stan_out$mean,
    sample_yi = est$yi
)
