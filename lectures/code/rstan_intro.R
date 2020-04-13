# example illustrating RStan

data <- read.table(here::here("lectures/beetle_data.txt"), header = F)
colnames(data) <- c('w','y','n')
w <- data$w
y <- data$y
n <- data$n
i <- length(n)

library(rstan)
# options(mc.cores = parallel::detectCores())  

fit1 <- stan(
    file = "lectures/code/beetles.stan",  # Stan program
    data = list(w = w,
                y = y,
                n = n,
                i = i),    # named list of data
    chains = 2,             # number of Markov chains
    warmup = 1000,          # number of warmup iterations per chain
    iter = 2000,            # total number of iterations per chain
    cores = 2,              # number of cores (could use one per chain)
    refresh = 0             # no progress shown
)

summary(fit1)
traceplot(fit1, pars = c("mu", "tau"), inc_warmup = TRUE, nrow = 2)

pairs(fit1, pars = c("mu", "tau", "lp__"), las = 1)

sampler_params <- get_sampler_params(fit1, inc_warmup = TRUE)
summary(do.call(rbind, sampler_params), digits = 2)
lapply(sampler_params, summary, digits = 2)

para_samples <- extract(fit1, pars=c("theta", "mu", "tau", "lp__"),permuted=FALSE, inc_warmup=FALSE)

