# function to write and save BUGS model in .txt file

rm(list=ls())
library(BRugs)


# Option 1: Use the 'WriteModel' function implemented in BRugs
lin_model <- function(){
    
    for(i in 1:N) {
        y[i] ~ dnorm(mu[i],y_prec)
        mu[i] <- beta[1]+beta[2]*x[i]
    }
    
    beta[1] ~ dnorm(beta_mean[1],beta_prec[1])
    beta[2] ~ dnorm(beta_mean[2],beta_prec[2])
    
    y_prec <- 1/(sigma*sigma)
    sigma ~ dunif(0.01,100)
    
}

dir <- 'D:/temp/'
modelfile <- file.path(dir,"linear_model.txt")
writeModel(lin_model,modelfile)


# Option 2: Use the base R command 'writeLines'
modelstring <- "
model {

    for(i in 1:N) {
        y[i] ~ dnorm(mu[i],y_prec)
        mu[i] <- beta[1]+beta[2]*x[i]
    }
    
    beta[1] ~ dnorm(beta_mean[1],beta_prec[1])
    beta[2] ~ dnorm(beta_mean[2],beta_prec[2])
    
    y_prec <- 1/(sigma*sigma)
    sigma ~ dunif(0.01,100)

}
"
dir <- 'D:/temp/'
modelfile <- file.path(dir,"linear_model.txt")
writeLines(modelstring,con=modelfile)


# Option 3: Directly write and save your model in a text editor
