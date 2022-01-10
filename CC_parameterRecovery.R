# seed RNG and load libraries
set.seed(1982)
library(R2jags) # does not work in my R probably needs to be ucloud
library(polspline)

# simulation information
ntrials <- 10
nagents <- 264
niterations <- 10
vals <- seq(0,20,1) #possible values to contribute - from 0 to 20 tokens

# empty data frame to fill with true parameters (randomly generated)
true <- c()
true$omega1 <- array(0,c(nagents,niterations))
true$lambda <- array(0,c(nagents,niterations))
true$gamma <- array(0,c(nagents,niterations))
true$pbeta <- array(0,c(nagents,niterations))

# empty data frame to fill with parameter estimates from jags model
MAP <- c()
MAP$omega1 <- array(0,c(nagents,niterations))
MAP$lambda <- array(0,c(nagents,niterations))
MAP$gamma <- array(0,c(nagents,niterations))
MAP$pbeta <- array(0,c(nagents,niterations))

# load simulation function into workspace
#setwd("C:/Users/au199986/Dropbox/Courses/F20/CognitiveModeling/Module5/CCmodel")
source("Decision-making/CC_fun.R")

# run and fit 1 iterations of the model
for (i in 1:niterations) {

  print(i) #tells you which iteration you're on
  
  #-----------------------------------------------------------
  #set free parameters - we can do inference on these
  parameters <- c()
  #initial beliefs about others average contribution. Set to max as simplification. Could have been free/estimated
  parameters$Gb1 <- c(rep(20, nagents)) 
  #initial weighting of beliefs about others contributionsrelative to own prefs. Higher number means > social influence
  parameters$omega1 <- runif(nagents,.1,1) 
  #decay rate in weighting of beliefs about others - prefs dominate over time following decay function
  parameters$lambda <- runif(nagents,.1,1)
  #parameter weighting of beliefs about what others will contribute, relative to observed contribution (learning rate)
  parameters$gamma <- runif(nagents,.1,1)
  #intercept of linear model relating preferred contributions to possible contribution values
  parameters$p0 <- c(rep(0, nagents)) 
  #slope of linear model relating preferred contributions to possible contribution values
  #this is capped at .7, in parameter recovery. values higher than this leave no room for beliefs because prefs are so high
  #should be free parameter in inference though
  parameters$pbeta <- runif(nagents,.1,.7) 
  #-----------------------------------------------------------
  
  #-----------------------------------------------------------
  # run simulation function and save results in "sims"
  sims <- CC_fun(nagents,ntrials,vals,parameters)
  # Sims has c and omegas the actual contributions and weighting of beliefs about others' contribution relative to prefs
  # aka simulation data
  # Here we can put in our data
  
  # save result of simulation as data
  c <- sims$c # simulated contributions
  Ga <- colMeans(c) #average group member contribution

  # save true parameters inputted into simulation
  true$omega1[,i] <- parameters$omega1 
  true$lambda[,i] <- parameters$lambda
  true$gamma[,i] <- parameters$gamma
  true$pbeta[,i] <- parameters$pbeta
  #-----------------------------------------------------------
    
  #-----------------------------------------------------------
  #prepare jags model for inference
  data <- list("ntrials", "nagents", "vals", "c","Ga") #data inputted into jags
  params <- c("omega1","lambda","gamma","p0","pbeta","c","omega") #parameters we'll track in jags

  # load and run jags model
  samples <- jags(data, inits=NULL, params,
       model.file ="Decision-making/CC_jags_parameter_recovery.txt",
       n.chains=3, n.iter=5000, n.burnin=1000, n.thin=1)
  
  # save maximum a posteriori (MAP) values for parameters from fitted model
  for (n in 1:nagents) {
    
    X <- samples$BUGSoutput$sims.list$omega1[,n]
    MAP$omega1[n,i] <-density(X)$x[which(density(X)$y==max(density(X)$y))]
  
    X <- samples$BUGSoutput$sims.list$lambda[,n]
    MAP$lambda[n,i] <-density(X)$x[which(density(X)$y==max(density(X)$y))]
    
    X <- samples$BUGSoutput$sims.list$gamma[,n]
    MAP$gamma[n,i] <-density(X)$x[which(density(X)$y==max(density(X)$y))]
    
    X <- samples$BUGSoutput$sims.list$pbeta[,n]
    MAP$pbeta[n,i] <-density(X)$x[which(density(X)$y==max(density(X)$y))]
    
  }
  #------------------------------------------------------------------------
  
}

#--------------------------------------------------------
# plot true vs fitted parameters
par(mfrow=c(4,3))

# top row - omega1 for agents 1-3
plot(true$omega1[1,],MAP$omega1[1,])
plot(true$omega1[2,],MAP$omega1[2,])
plot(true$omega1[3,],MAP$omega1[3,])

# 2nd row - lambda for agents 1-3
plot(true$lambda[1,],MAP$lambda[1,])
plot(true$lambda[2,],MAP$lambda[2,])
plot(true$lambda[3,],MAP$lambda[3,])

# 3rd row - gamma for agents 1-3
plot(true$gamma[1,],MAP$gamma[1,])
plot(true$gamma[2,],MAP$gamma[2,])
plot(true$gamma[3,],MAP$gamma[3,])

# bottom row - pbeta for agents 1-3
plot(true$pbeta[1,],MAP$pbeta[1,])
plot(true$pbeta[2,],MAP$pbeta[2,])
plot(true$pbeta[3,],MAP$pbeta[3,])

cor.test(true$pbeta[,],MAP$pbeta[,])

#
# You know the true parameters from the simulated data and then you try to run parameter recovery. 

# We need to fit the jags on the data we want to use



