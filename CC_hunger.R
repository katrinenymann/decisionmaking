#---SCRIPT TO RUN JAGS---

#install.packages("gitcreds")
#library(gitcreds)
#gitcreds_set()

library(R2jags) ; library(polspline) ; library(lmerTest) ; library(ggplot2)
library(tibble) ; library(dplyr) ; library(report)

# data from: https://zenodo.org/record/3764693

seeds <- list(13, 422, 1982, 1997, 2021) #5 different seeds

for (seed in seeds) {
  set.seed(seed)
  print("----------------------------------------")
  print(paste("NOW RUNNING WITH SEED ", seed))
  
  # load data
  a <- read.csv("Experiment_2.csv")
  
  # Preprocess data
  a$UniqueID=factor(paste(a$Experiment, a$Participant.Number))
  a$Previous.Rounds=a$Round-1
  np=subset(a, Punishment.Round=="N")
  
  # extract info about experiment from data
  ntrials <- length(unique(np$Round)) # number of rounds played (aka trials)
  nagents <- length(unique(np$UniqueID)) # number of participants
  vals <- seq(0,20,1) #possible values to contribute - from 0 to 20 tokens #@CHANGED TO 0 NOW
  
  # extract the relevant contribution values from the overall dataframe
  c <- matrix(np$Contribution, nrow = ntrials, ncol = nagents) # converting the contributions to a matrix to feed to JAGS
  c <- t(c) # transposing the matrix so the dimensions fit what we've told JAGS
  Ga <- matrix(np$Group.Contribution/4, nrow = ntrials, ncol = nagents) # specifying a "group average" for each participant
  Ga <- t(Ga)
  #Ga <- colMeans(c) # common group average for *all* participants
  
  # empty data frame to fill with parameter estimates from jags model
  MAP <- c()
  #MAP$omega1 <- array(0,nagents)
  #MAP$lambda <- array(0,nagents)
  #MAP$gamma <- array(0,nagents)
  MAP$pbeta <- array(0,nagents)
  
  #-----------------------------------------------------------

  #prepare jags model for inference
  data <- list("ntrials", "nagents", "vals", "c","Ga") #data inputted into jags
  params <- c("pbeta") #parameters we'll track in jags
  
  # load and run jags model
  samples <- jags(data, inits=NULL, params,
       model.file ="CC_jags.txt",
       n.chains=3, n.iter=5000, n.burnin=1000, n.thin=1)
  
  #SAVE OUTPUT for use in another session
  file_name1 <- paste("jags_", seed, ".Rdata", sep = "")
  save(samples, file = file_name1)
  mcmc13_0 <- as.mcmc(samples)
  file_name2 <- paste("mc_", seed, ".mcmc", sep = "")
  save(mcmc13_0, file = file_name2)
  
  # save maximum a posteriori (MAP) values for parameters from fitted model (see CC_jags.txt for more details)
  for (n in 1:nagents) {
    #X <- samples$BUGSoutput$sims.list$omega1[,n]
    #MAP$omega1[n] <-density(X)$x[which(density(X)$y==max(density(X)$y))]
  
    #X <- samples$BUGSoutput$sims.list$lambda[,n]
    #MAP$lambda[n] <-density(X)$x[which(density(X)$y==max(density(X)$y))]
    
    #X <- samples$BUGSoutput$sims.list$gamma[,n]
    #MAP$gamma[n] <-density(X)$x[which(density(X)$y==max(density(X)$y))]
    
    X <- samples$BUGSoutput$sims.list$pbeta[,n]
    MAP$pbeta[n] <-density(X)$x[which(density(X)$y==max(density(X)$y))] # this is the preference slope
  }
  
  # Add participant info
  MAP$UniqueID = unique(np$UniqueID)
  np_round1 = np %>% subset(Round == 1) %>% select(UniqueID, Condition, Gender, BreakfastToday, BreakfastUsually, HowHungry, Punishment.First, GroupHunger, Condition.Name)
  df = merge(MAP, np_round1, by = "UniqueID")
  
  #save MAP estimates
  file_name3 <- paste("MAPdf_", seed, ".csv", sep = "")
  write.csv(df,file = file_name3, row.names = FALSE) #this is then loaded in analysis script

}