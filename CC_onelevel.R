#LOOP
#seeds <- list(13, 422, 1982, 1997, 2021) #5 different ones
seeds <- list(13) #just running with one
for (seed in seeds) {
  set.seed(seed)
  print("----------------------------------------")
  print(paste("NOW RUNNING WITH SEED:", seed))
  
  # seed RNG and load libraries
  library(R2jags)
  library(polspline)
  
  # data from here: https://zenodo.org/record/3764693
  
  # set working directory and load data
  a <- read.csv("Experiment_2.csv")
  
  # preprocess data
  a$UniqueID=factor(paste(a$Experiment, a$Participant.Number))
  a$Previous.Rounds=a$Round-1
  np=subset(a, Punishment.Round=="N") #only look at no punishment
  
  # extract info about experiment from data
  ntrials <- length(unique(a$Round)) # number of rounds played (aka trials)
  nagents <- length(unique(a$Participant.Number))*2 # number of participants #
  vals <- seq(1,20,1) #possible values to contribute - from 0 to 20 tokens
  
  # extract the relevant contribution values from the overall data frame
  c <- matrix(np$Contribution, nrow = ntrials, ncol = nagents) # converting the contributions to a matrix to feed to JAGS
  c <- t(c) # transposing the matrix so the dimensions fit what we've told JAGS
  Ga <- matrix(np$Group.Contribution/4, nrow = ntrials, ncol = nagents) # specifying a "group average" for each participant
  Ga <- t(Ga)
  #Ga <- colMeans(c) # common group average for *all* participants
  
  # empty data frame to fill with parameter estimates from jags model
  MAP <- c()
  MAP$omega1 <- array(0,nagents)
  MAP$lambda <- array(0,nagents)
  MAP$gamma <- array(0,nagents)
  MAP$pbeta <- array(0,nagents)
  #-----------------------------------------------------------
  
  #-----------------------------------------------------------
  #prepare jags model for inference
  data <- list("ntrials", "nagents", "vals", "c","Ga") #data inputted into jags
  params <- c("omega1","lambda","gamma","p0","pbeta","c","omega") #parameters we'll track in jags
  # load and run jags model
  samples <- jags(data, inits=NULL, params,
                  model.file ="CC_jags - Copy_Andreas.txt",
                  n.chains=3, n.iter=5000, n.burnin=1000, n.thin=1)
  
  # save maximum a posteriori (MAP) values for parameters from fitted model (see CC_jags.txt for more details)
  for (n in 1:nagents) {
    
    X <- samples$BUGSoutput$sims.list$omega1[,n]
    MAP$omega1[n] <-density(X)$x[which(density(X)$y==max(density(X)$y))]
    
    X <- samples$BUGSoutput$sims.list$lambda[,n]
    MAP$lambda[n] <-density(X)$x[which(density(X)$y==max(density(X)$y))]
    
    X <- samples$BUGSoutput$sims.list$gamma[,n]
    MAP$gamma[n] <-density(X)$x[which(density(X)$y==max(density(X)$y))]
    
    X <- samples$BUGSoutput$sims.list$pbeta[,n]
    MAP$pbeta[n] <-density(X)$x[which(density(X)$y==max(density(X)$y))] # this is the slope
    
  }
  
  
  #@ida testing (omega vs omega1??? normally distributed???):
  test <- samples$BUGSoutput$sims.list$omega[,2,9] #shape: jags runs, participant, round
  hist(test)
  
  test <- samples$BUGSoutput$sims.list$pbeta[,5] #shape: jags runs, participant, round
  hist(test)
  density(test)$x[which(density(test)$y==max(density(test)$y))]
  
  dim(samples$BUGSoutput$sims.list$gamma)
  
  mcyay <- as.mcmc(samples)
  
  #file_name <- paste("MAP", seed, "csv", sep = ".")
  #write.csv(MAP,file_name, row.names = FALSE)
  
  #maps <- as.data.frame(MAP)
  
  
} #END of seed looping


#---PBETA
#(slope of linear model relating preferred contributions to possible contribution values)

# extracting a gender vector that aligns with the dimensions of MAP with only one entry per participant
# investigating gender effects
idx <- seq(1,length(b$Gender),10) #index looking at one trial per participant (and thereby also allows us to elicit a vector of gender)
#^^i skridt af 10 - tager bare første række for hver deltager ;)
sex <- b$Gender[idx] # extracting a gender vector that aligns with the dimensions of MAP with only one entry per participant
t.test(MAP$pbeta[sex=="F"],MAP$pbeta[sex=="M"])
boxplot(MAP$pbeta[sex=="F"],MAP$pbeta[sex=="M"])

#ida:
condition <- b$Condition[idx]
t.test(MAP$pbeta[condition=="0"],MAP$pbeta[condition=="1"])
boxplot(MAP$pbeta[condition=="0"],MAP$pbeta[condition=="1"])

#AVERAGE CONTRIBUTION
contr_avg <- rowMeans(c) # calculating the average contribution per participant

t.test(contr_avg[sex=="F"],contr_avg[sex=="M"])
boxplot(contr_avg[sex=="F"],contr_avg[sex=="M"])

t.test(contr_avg[condition=="0"],contr_avg[condition=="1"])
boxplot(contr_avg[condition=="0"],contr_avg[condition=="1"])


#værdien kommer ud som EN værdi per subjekt
#kan bygge lmer på hernede med hunger - og interaktioner osv...
#^^no difference in slope - and no difference in contribution!
#make the interactions we want... or e.g. look at hunger!

#session: he doesn't understand the variable
#Trial.Number: first 10 with punishment, then 10 without
#Round: antal runder inden for condition. Det er 10! (så baseret på deres data!)

#Andreas har ikke smidt p0 ind - hvis vi vil have intercept, skal vi selv gøre det

#do.call(rbind.data.frame, your_list)
write.csv(MAP,"MAP_values.csv", row.names = FALSE)