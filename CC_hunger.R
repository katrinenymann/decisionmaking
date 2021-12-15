# seed RNG and load libraries
set.seed(1982)
library(R2jags)
library(polspline)
library(lmerTest)
library(ggplot2)
library(tibble)

# data from here:
# https://zenodo.org/record/3764693

# set working directory and load data
#setwd("/Users/au183362/Documents/asst_prof/teaching/decision_making/E2021/classes/open_data/Fraser_Nettle/")
a <- read.csv("Experiment_2.csv", sep = ";")

# Preprocess data
a$UniqueID=factor(paste(a$Experiment, a$Participant.Number))
a$Previous.Rounds=a$Round-1
np=subset(a, Punishment.Round=="N")

# Manipulation checks
mc.data=subset(a, Round==1 & Punishment.Round=="N")
table(mc.data$BreakfastToday, mc.data$Condition.Name)
t.test(mc.data$HowHungry~mc.data$Condition.Name, var.equal=T)
# Significantly more hungry when assigned to the no breakfast group

# Their analysis
### Model 2 (no punishment game, round 1)####
m2 = lm(Contribution~Condition.Name, data=subset(a, Punishment.Round=="N" & Previous.Rounds==0))
summary(m2)
#### Model 3 (no punishment game, after round 1)#####
m3=lmer(Contribution~Lagged.Contribution+Lagged.MCO*Condition.Name + (1|UniqueGroup/UniqueID), data=subset(a, Punishment.Round=="N"), REML=F)
summary(m3)

# extract info about experiment from data
ntrials <- length(unique(np$Round)) # number of rounds played (aka trials)
nagents <- length(unique(np$UniqueID)) # number of participants
vals <- seq(1,20,1) #possible values to contribute - from 0 to 20 tokens

# extract the relevant contribution values from the overall dataframe
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
     model.file ="CC_jags.txt",
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

MAP <- read.csv("../Andreas/MAP.13.csv")

# Add participant info
MAP$UniqueID = unique(np$UniqueID)
np_round1 = np %>% subset(Round == 1) %>% select(UniqueID, Condition, Gender, BreakfastToday, BreakfastUsually, HowHungry, Punishment.First, GroupHunger, Condition.Name)
df = merge(MAP, np_round1, by = "UniqueID")

model_pbeta = lm(pbeta ~ BreakfastToday + BreakfastUsually + HowHungry + BreakfastToday*BreakfastUsually, df)
summary(model_pbeta)

#Gender effects
t.test(df$pbeta[df$Gender=="F"],df$pbeta[df$Gender=="M"])
boxplot(df$pbeta[df$Gender=="F"],df$pbeta[df$Gender=="M"])

df$contr_avg <- rowMeans(c) # calculating the average contribution per participant

model_pbeta = lm(pbeta ~ BreakfastToday + BreakfastUsually + HowHungry + BreakfastToday*BreakfastUsually, df)
model_omega1 = lm(omega1 ~ BreakfastToday + BreakfastUsually + HowHungry + BreakfastToday*BreakfastUsually, df)
model_lambda = lm(lambda ~ BreakfastToday + BreakfastUsually + HowHungry + BreakfastToday*BreakfastUsually, df)
model_gamma = lm(gamma ~ BreakfastToday + BreakfastUsually + HowHungry + BreakfastToday*BreakfastUsually, df)

summary(model_pbeta)
summary(model_omega1)
summary(model_lambda)
summary(model_gamma)

# plots
mydat <- tibble(x=0:20, 
                y = 0:20)

ggplot(mydat, aes(x=x, y=y)) + 
  geom_segment(df, mapping = aes(x = 0, xend = 20, y = 0, yend = pbeta*20, color = as.factor(Condition.Name))) +
  xlim(0, 20) +
  ylim(0, 20) +
  xlab("Belief about other' contribution") +
  ylab("My contribution") +
  labs(color = "Condition") +
  ggtitle("Slopes for the Conditional Cooperation")
