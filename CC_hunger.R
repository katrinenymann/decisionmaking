# set seed and load libraries
#Ida is back #this time actually me #hmm this time?
set.seed(13)

install.packages("gitcreds")
library(gitcreds)
gitcreds_set()

library(R2jags)
library(polspline)
library(lmerTest)
library(ggplot2)
library(tibble)
#install.packages("devtools")
#devtools::install_github("easystats/report")
library(report)

# data from here:
# https://zenodo.org/record/3764693

# set working directory and load data
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

samples = jags_output #load jags output from earlier run

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

p_load(tidyverse)

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

report(model_pbeta)

# plots
mydat <- tibble(x=0:20, 
                y = 0:20)

ggplot(mydat, aes(x=x, y=y)) + 
  geom_segment(df, mapping = aes(x = 0, xend = 20, y = 0, yend = pbeta*20, color = as.factor(BreakfastToday))) +
  xlim(0, 20) +
  ylim(0, 20) +
  xlab("Belief about other' contribution") +
  ylab("My contribution") +
  labs(color = "Condition") +
  ggtitle("Slopes for the Conditional Cooperation")

df_nb = df %>% subset(BreakfastToday == 0)
df_b = df %>% subset(BreakfastToday == 1)

ggplot(mydat, aes(x=x, y=y)) + 
  geom_segment(df_nb, mapping = aes(x = 0, xend = 20, y = 0, yend = mean(pbeta*20), color = as.factor(BreakfastToday))) +
  geom_segment(df_b, mapping = aes(x = 0, xend = 20, y = 0, yend = mean(pbeta*20), color = as.factor(BreakfastToday))) +
  xlim(0, 20) +
  ylim(0, 20) +
  xlab("Belief about other' contribution") +
  ylab("My contribution") +
  labs(color = "Condition") +
  ggtitle("Slopes for the Conditional Cooperation")

df_nb = df %>% subset(Condition == 0)
df_b = df %>% subset(Condition == 1)

ggplot(mydat, aes(x=x, y=y)) + 
  geom_segment(df_nb, mapping = aes(x = 0, xend = 20, y = 0, yend = mean(pbeta*20), color = as.factor(Condition))) +
  geom_segment(df_b, mapping = aes(x = 0, xend = 20, y = 0, yend = mean(pbeta*20), color = as.factor(Condition))) +
  xlim(0, 20) +
  ylim(0, 20) +
  xlab("Belief about other' contribution") +
  ylab("My contribution") +
  labs(color = "Condition") +
  ggtitle("Slopes for the Conditional Cooperation")

# Lambda: Are one's decay rate determined by breakfast eating? 
# In the paper they find the no breakfast group to be more influenced by the other group members 
# aka weigting you belief over preference. Aka smaller decay rate
model_lambda = lm(lambda ~ BreakfastToday + BreakfastUsually + HowHungry + BreakfastToday*BreakfastUsually, df)


df_nb = df %>% subset(BreakfastToday == 0)
df_b = df %>% subset(BreakfastToday == 1)

# Boxplot
ggplot(df, aes(x=as.factor(BreakfastToday), y=lambda)) + 
  geom_boxplot() +
  geom_point() + 
  xlab("Condition") +
  ylab("Value") +
  ggtitle("Initial weighting of beliefs per condition")

### 
# Omega1 initial weigting of beliefs

model_omega1 = lm(omega1 ~ BreakfastToday + BreakfastUsually + HowHungry + BreakfastToday*BreakfastUsually, df)
report(model_omega1)

# Boxplot
ggplot(df, aes(x=as.factor(BreakfastToday), y=omega1)) + 
  geom_boxplot() +
  geom_point() +
  xlab("Condition") +
  ylab("Value") +
  ggtitle("Gamma per condition")

### Gamma
# The higher the gamma, the more you weight your beliefs from last round
# The lower the more you weight what you observed in the last round
# Hungry people should have lower gamma since they weight more about what they observed / influence
model_gamma = lm(gamma ~ BreakfastToday + BreakfastUsually + HowHungry + BreakfastToday*BreakfastUsually, df)

report(model_gamma)

# Boxplot
ggplot(df, aes(x=as.factor(BreakfastToday), y=gamma)) + 
  geom_boxplot() +
  geom_point() +
  xlab("Condition") +
  ylab("Value") +
  ggtitle("Decay rate per condition")
