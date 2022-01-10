
CC_fun <- function(nagents,ntrials,vals,parameters) { 
  
  #free parameters - we can do inference on these
  Gb1 <- parameters$Gb1
  omega1 <- parameters$omega1 #initial weighting of beliefs about others contributions in choice of own contribution, relative to prefs
  lambda <- parameters$lambda #decay rate in weighting of beliefs about others - prefs/predictions dominate over time
  gamma <- parameters$gamma #parameter weighting of beliefs about what others will contribute, relative to observed contribution
  p0 <- parameters$p0 #intercept of linear model relating preferred contributions to possible contribution values
  pbeta <- parameters$pbeta #slope of linear model relating preferred contributions to possible contribution values
  
  #simulation arrays - to be filled   
  Ga <- array(0,c(ntrials)) #observed others' contribution (mean of all others)
  Gb <- array(0,c(nagents,ntrials)) #beliefs about others' contribution (mean of all others)
  p <- array(0,c(nagents,ntrials)) #prefered contribution on each trial (independent of beliefs)
  omega <- array(0,c(nagents,ntrials)) #weighting of beliefs about others' contribution relative to prefs
  c <- array(0,c(nagents,ntrials)) #actual contributions
  
  # agents preferences - assumed to be a linear function of possible values - linear function has p0 and pbeta as params
  pvals <- array(0,c(nagents,length(vals)))
  for (n in 1:nagents) {
    pvals[n,] <- p0[n] + (pbeta[n]*vals) #vector of preferred contributions for each possible value - assume linear relationship
    # this deviates from paper, which also has "triangular" realationships for some (i.e. increase co-operation up to a value then <)
  }
  
  # set omega as starting weighting for beliefs relative to preferences as parameter
  omega[,1] <- omega1
  
  # set starting beliefs about what others will contribute
  Gb[,1] <- Gb1
  
  # set average first trial contribution to average belief, assume full co-operation at outset. 
  # Reasonable simplification.
  c[,1] <- Gb1
  Ga[1] <- mean(Gb1)
  
  # run simulation
  for (t in 2:ntrials) {
    
    for (n in 1:nagents) {
    
      # update beliefs about what others contribute as a weighted average of prior beliefs and observed contributions
      #from page 549 "subjects belief in a given period is a weighted average of what he or she believed 
      #about others in the previous period and his or her observations of others' contributions" 
      Gb[n,t] <- (gamma[n]*(Gb[n,t-1]))+((1-gamma[n])*(Ga[t-1]))
      
      # determine what people "predict" or "prefer" to contribute, given the group contribution, outside
      # of, or independent of, their willingness to co-operate. Index preferences relative to believed contribution of others
      p[n,t] <- pvals[n,round(Gb[n,t])]
      
      # update relative weighting of beliefs about others contributions, using decay function
      #page 550 "in later periods, predicted (i.e. preferred) contibution becomes more important than belief"
      omega[n,t] <- omega[n,t-1]*(1-lambda[n]) #1 - lamba - because lambda is a decay rate
      # departure from paper - possible model innovation
      
      #page 550 subjects contribute a weighted average of predicted (i.e. preferred) contribution and belief
      c[n,t] <- ceiling((omega[n,t])*Gb[n,t] + ((1-omega[n,t])*p[n,t])) 
      
    }
    # recode average contribution as observed belief about how much each agent contributed Ga
    Ga[t] <- sum(c[,t])/nagents
    
  }
  
  #return contribution and weighting of beliefs/preferences as result of function.  
  result <- list(c=c,
                 omega=omega)
  
  return(result)

}
