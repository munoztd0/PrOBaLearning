#-------------------------------------------------------------------------%
#   SIMULATION and RUN OF THE PST TASK    %
#   David Munoz Tord
#-------------------------------------------------------------------------%
# PARAMETERS
# lrG      = learning rate gain [0 - 1]
# lrL      = learning rate loss [0 - 1]
# b        = softmax inverse temperature [0 - 10]


simulatedualPST <- function(alphaG, alphaL, beta, data){
  
  #tidy up
  N       <- length(data$type)
  #litle trick to get the types of trials from 12 34 56
  option1 <- data$type %/% 10  #modulo 12 = 1, 34 = 3, 56 = 5
  option2 <- data$type %% 10  #numerical division 12 = 2, 34 = 4, 56 = 6
  choice  <- rep(0, N); reward  <- rep(0, N)
  
  ev = rep(0.5, 6) #initialization
  pesim = c(); evsim = c()
  #loop through trials
  for (t in  1:N) {
    
    p =  exp(beta*ev[option1[t]])/(exp(beta*ev[option1[t]]) + (exp(beta*ev[option2[t]]))) #compute probability of chosing option 1
    
    
    choice[t] = ifelse(runif(1) < p, 1, 0) #probabilistic outcome
    
    
    if (choice[t] == 1) {
      co = option1[t]
      cof = option2[t]
    } else  {
      co = option2[t]
      cof = option1[t]
    }
    
    if (co == 1) {
      reward[t] = ifelse(sample(1:10, 1) <= 8, 1, 0)
    } else if (co == 2) {
      reward[t] = ifelse(sample(1:10, 1) <= 2, 1, 0)
    } else if (co == 3) {
      reward[t] = ifelse(sample(1:10, 1) <= 7, 1, 0)
    } else if (co == 4) {
      reward[t] = ifelse(sample(1:10, 1) <= 3, 1, 0)
    } else if (co == 5) {
      reward[t] = ifelse(sample(1:10, 1) <= 6, 1, 0)
    } else if (co == 6) {
      reward[t] = ifelse(sample(1:10, 1) <= 4, 1, 0)
    }
    
    pe = reward[t] - ev[co]  #prediction error
    pesim = rbind(pesim, pe)
    evsim = rbind(evsim, ev)
    alpha = ifelse(reward[t] > 0, alphaG, alphaL) # differential lr for loss and gain
    ev[co] = ev[co] + alpha * pe  #Q or expected value update
    #ev[cof] = 1-ev[co]
  }
  
  data$choice = choice
  data$reward = reward
  data$pe = pesim
  data$ev1 = evsim[,1]
  data$ev2 = evsim[,2]
  data$ev3 = evsim[,3]
  data$ev4 = evsim[,4]
  data$ev5 = evsim[,5]
  data$ev6 = evsim[,6]
  return(invisible(data))
}



PST_q_dual <- function(x, data){
  alphaG = x[1]
  alphaL = x[2]
  beta  = x[3]
  #tidy up
  N       <- length(data$type)
  #litle trick to get the types of trials from 12 34 56
  option1 <- data$type %/% 10  #modulo 12 = 1, 34 = 3, 56 = 5
  option2 <- data$type %% 10  #numerical division 12 = 2, 34 = 4, 56 = 6
  choice  <- data$choice
  reward  <- data$reward
  
  ev = rep(0, 6) #initialization
  lik = c()
  #loop through trials
  for (t in  1:N) {
    
    co = ifelse(choice[t] > 0, option1[t], option2[t])
    
    pe = reward[t] - ev[co]  #prediction error
    alpha = ifelse(reward[t] > 0, alphaG, alphaL) # differential lr for loss and gain
    ev[co] = ev[co] + alpha * pe  #Q or expected value update
    
    dQ = ev[option1[t]] - ev[option2[t]];
    phi = 1/(1+exp(-beta *dQ))
    
    #find the log likelihood for the choice made in each trial
    log_likelihood <- (choice[t] * log(phi)) + ((1-choice[t]) * log(1-phi))
    lik =  rbind(lik, log_likelihood)
  }
  NeglogLik <- -sum(lik)
  return(invisible(NeglogLik))
}
