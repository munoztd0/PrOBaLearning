## R code for FOR PROBA LEARNING TASK OBIWAN
# last modified on April 2020 by David MUNOZ TORD


# PRELIMINARY STUFF ----------------------------------------
#if there is any bug please run this line below once ant then rerun the script
#invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE))

#load packages
if(!require(pacman)) {
  install.packages("pacman")
  library(pacman)
}

pacman::p_load(tidyverse, plyr,dplyr,readr,rlist, ggpubr, pracma)



# SETUP ------------------------------------------------------------------

task = 'PBlearning'


# Set working directory #change here if the switchdrive is not on your home folder
analysis_path <- file.path('~/Desktop/SwitchDrive/OBIWAN/PROBA_LEARNING/Analysis') 
figures_path  <- file.path('~/Desktop/SwitchDrive/OBIWAN/PROBA_LEARNING/Figures') 

setwd(analysis_path)



# simulate behavior -------------------------------------------------------
set.seed(666)
source('~/Desktop/SwitchDrive/OBIWAN/PROBA_LEARNING/Analysis/simulatePST.R', echo=TRUE)

# create data structure ----------------------------------------------

s = 500 #number of subjects to simulate
simDATA = c(); paramA = c(); paramB = c()
for (i in 1:s) {
  beta = rgamma(1, shape =4, scale =0.5); paramB = rbind(paramB, beta) #rgamma(1, shape =4, scale =0.5)
  alpha = rand(); paramA = rbind(paramA, alpha)  #rbeta(1, shape1=5, shape2=1.5) 1-rbeta(1, shape1=5, shape2=1.5)
  ID = rep(i, each=90)
  type = sample(rep(c(12, 34, 56), each=30))
  data = cbind(ID, type); data = as_tibble(data)
  print(paste('simulate sub', i))
  d = simulatePST(alpha, beta, data)
  simDATA = rbind(simDATA, d)
  #print(paste('done sub', i))
  
}
d$trial = 1:length(d$ID)
data_long <- gather(d, option, ev, ev1:ev6, factor_key=TRUE)
data_long$option <- factor(data_long$option, levels = c("ev1", "ev3", "ev5", "ev2", "ev4", "ev6"))
# data_long %>%
#   ggplot(aes(x = trial, y = ev)) +
#   geom_line(size=0.5) +
#   #geom_line(aes(y=pe), color ='red',size=0.5) +
#   ylim(0,1) +
#   facet_wrap("option")

# Re-estimate model's parameters  ----------------------------------------------------------
source('~/Desktop/SwitchDrive/OBIWAN/PROBA_LEARNING/Analysis/PST_Q_learning.R', echo=TRUE)
set.seed(666); Nrep = 100; k = 2 # number of free parameteres
subj = unique(simDATA$ID)
alpha = c(); beta = c(); nll= c(); ID = c(); group = c(); trials = c()
LB = c(0, 0); UB = c(1, 10) # parameters lower and upper bounds
for (s in subj) {
  data = subset(simDATA, ID == s)
  param_rep = c(); nll_rep = c()
  for (i in  1:Nrep) {
    x0 = c(rand(), rgamma(1, shape =4, scale =0.5)); #different parameter initial values to avoid local minima
    f = fmincon(x0=x0,fn=PST_q, data = data, lb = LB, ub = UB) #optimize
    param_rep = rbind(param_rep, f$par); nll_rep = rbind(nll_rep, f$value)
  }
  pos = which.min(nll_rep)
  alpha = rbind(alpha,param_rep[pos,1]); beta = rbind(beta,param_rep[pos,2]); nll = rbind(nll,nll_rep[pos]); ID = rbind(ID,s); group = rbind(group, ifelse(s > 199, 'obese', 'lean')); trials = rbind(trials,length(data$ID))
  print(paste('done subj', s))
}

dfsim = cbind(ID, alpha, beta, nll, group, trials)
colnames(dfsim) = c('ID', 'alpha','beta', 'nll', 'group', 'trials'); 
dfsim = as_tibble(dfsim); dfsim$group =as.factor(dfsim$group); dfsim$ID =as.factor(dfsim$ID)
dfsim[] <- lapply(dfsim, function(x) {if(is.character(x)) as.numeric(as.character(x)) else x})  



# check correlation between true participant's parameters and their  --------

dfsim$alphaO =paramA; dfsim$betaO = paramB; dfsim = subset(dfsim, beta != 10); dfsim = subset(dfsim, alpha != 1) # removing maxed out valued

format_pval <- function(pval){
  pval <- scales::pvalue(pval, accuracy= 0.001, add_p = TRUE)
  gsub(pattern = "(=|<)", replacement = " \\1 ", x = pval)
}

p1 = ggscatter(dfsim, x = "alphaO", y = "alpha", 
               add = "reg.line", conf.int = T,
               add.params = list(color = "black", fill = "grey", size = 0.75), xlab = "Original", ylab = "Recovered", title = 'Parameter Recovery: \u03B1 (Learning Rate)', show.legend.text = FALSE ) + 
  stat_cor(aes(label = paste(..rr.label..,format_pval(..p..), sep = "~`,`~")), label.y = 0.01,label.x = 0.55) + theme(plot.title = element_text(hjust = 0.5)); p1

p2 = ggscatter(dfsim, x = "betaO", y = "beta", p.digits=2,p.accuracy=2,
               add = "reg.line", conf.int = T,
               add.params = list(color = "black", fill = "grey", size = 0.75), xlab = "Original", ylab = "Recovered", title = 'Parameter Recovery: \u03B2 (Choice Consistency)', show.legend.text = FALSE ) + 
  stat_cor(aes(label = paste(..rr.label.., format_pval(..p..), sep = "~`,`~"), p.digits=2,p.accuracy=2), label.y = 0.1,label.x = 4) + theme(plot.title = element_text(hjust = 0.5)); p2

cairo_pdf(file.path(figures_path,'Figure_alpha_full_sim.pdf'))
print(p1)
dev.off()


cairo_pdf(file.path(figures_path,'Figure_betafull_sim.pdf'))
print(p2)
dev.off()


#correlation between sim param fit
  p3 = ggscatter(dfsim, x = "alpha", y = "beta", 
               add = "reg.line", conf.int = T,cor.coef = F,
               add.params = list(color = "black", fill = "grey", size = 0.75), xlab = "Alpha", ylab = "Beta", title = 'Fit parameter correlation', show.legend.text = FALSE )  + theme(plot.title = element_text(hjust = 0.5)); p3


cairo_pdf(file.path(figures_path,'Figure_corfull_sim.pdf'))
print(p3)
dev.off()


  
# checks ------------------------------------------------------------------


#checking expected value
# d$trial = 1:length(d$ID)
# data_long <- gather(d, option, ev, ev1:ev6, factor_key=TRUE)
# data_long$option <- factor(data_long$option, levels = c("ev1", "ev3", "ev5", "ev2", "ev4", "ev6"))
# data_long %>%
#   ggplot(aes(x = trial, y = ev)) +
#   geom_line(size=0.5) +
#   #geom_line(aes(y=pe), color ='red',size=0.5) +
#   ylim(0,1) +
#   facet_wrap("option")