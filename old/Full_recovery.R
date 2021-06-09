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

pacman::p_load(tidyverse, plyr,dplyr,readr,rlist, ggpubr, NlcOptim,pracma, here, foreach, parallel, doSNOW, future, corrplot, RColorBrewer)



# SETUP ------------------------------------------------------------------

task = 'PBlearning'


# Set working directory #change here if the switchdrive is not on your home folder
analysis_path <- here::i_am()
  file.path('~/Desktop/SwitchDrive/OBIWAN/PROBA_LEARNING/Analysis') 
figures_path  <- file.path('~/Desktop/SwitchDrive/OBIWAN/PROBA_LEARNING/Figures') 

setwd(analysis_path)


options(mc.cores = parallel::detectCores())

cores=parallel::detectCores()
# cl <- future::makeCluster(cores[1]/2, outfile="") #for my machine
# #cl <- parallel::makeCluster(cores[1]/2, outfile="") # for centOS cluster
# registerDoSNOW(cl)

# simulate behavior -------------------------------------------------------
set.seed(666)
source('~/Desktop/SwitchDrive/OBIWAN/PROBA_LEARNING/Analysis/PST_dualQ_learning.R', echo=F)

# create data structure ----------------------------------------------

s = 50 #number of subjects to simulate
simDATA = c(); 

paramA = rand(s,1);  #rbeta(1, shape1=5, shape2=1.5) 1-rbeta(1, shape1=5, shape2=1.5)
paramB = rand(s,1);  #rbeta(1, shape1=5, shape2=1.5) 1-rbeta(1, shape1=5, shape2=1.5)
paramC = rgamma(s, shape =4, scale =0.5); #rgamma(1, shape =4, scale =0.5)


create <- function(i){ 
  alphaG = paramA[i]
  alphaL = paramB[i]
  beta = paramC[i]
  ID = rep(i, each=180)
  type = sample(rep(c(12, 34, 56), each=60))
  data = cbind(ID, type); data = as_tibble(data)
  d = simulatedualPST(alphaG, alphaL, beta, data)
  return(d)
}

pb <- txtProgressBar(min = 1, max = s, style = 3)

simDATA <- foreach(i = 1:s, .combine = 'rbind', .packages = c("pracma", "tidyverse")) %dopar% { 
  sol = create(i) 
  setTxtProgressBar(pb, i) 
  return(sol)
}

# d = subset(simDATA, ID%in%1:10)
# d$trial = rep(1:180, 10)
# data_long <- gather(d, option, ev, ev1:ev6, factor_key=TRUE)
# data_long$option <- factor(data_long$option, levels = c("ev1", "ev3", "ev5", "ev2", "ev4", "ev6"))
# 
# data_long %>%
#   ggplot(aes(x = trial, y = ev, color =as.factor(ID))) +
#   geom_line(size=0.5) +
#   #geom_line(aes(y=pe), color ='red',size=0.5) +
#   ylim(0,1) +
#   facet_wrap("option")

simDATA$subjID = simDATA$ID
modelSIM = hBayesDM::pst_gainloss_Q(simDATA, niter=5000, nwarmup=1000, nchain=3, ncore=8);

# Re-estimate model's parameters  ----------------------------------------------------------
set.seed(666); Nrep = 10; k = 2 # number of free parameteres
subj = unique(simDATA$ID)
dft = c(); LB = c(0, 0, 0); UB = c(1, 1, 10) # parameters lower and upper bounds

sim <- function(simDATA, s){  
  data = subset(simDATA, ID == s)
  param_rep = c(); nll_rep = c()
  for (i in  1:Nrep) {
    x0 = c(rand(),rand(), rgamma(1, shape =4, scale =0.5)); #different parameter initial values to avoid local minima
    f = fmincon(x0=x0,fn=PST_q_dual, data = data, lb = LB, ub = UB) #optimize
    param_rep = rbind(param_rep, f$par); nll_rep = rbind(nll_rep, f$value)
  }
  pos = which.min(nll_rep)
  dft$alphaG = param_rep[pos,1]; dft$alphaL = param_rep[pos,2]; dft$beta = param_rep[pos,3]; dft$nll = nll_rep[pos]; dft$ID = s;  dft$trials =length(data$ID)
  return(dft)
}


dfsimO = foreach(s = subj, .combine = 'rbind', .packages="pracma", .export ="dft", .errorhandling = 'pass') %dopar% {
  solu = tryCatch(sim(simDATA , s),
                 error = function(e) return(list(alphaG = 0, alphaL = 0, beta = 0, nll= 0, ID=s, trials=180)))
  setTxtProgressBar(pb, s) 
  return(solu)}


dfsim <- data.frame(matrix(unlist(dfsimO), ncol =size(t(dfsimO)), byrow = F))
colnames(dfsim) = c('alphaG','alphaL','beta', 'nll','ID', 'trials'); 
dfsim = as_tibble(dfsim);  dfsim$ID =as.factor(dfsim$ID)
#dfsim[] <- lapply(dfsim, function(x) {if(is.character(x)) as.numeric(as.character(x)) else x})  

save(dfsimO, "dfsim.RData")

# check correlation between true participant's parameters and their  --------

dfsim$alphaGO =c(paramA); dfsim$alphaLO = c(paramB); dfsim$betaO = paramC; dfsim = filter(dfsim, beta < 9.9999 & beta > 1e-3); dfsim = filter(dfsim, alphaG < 0.9999 & alphaG > 1e-3); dfsim = subset(dfsim, alphaL < 0.9999 &  alphaL > 1e-3)  # removing maxed out valued



p1 = ggscatter(dfsim, x = "alphaGO", y = "alphaG", 
               add = "reg.line", conf.int = T,
               add.params = list(color = "black", fill = "grey", size = 0.75), xlab = "Original", ylab = "Recovered", title = '\u03B1 Gain (Positive Learning Rate)', show.legend.text = FALSE )  + theme(plot.title = element_text(hjust = 0.5)); p1

p2 = ggscatter(dfsim, x = "alphaLO", y = "alphaL", p.digits=2,p.accuracy=2,
               add = "reg.line", conf.int = T,
               add.params = list(color = "black", fill = "grey", size = 0.75), xlab = "Original", ylab = "Recovered", title = '\u03B1 Loss (Negative Learning Rate)', show.legend.text = FALSE ) + theme(plot.title = element_text(hjust = 0.5)); p2

p3 = ggscatter(dfsim, x = "betaO", y = "beta", p.digits=2,p.accuracy=2,
               add = "reg.line", conf.int = T,
               add.params = list(color = "black", fill = "grey", size = 0.75), xlab = "Original", ylab = "Recovered", title = '\u03B2 (Choice Consistency)', show.legend.text = FALSE )  + theme(plot.title = element_text(hjust = 0.5)); p3

cairo_pdf(file.path(figures_path,'Figure_alphaG_full_sim.pdf'))
print(p1)
dev.off()


cairo_pdf(file.path(figures_path,'Figure_alphaL_full_sim.pdf'))
print(p2)
dev.off()


cairo_pdf(file.path(figures_path,'Figure_beta_full_sim.pdf'))
print(p3)
dev.off()


#correlation between param fit
corrplot(cor(dfsim[c(1:3)]), type="upper", order="hclust",
         col=brewer.pal(n=8, name="PuOr"))


#stats
test1 = summary(lm(alpha ~ alphaGO, data=dfsim)); test1$r.squared 
tesBF1 = BayesFactor::lmBF(alphaG ~ alphaGO, data=dfsim); tesBF1  

test2 = summary(lm(alphaL ~ alphaLO, data=dfsim)); test2$r.squared 
tesBF12= BayesFactor::lmBF(alphaL ~ alphaLO, data=dfsim); tesBF2  

test3 = summary(lm(beta ~ betaO, data=dfsim)); test3$r.squared 
tesBF3 = BayesFactor::lmBF(beta ~ betaO, data=dfsim); testBF3  

  
# checks ------------------------------------------------------------------

# format_pval <- function(pval){
#   pval <- scales::pvalue(pval, accuracy= 0.001, add_p = TRUE)
#   gsub(pattern = "(=|<)", replacement = " \\1 ", x = pval)
# }

# + stat_cor(aes(label = paste(..rr.label.., format_pval(..p..), sep = "~`,`~"), p.digits=2,p.accuracy=2), label.y = 0.1,label.x = 3)


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