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

pacman::p_load(tidyverse, plyr,dplyr,readr,rlist, ggpubr, pracma, here, foreach, parallel, corrplot, RColorBrewer, hBayesDM, BayesFactor)



# SETUP ------------------------------------------------------------------

task = 'PBlearning'


# Set working directory #change here if the switchdrive is not on your home folder
path <- here::i_am("Analysis/RECOVERY.R")
figures_path  <- here('Figures') 


options(mc.cores = parallel::detectCores(), warn=-1) #to mulithread

# simulate behavior -------------------------------------------------------
set.seed(666)
source( here('Analysis', 'PST_dualQ_learning.R') , echo=F)

# create data structure ----------------------------------------------

s = 500 #number of subjects to simulate
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


simDATA$subjID = simDATA$ID
modelSIM = hBayesDM::pst_gainloss_Q(simDATA, niter=5000, nwarmup=1000, nchain=4, ncore=8);

dfsim = modelSIM$allIndPars 
save(dfsim, file= "dfsim.RData")

# check correlation between original and recovered parameters  --------

dfsim$alphaGO =c(paramA); dfsim$alphaLO = c(paramB); dfsim$betaO = paramC; 


p1 = ggscatter(dfsim, x = "alphaGO", y = "alpha_pos", 
               add = "reg.line", conf.int = T,
               add.params = list(color = "black", fill = "grey", size = 0.75), xlab = "Original", ylab = "Recovered", title = '\u03B1 Gain (Positive Learning Rate)', show.legend.text = FALSE )  + theme(plot.title = element_text(hjust = 0.5))+  scale_y_continuous(breaks = c(seq.int(0,1, by = 0.25)), limits = c(0,1)) ; p1

p2 = ggscatter(dfsim, x = "alphaLO", y = "alpha_neg", p.digits=2,p.accuracy=2,
               add = "reg.line", conf.int = T,
               add.params = list(color = "black", fill = "grey", size = 0.75), xlab = "Original", ylab = "Recovered", title = '\u03B1 Loss (Negative Learning Rate)', show.legend.text = FALSE ) + theme(plot.title = element_text(hjust = 0.5))+  scale_y_continuous(breaks = c(seq.int(0,1, by = 0.25)), limits = c(0,1)); p2

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
test1 = summary(lm(alpha_pos ~ alphaGO, data=dfsim)); test1$r.squared 
tesBF1 = BayesFactor::lmBF(alpha_pos ~ alphaGO, data=dfsim); tesBF1  

test2 = summary(lm(alpha_neg ~ alphaLO, data=dfsim)); test2$r.squared 
tesBF12= BayesFactor::lmBF(alpha_neg ~ alphaLO, data=dfsim); tesBF2  

test3 = summary(lm(beta ~ betaO, data=dfsim)); test3$r.squared 
tesBF3 = BayesFactor::lmBF(beta ~ betaO, data=dfsim); testBF3  

#save.image(file= "SIM.RData")