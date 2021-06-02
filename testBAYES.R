

## R code for FOR PROBA LEARNING TASK OBIWAN
# last modified on April 2020 by David MUNOZ TORD
#invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE))
# PRELIMINARY STUFF ----------------------------------------
if(!require(pacman)) {
  install.packages("pacman")
  library(pacman)
}

if(!require(bmsR)) {
  install.packages("remotes")
  remotes::install_github("mattelisi/mlisi") # required for some dependancies
  remotes::install_github("mattelisi/bmsR")
  library(bmsR)
}


pacman::p_load(tidyverse, parallel, hBayesDM, ggpubr,viridis, tidyBF) 


pal = viridis::inferno(n=5) # specialy conceived for colorblindness
pal[6] = "#21908CFF" # add one

options(mc.cores = parallel::detectCores(),warn=-1) #to mulithread

# SETUP ------------------------------------------------------------------

task = 'PBlearning'


# Set working directory #change here if the switchdrive is not on your home folder
analysis_path <- file.path('~/Desktop/SwitchDrive/OBIWAN/PROBA_LEARNING/Analysis') 

setwd(analysis_path)

# open dataset 
full <- read_csv("~/Desktop/SwitchDrive/OBIWAN/PROBA_LEARNING/PBLearning.csv")




# Preprocess --------------------------------------------------------------

data  <- subset(full, Session == 1) #subset #only session one 
data  <- subset(data, Phase == 'proc1') #only learning phase

#factorize and rename
data$type = as.factor(revalue(data$imcor, c(A="AB", C="CD", E="EF")))
data$reward = revalue(data$feedback, c(Negatif=0, Positif=1))
data$side = revalue(data$Stim.RESP, c(x='L', n='R'))
data$subjID = data$Subject

#This loop is there to transform everything into one column "choice"
#this column takes a 1 if the action was to choose either A, C or E 
#and takes 0 if the response is either B, D or F (and that independently of the side)

data$choice = c(1:length(data$Trial)) #initialize variable
for (i in  1:length(data$Trial)) {
  if((data$side[i] == 'L')&(data$img[i] == 'A' || data$img[i] == 'C' || data$img[i] == 'E')) {
    data$choice[i] = 1
  } else if ((data$side[1] == 'R')&(data$imd[i] == 'A' || data$imd[i] == 'C' || data$imd[i] == 'E')) {
    data$choice[i] = 1}
  else {
    data$choice[i] = 0}
}

data$reward = as.numeric(data$reward)
data$type = revalue(data$type, c(AB=12, CD=34, EF=56))
data$type = as.numeric(as.character(data$type))


bs = ddply(data, .(Subject, imcor), summarise, acc = mean(Stim.ACC, na.rm = TRUE)) 

# Crtierium chose A at 65%, C at 60% and E at 50% and min 30 trials.
# bs_wide <- spread(bs, imcor, acc)
# bs_wide$pass = c(1:length(bs_wide$Subject)) #initialize variable
# for (i in  1:length(bs_wide$Subject)) {
#   if((bs_wide$A[i] >= 0.65) && (bs_wide$C[i] >=  0.60) && (bs_wide$E[i] >= 0.50 )) 
#   {bs_wide$pass[i] = 1} 
#   else {bs_wide$pass[i] = 0}
# }
# 
# data = merge(data, bs_wide[ , c("Subject", "pass")], by = "Subject", all.x=TRUE)
# length(unique(data$Subject))
# data = subset(data, pass == 1)
# length(unique(data$Subject))

count_trial = data %>% group_by(subjID) %>%tally()

dataclean <- select(data, c(subjID, type, choice, reward, Group))



model1 = hBayesDM::pst_gainloss_Q(dataclean, niter=1000, nwarmup=200, nchain=4, ncore=8); k1 = 3
model2 = hBayesDM::pst_Q(dataclean, niter=1000, nwarmup=200, nchain=4, ncore=8); k2 = 2

#load("~/OBIWAN/DERIVATIVES/BEHAV/PBL_OBIWAN_T0.RData") # if you dont want ot recompute and go directly to stats

log_ind1 = t(summarise_each(as_tibble(test1[["parVals"]][["log_lik"]]), funs(mean)))
log_ind2 = t(summarise_each(as_tibble(test2[["parVals"]][["log_lik"]]), funs(mean)))

AIC1 = -2*(-log_ind1) + 2*k1
AIC2 = -2*(-log_ind1) + 2*k2

BIC1 =  k1*log(count_trial$n) -2*(log_ind1) #Bayesian information criterion for each subj
BIC2 = k2*log(count_trial$n) -2*(log_ind2) #Bayesian information criterion for each subj

model_AIC = cbind(AIC1, AIC2)
model_BIC = cbind(BIC1, BIC2)

#Variational Bayesian Analysis (Daunizeau & Rigoux)
#a list with the posterior estimates of the Dirichlet parameters (alpha), the expected model frequencies (r), the exceedance probabilities (xp), the Bayesian Omnibus Risk (bor), and the protected exceedance probabilities (pxp). 
VBA = VB_bms(model_BIC) ; VBA  #model 1 is "better"


model_fit = printFit(test1, test2, ic = "waic") #Watanabe–Akaike information criterion (WAIC), generalized version of the Akaike information criterion 


test1$allIndPars$group = ifelse(test1$allIndPars$subjID > 199, 'obese', 'lean')

df = as_tibble(test1$allIndPars)


#%%%%%%%%%% alpha Gain
BF_alpha_pos = tidyBF::bf_ttest(df, group, alpha_pos, output = "dataframe", paired = F, iterations = 50000); 

fig_aplha_pos = 
  ggplot(df, aes(group, alpha_pos)) +
  geom_boxplot(aes(fill=group), alpha = 0.3, outlier.shape = 4) + 
  geom_point(position = position_jitter(width = 0.15, seed=123 ))+
  labs(caption  = expr(paste("BF"["10"],  " = ", !!format(round(BF_alpha_pos$bf10[1], digits=2), nsmall = 2), ", ", widehat(italic(delta))[" median"]^" posterior",  " = ", !!format(round(BF_alpha_pos$estimate[1], digits=2), nsmall = 2), ", CI"[" 95%"]^" HDI", " [", !!format(round(BF_alpha_pos$conf.low[1], digits=2), nsmall = 2), ", ", !!format(round(BF_alpha_pos$conf.high[1], digits=2), nsmall = 2), "]")), y = expression(alpha~"Gain (positive learning rate)" )) + 
  ggthemes::theme_fivethirtyeight(base_size = 32) + 
  theme(axis.title.y = element_text(size =32), panel.grid.major.x = element_blank() )  + scale_fill_manual(name = "",values=c("lean" = pal[1],"obese"=pal[6]), guide = 'none') + scale_x_discrete(labels=c("Lean", "Obese")) 

paste("cohen's d =" , BF_alpha_pos$estimate[2], "and 95%CI [", BF_alpha_pos$conf.low[2], ',', BF_alpha_pos$conf.high[2], ']')  #  and 95%CI

ttest_alpha_pos = t.test(df$alpha_pos[df$group == "lean"], df$alpha_pos[df$group == "obese"], paired=F, var.equal = F); ttest_alpha_pos



#%%%%%%%%%%%%%%%%%%%%% alpha Loss
BF_alpha_neg = tidyBF::bf_ttest(df, group, alpha_neg, output = "dataframe", paired = F, iterations = 50000); 

fig_aplha_neg = 
  ggplot(df, aes(group, alpha_neg)) +
  geom_boxplot(aes(fill=group), alpha = 0.3, outlier.shape = 4) + 
  geom_point(position = position_jitter(width = 0.15, seed=123 ))+
  labs(caption  = expr(paste("BF"["10"],  " = ", !!format(round(BF_alpha_neg$bf10[1], digits=2), nsmall = 2), ", ", widehat(italic(delta))[" median"]^" posterior",  " = ", !!format(round(BF_alpha_neg$estimate[1], digits=2), nsmall = 2), ", CI"[" 95%"]^" HDI", " [", !!format(round(BF_alpha_neg$conf.low[1], digits=2), nsmall = 2), ", ", !!format(round(BF_alpha_neg$conf.high[1], digits=2), nsmall = 2), "]")), y = expression(alpha~"Loss (negative learning rate)" )) + 
  ggthemes::theme_fivethirtyeight(base_size = 32) + 
  theme(axis.title.y = element_text(size =32), panel.grid.major.x = element_blank() )  + scale_fill_manual(name = "",values=c("lean" = pal[1],"obese"=pal[6]), guide = 'none') + scale_x_discrete(labels=c("Lean", "Obese")) 


paste("cohen's d =" , BF_alpha_neg$estimate[2], "and 95%CI [", BF_alpha_neg$conf.low[2], ',', BF_alpha_neg$conf.high[2], ']')  #  and 95%CI

ttest_alpha_neg = t.test(df$alpha_neg[df$group == "lean"], df$alpha_neg[df$group == "obese"], paired=F, var.equal = F); ttest_alpha_neg


BF_beta = tidyBF::bf_ttest(df, group, beta, output = "dataframe", paired = F, iterations = 50000); 


ggplot(df, aes(group, beta)) +
  geom_boxplot() + 
  geom_point(position = position_jitter(width = 0.15))+
  #geom_density() + # two-sample t-test results in an expression
  labs(subtitle = expr(paste("BF"["01"],  " = ", !!format(round(BF_alpha_pos$bf10[1], digits=2), nsmall = 2), ", ", widehat(italic(delta))["median"]^"posterior",  " = ", !!format(round(BF_beta$estimate[1], digits=2), nsmall = 2), ", CI"["95%"]^"HDI", " [", !!format(round(BF_beta$conf.low[1], digits=2), nsmall = 2), ", ", !!format(round(BF_beta$conf.high[1], digits=2), nsmall = 2), "]", ", ", italic("r")["Cauchy"]^"JZS", " = ", "0.71"))) + ggthemes::theme_fivethirtyeight()

lean = subset(dataclean, Group == 'C')
obese = subset(dataclean, Group == 'O')


group1_1 = hBayesDM::pst_gainloss_Q(lean, niter=10000, nwarmup=2000, nchain=4, ncore=8)

plot(group1_1, type="trace", fontSize=11) 
plot(group1_1)
rhat(group1_1)
plotInd(group1_1)

group2_1 = hBayesDM::pst_gainloss_Q(obese, niter=10000, nwarmup=2000, nchain=4, ncore=8)

plot(group2_1, type="trace", fontSize=11) 
plot(group2_1)
rhat(group2_1)
x = plotInd(group2_1, pars = c("beta"), fontSize=11); q <- ggplot_build(x)
q[["plot"]][["data"]][["params"]] = ""; p <- ggplot_gtable(q)


y = plotInd(group2_1, pars = c("mu_beta"))
figure1 = ggarrange(x, y, ncol = 2); figure1
plotInd(group2_1, pars = c("mu_alpha_pos","mu_alpha_neg"))
plotInd(group2_1, pars = c("mu_alpha_neg"))

a_1 = 1/(1 + exp(-a))

b_1 = 1/(1 + exp(-b))





group1_2 = hBayesDM::pst_Q(lean, niter=10000, nwarmup=2000, nchain=4, ncore=8)

plot(group1_2, type="trace", fontSize=11) 
plot(group1_2)
rhat(group1_2)
plotInd(group1_2)

group2_2 = hBayesDM::pst_Q(obese, niter=10000, nwarmup=2000, nchain=4, ncore=8)

plot(group2_2, type="trace", fontSize=11) 
plot(group2_2)
rhat(group2_2)
plotInd(group2_2)


    

printFit(group1_1, group1_2, ic = "looic") # pst_gainloss_Q best 
printFit(group2_1, group2_2, ic = "looic") # pst_gainloss_Q best



## After model fitting is complete for both groups,
## evaluate the group difference (e.g., on the 'pi' parameter) by examining the posterior distribution of group mean differences.

diffDist = group1_1$parVals$mu_beta - group2_1$parVals$mu_beta  # group1 - group2
HDIofMCMC(diffDist)  # Compute the 95% Highest Density Interval (HDI).
plotHDI(diffDist)    # plot the group mean differences

diffDist = group1_1$parVals$mu_alpha_pos - group2_1$parVals$mu_alpha_pos  # group1 - group2
HDIofMCMC(diffDist)  # Compute the 95% Highest Density Interval (HDI).
plotHDI(diffDist)    # plot the group mean differences

diffDist = group1_1$parVals$mu_alpha_neg - group2_1$parVals$mu_alpha_neg  # group1 - group2
HDIofMCMC(diffDist)  # Compute the 95% Highest Density Interval (HDI).
plotHDI(diffDist)    # plot the group mean differences


#%%%%%%%%%%%%%%%%%%%


diffDist = group1_2$parVals$mu_beta - group2_2$parVals$mu_beta  # group1 - group2
HDIofMCMC(diffDist)  # Compute the 95% Highest Density Interval (HDI).
plotHDI(diffDist)    # plot the group mean differences

diffDist = group1_2$parVals$mu_alpha - group2_2$parVals$mu_alpha  # group1 - group2
HDIofMCMC(diffDist)  # Compute the 95% Highest Density Interval (HDI).
plotHDI(diffDist)    # plot the group mean differences