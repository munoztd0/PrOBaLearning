## R code for FOR PROBA LEARNING TASK OBIWAN
# last modified on June 2020 by David MUNOZ TORD
#invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE))
# PRELIMINARY STUFF ----------------------------------------
if(!require(pacman)) {
  install.packages("pacman")
  library(pacman)
}

devtools::install_github('kongdd/Ipaper')

pacman::p_load(tidyverse, plyr, dplyr, parallel, Ipaper, hBayesDM, ggpubr, viridis, tidyBF, here, ggthemes, devtools, JWileymisc) 



pal = viridis::inferno(n=5) ; pal[6] = "#21908CFF"# specialy conceived for colorblindness # add one

options(mc.cores = parallel::detectCores(), warn= -1, readr.num_columns = 0) #to mulithread

# SETUP ------------------------------------------------------------------

task = 'PBlearning'


# Set working directory #change here if the switchdrive is not on your home folder
path <- here::i_am("Analysis/PrOBaLearning.R")
figures_path  <- here('Figures') 

# open dataset 
full <- read_csv(here("PBLearning.csv"))
info <- read_csv(here("info.csv"))


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

count_trial = data %>% group_by(subjID) %>%tally()

dataclean <- dplyr::select(data, c(subjID, type, choice, reward, Group))



### Demographics -------------------

bs = ddply(dataclean, .(subjID, Group), summarise, t = mean(Group, na.rm = TRUE)) 
demo = merge(bs, info, by.x = "subjID", by.y = "id"); demo$gender= as.factor(demo$gender)

egltable(c("age", "gender", "BMI_t1"),  g = "Group", data = demo, strict = T)

BMI = ddply(demo,~Group,summarise,mean=mean(BMI_t1),sd=sd(BMI_t1), min = min(BMI_t1), max = max(BMI_t1), n = length(unique(subjID))); BMI
AGE = ddply(demo,~Group,summarise,mean=mean(age),sd=sd(age), min = min(age), max = max(age), n = length(unique(subjID))); AGE


# MODELING ----------------------------------------------------------------


model1 = hBayesDM::pst_gainloss_Q(dataclean, niter=10000, nwarmup=2000, nchain=4, ncore=8);
model2 = hBayesDM::pst_Q(dataclean, niter=10000, nwarmup=2000, nchain=4, ncore=8); 

#or directly
#load(here("PST.RData")) # if you dont want ot recompute and go directly to stats

# evaluate models
model_fit = printFit(model1, model2, ic = "both"); model_fit #Watanabe–Akaike information criterion (WAIC), generalized version of the Akaike information criterion 
#model1 is better in both
# Model    LOOIC     WAIC L
# 1 pst_gainloss_Q 10433.68 10403.39  
# 2          pst_Q 10494.44 10454.70  

#check that MCMC samples are indeed well mixed and converged
plot(model1, type="trace", fontSize=11, ncols = 1)

#check posterior sample for hyperparameters
x =plotInd(model1, "mu_beta")+theme(axis.text.y = element_blank()); y= plotInd(model1, "mu_alpha_pos")+theme(axis.text.y = element_blank()); z =plotInd(model1, "mu_alpha_neg")+theme(axis.text.y = element_blank())
posterior = ggarrange(x, y, y, ncol = 3, labels = c("Beta", "Alpha+", "Alpha-"),  font.label = list(size = 32)); posterior 

model1$allIndPars$group = ifelse(model1$allIndPars$subjID > 199, 'obese', 'lean') #recreate groups

df = as_tibble(model1$allIndPars)


# STATS -------------------------------------------------------------------


#%%%%%%%%%% alpha Gain
BF_alpha_pos = tidyBF::bf_ttest(df, group, alpha_pos, output = "dataframe", paired = F, iterations = 50000); 

paste("cohen's d =" , BF_alpha_pos$estimate[2], "and 95%CI [", BF_alpha_pos$conf.low[2], ',', BF_alpha_pos$conf.high[2], ']') 

ttest_alpha_pos = t.test(df$alpha_pos[df$group == "lean"], df$alpha_pos[df$group == "obese"], paired=F, var.equal = F); ttest_alpha_pos



#%%%%%%%%%%%%%%%%%%%%% alpha Loss
BF_alpha_neg = tidyBF::bf_ttest(df, group, alpha_neg, output = "dataframe", paired = F, iterations = 50000); 


paste("cohen's d =" , BF_alpha_neg$estimate[2], "and 95%CI [", BF_alpha_neg$conf.low[2], ',', BF_alpha_neg$conf.high[2], ']')  

ttest_alpha_neg = t.test(df$alpha_neg[df$group == "lean"], df$alpha_neg[df$group == "obese"], paired=F, var.equal = F); ttest_alpha_neg

df$alpha_negZ = scale(df$alpha_neg)
dfout = filter(df, subjID != 109)# without huge outlier more than 3 SD

BF_alpha_neg2 = tidyBF::bf_ttest(dfout, group, alpha_neg, output = "dataframe", paired = F, iterations = 50000); 


paste("cohen's d =" , BF_alpha_neg2$estimate[2], "and 95%CI [", BF_alpha_neg2$conf.low[2], ',', BF_alpha_neg2$conf.high[2], ']')  

ttest_alpha_neg2 = t.test(dfout$alpha_neg[dfout$group == "lean"], dfout$alpha_neg[dfout$group == "obese"], paired=F, var.equal = F); ttest_alpha_neg2

#%%%%%%%%%%%%%%%%%%%%%%%%%%%% Beta

BF_beta = tidyBF::bf_ttest(df, group, beta, output = "dataframe", paired = F, iterations = 50000); 



paste("cohen's d =" , BF_beta$estimate[2], "and 95%CI [", BF_alpha_neg$conf.low[2], ',', BF_alpha_neg$conf.high[2], ']')  

ttest_beta = t.test(df$alpha_neg[df$group == "lean"], df$alpha_neg[df$group == "obese"], paired=F, var.equal = F); ttest_beta






# check posterior distrib -------------------------------------------------


# alternative way to evaluate group differences by examining the posterior distribution of group mean differences.
lean = dataclean[dataclean$Group == "C",]; obese = dataclean[dataclean$Group == "O",]; 
group1 = hBayesDM::pst_gainloss_Q(lean, niter=10000, nwarmup=2000, nchain=4, ncore=8)
group2 = hBayesDM::pst_gainloss_Q(obese, niter=10000, nwarmup=2000, nchain=4, ncore=8)

## After model fitting is complete for both groups

diffDistA = group1$parVals$mu_beta - group2$parVals$mu_beta  # group1 - group2
HDI_beta = HDIofMCMC(diffDistA)  # Compute the 95% Highest Density Interval (HDI).
e =plotHDI(diffDistA) +labs(x=expression(beta~"Gain: Lean - "~beta~"Gain: Obese" ))   # plot the group mean differences

diffDistB = group1$parVals$mu_alpha_pos - group2$parVals$mu_alpha_pos 
HDI_alpha_pos = HDIofMCMC(diffDistB) ; f= plotHDI(diffDistB) +labs(x=expression(alpha~"Gain: Lean - "~alpha~"Gain: Obese" ))

diffDistC = group1$parVals$mu_alpha_neg - group2$parVals$mu_alpha_neg  
HDI_alpha_neg = HDIofMCMC(diffDistC); g= plotHDI(diffDistC)  +labs(x=expression(alpha~"Loss: Lean - "~alpha~"Loss: Obese" ))  




# Create main plot -------------------------------------------------------------------------


fig_beta = 
  ggplot(df, aes(group, beta)) +
  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .3, size =2)+
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .3,linetype = "dashed") + 
  geom_boxplot2(aes(fill=group), width = 0.3, width.errorbar = 0.3,  fatten = NULL)  +
  geom_point(position = position_jitter(width = 0.15, seed=123 ))+
  labs(y = expression(beta~" (choice inverse temperature)" )) + 
  theme_fivethirtyeight(base_size = 28) + 
  theme(axis.title.y = element_text(size =28), panel.grid.major.x = element_blank(), plot.background = element_blank(), panel.background = element_blank())  + scale_fill_manual(name = "",values=c("lean" = alpha(pal[1], 1/3),"obese"=alpha(pal[6], 1/3)), guide = 'none') + scale_x_discrete(labels=c("Lean", "Obese"))


fig_alpha_pos = 
  ggplot(df, aes(group, alpha_pos)) +
  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .3, size =2)+
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .3,linetype = "dashed") + 
  geom_boxplot2(aes(fill=group), width = 0.3, width.errorbar = 0.3,  fatten = NULL)  +
  geom_point(position = position_jitter(width = 0.15, seed=123 ))+
  labs(y = expression(alpha~"Gain (positive learning rate)" )) + 
  ggthemes::theme_fivethirtyeight(base_size = 28) + 
  theme(axis.title.y = element_text(size =28), panel.grid.major.x = element_blank(), plot.background = element_blank(), panel.background = element_blank()) + scale_fill_manual(name = "",values=c("lean" = alpha(pal[1], 1/3),"obese"=alpha(pal[6], 1/3)), guide = 'none') + scale_x_discrete(labels=c("Lean", "Obese"))

fig_alpha_neg = 
  ggplot(dfout, aes(group, alpha_neg)) +
  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .3, size =2)+
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .3,linetype = "dashed") + 
  geom_boxplot2(aes(fill=group), width = 0.3, width.errorbar = 0.3,  fatten = NULL)  +
  geom_point(position = position_jitter(width = 0.15, seed=123 ))+
  labs(y = expression(alpha~"Loss (negative learning rate)" )) + 
  ggthemes::theme_fivethirtyeight(base_size = 28) + 
  theme(axis.title.y = element_text(size =28), panel.grid.major.x = element_blank(), plot.background = element_blank(), panel.background = element_blank()) + scale_fill_manual(name = "",values=c("lean" = alpha(pal[1], 1/3),"obese"=alpha(pal[6], 1/3)), guide = 'none') + scale_x_discrete(labels=c("Lean", "Obese")) 


figure1 = ggarrange(fig_beta, fig_alpha_pos, fig_alpha_neg, ncol = 3, labels = c("A", "B", "C"),  font.label = list(size = 32)); figure1 # + coord_fixed(ratio=1.5)

figure2 = ggarrange(e, f, g, ncol = 3, labels = c("A", "B", "C"),  font.label = list(size = 24)); figure2

cairo_pdf(file.path(figures_path,'Figure_main.pdf'))
print(figure1)
dev.off()

cairo_pdf(file.path(figures_path,'Figure_HDI.pdf'))
print(figure2)
dev.off()


# -------------------------------------------------------------------------
#betas
expr(paste("BF"["10"],  " = ", !!format(round(BF_beta$bf10[1], digits=2), nsmall = 2), ", ", widehat(italic(delta))[" median"]^" posterior",  " = ", !!format(round(median(diffDistA), digits=2), nsmall = 2), ", CI"[" 95%"]^" HDI", " [", !!format(round(HDI_beta[1], digits=2), nsmall = 2), ", ", !!format(round(HDI_beta[2], digits=2), nsmall = 2), "]"))

#alpha pos
expr(paste("BF"["10"],  " = ", !!format(round(BF_alpha_pos$bf10[1], digits=2), nsmall = 2), ", ", widehat(italic(delta))[" median"]^" posterior",  " = ", !!format(round(median(diffDistB), digits=2), nsmall = 2), ", CI"[" 95%"]^" HDI", " [", !!format(round(HDI_alpha_pos[1], digits=2), nsmall = 2), ", ", !!format(round(HDI_alpha_pos[2], digits=2), nsmall = 2), "]"))

#alpha neg2
expr(paste("BF"["10"],  " = ", !!format(round(BF_alpha_neg2$bf10[1], digits=2), nsmall = 2), ", ", widehat(italic(delta))[" median"]^" posterior",  " = ", !!format(round(median(diffDistC), digits=2), nsmall = 2), ", CI"[" 95%"]^" HDI", " [", !!format(round(HDI_alpha_neg[1], digits=2), nsmall = 2), ", ", !!format(round(HDI_alpha_neg[2], digits=2), nsmall = 2), "]"))


#save.image(file= "PST.RData")
