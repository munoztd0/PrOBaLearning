## R code for FOR PROBA LEARNING TASK OBIWAN
# last modified on April 2020 by David MUNOZ TORD
#invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE))
# PRELIMINARY STUFF ----------------------------------------
if(!require(pacman)) {
  install.packages("pacman")
  library(pacman)
}

if(!require(tidyBF)) {
  install_version("tidyBF", version = "0.3.0")
  library(tidyBF)
}

pacman::p_load(tidyverse, plyr, dplyr,readr, car, BayesFactor, sjmisc, effectsize) #whatchout to have tidyBF 0.3.0

# SETUP ------------------------------------------------------------------

task = 'PBlearning'


# Set working directory #change here if the switchdrive is not on your home folder
analysis_path <- file.path('~/Desktop/SwitchDrive/OBIWAN/PROBA_LEARNING/Analysis') 
figures_path  <- file.path('~/Desktop/SwitchDrive/OBIWAN/PROBA_LEARNING/Figures') 

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
bs_wide <- spread(bs, imcor, acc)
bs_wide$pass = c(1:length(bs_wide$Subject)) #initialize variable
for (i in  1:length(bs_wide$Subject)) {
  if((bs_wide$A[i] >= 0.65) && (bs_wide$C[i] >=  0.60) && (bs_wide$E[i] >= 0.50 )) 
    {bs_wide$pass[i] = 1} 
  else {bs_wide$pass[i] = 0}
}

data = merge(data, bs_wide[ , c("Subject", "pass")], by = "Subject", all.x=TRUE)

data = subset(data, pass == 1)

dataclean <- select(data, c(subjID, type, choice, reward, Group))
count_trial = dataclean %>% group_by(subjID, type) %>%tally()


# Time to acquire reach criterium -----------------------------------------------------

df = tally(group_by(dataclean, subjID, Group))
densityPlot(df$n)

source('~/Desktop/SwitchDrive/OBIWAN/PROBA_LEARNING/Analysis/cohen_d_ci.R', echo=F)

# unpaired Bayesian t-test
BF = tidyBF::bf_ttest(df, Group, n, output = "dataframe", paired = F, iterations = 50000); 
BF$bf10; BF$estimate; BF$conf.low; BF$conf.high
#[1] 0.3470477 [1] 8.649697 [1] -15.46683 [1] 33.81338
lean = subset(df, Group == 'C')
obese = subset(df, Group == 'O')

# unpaired frequentist t-test
Ttest = t.test(lean$n, obese$n); Ttest
# Welch Two Sample t-test
# 
# data:  lean$n and obese$n
# t = -0.72678, df = 34.595, p-value = 0.4722
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -38.57293  18.24174
# sample estimates:
#   mean of x mean of y 
# 69.73684  79.90244 
cohen_d_ci(lean$n, obese$n,  paired=F, var.equal=F, conf.level=0.95)
# Cohen's d_s Hedges' g_s 95 % CI lower bound 95 % CI upper bound
# Effect size  -0.2030256  -0.2003889  -0.7460149    0.3443384

# Plot --------------------------------------------------------------------
source('~/Desktop/SwitchDrive/OBIWAN/PROBA_LEARNING/Analysis/rainclouds.R', echo=F)

averaged_theme <- theme_bw(base_size = 32, base_family = "Helvetica")+
  theme(strip.text.x = element_text(size = 32, face = "bold"),
        strip.background = element_rect(color="white", fill="white", linetype="solid"),
        legend.position=c(.9,.9),
        legend.title  = element_text(size = 12),
        legend.text  = element_text(size = 10),
        legend.key.size = unit(0.2, "cm"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        panel.grid.major.x = element_blank() ,
        panel.grid.major.y = element_line(size=.2, color="lightgrey") ,
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size =  30),
        axis.line = element_line(size = 0.5),
        panel.border = element_blank())


pal = viridis::inferno(n=5) # specialy conceived for colorblindness
pal[6] = "#21908CFF" # add one

pp <- ggplot(df, aes(x = Group, y = n, 
                     fill = Group, color = Group)) +
  geom_flat_violin(scale = "count", trim = FALSE, alpha = .2, color = NA)+
  geom_point(aes(x = Group), alpha = .3, position = position_jitter(width = 0.05)) +
  geom_boxplot(width = 0.05 , alpha = 0.1)+
  ylab('Trials to achieve criterion')+ xlab('')+
  scale_y_continuous(expand = c(0, 0), breaks = c(seq.int(0,200, by = 50)), limits = c(0,210)) +
  scale_x_discrete(labels=c("Obese", "Lean")) +
  scale_fill_manual(values=c("C"= pal[1], "O"=  pal[6]), guide = 'none') +
  scale_color_manual(values=c("C"= pal[1], "O"=  pal[6]), guide = 'none') +
  theme_bw()


ppp <- pp + averaged_theme 
ppp

cairo_pdf(file.path(figures_path,'Figure_trialT0.pdf'))
print(ppp)
dev.off()


