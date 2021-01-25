# PRELIMINARY STUFF ----------------------------------------
if(!require(pacman)) {
  install.packages("pacman", "devtools")
  library(devtools)
  library(pacman)
}

if(!require(tidyBF)) {
  install_version("tidyBF", version = "0.3.0")
  library(tidyBF)
}

pacman::p_load(tidyverse, plyr,dplyr,readr, car, BayesFactor, sjmisc, effectsize, ggpubr) #whatchout to have tidyBF 0.3.0

# Set working directory #change here if the switchdrive is not on your home folder
analysis_path <- file.path('~/Desktop/SwitchDrive/OBIWAN/PROBA_LEARNING/Analysis') 
figures_path  <- file.path('~/Desktop/SwitchDrive/OBIWAN/PROBA_LEARNING/Figures') 

setwd(analysis_path)

# STATS -------------------------------------------------------------------

#load data
load("~/Desktop/SwitchDrive/OBIWAN/PROBA_LEARNING/PBL_OBIWAN_T0.RData")

source('~/Desktop/SwitchDrive/OBIWAN/PROBA_LEARNING/Analysis/cohen_d_ci.R', echo=F)

lean = subset(df, group == 'lean')
obese = subset(df, group == 'obese')

#### For alpha 

# unpaired Bayesian t-test
BF_alpha = tidyBF::bf_ttest(df, group, alpha, output = "dataframe", paired = F, iterations = 50000); 
BF_alpha$bf10; BF_alpha$estimate; BF_alpha$conf.low; BF_alpha$conf.high
# [1] 0.3172115 # [1] 0.01908016 # [1] -0.05750989 # [1] 0.09521447

# unpaired frequentist t-test
TtestA = t.test(lean$alpha, obese$alpha); TtestA
# Welch Two Sample t-test
# 
# data:  lean$alpha and obese$alpha
# t = -0.70366, df = 57.433, p-value = 0.4845
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.09224584  0.04426763
# sample estimates:
#   mean of x mean of y 
# 0.1718557 0.1958448

cohen_d_ci(lean$alpha, obese$alpha,  paired=F, var.equal=F, conf.level=0.95)
# Cohen's d_s Hedges' g_s 95 % CI lower bound 95 % CI upper bound
  # Effect size  -0.1561702   -0.154142    -0.7395484    0.3506508

#### For beta

# unpaired Bayesian t-test
BF_beta = tidyBF::bf_ttest(df, group, beta, output = "dataframe", paired = F, iterations = 50000); 
BF_beta$bf10; BF_beta$estimate; BF_beta$conf.low; BF_beta$conf.high
  # [1] 0.2784632 # [1] 0.038279 # [1] -1.179893 # [1] 1.205852

# unpaired frequentist t-test
TtestB = t.test(lean$beta, obese$beta); TtestB
# Welch Two Sample t-test
# 
# data:  lean$beta and obese$beta
# t = -0.025379, df = 36.212, p-value = 0.9799
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -1.334214  1.301228
# sample estimates:
#   mean of x mean of y 
# 7.598949  7.615442 

cohen_d_ci(lean$beta, obese$beta,  paired=F, var.equal=F, conf.level=0.95)
# Cohen's d_s  Hedges' g_s 95 % CI lower bound 95 % CI upper bound
# Effect size -0.006959211 -0.006868832          -0.5509596           0.5369334

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

pp <- ggplot(df, aes(x = group, y = alpha, 
                     fill = group, color = group)) +
  geom_flat_violin(scale = "count", trim = FALSE, alpha = .2, color = NA)+
  geom_point(aes(x = group), alpha = .3, position = position_jitter(width = 0.05)) +
  geom_boxplot(width = 0.05 , alpha = 0.1)+
  ylab('\u03B1 (Learning Rate)')+ xlab('')+
  scale_y_continuous(expand = c(0, 0), breaks = c(seq.int(0,1, by = 0.25)), limits = c(0,1)) +
  scale_x_discrete(labels=c("Obese", "Lean")) +
  scale_fill_manual(values=c("lean"= pal[1], "obese"=  pal[6]), guide = 'none') +
  scale_color_manual(values=c("lean"= pal[1], "obese"=  pal[6]), guide = 'none') +
  theme_bw()


ppp <- pp + averaged_theme 
ppp

cairo_pdf(file.path(figures_path,'Figure_alpha.pdf'))
print(ppp)
dev.off()


pp <- ggplot(df, aes(x = group, y = beta, 
                     fill = group, color = group)) +
  geom_flat_violin(scale = "count", trim = FALSE, alpha = .2, color = NA)+
  geom_point(aes(x = group), alpha = .3, position = position_jitter(width = 0.05)) +
  geom_boxplot(width = 0.05 , alpha = 0.1)+
  ylab('\u03B2 (Choice Consistency)')+ xlab('')+
  scale_y_continuous(expand = c(0, 0), breaks = c(seq.int(0,10, by = 2)), limits = c(0,12)) +
  scale_x_discrete(labels=c("Obese", "Lean")) +
  scale_fill_manual(values=c("lean"= pal[1], "obese"=  pal[6]), guide = 'none') +
  scale_color_manual(values=c("lean"= pal[1], "obese"=  pal[6]), guide = 'none') +
  theme_bw()


ppp <- pp + averaged_theme 
ppp

cairo_pdf(file.path(figures_path,'Figure_beta.pdf'))
print(ppp)
dev.off()




p1 = ggscatter(df, x = "alpha", y = "beta", p.digits=2,p.accuracy=2,
               add = "reg.line", conf.int = T, cor.coef = TRUE, 
               cor.coeff.args = list(method = "pearson", label.x = 0.6,label.y = 5, label.sep = "\n"),   add.params = list(color = "black", fill = "grey", size = 0.75), xlab = "Alpha", ylab = "Beta", title = 'Fit Parameter Correlation', show.legend.text = T ) ; p1

cairo_pdf(file.path(figures_path,'Figure_cor.pdf'))
print(p1)
dev.off()



