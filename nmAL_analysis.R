library(ggplot2)
library(dplyr)

# load in note matrices for each individual
load("new_note_matrices.RData")

# calculate normalized arc length = arc length / note duration
for (i in 1:length(redAL_nm_41_13_1)){
  nAL <- c()
  for (j in 1:nrow(redAL_nm_41_13_1[[i]])){
    temp <- redAL_nm_41_13_1[[i]]$arcL[j] / redAL_nm_41_13_1[[i]]$note.durs[j]
    nAL <- append(nAL,temp)
  }
  redAL_nm_41_13_1[[i]]$nAL <- nAL
}

for (i in 1:length(redAL_nm_49_12_1)){
  nAL <- c()
  for (j in 1:nrow(redAL_nm_49_12_1[[i]])){
    temp <- redAL_nm_49_12_1[[i]]$arcL[j] / redAL_nm_49_12_1[[i]]$note.durs[j]
    nAL <- append(nAL,temp)
  }
  redAL_nm_49_12_1[[i]]$nAL <- nAL
}

for (i in 1:length(redAL_nm_53_2_1)){
  nAL <- c()
  for (j in 1:nrow(redAL_nm_53_2_1[[i]])){
    temp <- redAL_nm_53_2_1[[i]]$arcL[j] / redAL_nm_53_2_1[[i]]$note.durs[j]
    nAL <- append(nAL,temp)
  }
  redAL_nm_53_2_1[[i]]$nAL <- nAL
}

for (i in 1:length(redAL_nm_54_2_2)){
  nAL <- c()
  for (j in 1:nrow(redAL_nm_54_2_2[[i]])){
    temp <- redAL_nm_54_2_2[[i]]$arcL[j] / redAL_nm_54_2_2[[i]]$note.durs[j]
    nAL <- append(nAL,temp)
  }
  redAL_nm_54_2_2[[i]]$nAL <- nAL
}

for (i in 1:length(redAL_nm_55_4_3)){
  nAL <- c()
  for (j in 1:nrow(redAL_nm_55_4_3[[i]])){
    temp <- redAL_nm_55_4_3[[i]]$arcL[j] / redAL_nm_55_4_3[[i]]$note.durs[j]
    nAL <- append(nAL,temp)
  }
  redAL_nm_55_4_3[[i]]$nAL <- nAL
}

for (i in 1:length(redAL_nm_55_5_2)){
  nAL <- c()
  for (j in 1:nrow(redAL_nm_55_5_2[[i]])){
    temp <- redAL_nm_55_5_2[[i]]$arcL[j] / redAL_nm_55_5_2[[i]]$note.durs[j]
    nAL <- append(nAL,temp)
  }
  redAL_nm_55_5_2[[i]]$nAL <- nAL
}

for (i in 1:length(redAL_nm_55_6_1)){
  nAL <- c()
  for (j in 1:nrow(redAL_nm_55_6_1[[i]])){
    temp <- redAL_nm_55_6_1[[i]]$arcL[j] / redAL_nm_55_6_1[[i]]$note.durs[j]
    nAL <- append(nAL,temp)
  }
  redAL_nm_55_6_1[[i]]$nAL <- nAL
}

for (i in 1:length(redAL_nm_56_1_3)){
  nAL <- c()
  for (j in 1:nrow(redAL_nm_56_1_3[[i]])){
    temp <- redAL_nm_56_1_3[[i]]$arcL[j] / redAL_nm_56_1_3[[i]]$note.durs[j]
    nAL <- append(nAL,temp)
  }
  redAL_nm_56_1_3[[i]]$nAL <- nAL
}

for (i in 1:length(redAL_nm_56_3_2)){
  nAL <- c()
  for (j in 1:nrow(redAL_nm_56_3_2[[i]])){
    temp <- redAL_nm_56_3_2[[i]]$arcL[j] / redAL_nm_56_3_2[[i]]$note.durs[j]
    nAL <- append(nAL,temp)
  }
  redAL_nm_56_3_2[[i]]$nAL <- nAL
}

for (i in 1:length(redAL_nm_63_1_1)){
  nAL <- c()
  for (j in 1:nrow(redAL_nm_63_1_1[[i]])){
    temp <- redAL_nm_63_1_1[[i]]$arcL[j] / redAL_nm_63_1_1[[i]]$note.durs[j]
    nAL <- append(nAL,temp)
  }
  redAL_nm_63_1_1[[i]]$nAL <- nAL
}

# Average normalized arc length per song for each individual
mean_nAL41 <- lapply(redAL_nm_41_13_1,function(x){ summarize(x,avg.AL = mean(nAL))}) %>%
  unlist(use.names = FALSE)
mean_curve_41_13_1 <- data.frame(ID = rep("41-13-1",length(mean_nAL41)), 
                                 Timepoint = rep(c("pre","post"),
                                                 times = c(lengths[1:2,3])),
                                 Treatment = rep("sham",times = length(redAL_nm_41_13_1)),
                                 mean_nAL = mean_nAL41)

mean_nAL49 <- lapply(redAL_nm_49_12_1,function(x){ summarize(x,avg.AL = mean(nAL))}) %>%
  unlist(use.names = FALSE)
mean_curve_49_12_1 <- data.frame(ID = rep("49-12-1",length(mean_nAL49)), 
                                 Timepoint = rep(c("pre","post"),
                                                 times = c(lengths[3:4,3])),
                                 Treatment = rep("ablation",times = length(redAL_nm_49_12_1)),
                                 mean_nAL = mean_nAL49)

mean_nAL53 <- lapply(redAL_nm_53_2_1,function(x){ summarize(x,avg.AL = mean(nAL))}) %>%
  unlist(use.names = FALSE)
mean_curve_53_2_1 <- data.frame(ID = rep("53-2-1",length(mean_nAL53)), 
                                Timepoint = rep(c("pre","post"),
                                                times = c(lengths[5:6,3])),
                                Treatment = rep("ablation",times = length(redAL_nm_53_2_1)),
                                mean_nAL = mean_nAL53)

mean_nAL54 <- lapply(redAL_nm_54_2_2,function(x){ summarize(x,avg.AL = mean(nAL))}) %>%
  unlist(use.names = FALSE)
mean_curve_54_2_2 <- data.frame(ID = rep("54-2-2",length(mean_nAL54)), 
                                Timepoint = rep(c("pre","post"),
                                                times = c(lengths[7:8,3])),
                                Treatment = rep("sham",times = length(redAL_nm_54_2_2)),
                                mean_nAL = mean_nAL54)

mean_nAL554 <- lapply(redAL_nm_55_4_3,function(x){ summarize(x,avg.AL = mean(nAL))}) %>%
  unlist(use.names = FALSE)
mean_curve_55_4_3 <- data.frame(ID = rep("55-4-3",length(mean_nAL554)), 
                                Timepoint = rep(c("pre","post"),
                                                times = c(lengths[9:10,3])),
                                Treatment = rep("ablation",times = length(redAL_nm_55_4_3)),
                                mean_nAL = mean_nAL554)

mean_nAL555 <- lapply(redAL_nm_55_5_2,function(x){ summarize(x,avg.AL = mean(nAL))}) %>%
  unlist(use.names = FALSE)
mean_curve_55_5_2 <- data.frame(ID = rep("55-5-2",length(mean_nAL555)), 
                                Timepoint = rep(c("pre","post"),
                                                times = c(lengths[11:12,3])),
                                Treatment = rep("ablation",times = length(redAL_nm_55_5_2)),
                                mean_nAL = mean_nAL555)

mean_nAL556 <- lapply(redAL_nm_55_6_1,function(x){ summarize(x,avg.AL = mean(nAL))}) %>%
  unlist(use.names = FALSE)
mean_curve_55_6_1 <- data.frame(ID = rep("55-6-1",length(mean_nAL556)), 
                                Timepoint = rep(c("pre","post"),
                                                times = c(lengths[13:14,3])),
                                Treatment = rep("sham",times = length(redAL_nm_55_6_1)),
                                mean_nAL = mean_nAL556)

mean_nAL561 <- lapply(redAL_nm_56_1_3,function(x){ summarize(x,avg.AL = mean(nAL))}) %>%
  unlist(use.names = FALSE)
mean_curve_56_1_3 <- data.frame(ID = rep("56-1-3",length(mean_nAL561)), 
                                Timepoint = rep(c("pre","post"),
                                                times = c(lengths[15:16,3])),
                                Treatment = rep("sham",times = length(redAL_nm_56_1_3)),
                                mean_nAL = mean_nAL561)

mean_nAL563 <- lapply(redAL_nm_56_3_2,function(x){ summarize(x,avg.AL = mean(nAL))}) %>%
  unlist(use.names = FALSE)
mean_curve_56_3_2 <- data.frame(ID = rep("56-3-2",length(mean_nAL563)), 
                                Timepoint = rep(c("pre","post"),
                                                times = c(lengths[17:18,3])),
                                Treatment = rep("ablation",times = length(redAL_nm_56_3_2)),
                                mean_nAL = mean_nAL563)

mean_nAL63 <- lapply(redAL_nm_63_1_1,function(x){ summarize(x,avg.AL = mean(nAL))}) %>%
  unlist(use.names = FALSE)
mean_curve_63_1_1 <- data.frame(ID = rep("63-1-1",length(mean_nAL63)), 
                                Timepoint = rep(c("pre","post"),
                                                times = c(lengths[19:20,3])),
                                Treatment = rep("sham",times = length(redAL_nm_63_1_1)),
                                mean_nAL = mean_nAL63)

# all individuals into one dataframe
curve_msrs <- rbind(mean_curve_41_13_1,mean_curve_49_12_1,mean_curve_53_2_1,
                    mean_curve_54_2_2,mean_curve_55_4_3,mean_curve_55_5_2,
                    mean_curve_55_6_1,mean_curve_56_1_3,mean_curve_56_3_2,
                    mean_curve_63_1_1)
curve_msrs$Timepoint <- factor(curve_msrs$Timepoint, levels = c("pre","post"))
curve_msrs$Treatment <- factor(curve_msrs$Treatment, levels = c("sham","ablation"))

######################################################################################
# Modeling effect of treatment and timepoint on average, normalized arc length
library(lme4)
library(lmtest)
library(MuMIn)
library(sjPlot)
library(effects)
library(AICcmodavg)

mod1 <- lmer(data=curve_msrs,mean_nAL ~ 1 + (1|ID),REML=FALSE) # only subject as random effect
mod2 <- lmer(data=curve_msrs,mean_nAL ~ Timepoint + (1|ID),REML=FALSE) # timepoint as fixed effect
mod3 <- lmer(data=curve_msrs,mean_nAL ~ Treatment + Timepoint + (1|ID),REML=FALSE) # treatment and timepoint as fixed effects
mod4 <- lmer(data=curve_msrs,mean_nAL ~ Treatment * Timepoint + (1|ID),REML=FALSE) # treatment, timepoint and interaction
mod5 <- lmer(data=curve_msrs,mean_nAL ~ Treatment * Timepoint + (Timepoint | ID), REML=FALSE) # treatment, timepoint, interaction = fixed effects. random slopes model - individual responds diff to timepoint

lrtest(mod1,mod2) # mod 2 better than mod 1
lrtest(mod2,mod3) # mod 3 not better than mod 2
lrtest(mod2,mod4) # mod 4 better than mod 2
lrtest(mod4,mod5) # mod 5 better than mod 4
summary(mod5) # mean_AL = 19124.1 + 2949.6(ablation) + 390(post-surgery) -7495.2(ablation:post)
r.squaredGLMM(mod5) # R2m - marginal - var explained by fixed factors, R2c- full model
# R2m = 0.3581598, R2c = 0.6131519

# which models to calculate AIC for
models <- list(mod1,mod2,mod3,mod4,mod5)
# model names
model.names <- c("null","tpt","trt.tpt","trt.tpt.int","trt.tpt.int.rslope")
# AIC for each model
aictab(cand.set=models,modnames = model.names)
# K: The number of parameters in the model.
# AICc: The AIC value of the model. The lowercase ‘c’ indicates that the AIC has been calculated from the AIC corrected for small sample sizes.
# Delta_AICc: The difference between the AIC of the best model compared to the current model being compared.
# AICcWt: The proportion of the total predictive power that can be found in the model.
# Cum.Wt: The cumulative sum of the AIC weights.
# LL: The log-likelihood of the model. This tells us how likely the model is, given the data we used.

# plot effect sizes for best model
sjPlot::plot_model(mod5,
                   show.values=TRUE,show.p=TRUE,colors="gs") + theme_classic()
sjPlot::tab_model(mod5)
effects_intx <- effects::effect(term="Treatment:Timepoint",mod=mod5)
summary(effects_intx)
plot_model(mod5,type="int",show.p=TRUE,colors="gs") + theme_classic() # predicted values of mean_nAL

######################################################################################################
# 55-4-3 is an outlier because it has a reduced PCA and TA in addition to CT. 
# Here, I redo the modelling with this sample removed to see if that changes the results

red_curve_msrs <- curve_msrs %>% filter(ID != "55-4-3")

rmod1 <- lmer(data=red_curve_msrs,mean_nAL ~ 1 + (1|ID),REML=FALSE) # only subject as random effect
rmod2 <- lmer(data=red_curve_msrs,mean_nAL ~ Timepoint + (1|ID),REML=FALSE) # timepoint as fixed effect
rmod3 <- lmer(data=red_curve_msrs,mean_nAL ~ Treatment + Timepoint + (1|ID),REML=FALSE) 
# trt + timepoint as fixed
rmod4 <- lmer(data=red_curve_msrs,mean_nAL ~ Treatment * Timepoint + (1|ID),REML=FALSE) 
rmod5 <- lmer(data=red_curve_msrs,mean_nAL ~ Treatment * Timepoint + (Timepoint | ID), REML=FALSE)

# adding trt,time interaction
lrtest(rmod1,rmod2) # mod 2 better than mod 1
lrtest(rmod2,rmod3) # mod 3 not better than mod 2
lrtest(rmod2,rmod4) # mod 4 better than mod 2
lrtest(rmod4,rmod5) # mod 5 better than mod 4
summary(rmod5) # mean_AL = 19124.1 + 2949.6(ablation) + 390(post-surgery) -7495.2(ablation:post)
r.squaredGLMM(rmod5) # R2m - marginal - var explained by fixed factors, R2c- full model
# R2m = 0.2890925, R2c = 0.5614361

# which models to calculate AIC
rmodels <- list(rmod1,rmod2,rmod3,rmod4,rmod5)
# model names
rmodel.names <- c("null","tpt","trt.tpt","trt.tpt.int","trt.tpt.int.rslope")
# AIC for each model
aictab(cand.set=rmodels,modnames = rmodel.names)
# K: The number of parameters in the model.
# AICc: The AIC value of the model. The lowercase ‘c’ indicates that the AIC has been calculated from the AIC corrected for small sample sizes.
# Delta_AICc: The difference between the AIC of the best model compared to the current model being compared.
# AICcWt: The proportion of the total predictive power that can be found in the model.
# Cum.Wt: The cumulative sum of the AIC weights.
# LL: The log-likelihood of the model. This tells us how likely the model is, given the data we used.

# plot effect sizes
sjPlot::plot_model(rmod5,
                   show.values=TRUE,show.p=TRUE,colors="gs") + theme_classic()
sjPlot::tab_model(rmod5)
effects_intx <- effects::effect(term="Treatment:Timepoint",mod=rmod5)
summary(effects_intx)
plot_model(rmod5,type="int",show.p=TRUE,colors="gs") + theme_classic() # predicted values of mean_nAL

######################################################################################################
# Pull in muscle measurement data and calculate volume estimates for muscle. 
# Measurements were done in imageJ (v. 2.3.0/1.53q) using the polygon tool.

library(tidyverse)
library(dplyr)
library(patchwork)
library(cowplot)
library(ggpubr)

file <- "muscle_ablation_muscle_measurements.csv"
measurements <- read_csv(file,col_names=TRUE)

# LTA = lateral thyroarytenoid muscle, MTA = medial thyroarytenoid muscle, TA = thyroarytenoid muscle
measurements <- measurements %>% mutate("TA_area" = `LTA area` + `MTA area`)
# calculate volume by taking the sum of each area measurement (in mm squared) multiplied by the distance between each section (thickness) = 0.02 mm
volumes <- measurements %>% group_by(ID) %>% summarize("CT_Volume" = sum(`CT area (mm^2)`*0.02,na.rm=TRUE),
                                                       "LCA_Volume" = sum(`LCA area`*0.02,na.rm=TRUE),
                                                       "LTA_Volume" = sum(`LTA area`*0.02,na.rm=TRUE),
                                                       "MTA_Volume" = sum(`MTA area`*0.02,na.rm=TRUE),
                                                       "TA_Volume" = sum(TA_area*0.02,na.rm=TRUE))
# add in information about which sample got which treatment
volumes$Treatment <- c("sham","ablation","ablation","sham","ablation","ablation","sham","sham","ablation","sham")
save(volumes,file="muscle_volumes.RData")

#############################################################################################
# Load in muscle volume data and see whether muscle volumes predict mean, normalized arc length

load("muscle_volumes.RData")
curve_diffs <- curve_msrs %>% group_by(ID, Timepoint) %>% 
  summarize('mean_mean_nAL' = mean(mean_nAL))
mean_nAL_diff <- curve_diffs[curve_diffs$Timepoint == "post",3] - curve_diffs[curve_diffs$Timepoint == "pre",3]

nALvMuscle <- cbind("ID" = volumes$ID, "Treatment" = volumes$Treatment,"CT_Volume" = volumes$CT_Volume, 
                    "LCA_Volume" = volumes$LCA_Volume, "TA_Volume" = volumes$TA_Volume,
                    "mnAL_diff" = mean_nAL_diff)
names(nALvMuscle) <- c("ID","Trt","CT","LCA","TA","mnALdiff")
nALvMuscle$Trt <- factor(nALvMuscle$Trt,levels = c("sham","ablation"))

# modeling to see if CT volume predicts difference in curve measures (pre to post surg)
# create a model for every combination of predictors
vmod1 <- lm(data=nALvMuscle,mnALdiff ~ 1) # intercept only model
summary(vmod1) # intercept significant p<0.01
vmod2 <- lm(data=nALvMuscle,mnALdiff ~ CT) # cricothyroid volume
summary(vmod2) # CT significant p<0.01
vmod3 <- lm(data=nALvMuscle,mnALdiff ~ LCA) # lateral cricoarytenoid volume
summary(vmod3) # LCA NOT significant p=0.0544
vmod4 <- lm(data=nALvMuscle,mnALdiff ~ TA) # thyroarytenoid volume
summary(vmod4) # TA NOT significant p=0.299
vmod5 <- lm(data=nALvMuscle,mnALdiff ~ CT + LCA) # cricothyroid and lateral cricoarytenoid volumes
summary(vmod5) # CT significant p < 0.05, LCA NOT significant p=0.381
vmod6 <- lm(data=nALvMuscle,mnALdiff ~ CT + TA) # cricothyroid and thyroarytenoid volumes
summary(vmod6) # CT significant p = 0.018, TA NOT significant p = 0.944
vmod7 <- lm(data=nALvMuscle,mnALdiff ~ CT + LCA + TA) # all three muscles
summary(vmod7) # CT significant p < 0.05, TA and LCA NOT significant (p=0.13 and 0.2 respectively)

# compare modules using a likelihood ratio test
lrtest(vmod1,vmod2) # Adding CT improvement over null model (~1)
lrtest(vmod2,vmod5) # Adding LCA does not improve model
lrtest(vmod6,vmod2) # Adding TA does not improve model
lrtest(vmod7,vmod2) # Adding both LCA and TA does not improve model
# best model is one that only includes CT as a predictor 

# which models to calculate AIC
vmodels <- list(vmod1,vmod2,vmod3,vmod4,vmod5,vmod6,vmod7)
# model names
vmodel.names <- c("null","CT","LCA","TA","CT.LCA","CT.TA","CT.LCA.TA")
# AIC for each model
aictab(cand.set=vmodels,modnames = vmodel.names)
# K: The number of parameters in the model.
# AICc: The AIC value of the model. The lowercase ‘c’ indicates that the AIC has been calculated from the AIC corrected for small sample sizes.
# Delta_AICc: The difference between the AIC of the best model compared to the current model being compared.
# AICcWt: The proportion of the total predictive power that can be found in the model.
# Cum.Wt: The cumulative sum of the AIC weights.
# LL: The log-likelihood of the model. This tells us how likely the model is, given the data we used.

# Boxplots of relationship between muscle volumes and mean, normalized arc length
library(ggpmisc)

fa <- ggplot(nALvMuscle, aes(x=CT,y=mnALdiff)) + geom_point(aes(col=Trt,pch=Trt),show.legend = FALSE,size=3) + 
  # geom_text(aes(label = ID),vjust=-0.1,hjust=-0.1) +
  scale_color_manual(values=c("#1F78B4","#33A02C")) + 
  xlab(expression("CT Volume (mm"^2*")")) + 
  ylab(expression(paste(Delta*" Frequency Modulation"))) +
  geom_hline(yintercept = 0, linetype = 2,col="gray") +
  stat_poly_line(col="black",alpha=0.3) + 
  stat_poly_eq(mapping=use_label(labels=c("R2","P"),aes(size=7))) + theme_classic() +
  theme(axis.line = element_line(linewidth = 1),
        axis.text = element_text(size=14),
        axis.title = element_text(size=18))
fa

fb <- ggplot(nALvMuscle, aes(x=LCA,y=mnALdiff)) + geom_point(aes(col=Trt),show.legend = FALSE) +
  scale_color_manual(values=c("#1F78B4","#33A02C")) + 
  xlab(expression("LCA Volume (mm"^2*")")) + ylab("post-pre mnAL") +
  stat_poly_line(col="black",alpha=0.3) + stat_poly_eq() + theme_classic()

fc <- ggplot(nALvMuscle, aes(x=TA,y=mnALdiff)) + geom_point(aes(col=Trt),show.legend = FALSE) +
  scale_color_manual(values=c("#1F78B4","#33A02C")) + 
  xlab(expression("TA Volume (mm"^2*")")) + ylab("post-pre mnAL") +
  stat_poly_line(col="black",alpha=0.3) + stat_poly_eq() + theme_classic()

right_side <- cowplot::plot_grid(fb,fc,align="v",ncol=1,label_size = 12)
cowplot::plot_grid(fa,right_side,align="hv",ncol=2,label_size = 12)

##################################################################################
# Relationship between muscle volume and mean, normalized arc length excluding 55-4-3
# remove data from 55-4-3 individual
red_nALvM <- nALvMuscle %>% filter(ID != "55-4-3")

ra <- ggplot(red_nALvM, aes(x=CT,y=mnALdiff)) + geom_point(aes(col=Trt),show.legend = FALSE) + 
  # geom_text(aes(label = ID),vjust=-0.1,hjust=-0.1) +
  scale_color_manual(values=c("#1F78B4","#33A02C")) + 
  xlab(expression("CT Volume (mm"^2*")")) +
  ylab("post-pre mnAL") +
  stat_poly_line(col="black",alpha=0.3) + stat_poly_eq() + theme_classic()

rb <- ggplot(red_nALvM, aes(x=LCA,y=mnALdiff)) + geom_point(aes(col=Trt),show.legend = FALSE) +
  scale_color_manual(values=c("#1F78B4","#33A02C")) + 
  xlab(expression("LCA Volume (mm"^2*")")) +
  ylab("post-pre mnAL") +
  stat_poly_line(col="black",alpha=0.3) + stat_poly_eq() + theme_classic()

rc <- ggplot(red_nALvM, aes(x=TA,y=mnALdiff)) + geom_point(aes(col=Trt),show.legend = FALSE) +
  scale_color_manual(values=c("#1F78B4","#33A02C")) + 
  xlab(expression("TA Volume (mm"^2*")")) +
  ylab("post-pre mnAL") +
  stat_poly_line(col="black",alpha=0.3) + stat_poly_eq() + theme_classic()

rright_side <- cowplot::plot_grid(rb,rc,align="v",ncol=1,label_size = 12)
cowplot::plot_grid(ra,rright_side,align="hv",ncol=2,label_size = 12)