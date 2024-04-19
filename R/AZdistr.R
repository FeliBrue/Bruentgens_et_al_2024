# Title: Nested active zone distribution analysis
#
# Description:
# This script performs a nested active zone distribution analysis. 
# It calculates the 3D distance of active zones within each nested bouton.
#
# Author: Daniel Parthier / Felicitas Brüntgens
# Date: 09.04.2024
# Last modified: 19.04.2024
# Version: 1.0


# Loading libraries ####
library(data.table)
library(readxl)
library(ggplot2)
library(ggbeeswarm)
library(lme4)
library(report)


#Loading the data ####
AZdistribution <- fread(file = "Data/AZdistribution.csv", sep = ";", stringsAsFactors = T)


# Data visualization ####
ggplot(data = AZdistribution, mapping = aes(y=AZdist, x = Group, colour=Group))+
  geom_line(data = AZdistribution[,.(AZdist=mean(AZdist)),by=.(Animal,Group)],
           inherit.aes = F, aes(y=AZdist, x = Group, group = Animal, colour=Group),
           linewidth = 1, alpha = 0.4)+
  geom_beeswarm(cex = 2)+
  scale_colour_manual(values = c("TKO Control" = "#e61c47",
                                 "TKO FSK" = "#FECC66",
                                 "WT Control" = "#45abff", 
                                 "WT FSK" = "#6AB560"))+
  scale_y_continuous(name = expression("AZ 3D Distance (nm)"), 
                     limits = c(0, 2000))+
  labs(x="")+
  theme_classic()+
  theme(legend.position = "none")


# QQ-plot to assess residuals ####
FullModel <- glmer(AZdist ~ Treatment + Genotype + (1|Animal) + (1|SampleID) + (1|BoutonID),
                   data = AZdistribution,
                   family = Gamma(link = "log"),
                   control = glmerControl(nAGQ0initStep = 0))

summary(FullModel)


FullModelLmer <- lmer(AZdist ~ Treatment * Genotype + (1|Animal) + (1|SampleID) + (1|BoutonID),
                      data = AZdistribution,
                      lmerControl(optimizer = "bobyqa"), REML = F)


AZdistribution[,`:=`(`Scaled Gaussian Residuals`=scale(residuals(FullModelLmer)),
             `Scaled Gamma Residuals`=scale(residuals(FullModel))),]
ResidualDT <- melt(AZdistribution, measure.vars = c("Scaled Gaussian Residuals", "Scaled Gamma Residuals"),
                   id.vars = c("Animal", "Treatment","Genotype"),
                   variable.name = "Distribution",
                   value.name = "Residuals")

ggplot(data = ResidualDT, aes(sample = Residuals))+
  stat_qq()+
  stat_qq_line()+
  facet_wrap(~Distribution)+
  scale_x_continuous(name = "Theoretical quantiles")+
  scale_y_continuous(name = "Scaled Residuals")+
  theme_classic()



# Build general linear models ####
FullModel <- glmer(AZdist ~ Treatment + Genotype + (1|Animal) + (1|SampleID) + (1|BoutonID),
                   data = AZdistribution,
                   family = Gamma(link = "log"),
                   control = glmerControl(nAGQ0initStep = 0))

summary(FullModel)


NullModel <- glmer(AZdist ~ 1 + (1|Animal) + (1|SampleID) + (1|BoutonID),
                   data = AZdistribution,
                   family = Gamma(link = "log"),
                   control = glmerControl(nAGQ0initStep = 0))

(anova(NullModel, FullModel))
