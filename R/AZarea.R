# Title: Nested active zone area analysis
#
# Description:
# This script performs a nested active zone area analysis. 
# It calculates the area of active zones within each nested bouton.
#
# Author: Daniel Parthier / Felicitas Brüntgens
# Date: 05.10.2023
# Last modified: 10.04.2024
# Version: 1.0

# an Daniel:
# Sollten wir hier konsequenterweise nach FullModel vs. NullModel aufhören?

# Load libraries ####
library(lme4)
library(ggplot2)
library(ggbeeswarm)
library(data.table)
library(emmeans)

# Loading the data ####
AZarea <- fread(file = "Data/AZarea_DataSet.csv", sep = ";", stringsAsFactors = T)

# Data visualization ####
ggplot(data = AZarea, mapping = aes(y=Active.zone.area, x = Condition, colour=Condition))+
  geom_line(data = AZarea[,.(Active.zone.area=mean(Active.zone.area)),by=.(Condition, Animal, Genotype)],
            inherit.aes = F, aes(y=Active.zone.area, x = Condition, group=Animal, colour=Genotype),
            linewidth = 1, alpha=0.4)+
  geom_beeswarm(cex = 2.5)+
  scale_colour_manual(values = c("WT" = "#45abff",
                                 "WT + FSK" = "#6AB560",
                                 "TKO" = "#e61c47",
                                 "TKO + FSK" = "#FECC66"))+
  scale_y_continuous(name = expression(paste("AZ Area ",  "n", m^2)))+
  labs(x="")+
  theme_classic()

# QQ-plot to assess residuals ####
FullModelLmer <- lmer(Active.zone.area ~ Genotype * Treatment + (1|BoutonID) + (1|Animal) + (1|SampleID),
                      data = AZarea)

FullModel <- glmer(Active.zone.area ~ Genotype * Treatment + (1|BoutonID) + (1|Animal) + (1|SampleID),
                   data = AZarea, family = Gamma(link = "log"), 
                   control = glmerControl(nAGQ0initStep = 0))

AZarea[,`:=`(`Scaled Gaussian Residuals`=scale(residuals(FullModelLmer)),
             `Scaled Gamma Residuals`=scale(residuals(FullModel))),]
ResidualDT <- melt(AZarea, measure.vars = c("Scaled Gaussian Residuals", "Scaled Gamma Residuals"),
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
FullModel <- glmer(Active.zone.area ~ Genotype * Treatment + (1|BoutonID) + (1|Animal) + (1|SampleID),
                   data = AZarea, family = Gamma(link = "log"), 
                   control = glmerControl(nAGQ0initStep = 0))

NullModel <- glmer(Active.zone.area ~ 1 + (1|BoutonID) + (1|Animal) + (1|SampleID),
                   data = AZarea, family = Gamma(link = "log"),
                   control = glmerControl(nAGQ0initStep = 0))
(anova(FullModel, NullModel))