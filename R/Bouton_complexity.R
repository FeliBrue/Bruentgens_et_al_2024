# Title: Nested bouton complexity analysis
#
# Description:
# This script performs a nested bouton complexity analysis. 
# It calculates the bouton complexity within each nested animal.
#
# Author: Daniel Parthier / Felicitas Brüntgens
# Date: 17.04.2024
# Last modified: 19.04.2024
# Version: 1.0

# Loading libraries ####
library(data.table)
library(readxl)
library(ggplot2)
library(ggbeeswarm)
library(lme4)


#Loading the data ####
Complexity <- fread(file = "Data/Bouton_complexity.csv", sep = ";", stringsAsFactors = T)


# Data visualization ####
ggplot(data = Complexity, mapping = aes(y= complexity, x = Group, colour=Group))+
  geom_line(data = Complexity[,.(complexity=mean(complexity)),by=.(Animal,Group)],
             inherit.aes = F, aes(y=complexity, x = Group, group = Animal, colour=Group),
             linewidth = 1, alpha = 0.4)+
  geom_beeswarm(cex = 2.5)+
  scale_colour_manual(values = c("WT_ctrl" = "#45abff", 
                                 "WT_fsk" = "#6AB560", 
                                 "TKO_ctrl" = "#e61c47", 
                                 "TKO_fsk" = "#FECC66"))+
  scale_y_continuous(name = expression("Complexity"), 
                     limits = c(0,8))+
  labs(x="")+
  theme_classic()+
  theme(legend.position = "none")


# QQ-plot to assess residuals ####
FullModelLmer <- lmer(formula = complexity ~ Genotype * treatment + (1|Animal) + (1|SampleID) + (1|BoutonID),
                      data = Complexity)

FullModel <- glmer(formula = complexity ~ Genotype * treatment + (1|Animal) + (1|SampleID) + (1|BoutonID),
                   data = Complexity,
                   family = Gamma(link = "log"),
                   control = glmerControl(nAGQ0initStep = 0))

Complexity[,`:=`(`Scaled Gaussian Residuals`=scale(residuals(FullModelLmer)),`Scaled Gamma Residuals`=scale(residuals(FullModel))),]
ResidualDT <- melt(Complexity, measure.vars = c("Scaled Gaussian Residuals", "Scaled Gamma Residuals"), variable.name = "Distribution", value.name = "Residuals")

ggplot(data = ResidualDT, aes(sample = Residuals))+
  stat_qq()+
  stat_qq_line(linetype="dashed")+
  facet_wrap(~Distribution)+
  scale_x_continuous(name = "Theoretical quantiles")+
  scale_y_continuous(name = "Scaled Residuals")+
  theme_classic()


# Build general linear models ####
FullModel <- glmer(formula = complexity ~ Genotype * treatment + (1|Animal) + (1|SampleID) + (1|BoutonID),
                    data = Complexity,
                    family = Gamma(link = "log"),
                    control = glmerControl(nAGQ0initStep = 0))

NullModel <- glmer(formula = complexity ~ 1 + (1|Animal) + (1|SampleID) + (1|BoutonID),
                   data = Complexity,
                   family = Gamma(link = "log"),
                   glmerControl(nAGQ0initStep = 0))
(anova(FullModel, NullModel))