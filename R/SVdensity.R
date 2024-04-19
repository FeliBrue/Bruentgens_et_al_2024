# Title: Nested SV density analysis
#
# Description:
# This script performs a nested SV density analysis.
# It calculates the density of synaptic vesicles within each nested bouton.
#
# Author: Daniel Parthier / Felicitas Brüntgens
# Date: 03.04.2024
# Last modified: 09.04.2024
# Version: 1.0


# Load libraries ####
library(data.table)
library(ggplot2)
library(ggbeeswarm)
library(lme4)

# Loading the data ####
SVDensity <- fread(file="Data/SVDensity_DataSet.csv", sep=";", stringsAsFactors = T)

# Data visualization ####
ggplot(data = SVDensity, mapping = aes(y=SVdensity, x = Group, colour=Group))+
  geom_beeswarm(cex = 2.5)+
  scale_colour_manual(values = c("WT" = "#45abff", "TKO" = "#e61c47"))+
  scale_y_continuous(name = expression("SV density " (SV/µm^3)),
                     limits = c(0, NA))+
  labs(x="")+
  theme_classic()

# QQ-plot to assess residuals ####
FullModelLmer <- lmer(formula = SVdensity ~ Genotype + (1|Animal),
                      data = SVDensity)
FullModel <- glmer(formula = SVdensity ~ Genotype + (1|Animal),
                   data = SVDensity,
                   family = Gamma(link = "log"),
                   control = glmerControl(nAGQ0initStep = 0))

SVDensity[,`:=`(`Scaled Gaussian Residuals`=scale(residuals(FullModelLmer)),`Scaled Gamma Residuals`=scale(residuals(FullModel))),]
ResidualDT <- melt(SVDensity, measure.vars = c("Scaled Gaussian Residuals", "Scaled Gamma Residuals"), variable.name = "Distribution", value.name = "Residuals")

ggplot(data = ResidualDT, aes(sample = Residuals))+
  stat_qq()+
  stat_qq_line(linetype="dashed")+
  facet_wrap(~Distribution)+
  scale_x_continuous(name = "Theoretical quantiles")+
  scale_y_continuous(name = "Scaled Residuals")+
  theme_classic()

# Build general linear models ####
FullModel <- glmer(formula = SVdensity ~ Genotype + (1|Animal),
                   data = SVDensity,
                   family = Gamma(link = "log"),
                   control = glmerControl(nAGQ0initStep = 0))

NullModel <- glmer(formula = SVdensity ~ 1 + (1|Animal),
                   data = SVDensity,
                   family = Gamma(link = "log"),
                   glmerControl(nAGQ0initStep = 0))
(AllFixedEffects <- anova(FullModel, NullModel))

