# Title: SV nested nearest neighbor analysis
#
# Description:
# This script performs a nested nearest neighbor analysis of synaptic vesicles.
#
# Author: Daniel Parthier
# Date: 19.04.2024
# Last modified: 19.04.2024
# Version: 1.0

# Load libraries ####
library(readxl)
library(data.table)
library(ggplot2)
library(ggbeeswarm)
library(lme4)

DataSet <- fread(input = "Data/MNND.csv", sep = ";", stringsAsFactors = TRUE)
DataSet[,Genotype:=factor(x = Genotype, levels = c("WT", "TKO")),]

# Data visualization ####
ggplot(data = DataSet[,.(Distance=mean(Distance)),by=.(BoutonID, SampleID, Genotype, Animal)], mapping = aes(y=Distance, x = Genotype, colour=Genotype))+
  geom_beeswarm(cex = 2.5)+
  scale_colour_manual(values = c("WT" = "#45abff",
                                 "TKO" = "#e61c47"))+
  scale_y_continuous(limits = c(0,125))+
  
  xlab("Genotype")+
  ylab("Mean Nearest Neighbor Distance (nm)")+
  theme_classic()+
  theme(legend.position = "none")

# QQ-plot to assess residuals ####
FullModelLmer <- lmer(formula = (Distance) ~ Genotype + (1|Animal) + (1|SampleID) + (1|BoutonID),
                      data = DataSet, REML = FALSE)
FullModel <- glmer(formula = Distance ~ Genotype + (1|Animal) + (1|SampleID) + (1|BoutonID),
                   data = DataSet,
                   family = Gamma(link = "log"),
                   control = glmerControl(nAGQ0initStep = 0, optimizer = "bobyqa"))

DataSet[,`:=`(Gaussian=scale(residuals(FullModelLmer)),Gamma=scale(residuals(FullModel))),]
ResidualDT <- melt(DataSet, measure.vars = c("Gaussian", "Gamma"), id.vars = c("BoutonID", "Animal", "Genotype"), variable.name = "Distribution", value.name = "Residuals")

ggplot(data = ResidualDT, aes(sample = Residuals))+
  stat_qq()+
  stat_qq_line()+
  facet_wrap(~Distribution)+
  scale_x_continuous(name = "Theoretical quantiles")+
  scale_y_continuous(name = "Scaled Residuals")+
  theme_classic()

# Build general linear models ####
NullModel <- glmer(formula = Distance ~ 1 + (1|Animal) + (1|SampleID) + (1|BoutonID),
                   data = DataSet,
                   family = Gamma(link = "log"),
                   control = glmerControl(nAGQ0initStep = 0, optimizer = "bobyqa"))
(anova(NullModel, FullModel))
