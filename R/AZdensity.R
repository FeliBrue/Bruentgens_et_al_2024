# Title: Nested active zone density analysis
#
# Description:
# This script performs a nested active zone density analysis. 
# It calculates the area of active zones within each nested bouton.
#
# Author: Daniel Parthier / Felicitas Brüntgens
# Date: 05.10.2023
# Last modified: 19.04.2024
# Version: 1.0


# Load libraries ####
library(lme4)
library(ggplot2)
library(ggbeeswarm)
library(data.table)
library(emmeans)


# Loading the data ####
AZdensity <- fread(file = "Data/AZdensity.csv", sep = ";", stringsAsFactors = T)

# Data visualization ####
ggplot(data = AZdensity, mapping = aes(y=AZdensity, x = Group, colour=Group))+
  geom_line(data = AZdensity[,.(AZdensity=mean(AZdensity)),by=.(Animal,Group)],
            inherit.aes = F, aes(y=AZdensity, x = Group, group = Animal, colour=Group),
            linewidth = 1, alpha = 0.4)+
  geom_beeswarm(cex = 2.5)+
  scale_colour_manual(values = c("WT" = "#45abff", 
                                 "WT + FSK" = "#6AB560", 
                                 "TKO" = "#e61c47", 
                                 "TKO + FSK" = "#FECC66"))+
  scale_y_continuous(name = expression("AZ Density " (AZ/µm^3)), 
                     limits = c(0,25))+
  labs(x="")+
  theme_classic()+
  theme(legend.position = "none")


# QQ-plot to assess residuals ####
FullModelLmer <- lmer(AZdensity ~ Genotype * Treatment + (1|Animal) + (1|SampleID),
                      data = AZdensity)

FullModel <- glmer(AZdensity ~ Genotype * Treatment + (1|Animal) + (1|SampleID),
                   data = AZdensity, family = Gamma(link = "log"), 
                   control = glmerControl(nAGQ0initStep = 0))

AZdensity[,`:=`(`Scaled Gaussian Residuals`=scale(residuals(FullModelLmer)),
             `Scaled Gamma Residuals`=scale(residuals(FullModel))),]
ResidualDT <- melt(AZdensity, measure.vars = c("Scaled Gaussian Residuals", "Scaled Gamma Residuals"),
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
FullModel <- glmer(AZdensity ~ Genotype * Treatment + (1|Animal) + (1|SampleID),
                   data = AZdensity, 
                   family = Gamma(link = "log"), 
                   control = glmerControl(nAGQ0initStep = 0))

NullModel <- glmer(AZdensity ~ 1 + (1|Animal) + (1|SampleID),
                   data = AZdensity, 
                   family = Gamma(link = "log"),
                   control = glmerControl(nAGQ0initStep = 0))

(anova(FullModel, NullModel))


## Treatment, Genotype and Interaction models ####
TreatmentNullModel <- glmer(formula = AZdensity ~ Genotype + (1|Animal) + (1|SampleID),
                            data = AZdensity,
                            family = Gamma(link = "log"),
                            control = glmerControl(nAGQ0initStep = 0))

(anova(FullModel, TreatmentNullModel))


GenotypeNullModel <- glmer(formula = AZdensity ~ Treatment + (1|Animal) + (1|SampleID), 
                           data = AZdensity, 
                           family = Gamma(link = "log"), 
                           control = glmerControl(nAGQ0initStep = 0))

(anova(FullModel, GenotypeNullModel))


InteractionNullModel <- glmer(formula = AZdensity ~ Treatment + Genotype + (1|Animal) + (1|SampleID),
                              data = AZdensity,
                              family = Gamma(link = "log"), 
                              control = glmerControl(nAGQ0initStep = 0))

(anova(FullModel, InteractionNullModel))


# Post-hoc analysis with Estimated Marginal Means ####
MarginalEffects <- emmeans(object = FullModel, specs =  ~ Genotype * Treatment, 
                           type="response", tran = "log")
MarginalEffects
ContrastMatrix <- diag(4)

GroupGenotype <- MarginalEffects@grid$Genotype
GroupTreatment <- MarginalEffects@grid$Treatment

### Contrast and adjustment for multiple testing (Benjamini-Hochberg) ####
MarginalEffectTest <- contrast(object = MarginalEffects, type = "response",
                                     method = list("WT / (WT + FSK)" = ContrastMatrix[GroupGenotype=="WT"&GroupTreatment=="Control",] - ContrastMatrix[GroupGenotype=="WT"&GroupTreatment=="FSK",],
                                                   "WT / TKO" = ContrastMatrix[GroupGenotype=="WT"&GroupTreatment=="Control",] - ContrastMatrix[GroupGenotype=="TKO"&GroupTreatment=="Control",],
                                                   "TKO / (TKO + FSK)" = ContrastMatrix[GroupGenotype=="TKO"&GroupTreatment=="Control",] - ContrastMatrix[GroupGenotype=="TKO"&GroupTreatment=="FSK",],
                                                   "(WT + FSK) / (TKO + FSK)" = ContrastMatrix[GroupGenotype=="WT"&GroupTreatment=="FSK",] - ContrastMatrix[GroupGenotype=="TKO"&GroupTreatment=="FSK",],
                                                   "TKO / (WT + FSK)" = ContrastMatrix[GroupGenotype=="TKO"&GroupTreatment=="Control",] - ContrastMatrix[GroupGenotype=="WT"&GroupTreatment=="FSK",]), adjust="BH")
print(MarginalEffectTest)

### Plot marginal effects ####
(MarginalPlot <- plot(MarginalEffectTest, colors = "grey10") +
   labs(title = "Estimated AZ Density Ratio")+
   ylab(label = "")+
   geom_vline(xintercept = 1, linetype="dashed")+
   xlab(label = "Ratio")+
   theme_classic())
