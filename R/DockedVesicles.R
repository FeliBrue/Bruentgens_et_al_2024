# Title: Nested docked vesicle density analysis
#
# Description:
# This script performs a nested docked vesicle density analysis.
# It calculates the density of docked vesicles per volume
# within each nested active zone.
#
# Author: Daniel Parthier / Felicitas Brüntgens
# Date: 04.04.2024
# Last modified: 10.04.2024
# Version: 1.0

# Load libraries ####
library(data.table)
library(ggplot2)
library(ggbeeswarm)
library(lme4)
library(emmeans)

# Loading the data ####
DataSetVolume <- fread(file = "Data/DockedVesicles_DataSetVolume.csv", sep = ";", stringsAsFactors = T)

# Data visualization ####
## SV density per bouton volume ####
ggplot(data = DataSetVolume, mapping = aes(y=Density, x = Condition, colour=Condition))+
  geom_line(data = DataSetVolume[,.(Density=mean(Density)),by=.(Condition, Animal, Genotype)],
            inherit.aes = F, aes(y=Density, x = Condition, group=Animal, colour=Genotype),
            linewidth = 1, alpha=0.4)+
  geom_beeswarm(cex = 2.5)+
  scale_colour_manual(values = c("WT" = "#45abff",
                                 "WT + FSK" = "#6AB560",
                                 "TKO" = "#e61c47",
                                 "TKO + FSK" = "#FECC66"))+
  scale_y_continuous(name = expression("Docked SV Density / µm3"),
                     limits = c(0, NA))+
  labs(x="")+
  theme_classic()+
  theme(legend.position = "none")


# QQ-plot to assess residuals ####
## SV density per bouton volume ####
FullModelVolume <- glmer(formula = Density ~ Genotype  * Treatment + (1|Animal),
                         data = DataSetVolume,
                         family = Gamma(link = "log"),
                         control = glmerControl(nAGQ0initStep = 0))

FullModelVolumeLmer <- lmer(formula = Density ~ Genotype  * Treatment + (1|Animal),
                         data = DataSetVolume)
DataSetVolume[,`:=`(`Scaled Gaussian Residuals`=scale(residuals(FullModelVolumeLmer)),`Scaled Gamma Residuals`=scale(residuals(FullModelVolume))),]
ResidualDTVolume <- melt(DataSetVolume, measure.vars = c("Scaled Gaussian Residuals", "Scaled Gamma Residuals"), id.vars = c("Date", "SampleID", "Bouton", "Animal", "Genotype"), variable.name = "Distribution", value.name = "Residuals")

ggplot(data = ResidualDTVolume, aes(sample = Residuals))+
  stat_qq()+
  stat_qq_line()+
  facet_wrap(~Distribution)+
  scale_x_continuous(name = "Theoretical quantiles")+
  scale_y_continuous(name = "Scaled Residuals")+
  theme_classic()

# Build general linear models ####
## SV density per bouton volume ####
FullModelVolume <- glmer(formula = Density ~ Genotype  * Treatment + (1|Animal),
                     data = DataSetVolume,
                     family = Gamma(link = "log"),
                     control = glmerControl(nAGQ0initStep = 0))
summary(FullModelVolume)


NullModelVolume <- glmer(formula = Density ~ 1 + (1|Animal),
                     data = DataSetVolume,
                     family = Gamma(link = "log"),
                     control = glmerControl(nAGQ0initStep = 0))
AllFixedEffectsVolume <- anova(FullModelVolume, NullModelVolume)
print(AllFixedEffectsVolume)


## Treatment, Genotype and Interaction models ####
TreatmentNullModelVolume <- glmer(formula = Density ~ Genotype + (1|Animal),
                     data = DataSetVolume,
                     family = Gamma(link = "log"),
                     control = glmerControl(nAGQ0initStep = 0))
TreatmentFixedEffectsVolume <- anova(FullModelVolume, TreatmentNullModelVolume)
print(TreatmentFixedEffectsVolume)

GenotypeNullModelVolume <- glmer(formula = Density ~ Treatment + (1|Animal),
                         data = DataSetVolume,
                         family = Gamma(link = "log"),
                         control = glmerControl(nAGQ0initStep = 0))
GenotypeFixedEffectsVolume <- anova(FullModelVolume, GenotypeNullModelVolume)
print(GenotypeFixedEffectsVolume)

InteractionNullModelVolume <- glmer(formula = Density ~ Genotype + Treatment + (1|Animal),
                          data = DataSetVolume,
                          family = Gamma(link = "log"),
                          control = glmerControl(nAGQ0initStep = 0))
InteractionFixedEffectsVolume <- anova(FullModelVolume, InteractionNullModelVolume)
print(InteractionFixedEffectsVolume)


# Post-hoc analysis with Estimated Marginal Means ####
## SV density per bouton volume ####
MarginalEffectsVolume <- emmeans(object = FullModelVolume, specs =  ~ Genotype * Treatment, type="response", tran = "log")
MarginalEffectsVolume
ContrastMatrix <- diag(4)

GroupGenotype <- MarginalEffectsVolume@grid$Genotype
GroupTreatment <- MarginalEffectsVolume@grid$Treatment

### Contrast and adjustment for multiple testing (Benjamini-Hochberg) ####
MarginalEffectTestVolume <- contrast(object = MarginalEffectsVolume, type = "response",
                               method = list("WT / (WT + FSK)" = ContrastMatrix[GroupGenotype=="WT"&GroupTreatment=="Control",] - ContrastMatrix[GroupGenotype=="WT"&GroupTreatment=="FSK",],
                                             "WT / TKO" = ContrastMatrix[GroupGenotype=="WT"&GroupTreatment=="Control",] - ContrastMatrix[GroupGenotype=="TKO"&GroupTreatment=="Control",],
                                             "TKO / (TKO + FSK)" = ContrastMatrix[GroupGenotype=="TKO"&GroupTreatment=="Control",] - ContrastMatrix[GroupGenotype=="TKO"&GroupTreatment=="FSK",],
                                             "(WT + FSK) / (TKO + FSK)" = ContrastMatrix[GroupGenotype=="WT"&GroupTreatment=="FSK",] - ContrastMatrix[GroupGenotype=="TKO"&GroupTreatment=="FSK",],
                                             "TKO / (WT + FSK)" = ContrastMatrix[GroupGenotype=="TKO"&GroupTreatment=="Control",] - ContrastMatrix[GroupGenotype=="WT"&GroupTreatment=="FSK",]), adjust="BH")
print(MarginalEffectTestVolume)

### Plot marginal effects ####
(MarginalPlotVolume <- plot(MarginalEffectTestVolume, colors = "grey10") +
   labs(title = "Estimated Docked SV Density Ratio")+
   ylab(label = "")+
   geom_vline(xintercept = 1, linetype="dashed")+
   xlab(label = "Ratio")+
   theme_classic())
