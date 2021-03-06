setwd("C:/ArthropodsArticle/2.0 _Finland_Analysis")

library(car)
library(MASS)
dvlmdvmdlvmdvlmdlmdl

This is a local test to assess stuff



Results <- cbind.data.frame(Variables, Alphas, Betas)
ResultsTable <- Results

#Removing the NA
Results7 <- Results[complete.cases(Results),]

Results4 <- as.data.frame(matrix(Results))
str(Results4)

Results2 <- as.list(cbind.data.frame(Variables, Alphas, Betas))
str(Results2)



#MAKING THE RESULTS EXPORTABLE INTO CSV
Results <- apply(Results, 2 , as.character)

#NAMING THE TRAIL SEGMENTS
rownames(Results) <- rownames(Alphas)

#PASSING RESULTS TO FILE
write.csv(Results, file = "C:/ArthropodsArticle/2.0 _Finland_Analysis/RESULTS2.csv")


Results3 <- read.csv("RESULTS.csv", header = TRUE)
Results8 <- read.csv("RESULTS2.csv", header = TRUE)
Results8
str(Results3)
str(Results8)

head(Results8)
summary(Results8)

hist(Results8$TAlphaAll, col="lightblue")

shapiro.test(Results8$TAlphaAll)

ResultsWithoutControls <- Results3[-c(5,6,11,12,17,18),-c(1:2)]



###########################################################################
#####                           GLMM                                  #####
###########################################################################

# Standard:     glmm <- lme(fruits~bodysize, random = ~1|site)

###
### ALPHAS GLMM###
###

##### All individuals################


#TAlphaAll

glmm_TAlphaAll_AB<-glmer(TAlphaAll ~ Dist_edge + Dist_trail_beginning + (1 | ForestID), data = Results8 , family = poisson)
summary(glmm_TAlphaAll_AB)

glmmT_AlphaAll_A<-glmer(TAlphaAll ~ Dist_edge + (1 | ForestID), data = Results8, family = poisson)
summary(glmmT_AlphaAll_A)

glmm_TAlphaAll_B<-glmer(TAlphaAll ~  Dist_trail_beginning + (1 | ForestID), data = Results8 , family = poisson)

summary(glmm_TAlphaAll_B)

glmm_TAlphaAll_C<-glmer(TAlphaAll ~  treatment + Dist_edge + (1 | ForestID), data = Results8 , family = poisson)
summary(glmm_TAlphaAll_C)

#FAlphaAll

glmm_FAlphaAll_AB <- glmer(FAlphaAll ~ Dist_edge + Dist_trail_beginning + (1 | ForestID), data = Results8 , family = gaussian)
summary(glmm_FAlphaAll_AB)

glmmT_AlphaAll_A<-glmer(FAlphaAll ~ Dist_edge + (1 | ForestID), data = Results8, family = gaussian)
summary(glmmT_AlphaAll_A)

glmm_FAlphaAll_B<-glmer(FAlphaAll ~  Dist_trail_beginning + (1 | ForestID), data = Results8 , family = gaussian)
summary(glmm_FAlphaAll_B)

glmm_FAlphaAll_C<-glmer(FAlphaAll ~  treatment+ (1 | ForestID), data = Results8 , family = gaussian)
# Não está a funcionar
summary(glmm_FAlphaAll_C)


##### Natives ####

#TAlphaAll

glmm_TAlphaNat_AB<-glmer(TAlphaNat ~ Dist_edge + Dist_trail_beginning + (1 | ForestID), data = Results8 , family = poisson)
summary(glmm_TAlphaNat_AB)

glmmT_AlphaNat_A<-glmer(TAlphaNat ~ Dist_edge + (1 | ForestID), data = Results8, family = poisson)
summary(glmmT_AlphaNat_A)

glmm_TAlphaNat_B<-glmer(TAlphaNat ~  Dist_trail_beginning + (1 | ForestID), data = Results8 , family = poisson)
summary(glmm_TAlphaNat_B)

glmm_TAlphaNat_c<-glmer(TAlphaNat ~  Dist_edge + treatment + (1 | ForestID), data = Results8 , family = poisson)
summary(glmm_TAlphaNat_c)

#FAlphaNat

glmm_FAlphaNat_AB <- glmer(FAlphaNat ~ Dist_edge + Dist_trail_beginning + (1 | ForestID), data = Results8 , family = gaussian)
summary(glmm_FAlphaNat_AB)

glmmT_AlphaNat_A<-glmer(FAlphaNat ~ Dist_edge + (1 | ForestID), data = Results8, family = gaussian)
summary(glmmT_AlphaNat_A)

glmm_FAlphaNat_B<-glmer(FAlphaNat ~  Dist_trail_beginning + (1 | ForestID), data = Results8 , family = gaussian)
summary(glmm_FAlphaNat_B)

glmm_FAlphaNat_C<-glmer(FAlphaNat ~  treatment + (1 | ForestID), data = Results8 , family = gaussian)
summary(glmm_FAlphaNat_C)

##### Non-Natives ################

#TAlphaAll

glmm_TAlphaNInd_AB<-glmer(TAlphaNInd ~ Dist_edge + Dist_trail_beginning + (1 | ForestID), data = Results8 , family = poisson)
summary(glmm_TAlphaNInd_AB)

glmmT_AlphaNInd_A<-glmer(TAlphaNInd ~ Dist_edge + (1 | ForestID), data = Results8, family = poisson)
summary(glmmT_AlphaNInd_A)

glmm_TAlphaNInd_B<-glmer(TAlphaNInd ~  Dist_trail_beginning + (1 | ForestID), data = Results8 , family = poisson)
summary(glmm_TAlphaNInd_B)


# O TC250 dá interação significativa
glmm_TAlphaNInd_C<-glmer(TAlphaNInd ~  treatment + (1 | ForestID), data = Results8 , family = poisson)
summary(glmm_TAlphaNInd_C)


#FAlphaNInd

glmm_FAlphaNInd_AB <- glmer(FAlphaNInd ~ Dist_edge + Dist_trail_beginning + (1 | ForestID), data = Results8 , family = gaussian)

summary(glmm_FAlphaNInd_AB)

glmmT_AlphaNInd_A<-glmer(FAlphaNInd ~ Dist_edge + (1 | ForestID), data = Results8, family = gaussian)
summary(glmmT_AlphaNInd_A)

glmm_FAlphaNInd_B<-glmer(FAlphaNInd ~  Dist_trail_beginning + (1 | ForestID), data = Results8 , family = gaussian)
summary(glmm_FAlphaNInd_B)

###
### BETAS GLMM###
###

##### All individuals################


#TBetaAll

glmm_BetaAllTotalVector_AB<-glmer(BetaAllTotalVector ~  treatment + (1 | ForestID), data = Results8 , family = binomial)
summary(glmm_BetaAllTotalVector_AB)

glmmT_BetaAllTotalVector_A<-glmer(BetaAllTotalVector ~ Dist_edge + (1 | ForestID), data = Results8, family = binomial)
summary(glmmT_BetaAllTotalVector_A)

glmm_BetaAllTotalVector_B<-glmer(BetaAllTotalVector ~  Dist_trail_beginning + (1 | ForestID), data = Results8 , family = binomial)
summary(glmm_BetaAllTotalVector_B)

glmm_BetaAllTotalVector_C<-glmer(BetaAllTotalVector ~  treatment + (1 | ForestID), data = Results8 , family = binomial)
summary(glmm_BetaAllTotalVector_C)

#FBetaAll

glmm_FBetaAll_AB <- glmer(BetaFuncAllTotal ~ Dist_edge + treatment  (1 | ForestID), data = Results8 , family = gaussian)
summary(glmm_FAlphaAll_AB)

glmmT_BetaAll_A<-glmer(BetaFuncAllTotal ~ treatment + (1 | ForestID), data = Results8, family = gaussian)
summary(glmmT_AlphaAll_A)

glmm_FBetaAll_B <- glmer(BetaFuncAllTotal ~ Dist_edge + (1 | ForestID), data = Results8 , family = gaussian)
summary(glmm_FAlphaAll_B)

##### natives################


#TBetaAll

glmm_BetaNatTotalVector_AB<-glmer(BetaNatTotalVector ~  treatment + (1 | ForestID), data = Results8 , family = binomial)
summary(glmm_BetaNatTotalVector_AB)

glmmT_BetaNatTotalVector_A<-glmer(BetaNatTotalVector ~ Dist_edge + (1 | ForestID), data = Results8, family = binomial)
summary(glmmT_BetaNatTotalVector_A)

glmm_BetaNatTotalVector_B<-glmer(BetaNatTotalVector ~  Dist_trail_beginning + (1 | ForestID), data = Results8 , family = binomial)
summary(glmm_BetaNatTotalVector_B)

glmm_BetaNatTotalVector_C<-glmer(BetaNatTotalVector ~  treatment + (1 | ForestID), data = Results8 , family = binomial)
summary(glmm_BetaNatTotalVector_C)

#FBetaNat

glmm_FBetaNat_AB <- glmer(BetaFuncNatTotal ~ Dist_edge + treatment  (1 | ForestID), data = Results8 , family = gaussian)
summary(glmm_FAlphaNat_AB)

glmmT_BetaNat_A<-glmer(BetaFuncNatTotal ~ treatment + (1 | ForestID), data = Results8, family = gaussian)
summary(glmmT_AlphaNat_A)

glmm_FBetaNat_B <- glmer(BetaFuncNatTotal ~ Dist_edge + (1 | ForestID), data = Results8 , family = gaussian)
summary(glmm_FAlphaNat_B)

