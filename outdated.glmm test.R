setwd("C:/ArthropodsArticle/2.0 _Finland_Analysis")

library(car)
library(MASS)
library(lme4)


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
Results8 <- read.csv("RESULTS2.csv", header = TRUE, sep=";")
Results8
str(Results3)
str(Results8)

ResultsWithoutControls <- Results3[-c(5,6,11,12,17,18),-c(1:2)]


###########################################################################
#####                           GLMM                                  #####
###########################################################################

glmm <- lme(fruits~bodysize, random = ~1|site)

###
### ALPHAS GLMM###
###

##### All individuals################


#TAlphaAll

glmm_TAlphaAll_AB<-glmer(TAlphaAll ~ Dist_edge + Dist_trail_beginning + (1 | ForestID), data = Results8 , family = poisson)
summary(glmm_TAlphaAll_AB)

glmmT_AlphaAll_A<-glmer(TAlphaAll ~ Dist_edge + (1 | ForestID), data = Results2, family = poisson)
summary(glmmT_AlphaAll_A)

glmm_TAlphaAll_B<-glmer(TAlphaAll ~  Dist_trail_beginning + (1 | ForestID), data = Results2 , family = poisson)

summary(glmm_TAlphaAll_B)

glmm_TAlphaAll_C<-glmer(TAlphaAll ~  treatment + Dist_edge + (1 | ForestID), data = Results2 , family = poisson)

summary(glmm_TAlphaAll_C)

#FAlphaAll

glmm_FAlphaAll_AB <- glmer(FAlphaAll ~ Dist_edge + Dist_trail_beginning + (1 | ForestID), data = Results2 , family = poisson)
summary(glmm_FAlphaAll_AB)

glmmT_AlphaAll_A<-glmer(FAlphaAll ~ Dist_edge + (1 | ForestID), data = Results2, family = poisson)
summary(glmmT_AlphaAll_A)

glmm_FAlphaAll_B<-glmer(FAlphaAll ~  Dist_trail_beginning + (1 | ForestID), data = Results2 , family = poisson)
summary(glmm_FAlphaAll_B)

glmm_FAlphaAll_C<-glmer(FAlphaAll ~  treatment+ (1 | ForestID), data = Results2 , family = poisson)
# Não está a funcionar
summary(glmm_FAlphaAll_C)


##### Natives ####

#TAlphaAll

glmm_TAlphaNat_AB<-glmer(TAlphaNat ~ Dist_edge + (1 | ForestID), data = Results2 , family = poisson)
summary(glmm_TAlphaNat_AB)

glmmT_AlphaNat_A<-glmer(TAlphaNat ~ Dist_edge + (1 | ForestID), data = Results2, family = poisson)
summary(glmmT_AlphaNat_A)



glmm_TAlphaNat_c<-glmer(TAlphaNat ~  Dist_edge + treatment + (1 | ForestID), data = Results2 , family = poisson)
summary(glmm_TAlphaNat_c)

#FAlphaNat

glmm_FAlphaNat_AB <- glmer(FAlphaNat ~ Dist_edge +  (1 | ForestID), data = Results8 , family = poisson)
summary(glmm_FAlphaNat_AB)

glmmT_AlphaNat_A<-glmer(FAlphaNat ~ Dist_edge + (1 | ForestID), data = Results2, family = poisson)
summary(glmmT_AlphaNat_A)

glmm_FAlphaNat_B<-glmer(FAlphaNat ~  Dist_trail_beginning + (1 | ForestID), data = Results2 , family = poisson)
summary(glmm_FAlphaNat_B)

glmm_FAlphaNat_C<-glmer(FAlphaNat ~  treatment + (1 | ForestID), data = Results2 , family = poisson)
summary(glmm_FAlphaNat_C)

##### Non-Natives ################

#TAlphaAll

#Quase significativo - 0,06
glmm_TAlphaNInd_AB<-glmer(TAlphaNInd ~ Dist_edge + Dist_trail_beginning + treatment + (1 | ForestID), data = Results2 , family = poisson)
summary(glmm_TAlphaNInd_AB)

glmmT_AlphaNInd_A<-glmer(TAlphaNInd ~ Dist_edge + (1 | ForestID), data = Results2, family = poisson)
summary(glmmT_AlphaNInd_A)

glmm_TAlphaNInd_B<-glmer(TAlphaNInd ~  Dist_trail_beginning + (1 | ForestID), data = Results2 , family = poisson)    
summary(glmm_TAlphaNInd_B)


# O TC250 dá interação significativa
glmm_TAlphaNInd_C<-glmer(TAlphaNInd ~  treatment + (1 | ForestID), data = Results2 , family = poisson)
summary(glmm_TAlphaNInd_C)


#FAlphaNInd

glmm_FAlphaNInd_AB <- glmer(FAlphaNInd ~ Dist_edge + Dist_trail_beginning + (1 | ForestID), data = Results2 , family = poisson)

summary(glmm_FAlphaNInd_AB)

glmmT_AlphaNInd_A<-glmer(FAlphaNInd ~ Dist_edge + (1 | ForestID), data = Results2, family = poisson)
summary(glmmT_AlphaNInd_A)

glmm_FAlphaNInd_B<-glmer(FAlphaNInd ~  treatment+ (1 | ForestID), data = Results2 , family = gaussian)


summary(glmm_FAlphaNInd_B)

###
### BETAS GLMM###
###

##### All individuals################


#TBetaAll

glmm_BetaAllTotalVector_AB<-glmer(BetaAllTotalVector ~  treatment + (1 | ForestID), data = Results8, family = binomial)
summary(glmm_BetaAllTotalVector_AB)

glmm_BetaAllTotalVector_A<-glmer(BetaAllTotalVector ~ Dist_edge + (1 | ForestID), data = Results8, family = binomial)
summary(glmm_BetaAllTotalVector_A)

glmm_BetaAllTotalVector_B<-glmer(BetaAllTotalVector ~  Dist_trail_beginning + (1 | ForestID), data = Results8 , family = binomial)
summary(glmm_BetaAllTotalVector_B)

glmm_BetaAllTotalVector_C<-glmer(BetaAllTotalVector ~  treatment + Dist_trail_beginning + (1 | ForestID), data = Results8 , family = binomial)
summary(glmm_BetaAllTotalVector_C)

#FAlphaAll

glmm_FAlphaAll_AB <- glmer(FAlphaAll ~ Dist_edge + Dist_trail_beginning + (1 | ForestID), data = Results8 , family = gaussian)
summary(glmm_FAlphaAll_AB)

glmmT_AlphaAll_A<-glmer(FAlphaAll ~ Dist_edge + (1 | ForestID), data = Results8, family = gaussian)
summary(glmmT_AlphaAll_A)

glmm_FAlphaAll_B<-glmer(FAlphaAll ~  Dist_trail_beginning + (1 | ForestID), data = Results8 , family = gaussian)
summary(glmm_FAlphaAll_B)

#
### BETAS GLMM###
###

##### All individuals################


#TBetaAll

glmm_BetaAllTotalVector_AB<-glmer(BetaAllTotalVector ~  treatment + (1 | ForestID), data = Results4 , family = poisson)
summary(glmm_BetaAllTotalVector_AB)

glmmT_BetaAllTotalVector_A<-glmer(BetaAllTotalVector ~ Dist_edge + (1 | ForestID), data = Results2, family = poisson)
summary(glmmT_BetaAllTotalVector_A)

glmm_BetaAllTotalVector_B<-glmer(BetaAllTotalVector ~  Dist_trail_beginning + (1 | ForestID), data = Results2 , family = poisson)
summary(glmm_BetaAllTotalVector_B)

glmm_BetaAllTotalVector_C<-glmer(BetaAllTotalVector ~  treatment + Dist_trail_beginning + (1 | ForestID), data = Results4 , family = poisson)

summary(glmm_BetaAllTotalVector_C)

#FAlphaAll

glmm_FAlphaAll_AB <- glmer(FAlphaAll ~ Dist_edge + Dist_trail_beginning + (1 | ForestID), data = Results8 , family = gaussian)
summary(glmm_FAlphaAll_AB)

glmmT_AlphaAll_A<-glmer(FAlphaAll ~ Dist_edge + (1 | ForestID), Results8 , family = gaussian)
summary(glmmT_AlphaAll_A)

glmm_FAlphaAll_B<-glmer(FAlphaAll ~  Dist_trail_beginning + (1 | ForestID), Results8 , family = gaussian)
summary(glmm_FAlphaAll_B)









########ARTHROPODS ARTICLE########
#I did not log transform abundance tables. Prepare table.
#Use PCA to pick apart correlated traits? 
#Use kdist procedure to check for correlation between traits
#Check if the trait values are scaled
#show suitability tests to make the glmm. Still, not sure if I did the necessary for glmm testing
# Do a null model for alphas and Betas?
#Reapply line 138  (function cor)to all dendrogram applications (see if the correlation value is high enough (>80 aprox))
#Square root or log transform the abundance matrix
#Use MDS to understand how functional and taxonomical traits vary along treatments (need more explaining)
#Testing without 50m sampling area
#




