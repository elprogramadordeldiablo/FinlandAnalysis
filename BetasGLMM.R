###
### BETAS GLMM###
###

##### All individuals################


#TBetaAll

glmm_BetaAllTotalVector_AB<-glmer(BetaAllTotalVector ~  treatment + (1 | ForestID), data = Results2 , family = binomial)
summary(glmm_BetaAllTotalVector_AB)

glmmT_BetaAllTotalVector_A<-glmer(BetaAllTotalVector ~ Dist_edge + (1 | ForestID), data = Results2, family = binomial)
summary(glmmT_BetaAllTotalVector_A)

glmm_BetaAllTotalVector_B<-glmer(BetaAllTotalVector ~  Dist_trail_beginning + (1 | ForestID), data = Results2 , family = binomial)
summary(glmm_BetaAllTotalVector_B)

glmm_BetaAllTotalVector_C<-glmer(BetaAllTotalVector ~  treatment + (1 | ForestID), data = Results2 , family = binomial)
summary(glmm_BetaAllTotalVector_C)

#FBetaAll

glmm_FBetaAll_AB <- glmer(BetaFuncAllTotalVector~ Dist_edge + treatment + (1 | ForestID), data = Results2 , family = binomial)
summary(glmm_FBetaAll_AB)

glmm_FBetaAll_A<-glmer(BetaFuncAllTotalVector ~ treatment + (1 | ForestID), data = Results2, family = binomial)
summary(glmm_FBetaAll_A)

glmm_FBetaAll_B <- glmer(BetaFuncAllTotalVector ~ Dist_edge + (1 | ForestID), data = Results2 , family = binomial)
summary(glmm_FBetaAll_B)

##### natives################


#TBetaAll

glmm_BetaNatTotalVector_AB<-glmer(BetaNatTotalVector ~  treatment + (1 | ForestID), data = Results2 , family = binomial)
summary(glmm_BetaNatTotalVector_AB)

glmmT_BetaNatTotalVector_A<-glmer(BetaNatTotalVector ~ Dist_edge + (1 | ForestID), data = Results2, family = binomial)
summary(glmmT_BetaNatTotalVector_A)

glmm_BetaNatTotalVector_B<-glmer(BetaNatTotalVector ~  Dist_trail_beginning + (1 | ForestID), data = Results2 , family = binomial)
summary(glmm_BetaNatTotalVector_B)

glmm_BetaNatTotalVector_C<-glmer(BetaNatTotalVector ~  treatment + (1 | ForestID), data = Results2 , family = binomial)
summary(glmm_BetaNatTotalVector_C)

#FAlphaNat

glmm_FBetaNat_AB <- glmer(BetaFuncNatTotalVector ~ Dist_edge + treatment + (1 | ForestID), data = Results2 , family = binomial)
summary(glmm_FBetaNat_AB)

glmmT_BetaNat_A<-glmer(BetaFuncNatTotalVector ~ treatment + (1 | ForestID), data = Results2, family = binomial)
summary(glmmT_AlphaNat_A)

glmm_FBetaNat_B <- glmer(BetaFuncNatTotalVector ~ Dist_edge + (1 | ForestID), data = Results2 , family = binomial)
summary(glmm_FBetaNat_B)
