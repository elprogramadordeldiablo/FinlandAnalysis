###
### ALPHAS GLMM###
###

##### All individuals################


#TAlphaAll

glmm_TAlphaAll_AB<-glmer(TAlphaAll ~ Dist_edge + Dist_trail_beginning + (1 | ForestID), data = Results2 , family = poisson)
summary(glmm_TAlphaAll_AB)

glmmT_AlphaAll_A<-glmer(TAlphaAll ~ Dist_edge + (1 | ForestID), data = Results2, family = poisson)
summary(glmmT_AlphaAll_A)

glmm_TAlphaAll_B<-glmer(TAlphaAll ~  Dist_trail_beginning + (1 | ForestID), data = Results2 , family = poisson)
summary(glmm_TAlphaAll_B)

glmm_TAlphaAll_C<-glmer(TAlphaAll ~  treatment + Dist_edge + (1 | ForestID), data = Results2 , family = poisson)
summary(glmm_TAlphaAll_C)

#FAlphaAll

glmm_FAlphaAll_AB <- glmer(FAlphaAll ~ Dist_edge + Dist_trail_beginning + (1 | ForestID), data = Results2 , family = gaussian)
summary(glmm_FAlphaAll_AB)

glmmT_AlphaAll_A<-glmer(FAlphaAll ~ Dist_edge + (1 | ForestID), data = Results2, family = gaussian)
summary(glmmT_AlphaAll_A)

glmm_FAlphaAll_B<-glmer(FAlphaAll ~  Dist_trail_beginning + (1 | ForestID), data = Results2 , family = gaussian)
summary(glmm_FAlphaAll_B)

glmm_FAlphaAll_C<-glmer(FAlphaAll ~  treatment+ (1 | ForestID), data = Results2 , family = gaussian)
summary(glmm_FAlphaAll_C)


##### Natives ####

#TAlphaAll

glmm_TAlphaNat_AB<-glmer(TAlphaNat ~ Dist_edge + Dist_trail_beginning + (1 | ForestID), data = Results2 , family = poisson)
summary(glmm_TAlphaNat_AB)

glmmT_AlphaNat_A<-glmer(TAlphaNat ~ Dist_edge + (1 | ForestID), data = Results2, family = poisson)
summary(glmmT_AlphaNat_A)

glmm_TAlphaNat_B<-glmer(TAlphaNat ~  Dist_trail_beginning + (1 | ForestID), data = Results2 , family = poisson)
summary(glmm_TAlphaNat_B)

glmm_TAlphaNat_c<-glmer(TAlphaNat ~  Dist_edge + treatment + (1 | ForestID), data = Results2 , family = poisson)
summary(glmm_TAlphaNat_c)

#FAlphaNat

glmm_FAlphaNat_AB <- glmer(FAlphaNat ~ Dist_edge + Dist_trail_beginning + (1 | ForestID), data = Results2 , family = gaussian)
summary(glmm_FAlphaNat_AB)

glmmT_AlphaNat_A<-glmer(FAlphaNat ~ Dist_edge + (1 | ForestID), data = Results2, family = gaussian)
summary(glmmT_AlphaNat_A)

glmm_FAlphaNat_B<-glmer(FAlphaNat ~  Dist_trail_beginning + (1 | ForestID), data = Results2 , family = gaussian)
summary(glmm_FAlphaNat_B)

glmm_FAlphaNat_C<-glmer(FAlphaNat ~  treatment + (1 | ForestID), data = Results2 , family = gaussian)
summary(glmm_FAlphaNat_C)

##### Non-Natives ################

#TAlphaAll

glmm_TAlphaNInd_AB<-glmer(TAlphaNInd ~ Dist_edge + Dist_trail_beginning + (1 | ForestID), data = Results2 , family = poisson)
summary(glmm_TAlphaNInd_AB)

glmmT_AlphaNInd_A<-glmer(TAlphaNInd ~ Dist_edge + (1 | ForestID), data = Results2, family = poisson)
summary(glmmT_AlphaNInd_A)



glmm_TAlphaNInd_B<-glmer(TAlphaNInd ~  Dist_trail_beginning + (1 | ForestID), data = Results2 , family = poisson)
summary(glmm_TAlphaNInd_B)


# O TC250 dá interação significativa
glmm_TAlphaNInd_C<-glmer(TAlphaNInd ~  treatment + (1 | ForestID), data = Results2 , family = poisson)
summary(glmm_TAlphaNInd_C)


#FAlphaNInd

glmm_FAlphaNInd_AB <- glmer(FAlphaNInd ~ Dist_edge + Dist_trail_beginning + (1 | ForestID), data = Results2 , family = gaussian)
summary(glmm_FAlphaNInd_AB)

glmmT_AlphaNInd_A<-glmer(FAlphaNInd ~ Dist_edge + (1 | ForestID), data = Results2, family = gaussian)
summary(glmmT_AlphaNInd_)A

glmm_FAlphaNInd_B<-glmer(FAlphaNInd ~  Dist_trail_beginning + (1 | ForestID), data = Results2 , family = gaussian)
summary(glmm_FAlphaNInd_B)