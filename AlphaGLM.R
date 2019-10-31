###
### ALPHAS GLMM###
###

##### All individuals################


# Testing several forms to evaluate the model ----

#TAlphaAll

glmm_TAlphaAll_AB <-glmer(TAlphaAll ~ Dist_edge + Dist_trail_beginning + (1 | ForestID), data = Results2 , 
                          family = poisson, control = glmerControl(optimizer="bobyqa", 
                                                                   optCtrl = list(maxfun = 1000000)))
summary(glmm_TAlphaAll_AB)

#plotting the residuals for an exploratory evaluation of the model. Will checkif they are a approximate constant,
# or have outliers. Will check for constant variation across the fitted range. 
plot(glmm_TAlphaAll_AB) 

# For generalized models it is often more useful to examine the residuals plotted on the link scale, ??,
# instead of the response scale

ggplot(data.frame(eta=predict(glmm_TAlphaAll_AB,type="link"),pearson=residuals(glmm_TAlphaAll_AB,type="pearson")),
       aes(x=eta,y=pearson)) +
  geom_point() +
  theme_bw()

#Checking linearity in each variable

ggplot(data.frame(x1="Results2$Dist_trail",pearson=residuals(glmm_TAlphaAll_AB,type="pearson")),
       aes(x=x1,y=pearson)) +
  geom_point() +
  theme_bw()

ggplot(data.frame(x2="Results2$Dist_edge",pearson=residuals(glmm_TAlphaAll_AB,type="pearson")),
       aes(x=x2,y=pearson)) +
  geom_point() +
  theme_bw()

ggplot(data.frame(x3="Results2$Dist_trail_beginning",pearson=residuals(glmm_TAlphaAll_AB,type="pearson")),
       aes(x=x3,y=pearson)) +
  geom_point() +
  theme_bw()

#Assessing for Independent variable independence

means <- aggregate(Results2[,c("Dist_trail","Dist_edge")],by=list(Results2$ForestID),FUN=mean)
lmcoefs <- summary(lm(TAlphaAll ~ Dist_trail + Dist_edge + ForestID, data=Results2))$coefficients[,"Estimate"]
means$effects <- c(0,lmcoefs[substr(names(lmcoefs),1,2) == "ForestID"])
means$effects <- means$effects - mean(means$effects)

cor(means[,c("Dist_trail","Dist_edge","effects")])


## Assessing the nomality of the residuals

qqnorm(residuals(glmm_TAlphaAll_AB)) # it seems fine

### analysing sensitivity to data 

ggplot(data.frame(lev=hatvalues(glmm_TAlphaAll_AB),pearson=residuals(glmm_TAlphaAll_AB,type="pearson")),
       aes(x=lev,y=pearson)) +
  geom_point() +
  theme_bw()
# ir responds with Warning message:
#   In hatvalues.merMod(glmm_TAlphaAll_AB) :
#   the hat matrix may not make sense for GLMMs

# Determining which of the observations have the highest leverage and displaying these observations. 
# Generation a new model without these observations and then comparing the coefficients 
# for the will all observations to this new model with some observations removed.

levId <- which(hatvalues(glmm_TAlphaAll_AB) >= .172)
# ##Result: Warning message:
# In hatvalues.merMod(glmm_TAlphaAll_AB) :
#   the hat matrix may not make sense for GLMMs

Results2[levId,c("TAlphaAll","Dist_trail","Dist_edge","ForestID")]

summary(Results2[,c("TAlphaAll","Dist_trail","Dist_edge")])

mmLev <- lmer(y ~ x1 + x2 + (1|g1), data=pbDat[-c(levId),])
mmLevCD <- data.frame(effect=fixef(mm),
                      change=(fixef(mmLev) - fixef(mm)),
                      se=sqrt(diag(vcov(mm)))
)
rownames(mmLevCD) <- names(fixef(mmLev))
mmLevCD$multiples <- abs(mmLevCD$change / mmLevCD$se)
mmLevCD



# The following tests ----
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 


glmmT_AlphaAll_A<-glmer(TAlphaAll ~ Dist_edge + (1 | ForestID), data = Results2, 
                        family = poisson, control = glmerControl(optimizer="bobyqa", 
                                                                 optCtrl = list(maxfun = 1000000)) )
summary(glmmT_AlphaAll_A)
plot(glmmT_AlphaAll_A)


glmm_TAlphaAll_B<-glmer(TAlphaAll ~  Dist_trail_beginning + (1 | ForestID), data = Results2 ,
                        family = poisson, control = glmerControl(optimizer="bobyqa", 
                                                                 optCtrl = list(maxfun = 1000000)))
summary(glmm_TAlphaAll_B)

glmm_TAlphaAll_C<-glmer(TAlphaAll ~  treatment + Dist_edge + (1 | ForestID), data = Results2 , 
                        family = poisson, control = glmerControl(optimizer="bobyqa", 
                                                                 optCtrl = list(maxfun = 1000000)))
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

glmm_TAlphaNInd_AB<-glmer(TAlphaNInd ~ Dist_edge + Dist_trail_beginning + (1 | ForestID), data = Results2 , 
                          family = poisson, control = glmerControl(optimizer="bobyqa", 
                                                                   optCtrl = list(maxfun = 1000000)))
summary(glmm_TAlphaNInd_AB)

glmmT_AlphaNInd_A<-glmer(TAlphaNInd ~ Dist_edge + (1 | ForestID), data = Results2, 
                         family = poisson, control = glmerControl(optimizer="bobyqa", 
                                                                  optCtrl = list(maxfun = 1000000)))
summary(glmmT_AlphaNInd_A)



glmm_TAlphaNInd_B<-glmer(TAlphaNInd ~  Dist_trail_beginning + (1 | ForestID), data = Results2 , 
                         family = poisson, control = glmerControl(optimizer="bobyqa", 
                                                                  optCtrl = list(maxfun = 1000000)))
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

