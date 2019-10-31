

glmA1<-glm(BetaAllTotalVector ~  treatment + Dist_edge, data = Results2, family = binomial)
summary(glmA1)

glmA2<-glm(BetaAllTotalVector ~  treatment, data = Results2, family = binomial)
summary(glmA2)

glmA3<-glm(BetaAllTotalVector ~ ForestID, data = Results2, family = binomial )
summary(glmA3)

glmA4<-glm(BetaAllTotalVector ~ Dist_edge, data = Results2, family = binomial)
summary(glmA4)

glmA5<-glm(BetaAllTotalVector ~ Dist_trail_beginning, data = Results2, family = binomial() )
summary(glmA5)

glmA6<-glm(BetaAllTotalVector ~ Dist_trail, data = Results2, family = binomial() )
summary(glmA6)

glmA7<-glm(BetaAllTotalVector ~ Dist_trail + Dist_trail_beginning + Dist_edge + ForestID + treatment, data = Results2, family = binomial() )
summary(glmA7)

glmA8<-glm(BetaAllTotalVector ~ Dist_trail + Dist_trail_beginning + ForestID + treatment, data = Results2, family = binomial() )
summary(glmA8)

glmA9<-glm(BetaAllTotalVector ~ Dist_trail + Dist_trail_beginning + Dist_edge + ForestID + treatment, data = Results2, family = binomial() )
summary(glmA9)

glmA10<-glm(BetaAllTotalVector ~ Dist_trail + Dist_trail_beginning + Dist_edge + ForestID + treatment, data = Results2, family = binomial)
summary(glmA10)




#Beta Functional

glmB1<-glm(BetaAllTotalVector ~  treatment + ForestID, data = Results2 )
summary(glmB1)

glmB2<-glm(BetaAllTotalVector ~  treatment, data = Results2 )
summary(glmB2)

glmB3<-glm(BetaAllTotalVector ~ ForestID, data = Results2 )
summary(glmB3)

glmB4<-glm(BetaAllTotalVector ~ Dist_edge, data = Results2 )
summary(glmB4)

glmB5<-glm(BetaAllTotalVector ~ Dist_trail_beginning, data = Results2 )
summary(glmB5)

glmB6<-glm(BetaAllTotalVector ~ Dist_trail, data = Results2 )
summary(glmB6)

