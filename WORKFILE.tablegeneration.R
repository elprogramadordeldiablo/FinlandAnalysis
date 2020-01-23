options(na.action = "na.fail")        
library(BAT)
library(readr)
library(FD) 
library(car)
library(MASS)
library(lme4)
library(here)       # to locate files
library(data.table) # to work with data
library(dplyr)      # to manage data
library(magrittr)   # to use the pipe operator %>% 
library(MuMIn)
library(glmmTMB)
library(bbmle)
library(performance)

Models.test = dredge(glmmTMB(abund.all ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data= Results2, family = poisson), extra =list(R2 = function(x)     {
  r.squaredGLMM(x, null = nullmodel)["delta ", ]})) 
summary(Models.test)
str(Models.test)
performance::r2(Models.test)
check_overdispersion(Models.test)
Models.test$R2.R2m
model.extract(Models.test)

Models1 <- dredge(glmmTMB(TAlphaAll ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data= Results2, family = "poisson"))
summary(Models1)

TalphaAll.test <- dredge(glmmTMB(TAlphaAll ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std +    (1 | ForestID), data= Results2,family = "poisson") )
str(TalphaAll.test)
model.output(TalphaAll.test)


model1.test <- glmmTMB(abund.all ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data= Results2, family = poisson)
summary(model1.test)
performance::r2(model1.test)
r.squaredGLMM(model1.test)
check_overdispersion(Model1.test)
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
overdisp_fun(model1.test)
p1 <- plot(model1.test,id=0.05,idLabels=~.obs)

str(model1.test)
model1.test$
model.extract(model1.test,"weights")
# (all.tax.btotal, all.tax.brich, all.tax.brepl, nat.tax.btotal, nat.tax.brich, nat.tax.brepl, nind.tax.btotal, nind.tax.brich, nind.tax.brepl, all.func.btotal, all.func.brich, 
# all.func.brepl, nat.func.btotal, nat.func.brich, nat.func.brepl, nind.func.btotal, nind.func.brich, nind.func.brepl)


fam3 <- c(NA, NA, NA, NA,NA, NA,NA, "poisson", "poisson", "poisson", "Gamma", "Gamma", "Gamma","nbinom1", "nbinom1", "nbinom1","beta_family", "beta_family","beta_family","beta_family", "beta_family", 
          "beta_family", "beta_family", "beta_family", "beta_family", "beta_family", "beta_family", "beta_family", "beta_family", "beta_family", "beta_family", "beta_family", "beta_family", 
          "beta_family", "beta_family", "beta_family", "beta_family")


str(fam3)
fam3

fam3backup <- c(NA, NA, NA, NA,NA, NA,NA, "poisson", "poisson", "poisson", "Gamma", "Gamma", "Gamma", "beta_family", "beta_family", 
                "beta_family", "beta_family", "beta_family", "beta_family", "beta_family", "beta_family", "beta_family", "beta_family", "beta_family", "beta_family", "beta_family", "beta_family", 
                "beta_family", "beta_family", "beta_family", "beta_family")
write.csv(fam3backup, file = here("results","ddistr.families.csv"), row.names = TRUE)

#### Teste de código para obter weights

aic.weights = c()
models.df2 = data.frame()
Models = model.extract() 
for(i in 8:37){
  newTable = Results2[,c(i,1,5,6,7)]
  colnames(newTable)[1] = "y"
  Models[[i]] = dredge(glmmTMB(y ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data= newTable, family = fam3[i]))
  names(Models)[[i]] = colnames(Results2[i])
  Models[[i]] = apply(Models[[i]], 2 , as.numeric)
  #############################

  w = c()
  for(k in 1:3){
    w[k] = sum(Models[[i]][!is.na(Models[[i]][,(k+2)]),10])
  }
  aic.weights = rbind(aic.weights, w)

    
  #############################
  Models[[i]] = cbind(rep(colnames(Results2)[i], nrow(Models[[i]])), Models[[i]])
  models.df2 = rbind(models.df, as.data.frame(Models[[i]]))

  #for(j in 2:ncol(models.df)){
  #  models.df[,j] <- as.numeric(models.df[,j])
  #}
}
colnames(aic.weights) = colnames(models.df)[3:5]
rownames(aic.weights) = unique(models.df[,1])
models.df2
aic.weights

#selecting delta <2

AICmodels <- filter(models.df, delta < 2)
AICmodels


str(models.df)
aic.weights

write.csv(resultspart1, file = here("results","dredge.models2.csv"), row.names = TRUE)
write.csv(AICmodels, file = here("results","aic.models2.csv"), row.names = TRUE)
write.csv(aic.weights, file = here("results","aic.weights.csv"), row.names = TRUE)

###### TESTE PARA OBTER R^2

models.df3 = data.frame()
Models3 = list() 
for(i in 8:37){
  newTable = Results2[,c(i,1,5,6,7)]
  colnames(newTable)[1] = "y"
  Models3[[i]] = dredge(glmer(y ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data= newTable, family = fam3[i]))
  names(Models3)[[i]] = colnames(Results2[i])
  Models3[[i]] = apply(Models3[[i]], 2 , as.numeric)
  Models3[[i]] = cbind(rep(colnames(Results2)[i], nrow(Models[[i]])), Models3[[i]])
  models.df3 = rbind(models.df, as.data.frame(Models3[[i]]))
  # for(j in 2:ncol(models.df)){
  #   models.df3[,j] <- as.numeric(models.df3[,j])
  # }
}
models.df3

Models

###### Backup di codigo original para criar a tabela
models.df = data.frame()
Models = list() 
for(i in 8:37){
  newTable = Results2[,c(i,1,5,6,7)]
  colnames(newTable)[1] = "y"
  Models[[i]] = dredge(glmmTMB(y ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data= newTable, family = fam3[i]))
  names(Models)[[i]] = colnames(Results2[i])
  Models[[i]] = apply(Models[[i]], 2 , as.numeric)
  Models[[i]] = cbind(rep(colnames(Results2)[i], nrow(Models[[i]])), Models[[i]])
  models.df = rbind(models.df, as.data.frame(Models[[i]]))
  for(j in 2:ncol(models.df)){
    models.df[,j] <- as.numeric(models.df[,j])
  }
}

Models
models.df

resultspart1 <- models.df
unlisted.models <- unlist(Models)
str(unlisted.models)
unlisted.models <- data.frame(unlisted.models)
