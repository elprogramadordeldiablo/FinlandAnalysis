
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

# (all.tax.btotal, all.tax.brich, all.tax.brepl, nat.tax.btotal, nat.tax.brich, nat.tax.brepl, nind.tax.btotal, nind.tax.brich, nind.tax.brepl, all.func.btotal, all.func.brich, 
# all.func.brepl, nat.func.btotal, nat.func.brich, nat.func.brepl, nind.func.btotal, nind.func.brich, nind.func.brepl)


fam3 <- c(NA, NA, NA, NA,NA, NA,NA, "poisson", "poisson", "poisson", "Gamma", "Gamma", "poisson", "beta_family", "beta_family", 
          "beta_family", "beta_family", "beta_family", "beta_family", "beta_family", "beta_family", "beta_family", "beta_family", "beta_family", "beta_family", "beta_family", "beta_family", 
          "beta_family", "beta_family", "beta_family", "beta_family")
fam3





#selecting delta <2

AICmodels <- filter(models.df, delta < 2)
AICmodels

Count <- colSums(!is.na(models.df[,4:6]))
Sum <- apply(!is.na(models.df[,4:6]), 2, function(x) sum(models.df$Values[x]))
data.frame(Count, Sum)

#### Teste de código para obter weights
ind.vars = data.frame()
aic.weights = data.frame()
myfun = function(x){sum(Models$weight[x])}
models.df = data.frame()
Models = list() 
for(i in 8:31){
  newTable = Results2[,c(i,1,5,6,7)]
  colnames(newTable)[1] = "y"
  Models[[i]] = dredge(glmmTMB(y ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data= newTable, family = fam3[i]))
  names(Models)[[i]] = colnames(Results2[i])
  Models[[i]] = apply(Models[[i]], 2 , as.numeric)
  #############################
  for(k in 3:5){
    ind.vars[[k]]= cbind(ind.vars, apply(!is.na(Models[,k]), 2, myfun))    ### Error in Models[, k] : incorrect number of dimensions
    aic.weights[[i]]= rbind(aic.weights, ind.vars[[k]])
    aic.weights[[i]]= cbind(colnames(Results2[i], aic.weights[[i]]))
    }
  #############################
  Models[[i]] = cbind(rep(colnames(Results2)[i], nrow(Models[[i]])), Models[[i]])
  models.df = rbind(models.df, as.data.frame(Models[[i]]))

  for(j in 2:ncol(models.df)){
    models.df[,j] <- as.numeric(models.df[,j])
  }
}
models.df
aic.weights

ind.vars




###### Backup di codigo original para criar a tabela
models.df = data.frame()
Models = list() 
for(i in 8:31){
  newTable = Results2[,c(i,1,5,6,7)]
  colnames(newTable)[1] = "y"
  Models[[i]] = dredge(glmmTMB(y ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data= newTable, family = fam3[i]))
  names(Models)[[i]] = colnames(Results2[i])
  Models[[i]] = apply(Models[[i]], 2 , as.numeric)
  Models[[i]] = cbind(rep(colnames(Results2)[i], nrow(Models[[i]])), Models[[i]])
  #test <- rbind(models.df, cbind(names(Models)[[i]],Models[[i]]))
  models.df = rbind(models.df, as.data.frame(Models[[i]]))
        #models.df = cbind(names(Models)[[i]],rbind(models.df, apply(Models[[i]], 2 , as.numeric)))
        #models.df = cbind(names(Models)[[i]],rbind(models.df, apply(Models[[i]], 2 , as.numeric)))
         #models.df2 =cbind(names(Models)[[i]], models.df)
  for(j in 2:ncol(models.df)){
    models.df[,j] <- as.numeric(models.df[,j])
              }
      }
models.df

Models
