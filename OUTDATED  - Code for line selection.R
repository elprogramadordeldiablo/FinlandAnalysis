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

sample1 = c(23, 45, 2, 0, NA, NA)
sample2 = c(45, 12, 3, 16, 45, 34)

data = rbind(sample1, sample2)
colnames(data) = c("distance", "plot", "l4m0", "l5m4", "l3m5", "l2m4")

for(i in 1:nrow(data)){
  basicInfo = data[i, 1:2]
  for(j in seq(3, 6, 2)){
    newSp = data[i, c(j, j+1)]
    if(!any(is.na(newSp)))
      newData = rbind(newData,c(basicInfo, colnames(data)[j], newSp))
  }
}

is.numeric(as.data.frame(newData)$V1)



#Working code


sample1 = c(23, 45, 2, 0, NA, NA)
sample2 = c(45, 12, 3, 16, 45, 34)

newData <- list() # este fui eu que acrescentei
data = rbind(sample1, sample2)
colnames(data) = c("distance", "plot", "l4m0", "l5m4", "l3m5", "l2m4")

 for(i in 1:nrow(data)){            # 
    basicInfo = data[i, 1:2]
      for(j in seq(3, 6, 2)){
           newSp = data[i, c(j, j+1)]
           if(!any(is.na(newSp)))
               newData = rbind(newData,c(basicInfo, colnames(data)[j], newSp))
           }
     }
 newData
 
 glmm.families <- glmm.families[-1,-2]
 glmm.families
fam <- read.csv2(here("data", "glmm.families.csv"), row.names= 1, header = TRUE, sep = ",", dec = ".")


ww <- (Results2[,10])
ww2 <- fam3[10]
ww2




fam3 <- c(NA, NA, NA, NA,NA, NA,NA, "poisson", "poisson", "poisson", "Gamma", "Gamma", "poisson", "beta_family", "beta_family", 
                "beta_family", "beta_family", "beta_family", "beta_family", "beta_family", "beta_family", "beta_family", "beta_family", "beta_family", "beta_family", "beta_family", "beta_family", 
                "beta_family", "beta_family", "beta_family", "beta_family")
fam3

#Working version!!!

Models = list() 
for(i in 8:31){
  cat(i)
  Models[[i]] = dredge(glmmTMB(Results2)[i] ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data= Results2, family = fam3[i])
  names(Models)[[i]] = colnames(Results2[i])
}


t = Results2[,i]
glmmTMB(t ~ Results2$Dist_edge_std + Results2$Dist_trail_std + Results2$Dist_trail_beginning_std + (1 | Results2$ForestID), family = gaussian)



## Testing expression separation

cols = list() 
modelnames = list()
for(i in 8:31){
  cat(i)
  cols[[i]] = colnames(Results2)[i]
  modelnames[[i]] = paste( cols[[i]]," ~ Dist_edge_std + $Dist_trail_std + $Dist_trail_beginning_std + (1 | $ForestID), data= Results2, family =" , fam3[i])
  #names(modelnames)[[i]] = colnames(Results2[i])
}


modelnames2 <- options(modelnames)
modelnames[[9]]
glmmTMB(modelnames[[9]])

modelnames

dredge(modelnames[[9]])


cols[[8]]
modelnames[[8]]
t = Results2[,i]
glmmTMB(t ~ Results2$Dist_edge_std + Results2$Dist_trail_std + Results2$Dist_trail_beginning_std + (1 | Results2$ForestID), family = gaussian)

dredge(glmmTMB(TAlphaAll  ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data= Results2, family = poisson))



gmteste <- colnames(Results2)[8]
gmteste

 
## Making separate datasets for each variable

require(MuMIn)
data(Cement)
d <- data.frame(Cement)
idx <- seq(11,13)
avgmod.95p <- list()
for (i in 1:length(idx)){
  d2 <- d[1:idx[i],]
  fm1 <- lm(y ~ ., data = d2)
  dd <- dredge(fm1, extra = c("R^2", F = function(x)
    summary(x)$fstatistic[[1]]))
  
  # 95% confidence set:
  confset.95p <- get.models(dd, cumsum(weight) <= .95)
  avgmod.95p[[i]] <- model.avg(confset.95p)
  
  
  ## Creating table
  
  Models = list() 
  for(i in 8:31){
    newTable = Results2[,c(i,1,5,6,7)]
    colnames(newTable)[1] = "y"
    Models[[i]] = dredge(glmmTMB(y ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data= newTable, family = fam3[i]))
    names(Models)[[i]] = colnames(Results2[i])
  }
  
  
  t = Results2[,i]
  glmmTMB(t ~ Results2$Dist_edge_std + Results2$Dist_trail_std + Results2$Dist_trail_beginning_std + (1 | Results2$ForestID), family = gaussian)
  
  
Models
str(Models)  
 Models
 
 Model9 <- Models[[9]]
 Model9
 str(Model9)
 models.df <- data.frame(t(matrix(unlist(Model9), nrow=length(Model9),  byrow=T)),stringsAsFactors = FALSE)
 
 str(models.df)
 
models.df2 <- apply(Models[[9]], 2 , as.numeric)
models.df2



 subListExtract(Models, name, simplify = FALSE, keep.names = TRUE)
 
valid.models <- Models[Models$AICc <= 2 ] 
  
valid.models
