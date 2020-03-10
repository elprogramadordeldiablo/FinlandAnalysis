

# Name: Rui Miguel Carvalho
# Date of creation: 7/2/2020



# Libraries -----------
library(BAT)
library(readr)
library(FD) 
library(alphahull)
library(hypervolume)
library(car)
library(MASS)
library(lme4)
#source(file = "HighstatLibV10.R") #cool tools to support
#library(factoextra) # Useful for PCA analysis
library(here)
library(data.table) # to work with data
library(dplyr)      # to manage data
library(magrittr)   # to use the pipe operator %>% 
library(MuMIn)
library(glmmTMB)
library(bbmle)
library(performance)
library(see)
library(pscl)
library(DHARMa)
library(rlang)
library(sjPlot)



# Load files -----------

#Ficheiro com as variável de distâncias com edge, trilhos e etc - para o GLMM

vars.veg <- read.csv2(here("data.veg","2017.ta.veg.var.csv"), row.names=1, header=TRUE,  stringsAsFactors = T, dec = ".", sep = ",")
str(vars.veg)
vars.veg$Dist_trail <- as.numeric(vars.veg$Dist_trail)
vars.veg$Dist_edge <- as.numeric(vars.veg$Dist_edge)
vars.veg$Dist_trail_beginning <- as.numeric(vars.veg$Dist_trail_beginning)

vars.veg$Dist_trail_std <- scale(vars.veg$Dist_trail, center = F)
vars.veg$Dist_edge_std <- scale(vars.veg$Dist_edge, center = F)
vars.veg$Dist_trail_beginning_std <- scale(vars.veg$Dist_trail_beginning, center = F)



#---- Ficheiro com as abundâncias por area de amostragem, para todas as amostras
all.area.2017 <- read.csv2(here("data.veg","2017.ta.all.csv"),row.names=1, header = TRUE, sep = ",")
str(all.area.2017)
species.veg <- colnames(all.area.2017)    ##Estava assim originalmente, mas parece-me trocado
sites.veg <- row.names(all.area.2017)
sites.veg

# 
Invasives <- read.csv2(here("data","Invasives.csv"),row.names=1, header = TRUE, sep = ",")
str(Invasives)
species.inv <- colnames(Invasives)    ##Estava assim originalmente, mas parece-me trocado
sites.inv <- row.names(Invasives)
sites.inv

#Taxonomical Alpha
Alpha_All.2017 <- alpha(all.area.2017)
colnames(Alpha_All.2017) = "Alpha_ll" 

Alpha_invasives <- alpha(Invasives)
colnames(Alpha_invasives) = "Alpha_inv" 

###### Taxonomical Beta
beta.all.ta.2017 <- beta(all.area.2017, abund=TRUE)
BetaInvasives2017 <- beta(Invasives, abund=TRUE)

#### Separating Total, Richness and Replacement  TAXONOMICAL betas into different data frames ----
betatotal.all.ta.2017 <- data.frame(as.matrix(beta.all.ta.2017[["Btotal"]]), row.names= sites.veg, sep= "tab")
colnames(betatotal.all.ta.2017) <- sites.veg

betarich.all.ta.2017 <- data.frame(as.matrix(beta.all.ta.2017[["Brich"]]), row.names= sites.veg,sep= "tab")
colnames(betarich.all.ta.2017) <- sites.veg

betarepl.all.ta.2017 <- data.frame(as.matrix(beta.all.ta.2017[["Brepl"]]),row.names= sites.veg, sep= "tab")
colnames(betarepl.all.ta.2017) <- sites.veg

# Beta partitions ----

betatotal.inv.2017 <- data.frame(as.matrix(BetaInvasives2017[["Btotal"]]), row.names= sites.inv, sep= "tab")
 colnames(betatotal.inv.2017) <- sites.inv

betarich.inv.2017 <- data.frame(as.matrix(BetaInvasives2017[["Brich"]]), row.names= sites.inv,sep= "tab")
colnames(betarich.inv.2017) <- sites.inv

betarepl.inv.2017 <- data.frame(as.matrix(BetaInvasives2017[["Brepl"]]),row.names= sites.inv, sep= "tab")
colnames(betarepl.inv.2017) <- sites.inv

# Separating the TAXONOMICAL beta values between the Control 250/Max and the other sampling areas from each trail ----

AA1v <- betatotal.all.ta.2017[1,1]
BB1v <- betatotal.all.ta.2017[5,c(2:5)]
CC1v <- betatotal.all.ta.2017[10,c(6:10)]
DD1v <- betatotal.all.ta.2017[14,c(11:14)]
all.tax.btotal.2017 <- c(AA1v, BB1v,  CC1, DD1v)

AA2v <- betarich.all.ta.2017[1,1]
BB2v <- betarich.all.ta.2017[5,c(2:5)]
CC2v <- betarich.all.ta.2017[10,c(6:10)]
DD2v <- betarich.all.ta.2017[14,c(11:14)]
all.tax.brich.2017 <- c(AA2v, BB2v,  CC2v, DD2v)

AA3v <- betarepl.all.ta.2017[1,1]
BB3v <- betarepl.all.ta.2017[5,c(2:5)]
CC3v <- betarepl.all.ta.2017[10,c(6:10)]
DD3v <- betarepl.all.ta.2017[14,c(11:14)]
all.tax.brepl.2017 <- c(AA3v, BB3v,  CC3v, DD3v)

#Invasives partitioning

AA1vi <- betatotal.inv.2017[1,1]
BB1vi <- betatotal.inv.2017[5,c(2:5)]
CC1vi <- betatotal.inv.2017[10,c(6:10)]
DD1vi <- betatotal.inv.2017[14,c(11:14)]
inv.tax.btotal.2017 <- c(AA1vi, BB1vi,  CC1, DD1vi)

AA2vi <- betarich.inv.2017[1,1]
BB2vi <- betarich.inv.2017[5,c(2:5)]
CC2vi <- betarich.inv.2017[10,c(6:10)]
DD2vi <- betarich.inv.2017[14,c(11:14)]
inv.tax.brich.2017 <- c(AA2vi, BB2vi,  CC2vi, DD2vi)

AA3vi <- betarepl.inv.2017[1,1]
BB3vi <- betarepl.inv.2017[5,c(2:5)]
CC3vi <- betarepl.inv.2017[10,c(6:10)]
DD3vi <- betarepl.inv.2017[14,c(11:14)]
inv.tax.brepl.2017 <- c(AA3vi, BB3vi,  CC3vi, DD3vi)
#
#### COMPILING ALL BETA INFORMATION INTO A TABLE, AND EXPORTING IT TO A FILE ----

betas.veg <- as.data.frame(t(rbind(all.tax.btotal.2017,all.tax.brich.2017,all.tax.brepl.2017,inv.tax.btotal.2017, inv.tax.brich.2017, inv.tax.brepl.2017 )))

betas.veg[betas.veg == 0] <- 0.001
betas.veg[betas.veg == 1] <- 0.999
str(betas.veg)

###########################################################################
#                               RESULTS                                ####
###########################################################################

Results.veg <- cbind.data.frame(vars.veg, Alpha_All.2017, Alpha_invasives, betas.veg)

#MAKING THE RESULTS EXPORTABLE INTO CSV
Results.veg <- apply(Results.veg, 2 , as.character, header=TRUE)

#NAMING THE TRAIL SEGMENTS
rownames(Results.veg) <- rownames(betas.veg)
Results.veg <- Results.veg[-1,]

#PASSING RESULTS TO FILE
write.csv(Results.veg, file = here("results","RESULTS.VEG.csv"), row.names = TRUE)

Results.veg.import <- read.csv2(here("results","RESULTS.VEG.CSV"), header=TRUE, row.names = 1,  stringsAsFactors = T,sep = ",", dec = ".")
Results.veg.import

###########################################################################
#                           Generating models                          ####
###########################################################################

fam.veg <- c(NA, NA, NA, NA,NA, NA,NA, "poisson","beta_family","beta_family", "beta_family", "beta_family","beta_family", "beta_family")

## Alpha, abundances and proportions

dredge.alpha <- dredge(glmmTMB(Alpha_ll ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data=Results.veg.import , family = "poisson"))
dredge.alpha <- data.frame(dredge.alpha)
write.csv(dredge.alpha, file = here("results","test.alpha.VEG.csv"), row.names = TRUE)


dredge.alpha.inv <- dredge(glmmTMB(Alpha_inv ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data=Results.veg.import , family = "poisson"))
dredge.alpha.inv <- data.frame(dredge.alpha.inv)
write.csv(dredge.alpha.inv, file = here("results","test.alpha.inv.csv"), row.names = TRUE)

## Betas
withoutcontrols.veg <- Results.veg.import[-c(1,4,19,13),]

#models.beta.veg = data.frame()
Models.veg = list() 
for(i in 10:15){
  newTable.veg = withoutcontrols.veg[,c(i,1,5,6,7)]
  colnames(newTable.veg)[1] = "y"
  Models.veg[[i]] = dredge(glmmTMB(y ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data= newTable.veg, family = fam.veg[i]))
  names(Models.veg)[[i]] = colnames(withoutcontrols.veg[i])
  Models.veg[[i]] = apply(Models.veg[[i]], 2 , as.numeric)
  Models.veg[[i]] = cbind(rep(colnames(withoutcontrols.veg)[i] , nrow(Models.veg[[i]])), Models.veg[[i]])
  colnames(models.beta.veg) = colnames(as.data.frame(Models.veg[[i]]))
  models.beta.veg = rbind(models.beta.veg, as.data.frame(Models.veg[[i]]))
  #for(j in 2:ncol(models.beta.veg)){
  #  models.beta.veg[,j] <- as.numeric(models.beta.veg[,j])
  #}
}

write.csv(models.beta.veg, file = here("results","models.beta.veg.csv"), row.names = TRUE)

Models.veg.total <- data.frame(Models.veg$all.tax.btotal.2017)
Models.veg.rich <- data.frame(Models.veg$all.tax.brich.2017)
Models.veg.repl <- data.frame(Models.veg$all.tax.brepl.2017)

dredge.veg <- rbind(Models.veg.total,Models.veg.rich,Models.veg.repl)
write.csv(dredge.veg, file = here("results","test.betas.VEG.csv"), row.names = TRUE)

####For Betas ----

aic.weights.betas.veg = data.frame()
models.df5 = data.frame()
Models.veg2 = list() 
for(i in 10:15){
  newTable = Results2[,c(i,1,5,6,7)]
  colnames(newTable)[1] = "y"
  Models.veg2[[i]] = dredge(glmmTMB(y ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data= newTable, family = fam3[i]))
  names(Models.veg2)[[i]] = colnames(Results.veg[i])
  Models.veg2[[i]] = apply(Models.veg2[[i]], 2 , as.numeric)
  #############################
  
  w = c()
  for(k in 1:3){
    w[k] = sum(Models.veg2[[i]][!is.na(Models.veg2[[i]][,(k+2)]),10])
  }
  aic.weights.betas.veg = rbind(aic.weights.betas.veg, w)
  
  #############################
  Models.veg2[[i]] = cbind(rep(colnames(Results.veg)[i], nrow(Models.veg2[[i]])), Models.veg2[[i]])
  models.df3 = rbind(models.df3, as.data.frame(Models.veg2[[i]]))
  
  #for(j in 2:ncol(models.df)){
  #  models.df[,j] <- as.numeric(models.df[,j])
  #}
}
colnames(aic.weights.betas) = colnames(models.df)[4:6]
rownames(aic.weights.betas) = unique(models.df2[,1])
models.df2
aic.weights.betas

model.weights <- rbind(aic.weights1,aic.weights.betas)

write.csv(model.weights, file = here("results","aic.weights.csv"), row.names = TRUE)

#selecting delta <2

AICmodels <- filter(dredge.models, delta < 2)
str(AICmodels)

write.csv(AICmodels, file = here("results","aic.models.csv"), row.names = TRUE)

###########################################################################
#                     Examining each model                            #####
###########################################################################

###### 1  - Crating functions to aid its interpertation

### Overdispersion function
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

## Making a function to extract all the necessary outputs out of each model

model.output <- function(model) {
  
  print(" ")
  print("CHECKING COLINEARITY")
  print(" ")
  #print(check_conversion(model))
  print(" ")
  print("CHECKING COLINEARITY")
  print(" ")
  print(check_collinearity(model))
  print(" ")
  print("TESTING FOR OVERDISPERSION - P-VALUES < 0.05 meand overdispersion")
  print(" ")
  #print(check_overdispersion(model))
  print(" ")
  print(overdisp_fun(model))
  print(" ")
  print("STANDARD R2")
  print(" ")
  print(performance::r2(model))
  print(" ")
  print("NAKAGAWA R2")
  print(" ")
  print(performance::r2_nakagawa(model))
  print(" ")
  print("MODEL SUMMARY")
  print(" ")
  print(summary(model))
  print(" ")
  print("RANDOM EFFECTS")
  print(" ")
  print(tab_model(model))
  print(" ")
  print("CHECKING MODEL GRAPHICALLY")
  print(" ")
  print(check_model(model))
  print(" ")
  print("CHECKING MODEL GRAPHICALLY")
  print(" ")
  print(check_model(model))
  
}

model.output.poisson <- function(model) {
  
  print(" ")
  print("CHECKING COLINEARITY")
  print(" ")
  #print(check_conversion(model))
  print(" ")
  print("CHECKING COLINEARITY")
  print(" ")
  print(check_collinearity(model))
  print(" ")
  print("TESTING FOR OVERDISPERSION - P-VALUES < 0.05 meand overdispersion")
  print(" ")
  print(check_overdispersion(model))
  print(" ")
  print(overdisp_fun(model))
  print(" ")
  print("STANDARD R2")
  print(" ")
  print(performance::r2(model))
  print(" ")
  print("NAKAGAWA R2")
  print(" ")
  print(performance::r2_nakagawa(model))
  print(" ")
  print("MODEL SUMMARY")
  print(" ")
  print(summary(model))
  print(" ")
  print("RANDOM EFFECTS")
  print(" ")
  print(tab_model(model))
  print(" ")
  print("CHECKING MODEL GRAPHICALLY")
  print(" ")
  print(check_model(model))
  print(" ")
  print("CHECKING MODEL GRAPHICALLY")
  print(" ")
  print(check_model(model))
  
}

## Obtaining the several R2 values
###############
###CRUDE SOLUTION FOR CALCULATING R 
##################3

version1 <- function(m) {
  lmfit <-  lm(model.response(model.frame(m)) ~ fitted(m))
  summary(lmfit)$r.squared
  
}

## ??20 (Xu 2003), which is almost the same, is based on comparing the residual variance of the full model against the 
##residual variance of a (fixed) intercept-only null model:

version2 <- function(m){
  
  1-var(residuals(m))/var(model.response(model.frame(m)))
  
}

## Another possibility is the squared correlation between the response variable and the predicted values:
version3 <- function(m){
  cor(model.response(model.frame(m)),predict(m,type="response"))^2
}

r2 <- function(model){
  
  print(" ")
  print("STANDARD R2")
  print(" ")
  print(performance::r2(model))
  print(" ")
  print(" ")
  print("NAKAGAWA")
  print(" ")
  print(performance::r2_nakagawa(model))
  print(" ")
  print("VERSION1")
  print(" ")
  print(version1(model))
  print(" ")
  print("VERSION2")
  print(" ")
  print(version2(model))
  # print(" ")
  # print("VERSION3")
  # print(" ")
  # print(version3(model))
} 

