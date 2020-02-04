library(see)
library(pscl)

# Loading the excel file with the formulas

model.formulas <- read.csv2(here("data","model.formulas.csv"), header=TRUE,  stringsAsFactors = T, dec = ".")
model.formulas


model.lala = list()
model.lala <- lapply(model.formulas[,1], model.output)

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
  print("CHECKING MODEL GRAPHICALLY")
  print(" ")
  print(check_model(model))
  print(" ")
  print("CHECKING MODEL GRAPHICALLY")
  print(" ")
  print(check_model(model))
  
}
return(model.output(TAlphaAll.glmm))
## Applying the extraction formula to all models

model.results = list()
for (i in unique(model.formulas[,"model.name"])) {
  aaaa <- return(model.output[[i]])
  model.results = c(model.results, aaaa)
}

model.results
str(model.results)



check_convergence(abund.nat.glmm.12)


# for(i in unique(Alpha_controls[,"Trail_Sampling.Area"])){
#   Alpha_controlsMatrix <- Alpha_controls[Alpha_controls[,"Trail_Sampling.Area"]== i, -1]  #definir a matriz a analisar
#   alphaaccumlist[[i]] <- (alpha.accum(Alpha_controlsMatrix))
#   alphaaccumlist[[i]]
#   
#   #plot(alphaaccumlist[[i]][,1], alphaaccumlist[[i]][,2], xlab="Samples", ylab="Individuals", )
#   plot(alphaaccumlist[[i]][,1], alphaaccumlist[[i]][,3], xlab="Samples", ylab="Individuals", 
#        ylim=c(0,50), main = i, col="blue", type = "p")
#   
#   lines(alphaaccumlist[[i]][,1], alphaaccumlist[[i]][,4])
#   lines(alphaaccumlist[[i]][,1], alphaaccumlist[[i]][,5], col="yellow")
#   lines(alphaaccumlist[[i]][,1], alphaaccumlist[[i]][,8], col="red")
#   lines(alphaaccumlist[[i]][,1], alphaaccumlist[[i]][,13], col="red", type = "p")
#   lines(alphaaccumlist[[i]][,1], alphaaccumlist[[i]][,17], col="green")
#   lines(alphaaccumlist[[i]][,1], alphaaccumlist[[i]][,19], col="green", type = "p")






TAlphaAll.glmm = glmmTMB(TAlphaAll ~    (1 | ForestID), data= Results2,family = "poisson") 
summary(TAlphaAll.glmm)
performance::r2(TAlphaAll.glmm)
overdisp_fun(model1.test)
model.output(TAlphaAll.glmm)
bb <- model.output(TAlphaAll.glmm)

lalala <- list()

lalala <- c(lalala, aa ,bb)


TAlphaAll.glmm.1 = glmmTMB(TAlphaAll ~    (1 | ForestID), data= Results2,family = "poisson") 
model.output(TAlphaAll.glmm.1) 

TAlphaNat.glmm.2 = glmmTMB(TAlphaNat ~    (1 | ForestID), data= Results2,family = poisson(link="log")) 
model.output(TAlphaNat.glmm.2) 
check_model(TAlphaNat.glmm.2)
TAlphaNInd.glmm.3 = glmmTMB(TAlphaNInd ~    (1 | ForestID), data= Results2,family = "binomial") 
model.output(TAlphaNInd.glmm.3) 

TAlphaNInd.glmm.4 = glmmTMB(TAlphaNInd ~  Dist_trail_beginning_std + (1 | ForestID), data= Results2,family = "poisson") 
TAlphaNInd.glmm.5 = glmmTMB(TAlphaNInd ~   Dist_trail_std +(1 | ForestID), data= Results2,family = "poisson") 
TAlphaNInd.glmm.6 = glmmTMB(TAlphaNInd ~ Dist_edge_std +  (1 | ForestID), data= Results2,family = "poisson") 
FAlphaAll.glmm.7 = glmmTMB(FAlphaAll ~    (1 | ForestID), data= Results2,family = "poisson") 
FAlphaAll.glmm.8 = glmmTMB(FAlphaAll ~  Dist_trail_beginning_std + (1 | ForestID), data= Results2,family = Gamma(link="log")) 

FAlphaNat.glmm.9 = glmmTMB(FAlphaNat ~    (1 | ForestID), data= Results2,family = "Gamma") 
FAlphaNInd.glmm.10 = glmmTMB(FAlphaNInd ~    (1 | ForestID), data= Results2,family = "Gamma") 


abund.all.glmm.11 = glmmTMB(abund.all ~ Dist_edge_std +Dist_trail_beginning_std +Dist_trail_std +(1 | ForestID), data= Results2,family = nbinom2(link = "log")) 
model.output(abund.all.glmm.11) 
check_model(abund.all.glmm.11)
performance::check_conversion(abund.all.glmm.11)
performance::check_convergence(abund.all.glmm.11)
check_distribution(abund.all.glmm.11)
plot(check_distribution(abund.all.glmm.11))


## testing different distr families
abund.all.glmm.11.nb2 = glmmTMB(abund.all ~ Dist_edge_std +Dist_trail_beginning_std + Dist_trail_std + (1 | ForestID), data= Results2,family = nbinom2(link = "log")) 
model.output(abund.all.glmm.11) 


abund.all.glmm.11.nb1 = glmmTMB(abund.all ~ Dist_edge_std +Dist_trail_beginning_std + Dist_trail_std + (1 | ForestID), data= Results2,family = nbinom1) 
model.output(abund.all.glmm.11) 


abund.all.glmm.11.poiss = glmmTMB(abund.all ~ Dist_edge_std +Dist_trail_beginning_std + Dist_trail_std + (1 | ForestID), data= Results2,family = "poisson") 
model.output(abund.all.glmm.11) 

odTest(abund.all.glmm.11)


check_overdispersion(abund.all.glmm.11)
######testing another function for glmm
abund.all.glmer2.11 = glmer(abund.all ~ Dist_edge_std +Dist_trail_beginning_std +Dist_trail_std +(1 | ForestID), data= Results2,family = "poisson") 
performance::check_convergence(abund.all.glmer2.11)

model.output(abund.all.glmer2.11) 
Residuals.8 <- resid(abund.all.glmer2.11)
check_collinearity(abund.all.glmer2.11)
classify_distribution(TAlphaAll)
R.squared.GLMM(abund.all.glmer.11)
performance::r2(FAlphaAll.glmm.8)
R.squared.GLMM(abund.all.glmer2.11)
performance::r2(abund.all.glmer2.11)

abund.nat.glmm.12 = glmmTMB(abund.nat ~ Dist_edge_std +Dist_trail_beginning_std +Dist_trail_std +(1 | ForestID), data= Results2,family = nbinom1) 
abund.nind.glmm.13 = glmmTMB(abund.nind ~ Dist_edge_std +  (1 | ForestID), data= Results2,family = nbinom1) 
abund.nind.glmm.14 = glmmTMB(abund.nind ~   Dist_trail_std +(1 | ForestID), data= Results2,family = "poisson") 
prop.Talpha.glmm.15 = glmmTMB(prop.Talpha ~    (1 | ForestID), data= Results2,family = "beta_family") 
prop.Talpha.glmm.16 = glmmTMB(prop.Talpha ~   Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
prop.Falpha.glmm.17 = glmmTMB(prop.Falpha ~    (1 | ForestID), data= Results2,family = "beta_family") 
prop.Falpha.glmm.18 = glmmTMB(prop.Falpha ~  Dist_trail_beginning_std + (1 | ForestID), data= Results2,family = "beta_family") 
prop.abund.glmm.19 = glmmTMB(prop.abund ~    (1 | ForestID), data= Results2,family = "beta_family") 
all.tax.btotal.glmm.20 = glmmTMB(all.tax.btotal ~   Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
all.tax.btotal.glmm.21 = glmmTMB(all.tax.btotal ~ Dist_edge_std + Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
all.tax.brich.glmm.22 = glmmTMB(all.tax.brich ~   Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
all.tax.brich.glmm.23 = glmmTMB(all.tax.brich ~ Dist_edge_std + Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
all.tax.brepl.glmm.24 = glmmTMB(all.tax.brepl ~   Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
all.tax.brepl.glmm.25 = glmmTMB(all.tax.brepl ~ Dist_edge_std + Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
nat.tax.btotal.glmm.26 = glmmTMB(nat.tax.btotal ~   Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
nat.tax.btotal.glmm.27 = glmmTMB(nat.tax.btotal ~ Dist_edge_std + Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
nat.tax.brich.glmm.28 = glmmTMB(nat.tax.brich ~    (1 | ForestID), data= Results2,family = "beta_family") 
nat.tax.brepl.glmm.29 = glmmTMB(nat.tax.brepl ~   Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
nat.tax.brepl.glmm.30 = glmmTMB(nat.tax.brepl ~ Dist_edge_std + Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
nind.tax.btotal.glmm.31 = glmmTMB(nind.tax.btotal ~ Dist_edge_std + Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
nind.tax.btotal.glmm.32 = glmmTMB(nind.tax.btotal ~    (1 | ForestID), data= Results2,family = "beta_family") 
nind.tax.btotal.glmm.33 = glmmTMB(nind.tax.btotal ~   Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
nind.tax.brich.glmm.34 = glmmTMB(nind.tax.brich ~ Dist_edge_std + Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
nind.tax.brepl.glmm.35 = glmmTMB(nind.tax.brepl ~    (1 | ForestID), data= Results2,family = "beta_family") 
nind.tax.brepl.glmm.36 = glmmTMB(nind.tax.brepl ~ Dist_edge_std +  (1 | ForestID), data= Results2,family = "beta_family") 
all.func.btotal.glmm.37 = glmmTMB(all.func.btotal ~   Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
all.func.btotal.glmm.38 = glmmTMB(all.func.btotal ~ Dist_edge_std + Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
all.func.brich.glmm.39 = glmmTMB(all.func.brich ~   Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
all.func.brich.glmm.40 = glmmTMB(all.func.brich ~ Dist_edge_std + Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
all.func.brepl.glmm.41 = glmmTMB(all.func.brepl ~   Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
all.func.brepl.glmm.42 = glmmTMB(all.func.brepl ~ Dist_edge_std + Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
nat.func.btotal.glmm.43 = glmmTMB(nat.func.btotal ~   Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
nat.func.btotal.glmm.44 = glmmTMB(nat.func.btotal ~ Dist_edge_std + Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
nat.func.brich.glmm.45 = glmmTMB(nat.func.brich ~   Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
nat.func.brich.glmm.46 = glmmTMB(nat.func.brich ~ Dist_edge_std + Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
nat.func.brepl.glmm.47 = glmmTMB(nat.func.brepl ~   Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
nind.func.btotal.glmm.48 = glmmTMB(nind.func.btotal ~ Dist_edge_std + Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
nind.func.brich.glmm.49 = glmmTMB(nind.func.brich ~ Dist_edge_std + Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
nind.func.brepl.glmm.50 = glmmTMB(nind.func.brepl ~    (1 | ForestID), data= Results2,family = "beta_family") 
nind.func.brepl.glmm.51 = glmmTMB(nind.func.brepl ~ Dist_edge_std +  (1 | ForestID), data= Results2,family = "beta_family") 
### glmer alternative
nind.func.brepl.glmer.51 = glmer(nind.func.brepl ~ Dist_edge_std +  (1 | ForestID), data= Results2,family = "beta_family", control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1000000))) 
nat.func.brich.glmer.45 = glmer(nat.func.brich ~   Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 

pdf('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA.pdf')
model.output.poisson(TAlphaAll.glmm.1) 
model.output.poisson(TAlphaNat.glmm.2) 
model.output.poisson(TAlphaNInd.glmm.3) 
model.output.poisson(TAlphaNInd.glmm.4) 
model.output.poisson(TAlphaNInd.glmm.5) 
model.output.poisson(TAlphaNInd.glmm.6) 
model.output(FAlphaAll.glmm.7) 
model.output(FAlphaAll.glmm.8) 
model.output(FAlphaNat.glmm.9) 
model.output(FAlphaNInd.glmm.10) 
model.output.poisson(abund.all.glmm.11) 
model.output(abund.nat.glmm.12) 
model.output(abund.nind.glmm.13) 
model.output(abund.nind.glmm.14) 
model.output(prop.Talpha.glmm.15) 
model.output(prop.Talpha.glmm.16) 
model.output(prop.Falpha.glmm.17) 
model.output(prop.Falpha.glmm.18) 
model.output(prop.abund.glmm.19) 
model.output(all.tax.btotal.glmm.20) 
model.output(all.tax.btotal.glmm.21) 
model.output(all.tax.brich.glmm.22) 
model.output(all.tax.brich.glmm.23) 
model.output(all.tax.brepl.glmm.24) 
model.output(all.tax.brepl.glmm.25) 
model.output(nat.tax.btotal.glmm.26) 
model.output(nat.tax.btotal.glmm.27) 
model.output(nat.tax.brich.glmm.28) 
model.output(nat.tax.brepl.glmm.29) 
model.output(nat.tax.brepl.glmm.30) 
model.output(nind.tax.btotal.glmm.31) 
model.output(nind.tax.btotal.glmm.32) 
model.output(nind.tax.btotal.glmm.33) 
model.output(nind.tax.brich.glmm.34) 
model.output(nind.tax.brepl.glmm.35) 
model.output(nind.tax.brepl.glmm.36) 
model.output(all.func.btotal.glmm.37) 
model.output(all.func.btotal.glmm.38) 
model.output(all.func.brich.glmm.39) 
model.output(all.func.brich.glmm.40) 
model.output(all.func.brepl.glmm.41) 
model.output(all.func.brepl.glmm.42) 
model.output(nat.func.btotal.glmm.43) 
model.output(nat.func.btotal.glmm.44) 
model.output(nat.func.brich.glmm.45) 
model.output(nat.func.brich.glmm.46) 
model.output(nat.func.brepl.glmm.47) 
model.output(nind.func.btotal.glmm.48) 
model.output(nind.func.brich.glmm.49) 
model.output(nind.func.brepl.glmm.50) 
model.output(nind.func.brepl.glmm.51) 
dev.off()



###https://www.r-bloggers.com/generalized-linear-models-understanding-the-link-function/
### está apresnetado algum código para a funç\ao predict