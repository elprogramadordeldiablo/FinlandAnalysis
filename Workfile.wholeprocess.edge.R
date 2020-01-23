####ABORTED - CANNOT TRAT A CONTINUOUSVARIABLE AS A RANDOM EFFECT###


# Name: Rui Miguel Carvalho
# Date of creation: 10/13/2012
# Date of last update: 
# 
# 

# Index --------------------

# CALC COMPLETENESS - Faz as curvas de acumulação para os controlos 250
#  SETTTING FILES
#  CALC_COMPLETENESS
#  CALC_TD_ALPHA
#  CALC_TD_BETA
#  CALC_FD_ALPHA
# CALC_FD_BETA
# Summing up Alpha values - 162
# Summing up Beta values - 336
# RESULTS - 427
# Data Exploration - 451
# GLMM - 478


# Libraries -----------------------------------------------
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



# Load files --------------------

#Ficheiro com as variável de distâncias com edge, trilhos e etc - para o GLMM

Variables <- read.csv2(here("data","GLMM_Variables.csv"), row.names=1, header=TRUE,  stringsAsFactors = T, dec = ".")
Variables$Dist_trail <- as.numeric(Variables$Dist_trail)
Variables$Dist_edge <- as.numeric(Variables$Dist_edge)
Variables$Dist_trail_beginning <- as.numeric(Variables$Dist_trail_beginning)

Variables$Dist_trail_std <- scale(Variables$Dist_trail, center = F)
Variables$Dist_edge_std <- scale(Variables$Dist_edge, center = F)
Variables$Dist_trail_beginning_std <- scale(Variables$Dist_trail_beginning, center = F)
str(Variables)

#Ficheiro com as variável de distâncias com edge, trilhos e etc - para o GLMM

abundances <- read.csv2(here("data","abundances.csv"), row.names=1, header=TRUE,  stringsAsFactors = T, dec = ".")

# Ficheiro com as abundâncias por amostra para os controlos 250/Max.
Alpha_controls <- read.csv2(here("data","Control250_Fin.csv"), header=TRUE,  stringsAsFactors = T, dec = ".")
str(Alpha_controls)


# Matrix with species (rows)x traits (cols) for all species in all plots


Traits  <- read.csv2(here("data","Traits_All_Fin.csv"), row.names=1, header=TRUE,  stringsAsFactors = T, dec = ".")
Traits <- Traits[,-c(1:3)] 
species <- rownames(Traits)
str(Traits)


# Matrix with species (rows)x traits (cols) for all natives species in all plots

TraitsNat  <- read.csv2(here("data","Traits_Nat_Fin.csv"), row.names=1, header=TRUE,  stringsAsFactors = T, dec = ".")
TraitsNat <- TraitsNat[,-c(1:3)] 
speciesNat <- rownames(TraitsNat)
str(Traits)

# Matrix with species (rows)x traits (cols) for all Non-Indigenous species in all plots

TraitsNInd  <- read.csv2(here("data","Traits_NInd_Fin.csv"), row.names=1, header=TRUE,  stringsAsFactors = T, dec = ".")
TraitsNInd <- TraitsNInd[,-c(1:3)] 
SpeciesNInd <- rownames(TraitsNInd)
str(TraitsNInd)


# Ficheiro com as abundâncias por Área de amostragem, para todas as amostras
SAAll <- read.csv2(here("data","All_Fin.csv"),row.names=1, header = TRUE)
colnames(SAAll) <- species   ##Estava assim originalmente, mas parece-me trocado
sites <- row.names(SAAll)
sites


# HEllinger transformation
#SAAll <- decostand(SAAll, "hellinger") 
#SAAll

# Ficheiro com as abundâncias por área de amostragem, para as espécies ind�???genas

SANat <- read.csv2(here("data","Nat_Fin.csv"),row.names=1, header = TRUE)

colnames(SANat) <- speciesNat 
# HEllinger transformation
#SANat <- decostand(SANat, "hellinger") 
#SANat


# Ficheiro com as abundâncias por área de amostragem, para as espécies não-ind�???genas
SANInd <- read.csv2(here("data","NInd_Fin.csv"),row.names=1, header = TRUE)
colnames(SANInd) <- SpeciesNInd

# HEllinger transformation
#SANInd <- decostand(SANInd, "hellinger") 
#SANInd

# Vector with weight trail and treatment names for MDS
MDSfile <- read.csv2(here("data","MDS_vectors.csv"))
trail <- MDSfile[,1]
trail
trail <- as.vector(trail)
treatment <- MDSfile[,2]
str(treatment)
treatment <- as.vector(treatment)
str(treatment)

# Vector with weight ratios for trails - used in cluster::daisy to compensate decomposition of one variable by
# multiple columns
Weightsfile <- read.csv2(here("data","Weight_Ratios_Traits.csv"), header = TRUE, dec=".")
weights <- Weightsfile[,-c(1:4)]
weights <- as.vector(weights)

Weightsfile2 <- read.csv2(here("data","Weight_Ratios_Traits.csv"), header = TRUE, dec=".")


fam <- read.csv2(here("results","RESULTS copy.csv"), header=TRUE, row.names = 1,  stringsAsFactors = T, sep = ";", dec = ".")
fam <- fam[1,]


fam2 <- read.csv2(here("data","glmm.families.csv"), header=TRUE, row.names = 1,  stringsAsFactors = T, sep = ",", dec = ".")
fam2


# Looking into the datasets -------------------------------

# Ver as estruturas dos dados (garantir os factors e os valores numericos) 
str(Variables)
str(Alpha_controls)
str(Traits)
str(TraitsNat)
str(TraitsNInd)
str(SAAll)
str(SANat)
str(SANInd)
str(weights)
str(treatment)
# CALCULATING COMPLETENESS FOR THE CONTROLS FOR BETA ANALYSYS ----

Alpha_controls <- read.csv2(here("data","Control250_Fin.csv"), header = TRUE)

alphaaccumlist <-list()
alphaaccumlist
plotlist <- list()
par(mfrow=c(2,3))
pdf('AccumCurvesControls.pdf')
for(i in unique(Alpha_controls[,"Trail_Sampling.Area"])){
  Alpha_controlsMatrix <- Alpha_controls[Alpha_controls[,"Trail_Sampling.Area"]== i, -1]  #definir a matriz a analisar
  alphaaccumlist[[i]] <- (alpha.accum(Alpha_controlsMatrix))
  alphaaccumlist[[i]]
  
  #plot(alphaaccumlist[[i]][,1], alphaaccumlist[[i]][,2], xlab="Samples", ylab="Individuals", )
  plot(alphaaccumlist[[i]][,1], alphaaccumlist[[i]][,3], xlab="Samples", ylab="Individuals", 
       ylim=c(0,50), main = i, col="blue", type = "p")
  
  lines(alphaaccumlist[[i]][,1], alphaaccumlist[[i]][,4])
  lines(alphaaccumlist[[i]][,1], alphaaccumlist[[i]][,5], col="yellow")
  lines(alphaaccumlist[[i]][,1], alphaaccumlist[[i]][,8], col="red")
  lines(alphaaccumlist[[i]][,1], alphaaccumlist[[i]][,13], col="red", type = "p")
  lines(alphaaccumlist[[i]][,1], alphaaccumlist[[i]][,17], col="green")
  lines(alphaaccumlist[[i]][,1], alphaaccumlist[[i]][,19], col="green", type = "p")
}
dev.off()

#CALC_COMPLETENESS ----

# Quero usar o fichero dos contrlos para calcular os estimadores, e da�??? obter completeness.
# sitescontr <- Alpha_controls[,1]
# species_contr <- Alpha_controls[,1]
# rownames(Alpha_controls) <-sitescontr
# colnames(Alpha_controls) <- species_contr
# Alpha_controls <- Alpha_controls[,-1]
# alpha.estimate(Alpha_controls)





###   Calculating functional elements ----

## All species tree and FD Alpha
dissimAll <- cluster::daisy(Traits, "gower", weights = weights)  # Calculating distances between species
#      https://www.rdocumentation.org/packages/cluster/versions/2.0.7-1/topics/daisy


tree <- hclust(dissimAll, "average")     
par(mfrow=c(1,1))
dev.copy(device = jpeg, filename = 'Tree.jpeg', width = 1000, height = 500)
plot(tree, hang = -1)  ## PQ O HANG= - 1?
dev.off()
#plot(as.dendrogram(tree))

# All natives tree and FD Alpha----
dissimNat <- cluster::daisy(TraitsNat, "gower", weights = weights)  # Calculating distances between species
treeNat <- hclust(dissimNat, "average")          # building the dendrogrma

cor(dissimNat, cophenetic(treeNat))

par(mfrow=c(1,1))
dev.copy(device = jpeg, filename = 'TreeNat.jpeg', width = 1000, height = 500)
plot(treeNat, hang = -1)  ## PQ O HANG= - 1?
dev.off()
#  plot(as.dendrogram(treeNat))
str(TraitsNInd)

# All Non-Ind tree and FD Alpha ----
dissimNInd <- cluster::daisy(TraitsNInd, "gower", weights = weights)  # Calculating distances between species
treeNInd <- hclust(dissimNInd, "average")          # building the 
par(mfrow=c(1,1))
dev.copy(device = jpeg, filename = 'TreeNInd.jpeg', width = 1000, height = 500)
plot(treeNInd, hang = -1)  ## PQ O HANG= - 1?
dev.off()

# ALPHAs for All species, Natives and Non-Indigenous ----

#Taxonomical Alpha
Alpha_All <- alpha(SAAll)
Alpha_Nat<- alpha(SANat)
Alpha_NInd<- alpha(SANInd)

#Functional Alpha
FDalphaAll <- alpha(SAAll, tree)
FDalphaNat <- alpha(SANat, treeNat)
FDalphaNInd <- alpha(SANInd, treeNInd)

Alphas <- as.data.frame(cbind(Alpha_All, Alpha_Nat, Alpha_NInd, FDalphaAll, FDalphaNat, FDalphaNInd))
colnames(Alphas) <- cbind("TAlphaAll", "TAlphaNat", "TAlphaNInd","FAlphaAll", "FAlphaNat", "FAlphaNInd" )
Alphas # Will be printed along with BETAS in the RESULTS file


#Removing absolute zeros from functional results

Alphas$FAlphaNInd[Alphas$FAlphaNInd == 0] <- 0.0001

##Abundances
Alphas$abund.all=rowSums(SAAll[,1:39])
Alphas$abund.nat=rowSums(SANat[,1:21])
Alphas$abund.nind=rowSums(SANInd[,1:18])

## Proportions
Alphas$prop.Talpha=(Alphas$TAlphaNat / Alphas$TAlphaAll)
Alphas$prop.Talpha[Alphas$prop.Talpha == 1] <- 0.9999
Alphas$prop.Talpha[Alphas$prop.Talpha == 0] <- 0.0001

Alphas$prop.Falpha=Alphas$FAlphaNat / Alphas$FAlphaAll
Alphas$prop.Falpha[Alphas$prop.Falpha == 1] <- 0.9999
Alphas$prop.Falpha[Alphas$prop.Falpha == 0] <- 0.0001

Alphas$prop.abund=Alphas$abund.nat / Alphas$abund.all
Alphas$prop.abund[Alphas$prop.abund == 1] <- 0.9999
Alphas$prop.abund[Alphas$prop.abund == 0] <- 0.0001

Alphas

# Calculating Betas for All species, Natives and Non-Indigenous ----
# Beta results come as a list (class(BetaAll) = list), so we transform them into data frame to export them into csv

# Taxonomical Beta
BetaAll <- beta(SAAll)
BetaNat <- beta(SANat)
BetaNInd <- beta(SANInd) 

#Functional Beta
BetaFuncAll <- beta(SAAll, tree, abund =  T)
BetaFuncNat <- beta(SANat, treeNat, abund = T) 
BetaFuncNInd <- beta(SANInd, treeNInd, abund= T)


## Separating Total, Richness and Replacement  TAXONOMICAL betas into different data frames ----
BetaAllTotal <- data.frame(as.matrix(BetaAll[["Btotal"]]), row.names= sites, sep= "tab")
colnames(BetaAllTotal) <- sites
BetaAllRich <- data.frame(as.matrix(BetaAll[["Brich"]]), sep= "tab")
BetaAllRepl <- data.frame(as.matrix(BetaAll[["Brepl"]]), sep= "tab")
BetaNatTotal <- data.frame(as.matrix(BetaNat[["Btotal"]]), sep= "tab")
BetaNatRich <- data.frame(as.matrix(BetaNat[["Brich"]]), sep= "tab")
BetaNatRepl <- data.frame(as.matrix(BetaNat[["Brepl"]]), sep= "tab")
BetaNIndTotal <- data.frame(as.matrix(BetaNInd[["Btotal"]]), sep= "tab")
BetaNIndRich <- data.frame(as.matrix(BetaNInd[["Brich"]]), sep= "tab")
BetaNIndRepl <- data.frame(as.matrix(BetaNInd[["Brepl"]]), sep= "tab")

## Separating Total, Richness and Replacement  FUNCTIONAL betas into different data frames
BetaFuncAllTotal <- data.frame(as.matrix(BetaFuncAll[["Btotal"]]), sep= "tab")
BetaFuncAllRich <- data.frame(as.matrix(BetaFuncAll[["Brich"]]), sep= "tab")
BetaFuncAllRepl <- data.frame(as.matrix(BetaFuncAll[["Brepl"]]), sep= "tab")
BetaFuncNatTotal <- data.frame(as.matrix(BetaFuncNat[["Btotal"]]), sep= "tab")
BetaFuncNatRich <- data.frame(as.matrix(BetaFuncNat[["Brich"]]), sep= "tab")
BetaFuncNatRepl <- data.frame(as.matrix(BetaFuncNat[["Brepl"]]), sep= "tab")
BetaFuncNIndTotal <- data.frame(as.matrix(BetaFuncNInd[["Btotal"]]), sep= "tab")
BetaFuncNIndRich <- data.frame(as.matrix(BetaFuncNInd[["Brich"]]), sep= "tab")
BetaFuncNIndRepl <- as.data.frame(as.matrix(BetaFuncNInd[["Brepl"]]), sep= "tab")

# Separating the TAXONOMICAL beta values between the Control 250/Max and the other sampling areas from each trail ----

A01 <- BetaAllTotal[1,1]
AA1 <- BetaAllTotal[6,c(2:6)]
BB1 <- BetaAllTotal[11,c(7:11)]
CC1 <- BetaAllTotal[17,c(12:17)]
DD1 <- BetaAllTotal[21,c(18:21)]
all.tax.btotal <- c(A01,AA1, BB1,  CC1, DD1)

A02 <- BetaAllRich[1,1]
AA2 <- BetaAllRich[6,c(2:6)]
BB2 <- BetaAllRich[11,c(7:11)]
CC2 <- BetaAllRich[17,c(12:17)]
DD2<- BetaAllRich[21,c(18:21)]
all.tax.brich <- c(A02,AA2, BB2,  CC2, DD2)

A03 <- BetaAllRepl[1,1]
AA3 <- BetaAllRepl[6,c(2:6)]
BB3 <- BetaAllRepl[11,c(7:11)]
CC3 <- BetaAllRepl[17,c(12:17)]
DD3 <- BetaAllRepl[21,c(18:21)]
all.tax.brepl <- c(A03,AA3, BB3, CC3, DD3)

A01 <- BetaNatTotal[1,1]
AA1 <- BetaNatTotal[6,c(2:6)]
BB1 <- BetaNatTotal[11,c(7:11)]
CC1 <- BetaNatTotal[17,c(12:17)]
DD1 <- BetaNatTotal[21,c(18:21)]
nat.tax.btotal <- c(A01,AA1, BB1,  CC1, DD1)

A02 <- BetaNatRich[1,1]
AA2 <- BetaNatRich[6,c(2:6)]
BB2 <- BetaNatRich[11,c(7:11)]
CC2 <- BetaNatRich[17,c(12:17)]
DD2<- BetaNatRich[21,c(18:21)]
nat.tax.brich <- c(A02,AA2, BB2,  CC2, DD2)

A03 <- BetaNatRepl[1,1]
AA3 <- BetaNatRepl[6,c(2:6)]
BB3 <- BetaNatRepl[11,c(7:11)]
CC3 <- BetaNatRepl[17,c(12:17)]
DD3 <- BetaNatRepl[21,c(18:21)]
nat.tax.brepl <- c(A03,AA3, BB3, CC3, DD3)

A01 <- BetaNIndTotal[1,1]
AA1 <- BetaNIndTotal[6,c(2:6)]
BB1 <- BetaNIndTotal[11,c(7:11)]
CC1 <- BetaNIndTotal[17,c(12:17)]
DD1 <- BetaNIndTotal[21,c(18:21)]
nind.tax.btotal <- c(A01,AA1, BB1,  CC1, DD1)

A02 <- BetaNIndRich[1,1]
AA2 <- BetaNIndRich[6,c(2:6)]
BB2 <- BetaNIndRich[11,c(7:11)]
CC2 <- BetaNIndRich[17,c(12:17)]
DD2<- BetaNIndRich[21,c(18:21)]
nind.tax.brich <- c(A02,AA2, BB2,  CC2, DD2)

A03 <- BetaNIndRepl[1,1]
AA3 <- BetaNIndRepl[6,c(2:6)]
BB3 <- BetaNIndRepl[11,c(7:11)]
CC3 <- BetaNIndRepl[17,c(12:17)]
DD3 <- BetaNIndRepl[21,c(18:21)]
nind.tax.brepl <- c(A03,AA3, BB3, CC3, DD3)

###### Separating the FUNCTIONAL beta values between the Control 250/Max and the other sampling areas from each trail ---


A01 <- BetaFuncAllTotal[1,1]
AA1 <- BetaFuncAllTotal[6,c(2:6)]
BB1 <- BetaFuncAllTotal[11,c(7:11)]
CC1 <- BetaFuncAllTotal[17,c(12:17)]
DD1 <- BetaFuncAllTotal[21,c(18:21)]
all.func.btotal <- c(A01,AA1, BB1,  CC1, DD1)

A02 <- BetaFuncAllRich[1,1]
AA2 <- BetaFuncAllRich[6,c(2:6)]
BB2 <- BetaFuncAllRich[11,c(7:11)]
CC2 <- BetaFuncAllRich[17,c(12:17)]
DD2<- BetaFuncAllRich[21,c(18:21)]
all.func.brich <- c(A02,AA2, BB2,  CC2, DD2)

A03 <- BetaFuncAllRepl[1,1]
AA3 <- BetaFuncAllRepl[6,c(2:6)]
BB3 <- BetaFuncAllRepl[11,c(7:11)]
CC3 <- BetaFuncAllRepl[17,c(12:17)]
DD3 <- BetaFuncAllRepl[21,c(18:21)]
all.func.brepl <- c(A03,AA3, BB3, CC3, DD3)

A01 <- BetaFuncNatTotal[1,1]
AA1 <- BetaFuncNatTotal[6,c(2:6)]
BB1 <- BetaFuncNatTotal[11,c(7:11)]
CC1 <- BetaFuncNatTotal[17,c(12:17)]
DD1 <- BetaFuncNatTotal[21,c(18:21)]
nat.func.btotal <- c(A01,AA1, BB1,  CC1, DD1)

A02 <- BetaFuncNatRich[1,1]
AA2 <- BetaFuncNatRich[6,c(2:6)]
BB2 <- BetaFuncNatRich[11,c(7:11)]
CC2 <- BetaFuncNatRich[17,c(12:17)]
DD2<- BetaFuncNatRich[21,c(18:21)]
nat.func.brich <- c(A02,AA2, BB2,  CC2, DD2)

A03 <- BetaFuncNatRepl[1,1]
AA3 <- BetaFuncNatRepl[6,c(2:6)]
BB3 <- BetaFuncNatRepl[11,c(7:11)]
CC3 <- BetaFuncNatRepl[17,c(12:17)]
DD3 <- BetaFuncNatRepl[21,c(18:21)]
nat.func.brepl <- c(A03,AA3, BB3, CC3, DD3)

A01 <- BetaFuncNIndTotal[1,1]
AA1 <- BetaFuncNIndTotal[6,c(2:6)]
BB1 <- BetaFuncNIndTotal[11,c(7:11)]
CC1 <- BetaFuncNIndTotal[17,c(12:17)]
DD1 <- BetaFuncNIndTotal[21,c(18:21)]
nind.func.btotal <- c(A01,AA1, BB1,  CC1, DD1)

A02 <- BetaFuncNIndRich[1,1]
AA2 <- BetaFuncNIndRich[6,c(2:6)]
BB2 <- BetaFuncNIndRich[11,c(7:11)]
CC2 <- BetaFuncNIndRich[17,c(12:17)]
DD2<- BetaFuncNIndRich[21,c(18:21)]
nind.func.brich <- c(A02,AA2, BB2,  CC2, DD2)

A03 <- BetaFuncNIndRepl[1,1]
AA3 <- BetaFuncNIndRepl[6,c(2:6)]
BB3 <- BetaFuncNIndRepl[11,c(7:11)]
CC3 <- BetaFuncNIndRepl[17,c(12:17)]
DD3 <- BetaFuncNIndRepl[21,c(18:21)]
nind.func.brepl <- c(A03,AA3, BB3, CC3, DD3)


#### COMPILING ALL BETA INFORMATION INTO A TABLE, AND EXPORTING IT TO A FILE

Betas <- as.data.frame( t(rbind(all.tax.btotal, all.tax.brich, all.tax.brepl, nat.tax.btotal, nat.tax.brich, nat.tax.brepl, nind.tax.btotal, 
                                nind.tax.brich, nind.tax.brepl, all.func.btotal, all.func.brich, 
                                all.func.brepl, nat.func.btotal, nat.func.brich, nat.func.brepl, nind.func.btotal, nind.func.brich, nind.func.brepl)))

str(Betas)
Betas[Betas == 0] <- 0.0001
Betas[Betas == 1] <- 0.9999
str(Betas)


## Exporting TAXONOMICAL beta results to .csv (for the record - it won't be usedin analysis, as the RESULTS 
##      file compiles all the  useable results )



write.csv(BetaAllTotal,file = here("results", "Control250_Fin.csv"))
write.csv(BetaAllRich,file = here("results", "BetaAllRich.csv"))
write.csv(BetaAllRepl,file = here("results", "BetaAllRepl.csv"))
write.csv(BetaNatTotal,file = here("results", "BetaNatTotal.csv"))
write.csv(BetaNatRich,file = here("results", "BetaNatRich.csv"))
write.csv(BetaNatRepl,file = here("results", "BetaNatRepl.csv"))
write.csv(BetaNIndTotal,file = here("results", "BetaNIndTotal.csv"))
write.csv(BetaNIndRich,file = here("results", "BetaNIndRich.csv"))
write.csv(BetaNIndRepl,file = here("results", "BetaNIndRepl.csv"))

## Exporting FUNCTIONAL beta results to .csv (for the record - it won't be usedin analysis, as the RESULTS 
#       file compiles all the  useable results )

write.csv(BetaFuncAllTotal,file = here("results", "BetaFuncAllTotal.csv"))
write.csv(BetaFuncAllRich,file = here("results", "BetaFuncAllRich.csv"))
write.csv(BetaFuncAllRepl,file = here("results", "BetaFuncAllRepl.csv"))
write.csv(BetaFuncNatTotal,file = here("results", "BetaFuncNatTotal.csv"))
write.csv(BetaFuncNatRich,file = here("results", "BetaFuncNatRich.csv"))
write.csv(BetaFuncNatRepl,file = here("results", "BetaFuncNatRepl.csv"))
write.csv(BetaFuncNIndTotal,file = here("results", "BetaFuncNIndTotal.csv"))
write.csv(BetaFuncNIndRich,file = here("results", "BetaFuncNIndRich.csv"))
write.csv(BetaFuncNIndRepl,file = here("results", "BetaFuncNIndRepl.csv"))

###########################################################################
#####                            RESULTS                              #####
###########################################################################

str(Variables)
str(Alphas)
str(Betas)
Results <- cbind.data.frame(Variables, Alphas, Betas)
str(Results)


#MAKING THE RESULTS EXPORTABLE INTO CSV
Results <- apply(Results, 2 , as.character)



#NAMING THE TRAIL SEGMENTS
rownames(Results) <- rownames(Alphas)

#PASSING RESULTS TO FILE
write.csv(Results, file = here("results","RESULTS.csv"), row.names = TRUE)

Results2 <- read.csv2(here("results","RESULTS.CSV"), header=TRUE, row.names = 1,  stringsAsFactors = T,sep = ",", dec = ".")
Results2

Results3 <- read.csv(here("results","RESULTS.csv"), header = TRUE, row.names = 1, sep = ",", dec = ".")
ResultsWithoutControls <- Results3[-c(5,6,10,11,16,17),-c(1:2)]

###########################################################################
#####                            Generating models                              #####
###########################################################################

fam3 <- c(NA, NA, NA, NA,NA, NA,NA, "poisson", "poisson", "poisson", "Gamma", "Gamma", "Gamma","nbinom1", "nbinom1", "nbinom1","beta_family", "beta_family","beta_family","beta_family", "beta_family", 
          "beta_family", "beta_family", "beta_family", "beta_family", "beta_family", "beta_family", "beta_family", "beta_family", "beta_family", "beta_family", "beta_family", "beta_family", 
          "beta_family", "beta_family", "beta_family", "beta_family")


models.df = data.frame()
Models = list() 
for(i in 8:37){
  newTable = Results2[,c(i,1,5,6,7)]
  colnames(newTable)[1] = "y"
  Models[[i]] = dredge(glmmTMB(y ~ Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data= newTable, family = fam3[i]))
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

write.csv(resultspart1, file = here("results","dredge.models2.csv"), row.names = TRUE)

###########################################################################
#####             Retrieving model weights                            #####
###########################################################################
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

write.csv(aic.weights, file = here("results","aic.weights.csv"), row.names = TRUE)

#selecting delta <2

AICmodels <- filter(models.df, delta < 2)
AICmodels

write.csv(aic.weights, file = here("results","aic.weights.csv"), row.names = TRUE)


###########################################################################
#####                 Examining each model                            #####
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

###### Applying the formula to the models AIC<2


TAlphaAll.glmm.1 = glmmTMB(TAlphaAll ~    (1 | ForestID), data= Results2,family = "poisson") 
model.output(TAlphaAll.glmm.1)  
TAlphaNat.glmm.2 = glmmTMB(TAlphaNat ~    (1 | ForestID), data= Results2,family = "poisson") 
model.output(TAlphaNat.glmm.2)  
TAlphaNInd.glmm.3 = glmmTMB(TAlphaNInd ~    (1 | ForestID), data= Results2,family = "poisson") 
model.output(TAlphaNInd.glmm.3)  
TAlphaNInd.glmm.4 = glmmTMB(TAlphaNInd ~  Dist_trail_beginning_std + (1 | ForestID), data= Results2,family = "poisson") 
model.output(TAlphaNInd.glmm.4)  
TAlphaNInd.glmm.5 = glmmTMB(TAlphaNInd ~   Dist_trail_std +(1 | ForestID), data= Results2,family = "poisson") 
model.output(TAlphaNInd.glmm.5)  
TAlphaNInd.glmm.6 = glmmTMB(TAlphaNInd ~ Dist_edge_std +  (1 | ForestID), data= Results2,family = "poisson") 
model.output(TAlphaNInd.glmm.6)  
FAlphaAll.glmm.7 = glmmTMB(FAlphaAll ~    (1 | ForestID), data= Results2,family = "Gamma") 
model.output(FAlphaAll.glmm.7)  
FAlphaAll.glmm.8 = glmmTMB(FAlphaAll ~  Dist_trail_beginning_std + (1 | ForestID), data= Results2,family = "Gamma") 
model.output(FAlphaAll.glmm.8)  
FAlphaNat.glmm.9 = glmmTMB(FAlphaNat ~    (1 | ForestID), data= Results2,family = "Gamma") 
model.output(FAlphaNat.glmm.9)  
FAlphaNInd.glmm.10 = glmmTMB(FAlphaNInd ~    (1 | ForestID), data= Results2,family = "poisson") 
model.output(FAlphaNInd.glmm.10)  
abund.all.glmm.11 = glmmTMB(abund.all ~ Dist_edge_std +Dist_trail_beginning_std +Dist_trail_std +(1 | ForestID), data= Results2,family = "nbinom1") 
model.output(abund.all.glmm.11)  
abund.nat.glmm.12 = glmmTMB(abund.nat ~ Dist_edge_std +Dist_trail_beginning_std +Dist_trail_std +(1 | ForestID), data= Results2,family = "nbinom1") 
model.output(abund.nat.glmm.12)  
abund.nind.glmm.13 = glmmTMB(abund.nind ~ Dist_edge_std +  (1 | ForestID), data= Results2,family = "nbinom1") 
model.output(abund.nind.glmm.13)  
abund.nind.glmm.14 = glmmTMB(abund.nind ~   Dist_trail_std +(1 | ForestID), data= Results2,family = "nbinom1") 
model.output(abund.nind.glmm.14)  
prop.Talpha.glmm.15 = glmmTMB(prop.Talpha ~    (1 | ForestID), data= Results2,family = "beta_family") 
model.output(prop.Talpha.glmm.15)  
prop.Talpha.glmm.16 = glmmTMB(prop.Talpha ~   Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
model.output(prop.Talpha.glmm.16)  
prop.Falpha.glmm.17 = glmmTMB(prop.Falpha ~    (1 | ForestID), data= Results2,family = "beta_family") 
model.output(prop.Falpha.glmm.17)  
prop.Falpha.glmm.18 = glmmTMB(prop.Falpha ~  Dist_trail_beginning_std + (1 | ForestID), data= Results2,family = "beta_family") 
model.output(prop.Falpha.glmm.18)  
prop.abund.glmm.19 = glmmTMB(prop.abund ~    (1 | ForestID), data= Results2,family = "beta_family") 
model.output(prop.abund.glmm.19)  
all.tax.btotal.glmm.20 = glmmTMB(all.tax.btotal ~   Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
model.output(all.tax.btotal.glmm.20)  
all.tax.btotal.glmm.21 = glmmTMB(all.tax.btotal ~ Dist_edge_std + Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
model.output(all.tax.btotal.glmm.21)  
all.tax.brich.glmm.22 = glmmTMB(all.tax.brich ~   Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
model.output(all.tax.brich.glmm.22)  
all.tax.brich.glmm.23 = glmmTMB(all.tax.brich ~ Dist_edge_std + Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
model.output(all.tax.brich.glmm.23)  
all.tax.brepl.glmm.24 = glmmTMB(all.tax.brepl ~   Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
model.output(all.tax.brepl.glmm.24)  
all.tax.brepl.glmm.25 = glmmTMB(all.tax.brepl ~ Dist_edge_std + Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
model.output(all.tax.brepl.glmm.25)  
nat.tax.btotal.glmm.26 = glmmTMB(nat.tax.btotal ~   Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
model.output(nat.tax.btotal.glmm.26)  
nat.tax.btotal.glmm.27 = glmmTMB(nat.tax.btotal ~ Dist_edge_std + Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
model.output(nat.tax.btotal.glmm.27)  
nat.tax.brich.glmm.28 = glmmTMB(nat.tax.brich ~    (1 | ForestID), data= Results2,family = "beta_family") 
model.output(nat.tax.brich.glmm.28)  
nat.tax.brepl.glmm.29 = glmmTMB(nat.tax.brepl ~   Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
model.output(nat.tax.brepl.glmm.29)  
nat.tax.brepl.glmm.30 = glmmTMB(nat.tax.brepl ~ Dist_edge_std + Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
model.output(nat.tax.brepl.glmm.30)  
nind.tax.btotal.glmm.31 = glmmTMB(nind.tax.btotal ~ Dist_edge_std + Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
model.output(nind.tax.btotal.glmm.31)  
nind.tax.btotal.glmm.32 = glmmTMB(nind.tax.btotal ~    (1 | ForestID), data= Results2,family = "beta_family") 
model.output(nind.tax.btotal.glmm.32)  
nind.tax.btotal.glmm.33 = glmmTMB(nind.tax.btotal ~   Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
model.output(nind.tax.btotal.glmm.33)  
nind.tax.brich.glmm.34 = glmmTMB(nind.tax.brich ~ Dist_edge_std + Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
model.output(nind.tax.brich.glmm.34)  
nind.tax.brepl.glmm.35 = glmmTMB(nind.tax.brepl ~    (1 | ForestID), data= Results2,family = "beta_family") 
model.output(nind.tax.brepl.glmm.35)  
nind.tax.brepl.glmm.36 = glmmTMB(nind.tax.brepl ~ Dist_edge_std +  (1 | ForestID), data= Results2,family = "beta_family") 
model.output(nind.tax.brepl.glmm.36)  
all.func.btotal.glmm.37 = glmmTMB(all.func.btotal ~   Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
model.output(all.func.btotal.glmm.37)  
all.func.btotal.glmm.38 = glmmTMB(all.func.btotal ~ Dist_edge_std + Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
model.output(all.func.btotal.glmm.38)  
all.func.brich.glmm.39 = glmmTMB(all.func.brich ~   Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
model.output(all.func.brich.glmm.39)  
all.func.brich.glmm.40 = glmmTMB(all.func.brich ~ Dist_edge_std + Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
model.output(all.func.brich.glmm.40)  
all.func.brepl.glmm.41 = glmmTMB(all.func.brepl ~   Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
model.output(all.func.brepl.glmm.41)  
all.func.brepl.glmm.42 = glmmTMB(all.func.brepl ~ Dist_edge_std + Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
model.output(all.func.brepl.glmm.42)  
nat.func.btotal.glmm.43 = glmmTMB(nat.func.btotal ~   Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
model.output(nat.func.btotal.glmm.43)  
nat.func.btotal.glmm.44 = glmmTMB(nat.func.btotal ~ Dist_edge_std + Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
model.output(nat.func.btotal.glmm.44)  
nat.func.brich.glmm.45 = glmmTMB(nat.func.brich ~   Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
model.output(nat.func.brich.glmm.45)  
nat.func.brich.glmm.46 = glmmTMB(nat.func.brich ~ Dist_edge_std + Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
model.output(nat.func.brich.glmm.46)  
nat.func.brepl.glmm.47 = glmmTMB(nat.func.brepl ~   Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
model.output(nat.func.brepl.glmm.47)  
nind.func.btotal.glmm.48 = glmmTMB(nind.func.btotal ~ Dist_edge_std + Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
model.output(nind.func.btotal.glmm.48)  
nind.func.brich.glmm.49 = glmmTMB(nind.func.brich ~ Dist_edge_std + Dist_trail_std +(1 | ForestID), data= Results2,family = "beta_family") 
model.output(nind.func.brich.glmm.49)  
nind.func.brepl.glmm.50 = glmmTMB(nind.func.brepl ~    (1 | ForestID), data= Results2,family = "beta_family") 
model.output(nind.func.brepl.glmm.50)  
nind.func.brepl.glmm.51 = glmmTMB(nind.func.brepl ~ Dist_edge_std +  (1 | ForestID), data= Results2,family = "beta_family") 
model.output(nind.func.brepl.glmm.51)  



