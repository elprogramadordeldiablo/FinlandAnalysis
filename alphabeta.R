
  
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

# Ficheiro com as abundâncias por área de amostragem, para as espécies indígenas

SANat <- read.csv2(here("data","Nat_Fin.csv"),row.names=1, header = TRUE)

colnames(SANat) <- speciesNat 
# HEllinger transformation
#SANat <- decostand(SANat, "hellinger") 
#SANat


# Ficheiro com as abundâncias por área de amostragem, para as espécies não-indígenas
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

# Quero usar o fichero dos contrlos para calcular os estimadores, e daí obter completeness.
sitescontr <- Alpha_controls[,1]
species_contr <- Alpha_controls[,1]
rownames(Alpha_controls) <-sitescontr
colnames(Alpha_controls) <- species_contr
Alpha_controls <- Alpha_controls[,-1]
alpha.estimate(Alpha_controls)



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

# sdfjskdfjhskdjf ----
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

# Calculating Betas for All species, Natives and Non-Indigenous ----
# Beta results come as a list (class(BetaAll) = list), so we transform them into data frame to export them into csv

# Taxonomical Beta
BetaAll <- beta(SAAll)
BetaNat <- beta(SANat)
BetaNat <- beta(SANInd) 

#Functional Beta
BetaFuncAll <- beta(SAAll, tree, abund =  T)
BetaFuncNat <- beta(SANat, treeNat, abund = T) 
BetaFuncNInd <- beta(SANInd, treeNInd, abund= T)


## Separating Total, Richness and Replacement  TAXONOMICAL betas into different data frames
BetaAllTotal <- data.frame(as.matrix(BetaAll[["Btotal"]]), row.names= sites, sep= "tab")
colnames(BetaAllTotal) <- sites
BetaAllRich <- data.frame(as.matrix(BetaAll[["Brich"]]), sep= "tab")
BetaAllRepl <- data.frame(as.matrix(BetaAll[["Brepl"]]), sep= "tab")
BetaNatTotal <- data.frame(as.matrix(BetaNat[["Btotal"]]), sep= "tab")
BetaNatRich <- data.frame(as.matrix(BetaNat[["Brich"]]), sep= "tab")
BetaNatRepl <- data.frame(as.matrix(BetaNat[["Brepl"]]), sep= "tab")
BetaNIndTotal <- data.frame(as.matrix(BetaNat[["Btotal"]]), sep= "tab")
BetaNIndRich <- data.frame(as.matrix(BetaNat[["Brich"]]), sep= "tab")
BetaNIndRepl <- data.frame(as.matrix(BetaNat[["Brepl"]]), sep= "tab")

## Separating Total, Richness and Replacement  FUNCTIONAL betas into different data frames
BetaFuncAllTotal <- data.frame(as.matrix(BetaFuncAll[["Btotal"]]), sep= "tab")
BetaFuncAllRich <- data.frame(as.matrix(BetaFuncAll[["Brich"]]), sep= "tab")
BetaFuncAllRepl <- data.frame(as.matrix(BetaFuncAll[["Brepl"]]), sep= "tab")
BetaFuncNatTotal <- data.frame(as.matrix(BetaFuncNat[["Btotal"]]), sep= "tab")
BetaFuncNatRich <- data.frame(as.matrix(BetaFuncNat[["Brich"]]), sep= "tab")
BetaFuncNatRepl <- data.frame(as.matrix(BetaFuncNat[["Brepl"]]), sep= "tab")
BetaFuncNIndTotal <- data.frame(as.matrix(BetaFuncNat[["Btotal"]]), sep= "tab")
BetaFuncNIndRich <- data.frame(as.matrix(BetaFuncNat[["Brich"]]), sep= "tab")
BetaFuncNIndRepl <- as.data.frame(as.matrix(BetaFuncNat[["Brepl"]]), sep= "tab")

# Separating the TAXONOMICAL beta values between the Control 250/Max and the other sampling areas from each trail ----

A01 <- BetaAllTotal[1,1]
AA1 <- BetaAllTotal[6,c(2:6)]
BB1 <- BetaAllTotal[11,c(7:11)]
CC1 <- BetaAllTotal[17,c(12:17)]
DD1 <- BetaAllTotal[21,c(18:21)]
BetaAllTotalVector <- c(A01,AA1, BB1,  CC1, DD1)

A02 <- BetaAllRich[1,1]
AA2 <- BetaAllRich[6,c(2:6)]
BB2 <- BetaAllRich[11,c(7:11)]
CC2 <- BetaAllRich[17,c(12:17)]
DD2<- BetaAllRich[21,c(18:21)]
BetaAllRichVector <- c(A02,AA2, BB2,  CC2, DD2)

A03 <- BetaAllRepl[1,1]
AA3 <- BetaAllRepl[6,c(2:6)]
BB3 <- BetaAllRepl[11,c(7:11)]
CC3 <- BetaAllRepl[17,c(12:17)]
DD3 <- BetaAllRepl[21,c(18:21)]
BetaAllReplVector <- c(A03,AA3, BB3, CC3, DD3)

A01 <- BetaNatTotal[1,1]
AA1 <- BetaNatTotal[6,c(2:6)]
BB1 <- BetaNatTotal[11,c(7:11)]
CC1 <- BetaNatTotal[17,c(12:17)]
DD1 <- BetaNatTotal[21,c(18:21)]
BetaNatTotalVector <- c(A01,AA1, BB1,  CC1, DD1)

A02 <- BetaNatRich[1,1]
AA2 <- BetaNatRich[6,c(2:6)]
BB2 <- BetaNatRich[11,c(7:11)]
CC2 <- BetaNatRich[17,c(12:17)]
DD2<- BetaNatRich[21,c(18:21)]
BetaNatRichVector <- c(A02,AA2, BB2,  CC2, DD2)

A03 <- BetaNatRepl[1,1]
AA3 <- BetaNatRepl[6,c(2:6)]
BB3 <- BetaNatRepl[11,c(7:11)]
CC3 <- BetaNatRepl[17,c(12:17)]
DD3 <- BetaNatRepl[21,c(18:21)]
BetaNatReplVector <- c(A03,AA3, BB3, CC3, DD3)

A01 <- BetaNIndTotal[1,1]
AA1 <- BetaNIndTotal[6,c(2:6)]
BB1 <- BetaNIndTotal[11,c(7:11)]
CC1 <- BetaNIndTotal[17,c(12:17)]
DD1 <- BetaNIndTotal[21,c(18:21)]
BetaNIndTotalVector <- c(A01,AA1, BB1,  CC1, DD1)

A02 <- BetaNIndRich[1,1]
AA2 <- BetaNIndRich[6,c(2:6)]
BB2 <- BetaNIndRich[11,c(7:11)]
CC2 <- BetaNIndRich[17,c(12:17)]
DD2<- BetaNIndRich[21,c(18:21)]
BetaNIndRichVector <- c(A02,AA2, BB2,  CC2, DD2)

A03 <- BetaNIndRepl[1,1]
AA3 <- BetaNIndRepl[6,c(2:6)]
BB3 <- BetaNIndRepl[11,c(7:11)]
CC3 <- BetaNIndRepl[17,c(12:17)]
DD3 <- BetaNIndRepl[21,c(18:21)]
BetaNIndReplVector <- c(A03,AA3, BB3, CC3, DD3)

## Separating the FUNCTIONAL beta values between the Control 250/Max and the other sampling areas from each trail


A01 <- BetaFuncAllTotal[1,1]
AA1 <- BetaFuncAllTotal[6,c(2:6)]
BB1 <- BetaFuncAllTotal[11,c(7:11)]
CC1 <- BetaFuncAllTotal[17,c(12:17)]
DD1 <- BetaFuncAllTotal[21,c(18:21)]
BetaFuncAllTotalVector <- c(A01,AA1, BB1,  CC1, DD1)

A02 <- BetaFuncAllRich[1,1]
AA2 <- BetaFuncAllRich[6,c(2:6)]
BB2 <- BetaFuncAllRich[11,c(7:11)]
CC2 <- BetaFuncAllRich[17,c(12:17)]
DD2<- BetaFuncAllRich[21,c(18:21)]
BetaFuncAllRichVector <- c(A02,AA2, BB2,  CC2, DD2)

A03 <- BetaFuncAllRepl[1,1]
AA3 <- BetaFuncAllRepl[6,c(2:6)]
BB3 <- BetaFuncAllRepl[11,c(7:11)]
CC3 <- BetaFuncAllRepl[17,c(12:17)]
DD3 <- BetaFuncAllRepl[21,c(18:21)]
BetaFuncAllReplVector <- c(A03,AA3, BB3, CC3, DD3)

A01 <- BetaFuncNatTotal[1,1]
AA1 <- BetaFuncNatTotal[6,c(2:6)]
BB1 <- BetaFuncNatTotal[11,c(7:11)]
CC1 <- BetaFuncNatTotal[17,c(12:17)]
DD1 <- BetaFuncNatTotal[21,c(18:21)]
BetaFuncNatTotalVector <- c(A01,AA1, BB1,  CC1, DD1)

A02 <- BetaFuncNatRich[1,1]
AA2 <- BetaFuncNatRich[6,c(2:6)]
BB2 <- BetaFuncNatRich[11,c(7:11)]
CC2 <- BetaFuncNatRich[17,c(12:17)]
DD2<- BetaFuncNatRich[21,c(18:21)]
BetaFuncNatRichVector <- c(A02,AA2, BB2,  CC2, DD2)

A03 <- BetaFuncNatRepl[1,1]
AA3 <- BetaFuncNatRepl[6,c(2:6)]
BB3 <- BetaFuncNatRepl[11,c(7:11)]
CC3 <- BetaFuncNatRepl[17,c(12:17)]
DD3 <- BetaFuncNatRepl[21,c(18:21)]
BetaFuncNatReplVector <- c(A03,AA3, BB3, CC3, DD3)

A01 <- BetaFuncNIndTotal[1,1]
AA1 <- BetaFuncNIndTotal[6,c(2:6)]
BB1 <- BetaFuncNIndTotal[11,c(7:11)]
CC1 <- BetaFuncNIndTotal[17,c(12:17)]
DD1 <- BetaFuncNIndTotal[21,c(18:21)]
BetaFuncNIndTotalVector <- c(A01,AA1, BB1,  CC1, DD1)

A02 <- BetaFuncNIndRich[1,1]
AA2 <- BetaFuncNIndRich[6,c(2:6)]
BB2 <- BetaFuncNIndRich[11,c(7:11)]
CC2 <- BetaFuncNIndRich[17,c(12:17)]
DD2<- BetaFuncNIndRich[21,c(18:21)]
BetaFuncNIndRichVector <- c(A02,AA2, BB2,  CC2, DD2)

A03 <- BetaFuncNIndRepl[1,1]
AA3 <- BetaFuncNIndRepl[6,c(2:6)]
BB3 <- BetaFuncNIndRepl[11,c(7:11)]
CC3 <- BetaFuncNIndRepl[17,c(12:17)]
DD3 <- BetaFuncNIndRepl[21,c(18:21)]
BetaFuncNIndReplVector <- c(A03,AA3, BB3, CC3, DD3)


#### COMPILING ALL BETA INFORMATION INTO A TABLE, AND EXPORTING IT TO A FILE

Betas <- as.data.frame( t(rbind(
  BetaAllTotalVector, BetaAllRichVector, BetaAllReplVector, 
  BetaFuncAllTotalVector, BetaFuncAllRichVector, BetaFuncAllReplVector,
  BetaNatTotalVector, BetaNatRichVector, BetaNatReplVector, 
  BetaFuncNatTotalVector, BetaFuncNatRichVector, BetaFuncNatReplVector,
  BetaNIndTotalVector, BetaNIndRichVector, BetaNIndReplVector,
  BetaFuncNIndTotalVector, BetaFuncNIndRichVector, BetaFuncNIndReplVector)))

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

