# Header ----

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



here

# Load files --------------------

#Ficheiro com as variável de distâncias com edge, trilhos e etc - para o GLMM
#Variables <- read.csv2("GLMM_Variables.csv", header = TRUE, stringsAsFactors = F, na.strings = "-")
#Variables <- Variables
#str(Variables)

# Nova versão, copiada do menu Import Dataset. opção readr:
Variables <- read_delim("GLMM_Variables.csv", ";", escape_double = FALSE, 
                             col_types = cols(Dist_edge = col_number(), 
                            Dist_trail = col_number(), Dist_trail_beginning = col_number()), 
                             trim_ws = TRUE)
View(Variables)

# Ficheiro com as abundâncias por amostra para os controlos 250/Max.
Alpha_controls <- read.csv2("Control250_Fin.csv", header = TRUE)
str(Alpha_controls)


# Matrix with species (rows)x traits (cols) for all species in all plots


Traits  <- read.csv2(here("Traits_All_Fin.csv"), row.names=1, header=TRUE,  stringsAsFactors = T, dec = ".")
Traits <- Traits[,-c(1:3)] 
str(Traits)


# Matrix with species (rows)x traits (cols) for all natives species in all plots

Traits  <- read.csv2(here("Traits_Nat_Fin.csv"), row.names=1, header=TRUE,  stringsAsFactors = T, dec = ".")
Traits <- Traits[,-c(1:3)] 
str(Traits)

# Matrix with species (rows)x traits (cols) for all Non-Indigenous species in all plots

Traits  <- read.csv2(here("Traits_NInd_Fin.csv"), row.names=1, header=TRUE,  stringsAsFactors = T, dec = ".")
Traits <- Traits[,-c(1:3)] 
str(Traits)


# Ficheiro com as abundâncias por Área de amostragem, para todas as amostras
SAAll <- read.csv2(here("All_Fin.csv"))
sites <- SAAll[,1]
SAAll <- SAAll[,-1] 
rownames(SAAll) <- sites
colnames(SAAll) <- species
# HEllinger transformation
SAAll <- decostand(SAAll, "hellinger") 
SAAll

# Ficheiro com as abundâncias por área de amostragem, para as espécies indígenas

SANat <- read.csv2(here("Nat_Fin.csv"))
sites <- SANat[,1]
SANat <- SANat[,-1] 
rownames(SANat) <- sites
colnames(SANat) <- species_nat
# HEllinger transformation
SANat <- decostand(SANat, "hellinger") 
SANat

# Ficheiro com as abundâncias por área de amostragem, para as espécies não-indígenas
SANInd <- read.csv2(here("NInd_Fin.csv"))
sites <- SANInd[,1]
SANInd <- SANInd[,-1] 
rownames(SANInd) <- sites
colnames(SANInd) <- species_nind

# HEllinger transformation
SANInd <- decostand(SANInd, "hellinger") 
SANInd

# Vector with weight trail and treatment names for MDS
MDSfile <- read.csv2("MDS_vectors.csv")
trail <- MDSfile[,1]
trail
trail <- as.vector(trail)
treatment <- MDSfile[,2]
str(treatment)
treatment <- as.vector(treatment)
str(treatment)

# Vector with weight ratios for trails - used in cluster::daisy to compensate decomposition of one variable by
# multiple columns
Weightsfile <- read.csv2("Weight_Ratios_Traits.csv", header = TRUE, dec=".")
weights <- Weightsfile[,-c(1:4)]
weights <- as.vector(weights)

# Looking into the datasets -------------------------------

# Ver as estruturas dos dados (garantir os factors e os valores numericos)
str(Variables)
str(data.controls)
str(Traits)
str(TraitsNat)
str(TraitsNInd)
str(SAAll)
str(SANat)
str(SANInd)
str(weights)
str(treatment)
# CALCULATING COMPLETENESS FOR THE CONTROLS FOR BETA ANALYSYS ----

Alpha_controls <- read.csv2("Control250_Fin.csv", header = TRUE)

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
BetaAllTotal <- as.data.frame(as.matrix(BetaAll[["Btotal"]]), sep= "tab")
BetaAllRich <- as.data.frame(as.matrix(BetaAll[["Brich"]]), sep= "tab")
BetaAllRepl <- as.data.frame(as.matrix(BetaAll[["Brepl"]]), sep= "tab")
BetaNatTotal <- as.data.frame(as.matrix(BetaNat[["Btotal"]]), sep= "tab")
BetaNatRich <- as.data.frame(as.matrix(BetaNat[["Brich"]]), sep= "tab")
BetaNatRepl <- as.data.frame(as.matrix(BetaNat[["Brepl"]]), sep= "tab")
BetaNIndTotal <- as.data.frame(as.matrix(BetaNat[["Btotal"]]), sep= "tab")
BetaNIndRich <- as.data.frame(as.matrix(BetaNat[["Brich"]]), sep= "tab")
BetaNIndRepl <- as.data.frame(as.matrix(BetaNat[["Brepl"]]), sep= "tab")

## Separating Total, Richness and Replacement  FUNCTIONAL betas into different data frames
BetaFuncAllTotal <- as.data.frame(as.matrix(BetaFuncAll[["Btotal"]]), sep= "tab")
BetaFuncAllRich <- as.data.frame(as.matrix(BetaFuncAll[["Brich"]]), sep= "tab")
BetaFuncAllRepl <- as.data.frame(as.matrix(BetaFuncAll[["Brepl"]]), sep= "tab")
BetaFuncNatTotal <- as.data.frame(as.matrix(BetaFuncNat[["Btotal"]]), sep= "tab")
BetaFuncNatRich <- as.data.frame(as.matrix(BetaFuncNat[["Brich"]]), sep= "tab")
BetaFuncNatRepl <- as.data.frame(as.matrix(BetaFuncNat[["Brepl"]]), sep= "tab")
BetaFuncNIndTotal <- as.data.frame(as.matrix(BetaFuncNat[["Btotal"]]), sep= "tab")
BetaFuncNIndRich <- as.data.frame(as.matrix(BetaFuncNat[["Brich"]]), sep= "tab")
BetaFuncNIndRepl <- as.data.frame(as.matrix(BetaFuncNat[["Brepl"]]), sep= "tab")

# Separating the TAXONOMICAL beta values between the Control 250/Max and the other sampling areas from each trail ----

A01 <- BetaAllTotal[1,1]
AA1 <- BetaAllTotal[6,c(2:6)]
BB1 <- BetaAllTotal[12,c(7:12)]
CC1 <- BetaAllTotal[18,c(13:18)]
DD1 <- BetaAllTotal[22,c(19:22)]
BetaAllTotalVector <- c(A01,AA1, BB1,  CC1, DD1)

A02 <- BetaAllRich[1,1]
AA2 <- BetaAllRich[6,c(2:6)]
BB2 <- BetaAllRich[12,c(7:12)]
CC2 <- BetaAllRich[18,c(13:18)]
DD2<- BetaAllRich[22,c(19:22)]
BetaAllRichVector <- c(A02,AA2, BB2,  CC2, DD2)

A03 <- BetaAllRepl[1,1]
AA3 <- BetaAllRepl[6,c(2:6)]
BB3 <- BetaAllRepl[12,c(7:12)]
CC3 <- BetaAllRepl[18,c(13:18)]
DD3 <- BetaAllRepl[22,c(19:22)]
BetaAllReplVector <- c(A03,AA3, BB3, CC3, DD3)

A01 <- BetaNatTotal[1,1]
AA1 <- BetaNatTotal[6,c(2:6)]
BB1 <- BetaNatTotal[12,c(7:12)]
CC1 <- BetaNatTotal[18,c(13:18)]
DD1 <- BetaNatTotal[22,c(19:22)]
BetaNatTotalVector <- c(A01,AA1, BB1,  CC1, DD1)

A02 <- BetaNatRich[1,1]
AA2 <- BetaNatRich[6,c(2:6)]
BB2 <- BetaNatRich[12,c(7:12)]
CC2 <- BetaNatRich[18,c(13:18)]
DD2<- BetaNatRich[22,c(19:22)]
BetaNatRichVector <- c(A02,AA2, BB2,  CC2, DD2)

A03 <- BetaNatRepl[1,1]
AA3 <- BetaNatRepl[6,c(2:6)]
BB3 <- BetaNatRepl[12,c(7:12)]
CC3 <- BetaNatRepl[18,c(13:18)]
DD3 <- BetaNatRepl[22,c(19:22)]
BetaNatReplVector <- c(A03,AA3, BB3, CC3, DD3)

A01 <- BetaNIndTotal[1,1]
AA1 <- BetaNIndTotal[6,c(2:6)]
BB1 <- BetaNIndTotal[12,c(7:12)]
CC1 <- BetaNIndTotal[18,c(13:18)]
DD1 <- BetaNIndTotal[22,c(19:22)]
BetaNIndTotalVector <- c(A01,AA1, BB1,  CC1, DD1)

A02 <- BetaNIndRich[1,1]
AA2 <- BetaNIndRich[6,c(2:6)]
BB2 <- BetaNIndRich[12,c(7:12)]
CC2 <- BetaNIndRich[18,c(13:18)]
DD2<- BetaNIndRich[22,c(19:22)]
BetaNIndRichVector <- c(A02,AA2, BB2,  CC2, DD2)

A03 <- BetaNIndRepl[1,1]
AA3 <- BetaNIndRepl[6,c(2:6)]
BB3 <- BetaNIndRepl[12,c(7:12)]
CC3 <- BetaNIndRepl[18,c(13:18)]
DD3 <- BetaNIndRepl[22,c(19:22)]
BetaNIndReplVector <- c(A03,AA3, BB3, CC3, DD3)

## Separating the FUNCTIONAL beta values between the Control 250/Max and the other sampling areas from each trail


A01 <- BetaFuncAllTotal[1,1]
AA1 <- BetaFuncAllTotal[6,c(2:6)]
BB1 <- BetaFuncAllTotal[12,c(7:12)]
CC1 <- BetaFuncAllTotal[18,c(13:18)]
DD1 <- BetaFuncAllTotal[22,c(19:22)]
BetaFuncAllTotalVector <- c(A01,AA1, BB1,  CC1, DD1)

A02 <- BetaFuncAllRich[1,1]
AA2 <- BetaFuncAllRich[6,c(2:6)]
BB2 <- BetaFuncAllRich[12,c(7:12)]
CC2 <- BetaFuncAllRich[18,c(13:18)]
DD2<- BetaFuncAllRich[22,c(19:22)]
BetaFuncAllRichVector <- c(A02,AA2, BB2,  CC2, DD2)

A03 <- BetaFuncAllRepl[1,1]
AA3 <- BetaFuncAllRepl[6,c(2:6)]
BB3 <- BetaFuncAllRepl[12,c(7:12)]
CC3 <- BetaFuncAllRepl[18,c(13:18)]
DD3 <- BetaFuncAllRepl[22,c(19:22)]
BetaFuncAllReplVector <- c(A03,AA3, BB3, CC3, DD3)

A01 <- BetaFuncNatTotal[1,1]
AA1 <- BetaFuncNatTotal[6,c(2:6)]
BB1 <- BetaFuncNatTotal[12,c(7:12)]
CC1 <- BetaFuncNatTotal[18,c(13:18)]
DD1 <- BetaFuncNatTotal[22,c(19:22)]
BetaFuncNatTotalVector <- c(A01,AA1, BB1,  CC1, DD1)

A02 <- BetaFuncNatRich[1,1]
AA2 <- BetaFuncNatRich[6,c(2:6)]
BB2 <- BetaFuncNatRich[12,c(7:12)]
CC2 <- BetaFuncNatRich[18,c(13:18)]
DD2<- BetaFuncNatRich[22,c(19:22)]
BetaFuncNatRichVector <- c(A02,AA2, BB2,  CC2, DD2)

A03 <- BetaFuncNatRepl[1,1]
AA3 <- BetaFuncNatRepl[6,c(2:6)]
BB3 <- BetaFuncNatRepl[12,c(7:12)]
CC3 <- BetaFuncNatRepl[18,c(13:18)]
DD3 <- BetaFuncNatRepl[22,c(19:22)]
BetaFuncNatReplVector <- c(A03,AA3, BB3, CC3, DD3)

A01 <- BetaFuncNIndTotal[1,1]
AA1 <- BetaFuncNIndTotal[6,c(2:6)]
BB1 <- BetaFuncNIndTotal[12,c(7:12)]
CC1 <- BetaFuncNIndTotal[18,c(13:18)]
DD1 <- BetaFuncNIndTotal[22,c(19:22)]
BetaFuncNIndTotalVector <- c(A01,AA1, BB1,  CC1, DD1)

A02 <- BetaFuncNIndRich[1,1]
AA2 <- BetaFuncNIndRich[6,c(2:6)]
BB2 <- BetaFuncNIndRich[12,c(7:12)]
CC2 <- BetaFuncNIndRich[18,c(13:18)]
DD2<- BetaFuncNIndRich[22,c(19:22)]
BetaFuncNIndRichVector <- c(A02,AA2, BB2,  CC2, DD2)

A03 <- BetaFuncNIndRepl[1,1]
AA3 <- BetaFuncNIndRepl[6,c(2:6)]
BB3 <- BetaFuncNIndRepl[12,c(7:12)]
CC3 <- BetaFuncNIndRepl[18,c(13:18)]
DD3 <- BetaFuncNIndRepl[22,c(19:22)]
BetaFuncNIndReplVector <- c(A03,AA3, BB3, CC3, DD3)


#### COMPILING ALL BETA INFORMATION INTO A TABLE, AND EXPORTING IT TO A FILE

Betas <- as.data.frame( t(rbind(
  BetaAllTotalVector, BetaAllRichVector, BetaAllReplVector, 
  BetaFuncAllTotalVector, BetaFuncAllRichVector, BetaFuncAllReplVector,
  BetaNatTotalVector, BetaNatRichVector, BetaNatReplVector, 
  BetaFuncNatTotalVector, BetaFuncNatRichVector, BetaFuncNatReplVector,
  BetaNIndTotalVector, BetaNIndRichVector, BetaNIndReplVector,
  BetaFuncNIndTotalVector, BetaFuncNIndRichVector, BetaFuncNIndReplVector)))
## Exporting TAXONOMICAL beta results to .csv (for the record - it won't be usedin analysis, as the RESULTS 
##      file compiles all the  useable results )

XXXYYY <- Betas
XXXYYY
sdfsdf<- Alphas
sdfsdf

write.csv(BetaAllTotal,file = "C:/ArthropodsArticle/2.0 _Finland_Analysis/BetaResults/BetaAllTotal.csv")
write.csv(BetaAllRich,file = "C:/ArthropodsArticle/2.0 _Finland_Analysis/BetaResults/BetaAllRich.csv")
write.csv(BetaAllRepl,file = "C:/ArthropodsArticle/2.0 _Finland_Analysis/BetaResults/BetaAllRepl.csv")
write.csv(BetaNatTotal,file = "C:/ArthropodsArticle/2.0 _Finland_Analysis/BetaResults/BetaNatTotal.csv")
write.csv(BetaNatRich,file = "C:/ArthropodsArticle/2.0 _Finland_Analysis/BetaResults/BetaNatRich.csv")
write.csv(BetaNatRepl,file = "C:/ArthropodsArticle/2.0 _Finland_Analysis/BetaResults/BetaNatRepl.csv")
write.csv(BetaNIndTotal,file = "C:/ArthropodsArticle/2.0 _Finland_Analysis/BetaResults/BetaNIndTotal.csv")
write.csv(BetaNIndRich,file = "C:/ArthropodsArticle/2.0 _Finland_Analysis/BetaResults/BetaNIndRich.csv")
write.csv(BetaNIndRepl,file = "C:/ArthropodsArticle/2.0 _Finland_Analysis/BetaResults/BetaNIndRepl.csv")

## Exporting FUNCTIONAL beta results to .csv (for the record - it won't be usedin analysis, as the RESULTS 
#       file compiles all the  useable results )

write.csv(BetaFuncAllTotal,file = "C:/ArthropodsArticle/2.0 _Finland_Analysis/BetaFuncResults/BetaFuncAllTotal.csv")
write.csv(BetaFuncAllRich,file = "C:/ArthropodsArticle/2.0 _Finland_Analysis/BetaFuncResults/BetaFuncAllRich.csv")
write.csv(BetaFuncAllRepl,file = "C:/ArthropodsArticle/2.0 _Finland_Analysis/BetaFuncResults/BetaFuncAllRepl.csv")
write.csv(BetaFuncNatTotal,file = "C:/ArthropodsArticle/2.0 _Finland_Analysis/BetaFuncResults/BetaFuncNatTotal.csv")
write.csv(BetaFuncNatRich,file = "C:/ArthropodsArticle/2.0 _Finland_Analysis/BetaFuncResults/BetaFuncNatRich.csv")
write.csv(BetaFuncNatRepl,file = "C:/ArthropodsArticle/2.0 _Finland_Analysis/BetaFuncResults/BetaFuncNatRepl.csv")
write.csv(BetaFuncNIndTotal,file = "C:/ArthropodsArticle/2.0 _Finland_Analysis/BetaFuncResults/BetaFuncNIndTotal.csv")
write.csv(BetaFuncNIndRich,file = "C:/ArthropodsArticle/2.0 _Finland_Analysis/BetaFuncResults/BetaFuncNIndRich.csv")
write.csv(BetaFuncNIndRepl,file = "C:/ArthropodsArticle/2.0 _Finland_Analysis/BetaFuncResults/BetaFuncNIndRepl.csv")


##  abundance per site
raref.td <- alpha(SAAll, raref=1) #An integer specifying the number of individuals for rarefaction 
#(individual based). If raref < 1 no rarefaction is made. If raref = 1 rarefaction is made by the minimum 
#abundance among all sites. If raref > 1 rarefaction is made by the abundance indicated. If not specified, 
#default is 0.
raref.td <- alpha (SAAll, raref=1, runs=1000)
raref.fd <- alpha (SAAll, tree, raref=1, runs=1000)

par (mfrow=c(1,2), cex = 1,5, mar=c(11,6,2,2))

dev.copy(device = jpeg, filename = ' TD and FD Rarefaction Boxplots.jpeg', width = 1000, height = 500)
boxplot(t(raref.td), names=sites, main="TD rarefaction", ylab=expression(alpha), las = 2)
boxplot(t(raref.fd), names=sites, main="FD rarefaction", ylab=expression(alpha), las = 2)
dev.off()



rarefNat.td <- alpha (SANat, raref=1, runs=1000)
rarefNat.fd <- alpha (SANat, treeNat, raref=1, runs=1000)

par (mfrow=c(1,2), cex = 1,5, mar=c(11,6,2,2))

dev.copy(device = jpeg, filename = ' TD and FD Nat Rarefaction Boxplots.jpeg', width = 1000, height = 500)
boxplot(t(rarefNat.td), names=sites, main="TD rarefaction", ylab=expression(alpha), las = 2)
boxplot(t(rarefNat.fd), names=sites, main="FD rarefaction", ylab=expression(alpha), las = 2)
dev.off()

rarefNInd.td <- alpha (SANInd, raref=1, runs=1000)
rarefNInd.fd <- alpha (SANInd, treeNInd, raref=1, runs=1000)

par (mfrow=c(1,2), cex = 1,5, mar=c(11,6,2,2))

dev.copy(device = jpeg , filename = ' TD and FD NInd Rarefaction Boxplots.jpeg', width = 1000, height = 500)
boxplot(t(rarefNInd.td), names=sites, main="TD rarefaction", ylab=expression(alpha), las = 2)
boxplot(t(rarefNInd.fd), names=sites, main="FD rarefaction", ylab=expression(alpha), las = 2)
dev.off()

## To assess functional diversity of each site
# matrix with all sites's abundances for each species
#make ahh collumn headers be the same


BAT::dispersion(SAAll, tree)


beta.multi(SAAll, tree)
# Average   Variance
# Btotal 0.4554385 0.05825376
# Brepl  0.2412576 0.03085853
# Brich  0.2141809 0.02739523



###########################################################################
#####                            RESULTS                              #####
###########################################################################


Results <- cbind.data.frame(Variables, Alphas, Betas)
str(Results)

Results4 <- as.data.frame(matrix(Results))
str(Results4)

Results2 <- as.list(cbind.data.frame(Variables, Alphas, Betas))
str(Results)
Results8 <- read.csv("RESULTS2.csv", header = TRUE)




#MAKING THE RESULTS EXPORTABLE INTO CSV
Results <- apply(Results, 2 , as.character)

#NAMING THE TRAIL SEGMENTS
rownames(Results) <- rownames(Alphas)

#PASSING RESULTS TO FILE
write.csv(Results, file = "C:/ArthropodsArticle/2.0 _Finland_Analysis/RESULTS.csv")


Results3 <- read.csv("RESULTS.csv", header = TRUE)
ResultsWithoutControls <- Results3[-c(5,6,11,12,17,18),-c(1:2)]

###########################################################################
#####                       DATA EXPLORATION                          #####
###########################################################################

str(Results)

# 1. Check outliers
str(Results)
ResultsNumeric <- Results[,3:29 ]
VariablesNumeric <- Variables[,3:5]


the_variables <- c("Dist_trail", "Dist_edge", "Dist_trail_beginning", "TAlphaAll", "TAlphaNat", "TAlphaNInd","FAlphaAll", "FAlphaNat", "FAlphaNInd",
                   "BetaAllTotalVector", "BetaAllRichVector", "BetaAllReplVector", 
                   "BetaFuncAllTotalVector", "BetaFuncAllRichVector", "BetaFuncAllReplVector",
                   "BetaNatTotalVector", "BetaNatRichVector", "BetaNatReplVector", 
                   "BetaFuncNatTotalVector", "BetaFuncNatRichVector", "BetaFuncNatReplVector",
                   "BetaNIndTotalVector", "BetaNIndRichVector", "BetaNIndReplVector",
                   "BetaFuncNIndTotalVector", "BetaFuncNIndRichVector", "BetaFuncNIndReplVector")
the_variables2 <- c("TAlphaAll", "TAlphaNat", "TAlphaNInd","FAlphaAll")

the_variables
the_variables2

#Não encontrei a função mydotplot
Mydotplot(Results3[,the_variables])
dev.off()
# Zero trouble
# Os meus dados para as não indígenas têm bastante zeros, mas não estou a perceber como posso
#avaliar se é um problema, e o que enho de fazer a seguir
# é aplicável para count data - neste caso, iríamos ao ficheiro de abundâncias que alimenta os resultados

#tentativa de aplicar a uma das variáveis dependentes
sum(Results$TAlphaAll == 0)  #Number of zeros
100 * sum(Results$TAlphaAll == 0) / nrow(ResultsWithoutControls)  #% of zeros

#COLINEARIDADE
the_independent_variables <- c("Dist_trail", "Dist_edge", "Dist_trail_beginning")

pairs(Results3[, the_independent_variables], 
      lower.panel = panel.cor)
## Não funciona:  "Error in pairs.default(Results[, the_variables], lower.panel = panel.cor) : 
##non-numeric argument to 'pairs"

# RELAÇÕES ENTRE VARIÁVEIS x E y
#

##Tive de retirar os tratamentos (linhas que tinham NA, para funcionar)
## ´Tens alguma sugestão para fazer multiplot, de forma apoder ver todas as variáveis dependentes?

the_variables3 <-c("Dist_trail_beginning", "Dist_edge")
Myxyplot(Results3, the_variables3, "FAlphaNat",
         MyYlab = "Species diversity")
par(mfrow = c(1,2))

boxplot(TAlphaAll ~ Dist_trail_beginning, data = Results3)
boxplot(TAlphaAll ~ Dist_edge, data = Results3)

###########################################################################
#####                           MDS                                  #####
###########################################################################
par(mfrow=c(1,1))

#£metaMDS(Results3, distance = "bray", k = 2, try = 20, trymax = 20)

#ordiellipse(mds, groups  = as.factor(treatment), label = T)
mds1 <- metaMDS(BetaAll$Btotal, k=2, trymax = 100)
plot(mds1)
ordiellipse(mds1, groups  = as.factor(trail), label = T)
ordiellipse(mds1, groups  = as.factor(treatment), label = T)

mds2 <- metaMDS(BetaAll$Brepl, k=3, trymax = 100)
plot(mds2)
ordiellipse(mds2, groups  = as.factor(trail), label = T)
ordiellipse(mds2, groups  = as.factor(treatment), label = T)

mds3 <- metaMDS(BetaNat$Btotal, k=3, trymax = 100)
plot(mds3)
ordiellipse(mds3, groups  = as.factor(trail), label = T)
ordiellipse(mds3, groups  = as.factor(treatment), label = T)

mds4 <- metaMDS(BetaNat$Brepl, k=3, trymax = 100)
plot(mds4)
ordiellipse(mds4, groups  = as.factor(trail), label = T)
ordiellipse(mds4, groups  = as.factor(treatment), label = T)

mds5 <- metaMDS(BetaFuncAll$Btotal, k=3, trymax = 100)
plot(mds5)
ordiellipse(mds5, groups  = as.factor(trail), label = T)
ordiellipse(mds5, groups  = as.factor(treatment), label = T)

mds6 <- metaMDS(BetaFuncAll$Brepl, k=3, trymax = 100)
plot(mds6)
ordiellipse(mds6, groups  = as.factor(trail), label = T)
ordiellipse(mds6, groups  = as.factor(treatment), label = T)

mds7 <- metaMDS(BetaFuncNat$Btotal, k=3, trymax = 100)
plot(mds7)
ordiellipse(mds7, groups  = as.factor(trail), label = T)
ordiellipse(mds7, groups  = as.factor(treatment), label = T)

mds8 <- metaMDS(BetaFuncNat$Brepl, k=3, trymax = 100)
plot(mds8)
ordiellipse(mds8, groups  = as.factor(trail), label = T)
ordiellipse(mds8, groups  = as.factor(treatment), label = T)



###########################################################################
#####                           GLMM                                  #####
###########################################################################

Results8 <- read.csv("RESULTS2.csv", header = TRUE)

###
### ALPHAS GLMM###
###

##### All individuals################


#TAlphaAll

glmm_TAlphaAll_AB<-glmer(TAlphaAll ~ Dist_edge + Dist_trail_beginning + (1 | ForestID), data = Results8 , family = poisson)
summary(glmm_TAlphaAll_AB)

glmmT_AlphaAll_A<-glmer(TAlphaAll ~ Dist_edge + (1 | ForestID), data = Results8, family = poisson)
summary(glmmT_AlphaAll_A)

glmm_TAlphaAll_B<-glmer(TAlphaAll ~  Dist_trail_beginning + (1 | ForestID), data = Results8 , family = poisson)
summary(glmm_TAlphaAll_B)

glmm_TAlphaAll_C<-glmer(TAlphaAll ~  treatment + Dist_edge + (1 | ForestID), data = Results8 , family = poisson)
summary(glmm_TAlphaAll_C)

#FAlphaAll

glmm_FAlphaAll_AB <- glmer(FAlphaAll ~ Dist_edge + Dist_trail_beginning + (1 | ForestID), data = Results8 , family = gaussian)
summary(glmm_FAlphaAll_AB)

glmmT_AlphaAll_A<-glmer(FAlphaAll ~ Dist_edge + (1 | ForestID), data = Results8, family = gaussian)
summary(glmmT_AlphaAll_A)

glmm_FAlphaAll_B<-glmer(FAlphaAll ~  Dist_trail_beginning + (1 | ForestID), data = Results8 , family = gaussian)
summary(glmm_FAlphaAll_B)

glmm_FAlphaAll_C<-glmer(FAlphaAll ~  treatment+ (1 | ForestID), data = Results8 , family = gaussian)
summary(glmm_FAlphaAll_C)


##### Natives ####

#TAlphaAll

glmm_TAlphaNat_AB<-glmer(TAlphaNat ~ Dist_edge + Dist_trail_beginning + (1 | ForestID), data = Results8 , family = poisson)
summary(glmm_TAlphaNat_AB)

glmmT_AlphaNat_A<-glmer(TAlphaNat ~ Dist_edge + (1 | ForestID), data = Results8, family = poisson)
summary(glmmT_AlphaNat_A)

glmm_TAlphaNat_B<-glmer(TAlphaNat ~  Dist_trail_beginning + (1 | ForestID), data = Results8 , family = poisson)
summary(glmm_TAlphaNat_B)

glmm_TAlphaNat_c<-glmer(TAlphaNat ~  Dist_edge + treatment + (1 | ForestID), data = Results8 , family = poisson)
summary(glmm_TAlphaNat_c)

#FAlphaNat

glmm_FAlphaNat_AB <- glmer(FAlphaNat ~ Dist_edge + Dist_trail_beginning + (1 | ForestID), data = Results8 , family = gaussian)
summary(glmm_FAlphaNat_AB)

glmmT_AlphaNat_A<-glmer(FAlphaNat ~ Dist_edge + (1 | ForestID), data = Results8, family = gaussian)
summary(glmmT_AlphaNat_A)

glmm_FAlphaNat_B<-glmer(FAlphaNat ~  Dist_trail_beginning + (1 | ForestID), data = Results8 , family = gaussian)
summary(glmm_FAlphaNat_B)

glmm_FAlphaNat_C<-glmer(FAlphaNat ~  treatment + (1 | ForestID), data = Results8 , family = gaussian)
summary(glmm_FAlphaNat_C)

##### Non-Natives ################

#TAlphaAll

glmm_TAlphaNInd_AB<-glmer(TAlphaNInd ~ Dist_edge + Dist_trail_beginning + (1 | ForestID), data = Results8 , family = poisson)
summary(glmm_TAlphaNInd_AB)

glmmT_AlphaNInd_A<-glmer(TAlphaNInd ~ Dist_edge + (1 | ForestID), data = Results8, family = poisson)
summary(glmmT_AlphaNInd_A)



glmm_TAlphaNInd_B<-glmer(TAlphaNInd ~  Dist_trail_beginning + (1 | ForestID), data = Results8 , family = poisson)
summary(glmm_TAlphaNInd_B)


# O TC250 dá interação significativa
glmm_TAlphaNInd_C<-glmer(TAlphaNInd ~  treatment + (1 | ForestID), data = Results8 , family = poisson)
summary(glmm_TAlphaNInd_C)


#FAlphaNInd

glmm_FAlphaNInd_AB <- glmer(FAlphaNInd ~ Dist_edge + Dist_trail_beginning + (1 | ForestID), data = Results8 , family = gaussian)
summary(glmm_FAlphaNInd_AB)

glmmT_AlphaNInd_A<-glmer(FAlphaNInd ~ Dist_edge + (1 | ForestID), data = Results8, family = gaussian)
summary(glmmT_AlphaNInd_)A

glmm_FAlphaNInd_B<-glmer(FAlphaNInd ~  Dist_trail_beginning + (1 | ForestID), data = Results8 , family = gaussian)
summary(glmm_FAlphaNInd_B)

###
### BETAS GLMM###
###

##### All individuals################


#TBetaAll

glmm_BetaAllTotalVector_AB<-glmer(BetaAllTotalVector ~  treatment + (1 | ForestID), data = Results8 , family = binomial)
summary(glmm_BetaAllTotalVector_AB)

glmmT_BetaAllTotalVector_A<-glmer(BetaAllTotalVector ~ Dist_edge + (1 | ForestID), data = Results8, family = binomial)
summary(glmmT_BetaAllTotalVector_A)

glmm_BetaAllTotalVector_B<-glmer(BetaAllTotalVector ~  Dist_trail_beginning + (1 | ForestID), data = Results8 , family = binomial)
summary(glmm_BetaAllTotalVector_B)

glmm_BetaAllTotalVector_C<-glmer(BetaAllTotalVector ~  treatment + (1 | ForestID), data = Results8 , family = binomial)
summary(glmm_BetaAllTotalVector_C)

#FBetaAll

glmm_FBetaAll_AB <- glmer(BetaFuncAllTotalVector~ Dist_edge + treatment + (1 | ForestID), data = Results8 , family = binomial)
summary(glmm_FBetaAll_AB)

glmm_FBetaAll_A<-glmer(BetaFuncAllTotalVector ~ treatment + (1 | ForestID), data = Results8, family = binomial)
summary(glmm_FBetaAll_A)

glmm_FBetaAll_B <- glmer(BetaFuncAllTotalVector ~ Dist_edge + (1 | ForestID), data = Results8 , family = binomial)
summary(glmm_FBetaAll_B)

##### natives################


#TBetaAll

glmm_BetaNatTotalVector_AB<-glmer(BetaNatTotalVector ~  treatment + (1 | ForestID), data = Results8 , family = binomial)
summary(glmm_BetaNatTotalVector_AB)

glmmT_BetaNatTotalVector_A<-glmer(BetaNatTotalVector ~ Dist_edge + (1 | ForestID), data = Results8, family = binomial)
summary(glmmT_BetaNatTotalVector_A)

glmm_BetaNatTotalVector_B<-glmer(BetaNatTotalVector ~  Dist_trail_beginning + (1 | ForestID), data = Results8 , family = binomial)
summary(glmm_BetaNatTotalVector_B)

glmm_BetaNatTotalVector_C<-glmer(BetaNatTotalVector ~  treatment + (1 | ForestID), data = Results8 , family = binomial)
summary(glmm_BetaNatTotalVector_C)

#FAlphaNat

glmm_FBetaNat_AB <- glmer(BetaFuncNatTotalVector ~ Dist_edge + treatment + (1 | ForestID), data = Results8 , family = binomial)
summary(glmm_FBetaNat_AB)

glmmT_BetaNat_A<-glmer(BetaFuncNatTotalVector ~ treatment + (1 | ForestID), data = Results8, family = binomial)
summary(glmmT_AlphaNat_A)

glmm_FBetaNat_B <- glmer(BetaFuncNatTotalVector ~ Dist_edge + (1 | ForestID), data = Results8 , family = binomial)
summary(glmm_FBetaNat_B)







