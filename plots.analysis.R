

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

plot.all <- read.csv2(here("data.veg","plots.alpha.csv"), row.names=1, header=TRUE,  stringsAsFactors = T, dec = ".", sep = ",")
plot.inv <- read.csv2(here("data.veg","plots.alphainv.csv"), row.names=1, header=TRUE,  stringsAsFactors = T, dec = ".", sep = ",")

plot1vars <- read.csv2(here("data.veg","plot1vars.csv"), row.names=1, header=TRUE,  stringsAsFactors = T, dec = ".", sep = ",")
str(plot1vars)
plot1vars$Dist_trail <- as.numeric(plot1vars$Dist_trail)
plot1vars$Dist_edge <- as.numeric(plot1vars$Dist_edge)
plot1vars$Dist_trail_beginning <- as.numeric(plot1vars$Dist_trail_beginning)

plot1vars$Dist_trail_std <- scale(plot1vars$Dist_trail, center = F)
plot1vars$Dist_edge_std <- scale(plot1vars$Dist_edge, center = F)
plot1vars$Dist_trail_beginning_std <- scale(plot1vars$Dist_trail_beginning, center = F)

write.csv(plot1vars, file = here("data.veg","plot.vars.csv"), row.names = TRUE)

plot2vars <- read.csv2(here("data.veg","plot2vars.csv"), row.names=1, header=TRUE,  stringsAsFactors = T, dec = ".", sep = ",")
plot2vars$Dist_trail <- as.numeric(plot2vars$Dist_trail)
plot2vars$Dist_edge <- as.numeric(plot2vars$Dist_edge)
plot2vars$Dist_trail_beginning <- as.numeric(plot2vars$Dist_trail_beginning)

plot2vars$Dist_trail_std <- scale(plot2vars$Dist_trail, center = F)
plot2vars$Dist_edge_std <- scale(plot2vars$Dist_edge, center = F)
plot2vars$Dist_trail_beginning_std <- scale(plot2vars$Dist_trail_beginning, center = F)

plot3vars <- read.csv2(here("data.veg","plot3vars.csv"), row.names=1, header=TRUE,  stringsAsFactors = T, dec = ".", sep = ",")
str(plot3vars)
plot3vars$Dist_trail <- as.numeric(plot3vars$Dist_trail)
plot3vars$Dist_edge <- as.numeric(plot3vars$Dist_edge)
plot3vars$Dist_trail_beginning <- as.numeric(plot3vars$Dist_trail_beginning)

plot3vars$Dist_trail_std <- scale(plot3vars$Dist_trail, center = F)
plot3vars$Dist_edge_std <- scale(plot3vars$Dist_edge, center = F)
plot3vars$Dist_trail_beginning_std <- scale(plot3vars$Dist_trail_beginning, center = F)



#---- Ficheiro com as abundâncias por area de amostragem, para todas as amostras

plot1.all <- read.csv2(here("data.veg","plot1.all.csv"),row.names=1, header = TRUE, sep = ",")
str(plot1.all)
species.veg <- colnames(plot1.all)    ##Estava assim originalmente, mas parece-me trocado
sites.veg <- row.names(plot1.all)
sites.veg

plot2.all <- read.csv2(here("data.veg","plot2.all.csv"),row.names=1, header = TRUE, sep = ",")
str(plot2.all)

plot3.all <- read.csv2(here("data.veg","plot3.all.csv"),row.names=1, header = TRUE, sep = ",")
str(plot3.all)

plot1.inv <- read.csv2(here("data.veg","plot1.inv.csv"),row.names=1, header = TRUE, sep = ",")
str(plot1.inv)
species.inv <- colnames(plot1.inv)    ##Estava assim originalmente, mas parece-me trocado
sites.inv <- row.names(plot1.inv)
sites.inv

plot2.inv <- read.csv2(here("data.veg","plot2.inv.csv"),row.names=1, header = TRUE, sep = ",")
str(plot2.inv)

plot3.inv <- read.csv2(here("data.veg","plot3.inv.csv"),row.names=1, header = TRUE, sep = ",")
str(plot3.inv)


#Taxonomical Alpha



Alpha.all.2017 <- alpha(plot.all, abund=TRUE)
Alpha.all.2019 <- alpha(plot.all[40:65,])
Alpha.inv.2017 <- alpha(plot.inv[1:39,])
Alpha.inv.2019 <- alpha(plot.inv[40:65,])

Alphas.plots  <- cbind (Alpha.all.2017, Alpha.all.2019, Alpha.inv.2017, Alpha.inv.2019, fill=NA)
write.csv(Alphas.plots, file = here("results","Alphas.plots.csv"), row.names = TRUE)



alpha.plot1.all <- alpha(plot1.all)
colnames(alpha.plot1.all) = "alpha.plot1.all" 
alpha.plot2.all <- alpha(plot2.all)
colnames(alpha.plot2.all) = "alpha.plot2.all" 
alpha.plot3.all <- alpha(plot3.all)
colnames(alpha.plot3.all) = "alpha.plot3.all" 

alpha.plot1.inv <- alpha(plot1.inv)
colnames(alpha.plot1.inv) = "alpha.plot1.inv" 
alpha.plot2.inv <- alpha(plot2.inv)
colnames(alpha.plot2.inv) = "alpha.plot2.inv" 
alpha.plot3.inv <- alpha(plot3.inv)
colnames(alpha.plot3.inv) = "alpha.plot3.inv" 

alphas.plot1 <- cbind(alpha.plot1.all,  alpha.plot1.inv )
alphas.plot2 <- cbind( alpha.plot2.all,  alpha.plot2.inv)
alphas.plot3 <- cbind(alpha.plot3.all, alpha.plot3.inv )





###### Taxonomical Beta ----

beta.plot1.all <- beta(plot1.all, abund=TRUE)
beta.plot2.all <- beta(plot2.all, abund=TRUE)
beta.plot3.all <- beta(plot3.all, abund=TRUE)


beta.plot1.inv <- beta(plot1.inv, abund=TRUE)
beta.plot2.inv <- beta(plot2.inv, abund=TRUE)
beta.plot3.inv <- beta(plot3.inv, abund=TRUE)


#### Separating Total, Richness and Replacement  TAXONOMICAL betas into different data frames ----
betatotal.plot1.all <- data.frame(as.matrix(beta.plot1.all[["Btotal"]]), row.names= sites.veg, sep= "tab")
colnames(betatotal.plot1.all) <- sites.veg
betarich.plot1.all <- data.frame(as.matrix(beta.plot1.all[["Brich"]]), row.names= sites.veg,sep= "tab")
colnames(betarich.plot1.all) <- sites.veg
betarepl.plot1.all <- data.frame(as.matrix(beta.plot1.all[["Brepl"]]),row.names= sites.veg, sep= "tab")
colnames(betarepl.plot1.all) <- sites.veg

betatotal.plot2.all <- data.frame(as.matrix(beta.plot2.all[["Btotal"]]), row.names= sites.veg, sep= "tab")
colnames(betatotal.plot2.all) <- sites.veg
betarich.plot2.all <- data.frame(as.matrix(beta.plot2.all[["Brich"]]), row.names= sites.veg,sep= "tab")
colnames(betarich.plot2.all) <- sites.veg
betarepl.plot2.all <- data.frame(as.matrix(beta.plot2.all[["Brepl"]]),row.names= sites.veg, sep= "tab")
colnames(betarepl.plot2.all) <- sites.veg

betatotal.plot3.all <- data.frame(as.matrix(beta.plot3.all[["Btotal"]]), row.names= sites.veg, sep= "tab")
colnames(betatotal.plot3.all) <- sites.veg
betarich.plot3.all <- data.frame(as.matrix(beta.plot3.all[["Brich"]]), row.names= sites.veg,sep= "tab")
colnames(betarich.plot3.all) <- sites.veg
betarepl.plot3.all <- data.frame(as.matrix(beta.plot3.all[["Brepl"]]),row.names= sites.veg, sep= "tab")
colnames(betarepl.plot3.all) <- sites.veg


betatotal.plot1.inv <- data.frame(as.matrix(beta.plot1.inv[["Btotal"]]), row.names= sites.veg, sep= "tab")
colnames(betatotal.plot1.inv) <- sites.veg
betarich.plot1.inv <- data.frame(as.matrix(beta.plot1.inv[["Brich"]]), row.names= sites.veg,sep= "tab")
colnames(betarich.plot1.inv) <- sites.veg
betarepl.plot1.inv <- data.frame(as.matrix(beta.plot1.inv[["Brepl"]]),row.names= sites.veg, sep= "tab")
colnames(betarepl.plot1.inv) <- sites.veg

betatotal.plot2.inv <- data.frame(as.matrix(beta.plot2.inv[["Btotal"]]), row.names= sites.veg, sep= "tab")
colnames(betatotal.plot2.inv) <- sites.veg
betarich.plot2.inv <- data.frame(as.matrix(beta.plot2.inv[["Brich"]]), row.names= sites.veg,sep= "tab")
colnames(betarich.plot2.inv) <- sites.veg
betarepl.plot2.inv <- data.frame(as.matrix(beta.plot2.inv[["Brepl"]]),row.names= sites.veg, sep= "tab")
colnames(betarepl.plot2.inv) <- sites.veg

betatotal.plot3.inv <- data.frame(as.matrix(beta.plot3.inv[["Btotal"]]), row.names= sites.veg, sep= "tab")
colnames(betatotal.plot3.inv) <- sites.veg
betarich.plot3.inv <- data.frame(as.matrix(beta.plot3.inv[["Brich"]]), row.names= sites.veg,sep= "tab")
colnames(betarich.plot3.inv) <- sites.veg
betarepl.plot3.inv <- data.frame(as.matrix(beta.plot3.inv[["Brepl"]]),row.names= sites.veg, sep= "tab")
colnames(betarepl.plot3.inv) <- sites.veg

# Separating the TAXONOMICAL beta values between the Control 250/Max and the other sampling areas from each trail ----

AA1v <- betatotal.plot1.all[1,1]
BB1v <- betatotal.plot1.all[5,c(2:5)]
CC1v <- betatotal.plot1.all[10,c(6:10)]
DD1v <- betatotal.plot1.all[15,c(11:15)]
AA1c <- betatotal.plot1.all[16,16]
BB1c <- betatotal.plot1.all[20,c(17:20)]
CC1c <- betatotal.plot1.all[25,c(21:25)]
plot1.all.btotal <- c(AA1v, BB1v,  CC1, DD1v, AA1c, BB1c, CC1c)

AA1v <- betarich.plot1.all[1,1]
BB1v <- betarich.plot1.all[5,c(2:5)]
CC1v <- betarich.plot1.all[10,c(6:10)]
DD1v <- betarich.plot1.all[15,c(11:15)]
AA1c <- betarich.plot1.all[16,16]
BB1c <- betarich.plot1.all[20,c(17:20)]
CC1c <- betarich.plot1.all[25,c(21:25)]
plot1.all.brich<- c(AA1v, BB1v,  CC1, DD1v, AA1c, BB1c, CC1c)

AA1v <- betarepl.plot1.all[1,1]
BB1v <- betarepl.plot1.all[5,c(2:5)]
CC1v <- betarepl.plot1.all[10,c(6:10)]
DD1v <- betarepl.plot1.all[15,c(11:15)]
AA1c <- betarepl.plot1.all[16,16]
BB1c <- betarepl.plot1.all[20,c(17:20)]
CC1c <- betarepl.plot1.all[25,c(21:25)]
plot1.all.brepl <- c(AA1v, BB1v,  CC1, DD1v, AA1c, BB1c, CC1c)

AA1v <- betatotal.plot2.all[1,1]
BB1v <- betatotal.plot2.all[5,c(2:5)]
CC1v <- betatotal.plot2.all[10,c(6:10)]
DD1v <- betatotal.plot2.all[15,c(11:15)]
AA1c <- betatotal.plot2.all[16,16]
BB1c <- betatotal.plot2.all[20,c(17:20)]
CC1c <- betatotal.plot2.all[25,c(21:25)]
plot2.all.btotal <- c(AA1v, BB1v,  CC1, DD1v, AA1c, BB1c, CC1c)

AA1v <- betarich.plot2.all[1,1]
BB1v <- betarich.plot2.all[5,c(2:5)]
CC1v <- betarich.plot2.all[10,c(6:10)]
DD1v <- betarich.plot2.all[15,c(11:15)]
AA1c <- betarich.plot2.all[16,16]
BB1c <- betarich.plot2.all[20,c(17:20)]
CC1c <- betarich.plot2.all[25,c(21:25)]
plot2.all.brich<- c(AA1v, BB1v,  CC1, DD1v, AA1c, BB1c, CC1c)

AA1v <- betarepl.plot2.all[1,1]
BB1v <- betarepl.plot2.all[5,c(2:5)]
CC1v <- betarepl.plot2.all[10,c(6:10)]
DD1v <- betarepl.plot2.all[15,c(11:15)]
AA1c <- betarepl.plot2.all[16,16]
BB1c <- betarepl.plot2.all[20,c(17:20)]
CC1c <- betarepl.plot2.all[25,c(21:25)]
plot2.all.brepl <- c(AA1v, BB1v,  CC1, DD1v, AA1c,BB1c,  CC1c)

AA1v <- betatotal.plot3.all[1,1]
BB1v <- betatotal.plot3.all[5,c(2:5)]
CC1v <- betatotal.plot3.all[10,c(6:10)]
DD1v <- betatotal.plot3.all[15,c(11:15)]
AA1c <- betatotal.plot3.all[16,16]
BB1c <- betatotal.plot3.all[20,c(17:20)]
CC1c <- betatotal.plot3.all[25,c(21:25)]
plot3.all.btotal <- c(AA1v, BB1v,  CC1, DD1v, AA1c,BB1c, CC1c)

AA1v <- betarich.plot3.all[1,1]
BB1v <- betarich.plot3.all[5,c(2:5)]
CC1v <- betarich.plot3.all[10,c(6:10)]
DD1v <- betarich.plot3.all[15,c(11:15)]
AA1c <- betarich.plot3.all[16,16]
BB1c <- betarich.plot3.all[20,c(17:20)]
CC1c <- betarich.plot3.all[25,c(21:25)]
plot3.all.brich<- c(AA1v, BB1v,  CC1, DD1v, AA1c, BB1c, CC1c)

AA1v <- betarepl.plot3.all[1,1]
BB1v <- betarepl.plot3.all[5,c(2:5)]
CC1v <- betarepl.plot3.all[10,c(6:10)]
DD1v <- betarepl.plot3.all[15,c(11:15)]
AA1c <- betarepl.plot3.all[16,16]
BB1c <- betarepl.plot3.all[20,c(17:20)]
CC1c <- betarepl.plot3.all[25,c(21:25)]
plot3.all.brepl <- c(AA1v, BB1v,  CC1, DD1v, AA1c, BB1c, CC1c)

#Invasives partitioning

AA1v <- betatotal.plot1.inv[1,1]
BB1v <- betatotal.plot1.inv[5,c(2:5)]
CC1v <- betatotal.plot1.inv[10,c(6:10)]
DD1v <- betatotal.plot1.inv[15,c(11:15)]
AA1c <- betatotal.plot1.inv[16,16]
BB1c <- betatotal.plot1.inv[20,c(17:20)]
CC1c <- betatotal.plot1.inv[25,c(21:25)]
plot1.inv.btotal <- c(AA1v, BB1v,  CC1, DD1v, AA1c,BB1c,  CC1c)

AA1v <- betarich.plot1.inv[1,1]
BB1v <- betarich.plot1.inv[5,c(2:5)]
CC1v <- betarich.plot1.inv[10,c(6:10)]
DD1v <- betarich.plot1.inv[15,c(11:15)]
AA1c <- betarich.plot1.inv[16,16]
BB1c <- betarich.plot1.inv[20,c(17:20)]
CC1c <- betarich.plot1.inv[25,c(21:25)]
plot1.inv.brich<- c(AA1v, BB1v,  CC1, DD1v, AA1c,BB1c,  CC1c)

AA1v <- betarepl.plot1.inv[1,1]
BB1v <- betarepl.plot1.inv[5,c(2:5)]
CC1v <- betarepl.plot1.inv[10,c(6:10)]
DD1v <- betarepl.plot1.inv[15,c(11:15)]
AA1c <- betarepl.plot1.inv[16,16]
BB1c <- betarepl.plot1.inv[20,c(17:20)]
CC1c <- betarepl.plot1.inv[25,c(21:25)]
plot1.inv.brepl <- c(AA1v, BB1v,  CC1, DD1v, AA1c, BB1c, CC1c)

AA1v <- betatotal.plot2.inv[1,1]
BB1v <- betatotal.plot2.inv[5,c(2:5)]
CC1v <- betatotal.plot2.inv[10,c(6:10)]
DD1v <- betatotal.plot2.inv[15,c(11:15)]
AA1c <- betatotal.plot2.inv[16,16]
BB1c <- betatotal.plot2.inv[20,c(17:20)]
CC1c <- betatotal.plot2.inv[25,c(21:25)]
plot2.inv.btotal <- c(AA1v, BB1v,  CC1, DD1v, AA1c, BB1c, CC1c)

AA1v <- betarich.plot2.inv[1,1]
BB1v <- betarich.plot2.inv[5,c(2:5)]
CC1v <- betarich.plot2.inv[10,c(6:10)]
DD1v <- betarich.plot2.inv[15,c(11:15)]
AA1c <- betarich.plot2.inv[16,16]
BB1c <- betarich.plot2.inv[20,c(17:20)]
CC1c <- betarich.plot2.inv[25,c(21:25)]
plot2.inv.brich<- c(AA1v, BB1v,  CC1, DD1v, AA1c, BB1c, CC1c)

AA1v <- betarepl.plot2.inv[1,1]
BB1v <- betarepl.plot2.inv[5,c(2:5)]
CC1v <- betarepl.plot2.inv[10,c(6:10)]
DD1v <- betarepl.plot2.inv[15,c(11:15)]
AA1c <- betarepl.plot2.inv[16,16]
BB1c <- betarepl.plot2.inv[20,c(17:20)]
CC1c <- betarepl.plot2.inv[25,c(21:25)]
plot2.inv.brepl <- c(AA1v, BB1v,  CC1, DD1v, AA1c, BB1c, CC1c)

AA1v <- betatotal.plot3.inv[1,1]
BB1v <- betatotal.plot3.inv[5,c(2:5)]
CC1v <- betatotal.plot3.inv[10,c(6:10)]
DD1v <- betatotal.plot3.inv[15,c(11:15)]
AA1c <- betatotal.plot3.inv[16,16]
BB1c <- betatotal.plot3.inv[20,c(17:20)]
CC1c <- betatotal.plot3.inv[25,c(21:25)]
plot3.inv.btotal <- c(AA1v, BB1v,  CC1, DD1v, AA1c, BB1c, CC1c)

AA1v <- betarich.plot3.inv[1,1]
BB1v <- betarich.plot3.inv[5,c(2:5)]
CC1v <- betarich.plot3.inv[10,c(6:10)]
DD1v <- betarich.plot3.inv[15,c(11:15)]
AA1c <- betarich.plot3.inv[16,16]
BB1c <- betarich.plot3.inv[20,c(17:20)]
CC1c <- betarich.plot3.inv[25,c(21:25)]
plot3.inv.brich<- c(AA1v, BB1v,  CC1, DD1v, AA1c, BB1c, CC1c)

AA1v <- betarepl.plot3.inv[1,1]
BB1v <- betarepl.plot3.inv[5,c(2:5)]
CC1v <- betarepl.plot3.inv[10,c(6:10)]
DD1v <- betarepl.plot3.inv[15,c(11:15)]
AA1c <- betarepl.plot3.inv[16,16]
BB1c <- betarepl.plot3.inv[20,c(17:20)]
CC1c <- betarepl.plot3.inv[25,c(21:25)]
plot3.inv.brepl <- c(AA1v, BB1v,  CC1, DD1v, AA1c, BB1c, CC1c)


## BETA TEMPORAL ############

guilmon <- betatotal.plot1.all[1,16]
ser0 <- betatotal.plot1.all[2,17]
ser50 <- betatotal.plot1.all[3,18]
ser250 <- betatotal.plot1.all[4,19]
sercon <- betatotal.plot1.all[5,20]
mn0 <- betatotal.plot1.all[6,21]
mn50 <- betatotal.plot1.all[7,22]
mn250 <- betatotal.plot1.all[8,23]
mnmax <- betatotal.plot1.all[9,24]
mncon <- betatotal.plot1.all[10,25]

beta.temp <- rbind (guilmon, ser0, ser50, ser250, sercon, mn0, mn50, mn250, mnmax, mncon )
rownames(beta.temp) = rbind ("guilmon", "ser0", "ser50", "ser250", "sercon", "mn0", "mn50", "mn250", "mnmax", "mncon" )

#
#### COMPILING ALL BETA INFORMATION INTO A TABLE, AND EXPORTING IT TO A FILE ----

betas.plot1 <- as.data.frame(t(rbind(plot1.all.btotal, plot1.all.brich, plot1.all.brepl,plot1.inv.btotal, plot1.inv.brich, plot1.inv.brepl)))
rownames(betas.plot1) = sites.veg
betas.plot2 <- as.data.frame(t(rbind(plot2.all.btotal, plot2.all.brich, plot2.all.brepl,plot2.inv.btotal, plot2.inv.brich, plot2.inv.brepl )))
rownames(betas.plot2) = sites.veg
betas.plot3 <- as.data.frame(t(rbind(plot3.all.btotal, plot3.all.brich, plot3.all.brepl,plot3.inv.btotal, plot3.inv.brich, plot3.inv.brepl  )))
rownames(betas.plot3) = sites.veg




rownames(betas.veg) = sites.veg
betas.veg[betas.veg == 0] <- 0.001
betas.veg[betas.veg == 1] <- 0.999
str(betas.veg)
str(sites.veg)
###########################################################################
#                               RESULTS                                ####
###########################################################################

Results.plot1 <- cbind.data.frame(plot1vars, alphas.plot1 , betas.plot1)

Results.plot2 <- cbind.data.frame(plot2vars, alphas.plot2, betas.plot2)

Results.plot3 <- cbind.data.frame(plot3vars, alphas.plot3, betas.plot3)




#MAKING THE RESULTS EXPORTABLE INTO CSV
Results.plot1 <- apply(Results.plot1, 2 , as.character, header=TRUE)
Results.plot2 <- apply(Results.plot2, 2 , as.character, header=TRUE)
Results.plot3 <- apply(Results.plot3, 2 , as.character, header=TRUE)

#NAMING THE TRAIL SEGMENTS
rownames(Results.plot1) <- rownames(betas.veg)
rownames(Results.plot2) <- rownames(betas.veg)
rownames(Results.plot3) <- rownames(betas.veg)
#PASSING RESULTS TO FILE
write.csv(Results.plot1, file = here("results","RESULTS.PLOT1.csv"), row.names = TRUE)
write.csv(Results.plot2, file = here("results","RESULTS.PLOT2.csv"), row.names = TRUE)
write.csv(Results.plot3, file = here("results","RESULTS.PLOT3.csv"), row.names = TRUE)



Results.plot1.import <- read.csv2(here("results","RESULTS.PLOT1.CSV"), header=TRUE, row.names = 1,  stringsAsFactors = T,sep = ",", dec = ".")
Results.plot2.import <- read.csv2(here("results","RESULTS.PLOT2.CSV"), header=TRUE, row.names = 1,  stringsAsFactors = T,sep = ",", dec = ".")
Results.plot3.import <- read.csv2(here("results","RESULTS.PLOT3.CSV"), header=TRUE, row.names = 1,  stringsAsFactors = T,sep = ",", dec = ".")


###########################################################################
#                           Generating models                           ###
###########################################################################

## Fake database

dummy1 <- read.csv2(here("data.veg","dummy3.csv"), header=TRUE, row.names = 1,  stringsAsFactors = T,sep = ",", dec = ".")

alpha.all <- dredge(glmmTMB(resp ~ dist_begin + dist_trail + (1 | ID), data= dummy1 , family = "poisson"))
alpha.all
alpha.all1 <- glmmTMB(resp ~ dist_begin + dist_trail + (1 | ID), data= dummy1 , family = "poisson")
summary(alpha.all1)
performance::r2(alpha.all1)
version2(alpha.all1)


alpha.all2 <- glm(resp ~ dist_begin , data= dummy1 , family = "poisson")
summary(alpha.all2)
version3(alpha.all2)
performance::r2(alpha.all2)

alpha.all3 <- glmmTMB(resp ~  dist_trail + (1 | ID), data= dummy , family = "poisson")
summary(alpha.all3)
performance::r2(alpha.all3)
version3(alpha.all3)

alpha.all4 <- glmmTMB(resp ~  dist_begin + (1 | ID), data= dummy , family = "poisson")
summary(alpha.all4)
performance::r2(alpha.all4)
version3(alpha.all4)






## Uploading results for vegetation

test1veg <- read.csv2(here("data.veg","test1veg.CSV"), header=TRUE, row.names = 1,  stringsAsFactors = T,sep = ",", dec = ".")
test1veg$percent.alfa[test1veg$percent.alfa == 0] <- 0.001
test1veg$percent.alfa[test1veg$percent.alfa == 1] <- 0.999
test1veg$percent.abund[test1veg$percent.abund == 0] <- 0.001
test1veg$percent.abund[test1veg$percent.abund == 1] <- 0.999


alpha.all <- dredge(glmmTMB(alpha.all ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data= test1veg , family = "poisson"))
alpha.all
alpha.all1 <- glmer(alpha.all ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std  + (1 | ForestID), data= test1veg , family = "poisson")
summary (alpha.all1)
version3(alpha.all1)
performance::r2(alpha.all1)
alpha.inv <- dredge(glmmTMB(alpha.inv ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data= test1veg , family = "poisson"))

percent.alfa <- dredge(glmmTMB(percent.alfa ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data= test1veg , family = "beta_family"))
abund.all <- dredge(glmmTMB(abund.all ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data= test1veg , family = "poisson"))

abund.all1 <- glmmTMB(abund.all ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data= test1veg , family = "poisson")
summary(abund.all1 )
performance::r2(alpha.all1)

abund.all2 <- glmmTMB(abund.all ~ Dist_trail_std  + (1 | ForestID), data= test1veg , family = "poisson")
summary(abund.all2 )
performance::r2(alpha.all2)



abund.inv <- dredge(glmmTMB(abund.inv ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data= test1veg , family = "poisson"))
percent.abund <- dredge(glmmTMB(percent.abund ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data= test1veg , family = "beta_family"))

#R2
Ralpha.all <- performance::r2(glmmTMB(alpha.all ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data= test1veg , family = "poisson"))
Ralpha.inv <- version2(glmmTMB(alpha.inv ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data= test1veg , family = "poisson"))
#Rpercent.alfa <- version2(glmmTMB(percent.alfa ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data= test1veg , family = "beta_family"))
Rabund.all <- version2(glmmTMB(abund.all ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data= test1veg , family = "poisson"))
Rabund.inv <- version2(glmmTMB(abund.inv ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data= test1veg , family = "poisson"))
#Rpercent.abund <- version2(glmmTMB(percent.abund ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data= test1veg , family = "beta_family"))


test2veg <- read.csv2(here("data.veg","test2veg.CSV"), header=TRUE, row.names = 1,  stringsAsFactors = T,sep = ",", dec = ".")

test2veg$delta.percent.alpha[test2veg$delta.percent.alpha == 0] <- 0.001
test2veg$delta.percent.alpha[test2veg$delta.percent.alpha == 1] <- 0.999
test2veg$delta.percent.abund[test2veg$delta.percent.abund == 0] <- 0.001
test2veg$delta.percent.abund[test2veg$delta.percent.abund == 1] <- 0.999

delta.alpha <- dredge(glmmTMB(delta.alpha ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data= test2veg , family = "gaussian"))
delta.abund <- dredge(glmmTMB(delta.abund ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data= test2veg , family = "gaussian"))
delta.percent.alpha <- dredge(glmmTMB(delta.percent.alpha ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data= test2veg , family = "gaussian"))
delta.percent.abund <- dredge(glmmTMB(delta.percent.abund ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data= test2veg , family = "gaussian"))

#Rdelta.alpha <- version2(glmmTMB(delta.alpha ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data= test2veg , family = "gaussian"))
#Rdelta.abund <- version2(glmmTMB(delta.abund ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data= test2veg , family = "gaussian"))
Rdelta.percent.alpha <- version2(glmmTMB(delta.percent.alpha ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data= test2veg , family = "gaussian"))
#Rdelta.percent.abund <- version2(glmmTMB(delta.percent.abund ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data= test2veg , family = "gaussian"))

fam.veg <- c(NA, NA, NA, NA,NA, NA,NA, "poisson", "poisson","beta_family","beta_family", "beta_family", "beta_family","beta_family", "beta_family")

## GENERATING MODELS FOR PLOT 1

# Alpha

dredge.alpha.plot1 <- dredge(glmmTMB(alpha.plot1.all ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data= Results.plot1.import , family = "poisson"))
dredge.alpha.plot2 <- dredge(glmmTMB(alpha.plot2.all ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data= Results.plot2.import , family = "poisson"))
dredge.alpha.plot3 <- dredge(glmmTMB(alpha.plot3.all ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data= Results.plot3.import , family = "poisson"))
dredge.alpha.plot1.inv <- dredge(glmmTMB(alpha.plot1.inv ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data= Results.plot1.import , family = "poisson"))
dredge.alpha.plot2.inv <- dredge(glmmTMB(alpha.plot2.inv ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data= Results.plot2.import , family = "poisson"))
dredge.alpha.plot3.inv <- dredge(glmmTMB(alpha.plot3.inv ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data= Results.plot3.import , family = "poisson"))

#dredge.alpha.inv <- data.frame(dredge.alpha.inv)
#write.csv(dredge.alpha.inv, file = here("results","test.alpha.inv.csv"), row.names = TRUE)

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