
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

