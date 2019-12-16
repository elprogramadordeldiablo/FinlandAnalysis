
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



gm1 <- glmmTMB(TAlphaAll ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data = Results2 , family = poisson)
dredgm1 <- dredge(gm1)
dredgm1







(all.tax.btotal, all.tax.brich, all.tax.brepl, nat.tax.btotal, nat.tax.brich, nat.tax.brepl, nind.tax.btotal, nind.tax.brich, nind.tax.brepl, all.func.btotal, all.func.brich, 
all.func.brepl, nat.func.btotal, nat.func.brich, nat.func.brepl, nind.func.btotal, nind.func.brich, nind.func.brepl)