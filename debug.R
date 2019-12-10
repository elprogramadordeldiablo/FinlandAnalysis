gm6 <- glmmTMB(FAlphaNInd ~ Dist_edge_std + Dist_trail_std + Dist_trail_beginning_std + (1 | ForestID), data = Results2 , family = Gamma)




dredge(gm6)

