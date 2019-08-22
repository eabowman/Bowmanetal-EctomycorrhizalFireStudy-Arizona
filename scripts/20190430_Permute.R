## Script created by Liz Bowman April 30, 2019
## for creating a stratified model for analyses using the Permute package

#========================================================================================#
# Load data and libraries----
#========================================================================================#
#install.packages('permute')
library(permute)

#========================================================================================#
# Permutations----
#========================================================================================#
# sites/trees

#--Lowest level (within plot)
within.em <- Within(type = 'free')

#--Plot/site level
stsp.matrix$Site <- as.factor(stsp.matrix$Site)
plts.em <- Plots(strata = stsp.matrix$Site, type = 'free')

#--Blocks
#range
stsp.matrix$Range <- as.factor(stsp.matrix$Range)
blck.range <- Blocks(strata = stsp.matrix$Range)

#fire history
stsp.matrix$Burn_status <- as.factor(stsp.matrix$Burn_status)
blck.burn <- Blocks(strata = stsp.matrix$Burn_status)

#<< permutation design for testing effect of fire history >> ---------------------------
how.range <- how(within = within.em, plots = plts.em, blocks = blck.range)

#<< permutation design for testing effect of range >> ----------------------------------
how.burn <- how(within = within.em, plots = plts.em, blocks = blck.burn)

#<< permutation design for testing effect of range and fire history>> ------------------
h1 <- how(within = within.em, plots = plts.em)

