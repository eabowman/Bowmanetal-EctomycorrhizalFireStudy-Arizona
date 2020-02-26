## R code written by Elizabeth A. Bowman May 22, 2019
## University of Arizona, School of Plant Sciences, eabowman@email.arizona.edu
## Analyses to evaluate whether any fungal OTUs were significantly associated 
## fire history or range

#========================================================================================#
# Load data ----
#========================================================================================#

#--Load data
stsp.matrix <- read.csv(paste0(dat.dir,'97%_SitexSpecies_TipAb.csv'), as.is = T)

#--Separate out by range
scm.stsp <- stsp.matrix[stsp.matrix$Range == 'santa.catalina',]
pm.stsp <- stsp.matrix[stsp.matrix$Range == 'pinaleno',]

#--Separate out by fire history
fa.stsp <- stsp.matrix[stsp.matrix$Burn_status == 'burned',]
fu.stsp <- stsp.matrix[stsp.matrix$Burn_status == 'unburned',]

#========================================================================================#
# Indicator species analysis----
#========================================================================================#

#<< All sites >>---- 
#--isolate OTUs with great than 4 ocurrences
em.indsp <- stsp.matrix[12:length(stsp.matrix)]
em.indsp <- em.indsp[, colSums(em.indsp) > 4]

#--pull out fire history for grouping
fire_comm <- stsp.matrix$Burn_status
#--pull out mountain range for grouping
range_comm <- stsp.matrix$Range
#--pull out both burn and range information
burn_range_comm <- stsp.matrix$burn_range

#--Fire history ----
#--Association between species patterns and high/low elevation sites
indsp.fire.em <- multipatt(em.indsp, fire_comm, func = "IndVal") # uses default func = "IndVal"
#--print summary
summary(indsp.fire.em)
#--Components A and B. 
summary(indsp.fire.em, indvalcomp = TRUE) 
# A = the probability that the surveyed site belongs to the target site group given the 
# fact that the species has been found.
# B = the probability of finding the species in sites belonging to the site group.

#--shows all species with lower indval, not just sig.
summary(indsp.fire.em, alpha = 1)
#--shows all spp. even those with highest indval that do not show up in other summaries
indsp.fire.em$sign 

#<< Range >> ----
#--Association between species patterns and high/low elevation sites
indsp.range.em <- multipatt(em.indsp, range_comm, func = "IndVal") # uses default func = "IndVal"
#--print summary
summary(indsp.range.em)
#--Components A and B. 
summary(indsp.range.em, indvalcomp = TRUE) 
# A = the probability that the surveyed site belongs to the target site group given the 
# fact that the species has been found.
# B = the probability of finding the species in sites belonging to the site group.

#--shows all species with lower indval, not just sig.
summary(indsp.range.em, alpha = 1)
#--shows all spp. even those with highest indval that do not show up in other summaries
indsp.range.em$sign 

#<< Fire history within Santa Catalina Mts. >>----
scm.indsp <- scm.stsp[12:length(scm.stsp)]
scm.indsp <- scm.indsp[, colSums(scm.indsp) > 4]

#--pull out fire history for grouping
fire_comm <- scm.stsp$Burn_status

#--Association between species patterns and high/low elevation sites
indsp.range.em <- multipatt(scm.indsp, fire_comm, func = "IndVal") # uses default func = "IndVal"
#--print summary
summary(indsp.range.em)
#--Components A and B. 
summary(indsp.range.em, indvalcomp = TRUE) 
# A = the probability that the surveyed site belongs to the target site group given the 
# fact that the species has been found.
# B = the probability of finding the species in sites belonging to the site group.

#--shows all species with lower indval, not just sig.
summary(indsp.range.em, alpha = 1)
#--shows all spp. even those with highest indval that do not show up in other summaries
indsp.range.em$sign 

#<< Fire history within Pinaleno Mts. >>----
pm.indsp <- pm.stsp[12:length(pm.stsp)]
pm.indsp <- pm.indsp[, colSums(pm.indsp) > 4]

#--pull out fire history for grouping
fire_comm <- pm.stsp$Burn_status

#--Association between species patterns and high/low elevation sites
indsp.range.em <- multipatt(pm.indsp, fire_comm, func = "IndVal") # uses default func = "IndVal"
#--print summary
summary(indsp.range.em)
#--Components A and B. 
summary(indsp.range.em, indvalcomp = TRUE) 
# A = the probability that the surveyed site belongs to the target site group given the 
# fact that the species has been found.
# B = the probability of finding the species in sites belonging to the site group.

#--shows all species with lower indval, not just sig.
summary(indsp.range.em, alpha = 1)
#--shows all spp. even those with highest indval that do not show up in other summaries
indsp.range.em$sign 
