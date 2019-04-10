## R code written by Elizabeth A. Bowman Mar. 6, 2018
## University of Arizona, School of Plant Sciences, eabowman@email.arizona.edu
## Analyses to evaluate whether any fungal OTUs were significantly associated 
## fire history or range

#========================================================================================#
  # Load data and libraries----
#========================================================================================#

#----------------------------------------------------------------------------------------#
# Load libraries----
#----------------------------------------------------------------------------------------#
library(ggplot2);library (tidyr);library (vegan);library (dplyr); library(indicspecies)

#----------------------------------------------------------------------------------------#
# set up paths to directories----
#----------------------------------------------------------------------------------------#
#--path to directory where the climate data downloaded from BIOCLIM site is stored
dat.dir <- "~/Documents/PhD/2_EM_Fire_effect/data/"
fig.dir <- '~/Documents/PhD/2_EM_Fire_effect/figures_output/'
res.dir <- "~/Documents/PhD/2_EM_Fire_effect/results_output/"

#----------------------------------------------------------------------------------------#
# Load data and clean up: Tree level----
#----------------------------------------------------------------------------------------#
#--Load data
all.data <- read.csv(paste0(dat.dir,'20170806_OTU_data.csv'), as.is = T)

#--Remove sequences not assigned to an OTU
otu.data <- all.data[!is.na(all.data$otu.97), ]

#--Remove sequences not assigned to Ponderosa host
otu.data <- otu.data[otu.data$Host == 'Ponderosa',]

#--Add an OTU count column (1 sequence = 1 frequency count)
otu.data['otu.count'] <- 1

#--Creates matrix, grouping OTUs by tree number and getting a sum of EM tip abundance 
# for each OTU
otu.data %>%
  select(Sample_name,Burn_status,Range,Site,Tree,otu.97,otu.count) %>%
  spread(otu.97,otu.count) %>%
  group_by(Tree,Site,Range,Burn_status) %>%
  select(matches ("otu*")) %>%
  summarize_each(funs(sum(., na.rm = TRUE))) %>%
  as.data.frame() -> stsp.matrix

#--Creates matrix, grouping OTUs by tree number and getting a sum of EM tip abundance 
# for each OTU
otu.data %>%
  select(Sample_name,Burn_status,Range,Site,Tree,otu.97,Tip_count) %>%
  spread(otu.97,Tip_count) %>%
  group_by(Tree,Site,Range,Burn_status) %>%
  select(matches ("otu*")) %>%
  summarize_each(funs(sum(., na.rm = TRUE))) %>%
  as.data.frame() -> abund.matrix

for(i in 1:nrow(stsp.matrix)){
  stsp.matrix[i, 'burn_range'] <- 
    paste(stsp.matrix[i, 'Burn_status'], stsp.matrix[i, 'Range'])
}

colnames(stsp.matrix)
stsp.matrix <- stsp.matrix[c(1:4, 121, 5:120)]

for(i in 1:nrow(abund.matrix)){
  abund.matrix[i, 'burn_range'] <- 
    paste(abund.matrix[i, 'Burn_status'], abund.matrix[i, 'Range'])
}

colnames(abund.matrix)
abund.matrix <- abund.matrix[c(1:4, 121, 5:120)]

#========================================================================================#
# Indicator species analysis----
#========================================================================================#

#----------------------------------------------------------------------------------------#
# OTU frequency based
#----------------------------------------------------------------------------------------#
#<< Overall EM, with singletons included >>-----------------------------------------------
#--isolate OTUs with great than 4 ocurrences
em.indsp <- stsp.matrix[6:length(stsp.matrix)]
em.indsp <- em.indsp[, colSums(em.indsp) > 4]

#--pull out fire history for grouping
fire_comm <- stsp.matrix$Burn_status
#--pull out mountain range for grouping
range_comm <- stsp.matrix$Range
#--pull out both burn and range information
burn_range_comm <- stsp.matrix$burn_range

#<< Fire history >> -----------------
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

#<< Range >> -----------------
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

#----------------------------------------------------------------------------------------#
# Abundance
#----------------------------------------------------------------------------------------#
#--isolate OTUs with great than 4 ocurrences
em.indsp <- stsp.matrix[5:length(abund.matrix)]
em.indsp <- em.indsp[, colSums(em.indsp) > 4]

#--pull out fire history for grouping
fire_comm = abund.matrix$Burn_status
#--pull out mountain range for grouping
range_comm = abund.matrix$Range

#<< Fire history >> -----------------
#--Association between species patterns and high/low elevation sites
indsp.ab.fire.em <- multipatt(em.indsp, fire_comm, func = "IndVal") # uses default func = "IndVal"
#--print summary
summary(indsp.ab.fire.em)
#--Components A and B. 
summary(indsp.ab.fire.em, indvalcomp = TRUE) 
# A = the probability that the surveyed site belongs to the target site group given the 
# fact that the species has been found.
# B = the probability of finding the species in sites belonging to the site group.

#--shows all species with lower indval, not just sig.
summary(indsp.ab.fire.em, alpha = 1)
#--shows all spp. even those with highest indval that do not show up in other summaries
indsp.ab.fire.em$sign 

#<< Range >> -----------------
#--Association between species patterns and high/low elevation sites
indsp.ab.range.em <- multipatt(em.indsp, range_comm, func = "IndVal") # uses default func = "IndVal"
#--print summary
summary(indsp.ab.range.em)
#--Components A and B. 
summary(indsp.ab.range.em, indvalcomp = TRUE) 
# A = the probability that the surveyed site belongs to the target site group given the 
# fact that the species has been found.
# B = the probability of finding the species in sites belonging to the site group.

#--shows all species with lower indval, not just sig.
summary(indsp.ab.range.em, alpha = 1)
#--shows all spp. even those with highest indval that do not show up in other summaries
indsp.ab.range.em$sign 
