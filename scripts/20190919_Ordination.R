## Script created by Liz Bowman June 3, 2017
## for analyzing community differences between ranges and burn status of sites

#========================================================================================#
# Load data and libraries----
#========================================================================================#

#----------------------------------------------------------------------------------------#
# Load libraries----
#----------------------------------------------------------------------------------------#
library(ggplot2);library (tidyr);library (vegan);library (dplyr)

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
  select(Sample_name,Burn_status,Range,Site,Tree,otu.97,otu.count,Tip_count) %>%
  spread(otu.97,Tip_count) %>%
  group_by(Tree,Site,Range,Burn_status) %>%
  select(matches ("otu*")) %>%
  summarize_each(funs(sum(., na.rm = TRUE))) %>%
  as.data.frame() -> stsp.matrix

#--climate data: add to stsp.matrix
clim.data <- read.csv(paste0(dat.dir, 'climate_data.csv'))
for(i in unique(stsp.matrix$Site)) {
  stsp.matrix[stsp.matrix$Site == i, 'prec'] <- clim.data[clim.data$site == i, 'prec']
  stsp.matrix[stsp.matrix$Site == i, 'Tmax'] <-
    clim.data[clim.data$site == i, 'Tmax']
}

#--add forest type
stsp.matrix$forest <- NA
stsp.matrix[stsp.matrix$site %in% c('sc1','sc2','sc3','sc4','p2'), 'forest'] <- 'pine oak'
stsp.matrix[stsp.matrix$site %in% c('p3','p4'), 'forest'] <- 'pine df'
stsp.matrix[stsp.matrix$site == 'p1', 'forest'] <- 'pine'

#--reorder matrix
colnames(stsp.matrix)
stsp.matrix <- stsp.matrix[c(1:4,121:122,130,5:120)]

#----------------------------------------------------------------------------------------#
# Taxonomic data
#----------------------------------------------------------------------------------------#

tax.data <- read.csv(paste0(dat.dir,'20170806_OTU_data.csv'), as.is = T)
tax.data <- tax.data[tax.data$Host == 'Ponderosa',]

#--Add column for range burn combo
for (r in unique(tax.data$Range)){
  for (b in unique(tax.data$Burn_status)){
    tax.data[tax.data$Range == r & tax.data$Burn_status == b, 'rangeburn'] <- paste(r,b)
  }
}

#----------------------------------------------------------------------------------------#
# Create result tables
#----------------------------------------------------------------------------------------#
#--Table of Anosim results
anosim.res <- c("jaccard.all", "morisita.all",
                'jaccard.p', 'morisita.p',
                'jaccard.scm', 'morisita.scm',
                'morisita.fa', 'jaccard.fa',
                'morisita.fu', 'jaccard.fu')
anosim.res <- data.frame(anosim.res)
#--Table of PERMANOVA results
permanova.res <- data.frame(c('jaccard','morisita'))
colnames(permanova.res) <- "test"

#--outlier from Santa Catalina Mts. (LB056)***
stsp.matrix <- stsp.matrix[!stsp.matrix$Tree == 'LB056',]

#----------------------------------------------------------------------------------------#
# Create distinct tables for each burn history/range
#----------------------------------------------------------------------------------------#
FA.matrix <- stsp.matrix[stsp.matrix$Burn_status == 'burned',]
FU.matrix <- stsp.matrix[stsp.matrix$Burn_status == 'unburned',]

scm.matrix <- stsp.matrix[stsp.matrix$Range == 'santa.catalina',]
pm.matrix <- stsp.matrix[stsp.matrix$Range == 'pinaleno',]

#========================================================================================#
# Jaccard based dissimilarity index: Pinaleno----
#========================================================================================#

#--Remove outlier, NF16
pm.matrix <- pm.matrix[!pm.matrix$Tree == 'NF16',]

#--isolate otu data
comm.matrix <- pm.matrix[9:length(pm.matrix)]

#--comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 2]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 1, ] # remove rows with sums of 0
jaccard.matrix <- pm.matrix[row.names(comm.matrix),]

#--distance matrix using jaccard index
comm.dist.jaccard <- vegdist(comm.matrix, method = "jaccard", binary = TRUE)

#--NMDS analysis
jaccard.otu <- metaMDS(comm.dist.jaccard, dist = "bray", permutations = 999,
                       try = 100, trymax = 1000)

#--Stress
jaccard.otu$stress

#--add stress of NMDS to results table
anosim.res[which(anosim.res$anosim.res =="jaccard.p"), "stress.nmds"] <-
  jaccard.otu$stress

#--format and output NMDS plots to figure folder
# jpeg(filename = 'figures_output/NMDS_pinaleno_Jaccard.jpeg',
#      width = 700, height = 600,
#      quality = 100)
# par(mfrow = c(1,1), "mar"=c(6, 5, 5, 3))
# 
# #--Plot NMDS of EM community based on Jaccard index and OTU abundance
# plot(jaccard.otu, display = "sites", type = "n", cex.lab = 2,
#      cex.axis = 1.5, xlab = 'Axis 1', ylab = 'Axis 2')
# # colors for points
# color.vec <- data.frame(color = rep(NA, nrow(jaccard.matrix)),
#                         p.group = jaccard.matrix$Burn_status)
# color.vec[color.vec$p.group == 'burned', 'color'] <- 'black'
# color.vec[color.vec$p.group == 'unburned', 'color'] <- 'darkgreen'
# 
# color.vec[color.vec$p.group == 'burned', 'shape'] <- 15
# color.vec[color.vec$p.group == 'unburned', 'shape'] <- 15
# 
# #ordipointlabel(jaccard.otu, display = "sites")
# #--isolate points for pinaleno mts.
# burned <- row.names(jaccard.matrix[jaccard.matrix$Burn_status == 'burned',])
# unburned <- row.names(jaccard.matrix[jaccard.matrix$Burn_status == 'unburned',])
# 
# #--isolate points using rows isolated above
# points(jaccard.otu$points[burned,1:2], display = "sites", cex = 3,
#        pch = color.vec[color.vec$p.group == 'burned',
#                        'shape'],
#        col = color.vec[color.vec$p.group == 'burned',
#                        'color'],
#        bg = color.vec[color.vec$p.group == 'burned',
#                       'color'])
# points(jaccard.otu$points[unburned,1:2], display = "sites", cex = 3,
#        pch = color.vec[color.vec$p.group == 'unburned',
#                        'shape'],
#        col = color.vec[color.vec$p.group == 'unburned',
#                        'color'],
#        bg = color.vec[color.vec$p.group == 'unburned',
#                       'color'])
# # legend("bottomright", legend = c('FA',
# #                              'FU'), bty = "n",
# #        col = c('black', 'darkgreen'),
# #        pch = c(15,15),
# #        pt.bg =  c('black','darkgreen'), cex = 2)
# 
# #--Ordihull variations by range, burn, and both burn and range
# #ordihull(jaccard.otu, groups = color.vec$p.group)
# ordiellipse(jaccard.otu, groups = jaccard.matrix$Burn_status,
#             col = c('black', 'black'),
#             kind = 'ehull')
# #ordihull(jaccard.otu, groups = stsp.matrix$Range) # just mt. range
# #--Overlays
# fit <- envfit(jaccard.otu ~ prec, jaccard.matrix)
# fit
# # plot(fit, type = 'n')
# 
# dev.off()

#--BetaDisper
#--a multivariate analogue of Levene's test for homogeneity of variance
betadisper <- betadisper(comm.dist.jaccard, group = jaccard.matrix$Burn_status)
jaccard.betadisper <- anova(betadisper)
jaccard.betadisper

#--add Betadisper results to table
anosim.res[which(anosim.res$anosim.res =="jaccard.p"), "F.betadisper"] <-
  jaccard.betadisper$`F value`[1]
anosim.res[which(anosim.res$anosim.res =="jaccard.p"), "df.betadisper.1"] <-
  jaccard.betadisper$Df[1]
anosim.res[which(anosim.res$anosim.res =="jaccard.p"), "df.betadisper.2"] <-
  jaccard.betadisper$Df[2]
anosim.res[which(anosim.res$anosim.res =="jaccard.p"), "p.betadisper"] <-
  jaccard.betadisper$`Pr(>F)`[1]

#<< ANOSIM >>-----------------------------------------------------------------------------
jaccard.anosim <- anosim(comm.matrix, grouping = jaccard.matrix$Burn_status,
                         distance = "jaccard")
jaccard.anosim

#--Add results to data frame
anosim.res[which(anosim.res$anosim.res =="jaccard.p"), "r"] <- jaccard.anosim$statistic
anosim.res[which(anosim.res$anosim.res =="jaccard.p"), "p"] <- jaccard.anosim$signif

#========================================================================================#
# Morisita based dissimilarity index: Pinaleno----
#========================================================================================#

# #--comment to not remove outlier (NF19)
# stsp.matrix <- stsp.matrix[!stsp.matrix$Tree %in% c('NF19'),]

#--isolate site X species matrix only without metadata
comm.matrix <- pm.matrix[9:length(pm.matrix)]

#--comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 2]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 1, ] # remove rows with sums of 0
morisita.matrix <- pm.matrix[row.names(comm.matrix),]

#--distance matrix using jaccard index
comm.dist.morisita <- vegdist(comm.matrix, method = "horn", binary = F)

#--NMDS analysis
morisita.otu <- metaMDS(comm.dist.morisita, dist = "bray", permutations = 999,
                        try = 100, trymax = 1000)

#--add stress of NMDS to results table
anosim.res[which(anosim.res$anosim.res =="morisita.p"), "stress.nmds"] <-
  morisita.otu$stress

#--format and output NMDS plots to figure folder
# jpeg(filename = 'figures_output/NMDS_pinaleno_MorisitaHorn.jpeg',
#      width = 700, height = 600,
#      quality = 100)
# par(mfrow = c(1,1), "mar"=c(6, 5, 5, 3))
# 
# #--Plot NMDS of EM community based on Jaccard index and OTU abundance
# plot(morisita.otu, display = "sites", type = "n", cex.lab = 2,
#      cex.axis = 1.5, xlab = 'Axis 1', ylab = 'Axis 2')
# 
# # colors for points
# color.vec <- data.frame(color = rep(NA, nrow(morisita.matrix)),
#                         p.group = morisita.matrix$Burn_status)
# color.vec[color.vec$p.group == 'burned', 'color'] <- 'black'
# color.vec[color.vec$p.group == 'unburned', 'color'] <- 'darkgreen'
# 
# color.vec[color.vec$p.group == 'burned', 'shape'] <- 15
# color.vec[color.vec$p.group == 'unburned', 'shape'] <- 15
# 
# #ordipointlabel(jaccard.otu, display = "sites")
# #--isolate points for pinaleno mts.
# burned <- row.names(morisita.matrix[jaccard.matrix$Burn_status == 'burned',])
# unburned <- row.names(morisita.matrix[jaccard.matrix$Burn_status == 'unburned',])
# 
# #--isolate points using rows isolated above
# points(morisita.otu$points[burned,1:2], display = "sites", cex = 3,
#        pch = color.vec[color.vec$p.group == 'burned',
#                        'shape'],
#        col = color.vec[color.vec$p.group == 'burned',
#                        'color'],
#        bg = color.vec[color.vec$p.group == 'burned',
#                       'color'])
# points(morisita.otu$points[unburned,1:2], display = "sites", cex = 3,
#        pch = color.vec[color.vec$p.group == 'unburned',
#                        'shape'],
#        col = color.vec[color.vec$p.group == 'unburned',
#                        'color'],
#        bg = color.vec[color.vec$p.group == 'unburned',
#                       'color'])
# # legend("bottomleft", legend = c('FA',
# #                              'FU'), bty = "n",
# #        col = c('black', 'darkgreen'),
# #        pch = c(16,16),
# #        pt.bg =  c('black','darkgreen'), cex = 2)
# 
# #--Ordihull variations by range, burn, and both burn and range
# ordiellipse(morisita.otu, groups = morisita.matrix$Burn_status,
#             col = c('black', 'black'),
#             kind = 'ehull')
# 
# dev.off()

#--BetaDisper
#--a multivariate analogue of Levene's test for homogeneity of variance
betadisper <- betadisper(comm.dist.morisita, group = morisita.matrix$Burn_status)
morisita.betadisper <- anova(betadisper)
morisita.betadisper

#--add Betadisper results to table
anosim.res[which(anosim.res$anosim.res =="morisita.p"), "F.betadisper"] <-
  morisita.betadisper$`F value`[1]
anosim.res[which(anosim.res$anosim.res =="morisita.p"), "df.betadisper.1"] <-
  morisita.betadisper$Df[1]
anosim.res[which(anosim.res$anosim.res =="morisita.p"), "df.betadisper.2"] <-
  morisita.betadisper$Df[2]
anosim.res[which(anosim.res$anosim.res =="morisita.p"), "p.betadisper"] <-
  morisita.betadisper$`Pr(>F)`[1]

#<< ANOSIM >>-----------------------------------------------------------------------------
morisita.anosim <- anosim(comm.matrix, grouping = morisita.matrix$Burn_status,
                         distance = "horn")
morisita.anosim

#--Add results to data frame
anosim.res[which(anosim.res$anosim.res =="morisita.p"), "r"] <- morisita.anosim$statistic
anosim.res[which(anosim.res$anosim.res =="morisita.p"), "p"] <- morisita.anosim$signif

#========================================================================================#
# Jaccard based dissimilarity index: Santa Catalina----
#========================================================================================#

#--isolate otu data
comm.matrix <- scm.matrix[9:length(scm.matrix)]

#--comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 2]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 1, ] # remove rows with sums of 0
jaccard.matrix <- scm.matrix[row.names(comm.matrix),]

#--distance matrix using jaccard index
comm.dist.jaccard <- vegdist(comm.matrix, method = "jaccard", binary = TRUE)

#--NMDS analysis
jaccard.otu <- metaMDS(comm.dist.jaccard, dist = "bray", permutations = 999,
                       try = 100, trymax = 1000)

#--Stress
jaccard.otu$stress

#--add stress of NMDS to results table
anosim.res[which(anosim.res$anosim.res =="jaccard.scm"), "stress.nmds"] <-
  jaccard.otu$stress

#--format and output NMDS plots to figure folder
# jpeg(filename = 'figures_output/NMDS_SC_Jaccard.jpeg',
#      width = 700, height = 600,
#      quality = 100)
# par(mfrow = c(1,1), "mar"=c(6, 5, 5, 3))
# 
# #--Plot NMDS of EM community based on Jaccard index and OTU abundance
# plot(jaccard.otu, display = "sites", type = "n", cex.lab = 2,
#      cex.axis = 1.5, xlab = 'Axis 1', ylab = 'Axis 2')
# 
# # colors for points
# color.vec <- data.frame(color = rep(NA, nrow(jaccard.matrix)),
#                         p.group = jaccard.matrix$Burn_status)
# color.vec[color.vec$p.group == 'burned', 'color'] <- 'black'
# color.vec[color.vec$p.group == 'unburned', 'color'] <- 'darkgreen'
# 
# color.vec[color.vec$p.group == 'burned', 'shape'] <- 16
# color.vec[color.vec$p.group == 'unburned', 'shape'] <- 16
# 
# #ordipointlabel(jaccard.otu, display = "sites")
# #--isolate points for pinaleno mts.
# burned <- row.names(jaccard.matrix[jaccard.matrix$Burn_status == 'burned',])
# unburned <- row.names(jaccard.matrix[!jaccard.matrix$Burn_status == 'burned',])
# 
# #--isolate points using rows isolated above
# points(jaccard.otu$points[burned,1:2], display = "sites", cex = 3,
#        pch = color.vec[color.vec$p.group == 'burned',
#                        'shape'],
#        col = color.vec[color.vec$p.group == 'burned',
#                        'color'],
#        bg = color.vec[color.vec$p.group == 'burned',
#                       'color'])
# points(jaccard.otu$points[unburned,1:2], display = "sites", cex = 3,
#        pch = color.vec[color.vec$p.group == 'unburned',
#                        'shape'],
#        col = color.vec[color.vec$p.group == 'unburned',
#                        'color'],
#        bg = color.vec[color.vec$p.group == 'unburned',
#                       'color'])
# # legend("topright", legend = c('FA','FU'), bty = "n",
# #        col = c('black','darkgreen'),
# #        pch = c(16,16),
# #        pt.bg =  c('black','darkgreen'), cex = 2)
# 
# #--Ordihull variations by range, burn, and both burn and range
# #ordihull(jaccard.otu, groups = color.vec$p.group)
# ordiellipse(jaccard.otu, groups = jaccard.matrix$Burn_status,
#             col = c('black','black'),
#             kind = 'ehull')
# #ordihull(jaccard.otu, groups = stsp.matrix$Range) # just mt. range
# #--Overlays
# # fit <- envfit(jaccard.otu ~ po4.p.ppm + temp.warmest.quarter +
# #                 pH.su, jaccard.matrix)
# # fit
# # plot(fit, type = 'n')
# 
# dev.off()

#--BetaDisper
#--a multivariate analogue of Levene's test for homogeneity of variance
betadisper <- betadisper(comm.dist.jaccard, group = jaccard.matrix$Burn_status)
jaccard.betadisper <- anova(betadisper)
jaccard.betadisper

#--add Betadisper results to table
anosim.res[which(anosim.res$anosim.res =="jaccard.scm"), "F.betadisper"] <-
  jaccard.betadisper$`F value`[1]
anosim.res[which(anosim.res$anosim.res =="jaccard.scm"), "df.betadisper.1"] <-
  jaccard.betadisper$Df[1]
anosim.res[which(anosim.res$anosim.res =="jaccard.scm"), "df.betadisper.2"] <-
  jaccard.betadisper$Df[2]
anosim.res[which(anosim.res$anosim.res =="jaccard.scm"), "p.betadisper"] <-
  jaccard.betadisper$`Pr(>F)`[1]

#<< ANOSIM >>-----------------------------------------------------------------------------
jaccard.anosim <- anosim(comm.matrix, grouping = jaccard.matrix$Burn_status,
                         distance = "jaccard")
jaccard.anosim

#--Add results to data frame
anosim.res[which(anosim.res$anosim.res =="jaccard.scm"), "r"] <- jaccard.anosim$statistic
anosim.res[which(anosim.res$anosim.res =="jaccard.scm"), "p"] <- jaccard.anosim$signif

#========================================================================================#
# Morisita based dissimilarity index: Santa Catalina----
#========================================================================================#

# #--comment to not remove outlier (NF19)
# stsp.matrix <- stsp.matrix[!stsp.matrix$Tree %in% c('NF19'),]

#--isolate site X species matrix only without metadata
comm.matrix <- scm.matrix[9:length(scm.matrix)]

#--comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 2]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 1, ] # remove rows with sums of 0
morisita.matrix <- scm.matrix[row.names(comm.matrix),]

#--distance matrix using jaccard index
comm.dist.morisita <- vegdist(comm.matrix, method = "horn", binary = F)

#--NMDS analysis
morisita.otu <- metaMDS(comm.dist.morisita, dist = "bray", permutations = 999,
                        try = 100, trymax = 1000)

#--add stress of NMDS to results table
anosim.res[which(anosim.res$anosim.res =="morisita.scm"), "stress.nmds"] <-
  morisita.otu$stress

#--format and output NMDS plots to figure folder
# jpeg(filename = 'figures_output/NMDS_SC_MorisitaHorn.jpeg',
#      width = 700, height = 600,
#      quality = 100)
# par(mfrow = c(1,1), "mar"=c(6, 5, 5, 3))
# 
# #--Plot NMDS of EM community based on Jaccard index and OTU abundance
# plot(morisita.otu, display = "sites", type = "n", cex.lab = 2,
#      cex.axis = 1.5, xlab = 'Axis 1', ylab = 'Axis 2')
# # colors for points
# color.vec <- data.frame (color = rep (NA, nrow(morisita.matrix)),
#                          p.group = morisita.matrix$Burn_status)
# color.vec[color.vec$p.group == 'burned', 'color'] <- 'black'
# color.vec[color.vec$p.group == 'unburned', 'color'] <- 'darkgreen'
# color.vec[color.vec$p.group == 'burned', 'shape'] <- 16
# color.vec[color.vec$p.group == 'unburned', 'shape'] <- 16
# 
# #ordipointlabel(morisita.otu, display = "sites")
# 
# #--isolate points for pinaleno mts.
# burned <- row.names(morisita.matrix[morisita.matrix$Burn_status == 'burned',])
# unburned <- row.names(morisita.matrix[!morisita.matrix$Burn_status == 'burned',])
# 
# #--isolate points using rows isolated above
# points(morisita.otu$points[burned,1:2], display = "sites", cex = 3,
#        pch = color.vec[color.vec$p.group == 'burned',
#                        'shape'],
#        col = color.vec[color.vec$p.group == 'burned',
#                        'color'],
#        bg = color.vec[color.vec$p.group == 'burned',
#                       'color'])
# points(morisita.otu$points[unburned,1:2], display = "sites", cex = 3,
#        pch = color.vec[color.vec$p.group == 'unburned',
#                        'shape'],
#        col = color.vec[color.vec$p.group == 'unburned',
#                        'color'],
#        bg = color.vec[color.vec$p.group == 'unburned',
#                       'color'])
# # legend("bottomright", legend = c('FA','FU'), bty = "n",
# #        col = c('black','darkgreen'),
# #        pch = c(16,16),
# #        pt.bg =  c('black','darkgreen'), cex = 2)
# 
# #--Ordihull variations by range, burn, and both burn and range
# #ordihull(morisita.otu, groups = color.vec$p.group,
# #         col = c('black','green4','black','green4')),
# ordiellipse(morisita.otu, groups = morisita.matrix$Burn_status,
#             col = c('black','black'),
#             kind = 'ehull') # just burn status
# #ordihull(morisita.otu, groups = stsp.matrix$Range,
# #         col = c('black','green4','black','green4')) # just mt. range
# 
# dev.off()

#--BetaDisper
#--a multivariate analogue of Levene's test for homogeneity of variance
betadisper <- betadisper(comm.dist.morisita, group = morisita.matrix$Burn_status)
morisita.betadisper <- anova(betadisper)
morisita.betadisper

#--add Betadisper results to table
anosim.res[which(anosim.res$anosim.res =="morisita.scm"), "F.betadisper"] <-
  morisita.betadisper$`F value`[1]
anosim.res[which(anosim.res$anosim.res =="morisita.scm"), "df.betadisper.1"] <-
  morisita.betadisper$Df[1]
anosim.res[which(anosim.res$anosim.res =="morisita.scm"), "df.betadisper.2"] <-
  morisita.betadisper$Df[2]
anosim.res[which(anosim.res$anosim.res =="morisita.scm"), "p.betadisper"] <-
  morisita.betadisper$`Pr(>F)`[1]

#<< PERMANOVA >>------------------------------------------------------------------------
morisita.adonis.anosim <- adonis(comm.dist.morisita ~ Burn_status,
                                data = morisita.matrix, permutations = 1000)
morisita.adonis.anosim

#--Burn f.model, r2, p-value
anosim.res[which(anosim.res$anosim.res == "morisita.scm"), "F.model.burn"] <-
  morisita.adonis.anosim$aov.tab$F.Model[1]
anosim.res[which(anosim.res$anosim.res == "morisita.scm"), "r2.burn"] <-
  morisita.adonis.anosim$aov.tab$R2[1]
anosim.res[which(anosim.res$anosim.res == "morisita.scm"), "p.burn"] <-
  morisita.adonis.anosim$aov.tab$`Pr(>F)`[1]


#========================================================================================#
# Jaccard based dissimilarity index: Fire unaffected----
#========================================================================================#

#--isolate site X species matrix only without metadata
comm.matrix <- FU.matrix[9:length(FU.matrix)]

#--comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 2]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 1, ] # remove rows with sums of 0
jaccard.matrix <- FU.matrix[row.names(comm.matrix),]

#--distance matrix using jaccard index
comm.dist.jaccard <- vegdist(comm.matrix, method = "jaccard", binary = T)

#--NMDS analysis
jaccard.otu <- metaMDS(comm.dist.jaccard, dist = "bray", permutations = 999,
                        try = 100, trymax = 1000)

#--add stress of NMDS to results table
anosim.res[which(anosim.res$anosim.res =="jaccard.fu"), "stress.nmds"] <-
  jaccard.otu$stress

# jpeg(filename = 'figures_output/NMDS_Jaccard_all_FU.jpeg',
#      width = 700, height = 600,
#      quality = 100)
# par(mfrow = c(1,1), "mar"=c(6, 5, 5, 3))
# 
# #--Plot NMDS of EM community based on Jaccard index and OTU abundance
# plot(jaccard.otu, display = "sites", type = "n", cex.lab = 2,
#      cex.axis = 1.5, xlab = 'Axis 1',ylab = 'Axis 2')
# 
# # colors for points
# color.vec <- data.frame(shape = rep(NA, nrow(jaccard.matrix)),
#                         p.group = jaccard.matrix$Range)
# color.vec[color.vec$p.group == 'santa.catalina', 'shape'] <- 16
# color.vec[color.vec$p.group == 'pinaleno', 'shape'] <- 15
# 
# #ordipointlabel(morisita.otu, display = "sites")
# #--isolate points for pinaleno mts.
# pinaleno.unburn <- row.names(jaccard.matrix[jaccard.matrix$Range == 'pinaleno' & 
#                                               jaccard.matrix$Burn_status == 'unburned',])
# santa.catalina.unburn <- row.names(jaccard.matrix[jaccard.matrix$Range == 'santa.catalina' & 
#                                                     jaccard.matrix$Burn_status == 'unburned',])
# 
# #--isolate points using rows isolated above
# points(jaccard.otu$points[pinaleno.unburn,1:2], display = "sites", cex = 3,
#        pch = color.vec[color.vec$p.group == 'pinaleno' | 
#                          color.vec$p.group == 'pinaleno',
#                        'shape'],
#        col = 'darkgreen',
#        bg = 'darkgreen')
# points(jaccard.otu$points[santa.catalina.unburn,1:2], display = "sites", cex = 3,
#        pch = color.vec[color.vec$p.group == 'santa.catalina' | 
#                          color.vec$p.group == 'santa.catalina',
#                        'shape'],
#        col = 'darkgreen',
#        bg = 'darkgreen')
# # legend("topright", legend = c('Pinaleno Mts., burned site',
# #                               'Pinaleno Mts., unburned site',
# #                               'Santa Catalina Mts., burned site',
# #                               'Santa Catalina Mts., unburned site'), bty = "n",
# #         col = c('black','green4','black','green4'),
# #         pch = c(6,6,20,20),
# #         pt.bg =  c('black','green4','black','green4'), cex = 1.5)
# 
# # legend("bottomleft", legend = c('Pinaleno Mts.',
# #                                 'Santa Catalina Mts.'), bty = "n",
# #        col = c('black', 'black'),
# #        pch = c(15,16),
# #        pt.bg =  c('black','black'), cex = 2)
# # 
# # legend("bottomright", legend = c('FA',
# #                                  'FU'), bty = "n",
# #        col = c('black', 'darkgreen'),
# #        pch = c(16,16),
# #        pt.bg =  c('black', 'darkgreen'), cex = 2)
# 
# #--Ordihull variations by range, burn, and both burn and range
# for(i in 1:nrow(jaccard.matrix)){
#   row.i <- paste(jaccard.matrix[i, 'Range'],
#                  jaccard.matrix[i, 'Burn_status'])
#   jaccard.matrix[i,'range_burn'] <- row.i
# }
# ordiellipse(jaccard.otu,
#             groups = jaccard.matrix$range_burn,
#             col = c('black','black'),
#             kind = 'ehull')
# 
# dev.off()

#--BetaDisper
#--a multivariate analogue of Levene's test for homogeneity of variance
betadisper <- betadisper(comm.dist.jaccard, group = FU.matrix$Range)
jaccard.betadisper <- anova(betadisper)
jaccard.betadisper

#--add Betadisper results to table
anosim.res[which(anosim.res$anosim.res =="jaccard.fu"), "F.betadisper"] <-
  jaccard.betadisper$`F value`[1]
anosim.res[which(anosim.res$anosim.res =="jaccard.fu"), "df.betadisper.1"] <-
  jaccard.betadisper$Df[1]
anosim.res[which(anosim.res$anosim.res =="jaccard.fu"), "df.betadisper.2"] <-
  jaccard.betadisper$Df[2]
anosim.res[which(anosim.res$anosim.res =="jaccard.fu"), "p.betadisper"] <-
  jaccard.betadisper$`Pr(>F)`[1]

#<< PERMANOVA >>------------------------------------------------------------------------
jaccard.adonis.anosim <- adonis(comm.dist.jaccard ~ Range,
                                 data = jaccard.matrix, permutations = 1000)
jaccard.adonis.anosim

#--Burn f.model, r2, p-value
anosim.res[which(anosim.res$anosim.res == "jaccard.fu"), "F.model.burn"] <-
  jaccard.adonis.anosim$aov.tab$F.Model[1]
anosim.res[which(anosim.res$anosim.res == "jaccard.fu"), "r2.burn"] <-
  jaccard.adonis.anosim$aov.tab$R2[1]
anosim.res[which(anosim.res$anosim.res == "jaccard.fu"), "p.burn"] <-
  jaccard.adonis.anosim$aov.tab$`Pr(>F)`[1]

#========================================================================================#
# Morisita based dissimilarity index: Fire unaffected----
#========================================================================================#

#--isolate site X species matrix only without metadata
comm.matrix <- FU.matrix[9:length(FU.matrix)]

#--comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 2]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 1, ] # remove rows with sums of 0
morisita.matrix <- FU.matrix[row.names(comm.matrix),]

#--distance matrix using jaccard index
comm.dist.morisita <- vegdist(comm.matrix, method = "horn", binary = F)

#--NMDS analysis
morisita.otu <- metaMDS(comm.dist.morisita, dist = "bray", permutations = 999,
                        try = 100, trymax = 1000)

#--add stress of NMDS to results table
anosim.res[which(anosim.res$anosim.res =="morisita.fu"), "stress.nmds"] <-
  morisita.otu$stress

jpeg(filename = 'figures_output/NMDS_MorisitaHorn_all_FU.jpeg',
     width = 700, height = 600,
     quality = 100)
par(mfrow = c(1,1), "mar"=c(6, 5, 5, 3))

#--Plot NMDS of EM community based on Jaccard index and OTU abundance
plot(morisita.otu, display = "sites", type = "n", cex.lab = 2,
     cex.axis = 1.5, xlab = 'Axis 1',ylab = 'Axis 2',
     ylim = c(-0.5, 0.5))

# colors for points
color.vec <- data.frame(shape = rep(NA, nrow(morisita.matrix)),
                        p.group = morisita.matrix$Range)
color.vec[color.vec$p.group == 'santa.catalina', 'shape'] <- 16
color.vec[color.vec$p.group == 'pinaleno', 'shape'] <- 15

#ordipointlabel(morisita.otu, display = "sites")
#--isolate points for pinaleno mts.
pinaleno.unburn <- row.names(morisita.matrix[morisita.matrix$Range == 'pinaleno' &
                                               morisita.matrix$Burn_status == 'unburned',])
santa.catalina.unburn <- row.names(morisita.matrix[morisita.matrix$Range == 'santa.catalina' &
                                                     morisita.matrix$Burn_status == 'unburned',])

#--isolate points using rows isolated above
points(morisita.otu$points[pinaleno.unburn,1:2], display = "sites", cex = 3,
       pch = color.vec[color.vec$p.group == 'pinaleno' |
                         color.vec$p.group == 'pinaleno',
                       'shape'],
       col = 'darkgreen',
       bg = 'darkgreen')
points(morisita.otu$points[santa.catalina.unburn,1:2], display = "sites", cex = 3,
       pch = color.vec[color.vec$p.group == 'santa.catalina' |
                         color.vec$p.group == 'santa.catalina',
                       'shape'],
       col = 'darkgreen',
       bg = 'darkgreen')
# legend("topright", legend = c('Pinaleno Mts., burned site',
#                               'Pinaleno Mts., unburned site',
#                               'Santa Catalina Mts., burned site',
#                               'Santa Catalina Mts., unburned site'), bty = "n",
#         col = c('black','green4','black','green4'),
#         pch = c(6,6,20,20),
#         pt.bg =  c('black','green4','black','green4'), cex = 1.5)

# legend("bottomleft", legend = c('Pinaleno Mts.',
#                                 'Santa Catalina Mts.'), bty = "n",
#        col = c('black', 'black'),
#        pch = c(15,16),
#        pt.bg =  c('black','black'), cex = 2)
#
# legend("bottomright", legend = c('FA',
#                                  'FU'), bty = "n",
#        col = c('black', 'darkgreen'),
#        pch = c(16,16),
#        pt.bg =  c('black', 'darkgreen'), cex = 2)

#--Ordihull variations by range, burn, and both burn and range
for(i in 1:nrow(morisita.matrix)){
  row.i <- paste(morisita.matrix[i, 'Range'],
                 morisita.matrix[i, 'Burn_status'])
  morisita.matrix[i,'range_burn'] <- row.i
}
ordiellipse(morisita.otu,
            groups = morisita.matrix$range_burn,
            col = c('black','black'),
            kind = 'ehull')

dev.off()


#--BetaDisper
#--a multivariate analogue of Levene's test for homogeneity of variance
betadisper <- betadisper(comm.dist.morisita, group = FU.matrix$Range)
morisita.betadisper <- anova(betadisper)
morisita.betadisper

#--add Betadisper results to table
anosim.res[which(anosim.res$anosim.res =="morisita.fu"), "F.betadisper"] <-
  morisita.betadisper$`F value`[1]
anosim.res[which(anosim.res$anosim.res =="morisita.fu"), "df.betadisper.1"] <-
  morisita.betadisper$Df[1]
anosim.res[which(anosim.res$anosim.res =="morisita.fu"), "df.betadisper.2"] <-
  morisita.betadisper$Df[2]
anosim.res[which(anosim.res$anosim.res =="morisita.fu"), "p.betadisper"] <-
  morisita.betadisper$`Pr(>F)`[1]

#<< PERMANOVA >>------------------------------------------------------------------------
morisita.adonis.anosim <- adonis(comm.dist.morisita ~ Range,
                                data = morisita.matrix, permutations = 1000)
morisita.adonis.anosim

#--Burn f.model, r2, p-value
anosim.res[which(anosim.res$anosim.res == "morisita.fu"), "F.model.burn"] <-
  morisita.adonis.anosim$aov.tab$F.Model[1]
anosim.res[which(anosim.res$anosim.res == "morisita.fu"), "r2.burn"] <-
  morisita.adonis.anosim$aov.tab$R2[1]
anosim.res[which(anosim.res$anosim.res == "morisita.fu"), "p.burn"] <-
  morisita.adonis.anosim$aov.tab$`Pr(>F)`[1]


#========================================================================================#
# Jaccard based dissimilarity index: Fire affected----
#========================================================================================#

#--isolate site X species matrix only without metadata
comm.matrix <- FA.matrix[9:length(FA.matrix)]

#--comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 2]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 1, ] # remove rows with sums of 0
jaccard.matrix <- FA.matrix[row.names(comm.matrix),]

#--distance matrix using jaccard index
comm.dist.jaccard <- vegdist(comm.matrix, method = "jaccard", binary = T)

#--NMDS analysis
jaccard.otu <- metaMDS(comm.dist.jaccard, dist = "bray", permutations = 999,
                        try = 100, trymax = 1000)

#--add stress of NMDS to results table
anosim.res[which(anosim.res$anosim.res =="jaccard.fa"), "stress.nmds"] <-
  jaccard.otu$stress

# jpeg(filename = 'figures_output/NMDS_Jaccard_all_FA.jpeg',
#      width = 700, height = 600,
#      quality = 100)
# par(mfrow = c(1,1), "mar"=c(6, 5, 5, 3))
# 
# #--Plot NMDS of EM community based on Jaccard index and OTU abundance
# plot(jaccard.otu, display = "sites", type = "n", cex.lab = 2,
#      cex.axis = 1.5, xlab = 'Axis 1',ylab = 'Axis 2')
# 
# # colors for points
# color.vec <- data.frame(shape = rep(NA, nrow(jaccard.matrix)),
#                         p.group = jaccard.matrix$Range)
# color.vec[color.vec$p.group == 'santa.catalina', 'shape'] <- 16
# color.vec[color.vec$p.group == 'pinaleno', 'shape'] <- 15
# 
# #ordipointlabel(jaccard.otu, display = "sites")
# #--isolate points for pinaleno mts.
# pinaleno.burn <- row.names(jaccard.matrix[jaccard.matrix$Range == 'pinaleno' & 
#                                              jaccard.matrix$Burn_status == 'burned',])
# santa.catalina.burn <- row.names(jaccard.matrix[jaccard.matrix$Range == 'santa.catalina' & 
#                                                   jaccard.matrix$Burn_status == 'burned',])
# 
# #--isolate points using rows isolated above
# points(jaccard.otu$points[pinaleno.burn,1:2], display = "sites", cex = 3,
#        pch = color.vec[color.vec$p.group == 'pinaleno' | 
#                          color.vec$p.group == 'pinaleno',
#                        'shape'],
#        col = 'black',
#        bg = 'black')
# 
# points(jaccard.otu$points[santa.catalina.burn,1:2], display = "sites", cex = 3,
#        pch = color.vec[color.vec$p.group == 'santa.catalina' | 
#                          color.vec$p.group == 'santa.catalina',
#                        'shape'],
#        col = 'black',
#        bg = 'black')
# 
# # legend("topright", legend = c('Pinaleno Mts., burned site',
# #                               'Pinaleno Mts., unburned site',
# #                               'Santa Catalina Mts., burned site',
# #                               'Santa Catalina Mts., unburned site'), bty = "n",
# #         col = c('black','green4','black','green4'),
# #         pch = c(6,6,20,20),
# #         pt.bg =  c('black','green4','black','green4'), cex = 1.5)
# 
# # legend("bottomleft", legend = c('Pinaleno Mts.',
# #                                 'Santa Catalina Mts.'), bty = "n",
# #        col = c('black', 'black'),
# #        pch = c(15,16),
# #        pt.bg =  c('black','black'), cex = 2)
# # 
# # legend("bottomright", legend = c('FA',
# #                                  'FU'), bty = "n",
# #        col = c('black', 'darkgreen'),
# #        pch = c(16,16),
# #        pt.bg =  c('black', 'darkgreen'), cex = 2)
# 
# #--Ordihull variations by range, burn, and both burn and range
# for(i in 1:nrow(jaccard.matrix)){
#   row.i <- paste(jaccard.matrix[i, 'Range'],
#                  jaccard.matrix[i, 'Burn_status'])
#   jaccard.matrix[i,'range_burn'] <- row.i
# }
# ordiellipse(jaccard.otu,
#             groups = jaccard.matrix$range_burn,
#             col = c('black','black'),
#             kind = 'ehull')
# 
# dev.off()

#--BetaDisper
#--a multivariate analogue of Levene's test for homogeneity of variance
betadisper <- betadisper(comm.dist.jaccard, group = FA.matrix$Range)
jaccard.betadisper <- anova(betadisper)
jaccard.betadisper

#--add Betadisper results to table
anosim.res[which(anosim.res$anosim.res =="jaccard.fa"), "F.betadisper"] <-
  jaccard.betadisper$`F value`[1]
anosim.res[which(anosim.res$anosim.res =="jaccard.fa"), "df.betadisper.1"] <-
  jaccard.betadisper$Df[1]
anosim.res[which(anosim.res$anosim.res =="jaccard.fa"), "df.betadisper.2"] <-
  jaccard.betadisper$Df[2]
anosim.res[which(anosim.res$anosim.res =="jaccard.fa"), "p.betadisper"] <-
  jaccard.betadisper$`Pr(>F)`[1]

#<< PERMANOVA >>------------------------------------------------------------------------
jaccard.adonis.anosim <- adonis(comm.dist.jaccard ~ Range,
                                 data = jaccard.matrix, permutations = 1000)
jaccard.adonis.anosim

#--Burn f.model, r2, p-value
anosim.res[which(anosim.res$anosim.res == "jaccard.fa"), "F.model.burn"] <-
  jaccard.adonis.anosim$aov.tab$F.Model[1]
anosim.res[which(anosim.res$anosim.res == "jaccard.fa"), "r2.burn"] <-
  jaccard.adonis.anosim$aov.tab$R2[1]
anosim.res[which(anosim.res$anosim.res == "jaccard.fa"), "p.burn"] <-
  jaccard.adonis.anosim$aov.tab$`Pr(>F)`[1]


#========================================================================================#
# Morisita based dissimilarity index: Fire affected----
#========================================================================================#

#--isolate site X species matrix only without metadata
comm.matrix <- FA.matrix[9:length(FA.matrix)]

#--comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 2]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 1, ] # remove rows with sums of 0
morisita.matrix <- FA.matrix[row.names(comm.matrix),]

#--distance matrix using jaccard index
comm.dist.morisita <- vegdist(comm.matrix, method = "horn", binary = F)

#--NMDS analysis
morisita.otu <- metaMDS(comm.dist.morisita, dist = "bray", permutations = 999,
                        try = 100, trymax = 1000)

#--add stress of NMDS to results table
anosim.res[which(anosim.res$anosim.res =="morisita.fa"), "stress.nmds"] <-
  morisita.otu$stress


#<< Base r plot >> ----
jpeg(filename = 'figures_output/NMDS_MorisitaHorn_all_FA.jpeg',
     width = 700, height = 600,
     quality = 100)
par(mfrow = c(1,1), "mar"=c(6, 5, 5, 3))

#--Plot NMDS of EM community based on Jaccard index and OTU abundance
plot(morisita.otu, display = "sites", type = "n", cex.lab = 2,
     cex.axis = 1.5, xlab = 'Axis 1',ylab = 'Axis 2',
     ylim = c(-0.5, 0.5))

# colors for points
color.vec <- data.frame(shape = rep(NA, nrow(morisita.matrix)),
                        p.group = morisita.matrix$Range)
color.vec[color.vec$p.group == 'santa.catalina', 'shape'] <- 16
color.vec[color.vec$p.group == 'pinaleno', 'shape'] <- 15

#ordipointlabel(jaccard.otu, display = "sites")
#--isolate points for pinaleno mts.
pinaleno.burn <- row.names(morisita.matrix[morisita.matrix$Range == 'pinaleno' &
                                             morisita.matrix$Burn_status == 'burned',])
santa.catalina.burn <- row.names(morisita.matrix[morisita.matrix$Range == 'santa.catalina' &
                                                   morisita.matrix$Burn_status == 'burned',])

#--isolate points using rows isolated above
points(morisita.otu$points[pinaleno.burn,1:2], display = "sites", cex = 3,
       pch = color.vec[color.vec$p.group == 'pinaleno' |
                         color.vec$p.group == 'pinaleno',
                       'shape'],
       col = 'black',
       bg = 'black')

points(morisita.otu$points[santa.catalina.burn,1:2], display = "sites", cex = 3,
       pch = color.vec[color.vec$p.group == 'santa.catalina' |
                         color.vec$p.group == 'santa.catalina',
                       'shape'],
       col = 'black',
       bg = 'black')

# legend("topright", legend = c('Pinaleno Mts., burned site',
#                               'Pinaleno Mts., unburned site',
#                               'Santa Catalina Mts., burned site',
#                               'Santa Catalina Mts., unburned site'), bty = "n",
#         col = c('black','green4','black','green4'),
#         pch = c(6,6,20,20),
#         pt.bg =  c('black','green4','black','green4'), cex = 1.5)

# legend("bottomleft", legend = c('Pinaleno Mts.',
#                                 'Santa Catalina Mts.'), bty = "n",
#        col = c('black', 'black'),
#        pch = c(15,16),
#        pt.bg =  c('black','black'), cex = 2)
#
# legend("bottomright", legend = c('FA',
#                                  'FU'), bty = "n",
#        col = c('black', 'darkgreen'),
#        pch = c(16,16),
#        pt.bg =  c('black', 'darkgreen'), cex = 2)

#--Ordihull variations by range, burn, and both burn and range
for(i in 1:nrow(morisita.matrix)){
  row.i <- paste(morisita.matrix[i, 'Range'],
                 morisita.matrix[i, 'Burn_status'])
  morisita.matrix[i,'range_burn'] <- row.i
}
ordiellipse(morisita.otu,
            groups = morisita.matrix$range_burn,
            col = c('black','black'),
            kind = 'ehull')

dev.off()

#<< GGplot plot: Fig. 3B and inset (taxonomic data) >> ----
#--Create dataframe
data.scores <- as.data.frame(scores(morisita.otu))
data.scores$site <- morisita.matrix$Site
data.scores$tree <- morisita.matrix$Tree
data.scores$Range <- morisita.matrix$Range
data.scores[data.scores$Range == 'pinaleno','Range'] <- 'Pinaleno Mts.'
data.scores[data.scores$Range == 'santa.catalina','Range'] <- 'Santa Catalina Mts.'

fa <- ggplot() + 
  geom_point(data = data.scores,aes(x = NMDS1,
                                 y = NMDS2,
                                 shape=Range),size=4) + # add the point markers
  coord_equal() +
  theme_bw() + 
  scale_color_manual(values='black') +
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=28,margin = margin(t = 30)), # remove x-axis labels
        axis.title.y = element_text(size=28,margin = margin(r = 30)), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())

#--Taxonomic data --> inset
#Isolated burned data
burn.tax <- tax.data[tax.data$Burn_status == 'burned',]
burn.tax <- burn.tax[!burn.tax$Taxonomy_class %in% c('Saccharomycetes',
                                                     'Eurotiomycetes',NA),]
#Change levels of Burn_status and Range columns
burn.tax[burn.tax$Range == 'santa.catalina', 'Range'] <- 'Santa Catalina Mts.'
burn.tax[burn.tax$Range == 'pinaleno', 'Range'] <- 'Pinaleno Mts.'
burn.tax$Range <-  factor(burn.tax$Range, 
                          levels = c('Santa Catalina Mts.', 'Pinaleno Mts.'))

#Bar graph
burn.class <- ggplot(data = burn.tax, 
                     aes(x = Range,
                         fill = Taxonomy_class)) + 
  geom_bar(position = "fill") + 
  ylab("Proportion of \n sequences per class") +
  theme_bw() +
  xlab(element_blank()) +
  #ggtitle("Proportion of Classes by Topography") +
  scale_fill_brewer(palette = 6,
                    direction = -1) +
  guides (fill=guide_legend(title=NULL)) +
  theme(legend.position='right',
        # axis.title.x = element_text(margin = margin(t = 30)),
        # axis.title.y = element_text(margin = margin(r = 30)),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        axis.text = element_text(size=20, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

vp <- viewport(width = 0.4, height = 0.4, x = 1,
               y = unit(0.7, "lines"), just = c("right",
                                                "bottom"))

full <- function(){
  print(fa)
  print(burn.class, vp = vp)
}
full()

#--BetaDisper
#--a multivariate analogue of Levene's test for homogeneity of variance
betadisper <- betadisper(comm.dist.morisita, group = FA.matrix$Range)
morisita.betadisper <- anova(betadisper)
morisita.betadisper

#--add Betadisper results to table
anosim.res[which(anosim.res$anosim.res =="morisita.fa"), "F.betadisper"] <-
  morisita.betadisper$`F value`[1]
anosim.res[which(anosim.res$anosim.res =="morisita.fa"), "df.betadisper.1"] <-
  morisita.betadisper$Df[1]
anosim.res[which(anosim.res$anosim.res =="morisita.fa"), "df.betadisper.2"] <-
  morisita.betadisper$Df[2]
anosim.res[which(anosim.res$anosim.res =="morisita.fa"), "p.betadisper"] <-
  morisita.betadisper$`Pr(>F)`[1]

#<< ANOSIM >>-----------------------------------------------------------------------------
morisita.anosim <- anosim(comm.matrix, grouping = FA.matrix$Range,
                          distance = "horn")
morisita.anosim

#--Add results to data frame
anosim.res[which(anosim.res$anosim.res =="morisita.fa"), "r"] <- morisita.anosim$statistic
anosim.res[which(anosim.res$anosim.res =="morisita.fa"), "p"] <- morisita.anosim$signif

#--export anosim table
write.csv(anosim.res, paste0(res.dir, "ANOSIM_SequenceBased.csv"),
          row.names = F)

# #========================================================================================#
# # Venn diagram of PERMANOVA restuls----
# #========================================================================================#
# 
# #install.packages('VennDiagram')
# library(VennDiagram)
# 
# #<< Jaccard >>-----------------------------------------------------------------------------
# #fire history 0.04410 r-squared
# #range 0.07065 r-squared
# #fire history * range 0.09072 r-squared
# 
# jpeg(filename = paste0('/Volumes/Cenococcum/PhD/Dissertation/Chpt.2/Figures/',
#                        'Fig3d_jaccard.jpeg'), width = 700, height = 600,
#      quality = 100)
# 
# grid.newpage()
# draw.pairwise.venn(area1 = 13.482,
#                    area2 = 16.137,
#                    cross.area = 9.072,
#                    category = c("Fire history","Range"),
#                    lty = rep("blank", 2), col = 'black',fill = c("grey", "light blue"), 
#                    alpha = rep(0.5, 2), cat.pos = c(200, 160), cat.dist = 0.05,
#                    cex = 1.5, cat.cex = 2)
# 
# dev.off()
# 
# #<< Morisita-horn >>--------------------------------------------------------------------
# #fire history 0.06127 r-squared
# #range 0.06905 r-squared
# #fire history * range 0.08595 r-squared
# 
# jpeg(filename = paste0('/Volumes/Cenococcum/PhD/Dissertation/Chpt.2/Figures/',
#                        'Fig3d_morisitahorn.jpeg'), width = 700, height = 600,
#      quality = 100)
# 
# grid.newpage()
# draw.pairwise.venn(area1 = 14.722,
#                    area2 = 15.500,
#                    cross.area = 8.595,
#                    category = c("Fire history","Range"),
#                    lty = rep("blank", 2), col = 'black',fill = c("grey", "light blue"), 
#                    alpha = rep(0.5, 2), cat.pos = c(200, 160), cat.dist = 0.05,
#                    cex = 1.5, cat.cex = 2)
# 
# dev.off()