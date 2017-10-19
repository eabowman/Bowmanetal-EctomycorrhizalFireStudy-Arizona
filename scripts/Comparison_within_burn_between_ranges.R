## Script created by Liz Bowman June 3, 2017
## for analyzing whether burned/unburned sites are similar between ranges
## if they are then they can be used as a burned/unburned group (2 groups instead of 4)

#=========================================================================================
# Load data and libraries----
#=========================================================================================

#-----------------------------------------------------------------------------------------
# Load libraries----
#-----------------------------------------------------------------------------------------
library(ggplot2);library (tidyr);library (vegan);library (dplyr)

#-----------------------------------------------------------------------------------------
# set up paths to directories----
#-----------------------------------------------------------------------------------------
#--path to directory where the climate data downloaded from BIOCLIM site is stored
dat.dir <- "~/Documents/PhD/3_EM_Fire_effect/data/"
fig.dir <- '~/Documents/PhD/3_EM_Fire_effect/figures/'
res.dir <- "~/Documents/PhD/3_EM_Fire_effect/results/"

#-----------------------------------------------------------------------------------------
# Load data and clean up----
#-----------------------------------------------------------------------------------------
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
  select (Sample_name,Burn_status,Range,Site,Tree,otu.97,otu.count) %>%
  spread (otu.97,otu.count) %>%
  group_by(Tree,Site,Range,Burn_status) %>%
  select (matches ("otu.*")) %>%
  summarize_each (funs (sum (., na.rm = TRUE))) %>%
  as.data.frame () -> stsp.matrix

#--climate data: add to stsp.matrix
clim.data <- read.csv(paste0(dat.dir, 'climate_data.csv'))
for(i in unique(stsp.matrix$Site)) {
  stsp.matrix[stsp.matrix$Site == i, 'prec'] <- clim.data[clim.data$site == i, 'prec']
  stsp.matrix[stsp.matrix$Site == i, 'temp.warmest.quarter'] <-
    clim.data[clim.data$site == i, 'BIO10_red']
}

#<< Per site mean of significant soil data (pH, PO4, SO4, and B) >> ----------------------
soil.data <- read.csv(paste0(dat.dir, 'soil_data.csv'),as.is = T)
soil.sig <- c('site','pH.su','po4.p.ppm','so4.s.ppm','b.ppm')
soil.data <- soil.data[colnames(soil.data)%in% soil.sig]
for(s in unique(soil.data$site)){
  for(f in soil.sig[2:5]){
    stsp.matrix[stsp.matrix$Site == s, f] <- mean(soil.data[soil.data$site == s, f])
  }
}

#--reorder matrix
stsp.matrix <- stsp.matrix[c(1:4,121:126,5:120)]

#<< separate data tables into burned and unburned tables >> ----------------------
burned.matrix <- stsp.matrix[stsp.matrix$Burn_status == 'burned',]
unburned.matrix <- stsp.matrix[!stsp.matrix$Burn_status == 'burned',]

#-----------------------------------------------------------------------------------------
# Create result tables
#-----------------------------------------------------------------------------------------
#--Table of Anosim results
anosim.res <- c("jaccard", "morisita")
anosim.res <- data.frame(anosim.res)
#--Table of PERMANOVA results
permanova.res <- data.frame(c('jaccard','morisita'))
colnames(permanova.res) <- "test"

#=========================================================================================
# Jaccard based dissimilarity index: Burned sites----
#=========================================================================================

#--isolate otu data
burned.comm.matrix <- burned.matrix[11:length(burned.matrix)]

#--comment to include singletons
burned.comm.matrix <- burned.comm.matrix[colSums(burned.comm.matrix) > 1]
burned.comm.matrix <- burned.comm.matrix[rowSums(burned.comm.matrix) > 0, ] # remove rows with sums of 0
burned.matrix <- burned.matrix[row.names(burned.comm.matrix),]

#--distance matrix using jaccard index
burn.dist.jaccard <- vegdist(burned.comm.matrix, method = "jaccard", binary = TRUE)

#--NMDS analysis
jaccard.otu <- metaMDS(burn.dist.jaccard, dist = "bray", permutations = 999,
                       try = 100, trymax = 1000)

#--Stress
jaccard.otu$stress

#--add stress of NMDS to results table
anosim.res[which(anosim.res$anosim.res =="jaccard"), "stress.nmds"] <-
  jaccard.otu$stress

#--format and output NMDS plots to figure folder
jpeg(filename = paste0(fig.dir, 'burned_NDMS.jpeg'), width = 1200, height = 450,
     quality = 100)
par(mfrow = c(1,2), "mar"=c(6, 5, 5, 3))

#--Plot NMDS of EM community based on Jaccard index and OTU abundance
plot(jaccard.otu, display = "sites", type = "n", cex.lab = 1.5,
     cex.axis = 1.5, yaxt = "n")
axis (2, at = seq (-0.4, 0.4, by = 0.2), cex.axis = 1.5, las = 2)
# colors for points
color.vec <- data.frame(color = rep (NA, nrow(burned.matrix)),
                         p.group = burned.matrix$Range)
color.vec[color.vec$p.group == 'santa.catalina', 'color'] <- 'grey'
color.vec[color.vec$p.group == 'pinaleno', 'color'] <- 'black'

#ordipointlabel(jaccard.otu, display = "sites")

#--isolate points using rows isolated above
points(jaccard.otu, display = "sites", cex = 1.5,
       col = color.vec$color,
       bg = color.vec$color)

legend("bottomleft", legend = levels(color.vec$p.group), bty = "n",
       col = c('black','grey'),
       pch = c(6,6,20,20),
       pt.bg =  c('black','grey'), cex = 1)
#--Ordihull variations by range, burn, and both burn and range
#ordihull(jaccard.otu, groups = color.vec$p.group)
ordihull(jaccard.otu, groups = color.vec$p.group,
         col = c('black','grey'))
#ordihull(jaccard.otu, groups = stsp.matrix$Range) # just mt. range

#--BetaDisper
#--a multivariate analogue of Levene's test for homogeneity of variance
betadisper <- betadisper(burn.dist.jaccard, group = color.vec$p.group)
jaccard.betadisper <- anova(betadisper)

#--add Betadisper results to table
anosim.res[which(anosim.res$anosim.res =="jaccard"), "F.betadisper"] <-
  jaccard.betadisper$`F value`[1]
anosim.res[which(anosim.res$anosim.res =="jaccard"), "df.betadisper.1"] <-
  jaccard.betadisper$Df[1]
anosim.res[which(anosim.res$anosim.res =="jaccard"), "df.betadisper.2"] <-
  jaccard.betadisper$Df[2]
anosim.res[which(anosim.res$anosim.res =="jaccard"), "p.betadisper"] <-
  jaccard.betadisper$`Pr(>F)`[1]

#<< ANOSIM >>-----------------------------------------------------------------------------
jaccard.anosim <- anosim(burned.comm.matrix, grouping = color.vec$p.group,
                         distance = "jaccard")
jaccard.anosim

#--Add results to data frame
anosim.res[which(anosim.res$anosim.res =="jaccard"), "r"] <- jaccard.anosim$statistic
anosim.res[which(anosim.res$anosim.res =="jaccard"), "p"] <- jaccard.anosim$signif

#<< PERMANOVA >>--------------------------------------------------------------------------
#--adonis
jaccard.adonis <- adonis(formula = burn.dist.jaccard ~ Range * po4.p.ppm *
                         temp.warmest.quarter, data = burned.matrix, permutations = 1000)
jaccard.adonis

#--add results to data.frame
#--Range f.model, r2, p-value
permanova.res[which(permanova.res$test == "jaccard"), "F.model.range"] <-
 jaccard.adonis$aov.tab$F.Model[2]
permanova.res[which(permanova.res$test == "jaccard"), "r2.range"] <-
 jaccard.adonis$aov.tab$R2[2]
permanova.res[which(permanova.res$test == "jaccard"), "p.range"] <-
 jaccard.adonis$aov.tab$`Pr(>F)`[2]

#--Phosphate f.model, r2, p-value
permanova.res[which(permanova.res$test == "jaccard"), "F.model.po4"] <- 
  jaccard.adonis$aov.tab$F.Model[2]
permanova.res[which(permanova.res$test == "jaccard"), "r2.po4"] <- 
  jaccard.adonis$aov.tab$R2[2]
permanova.res[which(permanova.res$test == "jaccard"), "p.po4"] <- 
  jaccard.adonis$aov.tab$`Pr(>F)`[2]

#--prec f.model, r2, p-value
permanova.res[which(permanova.res$test == "jaccard"), "F.model.temp"] <- 
  jaccard.adonis$aov.tab$F.Model[3]
permanova.res[which(permanova.res$test == "jaccard"), "r2.temp"] <- 
  jaccard.adonis$aov.tab$R2[3]
permanova.res[which(permanova.res$test == "jaccard"), "p.temp"] <- 
  jaccard.adonis$aov.tab$`Pr(>F)`[3]

#=========================================================================================
# Morisita based dissimilarity index: Burned sites----
#=========================================================================================
#--isolate otu data
burned.comm.matrix <- burned.matrix[11:length(burned.matrix)]

#--comment to include singletons
burned.comm.matrix <- burned.comm.matrix[colSums(burned.comm.matrix) > 1]
burned.comm.matrix <- burned.comm.matrix[rowSums(burned.comm.matrix) > 0, ] # remove rows with sums of 0
burned.matrix <- burned.matrix[row.names(burned.comm.matrix),]

#--distance matrix using jaccard index
burn.dist.morisita <- vegdist(burned.comm.matrix, method = "morisita")

#--NMDS analysis
morisita.otu <- metaMDS(burn.dist.morisita, dist = "bray", permutations = 999,
                       try = 100, trymax = 1000)

#--Stress
morisita.otu$stress

#--add stress of NMDS to results table
anosim.res[which(anosim.res$anosim.res =="morisita"), "stress.nmds"] <-
  morisita.otu$stress

#--Plot NMDS of EM community based on Jaccard index and OTU abundance
plot(morisita.otu, display = "sites", type = "n", cex.lab = 1.5,
     cex.axis = 1.5, yaxt = "n")
axis (2, at = seq (-0.4, 0.4, by = 0.2), cex.axis = 1.5, las = 2)
# colors for points
color.vec <- data.frame(color = rep (NA, nrow(burned.matrix)),
                        p.group = burned.matrix$Range)
color.vec[color.vec$p.group == 'santa.catalina', 'color'] <- 'grey'
color.vec[color.vec$p.group == 'pinaleno', 'color'] <- 'black'

#ordipointlabel(jaccard.otu, display = "sites")

#--isolate points using rows isolated above
points(morisita.otu, display = "sites", cex = 1.5,
       col = color.vec$color,
       bg = color.vec$color)

legend("bottomleft", legend = levels(color.vec$p.group), bty = "n",
       col = c('black','grey'),
       pch = c(6,6,20,20),
       pt.bg =  c('black','grey'), cex = 1)
#--Ordihull variations by range, burn, and both burn and range
#ordihull(morisita.otu, groups = color.vec$p.group)
ordihull(morisita.otu, groups = color.vec$p.group,
         col = c('black','grey'))
#ordihull(morisita.otu, groups = stsp.matrix$Range) # just mt. range

#--BetaDisper
#--a multivariate analogue of Levene's test for homogeneity of variance
betadisper <- betadisper(burn.dist.morisita, group = color.vec$p.group)
morisita.betadisper <- anova(betadisper)

#--add Betadisper results to table
anosim.res[which(anosim.res$anosim.res =="morisita"), "F.betadisper"] <-
  morisita.betadisper$`F value`[1]
anosim.res[which(anosim.res$anosim.res =="morisita"), "df.betadisper.1"] <-
  morisita.betadisper$Df[1]
anosim.res[which(anosim.res$anosim.res =="morisita"), "df.betadisper.2"] <-
  morisita.betadisper$Df[2]
anosim.res[which(anosim.res$anosim.res =="morisita"), "p.betadisper"] <-
  morisita.betadisper$`Pr(>F)`[1]

#<< ANOSIM >>-----------------------------------------------------------------------------
morisita.anosim <- anosim(burned.comm.matrix, grouping = color.vec$p.group,
                         distance = "morisita")
morisita.anosim

#--Add results to data frame
anosim.res[which(anosim.res$anosim.res =="morisita"), "r"] <- morisita.anosim$statistic
anosim.res[which(anosim.res$anosim.res =="morisita"), "p"] <- morisita.anosim$signif

#<< PERMANOVA >>--------------------------------------------------------------------------
#--adonis
morisita.adonis <- adonis(formula = burn.dist.morisita ~ Range * po4.p.ppm *
                           temp.warmest.quarter, data = burned.matrix, permutations = 1000)
morisita.adonis

#--add results to data.frame
#--Range f.model, r2, p-value
permanova.res[which(permanova.res$test == "morisita"), "F.model.range"] <-
 morisita.adonis$aov.tab$F.Model[2]
permanova.res[which(permanova.res$test == "morisita"), "r2.range"] <-
 morisita.adonis$aov.tab$R2[2]
permanova.res[which(permanova.res$test == "morisita"), "p.range"] <-
 morisita.adonis$aov.tab$`Pr(>F)`[2]

#--Phosphate f.model, r2, p-value
permanova.res[which(permanova.res$test == "morisita"), "F.model.po4"] <- 
  morisita.adonis$aov.tab$F.Model[2]
permanova.res[which(permanova.res$test == "morisita"), "r2.po4"] <- 
  morisita.adonis$aov.tab$R2[2]
permanova.res[which(permanova.res$test == "morisita"), "p.po4"] <- 
  morisita.adonis$aov.tab$`Pr(>F)`[2]

#--prec f.model, r2, p-value
permanova.res[which(permanova.res$test == "morisita"), "F.model.temp"] <- 
  morisita.adonis$aov.tab$F.Model[3]
permanova.res[which(permanova.res$test == "morisita"), "r2.temp"] <- 
  morisita.adonis$aov.tab$R2[3]
permanova.res[which(permanova.res$test == "morisita"), "p.temp"] <- 
  morisita.adonis$aov.tab$`Pr(>F)`[3]

write.csv(anosim.res, file = paste0(res.dir,'Burned_anosim.csv'), row.names = F)
write.csv(permanova.res,  file = paste0(res.dir,'Burned_permanova.csv'), row.names = F)

dev.off()

rm(anosim.res, permanova.res)

#=========================================================================================
# Jaccard based dissimilarity index: Unburned sites----
#=========================================================================================

#-----------------------------------------------------------------------------------------
# Create result tables
#-----------------------------------------------------------------------------------------
#--Table of Anosim results
anosim.res <- c("jaccard", "morisita")
anosim.res <- data.frame(anosim.res)
#--Table of PERMANOVA results
permanova.res <- data.frame(c('jaccard','morisita'))
colnames(permanova.res) <- "test"


#--isolate otu data
unburned.comm.matrix <- unburned.matrix[11:length(unburned.matrix)]

#--comment to include singletons
unburned.comm.matrix <- unburned.comm.matrix[colSums(unburned.comm.matrix) > 1]
unburned.comm.matrix <- unburned.comm.matrix[rowSums(unburned.comm.matrix) > 0, ] # remove rows with sums of 0
unburned.matrix <- unburned.matrix[row.names(unburned.comm.matrix),]

#--distance matrix using jaccard index
unburn.dist.jaccard <- vegdist(unburned.comm.matrix, method = "jaccard", binary = TRUE)

#--NMDS analysis
jaccard.otu <- metaMDS(unburn.dist.jaccard, dist = "bray", permutations = 999,
                       try = 100, trymax = 1000)

#--Stress
jaccard.otu$stress

#--add stress of NMDS to results table
anosim.res[which(anosim.res$anosim.res =="jaccard"), "stress.nmds"] <-
  jaccard.otu$stress

#--format and output NMDS plots to figure folder
jpeg(filename = paste0(fig.dir, 'unburned_NDMS.jpeg'), width = 1200, height = 450,
     quality = 100)
par(mfrow = c(1,2), "mar"=c(6, 5, 5, 3))

#--Plot NMDS of EM community based on Jaccard index and OTU abundance
plot(jaccard.otu, display = "sites", type = "n", cex.lab = 1.5,
     cex.axis = 1.5, yaxt = "n")
axis (2, at = seq (-0.4, 0.4, by = 0.2), cex.axis = 1.5, las = 2)
# colors for points
color.vec <- data.frame(color = rep (NA, nrow(unburned.matrix)),
                        p.group = unburned.matrix$Range)
color.vec[color.vec$p.group == 'santa.catalina', 'color'] <- 'grey'
color.vec[color.vec$p.group == 'pinaleno', 'color'] <- 'black'

#ordipointlabel(jaccard.otu, display = "sites")

#--isolate points using rows isolated above
points(jaccard.otu, display = "sites", cex = 1.5,
       col = color.vec$color,
       bg = color.vec$color)

legend("bottomleft", legend = levels(color.vec$p.group), bty = "n",
       col = c('black','grey'),
       pch = c(6,6,20,20),
       pt.bg =  c('black','grey'), cex = 1)
#--Ordihull variations by range, burn, and both burn and range
#ordihull(jaccard.otu, groups = color.vec$p.group)
ordihull(jaccard.otu, groups = color.vec$p.group,
         col = c('black','grey'))
#ordihull(jaccard.otu, groups = stsp.matrix$Range) # just mt. range

#--BetaDisper
#--a multivariate analogue of Levene's test for homogeneity of variance
betadisper <- betadisper(unburn.dist.jaccard, group = color.vec$p.group)
jaccard.betadisper <- anova(betadisper)

#--add Betadisper results to table
anosim.res[which(anosim.res$anosim.res =="jaccard"), "F.betadisper"] <-
  jaccard.betadisper$`F value`[1]
anosim.res[which(anosim.res$anosim.res =="jaccard"), "df.betadisper.1"] <-
  jaccard.betadisper$Df[1]
anosim.res[which(anosim.res$anosim.res =="jaccard"), "df.betadisper.2"] <-
  jaccard.betadisper$Df[2]
anosim.res[which(anosim.res$anosim.res =="jaccard"), "p.betadisper"] <-
  jaccard.betadisper$`Pr(>F)`[1]

#<< ANOSIM >>-----------------------------------------------------------------------------
jaccard.anosim <- anosim(unburned.comm.matrix, grouping = color.vec$p.group,
                         distance = "jaccard")
jaccard.anosim

#--Add results to data frame
anosim.res[which(anosim.res$anosim.res =="jaccard"), "r"] <- jaccard.anosim$statistic
anosim.res[which(anosim.res$anosim.res =="jaccard"), "p"] <- jaccard.anosim$signif

#<< PERMANOVA >>--------------------------------------------------------------------------
#--adonis
jaccard.adonis <- adonis(formula = unburn.dist.jaccard ~ Range * po4.p.ppm *
                         temp.warmest.quarter, data = unburned.matrix,
                         permutations = 1000)
jaccard.adonis

#--add results to data.frame
#--Range f.model, r2, p-value
permanova.res[which(permanova.res$test == "jaccard"), "F.model.range"] <-
 jaccard.adonis$aov.tab$F.Model[2]
permanova.res[which(permanova.res$test == "jaccard"), "r2.range"] <-
 jaccard.adonis$aov.tab$R2[2]
permanova.res[which(permanova.res$test == "jaccard"), "p.range"] <-
 jaccard.adonis$aov.tab$`Pr(>F)`[2]

#--Phosphate f.model, r2, p-value
permanova.res[which(permanova.res$test == "jaccard"), "F.model.po4"] <- 
  jaccard.adonis$aov.tab$F.Model[2]
permanova.res[which(permanova.res$test == "jaccard"), "r2.po4"] <- 
  jaccard.adonis$aov.tab$R2[2]
permanova.res[which(permanova.res$test == "jaccard"), "p.po4"] <- 
  jaccard.adonis$aov.tab$`Pr(>F)`[2]

#--prec f.model, r2, p-value
permanova.res[which(permanova.res$test == "jaccard"), "F.model.temp"] <- 
  jaccard.adonis$aov.tab$F.Model[3]
permanova.res[which(permanova.res$test == "jaccard"), "r2.temp"] <- 
  jaccard.adonis$aov.tab$R2[3]
permanova.res[which(permanova.res$test == "jaccard"), "p.temp"] <- 
  jaccard.adonis$aov.tab$`Pr(>F)`[3]

#=========================================================================================
# Morisita based dissimilarity index: Unburned sites----
#=========================================================================================
#--isolate otu data
unburned.comm.matrix <- unburned.matrix[11:length(unburned.matrix)]

#--comment to include singletons
unburned.comm.matrix <- unburned.comm.matrix[colSums(unburned.comm.matrix) > 1]
unburned.comm.matrix <- unburned.comm.matrix[rowSums(unburned.comm.matrix) > 0, ] # remove rows with sums of 0
burned.matrix <- unburned.matrix[row.names(unburned.comm.matrix),]

#--distance matrix using jaccard index
unburn.dist.morisita <- vegdist(unburned.comm.matrix, method = "horn")

#--NMDS analysis
morisita.otu <- metaMDS(unburn.dist.morisita, dist = "bray", permutations = 999,
                        try = 100, trymax = 1000)

#--Stress
morisita.otu$stress

#--add stress of NMDS to results table
anosim.res[which(anosim.res$anosim.res =="morisita"), "stress.nmds"] <-
  morisita.otu$stress

#--Plot NMDS of EM community based on Jaccard index and OTU abundance
plot(morisita.otu, display = "sites", type = "n", cex.lab = 1.5,
     cex.axis = 1.5, yaxt = "n")
axis (2, at = seq (-0.4, 0.4, by = 0.2), cex.axis = 1.5, las = 2)
# colors for points
color.vec <- data.frame(color = rep (NA, nrow(unburned.matrix)),
                        p.group = unburned.matrix$Range)
color.vec[color.vec$p.group == 'santa.catalina', 'color'] <- 'grey'
color.vec[color.vec$p.group == 'pinaleno', 'color'] <- 'black'

#ordipointlabel(jaccard.otu, display = "sites")

#--isolate points using rows isolated above
points(morisita.otu, display = "sites", cex = 1.5,
       col = color.vec$color,
       bg = color.vec$color)

legend("bottomleft", legend = levels(color.vec$p.group), bty = "n",
       col = c('black','grey'),
       pch = c(6,6,20,20),
       pt.bg =  c('black','grey'), cex = 1)
#--Ordihull variations by range, burn, and both burn and range
#ordihull(morisita.otu, groups = color.vec$p.group)
ordihull(morisita.otu, groups = color.vec$p.group,
         col = c('black','grey'))
#ordihull(morisita.otu, groups = stsp.matrix$Range) # just mt. range

#--BetaDisper
#--a multivariate analogue of Levene's test for homogeneity of variance
betadisper <- betadisper(unburn.dist.morisita, group = color.vec$p.group)
morisita.betadisper <- anova(betadisper)

#--add Betadisper results to table
anosim.res[which(anosim.res$anosim.res =="morisita"), "F.betadisper"] <-
  morisita.betadisper$`F value`[1]
anosim.res[which(anosim.res$anosim.res =="morisita"), "df.betadisper.1"] <-
  morisita.betadisper$Df[1]
anosim.res[which(anosim.res$anosim.res =="morisita"), "df.betadisper.2"] <-
  morisita.betadisper$Df[2]
anosim.res[which(anosim.res$anosim.res =="morisita"), "p.betadisper"] <-
  morisita.betadisper$`Pr(>F)`[1]

#<< ANOSIM >>-----------------------------------------------------------------------------
morisita.anosim <- anosim(unburned.comm.matrix, grouping = color.vec$p.group,
                          distance = "horn")
morisita.anosim

#--Add results to data frame
anosim.res[which(anosim.res$anosim.res =="morisita"), "r"] <- morisita.anosim$statistic
anosim.res[which(anosim.res$anosim.res =="morisita"), "p"] <- morisita.anosim$signif

#<< PERMANOVA >>--------------------------------------------------------------------------
#--adonis
morisita.adonis <- adonis(formula = unburn.dist.morisita ~ Range * po4.p.ppm *
                            temp.warmest.quarter, data = burned.matrix, permutations = 1000)
morisita.adonis

#--add results to data.frame
#--Range f.model, r2, p-value
permanova.res[which(permanova.res$test == "morisita"), "F.model.range"] <-
  morisita.adonis$aov.tab$F.Model[2]
permanova.res[which(permanova.res$test == "morisita"), "r2.range"] <-
  morisita.adonis$aov.tab$R2[2]
permanova.res[which(permanova.res$test == "morisita"), "p.range"] <-
  morisita.adonis$aov.tab$`Pr(>F)`[2]

#--Phosphate f.model, r2, p-value
permanova.res[which(permanova.res$test == "morisita"), "F.model.po4"] <- 
  morisita.adonis$aov.tab$F.Model[2]
permanova.res[which(permanova.res$test == "morisita"), "r2.po4"] <- 
  morisita.adonis$aov.tab$R2[2]
permanova.res[which(permanova.res$test == "morisita"), "p.po4"] <- 
  morisita.adonis$aov.tab$`Pr(>F)`[2]

#--prec f.model, r2, p-value
permanova.res[which(permanova.res$test == "morisita"), "F.model.temp"] <- 
  morisita.adonis$aov.tab$F.Model[3]
permanova.res[which(permanova.res$test == "morisita"), "r2.temp"] <- 
  morisita.adonis$aov.tab$R2[3]
permanova.res[which(permanova.res$test == "morisita"), "p.temp"] <- 
  morisita.adonis$aov.tab$`Pr(>F)`[3]

write.csv(anosim.res, file = paste0(res.dir,'Unburned_anosim.csv'), row.names = F)
write.csv(permanova.res,  file = paste0(res.dir,'Unburned_permanova.csv'), row.names = F)

dev.off()
