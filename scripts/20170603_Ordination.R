## Script created by Liz Bowman June 3, 2017
## for analyzing community differences between ranges and burn status of sites

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
dat.dir <- "~/Documents/PhD/2_EM_Fire_effect/data/"
fig.dir <- '~/Documents/PhD/2_EM_Fire_effect/figures/'
res.dir <- "~/Documents/PhD/2_EM_Fire_effect/results/"

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
  select(Sample_name,Burn_status,Range,Site,Tree,otu.97,otu.count) %>%
  spread(otu.97,otu.count) %>%
  group_by(Tree,Site,Range,Burn_status) %>%
  select(matches ("otu*")) %>%
  summarize_each(funs(sum(., na.rm = TRUE))) %>%
  as.data.frame() -> stsp.matrix

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
# Jaccard based dissimilarity index----
#=========================================================================================

#--comment to not remove outlier (LB056)
stsp.matrix <- stsp.matrix[!stsp.matrix$Tree %in% c('LB056'),]

#--isolate otu data
comm.matrix <- stsp.matrix[11:length(stsp.matrix)]

#--comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) > 4]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 0, ] # remove rows with sums of 0
jaccard.matrix <- stsp.matrix[row.names(comm.matrix),]

#--distance matrix using jaccard index
comm.dist.jaccard <- vegdist(comm.matrix, method = "jaccard", binary = TRUE)

#--NMDS analysis
jaccard.otu <- metaMDS(comm.dist.jaccard, dist = "bray", permutations = 999,
                       try = 100, trymax = 1000)

#--Stress
jaccard.otu$stress

#--add stress of NMDS to results table
anosim.res[which(anosim.res$anosim.res =="jaccard"), "stress.nmds"] <-
  jaccard.otu$stress

#--format and output NMDS plots to figure folder
jpeg(filename = paste0(fig.dir, 'NDMS.jpeg'), width = 1200, height = 450,
     quality = 100)
par(mfrow = c(1,2), "mar"=c(6, 5, 5, 3))

#--Plot NMDS of EM community based on Jaccard index and OTU abundance
plot(jaccard.otu, display = "sites", type = "n", cex.lab = 1.5,
     cex.axis = 1.5, yaxt = "n")
axis (2, at = seq (-0.4, 0.4, by = 0.2), cex.axis = 1.5, las = 2)
# colors for points
color.vec <- data.frame(color = rep(NA, nrow(jaccard.matrix)),
                         p.group = paste(jaccard.matrix$Range,
                                         jaccard.matrix$Burn_status))
color.vec[color.vec$p.group == 'santa.catalina burned', 'color'] <- 'black'
color.vec[color.vec$p.group == 'pinaleno burned', 'color'] <- 'black'
color.vec[color.vec$p.group == 'santa.catalina unburned', 'color'] <- 'green4'
color.vec[color.vec$p.group == 'pinaleno unburned', 'color'] <- 'green4'

color.vec[color.vec$p.group == 'santa.catalina burned', 'shape'] <- 20
color.vec[color.vec$p.group == 'pinaleno burned', 'shape'] <- 6
color.vec[color.vec$p.group == 'santa.catalina unburned', 'shape'] <- 20
color.vec[color.vec$p.group == 'pinaleno unburned', 'shape'] <- 6

#ordipointlabel(jaccard.otu, display = "sites")
#--isolate points for pinaleno mts.
pinaleno <- row.names(jaccard.matrix[jaccard.matrix$Range == 'pinaleno',])
santa.catalina <- row.names(jaccard.matrix[jaccard.matrix$Range == 'santa.catalina',])

#--isolate points using rows isolated above
points(jaccard.otu$points[pinaleno,1:2], display = "sites", cex = 1.5,
       pch = color.vec[color.vec$p.group == 'pinaleno burned' | 
                         color.vec$p.group == 'pinaleno unburned',
                       'shape'],
       col = color.vec[color.vec$p.group == 'pinaleno burned' | 
                         color.vec$p.group == 'pinaleno unburned',
                       'color'],
       bg = color.vec[color.vec$p.group == 'pinaleno burned' | 
                        color.vec$p.group == 'pinaleno unburned',
                      'color'])
points(jaccard.otu$points[santa.catalina,1:2], display = "sites", cex = 1.5,
       pch = color.vec[color.vec$p.group == 'santa.catalina burned' | 
                         color.vec$p.group == 'santa.catalina unburned',
                       'shape'],
       col = color.vec[color.vec$p.group == 'santa.catalina burned' |
                         color.vec$p.group == 'santa.catalina unburned',
                       'color'],
       bg = color.vec[color.vec$p.group == 'santa.catalina burned' |
                        color.vec$p.group == 'santa.catalina unburned',
                      'color'])
legend("topleft", legend = levels(color.vec$p.group), bty = "n",
       col = c('black','green4','black','green4'),
       pch = c(6,6,20,20),
       pt.bg =  c('black','green4','black','green4'), cex = 1)
#--Ordihull variations by range, burn, and both burn and range
#ordihull(jaccard.otu, groups = color.vec$p.group)
ordihull(jaccard.otu, groups = color.vec$p.group,
         col = c('black','green4','black','green4'))
#ordihull(jaccard.otu, groups = stsp.matrix$Range) # just mt. range

#--BetaDisper
#--a multivariate analogue of Levene's test for homogeneity of variance
betadisper <- betadisper(comm.dist.jaccard, group = color.vec$p.group)
jaccard.betadisper <- anova(betadisper)
jaccard.betadisper

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
jaccard.anosim <- anosim(comm.matrix, grouping = color.vec$p.group,
                               distance = "jaccard")
jaccard.anosim

#--Add results to data frame
anosim.res[which(anosim.res$anosim.res =="jaccard"), "r"] <- jaccard.anosim$statistic
anosim.res[which(anosim.res$anosim.res =="jaccard"), "p"] <- jaccard.anosim$signif

#<< PERMANOVA >>--------------------------------------------------------------------------
#--adonis
jaccard.adonis <- adonis(formula = comm.dist.jaccard ~ Burn_status * po4.p.ppm *
                  temp.warmest.quarter, data = jaccard.matrix, permutations = 1000)
jaccard.adonis

#--add results to data.frame
#--Burn status f.model, r2, p-value
permanova.res[which(permanova.res$test == "jaccard"), "F.model.burn.status"] <- 
  jaccard.adonis$aov.tab$F.Model[1]
permanova.res[which(permanova.res$test == "jaccard"), "r2.burn.status"] <- 
  jaccard.adonis$aov.tab$R2[1]
permanova.res[which(permanova.res$test == "jaccard"), "p.burn.status"] <- 
  jaccard.adonis$aov.tab$`Pr(>F)`[1]

#--Range f.model, r2, p-value
#permanova.res[which(permanova.res$test == "jaccard"), "F.model.range"] <- 
#  jaccard.adonis$aov.tab$F.Model[2]
#permanova.res[which(permanova.res$test == "jaccard"), "r2.range"] <- 
#  jaccard.adonis$aov.tab$R2[2]
#permanova.res[which(permanova.res$test == "jaccard"), "p.range"] <- 
#  jaccard.adonis$aov.tab$`Pr(>F)`[2]

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
# Morisita based dissimilarity index----
#=========================================================================================

# #--comment to not remove outlier (NF19)
# stsp.matrix <- stsp.matrix[!stsp.matrix$Tree %in% c('NF19'),]

#--isolate site X species matrix only without metadata
comm.matrix <- stsp.matrix[11:length(stsp.matrix)]

#--comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) > 1]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 0, ] # remove rows with sums of 0
morisita.matrix <- stsp.matrix[row.names(comm.matrix),]

#--distance matrix using jaccard index
comm.dist.morisita <- vegdist(comm.matrix, method = "horn", binary = F)

#--NMDS analysis
morisita.otu <- metaMDS(comm.dist.morisita, dist = "bray", permutations = 999,
                       try = 100, trymax = 1000)

#--Stress
stressplot(morisita.otu)
morisita.otu$stress

#--add stress of NMDS to results table
anosim.res[which(anosim.res$anosim.res =="morisita"), "stress.nmds"] <-
  morisita.otu$stress

#--Plot NMDS of EM community based on Jaccard index and OTU abundance
plot(morisita.otu, display = "sites", type = "n", cex.lab = 1.5,
     cex.axis = 1.5, yaxt = "n")
axis (2, at = seq (-0.4, 0.4, by = 0.2), cex.axis = 1.5, las = 2)
# colors for points
color.vec <- data.frame (color = rep (NA, nrow(morisita.matrix)),
                         p.group = paste(morisita.matrix$Range,
                                         morisita.matrix$Burn_status))
color.vec[color.vec$p.group == 'santa.catalina burned', 'color'] <- 'black'
color.vec[color.vec$p.group == 'pinaleno burned', 'color'] <- 'black'
color.vec[color.vec$p.group == 'santa.catalina unburned', 'color'] <- 'green4'
color.vec[color.vec$p.group == 'pinaleno unburned', 'color'] <- 'green4'

color.vec[color.vec$p.group == 'santa.catalina burned', 'shape'] <- 20
color.vec[color.vec$p.group == 'pinaleno burned', 'shape'] <- 6
color.vec[color.vec$p.group == 'santa.catalina unburned', 'shape'] <- 20
color.vec[color.vec$p.group == 'pinaleno unburned', 'shape'] <- 6

#ordipointlabel(morisita.otu, display = "sites")

#--isolate points for pinaleno mts.
pinaleno <- row.names(morisita.matrix[morisita.matrix$Range == 'pinaleno',])
santa.catalina <- row.names(morisita.matrix[morisita.matrix$Range == 'santa.catalina',])
#--isolate points using rows isolated above
points(morisita.otu$points[pinaleno,1:2], display = "sites", cex = 1.5,
       pch = color.vec[color.vec$p.group == 'pinaleno burned' | 
                         color.vec$p.group == 'pinaleno unburned',
                       'shape'],
       col = color.vec[color.vec$p.group == 'pinaleno burned' | 
                         color.vec$p.group == 'pinaleno unburned',
                       'color'],
       bg = color.vec[color.vec$p.group == 'pinaleno burned' | 
                        color.vec$p.group == 'pinaleno unburned',
                      'color'])
points(morisita.otu$points[santa.catalina,1:2], display = "sites", cex = 1.5,
       pch = color.vec[color.vec$p.group == 'santa.catalina burned' | 
                         color.vec$p.group == 'santa.catalina unburned',
                       'shape'],
       col = color.vec[color.vec$p.group == 'santa.catalina burned' |
                         color.vec$p.group == 'santa.catalina unburned',
                       'color'],
       bg = color.vec[color.vec$p.group == 'santa.catalina burned' |
                        color.vec$p.group == 'santa.catalina unburned',
                      'color'])
# legend("bottomright", legend = levels(color.vec$p.group), bty = "n",
#        col = c('black','green4','black','green4'),
#        pch = c(6,6,20,20),
#        pt.bg =  c('black','green4','black','green4'), cex = 1)

#--Ordihull variations by range, burn, and both burn and range
#ordihull(morisita.otu, groups = color.vec$p.group,
#         col = c('black','green4','black','green4')),
ordihull(morisita.otu, groups = color.vec$p.group,
         col = c('black','green4','black','green4')) # just burn status
#ordihull(morisita.otu, groups = stsp.matrix$Range,
#         col = c('black','green4','black','green4')) # just mt. range

dev.off()

#--BetaDisper
#--a multivariate analogue of Levene's test for homogeneity of variance
betadisper <- betadisper(comm.dist.morisita, group = color.vec$p.group)
morisita.betadisper <- anova(betadisper)
morisita.betadisper

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
morisita.anosim <- anosim(comm.matrix, grouping = color.vec$p.group,
                         distance = "horn")
morisita.anosim

#--Add results to data frame
anosim.res[which(anosim.res$anosim.res =="morisita"), "r"] <- morisita.anosim$statistic
anosim.res[which(anosim.res$anosim.res =="morisita"), "p"] <- morisita.anosim$signif

#<< PERMANOVA >>--------------------------------------------------------------------------
#--adonis
morisita.adonis <- adonis(formula = comm.dist.morisita ~ Burn_status * po4.p.ppm *
                           temp.warmest.quarter, data = morisita.matrix, permutations = 1000)
morisita.adonis

#--add results to data.frame
#--Burn status f.model, r2, p-value
permanova.res[which(permanova.res$test == "morisita"), "F.model.burn.status"] <- 
  morisita.adonis$aov.tab$F.Model[1]
permanova.res[which(permanova.res$test == "morisita"), "r2.burn.status"] <- 
  morisita.adonis$aov.tab$R2[1]
permanova.res[which(permanova.res$test == "morisita"), "p.burn.status"] <- 
  morisita.adonis$aov.tab$`Pr(>F)`[1]
#--Range f.model, r2, p-value
#permanova.res[which(permanova.res$test == "morisita"), "F.model.range"] <- 
#  morisita.adonis$aov.tab$F.Model[2]
#permanova.res[which(permanova.res$test == "morisita"), "r2.range"] <- 
#  morisita.adonis$aov.tab$R2[2]
#permanova.res[which(permanova.res$test == "morisita"), "p.range"] <- 
#  morisita.adonis$aov.tab$`Pr(>F)`[2]

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

#--export anosim table
write.csv(anosim.res, paste0(res.dir, "comm_comp_anosim_res.csv"),
          row.names = F)

#--output permanova results
write.csv(permanova.res, paste0(res.dir, "comm_comp_permanova_res.csv"),
          row.names = F)
