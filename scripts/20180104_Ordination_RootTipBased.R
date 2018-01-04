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
# Load data and clean up----
#----------------------------------------------------------------------------------------#
#--Load data
all.data <- read.csv(paste0(dat.dir,'20170806_OTU_data.csv'), as.is = T)

#--Remove sequences not assigned to an OTU
otu.data <- all.data[!is.na(all.data$otu.97), ]

#--Remove sequences not assigned to Ponderosa host
otu.data <- otu.data[otu.data$Host == 'Ponderosa',]

#--Creates matrix, grouping OTUs by tree number and getting a sum of EM tip abundance 
# for each OTU
otu.data %>%
  select(Sample_name,Burn_status,Range,Site,Tree,otu.97,Tip_count) %>%
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

#<< Isolate soil data for PCA >> ----------------------
soil.data <- read.csv(paste0(dat.dir, 'soil_data.csv'),as.is = T)
# soil.sig <- c('site','pH.su','po4.p.ppm','so4.s.ppm','b.ppm','mg.ppm',
#               'k.ppm','no3.n.ppm')
all.soil <- c('site','pH.su','po4.p.ppm','so4.s.ppm','b.ppm','mg.ppm',
              'k.ppm','no3.n.ppm','EC.ds.m','ca.ppm','na.ppm','zn.ppm','fe.ppm','mn.ppm',
              'cu.ppm','ni.ppm','cec.meq.100g')
soil.data.sig <- soil.data[colnames(soil.data)%in% all.soil]
soil.data.sig[soil.data.sig$no3.n.ppm == '<1.0', 'no3.n.ppm'] <- 0
soil.data.sig$no3.n.ppm <- as.numeric(soil.data.sig$no3.n.ppm)
for(s in unique(soil.data$site)){
  for(f in all.soil[2:length(all.soil)]){
    stsp.matrix[stsp.matrix$Site == s, f] <- mean(soil.data.sig[soil.data.sig$site == s, f])
  }
}

#--reorder matrix
colnames(stsp.matrix)
soil.perm <- stsp.matrix[c(1:4,123:138)]
stsp.matrix <- stsp.matrix[c(1:4,121:122,5:120)]

#----------------------------------------------------------------------------------------#
# PCA of soil characteristics
#----------------------------------------------------------------------------------------#
soil.perma <- soil.perm[5:length(soil.perm)]
soil.pca <- prcomp(soil.perma,
                   center = TRUE,
                   scale. = TRUE) 
#print(soil.pca)
#plot(soil.pca, type = 'l')
summary(soil.pca)
soil.eigenvector <- scores(soil.pca, choices = 1)
soil.eigenvector <- data.frame(matrix(unlist(soil.eigenvector), nrow=41, byrow=T))
stsp.matrix$soil.pca <- soil.eigenvector$matrix.unlist.soil.eigenvector...nrow...41..byrow...T.

stsp.matrix <- stsp.matrix[c(1:6,123,7:122)]

#----------------------------------------------------------------------------------------#
# Create result tables
#----------------------------------------------------------------------------------------#
#--Table of Anosim results
anosim.res <- c("jaccard", "morisita")
anosim.res <- data.frame(anosim.res)
#--Table of PERMANOVA results
permanova.res <- data.frame(c('jaccard','morisita'))
colnames(permanova.res) <- "test"

#----------------------------------------------------------------------------------------#
# Create distinct tables for burned sites only and unburned only
#----------------------------------------------------------------------------------------#
burned.matrix <- stsp.matrix[stsp.matrix$Burn_status == 'burned',]
unburned.matrix <- stsp.matrix[!stsp.matrix$Burn_status == 'burned',]

#========================================================================================#
# Jaccard based dissimilarity index----
#========================================================================================#

#--comment to not remove outlier (LB056)
#jaccard.matrix <- stsp.matrix[!stsp.matrix$Tree %in% c('LB056'),]

#--isolate otu data
jaccard.matrix <- stsp.matrix
comm.matrix <- jaccard.matrix[8:length(jaccard.matrix)]

#--comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 4]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 0, ] # remove rows with sums of 0
jaccard.matrix <- jaccard.matrix[row.names(comm.matrix),]

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
jpeg(filename = paste0(fig.dir, 'NDMS_RootTipBased.jpeg'),
     width = 1200, height = 450,
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
legend("topright", legend = levels(color.vec$p.group), bty = "n",
       col = c('black','green4','black','green4'),
       pch = c(6,6,20,20),
       pt.bg =  c('black','green4','black','green4'), cex = 1)

#--Ordihull variations by range, burn, and both burn and range
ordihull(jaccard.otu, groups = color.vec$p.group,
         col = c('black','green4','black','green4'))
#ordiellipse(jaccard.otu, groups = color.vec$p.group,
#         col = c('black','green4','black','green4'),
#         kind = 'ehull')
#ordihull(jaccard.otu, groups = stsp.matrix$Range) # just mt. range
#--Overlays
fit <- envfit(jaccard.otu ~ soil.pca + Tmax + prec, jaccard.matrix)
fit
plot(fit, type = 'n')

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
jaccard.adonis <- adonis(formula = comm.dist.jaccard ~ Burn_status * Range * Tmax,
                         data = jaccard.matrix, permutations = 1000)
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
permanova.res[which(permanova.res$test == "jaccard"), "F.model.range"] <-
  jaccard.adonis$aov.tab$F.Model[2]
permanova.res[which(permanova.res$test == "jaccard"), "r2.range"] <-
  jaccard.adonis$aov.tab$R2[2]
permanova.res[which(permanova.res$test == "jaccard"), "p.range"] <-
  jaccard.adonis$aov.tab$`Pr(>F)`[2]

#--Tmax f.model, r2, p-value
permanova.res[which(permanova.res$test == "jaccard"), "F.model.temp"] <- 
  jaccard.adonis$aov.tab$F.Model[3]
permanova.res[which(permanova.res$test == "jaccard"), "r2.temp"] <- 
  jaccard.adonis$aov.tab$R2[3]
permanova.res[which(permanova.res$test == "jaccard"), "p.temp"] <- 
  jaccard.adonis$aov.tab$`Pr(>F)`[3]

#========================================================================================#
# Morisita based dissimilarity index----
#========================================================================================#

#--comment to not remove outlier (NF19)
#morisitia.matrix <- stsp.matrix[!stsp.matrix$Tree %in% c('NF19'),]

#--isolate site X species matrix only without metadata
morisita.matrix <- stsp.matrix
comm.matrix <- morisita.matrix[8:length(morisita.matrix)]

#--comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 4]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 0, ] # remove rows with sums of 0
morisita.matrix <- stsp.matrix[row.names(comm.matrix),]

#--distance matrix using jaccard index
comm.dist.morisita <- vegdist(comm.matrix, method = "horn", binary = F)

#--NMDS analysis
morisita.otu <- metaMDS(comm.dist.morisita, dist = "bray", permutations = 999,
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
ordihull(morisita.otu, groups = color.vec$p.group,
         col = c('black','green4','black','green4'))

#--Overlays
fit <- envfit(morisita.otu ~ soil.pca + Tmax + prec, morisita.matrix)
fit
plot(fit, type = 'n')

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

#<< PERMANOVA >>-----------------------------------------------------------------------------
morisita.anosim <- adonis(comm.matrix ~ Burn_status * Range, data = morisita.matrix,
                          distance = "horn", by = NULL)
morisita.anosim

#--Add results to data frame
anosim.res[which(anosim.res$anosim.res =="morisita"), "r"] <- morisita.anosim$statistic
anosim.res[which(anosim.res$anosim.res =="morisita"), "p"] <- morisita.anosim$signif

#<< PERMANOVA >>--------------------------------------------------------------------------
#--adonis
morisita.adonis <- adonis(formula = comm.dist.morisita ~ Burn_status * Range * Tmax,
                          data = morisita.matrix, permutations = 1000)
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
permanova.res[which(permanova.res$test == "morisita"), "F.model.range"] <-
  morisita.adonis$aov.tab$F.Model[2]
permanova.res[which(permanova.res$test == "morisita"), "r2.range"] <-
  morisita.adonis$aov.tab$R2[2]
permanova.res[which(permanova.res$test == "morisita"), "p.range"] <-
  morisita.adonis$aov.tab$`Pr(>F)`[2]

#--Tmax f.model, r2, p-value
permanova.res[which(permanova.res$test == "morisita"), "F.model.temp"] <- 
  morisita.adonis$aov.tab$F.Model[3]
permanova.res[which(permanova.res$test == "morisita"), "r2.temp"] <- 
  morisita.adonis$aov.tab$R2[3]
permanova.res[which(permanova.res$test == "morisita"), "p.temp"] <- 
  morisita.adonis$aov.tab$`Pr(>F)`[3]

#--export anosim table
write.csv(anosim.res, paste0(res.dir, "ANOSIM_RootTipBased.csv"),
          row.names = F)

#--output permanova results
write.csv(permanova.res, paste0(res.dir, "Permanova_RootTipBaed.csv"),
          row.names = F)


#========================================================================================#
# Jaccard based dissimilarity index: Burned----
#========================================================================================#

#--isolate otu data
comm.matrix <- burned.matrix[11:length(burned.matrix)]

#--comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 4]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 0, ] # remove rows with sums of 0
jaccard.matrix <- burned.matrix[row.names(comm.matrix),]

#--distance matrix using jaccard index
comm.dist.jaccard <- vegdist(comm.matrix, method = "jaccard", binary = TRUE)

#--NMDS analysis
jaccard.otu <- metaMDS(comm.dist.jaccard, dist = "bray", permutations = 999,
                       try = 100, trymax = 1000)

#--Stress
jaccard.otu$stress

#--format and output NMDS plots to figure folder
jpeg(filename = paste0(fig.dir, 'NDMS_burned.jpeg'), width = 1200, height = 450,
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

color.vec[color.vec$p.group == 'santa.catalina burned', 'shape'] <- 20
color.vec[color.vec$p.group == 'pinaleno burned', 'shape'] <- 6

#ordipointlabel(jaccard.otu, display = "sites")
#--isolate points for pinaleno mts.
pinaleno <- row.names(jaccard.matrix[jaccard.matrix$Range == 'pinaleno',])
santa.catalina <- row.names(jaccard.matrix[jaccard.matrix$Range == 'santa.catalina',])

#--isolate points using rows isolated above
points(jaccard.otu$points[pinaleno,1:2], display = "sites", cex = 1.5,
       pch = color.vec[color.vec$p.group == 'pinaleno burned',
                       'shape'],
       col = color.vec[color.vec$p.group == 'pinaleno burned',
                       'color'],
       bg = color.vec[color.vec$p.group == 'pinaleno burned',
                      'color'])
points(jaccard.otu$points[santa.catalina,1:2], display = "sites", cex = 1.5,
       pch = color.vec[color.vec$p.group == 'santa.catalina burned',
                       'shape'],
       col = color.vec[color.vec$p.group == 'santa.catalina burned',
                       'color'],
       bg = color.vec[color.vec$p.group == 'santa.catalina burned',
                      'color'])
legend("topright", legend = levels(color.vec$p.group), bty = "n",
       col = c('black','black'),
       pch = c(6,20),
       pt.bg =  c('black','black'), cex = 1)

#--Ordihull variations by range, burn, and both burn and range
#ordihull(jaccard.otu, groups = color.vec$p.group)
ordiellipse(jaccard.otu, groups = color.vec$p.group,
            col = 'black',
            kind = 'ehull')
#ordihull(jaccard.otu, groups = stsp.matrix$Range) # just mt. range
#--Overlays
# fit <- envfit(jaccard.otu ~ po4.p.ppm + temp.warmest.quarter +
#                 pH.su, jaccard.matrix)
# fit
# plot(fit, type = 'n')

#--BetaDisper
#--a multivariate analogue of Levene's test for homogeneity of variance
betadisper <- betadisper(comm.dist.jaccard, group = color.vec$p.group)
jaccard.betadisper <- anova(betadisper)
jaccard.betadisper

#<< ANOSIM >>-----------------------------------------------------------------------------
jaccard.anosim <- anosim(comm.matrix, grouping = ,
                         distance = "jaccard")
jaccard.anosim

#========================================================================================#
# Morisita based dissimilarity index: Burned----
#========================================================================================#

# #--comment to not remove outlier (NF19)
# stsp.matrix <- stsp.matrix[!stsp.matrix$Tree %in% c('NF19'),]

#--isolate site X species matrix only without metadata
comm.matrix <- burned.matrix[11:length(burned.matrix)]

#--comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 4]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 0, ] # remove rows with sums of 0
morisita.matrix <- burned.matrix[row.names(comm.matrix),]

#--distance matrix using jaccard index
comm.dist.morisita <- vegdist(comm.matrix, method = "horn", binary = F)

#--NMDS analysis
morisita.otu <- metaMDS(comm.dist.morisita, dist = "bray", permutations = 999,
                        try = 100, trymax = 1000)

#--Stress
stressplot(morisita.otu)
morisita.otu$stress


#--Plot NMDS of EM community based on Jaccard index and OTU abundance
plot(morisita.otu, display = "sites", type = "n", cex.lab = 1.5,
     cex.axis = 1.5)
# colors for points
color.vec <- data.frame (color = rep (NA, nrow(morisita.matrix)),
                         p.group = paste(morisita.matrix$Range,
                                         morisita.matrix$Burn_status))
color.vec[color.vec$p.group == 'santa.catalina burned', 'color'] <- 'black'
color.vec[color.vec$p.group == 'pinaleno burned', 'color'] <- 'black'
color.vec[color.vec$p.group == 'santa.catalina burned', 'shape'] <- 20
color.vec[color.vec$p.group == 'pinaleno burned', 'shape'] <- 6

#ordipointlabel(morisita.otu, display = "sites")

#--isolate points for pinaleno mts.
pinaleno <- row.names(morisita.matrix[morisita.matrix$Range == 'pinaleno',])
santa.catalina <- row.names(morisita.matrix[morisita.matrix$Range == 'santa.catalina',])
#--isolate points using rows isolated above
points(morisita.otu$points[pinaleno,1:2], display = "sites", cex = 1.5,
       pch = color.vec[color.vec$p.group == 'pinaleno burned',
                       'shape'],
       col = color.vec[color.vec$p.group == 'pinaleno burned',
                       'color'],
       bg = color.vec[color.vec$p.group == 'pinaleno burned',
                      'color'])
points(morisita.otu$points[santa.catalina,1:2], display = "sites", cex = 1.5,
       pch = color.vec[color.vec$p.group == 'santa.catalina burned',
                       'shape'],
       col = color.vec[color.vec$p.group == 'santa.catalina burned',
                       'color'],
       bg = color.vec[color.vec$p.group == 'santa.catalina burned',
                      'color'])
legend("bottomleft", legend = levels(color.vec$p.group), bty = "n",
       col = c('black','black'),
       pch = c(6,20),
       pt.bg =  c('black','black'), cex = 1)

#--Ordihull variations by range, burn, and both burn and range
#ordihull(morisita.otu, groups = color.vec$p.group,
#         col = c('black','green4','black','green4')),
ordiellipse(morisita.otu, groups = color.vec$p.group,
            col = 'black',
            kind = 'ehull') # just burn status
#ordihull(morisita.otu, groups = stsp.matrix$Range,
#         col = c('black','green4','black','green4')) # just mt. range


#--BetaDisper
#--a multivariate analogue of Levene's test for homogeneity of variance
betadisper <- betadisper(comm.dist.morisita, group = color.vec$p.group)
morisita.betadisper <- anova(betadisper)
morisita.betadisper

#<< ANOSIM >>-----------------------------------------------------------------------------
morisita.anosim <- anosim(comm.matrix, grouping = color.vec$p.group,
                          distance = "horn")
morisita.anosim

#dev.off()
