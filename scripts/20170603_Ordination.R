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
# PCA of soil characteristics: by tree
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
anosim.res <- c("jaccard.all", "morisita.all",
                'jaccard.p', 'morisita.p',
                'jaccard.scm', 'morisita.scm')
anosim.res <- data.frame(anosim.res)
#--Table of PERMANOVA results
permanova.res <- data.frame(c('jaccard','morisita'))
colnames(permanova.res) <- "test"


#--outlier from Santa Catalina Mts. (LB056)
stsp.matrix <- stsp.matrix[!stsp.matrix$Tree == 'LB056',]

#----------------------------------------------------------------------------------------#
# Create distinct tables for each range
#----------------------------------------------------------------------------------------#
pinaleno.matrix <- stsp.matrix[stsp.matrix$Range == 'pinaleno',]
santa.matrix <- stsp.matrix[stsp.matrix$Range == 'santa.catalina',]

#========================================================================================#
# Jaccard based dissimilarity index: tree level----
#========================================================================================#

#--comment to not remove outliers
jaccard.matrix <- stsp.matrix[!stsp.matrix$Tree %in% c('NF16', 'NF19'),]

#--isolate otu data
comm.matrix <- jaccard.matrix[8:length(jaccard.matrix)]

#--comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 4]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 0, ] # remove rows with sums of 0
jaccard.matrix <- jaccard.matrix[row.names(comm.matrix),]

#--distance matrix using jaccard index
comm.dist.jaccard <- vegdist(comm.matrix, method = "jaccard", binary = TRUE)

#--NMDS analysis
jaccard.otu <- metaMDS(comm.dist.jaccard, dist = "bray", permutations = 999,
                       try = 500, trymax = 1000)

#--Stress
jaccard.otu$stress

#--add stress of NMDS to results table
anosim.res[which(anosim.res$anosim.res =="jaccard.all"), "stress.nmds"] <-
  jaccard.otu$stress

#--format and output NMDS plots to figure folder
jpeg(filename = paste0('/Volumes/Cenococcum/PhD/Dissertation/Chpt.2/Figures/',
                       'Fig3c_jaccard.jpeg'),
     width = 700, height = 600,
     quality = 100)
par(mfrow = c(1,1), "mar"=c(6, 5, 5, 3))

#--Plot NMDS of EM community based on Jaccard index and OTU abundance
plot(jaccard.otu, display = "sites", type = "n", cex.lab = 2.5,
     cex.axis = 2.5, xlab = 'Axis 1',ylab =  'Axis 2')
# colors for points
color.vec <- data.frame(color = rep(NA, nrow(jaccard.matrix)),
                         p.group = jaccard.matrix$Range)
color.vec[color.vec$p.group == 'santa.catalina', 'shape'] <- 16
color.vec[color.vec$p.group == 'pinaleno', 'shape'] <- 15

#ordipointlabel(jaccard.otu, display = "sites")
#--isolate points for pinaleno mts.
pinaleno.burn <- row.names(jaccard.matrix[jaccard.matrix$Range == 'pinaleno' & 
                                            jaccard.matrix$Burn_status == 'burned',])
pinaleno.unburn <- row.names(jaccard.matrix[jaccard.matrix$Range == 'pinaleno' & 
                                              jaccard.matrix$Burn_status == 'unburned',])
santa.catalina.burn <- row.names(jaccard.matrix[jaccard.matrix$Range == 'santa.catalina' & 
                                                  jaccard.matrix$Burn_status == 'burned',])
santa.catalina.unburn <- row.names(jaccard.matrix[jaccard.matrix$Range == 'santa.catalina' & 
                                                    jaccard.matrix$Burn_status == 'unburned',])

#--isolate points using rows isolated above
points(jaccard.otu$points[pinaleno.burn,1:2], display = "sites", cex = 3,
       pch = color.vec[color.vec$p.group == 'pinaleno' | 
                         color.vec$p.group == 'pinaleno',
                       'shape'],
       col = 'black',
       bg = 'black')
points(jaccard.otu$points[pinaleno.unburn,1:2], display = "sites", cex = 3,
       pch = color.vec[color.vec$p.group == 'pinaleno' | 
                         color.vec$p.group == 'pinaleno',
                       'shape'],
       col = 'darkgreen',
       bg = 'darkgreen')
points(jaccard.otu$points[santa.catalina.burn,1:2], display = "sites", cex = 3,
       pch = color.vec[color.vec$p.group == 'santa.catalina' | 
                         color.vec$p.group == 'santa.catalina',
                       'shape'],
       col = 'black',
       bg = 'black')
points(jaccard.otu$points[santa.catalina.unburn,1:2], display = "sites", cex = 3,
       pch = color.vec[color.vec$p.group == 'santa.catalina' | 
                         color.vec$p.group == 'santa.catalina',
                       'shape'],
       col = 'darkgreen',
       bg = 'darkgreen')
# legend
legend("bottomleft", legend = c('Pinaleno Mts.',
                             'Santa Catalina Mts.'), bty = "n",
       col = c('black', 'black'),
       pch = c(15,16),
       pt.bg =  c('black','black'), cex = 2)

legend("topleft", legend = c('Burned',
                             'Unburned'), bty = "n",
       col = c('black', 'darkgreen'),
       pch = c(16,16),
       pt.bg =  c('black', 'darkgreen'), cex = 2)

#--Ordihull variations by range, burn, and both burn and range
for(i in 1:nrow(jaccard.matrix)){
  row.i <- paste(jaccard.matrix[i, 'Range'],
                 jaccard.matrix[i, 'Burn_status'])
  jaccard.matrix[i,'range_burn'] <- row.i
}
ordiellipse(jaccard.otu,
            groups = jaccard.matrix$range_burn,
            col = c('black','black'),
            kind = 'ehull')

#ordihull(jaccard.otu, groups = stsp.matrix$Range) # just mt. range
#--Overlays
fit <- envfit(jaccard.otu ~ soil.pca + Tmax + prec, jaccard.matrix)
fit
#plot(fit, type = 'n')

dev.off()

#--BetaDisper
#--a multivariate analogue of Levene's test for homogeneity of variance
betadisper <- betadisper(comm.dist.jaccard, group = jaccard.matrix$Burn_status)
jaccard.betadisper <- anova(betadisper)
jaccard.betadisper

#--add Betadisper results to table
anosim.res[which(anosim.res$anosim.res =="jaccard.all"), "F.betadisper"] <-
  jaccard.betadisper$`F value`[1]
anosim.res[which(anosim.res$anosim.res =="jaccard.all"), "df.betadisper.1"] <-
  jaccard.betadisper$Df[1]
anosim.res[which(anosim.res$anosim.res =="jaccard.all"), "df.betadisper.2"] <-
  jaccard.betadisper$Df[2]
anosim.res[which(anosim.res$anosim.res =="jaccard.all"), "p.betadisper"] <-
  jaccard.betadisper$`Pr(>F)`[1]

#<< ANOSIM >>---------------------------------------------------------------
jaccard.anosim <- anosim(comm.matrix, grouping = jaccard.matrix$Burn_status,
                         distance = "jaccard")
jaccard.anosim

#--Add results to data frame
anosim.res[which(anosim.res$anosim.res =="jaccard.all"), "r"] <- jaccard.anosim$statistic
anosim.res[which(anosim.res$anosim.res =="jaccard.all"), "p"] <- jaccard.anosim$signif

#<< PERMANOVA >>--------------------------------------------------------------------------
#--adonis
jaccard.adonis <- adonis(formula = comm.dist.jaccard ~ Burn_status*Range,
                          data = jaccard.matrix, permutations = 1000)
jaccard.adonis

#--add results to data.frame
#--Burn status f.model, r2, p-value
permanova.res[which(permanova.res$test == "jaccard"), "F.model.burn"] <-
  jaccard.adonis$aov.tab$F.Model[1]
permanova.res[which(permanova.res$test == "jaccard"), "r2.burn"] <-
  jaccard.adonis$aov.tab$R2[1]
permanova.res[which(permanova.res$test == "jaccard"), "p.burn"] <-
  jaccard.adonis$aov.tab$`Pr(>F)`[1]

#--Range f.model, r2, p-value
permanova.res[which(permanova.res$test == "jaccard"), "F.model.range"] <-
  jaccard.adonis$aov.tab$F.Model[2]
permanova.res[which(permanova.res$test == "jaccard"), "r2.range"] <-
  jaccard.adonis$aov.tab$R2[2]
permanova.res[which(permanova.res$test == "jaccard"), "p.range"] <-
  jaccard.adonis$aov.tab$`Pr(>F)`[2]

#--Burn/range f.model, r2, p-value
permanova.res[which(permanova.res$test == "jaccard"), "F.model.burn/range"] <- 
  jaccard.adonis$aov.tab$F.Model[3]
permanova.res[which(permanova.res$test == "jaccard"), "r2.burn/range"] <- 
  jaccard.adonis$aov.tab$R2[3]
permanova.res[which(permanova.res$test == "jaccard"), "p.burn/range"] <- 
  jaccard.adonis$aov.tab$`Pr(>F)`[3]

#========================================================================================#
# Morisita based dissimilarity index: tree level----
#========================================================================================#

#--isolate site X species matrix only without metadata
comm.matrix <- stsp.matrix[8:length(stsp.matrix)]

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
anosim.res[which(anosim.res$anosim.res =="morisita.all"), "stress.nmds"] <-
  morisita.otu$stress

jpeg(filename = paste0('/Volumes/Cenococcum/PhD/Dissertation/Chpt.2/Figures/',
                       'Fig3c_morisitahorn.jpeg'), width = 700, height = 600,
     quality = 100)
par(mfrow = c(1,1), "mar"=c(6, 5, 5, 3))

#--Plot NMDS of EM community based on Jaccard index and OTU abundance
plot(morisita.otu, display = "sites", type = "n", cex.lab = 2,
     cex.axis = 1.5, xlab = 'Axis 1',ylab = 'Axis 2')

# colors for points
color.vec <- data.frame(shape = rep(NA, nrow(morisita.matrix)),
                        p.group = morisita.matrix$Range)
color.vec[color.vec$p.group == 'santa.catalina', 'shape'] <- 16
color.vec[color.vec$p.group == 'pinaleno', 'shape'] <- 15

#ordipointlabel(jaccard.otu, display = "sites")
#--isolate points for pinaleno mts.
pinaleno.burn <- row.names(morisita.matrix[morisita.matrix$Range == 'pinaleno' & 
                                             morisita.matrix$Burn_status == 'burned',])
pinaleno.unburn <- row.names(morisita.matrix[morisita.matrix$Range == 'pinaleno' & 
                                              morisita.matrix$Burn_status == 'unburned',])
santa.catalina.burn <- row.names(morisita.matrix[morisita.matrix$Range == 'santa.catalina' & 
                                                  morisita.matrix$Burn_status == 'burned',])
santa.catalina.unburn <- row.names(morisita.matrix[morisita.matrix$Range == 'santa.catalina' & 
                                                    morisita.matrix$Burn_status == 'unburned',])

#--isolate points using rows isolated above
points(morisita.otu$points[pinaleno.burn,1:2], display = "sites", cex = 3,
       pch = color.vec[color.vec$p.group == 'pinaleno' | 
                         color.vec$p.group == 'pinaleno',
                       'shape'],
       col = 'black',
       bg = 'black')
points(morisita.otu$points[pinaleno.unburn,1:2], display = "sites", cex = 3,
       pch = color.vec[color.vec$p.group == 'pinaleno' | 
                         color.vec$p.group == 'pinaleno',
                       'shape'],
       col = 'darkgreen',
       bg = 'darkgreen')
points(morisita.otu$points[santa.catalina.burn,1:2], display = "sites", cex = 3,
       pch = color.vec[color.vec$p.group == 'santa.catalina' | 
                         color.vec$p.group == 'santa.catalina',
                       'shape'],
       col = 'black',
       bg = 'black')
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

legend("bottomleft", legend = c('Pinaleno Mts.',
                                'Santa Catalina Mts.'), bty = "n",
       col = c('black', 'black'),
       pch = c(15,16),
       pt.bg =  c('black','black'), cex = 2)

legend("bottomright", legend = c('Burned',
                                 'Unburned'), bty = "n",
       col = c('black', 'darkgreen'),
       pch = c(16,16),
       pt.bg =  c('black', 'darkgreen'), cex = 2)

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
#--Overlays
fit <- envfit(morisita.otu ~ soil.pca + Tmax + prec, morisita.matrix)
fit
#plot(fit, type = 'n')

dev.off()

#--BetaDisper
#--a multivariate analogue of Levene's test for homogeneity of variance
betadisper <- betadisper(comm.dist.morisita, group = morisita.matrix$Burn_status)
morisita.betadisper <- anova(betadisper)
morisita.betadisper

#--add Betadisper results to table
anosim.res[which(anosim.res$anosim.res =="morisita.all"), "F.betadisper"] <-
  morisita.betadisper$`F value`[1]
anosim.res[which(anosim.res$anosim.res =="morisita.all"), "df.betadisper.1"] <-
  morisita.betadisper$Df[1]
anosim.res[which(anosim.res$anosim.res =="morisita.all"), "df.betadisper.2"] <-
  morisita.betadisper$Df[2]
anosim.res[which(anosim.res$anosim.res =="morisita.all"), "p.betadisper"] <-
  morisita.betadisper$`Pr(>F)`[1]

#<< PERMANOVA >>------------------------------------------------------------------------
morisita.adonis.anosim <- adonis(comm.dist.morisita ~ Burn_status,
                                 data = morisita.matrix, permutations = 1000)
morisita.adonis.anosim

#--Burn f.model, r2, p-value
permanova.res[which(permanova.res$test == "morisita"), "F.model.burn"] <-
  morisita.adonis.anosim$aov.tab$F.Model[1]
permanova.res[which(permanova.res$test == "morisita"), "r2.burn"] <-
  morisita.adonis.anosim$aov.tab$R2[1]
permanova.res[which(permanova.res$test == "morisita"), "p.burn"] <-
  morisita.adonis.anosim$aov.tab$`Pr(>F)`[1]

#<< PERMANOVA >>--------------------------------------------------------------------------
#--adonis
morisita.adonis <- adonis(formula = comm.dist.morisita ~ Burn_status*Range,
                          data = morisita.matrix, permutations = 1000,
                          strata = morisita.matrix$Range)
morisita.adonis

#--add results to data.frame
#--Burn status f.model, r2, p-value
permanova.res[which(permanova.res$test == "morisita"), "F.model.burn"] <-
  morisita.adonis$aov.tab$F.Model[1]
permanova.res[which(permanova.res$test == "morisita"), "r2.burn"] <-
  morisita.adonis$aov.tab$R2[1]
permanova.res[which(permanova.res$test == "morisita"), "p.burn"] <-
  morisita.adonis$aov.tab$`Pr(>F)`[1]

#--Range f.model, r2, p-value
permanova.res[which(permanova.res$test == "morisita"), "F.model.range"] <-
 morisita.adonis$aov.tab$F.Model[2]
permanova.res[which(permanova.res$test == "morisita"), "r2.range"] <-
 morisita.adonis$aov.tab$R2[2]
permanova.res[which(permanova.res$test == "morisita"), "p.range"] <-
 morisita.adonis$aov.tab$`Pr(>F)`[2]

#--Burn/range f.model, r2, p-value
permanova.res[which(permanova.res$test == "morisita"), "F.model.burn/range"] <- 
  morisita.adonis$aov.tab$F.Model[3]
permanova.res[which(permanova.res$test == "morisita"), "r2.burn/range"] <- 
  morisita.adonis$aov.tab$R2[3]
permanova.res[which(permanova.res$test == "morisita"), "p.burn/range"] <- 
  morisita.adonis$aov.tab$`Pr(>F)`[3]


#--output permanova results
write.csv(permanova.res, paste0(res.dir, "Permanova_range_SequenceBased.csv"),
          row.names = F)

#========================================================================================#
# Jaccard based dissimilarity index: Pinaleno----
#========================================================================================#

#--isolate otu data
comm.matrix <- pinaleno.matrix[8:length(pinaleno.matrix)]

#--comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 4]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 0, ] # remove rows with sums of 0
jaccard.matrix <- pinaleno.matrix[row.names(comm.matrix),]

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
jpeg(filename = paste0('/Volumes/Cenococcum/PhD/Dissertation/Chpt.2/Figures/',
                       'Fig3a_jaccard.jpeg'), width = 700, height = 600,
     quality = 100)
par(mfrow = c(1,1), "mar"=c(6, 5, 5, 3))

#--Plot NMDS of EM community based on Jaccard index and OTU abundance
plot(jaccard.otu, display = "sites", type = "n", cex.lab = 2,
     cex.axis = 1.5, xlab = 'Axis 1', ylab = 'Axis 2')
# colors for points
color.vec <- data.frame(color = rep(NA, nrow(jaccard.matrix)),
                        p.group = jaccard.matrix$Burn_status)
color.vec[color.vec$p.group == 'burned', 'color'] <- 'black'
color.vec[color.vec$p.group == 'unburned', 'color'] <- 'darkgreen'

color.vec[color.vec$p.group == 'burned', 'shape'] <- 15
color.vec[color.vec$p.group == 'unburned', 'shape'] <- 15

#ordipointlabel(jaccard.otu, display = "sites")
#--isolate points for pinaleno mts.
burned <- row.names(jaccard.matrix[jaccard.matrix$Burn_status == 'burned',])
unburned <- row.names(jaccard.matrix[jaccard.matrix$Burn_status == 'unburned',])

#--isolate points using rows isolated above
points(jaccard.otu$points[burned,1:2], display = "sites", cex = 3,
       pch = color.vec[color.vec$p.group == 'burned',
                       'shape'],
       col = color.vec[color.vec$p.group == 'burned',
                       'color'],
       bg = color.vec[color.vec$p.group == 'burned',
                      'color'])
points(jaccard.otu$points[unburned,1:2], display = "sites", cex = 3,
       pch = color.vec[color.vec$p.group == 'unburned',
                       'shape'],
       col = color.vec[color.vec$p.group == 'unburned',
                       'color'],
       bg = color.vec[color.vec$p.group == 'unburned',
                      'color'])
legend("topright", legend = c('Burned',
                             'Unburned'), bty = "n",
       col = c('black', 'darkgreen'),
       pch = c(15,15),
       pt.bg =  c('black','darkgreen'), cex = 2)

#--Ordihull variations by range, burn, and both burn and range
#ordihull(jaccard.otu, groups = color.vec$p.group)
ordiellipse(jaccard.otu, groups = jaccard.matrix$Burn_status,
            col = c('black', 'black'),
            kind = 'ehull')
#ordihull(jaccard.otu, groups = stsp.matrix$Range) # just mt. range
#--Overlays
# fit <- envfit(jaccard.otu ~ po4.p.ppm + temp.warmest.quarter +
#                 pH.su, jaccard.matrix)
# fit
# plot(fit, type = 'n')

dev.off()

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
comm.matrix <- pinaleno.matrix[8:length(pinaleno.matrix)]

#--comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 4]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 0, ] # remove rows with sums of 0
morisita.matrix <- pinaleno.matrix[row.names(comm.matrix),]

#--distance matrix using jaccard index
comm.dist.morisita <- vegdist(comm.matrix, method = "horn", binary = F)

#--NMDS analysis
morisita.otu <- metaMDS(comm.dist.morisita, dist = "bray", permutations = 999,
                        try = 100, trymax = 1000)

#--add stress of NMDS to results table
anosim.res[which(anosim.res$anosim.res =="morisita.p"), "stress.nmds"] <-
  morisita.otu$stress

#--format and output NMDS plots to figure folder
jpeg(filename = paste0('/Volumes/Cenococcum/PhD/Dissertation/Chpt.2/Figures/',
                       'Fig3a-morisitahorn.jpeg'), width = 700, height = 600,
     quality = 100)
par(mfrow = c(1,1), "mar"=c(6, 5, 5, 3))

#--Plot NMDS of EM community based on Jaccard index and OTU abundance
plot(morisita.otu, display = "sites", type = "n", cex.lab = 2,
     cex.axis = 1.5, xlab = 'Axis 1', ylab = 'Axis 2')

# colors for points
color.vec <- data.frame(color = rep(NA, nrow(morisita.matrix)),
                        p.group = morisita.matrix$Burn_status)
color.vec[color.vec$p.group == 'burned', 'color'] <- 'black'
color.vec[color.vec$p.group == 'unburned', 'color'] <- 'darkgreen'

color.vec[color.vec$p.group == 'burned', 'shape'] <- 15
color.vec[color.vec$p.group == 'unburned', 'shape'] <- 15

#ordipointlabel(jaccard.otu, display = "sites")
#--isolate points for pinaleno mts.
burned <- row.names(morisita.matrix[jaccard.matrix$Burn_status == 'burned',])
unburned <- row.names(morisita.matrix[jaccard.matrix$Burn_status == 'unburned',])

#--isolate points using rows isolated above
points(morisita.otu$points[burned,1:2], display = "sites", cex = 3,
       pch = color.vec[color.vec$p.group == 'burned',
                       'shape'],
       col = color.vec[color.vec$p.group == 'burned',
                       'color'],
       bg = color.vec[color.vec$p.group == 'burned',
                      'color'])
points(morisita.otu$points[unburned,1:2], display = "sites", cex = 3,
       pch = color.vec[color.vec$p.group == 'unburned',
                       'shape'],
       col = color.vec[color.vec$p.group == 'unburned',
                       'color'],
       bg = color.vec[color.vec$p.group == 'unburned',
                      'color'])
legend("bottomleft", legend = c('Burned',
                             'Unburned'), bty = "n",
       col = c('black', 'darkgreen'),
       pch = c(16,16),
       pt.bg =  c('black','darkgreen'), cex = 2)

#--Ordihull variations by range, burn, and both burn and range
ordiellipse(morisita.otu, groups = morisita.matrix$Burn_status,
            col = c('black', 'black'),
            kind = 'ehull')

dev.off()

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
comm.matrix <- santa.matrix[8:length(santa.matrix)]

#--comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 4]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 0, ] # remove rows with sums of 0
jaccard.matrix <- santa.matrix[row.names(comm.matrix),]

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
jpeg(filename = paste0('/Volumes/Cenococcum/PhD/Dissertation/Chpt.2/Figures/',
                       'Fig3b_jaccard.jpeg'), width = 700, height = 600,
     quality = 100)
par(mfrow = c(1,1), "mar"=c(6, 5, 5, 3))

#--Plot NMDS of EM community based on Jaccard index and OTU abundance
plot(jaccard.otu, display = "sites", type = "n", cex.lab = 2,
     cex.axis = 1.5, xlab = 'Axis 1', ylab = 'Axis 2')

# colors for points
color.vec <- data.frame(color = rep(NA, nrow(jaccard.matrix)),
                        p.group = jaccard.matrix$Burn_status)
color.vec[color.vec$p.group == 'burned', 'color'] <- 'black'
color.vec[color.vec$p.group == 'unburned', 'color'] <- 'darkgreen'

color.vec[color.vec$p.group == 'burned', 'shape'] <- 16
color.vec[color.vec$p.group == 'unburned', 'shape'] <- 16

#ordipointlabel(jaccard.otu, display = "sites")
#--isolate points for pinaleno mts.
burned <- row.names(jaccard.matrix[jaccard.matrix$Burn_status == 'burned',])
unburned <- row.names(jaccard.matrix[!jaccard.matrix$Burn_status == 'burned',])

#--isolate points using rows isolated above
points(jaccard.otu$points[burned,1:2], display = "sites", cex = 3,
       pch = color.vec[color.vec$p.group == 'burned',
                       'shape'],
       col = color.vec[color.vec$p.group == 'burned',
                       'color'],
       bg = color.vec[color.vec$p.group == 'burned',
                      'color'])
points(jaccard.otu$points[unburned,1:2], display = "sites", cex = 3,
       pch = color.vec[color.vec$p.group == 'unburned',
                       'shape'],
       col = color.vec[color.vec$p.group == 'unburned',
                       'color'],
       bg = color.vec[color.vec$p.group == 'unburned',
                      'color'])
legend("topright", legend = c('Burned','Unburned'), bty = "n",
       col = c('black','darkgreen'),
       pch = c(16,16),
       pt.bg =  c('black','darkgreen'), cex = 2)

#--Ordihull variations by range, burn, and both burn and range
#ordihull(jaccard.otu, groups = color.vec$p.group)
ordiellipse(jaccard.otu, groups = jaccard.matrix$Burn_status,
            col = c('black','black'),
            kind = 'ehull')
#ordihull(jaccard.otu, groups = stsp.matrix$Range) # just mt. range
#--Overlays
# fit <- envfit(jaccard.otu ~ po4.p.ppm + temp.warmest.quarter +
#                 pH.su, jaccard.matrix)
# fit
# plot(fit, type = 'n')

dev.off()

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
comm.matrix <- santa.matrix[8:length(santa.matrix)]

#--comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 4]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 0, ] # remove rows with sums of 0
morisita.matrix <- santa.matrix[row.names(comm.matrix),]

#--distance matrix using jaccard index
comm.dist.morisita <- vegdist(comm.matrix, method = "horn", binary = F)

#--NMDS analysis
morisita.otu <- metaMDS(comm.dist.morisita, dist = "bray", permutations = 999,
                        try = 100, trymax = 1000)

#--add stress of NMDS to results table
anosim.res[which(anosim.res$anosim.res =="morisita.scm"), "stress.nmds"] <-
  morisita.otu$stress

#--format and output NMDS plots to figure folder
jpeg(filename = paste0('/Volumes/Cenococcum/PhD/Dissertation/Chpt.2/Figures/',
                       'Fig3b_morisitahorn.jpeg'), width = 700, height = 600,
     quality = 100)
par(mfrow = c(1,1), "mar"=c(6, 5, 5, 3))

#--Plot NMDS of EM community based on Jaccard index and OTU abundance
plot(morisita.otu, display = "sites", type = "n", cex.lab = 2,
     cex.axis = 1.5, xlab = 'Axis 1', ylab = 'Axis 2')
# colors for points
color.vec <- data.frame (color = rep (NA, nrow(morisita.matrix)),
                         p.group = morisita.matrix$Burn_status)
color.vec[color.vec$p.group == 'burned', 'color'] <- 'black'
color.vec[color.vec$p.group == 'unburned', 'color'] <- 'darkgreen'
color.vec[color.vec$p.group == 'burned', 'shape'] <- 16
color.vec[color.vec$p.group == 'unburned', 'shape'] <- 16

#ordipointlabel(morisita.otu, display = "sites")

#--isolate points for pinaleno mts.
burned <- row.names(morisita.matrix[morisita.matrix$Burn_status == 'burned',])
unburned <- row.names(morisita.matrix[!morisita.matrix$Burn_status == 'burned',])

#--isolate points using rows isolated above
points(morisita.otu$points[burned,1:2], display = "sites", cex = 3,
       pch = color.vec[color.vec$p.group == 'burned',
                       'shape'],
       col = color.vec[color.vec$p.group == 'burned',
                       'color'],
       bg = color.vec[color.vec$p.group == 'burned',
                      'color'])
points(morisita.otu$points[unburned,1:2], display = "sites", cex = 3,
       pch = color.vec[color.vec$p.group == 'unburned',
                       'shape'],
       col = color.vec[color.vec$p.group == 'unburned',
                       'color'],
       bg = color.vec[color.vec$p.group == 'unburned',
                      'color'])
legend("bottomright", legend = c('Burned','Unburned'), bty = "n",
       col = c('black','darkgreen'),
       pch = c(16,16),
       pt.bg =  c('black','darkgreen'), cex = 2)

#--Ordihull variations by range, burn, and both burn and range
#ordihull(morisita.otu, groups = color.vec$p.group,
#         col = c('black','green4','black','green4')),
ordiellipse(morisita.otu, groups = morisita.matrix$Burn_status,
            col = c('black','black'),
            kind = 'ehull') # just burn status
#ordihull(morisita.otu, groups = stsp.matrix$Range,
#         col = c('black','green4','black','green4')) # just mt. range

dev.off()

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

#<< ANOSIM >>-----------------------------------------------------------------------------
morisita.anosim <- anosim(comm.matrix, grouping = morisita.matrix$Burn_status,
                          distance = "horn")
morisita.anosim

#--Add results to data frame
anosim.res[which(anosim.res$anosim.res =="morisita.scm"), "r"] <- morisita.anosim$statistic
anosim.res[which(anosim.res$anosim.res =="morisita.scm"), "p"] <- morisita.anosim$signif

#--export anosim table
write.csv(anosim.res, paste0(res.dir, "ANOSIM_SequenceBased.csv"),
          row.names = F)

#========================================================================================#
# Venn diagram of PERMANOVA restuls----
#========================================================================================#

#install.packages('VennDiagram')
library(VennDiagram)

#<< Jaccard >>-----------------------------------------------------------------------------
#fire history 0.04410 r-squared
#range 0.07065 r-squared
#fire history * range 0.09072 r-squared

jpeg(filename = paste0('/Volumes/Cenococcum/PhD/Dissertation/Chpt.2/Figures/',
                       'Fig3d_jaccard.jpeg'), width = 700, height = 600,
     quality = 100)

grid.newpage()
draw.pairwise.venn(area1 = 13.482,
                   area2 = 16.137,
                   cross.area = 9.072,
                   category = c("Fire history","Range"),
                   lty = rep("blank", 2), col = 'black',fill = c("grey", "light blue"), 
                   alpha = rep(0.5, 2), cat.pos = c(200, 160), cat.dist = 0.05,
                   cex = 1.5, cat.cex = 2)

dev.off()

#<< Morisita-horn >>--------------------------------------------------------------------
#fire history 0.06127 r-squared
#range 0.06905 r-squared
#fire history * range 0.08595 r-squared

jpeg(filename = paste0('/Volumes/Cenococcum/PhD/Dissertation/Chpt.2/Figures/',
                       'Fig3d_morisitahorn.jpeg'), width = 700, height = 600,
     quality = 100)

grid.newpage()
draw.pairwise.venn(area1 = 14.722,
                   area2 = 15.500,
                   cross.area = 8.595,
                   category = c("Fire history","Range"),
                   lty = rep("blank", 2), col = 'black',fill = c("grey", "light blue"), 
                   alpha = rep(0.5, 2), cat.pos = c(200, 160), cat.dist = 0.05,
                   cex = 1.5, cat.cex = 2)

dev.off()