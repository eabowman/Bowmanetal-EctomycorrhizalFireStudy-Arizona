## Script created by Liz Bowman June 3, 2017
## for analyzing community differences between ranges and burn status of sites

#========================================================================================#
# Load data ----
#========================================================================================#

#----------------------------------------------------------------------------------------#
# Load data and clean up: Tree level----
#----------------------------------------------------------------------------------------#
stsp.matrix <- read.csv(paste0(dat.dir,'97%_SitexSpecies_TipAb.csv'), as.is = T)

# << Create distinct tables for each burn history/range >> --
FA.matrix <- stsp.matrix[stsp.matrix$Burn_status == 'burned',]
FU.matrix <- stsp.matrix[stsp.matrix$Burn_status == 'unburned',]

scm.matrix <- stsp.matrix[stsp.matrix$Range == 'santa.catalina',]
pm.matrix <- stsp.matrix[stsp.matrix$Range == 'pinaleno',]

#----------------------------------------------------------------------------------------#
# Create result tables
#----------------------------------------------------------------------------------------#
#--Table of Anosim results
anosim.res <- c('jaccard.overall','morisita.overall',
                'jaccard.p', 'morisita.p',
                'jaccard.scm', 'morisita.scm')
anosim.res <- data.frame(anosim.res)
#--Table of PERMANOVA results
permanova.res <- data.frame(c('jaccard','morisita'))
colnames(permanova.res) <- "test"

#========================================================================================#
# Jaccard based dissimilarity index: Overall----
#========================================================================================#

#--remove outliers
jaccard.matrix <- stsp.matrix[!stsp.matrix$Tree %in% c('NF19','NF16'),]
#jaccard.matrix <- stsp.matrix

#--isolate otu data
comm.matrix <- jaccard.matrix[12:length(jaccard.matrix)]

#--remove singletons; comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 4]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 1,]
jaccard.matrix <- jaccard.matrix[rownames(comm.matrix),]

#--distance matrix using jaccard index
comm.dist.jaccard <- vegdist(comm.matrix, method = 'jaccard', binary = F, na.rm = T)

#--NMDS analysis
jaccard.otu <- metaMDS(comm.dist.jaccard, dist = 'bray', permutations = 999,
                       try = 100)
jaccard.otu$stress

#--BetaDisper
#--a multivariate analogue of Levene's test for homogeneity of variance
betadisper <- betadisper(comm.dist.jaccard, group = jaccard.matrix$Burn_status)
jaccard.betadisper <- anova(betadisper)
jaccard.betadisper

#--add Betadisper results to table
anosim.res[which(anosim.res$anosim.res =="jaccard.overall"), "F.betadisper"] <-
  jaccard.betadisper$`F value`[1]
anosim.res[which(anosim.res$anosim.res =="jaccard.overall"), "df.betadisper.1"] <-
  jaccard.betadisper$Df[1]
anosim.res[which(anosim.res$anosim.res =="jaccard.overall"), "df.betadisper.2"] <-
  jaccard.betadisper$Df[2]
anosim.res[which(anosim.res$anosim.res =="jaccard.overall"), "p.betadisper"] <-
  jaccard.betadisper$`Pr(>F)`[1]

#<< PERMANOVA >>-----------------------------------------------------------------------------
jaccard.adonis <- adonis(comm.dist.jaccard ~ Burn_status * Range, data = jaccard.matrix)
jaccard.adonis

#--Burn f.model, r2, p-value
anosim.res[which(anosim.res$anosim.res == "jaccard.overall"), "F.model.burn"] <-
  morisita.adonis.anosim$aov.tab$F.Model[1]
anosim.res[which(anosim.res$anosim.res == "jaccard.overall"), "r2.burn"] <-
  morisita.adonis.anosim$aov.tab$R2[1]
anosim.res[which(anosim.res$anosim.res == "jaccard.overall"), "p.burn"] <-
  morisita.adonis.anosim$aov.tab$`Pr(>F)`[1]

#--Base R plot
#--format and output NMDS plots to figure folder
jpeg(filename = 'figures/NMDS_overall_Jaccard.jpeg',
     width = 700, height = 600,
     quality = 100)
par(mfrow = c(1,1), "mar"=c(6, 5, 5, 3))
#par(mar=c(7.1, 4.1, 4.1, 12.1), xpd=TRUE)

#--Plot NMDS of EM community based on Jaccard index and OTU abundance
plot(jaccard.otu, display = "sites", type = "n", cex.lab = 2,
     cex.axis = 1.5, xlab = 'Axis 1', ylab = 'Axis 2')
# color and shape for points
color.vec <- data.frame(color = rep(NA, nrow(jaccard.matrix)),
                        shape = rep(NA, nrow(jaccard.matrix)),
                        fire.group = jaccard.matrix$Burn_status,
                        range.group = jaccard.matrix$Range)
# populate color
color.vec[color.vec$fire.group == 'burned' & color.vec$range.group == 'pinaleno',
          'color'] <- 'black'
color.vec[color.vec$fire.group == 'unburned' & color.vec$range.group == 'pinaleno',
          'color'] <- 'darkgreen'
color.vec[color.vec$fire.group == 'burned' & color.vec$range.group == 'santa.catalina',
          'color'] <- 'darkgrey'
color.vec[color.vec$fire.group == 'unburned' & color.vec$range.group == 'santa.catalina',
          'color'] <- 'lightgreen'
# populate shape
color.vec[color.vec$range.group == 'pinaleno', 'shape'] <- 19
color.vec[color.vec$range.group == 'santa.catalina', 'shape'] <- 18

#ordipointlabel(jaccard.otu, display = "sites")

#--add points to plot
points(jaccard.otu$points, display = "sites", cex = 3,
       pch = color.vec$shape,
       col = color.vec$color,
       bg = color.vec$color)

#--Ordihull variations by range, burn, and both burn and range
ordihull(jaccard.otu, groups = paste0(color.vec$fire.group, color.vec$range.group),
         draw = 'polygon')
dev.off()

#--Legend
jpeg(filename = 'figures/Legend.jpeg',
     width = 700, height = 600,
     quality = 100)
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend('topleft', 
       legend=c('Santa Catalina Mts. burned',
                'Santa Catalina Mts. unburned',
                'Pinaleno Mts. burned',
                'Pinaleno Mts. unburned'),
       col = c('darkgrey','lightgreen','black','darkgreen'),
       pch = c(18,18,19,19),
       pt.cex = 3, cex = 1.5, bty = 'n')
dev.off()

#<< Pairwise >>-----------------------------------------------------------------------------
jac.overall <- read.csv('data/JaccardOverall_firehistory.csv', as.is = T)

# rownames(stsp.matrix) <- stsp.matrix$Tree
# #--remove outliers
# jaccard.matrix <- stsp.matrix[!stsp.matrix$Tree %in% c('NF19','NF16'),]
# 
# #--isolate otu data
# comm.matrix <- jaccard.matrix[12:length(jaccard.matrix)]
# 
# #--remove singletons; comment to include singletons
# comm.matrix <- comm.matrix[colSums(comm.matrix) >= 4]
# comm.matrix <- comm.matrix[rowSums(comm.matrix) > 1,]
# jaccard.matrix <- jaccard.matrix[rownames(comm.matrix),]
# 
# #--distance matrix using jaccard index
# comm.dist.jaccard <- vegdist(comm.matrix, method = 'jaccard', binary = F, na.rm = T)
# comm.dist.jaccard <- as.data.frame(as.matrix(comm.dist.jaccard))
# comm.dist.jaccard <- read.csv('data_output/JaccardOverall_firehistory_dist.csv', row.names = 1)
# comm.dist.jaccard$tree1 <- rownames(comm.dist.jaccard)
# 
# #--Make dataframe
# jac.overall <- gather(comm.dist.jaccard,'tree2', 'dissimilarity', -tree1)
# jac.overall <- jac.overall[!is.na(jac.overall$dissimilarity),]
# #--Add fire history data
# for(i in jac.overall$tree1){
#   for(t in jac.overall$tree2){
#     jac.overall[jac.overall$tree1 == i, 'fire.history.1'] <- jaccard.matrix[jaccard.matrix$Tree == i, 'Burn_status']
#     jac.overall[jac.overall$tree2 == t, 'fire.history.2'] <- jaccard.matrix[jaccard.matrix$Tree == t, 'Burn_status']
#     jac.overall[jac.overall$tree1 == i, 'range.1'] <- jaccard.matrix[jaccard.matrix$Tree == i, 'Range']
#     jac.overall[jac.overall$tree2 == t, 'range.2'] <- jaccard.matrix[jaccard.matrix$Tree == t, 'Range']
#   }
# }
# #--assess if within same or different fire history
# for(i in 1:nrow(jac.overall)){
#   if(jac.overall[i, 'fire.history.1'] == jac.overall[i, 'fire.history.2']){
#     jac.overall[i, 'comp'] <- 'Same'
#   } else{jac.overall[i, 'comp'] <- 'Different'}
# }

# #--assess if within same or different range origin
# for(i in 1:nrow(jac.overall)){
#   if(jac.overall[i, 'range.1'] == jac.overall[i, 'range.2']){
#     jac.overall[i, 'comp.range'] <- 'Same'
#   } else{jac.overall[i, 'comp.range'] <- 'Different'}
# }

#--Remove self-comparisons (dissimilarity = 0) and outliers
jac.overall <- jac.overall[jac.overall$dissimilarity > 0, ]
# same.mor.scm <- mor.scm[mor.scm$comp == 'Same',] 
# same.mor.scm <- same.mor.scm[same.mor.scm$dissimilarity > 0.4,]
# diff.mor.scm <- mor.scm[mor.scm$comp == 'Different',]
# diff.mor.scm <- diff.mor.scm[diff.mor.scm$dissimilarity > 0.57,]
# mor.scm.out <- bind_rows(same.mor.scm, diff.mor.scm)

#--Logit transform data
jac.overall$logit.dis <- logit(jac.overall$dissimilarity)

#--wilcox test
wilcox.test(dissimilarity ~ comp.fire, data = jac.overall)

ggplot(jac.overall, mapping = aes(x = comp.fire,
                                  y = logit.dis)) +
  geom_boxplot() +
  xlab('Fire history') +
  ylab('Logit Jaccard dissimilarity') +
  theme_classic()

#========================================================================================#
# Morisita based dissimilarity index: Overall----
#========================================================================================#

#--remove outliers
morisita.matrix <- stsp.matrix[!stsp.matrix$Tree %in% c('NF16','NF19'),]
#morisita.matrix <- stsp.matrix

#--isolate otu data
comm.matrix <- morisita.matrix[12:length(morisita.matrix)]

#--remove singletons; comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 2]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 1,]
morisita.matrix <- morisita.matrix[rownames(comm.matrix),]

#--distance matrix using jaccard index
comm.dist.morisita <- vegdist(comm.matrix, method = 'horn', binary = F, na.rm = T)

#--NMDS analysis
morisita.otu <- metaMDS(comm.dist.morisita, dist = 'bray', permutations = 999,
                       try = 100)
morisita.otu$stress

#--BetaDisper
#--a multivariate analogue of Levene's test for homogeneity of variance
betadisper <- betadisper(comm.dist.morisita, group = morisita.matrix$Burn_status)
jaccard.betadisper <- anova(betadisper)
jaccard.betadisper

#--add Betadisper results to table
anosim.res[which(anosim.res$anosim.res =="jaccard.overall"), "F.betadisper"] <-
  jaccard.betadisper$`F value`[1]
anosim.res[which(anosim.res$anosim.res =="jaccard.overall"), "df.betadisper.1"] <-
  jaccard.betadisper$Df[1]
anosim.res[which(anosim.res$anosim.res =="jaccard.overall"), "df.betadisper.2"] <-
  jaccard.betadisper$Df[2]
anosim.res[which(anosim.res$anosim.res =="jaccard.overall"), "p.betadisper"] <-
  jaccard.betadisper$`Pr(>F)`[1]

#<< PERMANOVA >>-----------------------------------------------------------------------------
morisita.adonis <- adonis(comm.dist.morisita ~ Burn_status * Range, data = morisita.matrix)
morisita.adonis

#--Burn f.model, r2, p-value
anosim.res[which(anosim.res$anosim.res == "morisita.overall"), "F.model.burn"] <-
  morisita.adonis.anosim$aov.tab$F.Model[1]
anosim.res[which(anosim.res$anosim.res == "morisita.overall"), "r2.burn"] <-
  morisita.adonis.anosim$aov.tab$R2[1]
anosim.res[which(anosim.res$anosim.res == "morisita.overall"), "p.burn"] <-
  morisita.adonis.anosim$aov.tab$`Pr(>F)`[1]


#--Base R plot
#--format and output NMDS plots to figure folder
jpeg(filename = 'figures/NMDS_overall_MorisitaHorn.jpeg',
     width = 700, height = 600,
     quality = 100)
par(mfrow = c(1,1), "mar"=c(6, 5, 5, 3))
#par(mar=c(7.1, 4.1, 4.1, 12.1), xpd=TRUE)

#--Plot NMDS of EM community based on Jaccard index and OTU abundance
plot(morisita.otu, display = "sites", type = "n", cex.lab = 2,
     cex.axis = 1.5, xlab = 'Axis 1', ylab = 'Axis 2')
# color and shape for points
color.vec <- data.frame(color = rep(NA, nrow(morisita.matrix)),
                        shape = rep(NA, nrow(morisita.matrix)),
                        fire.group = morisita.matrix$Burn_status,
                        range.group = morisita.matrix$Range)
# populate color
color.vec[color.vec$fire.group == 'burned' & color.vec$range.group == 'pinaleno',
          'color'] <- 'black'
color.vec[color.vec$fire.group == 'unburned' & color.vec$range.group == 'pinaleno',
          'color'] <- 'darkgreen'
color.vec[color.vec$fire.group == 'burned' & color.vec$range.group == 'santa.catalina',
          'color'] <- 'darkgrey'
color.vec[color.vec$fire.group == 'unburned' & color.vec$range.group == 'santa.catalina',
          'color'] <- 'lightgreen'
# populate shape
color.vec[color.vec$fire.group == 'burned', 'shape'] <- 19
color.vec[color.vec$fire.group == 'unburned', 'shape'] <- 18

#ordipointlabel(morisita.otu, display = "sites")

#--add points to plot
points(morisita.otu$points, display = "sites", cex = 3,
       pch = color.vec$shape,
       col = color.vec$color,
       bg = color.vec$color)

# #--Legend
# legend("topright", inset=c(-1,0), 
#        legend=c('Santa Catalina Mts. burned',
#                 'Santa Catalina Mts. unburned',
#                 'Pinaleno Mts. burned',
#                 'Pinaleno Mts. unburned'),
#        col = c('darkgrey','lightgreen','black','darkgreen'),
#        pch = c(19,18,19,18),
#        cex = 1)

#--Ordihull variations by range, burn, and both burn and range
ordihull(morisita.otu, groups = paste0(color.vec$fire.group, color.vec$range.group),
         draw = 'polygon')
dev.off()

#<< Pairwise >>-----------------------------------------------------------------------------
horn.overall <- read.csv('data/MorisitaHornOverall_firehistory.csv', as.is = T)
# rownames(stsp.matrix) <- stsp.matrix$Tree
# #--remove outliers
# morisita.matrix <- stsp.matrix[!stsp.matrix$Tree %in% c('NF19','NF16'),]
# 
# #--isolate otu data
# comm.matrix <- morisita.matrix[12:length(morisita.matrix)]
# 
# #--remove singletons; comment to include singletons
# comm.matrix <- comm.matrix[colSums(comm.matrix) >= 4]
# comm.matrix <- comm.matrix[rowSums(comm.matrix) > 1,]
# morisita.matrix <- morisita.matrix[rownames(comm.matrix),]
# 
# #--distance matrix using jaccard index
# comm.dist.horn <- vegdist(comm.matrix, method = 'horn', binary = F, na.rm = T)
# comm.dist.horn <- as.data.frame(as.matrix(comm.dist.horn))
# comm.dist.horn <- read.csv('./data_output/MorisitaHornOverall_firehistory_dist.csv', row.names = 1)
# comm.dist.horn$tree1 <- rownames(comm.dist.horn)
# 
# #--Make dataframe
# horn.overall <- gather(comm.dist.horn,'tree2', 'dissimilarity', -tree1)
# horn.overall <- horn.overall[!is.na(horn.overall$dissimilarity),]
#--Add fire history data
# for(i in horn.overall$tree1){
#   for(t in horn.overall$tree2){
#     horn.overall[horn.overall$tree1 == i, 'fire.history.1'] <- morisita.matrix[morisita.matrix$Tree == i, 'Burn_status']
#     horn.overall[horn.overall$tree2 == t, 'fire.history.2'] <- morisita.matrix[morisita.matrix$Tree == t, 'Burn_status']
#     horn.overall[horn.overall$tree1 == i, 'range.1'] <- morisita.matrix[morisita.matrix$Tree == i, 'Range']
#     horn.overall[horn.overall$tree2 == t, 'range.2'] <- morisita.matrix[morisita.matrix$Tree == t, 'Range']
#   }
# }
# #--assess if within same or different fire history
# for(i in 1:nrow(horn.overall)){
#   if(horn.overall[i, 'fire.history.1'] == horn.overall[i, 'fire.history.2']){
#     horn.overall[i, 'comp'] <- 'Same'
#   } else{horn.overall[i, 'comp'] <- 'Different'}
# }

# #--assess if within same or different range origin
# for(i in 1:nrow(horn.overall)){
#   if(horn.overall[i, 'range.1'] == horn.overall[i, 'range.2']){
#     horn.overall[i, 'comp.range'] <- 'Same'
#   } else{horn.overall[i, 'comp.range'] <- 'Different'}
# }

#--logit transform
horn.overall$logit.dis <- logit(horn.overall$dissimilarity)
#--Remove self-comparisons (dissimilarity = 0) and outliers
horn.overall <- horn.overall[horn.overall$dissimilarity > 0, ]
horn.overall <- horn.overall[horn.overall$dissimilarity < 0.9997,]

#--wilcox test
t.test(logit.dis ~ comp.fire, data = horn.overall)

ggplot(horn.overall, mapping = aes(x = comp.fire,
                                  y = logit.dis)) +
  geom_boxplot() +
  xlab('Fire history') +
  ylab('Jaccard dissimilarity') +
  theme_classic()

#========================================================================================#
# Jaccard based dissimilarity index: Pinaleno----
#========================================================================================#

#--Remove outlier, NF16
pm.matrix <- pm.matrix[!pm.matrix$Tree %in% c('NF19','NF16'),]

#--isolate otu data
comm.matrix <- pm.matrix[12:length(pm.matrix)]

#--comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 2]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 1, ] # remove rows with sums of 0
jaccard.matrix <- pm.matrix[row.names(comm.matrix),]

#--distance matrix using jaccard index
comm.dist.jaccard <- vegdist(comm.matrix, method = "jaccard", binary = F, na.rm = T)

#--NMDS analysis
jaccard.otu <- metaMDS(comm.dist.jaccard, dist = "bray", permutations = 999,
                       try = 100, trymax = 1000)

#--Stress
jaccard.otu$stress

#--add stress of NMDS to results table
anosim.res[which(anosim.res$anosim.res =="jaccard.p"), "stress.nmds"] <-
  jaccard.otu$stress


#--Base R plot
#--format and output NMDS plots to figure folder
jpeg(filename = 'figures/NMDS_pinaleno_Jaccard.jpeg',
     width = 700, height = 600,
     quality = 100)
par(mfrow = c(1,1), "mar"=c(6, 5, 5, 3))

#--Plot NMDS of EM community based on Jaccard index and OTU abundance
plot(jaccard.otu, display = "sites", type = "n", cex.lab = 2,
     cex.axis = 1.5, xlab = 'Axis 1', ylab = 'Axis 2')
# colors for points
color.vec <- data.frame(color = rep(NA, nrow(jaccard.matrix)),
                        p.group = jaccard.matrix$Burn_status)
color.vec[color.vec$p.group == 'burned', 'color'] <- 'black'
color.vec[color.vec$p.group == 'unburned', 'color'] <- 'black'

color.vec[color.vec$p.group == 'burned', 'shape'] <- 15
color.vec[color.vec$p.group == 'unburned', 'shape'] <- 0

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
legend("bottomleft", legend = c('Burned',
                             'Unburned'), bty = "n",
       col = c('black', 'black'),
       pch = c(15,0),
       pt.bg =  c('black','black'), cex = 2)

#--Ordihull variations by range, burn, and both burn and range
#ordihull(jaccard.otu, groups = color.vec$p.group)
ordiellipse(jaccard.otu, groups = jaccard.matrix$Burn_status,
            col = c('black', 'black'),
            kind = 'ehull')
#ordihull(jaccard.otu, groups = stsp.matrix$Range) # just mt. range
#--Overlays
# fit <- envfit(jaccard.otu ~ prec, jaccard.matrix)
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

#<< Pairwise >>-----------------------------------------------------------------------------
#--Remove outlier, NF16
pm.matrix <- pm.matrix[!pm.matrix$Tree %in% c('NF19','NF16'),]
rownames(pm.matrix) <- pm.matrix$Tree

#--isolate otu data
comm.matrix <- pm.matrix[12:length(pm.matrix)]

#--comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 4]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 1, ] # remove rows with sums of 0
jaccard.matrix <- pm.matrix[row.names(comm.matrix),]

#--distance matrix using jaccard index
comm.dist.jaccard <- vegdist(comm.matrix, method = "jaccard", binary = F, upper = T)
comm.dist.jaccard <- as.data.frame(as.matrix(comm.dist.jaccard))
comm.dist.jaccard$tree1 <- rownames(comm.dist.jaccard)

#--Make dataframe
jac.pin <- gather(comm.dist.jaccard,'tree2', 'dissimilarity', -tree1)
test <- distinct(jac.pin)
#--Add fire history data
for(i in jac.pin$tree1){
  for(t in jac.pin$tree2){
    jac.pin[jac.pin$tree1 == i, 'fire.history.1'] <- jaccard.matrix[jaccard.matrix$Tree == i, 'Burn_status']
    jac.pin[jac.pin$tree2 == t, 'fire.history.2'] <- jaccard.matrix[jaccard.matrix$Tree == t, 'Burn_status']
    }
}
#--assess if within same or different fire history
for(i in 1:nrow(jac.pin)){
  if(jac.pin[i, 'fire.history.1'] == jac.pin[i, 'fire.history.2']){
    jac.pin[i, 'comp'] <- 'Same'
  } else{jac.pin[i, 'comp'] <- 'Different'}
}

#--Remove self-comparisons (dissimilarity = 0) and outliers
jac.pin <- jac.pin[jac.pin$dissimilarity > 0, ]
jac.pin.out <- jac.pin[!jac.pin$dissimilarity < 0.92,]
#jac.pin.out <- jac.pin[jac.pin$comp == 'Different' & jac.pin$dissimilarity < 0.96,]

#--wilcox test
wilcox.test(dissimilarity ~ comp, data = jac.pin.out)

ggplot(jac.pin.out, mapping = aes(x = comp,
                              y = dissimilarity)) +
  geom_boxplot() +
  xlab('Fire history') +
  ylab('Jaccard dissimilarity') +
  theme_classic()


#========================================================================================#
# Morisita based dissimilarity index: Pinaleno----
#========================================================================================#

#--isolate site X species matrix only without metadata
comm.matrix <- pm.matrix[12:length(pm.matrix)]

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
jpeg(filename = 'figures/NMDS_pinaleno_MorisitaHorn.jpeg',
     width = 700, height = 600,
     quality = 100)
par(mfrow = c(1,1), "mar"=c(6, 5, 5, 3))

#--Plot NMDS of EM community based on Jaccard index and OTU abundance
plot(morisita.otu, display = "sites", type = "n", cex.lab = 2,
     cex.axis = 1.5, xlab = 'Axis 1', ylab = 'Axis 2')

# colors for points
color.vec <- data.frame(color = rep(NA, nrow(morisita.matrix)),
                        p.group = morisita.matrix$Burn_status)
color.vec[color.vec$p.group == 'burned', 'color'] <- 'black'
color.vec[color.vec$p.group == 'unburned', 'color'] <- 'black'

color.vec[color.vec$p.group == 'burned', 'shape'] <- 15
color.vec[color.vec$p.group == 'unburned', 'shape'] <- 0

#ordipointlabel(jaccard.otu, display = "sites")
#--isolate points for pinaleno mts.
burned <- row.names(morisita.matrix[morisita.matrix$Burn_status == 'burned',])
unburned <- row.names(morisita.matrix[morisita.matrix$Burn_status == 'unburned',])

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
       col = c('black', 'black'),
       pch = c(15,0),
       pt.bg =  c('black','black'), cex = 2)

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

#<< Pairwise >>-----------------------------------------------------------------------------
#--Remove outlier, NF16
pm.matrix <- pm.matrix[!pm.matrix$Tree %in% c('NF19','NF16'),]
rownames(pm.matrix) <- pm.matrix$Tree

#--isolate otu data
comm.matrix <- pm.matrix[12:length(pm.matrix)]

#--comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 4]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 1, ] # remove rows with sums of 0
morisita.matrix <- pm.matrix[row.names(comm.matrix),]

#--distance matrix using jaccard index
comm.dist.morisita <- vegdist(comm.matrix, method = "horn", binary = F, upper = T)
comm.dist.morisita <- as.data.frame(as.matrix(comm.dist.morisita))
comm.dist.morisita$tree1 <- rownames(comm.dist.morisita)

#--Make dataframe
mor.pin <- gather(comm.dist.morisita,'tree2', 'dissimilarity', -tree1)
#--Add fire history data
for(i in mor.pin$tree1){
  for(t in mor.pin$tree2){
    mor.pin[mor.pin$tree1 == i, 'fire.history.1'] <- morisita.matrix[morisita.matrix$Tree == i, 'Burn_status']
    mor.pin[mor.pin$tree2 == t, 'fire.history.2'] <- morisita.matrix[morisita.matrix$Tree == t, 'Burn_status']
  }
}
#--assess if within same or different fire history
for(i in 1:nrow(mor.pin)){
  if(mor.pin[i, 'fire.history.1'] == mor.pin[i, 'fire.history.2']){
    mor.pin[i, 'comp'] <- 'Same'
  } else{mor.pin[i, 'comp'] <- 'Different'}
}

#--Remove self-comparisons (dissimilarity = 0) and outliers
mor.pin <- mor.pin[mor.pin$dissimilarity > 0, ]
mor.pin.out <- mor.pin[!mor.pin$dissimilarity < 0.85,]
#jac.pin.out <- jac.pin[jac.pin$comp == 'Different' & jac.pin$dissimilarity < 0.96,]

#--wilcox test
wilcox.test(dissimilarity ~ comp, data = mor.pin.out)

ggplot(mor.pin.out, mapping = aes(x = comp,
                                  y = dissimilarity)) +
  geom_boxplot() +
  xlab('Fire history') +
  ylab('Morisita-Horn dissimilarity') +
  theme_classic()

#========================================================================================#
# Jaccard based dissimilarity index: Santa Catalina----
#========================================================================================#
#--Remove outlier
scm.matrix <- scm.matrix[!scm.matrix$Tree == 'LB056',]
#--isolate otu data
comm.matrix <- scm.matrix[12:length(scm.matrix)]

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
jpeg(filename = 'figures/NMDS_SC_Jaccard.jpeg',
     width = 700, height = 600,
     quality = 100)
par(mfrow = c(1,1), "mar"=c(6, 5, 5, 3))

#--Plot NMDS of EM community based on Jaccard index and OTU abundance
plot(jaccard.otu, display = "sites", type = "n", cex.lab = 2,
     cex.axis = 1.5, xlab = 'Axis 1', ylab = 'Axis 2')

# colors for points
color.vec <- data.frame(color = rep(NA, nrow(jaccard.matrix)),
                        p.group = jaccard.matrix$Burn_status)
color.vec[color.vec$p.group == 'burned', 'color'] <- 'black'
color.vec[color.vec$p.group == 'unburned', 'color'] <- 'black'

color.vec[color.vec$p.group == 'burned', 'shape'] <- 16
color.vec[color.vec$p.group == 'unburned', 'shape'] <- 1

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
legend("bottomleft", legend = c('Burned','Unburned'), bty = "n",
       col = c('black','black'),
       pch = c(16,1),
       pt.bg =  c('black','black'), cex = 2)

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

#<< Pairwise >>-----------------------------------------------------------------------------
rownames(scm.matrix) <- scm.matrix$Tree

#--isolate otu data
comm.matrix <- scm.matrix[12:length(scm.matrix)]

#--comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 4]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 1, ] # remove rows with sums of 0
jaccard.matrix <- scm.matrix[row.names(comm.matrix),]

#--distance matrix using jaccard index
comm.dist.jaccard <- vegdist(comm.matrix, method = "jaccard", binary = F, upper = T)
comm.dist.jaccard <- as.data.frame(as.matrix(comm.dist.jaccard))
comm.dist.jaccard$tree1 <- rownames(comm.dist.jaccard)

#--Make dataframe
jac.scm <- gather(comm.dist.jaccard,'tree2', 'dissimilarity', -tree1)
#--Add fire history data
for(i in jac.scm$tree1){
  for(t in jac.scm$tree2){
    jac.scm[jac.scm$tree1 == i, 'fire.history.1'] <- jaccard.matrix[jaccard.matrix$Tree == i, 'Burn_status']
    jac.scm[jac.scm$tree2 == t, 'fire.history.2'] <- jaccard.matrix[jaccard.matrix$Tree == t, 'Burn_status']
  }
}
#--assess if within same or different fire history
for(i in 1:nrow(jac.scm)){
  if(jac.scm[i, 'fire.history.1'] == jac.scm[i, 'fire.history.2']){
    jac.scm[i, 'comp'] <- 'Same'
  } else{jac.scm[i, 'comp'] <- 'Different'}
}

#--Remove self-comparisons (dissimilarity = 0) and outliers
jac.scm <- jac.scm[jac.scm$dissimilarity > 0, ]
jac.scm.out <- jac.scm[!jac.scm$dissimilarity < 0.83,]
#jac.pin.out <- jac.pin[jac.pin$comp == 'Different' & jac.pin$dissimilarity < 0.96,]

#--wilcox test
wilcox.test(dissimilarity ~ comp, data = jac.scm.out)

ggplot(jac.scm.out, mapping = aes(x = comp,
                                  y = dissimilarity)) +
  geom_boxplot() +
  xlab('Fire history') +
  ylab('Jaccard dissimilarity') +
  theme_classic()

#========================================================================================#
# Morisita based dissimilarity index: Santa Catalina----
#========================================================================================#

#--isolate site X species matrix only without metadata
comm.matrix <- scm.matrix[12:length(scm.matrix)]

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
jpeg(filename = 'figures/NMDS_SC_MorisitaHorn.jpeg',
     width = 700, height = 600,
     quality = 100)
par(mfrow = c(1,1), "mar"=c(6, 5, 5, 3))

#--Plot NMDS of EM community based on Jaccard index and OTU abundance
plot(morisita.otu, display = "sites", type = "n", cex.lab = 2,
     cex.axis = 1.5, xlab = 'Axis 1', ylab = 'Axis 2')
# colors for points
color.vec <- data.frame (color = rep (NA, nrow(morisita.matrix)),
                         p.group = morisita.matrix$Burn_status)
color.vec[color.vec$p.group == 'burned', 'color'] <- 'black'
color.vec[color.vec$p.group == 'unburned', 'color'] <- 'black'
color.vec[color.vec$p.group == 'burned', 'shape'] <- 16
color.vec[color.vec$p.group == 'unburned', 'shape'] <- 1

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
       col = c('black','black'),
       pch = c(16,1),
       pt.bg =  c('black','black'), cex = 2)

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

#<< Pairwise >>-----------------------------------------------------------------------------
rownames(scm.matrix) <- scm.matrix$Tree

#--isolate otu data
comm.matrix <- scm.matrix[12:length(scm.matrix)]

#--comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 4]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 1, ] # remove rows with sums of 0
morisita.matrix <- scm.matrix[row.names(comm.matrix),]

#--distance matrix using jaccard index
comm.dist.morisita <- vegdist(comm.matrix, method = "horn", binary = F, upper = T)
comm.dist.morisita <- as.data.frame(as.matrix(comm.dist.morisita))
comm.dist.morisita$tree1 <- rownames(comm.dist.morisita)

#--Make dataframe
mor.scm <- gather(comm.dist.morisita,'tree2', 'dissimilarity', -tree1)
#--Add fire history data
for(i in mor.scm$tree1){
  for(t in mor.scm$tree2){
    mor.scm[mor.scm$tree1 == i, 'fire.history.1'] <- morisita.matrix[jaccard.matrix$Tree == i, 'Burn_status']
    mor.scm[mor.scm$tree2 == t, 'fire.history.2'] <- morisita.matrix[jaccard.matrix$Tree == t, 'Burn_status']
  }
}
#--assess if within same or different fire history
for(i in 1:nrow(mor.scm)){
  if(mor.scm[i, 'fire.history.1'] == mor.scm[i, 'fire.history.2']){
    mor.scm[i, 'comp'] <- 'Same'
  } else{mor.scm[i, 'comp'] <- 'Different'}
}

#--Remove self-comparisons (dissimilarity = 0) and outliers
mor.scm <- mor.scm[mor.scm$dissimilarity > 0, ]
same.mor.scm <- mor.scm[mor.scm$comp == 'Same',] 
same.mor.scm <- same.mor.scm[same.mor.scm$dissimilarity > 0.4,]
diff.mor.scm <- mor.scm[mor.scm$comp == 'Different',]
diff.mor.scm <- diff.mor.scm[diff.mor.scm$dissimilarity > 0.57,]
mor.scm.out <- bind_rows(same.mor.scm, diff.mor.scm)

#--wilcox test
wilcox.test(dissimilarity ~ comp, data = mor.scm.out)

ggplot(mor.scm.out, mapping = aes(x = comp,
                                  y = dissimilarity)) +
  geom_boxplot() +
  xlab('Fire history') +
  ylab('Morisita-Horn dissimilarity') +
  theme_classic()
