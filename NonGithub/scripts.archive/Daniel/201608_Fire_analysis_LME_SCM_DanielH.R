# Script created by Liz Bowman on March 15, 2016
# for analysis of effects of fire on EM communities in the
# Santa Catalina Mts.
# eabowman@email.arizona.edu

#=========================================================================================
# Load libraries and data
#=========================================================================================
# << Libraries >> ------------------------------------------------------------------------
#install.packages('tidyr')
#install.packages('vegan')
#install.packages('dplyr')
#install.packages('ggplot2')
library (tidyr)
library (vegan)
library (dplyr)
library (ggplot2)

# << data paths >> -----------------------------------------------------------------------
dat.dir <- '~/Documents/PhD/3_EM_Fire_effect/data/'
fig.dir <- '~/Documents/PhD/3_EM_Fire_effect/figures/'
res.dir <- '~/Documents/PhD/3_EM_Fire_effect/results/'

# << Load data >> ------------------------------------------------------------------------
#--Burned site data
burned.seq <- read.csv (paste0 (dat.dir, 'Burned_data.csv'), as.is = T)
burned.sites <- read.csv (paste0 (dat.dir, 'Burned_sites.csv'), as.is = T)
burned.morpho <- read.csv (paste0 (dat.dir, 'Burn_sites_morphotypes.csv'), as.is = T)
burned.stats <- read.csv (paste0 (dat.dir, 'Burned_stats_core.csv'), as.is = T)
#--Unburned site data
unburned.seq <- read.csv (paste0 (dat.dir, 'Unburned_data.csv'), as.is = T)
unburned.stats <- read.csv (paste0 (dat.dir, 'Unburned_stats.csv'), as.is = T)
#--Soil
soil <- read.csv (paste0 (dat.dir, 'Soil.csv'), as.is = T)
#--OTU_clustering
otu.all <- read.csv(paste0 (dat.dir, 'All_OTU.csv'), as.is = T)

#=========================================================================================
# Clean-up data
#=========================================================================================

#--changed core to microsite and c to uphill, b to parallel, a to downhill
colnames (burned.seq) [7] <- 'microsite'
burned.seq[burned.seq$microsite == 'A', 'microsite'] <- 'Downhill'
burned.seq[burned.seq$microsite == 'B', 'microsite'] <- 'Parallel'
burned.seq[burned.seq$microsite == 'C', 'microsite'] <- 'Uphill'

#--Add unburned tree number to OTU.all
seq <- unburned.seq$sequencing_number
for (i in seq) {
  otu.all[otu.all$Sampling_number == i, 'tree_number'] <-
    unburned.seq[unburned.seq$sequencing_number == i, 'tree_number']
  otu.all[otu.all$Sampling_number == i, 'elevation'] <-
    unburned.seq[unburned.seq$sequencing_number == i, 'elevation']
  otu.all[otu.all$Sampling_number == i, 'sequence'] <-
    unburned.seq[unburned.seq$sequencing_number == i, 'sequence']
  otu.all[otu.all$Sampling_number == i, 'microsite'] <-
    unburned.seq[unburned.seq$sequencing_number == i, 'microsite']
}
#--Add burned tree number to OTU.all
seq <- burned.seq$sample
for (i in seq) {
  otu.all [otu.all$Sampling_number == i, 'tree_number'] <-
    burned.seq [burned.seq$sample == i, 'tree_number']
  otu.all [otu.all$Sampling_number == i, 'elevation'] <-
    burned.seq [burned.seq$sample == i, 'elevation']
  otu.all [otu.all$Sampling_number == i, 'sequence'] <-
    burned.seq [burned.seq$sample == i, 'sequence']
}

#write.csv(otu.all, paste0 (dat.dir,'otu_all_info.csv')

#--add burned and unburned data to otu.all
burn <- c(2183.6, 2181.8)
unburn <- c(2170, 2201)
otu.all [otu.all$elevation %in% burn, 'group'] <- 'Burned'
otu.all [otu.all$elevation %in% unburn, 'group'] <- 'Unburned'

#--add column with OTU count
otu.all$otu.count <- 1
otu.all %>%
  select (OTU.95,Sampling_number,otu.count,microsite,tree_number,elevation,group) %>%
  spread (OTU.95, otu.count) %>%
  group_by (tree_number, elevation, group) %>%
  select (matches ("OTU.*")) %>%
  summarize_each (funs (sum (., na.rm = TRUE))) %>%
  as.data.frame () -> otu.data

#=========================================================================================
# Preliminary
#=========================================================================================

# << Compare burned abundance based on dry root weight and total tip count >> ------------
tree <- burned.stats$root.core [c(1, 3:6, 8:length(burned.stats$root.core))]
for (i in tree) {
  burned.stats [burned.stats$root.core == i, 'ab.dry'] <-
    burned.stats [burned.stats$root.core == i, 'colonized.root.tips'] /
    burned.stats [burned.stats$root.core == i, 'dry.root.weight']
  burned.stats [burned.stats$root.core == i, 'ab.total.rt.tip'] <-
    burned.stats [burned.stats$root.core == i, 'colonized.root.tips'] /
    burned.stats [burned.stats$root.core == i, 'total.root.tips']
}

#--check for correlation between two different methods to calculate abundance
#x <- lm (burned.stats$ab.total.rt.tip ~ burned.stats$ab.dry)
#summary (x)
#plot (burned.stats$ab.dry, burned.stats$ab.total.rt.tip)
#abline(x)

#=========================================================================================
# Soil
#=========================================================================================

#--add Burned Unburned to environmental data
soil$site [1:6] <- 'Burned'
soil$site [7:12] <- 'Unburned'

#--replace <1.0 with 0
soil [soil$no3.n.ppm == '<1.0', 'no3.n.ppm'] <- 0

#--make all numeric data numeric
soil$no3.n.ppm <- as.numeric(soil$no3.n.ppm)
soil$ca.ppm <- as.numeric(soil$ca.ppm)
soil$mg.ppm <- as.numeric(soil$mg.ppm)
soil$k.ppm <- as.numeric(soil$k.ppm)
soil$fe.ppm <- as.numeric(soil$fe.ppm)
soil$mn.ppm <- as.numeric(soil$mn.ppm)
soil$site <- as.factor (soil$site)

#--make results data frame
soil.results <- data.frame (soil_trait = colnames(soil[3:19]), stat = NA, df = NA, p = NA)

#--pH vs burn status; not sig
boxplot(soil$pH.su ~ soil$site)
pH.soil <- kruskal.test(soil$pH.su ~ soil$site)
soil.results[soil.results$soil_trait == 'pH.su', 'stat'] <- pH.soil$statistic
soil.results[soil.results$soil_trait == 'pH.su', 'df'] <- pH.soil$parameter
soil.results[soil.results$soil_trait == 'pH.su', 'p'] <- pH.soil$p.value

#--electrical conductivity vs burn status; not sig
boxplot(soil$EC.ds.m ~ soil$site)
EC.soil <- kruskal.test(soil$EC.ds.m ~ soil$site)
soil.results[soil.results$soil_trait == 'EC.ds.m', 'stat'] <- EC.soil$statistic
soil.results[soil.results$soil_trait == 'EC.ds.m', 'df'] <- EC.soil$parameter
soil.results[soil.results$soil_trait == 'EC.ds.m', 'p'] <- EC.soil$p.value

#--Ca vs burn status; not sig
boxplot(soil$ca.ppm ~ soil$site)
ca.soil <- kruskal.test(soil$ca.ppm ~ soil$site)
soil.results[soil.results$soil_trait == 'ca.ppm', 'stat'] <- ca.soil$statistic
soil.results[soil.results$soil_trait == 'ca.ppm', 'df'] <- ca.soil$parameter
soil.results[soil.results$soil_trait == 'ca.ppm', 'p'] <- ca.soil$p.value

#--Mg vs burn status; not sig
boxplot(soil$mg.ppm ~ soil$site)
mg.soil <- kruskal.test(soil$mg.ppm ~ soil$site)
soil.results[soil.results$soil_trait == 'mg.ppm', 'stat'] <- mg.soil$statistic
soil.results[soil.results$soil_trait == 'mg.ppm', 'df'] <- mg.soil$parameter
soil.results[soil.results$soil_trait == 'mg.ppm', 'p'] <- mg.soil$p.value

#--Na vs burn status; not sig
boxplot(soil$na.ppm ~ soil$site)
na.soil <- kruskal.test(soil$na.ppm ~ soil$site)
soil.results[soil.results$soil_trait == 'na.ppm', 'stat'] <- na.soil$statistic
soil.results[soil.results$soil_trait == 'na.ppm', 'df'] <- na.soil$parameter
soil.results[soil.results$soil_trait == 'na.ppm', 'p'] <- na.soil$p.value

#--K vs burn status; not sig
boxplot(soil$k.ppm ~ soil$site)
k.soil <- kruskal.test(soil$k.ppm ~ soil$site)
soil.results[soil.results$soil_trait == 'k.ppm', 'stat'] <- k.soil$statistic
soil.results[soil.results$soil_trait == 'k.ppm', 'df'] <- k.soil$parameter
soil.results[soil.results$soil_trait == 'k.ppm', 'p'] <- k.soil$p.value

#--Zn vs burn status; not sig
boxplot(soil$zn.ppm ~ soil$site)
zn.soil <- kruskal.test(soil$zn.ppm ~ soil$site)
soil.results[soil.results$soil_trait == 'zn.ppm', 'stat'] <- zn.soil$statistic
soil.results[soil.results$soil_trait == 'zn.ppm', 'df'] <- zn.soil$parameter
soil.results[soil.results$soil_trait == 'zn.ppm', 'p'] <- zn.soil$p.value

#--Fe vs burn status; not sig
boxplot(soil$fe.ppm ~ soil$site)
fe.soil <- kruskal.test(soil$fe.ppm ~ soil$site)
soil.results[soil.results$soil_trait == 'fe.ppm', 'stat'] <- fe.soil$statistic
soil.results[soil.results$soil_trait == 'fe.ppm', 'df'] <- fe.soil$parameter
soil.results[soil.results$soil_trait == 'fe.ppm', 'p'] <- fe.soil$p.value

#--Mn vs burn status; not sig
boxplot(soil$mn.ppm ~ soil$site)
mn.soil <- kruskal.test(soil$mn.ppm ~ soil$site)
soil.results[soil.results$soil_trait == 'mn.ppm', 'stat'] <- mn.soil$statistic
soil.results[soil.results$soil_trait == 'mn.ppm', 'df'] <- mn.soil$parameter
soil.results[soil.results$soil_trait == 'mn.ppm', 'p'] <- mn.soil$p.value

#--Cu vs burn status; not sig
boxplot(soil$cu.ppm ~ soil$site)
cu.soil <- kruskal.test(soil$cu.ppm ~ soil$site)
soil.results[soil.results$soil_trait == 'cu.ppm', 'stat'] <- cu.soil$statistic
soil.results[soil.results$soil_trait == 'cu.ppm', 'df'] <- cu.soil$parameter
soil.results[soil.results$soil_trait == 'cu.ppm', 'p'] <- cu.soil$p.value

#--Ni vs burn status; not sig
boxplot(soil$ni.ppm ~ soil$site)
ni.soil <- kruskal.test(soil$ni.ppm ~ soil$site)
soil.results[soil.results$soil_trait == 'ni.ppm', 'stat'] <- ni.soil$statistic
soil.results[soil.results$soil_trait == 'ni.ppm', 'df'] <- ni.soil$parameter
soil.results[soil.results$soil_trait == 'ni.ppm', 'p'] <- ni.soil$p.value

#--remove LB020 because it is an outlier (see below)
soil.no3 <- soil [ !soil$no3.n.ppm == 3, ]
#--NO3.N vs burn status *
boxplot(soil.no3$no3.n.ppm ~ soil.no3$site, ylab = 'PPM')
no3.soil <- kruskal.test(soil.no3$no3.n.ppm ~ soil.no3$site)
#--check outliers
#--unburn outlier above is 3 (outside 2sd); tree LB020 NEEDS TO BE REMOVED
#unburn.no3.mean <- mean (soil [soil$site == 'Unburned', 'no3.n.ppm'])
#unburn.no3.sd <- sd (soil [soil$site == 'Unburned', 'no3.n.ppm'])
#unburn.no3.mean + 2*unburn.no3.sd
#--add stats to results table
soil.results[soil.results$soil_trait == 'no3.n.ppm', 'stat'] <- no3.soil$statistic
soil.results[soil.results$soil_trait == 'no3.n.ppm', 'df'] <- no3.soil$parameter
soil.results[soil.results$soil_trait == 'no3.n.ppm', 'p'] <- no3.soil$p.value

#--PO4.p vs burn status *
boxplot(soil$po4.p.ppm ~ soil$site, ylab = 'PPM')
po4.soil <- kruskal.test(soil$po4.p.ppm ~ soil$site)
#--check outliers
#--burn outlier is 13 (not outside 2 sd)
#burn.po4.mean <- mean (soil [soil$site == 'Burned', 'po4.p.ppm'])
#burn.po4.sd <- sd (soil [soil$site == 'Burned', 'po4.p.ppm'])
#burn.po4.mean + 2*burn.po4.sd
#--unburn outlier above is 14 (not outside 2sd); below is 4.0 (not outside 2sd)
#unburn.po4.mean <- mean (soil [soil$site == 'Unburned', 'po4.p.ppm'])
#unburn.po4.sd <- sd (soil [soil$site == 'Unburned', 'po4.p.ppm'])
#unburn.po4.mean + 2*unburn.po4.sd
#unburn.po4.mean - 2*unburn.po4.sd
#--add stats to results table
soil.results[soil.results$soil_trait == 'po4.p.ppm', 'stat'] <- po4.soil$statistic
soil.results[soil.results$soil_trait == 'po4.p.ppm', 'df'] <- po4.soil$parameter
soil.results[soil.results$soil_trait == 'po4.p.ppm', 'p'] <- po4.soil$p.value

#--remove outliers (LB020 and F9)
soil.so4 <- soil [!soil$tree_number %in% c('F9','LB020'), ]
#--SO4.S vs burn status *
boxplot(soil.so4$so4.s.ppm ~ soil.so4$site, ylab = 'PPM')
so4.soil <- kruskal.test(soil.so4$so4.s.ppm ~ soil.so4$site)
#--check outliers
#--burn outlier is 7.2 which should be remove; remove tree F9
#burn.so4.mean <- mean (soil [soil$site == 'Burned', 'so4.s.ppm'])
#burn.so4.sd <- sd (soil [soil$site == 'Burned', 'so4.s.ppm'])
#burn.so4.mean + 2*burn.so4.sd
#--unburn outlier above is 21 (not outside 2sd); below is 4.8, remove tree LB020
#unburn.so4.mean <- mean (soil [soil$site == 'Unburned', 'so4.s.ppm'])
#unburn.so4.sd <- sd (soil [soil$site == 'Unburned', 'so4.s.ppm'])
#unburn.so4.mean + 2*unburn.so4.sd
#unburn.so4.mean - 2*unburn.so4.sd
#--add stats to results table
soil.results[soil.results$soil_trait == 'so4.s.ppm', 'stat'] <- so4.soil$statistic
soil.results[soil.results$soil_trait == 'so4.s.ppm', 'df'] <- so4.soil$parameter
soil.results[soil.results$soil_trait == 'so4.s.ppm', 'p'] <- so4.soil$p.value

#--B vs burn status *
boxplot(soil$b.ppm ~ soil$site, ylab = 'PPM')
b.soil <- kruskal.test(soil$b.ppm ~ soil$site)
#--add stats to results table
soil.results[soil.results$soil_trait == 'b.ppm', 'stat'] <- b.soil$statistic
soil.results[soil.results$soil_trait == 'b.ppm', 'df'] <- b.soil$parameter
soil.results[soil.results$soil_trait == 'b.ppm', 'p'] <- b.soil$p.value

#--ESP vs burn status; not sig
boxplot(soil$esp ~ soil$site)
esp.soil <- kruskal.test(soil$esp ~ soil$site)
#--add stats to results table
soil.results[soil.results$soil_trait == 'esp', 'stat'] <- esp.soil$statistic
soil.results[soil.results$soil_trait == 'esp', 'df'] <- esp.soil$parameter
soil.results[soil.results$soil_trait == 'esp', 'p'] <- esp.soil$p.value

#--CEC vs burn status; not sig
boxplot(soil$cec.meq.100g ~ soil$site)
cec.soil <- kruskal.test(soil$cec.meq.100g ~ soil$site)
#--add stats to results table
soil.results[soil.results$soil_trait == 'cec.meq.100g', 'stat'] <- cec.soil$statistic
soil.results[soil.results$soil_trait == 'cec.meq.100g', 'df'] <- cec.soil$parameter
soil.results[soil.results$soil_trait == 'cec.meq.100g', 'p'] <- cec.soil$p.value

#--write out results table
write.csv (soil.results, paste0(res.dir, 'Soil_results.csv'))

# << check for correlation between soil traits (B, PO4, SO4, NO3) >> ---------------------
#--SO4 and NO3 are correlated
plot (soil$so4.s.ppm, soil$no3.n.ppm)
z <- lm (no3.n.ppm ~ so4.s.ppm, data = soil)
abline (z)
anova (z)

#--SO4 and PO4 are correlated
plot (soil$so4.s.ppm, soil$po4.p.ppm)
z <- lm (po4.p.ppm ~ so4.s.ppm, data = soil)
abline (z)
anova (z)

#--SO4 and B are correlated
plot (soil$so4.s.ppm, soil$b.ppm)
z <- lm (b.ppm ~ so4.s.ppm, data = soil)
abline (z)
anova (z)

#--NO3 and PO4 are NOT correlated
plot (soil$no3.n.ppm, soil$po4.p.ppm)
z <- lm (po4.p.ppm ~no3.n.ppm, data = soil)
abline (z)
anova (z)

#--NO3 and B are NOT correlated
plot (soil$no3.n.ppm, soil$b.ppm)
z <- lm (b.ppm ~ no3.n.ppm, data = soil)
abline (z)
anova (z)

#--B and PO4 are NOT correlated
plot(soil$po4.p.ppm, soil$b.ppm)
z <- lm (b.ppm ~ po4.p.ppm, data = soil)
abline (z)
anova (z)

#--Therefore use PO4, B, and NO3 for PERMANOVA analyses

#=========================================================================================
# Abundance
#=========================================================================================
#--Be sure to run preliminary section before this section

# << Check outliers >> -------------------------------------------------------------------
#--data before outliers were removed
ab.out <- data.frame (site = c(rep('Burned', 28), rep ('Unburned', 11)), abundance = NA)
ab.out[1:10, 'abundance'] <- burned.stats [c(1:10), 'ab.dry']
ab.out[11:21, 'abundance'] <- unburned.stats [c(1:11), 'abund.em']
#--burn sites
burn.mean <- mean (ab.out [ab.out$site == 'Burned', 'abundance'])
burn.sd <- sd (ab.out [ab.out$site == 'Burned', 'abundance'])
#--outlier for burned sites is 65
burn.mean + 2*burn.sd #--outlier for burned can be removed
#--unburned sites
unburn.mean <- mean (ab.out [ab.out$site == 'Unburned', 'abundance'])
unburn.sd <- sd (ab.out [ab.out$site == 'Unburned', 'abundance'])
#--outliers for unburned sites is 81 and 99
unburn.mean + 2*unburn.sd #--99 can be removed

#--check for normaily
shapiro.test(ab [ab$site == 'Burned', 'abundance'])
qqnorm(ab [ab$site == 'Burned', 'abundance'])
shapiro.test (ab [ab$site == 'Unburned', 'abundance'])
qqnorm(ab [ab$site == 'Unburned', 'abundance'])

# << Analyses w/o outliers >> ------------------------------------------------------------

#--make single file with burned and unburned abundance data
#--outliers were removed, see above
ab <- data.frame (site = c(rep('Burned', 9), rep ('Unburned', 10)), abundance = NA)
ab[1:9, 'abundance'] <- burned.stats [c(1:6, 8:10), 'ab.dry']
ab[10:19, 'abundance'] <- unburned.stats [c(1:9, 11), 'abund.em']

#--plot of burned versus unburned sites
#--check distrubtion; outliers are they 2 sd away from the mean
plot (ab$site, ab$abundance, ylab = 'Abundance')
#--data is normally distributed, see below
ab.ttest <- t.test (abundance ~ site, data = ab)

#--make results table
abund.res <- data.frame(t.stat = ab.ttest$statistic,
                        df = ab.ttest$parameter,
                        p = ab.ttest$p.value)

#--write out results table
write.csv (abund.res, paste0 (res.dir, 'Abundance_ttest_results.csv'))

#=========================================================================================
# Diversity
#=========================================================================================
#--isolate otu data
comm.matrix <- otu.data[4:59]
#--create data frame for fisher's alpha data
fish <- data.frame (tree_number = otu.data$tree_number, elevation = otu.data$elevation,
                    group = otu.data$group)
#--calculate fisher's alpha
fish$f.a. <- fisher.alpha(comm.matrix)

#--check for normality
shapiro.test(fish$f.a [fish$group == 'Burned'])
qqnorm(fish$f.a. [fish$group == 'Burned']) #--remove F9
shapiro.test(fish$f.a [fish$group == 'Unburned'])
qqnorm(fish$f.a. [fish$group == 'Unburned']) #--remove LB021

#--remove outliers (see above)
outliers <- c('F9', 'LB021', 'F4')
fish <- fish [!fish$tree_number %in% outliers, ]

#--plot of burned versus unburned sites
#--check distrubtion; outliers are they 2 sd away from the mean
plot (fish$group, fish$f.a., ylab = 'Fishers alpha')
#--data is normally distributed, see below
div.ttest <- t.test (f.a. ~ group, data = fish)

#--make results table
div.res <- data.frame(t.stat = div.ttest$statistic,
                      df = div.ttest$parameter,
                      p = div.ttest$p.value)

#--write out results table
write.csv (div.res, paste0 (res.dir, 'Diversity_ttest_results.csv'))

#=========================================================================================
# Community Structure: NMDS and ANOSIM*******Permanova/remove outlier ****Betadisper
#=========================================================================================
#--Table of Anosim results
anosim.res <- data.frame(test = c('morisita', 'jaccard'), p = NA, r = NA)
#--isolate otu data
comm.matrix <- otu.data[4:60]
#--Remove otu with singletons
comm.matrix <- comm.matrix[which (colSums(comm.matrix) >= 2)]

#--remove outlier LB056
#comm.matrix <- otu.data [-17,4:59]
#otu.dat.out <- otu.data[-17,]

#<< MORISTA >> ---------------------------------------------------------------------------
#--check for empty columns and rows
colSums(comm.matrix)
rowSums(comm.matrix)
#--distance matrix using MOrisita index
comm.dist.horn <- vegdist(comm.matrix, method = "morisita", binary = F)

#--Check variance of distance matrix; p < 0.05 means that data is not normally distributed
#--can't use this data for anosim; must use non=parametric analyses (e.g. permanova)
anova(betadisper(comm.dist.horn, otu.data$group))

#--NMDS analysis
horn.abund <- metaMDS(comm.dist.horn, dist = "bray", permutations = 999)

#--create output file
jpeg(filename = paste0 (fig.dir,"NMDS_mor_jac.jpeg"), width = 1100, height = 600)
#--Combine into one figure
par(mfrow = c(1,2), "mar"=c(6, 5, 5, 3))

#--Plot NMDS of EM community based on Jaccard index and OTU abundance
plot(horn.abund, display = "sites", type = "n", cex.lab = 1.5,
     cex.axis = 1.5, yaxt = "n")
axis (2, at = seq (-0.4, 0.4, by = 0.2), cex.axis = 1.5, las = 2)
# colors for points
color.vec <- data.frame (color = rep(NA, nrow(comm.matrix)),
                         p.group = otu.data$group)
color.vec <- sapply (color.vec$p.group, function (x) 
  if (x == 'Burned') {color.vec = 'red'} else {color.vec = 'black'})
#ordipointlabel(horn.abund, display = "sites")
points(horn.abund, display = "sites", cex = 2, pch = 20,
       col = color.vec,
       bg = color.vec)
#legend("topright", legend = levels(otu.env.data$group), bty = "n",
#       col = color.vec, pch = 21, pt.bg = color.vec, cex = 2)
ordihull(horn.abund, groups = otu.data[-53,'group'])

#--ANOSIM: Morisita
horn.anosim <- anosim (comm.matrix, grouping =  otu.data$group, distance = "morisita")
#--Add results to data frame
anosim.res [ which (anosim.res$test =="morisita"), "r"] <- horn.anosim$statistic
anosim.res [ which (anosim.res$test =="morisita"), "p"] <- horn.anosim$signif

#<< Jaccard >> ---------------------------------------------------------------------------
#--isolate otu data
comm.matrix <- otu.data[5:61]
#--Remove otu with singletons
comm.matrix <- comm.matrix[which (colSums(comm.matrix) >= 2)]
#--LB059 (row 53) is empty, remove
comm.matrix <- comm.matrix[-53,]

#--distance matrix using jaccard index
comm.dist.jaccard <- vegdist(comm.matrix, method = "jaccard", binary = TRUE)

#--Check variance of distance matrix; p < 0.05 means that data is not normally distributed
#--can't use this data for anosim; must use non=parametric analyses (e.g. permanova)
anova(betadisper(comm.dist.jaccard,  otu.data[-53,'group']))

#--NMDS analysis
jaccard.abund <- metaMDS(comm.dist.jaccard, dist = "bray", permutations = 999)

#--Plot NMDS of EM community based on Jaccard index and OTU abundance
plot(jaccard.abund, display = "sites", type = "n", cex.lab = 1.5,
     cex.axis = 1.5, yaxt = "n")
axis (2, at = seq (-0.4, 0.4, by = 0.2), cex.axis = 1.5, las = 2)
# colors for points
color.vec <- data.frame (color = rep (NA, nrow(comm.matrix)),
                         p.group =  otu.data[-53,'group'])
color.vec <- sapply (color.vec$p.group, function (x) 
  if (x == 'Burned') {color.vec = 'black'} else {color.vec = 'green'})
#ordipointlabel(morisita.abund, display = "sites")
points(jaccard.abund, display = "sites", cex = 2, pch = 20,
       col = color.vec,
       bg = color.vec)
#legend("topright", legend = levels(otu.env.data$group), bty = "n",
#       col = color.vec, pch = 21, pt.bg = color.vec, cex = 2)
ordihull(jaccard.abund, groups =  otu.data[-53,'group'])

#--close figure path
dev.off()

#--ANOSIM: Jaccard
jaccard.anosim <- anosim (comm.matrix, grouping = otu.data$group, distance = "jaccard")
#--Add results to data frame
anosim.res [ which (anosim.res$test =="jaccard"), "r"] <- jaccard.anosim$statistic
anosim.res [ which (anosim.res$test =="jaccard"), "p"] <- jaccard.anosim$signif

#--write out results
write.csv (anosim.res, paste0 (res.dir, 'Anosim_results.csv'))

#=========================================================================================
# Community Structure: PERMANOVA
#=========================================================================================
#--Add soil data to otu.data file
elevation <- unique (soil$elevation)
for ( i in elevation) {
  mean.po4 <- soil [soil$elevation == i, 'po4.p.ppm'] 
}

#--distance matrix using MOrisita index
comm.dist.morisita <- vegdist(comm.matrix, method = "morisita", binary = F)

#--adonis
morisita.perma <- adonis(comm.dist.morisita ~ group, data = otu.data)

#--make results table
permanova <- data.frame (df = morisita.perma$aov.tab$Df,
                         f.model = morisita.perma$aov.tab$F.Model,
                         r = morisita.perma$aov.tab$R2, 
                         p = morisita.perma$aov.tab$`Pr(>F)`)

#--write otu results
write.csv (permanova, paste0 (res.dir, 'permanova_results_morisita.csv'))

#=========================================================================================
# Taxonomy
#=========================================================================================
#--Make results table
chi <- data.frame (chi = NA, df = NA, p = NA)

#--Make count table of em classes
tax.em <- table(otu.all$group,
                otu.all$taxonomy_class)
#--Remove taxa with less than 1 occurence
#rare.em <- "Fungi incertae sedis"
#elvgr.em <- elvgr.em [ , which (!colnames(elvgr.em) %in% rare.em)]

#--Chi square test
em.chi <- chisq.test(tax.em)

#--Add results to results table; chi-square, degrees of freedom, and the p-value
chi$chi <- em.chi$statistic[[1]]
chi$df <- em.chi$parameter[[1]]
chi$p <- em.chi$p.value[[1]]
#--writ out results table
write.csv (chi, paste0(res.dir, 'taxonomy_results.csv'))

#--Remove taxa with less than 1 occurence
#rare.em <- c("Fungi incertae sedis", "Eurotiomycetes", "Sordariomycetes")
#em.chi <- em.data [which (!em.data$taxonomy_class %in% rare.em), ]
#--Bar graph
ggplot(data = otu.all, 
       aes(x = group,
           fill = taxonomy_class)) + 
  geom_bar(position = "fill") + 
  ylab("Proportion of sequences per class") +
  xlab("Fire status") +
  #ggtitle("Proportion of Classes by Topography") +
  scale_fill_brewer(palette = "BrBG") +
  guides (fill=guide_legend(title=NULL)) +
  theme_classic(base_size = 12)

#--export graph
ggsave ('taxonomy.jpeg', plot = last_plot(), width = 5, height = 5, device = "jpeg",
        path = fig.dir, dpi = 1000)
