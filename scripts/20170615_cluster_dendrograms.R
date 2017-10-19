## Script created by Liz Bowman June 3, 2017
## for analyzing community differences between ranges and burn status of sites

#=========================================================================================
# Load data and libraries----
#=========================================================================================

#-----------------------------------------------------------------------------------------
# Load libraries----
#-----------------------------------------------------------------------------------------
library(ggplot2);library(tidyr);library(vegan);library(dplyr)

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

#--Creates matrix, grouping OTUs by tree number
otu.data %>%
  select (Sample_name,Burn_status,Range,Site,Tree,otu.97,otu.count) %>%
  spread (otu.97,otu.count) %>%
  group_by(Tree,Site,Range,Burn_status) %>%
  select (matches ("otu.*")) %>%
  summarize_each (funs (sum (., na.rm = TRUE))) %>%
  as.data.frame () -> stsp.matrix

#--climate data: add to stsp.matrix
clim.data <- read.csv(paste0(dat.dir, 'climate_data.csv'), as.is = T)
for(i in unique(stsp.matrix$Site)) {
  stsp.matrix[stsp.matrix$Site == i, 'prec'] <- clim.data[clim.data$site == i, 'prec']
  stsp.matrix[stsp.matrix$Site == i, 'temp.warmest.quarter'] <-
    clim.data[clim.data$site == i, 'BIO11_red']
}
#--reorder matrix
stsp.matrix <- stsp.matrix[c(1:4,121,122,5:120)]

#--Creates matrix, grouping OTUs by site, range, and burn status
otu.data %>%
  select (Sample_name,Burn_status,Range,Site,Tree,otu.97,otu.count) %>%
  spread (otu.97,otu.count) %>%
  group_by(Site,Range,Burn_status) %>%
  select (matches ("otu.*")) %>%
  summarize_each (funs (sum (., na.rm = TRUE))) %>%
  as.data.frame () -> rangeburn.matrix

#--Soil data
soil.data <- read.csv(paste0(dat.dir, 'soil_data.csv'),as.is = T)

#--change < 1.0 in the nitrate column to 0. 
soil.data[soil.data$no3.n.ppm == '<1.0', 'no3.n.ppm'] <- 0
soil.data$no3.n.ppm <- as.numeric(soil.data$no3.n.ppm)

#--add range and burn information to soil data data frame
scm <- c('sc1','sc2','sc3','sc4')
burned <- c('sc3','sc4','p3','p4')
soil.data[soil.data$site %in% scm, 'range'] <- 'santa.catalina'
soil.data[!soil.data$site %in% scm, 'range'] <- 'pinaleno'
soil.data[soil.data$site %in% burned, 'burn_status'] <- 'burned'
soil.data[!soil.data$site %in% burned, 'burn_status'] <- 'unburned'

#--rearrange data frame
soil.data <- soil.data[c(1:3,21:22,4:20)]

#=========================================================================================
# cluster matrix: fungal community----
#=========================================================================================

#-----------------------------------------------------------------------------------------
# OTUs grouped by site----
#-----------------------------------------------------------------------------------------
morisita.dist <- vegdist(rangeburn.matrix[4:length(rangeburn.matrix)], method = 'horn',
                         binary = F, na.rm = T)
jaccard.dist <- vegdist(rangeburn.matrix[4:length(rangeburn.matrix)], method = 'jaccard')

morisita.clust <- hclust(morisita.dist)
jaccard.clust <- hclust(jaccard.dist)

plot(morisita.clust,labels = paste(rangeburn.matrix$Range, rangeburn.matrix$Burn_status))
plot(jaccard.clust,labels = paste(rangeburn.matrix$Range, rangeburn.matrix$Burn_status))

#-----------------------------------------------------------------------------------------
# OTUs grouped by tree----
#-----------------------------------------------------------------------------------------

# morisita.dist <- vegdist(stsp.matrix[7:length(stsp.matrix)], method = 'horn',
#                          binary = F, na.rm = T)
# jaccard.dist <- vegdist(stsp.matrix[7:length(stsp.matrix)], method = 'jaccard')
# 
# morisita.clust <- hclust(morisita.dist)
# jaccard.clust <- hclust(jaccard.dist)
# 
# plot(morisita.clust,labels = paste(stsp.matrix$Range, stsp.matrix$Burn_status))
# plot(jaccard.clust,labels = paste(stsp.matrix$Range, stsp.matrix$Burn_status))

#=========================================================================================
# Cluster matrix: Soil----
#=========================================================================================

#<< Mean soil characteristics per site >>-------------------------------------------------
soil.mean <- data.frame(site = unique(soil.data$site))
soil.mean[soil.mean$site %in% scm, 'range'] <- 'santa.catalina'
soil.mean[!soil.mean$site %in% scm, 'range'] <- 'pinaleno'
soil.mean[soil.mean$site %in% burned, 'burn_status'] <- 'burned'
soil.mean[!soil.mean$site %in% burned, 'burn_status'] <- 'unburned'

for(s in unique(soil.data$site)) {
  for(c in colnames(soil.data[6:length(soil.data)])) {
    soil.mean[soil.mean$site == s, c] <- mean(soil.data[soil.data$site == s, c])
  }
}

soil.mean.dist <- vegdist(soil.mean[4:length(soil.mean)], method = 'euclidean')

soil.mean.clust <- hclust(soil.mean.dist)

plot(soil.mean.clust,
     labels = soil.mean$site)

#<< Soil characteristics >>---------------------------------------------------------------
soil.dist <- vegdist(soil.data[4:length(soil.data)], method = 'euclidean')

soil.clust <- hclust(soil.dist)

plot(soil.clust,labels = soil.data$site)

#=========================================================================================
# Cluster matrix: Weather----
#=========================================================================================

clim.dist <- vegdist(stsp.matrix[5:6], method = 'euclidean')

clim.clust <- hclust(clim.dist)

plot(clim.clust, labels = paste(stsp.matrix$Range, stsp.matrix$Burn_status))

