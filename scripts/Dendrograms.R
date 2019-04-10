## Script created by Liz Bowman June 9, 2017
## for analyzing taxonomic differences using a heat map; as a complement to
## 20170603_Taxonomic_analyses.R scipt output

#========================================================================================#
# Load data and libraries----
#========================================================================================#

#----------------------------------------------------------------------------------------#
# Load libraries----
#----------------------------------------------------------------------------------------#
library(ggplot2)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(vegan)
#install.packages('dendextend')
library(dendextend)

#----------------------------------------------------------------------------------------#
# set up paths to directories----
#----------------------------------------------------------------------------------------#
#--path to directory where the climate data downloaded from BIOCLIM site is stored
dat.dir <- "~/Documents/PhD/2_EM_Fire_effect/data/"
fig.dir <- '~/Documents/PhD/2_EM_Fire_effect/figures_output/'
res.dir <- "~/Documents/PhD/2_EM_Fire_effect/results/"

#----------------------------------------------------------------------------------------#
# Load data: shortened----
#----------------------------------------------------------------------------------------#

#<< Taxonomic/OTU data >>----------------------------------------------

#--taxonomic data
tax.data <- read.csv(paste0(dat.dir,'20170806_OTU_data.csv'), as.is = T)

#--filter out doug fir hosts
tax.data <- tax.data[tax.data$Host == 'Ponderosa',]

#--Load data
all.data <- read.csv(paste0(dat.dir,'20170806_OTU_data.csv'), as.is = T)

#--Remove sequences not assigned to an OTU
otu.data <- all.data[!is.na(all.data$otu.97), ]

#--Remove sequences not assigned to Ponderosa host
otu.data <- otu.data[otu.data$Host == 'Ponderosa',]

#--Add an OTU count column (1 sequence = 1 frequency count)
otu.data['otu.count'] <- 1

#--Creates matrix, grouping OTUs by burn status and range
otu.data %>%
  select (Sample_name,Burn_status,Range,Site,Tree,otu.97,otu.count) %>%
  spread (otu.97,otu.count) %>%
  group_by(Site,Range, Burn_status) %>%
  select (matches ("otu.*")) %>%
  summarize_each (funs (sum (., na.rm = TRUE))) %>%
  as.data.frame () -> site.matrix

#--add concatenated range_burn column
for(i in 1:nrow(site.matrix)) {
  site.matrix[i,'range.burn'] <-
    paste(site.matrix[i,'Range'], site.matrix[i, 'Burn_status'])
}
site.matrix <- site.matrix[c(1:3,120, 4:119)]

#--Creates matrix, grouping OTUs by site, range, and burn status
otu.data %>%
  select (Sample_name,Burn_status,Range,Site,Tree,otu.97,otu.count) %>%
  spread (otu.97,otu.count) %>%
  group_by(Site,Range,Burn_status) %>%
  select (matches ("otu.*")) %>%
  summarize_each (funs (sum (., na.rm = TRUE))) %>%
  as.data.frame () -> rangeburn.matrix

#<< Mean climate characteristics per site >>----------------------------------------------

#--climate data
clim.data <- read.csv(paste0(dat.dir, 'climate_data.csv'))
#--add range data
scm <- c('sc1','sc2','sc3','sc4')
clim.data[clim.data$site %in% scm, 'range'] <- 'santa.catalina'
clim.data[!clim.data$site %in% scm, 'range'] <- 'pinaleno'

#--create separate data frame for means
clim.mean <- data.frame(site = unique(clim.data$site))
burned <- c('sc3','sc4','p3','p4')
clim.mean[clim.mean$site %in% scm, 'range'] <- 'santa.catalina'
clim.mean[!clim.mean$site %in% scm, 'range'] <- 'pinaleno'
clim.mean[clim.mean$site %in% burned, 'burn_status'] <- 'burned'
clim.mean[!clim.mean$site %in% burned, 'burn_status'] <- 'unburned'

for(s in unique(clim.mean$site)) {
  clim.mean[clim.mean$site == s, 'prec'] <-
    mean(clim.data[clim.data$site == s, 'prec'])
  clim.mean[clim.mean$site == s, 'Tmax'] <-
    mean(clim.data[clim.data$site == s, 'Tmax'])
}


#<< Mean soil characteristics per site >>-------------------------------------------------
soil.data <- read.csv(paste0(dat.dir, 'soil_data.csv'), as.is = T, header = T)
#--create separate data frame for means
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

#-----------------------------------------------------------------------------------------
# Load data: per tree----
#-----------------------------------------------------------------------------------------

#<< Taxonomic/OTU data >>----------------------------------------------

#--taxonomic data
tax.data <- read.csv(paste0(dat.dir,'20170806_OTU_data.csv'), as.is = T)

#--filter out doug fir hosts
tax.data <- tax.data[tax.data$Host == 'Ponderosa',]

#--Load data
all.data <- read.csv(paste0(dat.dir,'20170806_OTU_data.csv'), as.is = T)

#--Remove sequences not assigned to an OTU
otu.data <- all.data[!is.na(all.data$otu.97), ]

#--Remove sequences not assigned to Ponderosa host
otu.data <- otu.data[otu.data$Host == 'Ponderosa',]

#--Add an OTU count column (1 sequence = 1 frequency count)
otu.data['otu.count'] <- 1

#--Creates matrix, grouping OTUs by burn status and range
otu.data %>%
  select (Sample_name,Burn_status,Range,Site,Tree,otu.97,otu.count) %>%
  spread (otu.97,otu.count) %>%
  group_by(Tree, Site, Range, Burn_status) %>%
  select (matches ("otu.*")) %>%
  summarize_each (funs (sum (., na.rm = TRUE))) %>%
  as.data.frame () -> stsp.matrix

#--add concatenated range_burn column
for(i in 1:nrow(stsp.matrix)) {
  stsp.matrix[i,'range.burn'] <-
    paste(stsp.matrix[i,'Range'], stsp.matrix[i, 'Burn_status'])
}

#--rearrange stsp.matrix
stsp.matrix <- stsp.matrix[c(1:4,121,5:120)]

#<< Mean climate characteristics per site >>----------------------------------------------

#--climate data
clim.data <- read.csv(paste0(dat.dir, 'climate_data.csv'))
#--add range data
scm <- c('sc1','sc2','sc3','sc4')
clim.data[clim.data$site %in% scm, 'range'] <- 'santa.catalina'
clim.data[!clim.data$site %in% scm, 'range'] <- 'pinaleno'

#<< Mean soil characteristics per site >>-------------------------------------------------
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

# #--rearrange data frame
soil.data <- soil.data[c(1:3,21:22,4:20)]

#========================================================================================#
# Hierarchical cluster analysis----
#========================================================================================#
#----------------------------------------------------------------------------------------#
# OTUs grouped by site----
#----------------------------------------------------------------------------------------#
morisita.dist <- vegdist(rangeburn.matrix[4:length(rangeburn.matrix)], method = 'horn',
                         binary = F, na.rm = T)
jaccard.dist <- vegdist(rangeburn.matrix[4:length(rangeburn.matrix)], method = 'jaccard')

#========================================================================================#
# Cluster analyses with site descriptions----
#========================================================================================#
rownames(rangeburn.matrix) <- rangeburn.matrix$Site
rangeburn.matrix[4:length(rangeburn.matrix)] %>%
  vegdist(method = 'horn') %>%
  hclust(method = 'average') %>%
  as.dendrogram -> ab.dend
plot(ab.dend)
labels(ab.dend)
ab.dend %>%
  set("labels", c('SCM unburned','SCM burned','SCM burned','PM burned',
                  'PM burned','SCM unburned','PM unburned','PM unburned')) %>% 
  plot


#----------------------------------------------------------------------------------------#
# By site: Jaccard-----
#----------------------------------------------------------------------------------------#

jpeg(file = paste0(fig.dir,'Cluster_abundance_site_Jaccard.jpeg'),
     width = 800, height = 1200,
     quality = 100)

#divide your plotting device in 3 plotting regions with different dimensions
layout(matrix(c(1,2,3), 1, 3, byrow = TRUE), widths = c(4,2,2)) 
#--abundance
matrix.abu <- as.matrix(site.matrix[,5:119])
rownames(matrix.abu) <- site.matrix[,4]

#--presence/absence
matrix.pa <- as.matrix(site.matrix[,5:119])
rownames(matrix.pa) <- site.matrix[,4]

#Calculate distances among communities using the 
dist.abu <- vegdist(matrix.abu, method="horn") #to calculate the distances 
dist.pres <- vegdist(matrix.pa, method="jaccard")

#Hierarchical cluster analysis
ca.abu <- hclust(dist.abu,"average") 
ca.pres <- hclust(dist.pres, "average")

#Abundances
plot(as.dendrogram(ca.pres), horiz=T, hang = -1, 
     leaflab = 'none', cex = 3, lwd = 2)
title(xlab="Dissimilarity (Jaccard)", cex.lab=2)
mtext("Jaccard", cex = 3)

# plot(as.dendrogram(ca.pres), horiz=T)
# title(xlab="Dissimilarity (Jaccard)", cex.lab=2)
# mtext("Presence/Absence", cex = 1)

#--data frame for metadata
meta.matrix <- data.frame(range.burn = rownames(matrix.pa),
           range = rangeburn.matrix$Range,
           burn_status = rangeburn.matrix$Burn_status)

# #--Morisita-horn labels
# hm <- data.frame(rcol = c(1,1,1,1,2,2,2,2),
#            bcol = c(3,4,3,4,4,4,3,3))
# col.range <- brewer.pal(2,'Greens')
# col.burn <- c("#FDAE6B" ,"#E6550D")
# #<< Draw the colors: Morisita-horn >> ----------------------
# par(mai=c(.68,.1,.48,.1))
# image(z=t(hm), yaxt="n", xaxt="n",
#       col = c(col.range, col.burn)) 

#--Jaccard labels
hj <- data.frame(rcol = c(2,2,1,2,2,1,1,1),
                 bcol = c(3,3,3,4,4,4,3,4))

labels(as.dendrogram(ca.pres))
col.range <- brewer.pal(2,'Greens')
col.burn <- c("#FDAE6B" ,"#E6550D")
#<< Draw the colors: Jaccard >> ----------------------
par(mai=c(.95,.1,.95,.1))
image(z=t(hj), yaxt="n", xaxt="n",
      col = c(col.range, col.burn)) 

#<< draw the color code >> ------------------
par(mai=c(.7,.1,.7,.1))
plot(c(0,6), c(-5,24), type = "n", xlab = "", ylab = "", bty="n", yaxt="n",xaxt="n")

#--Range -------------------------
rect(xleft=c(1,1,1,1), 
     xright=c(1.6,1.6), 
     ybottom=c(15,16), 
     ytop=c(16,17), 
     col=rev(col.range), border="black")
text(x=1.3,y=18+0.5, labels="Range", font=2, cex=2)
text(1.7,16+0.3, labels="Pinaleno Mts.", font=3, cex=2, pos=4)
text(1.7,15+0.3, labels="Santa Catalina Mts.", font=3, cex=2, pos=4)

#--Burn_status
rect(xleft=c(1,1), 
     xright=c(1.6,1.6),
     ybottom=c(8,9),
     ytop=c(9,10), 
     col=rev(col.burn), border="black")
text(x=1.3,y=11+0.5, labels="Burn status", font=2, cex=2)
text(1.7,8+0.3, labels="Burned", font=1, cex=2, pos=4)
text(1.7,9+0.3, labels="Unburned", font=1, cex=2, pos=4)

dev.off()


#----------------------------------------------------------------------------------------#
# By site: Morisita-Horn-----
#----------------------------------------------------------------------------------------#

jpeg(file = paste0(fig.dir,'Cluster_abundance_site_MorisitaHorn.jpeg'),
     width = 800, height = 1200,
     quality = 100)

#divide your plotting device in 3 plotting regions with different dimensions
layout(matrix(c(1,2,3), 1, 3, byrow = TRUE), widths = c(4,2,2)) 
#--abundance
matrix.abu <- as.matrix(site.matrix[,5:119])
rownames(matrix.abu) <- site.matrix[,4]

#--presence/absence
matrix.pa <- as.matrix(site.matrix[,5:119])
rownames(matrix.pa) <- site.matrix[,4]

#Calculate distances among communities using the 
dist.abu <- vegdist(matrix.abu, method="horn") #to calculate the distances 
dist.pres <- vegdist(matrix.pa, method="jaccard")

#Hierarchical cluster analysis
ca.abu <- hclust(dist.abu,"average") 
ca.pres <- hclust(dist.pres, "average")

#Abundances
plot(as.dendrogram(ca.abu), horiz=T, hang = -1, 
     leaflab = 'none', cex = 5, lwd = 3)
title(xlab="Dissimilarity (Morisita-Horn)", cex.lab=2)
mtext("Morisita-Horn", cex = 3)

# plot(as.dendrogram(ca.pres), horiz=T)
# title(xlab="Dissimilarity (Jaccard)", cex.lab=2)
# mtext("Presence/Absence", cex = 1)

#--data frame for metadata
meta.matrix <- data.frame(range.burn = rownames(matrix.abu),
                          range = rangeburn.matrix$Range,
                          burn_status = rangeburn.matrix$Burn_status)

# #--Morisita-horn labels
# hm <- data.frame(rcol = c(1,1,1,1,2,2,2,2),
#            bcol = c(3,4,3,4,4,4,3,3))
# col.range <- brewer.pal(2,'Greens')
# col.burn <- c("#FDAE6B" ,"#E6550D")
# #<< Draw the colors: Morisita-horn >> ----------------------
# par(mai=c(.68,.1,.48,.1))
# image(z=t(hm), yaxt="n", xaxt="n",
#       col = c(col.range, col.burn)) 

#--Jaccard labels
hj <- data.frame(rcol = c(2,2,1,2,2,1,1,1),
                 bcol = c(3,3,3,4,4,4,3,4))

labels(as.dendrogram(ca.abu))
col.range <- brewer.pal(2,'Greens')
col.burn <- c("#FDAE6B" ,"#E6550D")
#<< Draw the colors: Jaccard >> ----------------------
par(mai=c(.95,.1,.95,.1))
image(z=t(hj), yaxt="n", xaxt="n",
      col = c(col.range, col.burn)) 

#<< draw the color code >> ------------------
par(mai=c(.7,.1,.7,.1))
plot(c(0,6), c(-5,24), type = "n", xlab = "", ylab = "", bty="n", yaxt="n",xaxt="n")

#--Range -------------------------
rect(xleft=c(1,1,1,1), 
     xright=c(1.6,1.6), 
     ybottom=c(15,16), 
     ytop=c(16,17), 
     col=rev(col.range), border="black")
text(x=1.3,y=18+0.5, labels="Range", font=2, cex=2)
text(1.7,16+0.3, labels="Pinaleno Mts.", font=3, cex=2, pos=4)
text(1.7,15+0.3, labels="Santa Catalina Mts.", font=3, cex=2, pos=4)

#--Burn_status
rect(xleft=c(1,1), 
     xright=c(1.6,1.6),
     ybottom=c(8,9),
     ytop=c(9,10), 
     col=rev(col.burn), border="black")
text(x=1.3,y=11+0.5, labels="Burn status", font=2, cex=2)
text(1.7,8+0.3, labels="Burned", font=1, cex=2, pos=4)
text(1.7,9+0.3, labels="Unburned", font=1, cex=2, pos=4)

dev.off()

#-----------------------------------------------------------------------------------------
# by tree-----
#-----------------------------------------------------------------------------------------

pdf(file=paste0(fig.dir,'Cluster_abundance_tree.pdf'),
    paper = "special", width = 6, height = 8, family="Times")

#divide your plotting device in 3 plotting regions with different dimensions
layout(matrix(c(1,2,3), 1, 3, byrow = TRUE), widths = c(4,2,2)) 

#--abundance
#--comment to include singletons
abu.meta <- stsp.matrix
matrix.abu <- abu.meta[6:length(abu.meta)]
matrix.abu <- matrix.abu[colSums(matrix.abu) >= 2]
matrix.abu <- matrix.abu[rowSums(matrix.abu) > 0, ] 
abu.meta <- abu.meta[row.names(matrix.abu),]
rownames(matrix.abu) <- abu.meta$Tree
matrix.abu <- as.matrix(matrix.abu)

#--presence/absence
pa.meta <- stsp.matrix
matrix.pa <- pa.meta[6:length(pa.meta)]
matrix.pa <- matrix.pa[colSums(matrix.pa) >= 2]
matrix.pa <- matrix.pa[rowSums(matrix.pa) > 0, ] 
pa.meta <- pa.meta[row.names(matrix.pa),]

#Calculate distances among communities using the 
dist.abu <- vegdist(matrix.abu, method="horn") #to calculate the distances 
dist.pres <- vegdist(matrix.pa, method="jaccard")

#Hierarchical cluster analysis
ca.abu <- hclust(dist.abu,"average") 
ca.pres <- hclust(dist.pres, "average")

#Abundances
plot(as.dendrogram(ca.abu), horiz=T, leaflab = 'none')
title(xlab="Dissimilarity (Morisita-Horn)", cex.lab=1)
mtext("Abundance", cex = 1)

# plot(as.dendrogram(ca.pres), horiz=T)
# title(xlab="Dissimilarity (Jaccard)", cex.lab=2)
# mtext("Presence/Absence", cex = 1)

#--data frame for metadata
labels <- c('NF20','NF17','NF18','NF12','NF11','LB060','F5','F17','F14','LB059','F8',
            'F2','F1','LB021','LB019','F11','LB055','LB020','F7','F10','F6','NF15',
            'F9','F4','LB023','F3','NF13','LB057','F16','F12','LB058','F19','F13',
            'NF14','F18','F15','LB056','LB024','F20','NF16','NF19')
clust.labels <- NA
for(i in labels){
  lab.i <- abu.meta[abu.meta$Tree == i, 'range.burn']
  clust.labels <- c(clust.labels, lab.i)
}
clust.labels <- clust.labels[2:length(clust.labels)]

rcol <- clust.labels
rcol <- replace(rcol,rcol == 'santa.catalina burned' |
                  rcol == 'santa.catalina unburned', 1)
rcol <- replace(rcol,rcol == 'pinaleno burned' | rcol == 'pinaleno unburned', 2)
rcol <- as.numeric(rcol)

bcol <- clust.labels
bcol <- replace(bcol,bcol == 'pinaleno burned' | bcol == 'santa.catalina burned', 4)
bcol <- replace(bcol,bcol == 'pinaleno unburned' | bcol == 'santa.catalina unburned', 3)
bcol <- as.numeric(bcol)
  
hm <- data.frame(rcol = rcol,
                 bcol = bcol)

col.range <- brewer.pal(2,'Greens')
col.burn <- c("#FDAE6B" ,"#E6550D")

#<< Draw the colors >> ----------------------
par(mai=c(0.85,0.1,0.85,0.1))
image(z=t(hm), yaxt="n", xaxt="n",
      col = c(col.range, col.burn)) 

#<< draw the color code >> ------------------
par(mai=c(.7,.1,.7,.1))
plot(c(0,6), c(-5,24), type = "n", xlab = "", ylab = "", bty="n", yaxt="n",xaxt="n")

#--Range
rect(xleft=c(1,1,1,1), 
     xright=c(1.6,1.6), 
     ybottom=c(15,16), 
     ytop=c(16,17), 
     col=rev(col.range), border="black")
text(x=1.3,y=18+0.5, labels="Range", font=2, cex=.8)
text(1.7,16+0.3, labels="Pinaleno Mts.", font=1, cex=.7, pos=4)
text(1.7,15+0.3, labels="Santa Catalina Mts.", font=1, cex=.7, pos=4)

#--Burn_status
rect(xleft=c(1,1), 
     xright=c(1.6,1.6),
     ybottom=c(8,9),
     ytop=c(9,10), 
     col=rev(col.burn), border="black")
text(x=1.3,y=11+0.5, labels="Burn status", font=2, cex=.8)
text(1.7,8+0.3, labels="Burned", font=1, cex=.7, pos=4)
text(1.7,9+0.3, labels="Unburned", font=1, cex=.7, pos=4)

dev.off()

#<< Pairwise comparison plot >> ----------------------
source('~/Documents/PhD/3_Sky_islands-FE/scripts/functions.R')
#--Range information
range.anosim <- anosim(stsp.matrix[6:length(stsp.matrix)], 
                       grouping = stsp.matrix$Range,
                       permutations = 999,
                       distance = 'horn')
plot.anosim(range.anosim)

#--Burn information
burn.anosim <- anosim(stsp.matrix[6:length(stsp.matrix)], 
                      grouping = stsp.matrix$Burn_status,
                      permutations = 999,
                      distance = 'horn')
plot.anosim(burn.anosim)

#--Range and burn information
rangeburn.anosim <- anosim(stsp.matrix[6:length(stsp.matrix)], 
                      grouping = stsp.matrix$range.burn,
                      permutations = 999,
                      distance = 'horn')
plot.anosim(rangeburn.anosim)

#-----------------------------------------------------------------------------------------
# Soil-----
#-----------------------------------------------------------------------------------------

pdf(file=paste0(fig.dir,'Cluster_soil.pdf'),
    paper = "special", width = 6, height = 8, family="Times")

#divide your plotting device in 3 plotting regions with different dimensions
layout(matrix(c(1,2,3), 1, 3, byrow = TRUE), widths = c(4,2,2)) 
#--soil data isolated
matrix.soil <- as.matrix(soil.data[,6:length(soil.data)])
rownames(matrix.soil) <- soil.data[,6]

#Calculate distances among communities using the 
dist.soil <- vegdist(matrix.soil, method="euclidean")

#Hierarchical cluster analysis
ca.soil <- hclust(dist.soil,"average") 

#Abundances
plot(as.dendrogram(ca.soil), horiz=T, leaflab = 'none', cex = 2,
     hang = -1)
title(xlab="Dissimilarity (Euclidean)", cex.lab=2)
mtext("Soil", cex = 2)

# plot(as.dendrogram(ca.pres), horiz=T)
# title(xlab="Dissimilarity (Jaccard)", cex.lab=2)
# mtext("Presence/Absence", cex = 1)

#--data frame for metadata
meta.matrix <- data.frame(range.burn = rownames(matrix.soil),
                          range = soil.data$range,
                          burn_status = soil.data$burn_status)
clust.labels <- c('pinaleno burned','pinaleno burned','pinaleno unburned',
                  'pinaleno burned', 'pinaleno unburned', 'pinaleno burned',
                  'pinaleno unburned','santa.catalina burned','santa.catalina unburned',
                  'santa.catalina burned','santa.catalina unburned',
                  'santa.catalina unburned','santa.catalina unburned',
                  'santa.catalina burned','santa.catalina burned',
                  'santa.catalina unburned','pinaleno unburned','santa.catalina burned',
                  'pinaleno unburned', 'pinaleno burned','pinaleno burned',
                  'santa.catalina burned','santa.catalina unburned','pinaleno unburned')
rcol <- rev(clust.labels)
rcol <- replace(rcol,rcol == 'santa.catalina burned' |
                  rcol == 'santa.catalina unburned', 1)
rcol <- replace(rcol,rcol == 'pinaleno burned' | rcol == 'pinaleno unburned', 2)
rcol <- as.numeric(rcol)

bcol <- rev(clust.labels)
bcol <- replace(bcol,bcol == 'pinaleno burned' | bcol == 'santa.catalina burned', 4)
bcol <- replace(bcol,bcol == 'pinaleno unburned' | bcol == 'santa.catalina unburned', 3)
bcol <- as.numeric(bcol)

hm <- data.frame(rcol = rcol,
                 bcol = bcol)

col.range <- brewer.pal(2,'Greens')
col.burn <- c("#FDAE6B" ,"#E6550D")

#<< Draw the colors >> ----------------------
par(mai=c(.86,.1,.7,.1))
image(z=t(hm), yaxt="n", xaxt="n",
      col = c(col.range, col.burn)) 

#<< draw the color code >> ------------------
par(mai=c(.7,.1,.7,.1))
plot(c(0,6), c(-5,24), type = "n", xlab = "", ylab = "", bty="n", yaxt="n",xaxt="n")

#--Range
rect(xleft=c(1,1,1,1), 
     xright=c(1.6,1.6), 
     ybottom=c(15,16), 
     ytop=c(16,17), 
     col=rev(col.range), border="black")
text(x=1.3,y=18+0.5, labels="Range", cex=2)
text(1.7,16+0.3, labels="Pinaleno Mts.", cex=1, pos=4)
text(1.7,15+0.3, labels="Santa Catalina Mts.", cex=1, pos=4)

#--Burn_status
rect(xleft=c(1,1), 
     xright=c(1.6,1.6),
     ybottom=c(8,9),
     ytop=c(9,10), 
     col=rev(col.burn), border="black")
text(x=1.3,y=11+0.5, labels="Burn \n status", cex=2)
text(1.7,8+0.3, labels="Burned", cex=1, pos=4)
text(1.7,9+0.3, labels="Unburned", cex=1, pos=4)

dev.off()
