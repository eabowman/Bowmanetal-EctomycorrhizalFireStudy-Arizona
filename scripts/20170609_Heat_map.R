## Script created by Liz Bowman June 9, 2017
## for analyzing taxonomic differences using a heat map; as a complement to
## 20170603_Taxonomic_analyses.R scipt output

#=========================================================================================
# Load data and libraries----
#=========================================================================================

#-----------------------------------------------------------------------------------------
# Load libraries----
#-----------------------------------------------------------------------------------------
library(ggplot2)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(vegan)

#-----------------------------------------------------------------------------------------
# set up paths to directories----
#-----------------------------------------------------------------------------------------
#--path to directory where the climate data downloaded from BIOCLIM site is stored
dat.dir <- "~/Documents/PhD/3_EM_Fire_effect/data/"
fig.dir <- '~/Documents/PhD/3_EM_Fire_effect/figures/'
res.dir <- "~/Documents/PhD/3_EM_Fire_effect/results/"

#-----------------------------------------------------------------------------------------
# Load data: shortened----
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
  clim.mean[clim.mean$site == s, 'temp.warmest.quarter'] <-
    mean(clim.data[clim.data$site == s, 'BIO10_red'])
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

#=========================================================================================
# Heat Map Basic----
#=========================================================================================
#-----------------------------------------------------------------------------------------
# Create data table with x, y, and z axes defined-----
#-----------------------------------------------------------------------------------------

#--check the structure of the data file
str(tax.data)

#--remove unclassified ones
tax.data <- tax.data[!tax.data$Taxonomy_class == '',]

#--table 
tax.table <- as.data.frame(table(tax.data$Site, tax.data$Taxonomy_class))
colnames(tax.table) <- c('site','taxonomy','count')
for (i in unique(tax.table$site)) {
  tax.table[tax.table$site == i, 'range'] <- unique(tax.data[tax.data$Site == i, 'Range'])
  tax.table[tax.table$site == i, 'burn_status'] <-
    unique(tax.data[tax.data$Site == i, 'Burn_status'])
}

tax.table

#--set data up for heat map with a column for the x, y, and z axes
#tax.table %>%
#  select(burn_status,range,site,taxonomy,count) %>%
#  spread(key = taxonomy, value = count) -> tax.long
#--samples names, class, and abundance are axes x, y, and z, respectively
#tax.long

#-----------------------------------------------------------------------------------------
# Plot-----
#-----------------------------------------------------------------------------------------

#--Remove rare species (Saccharomycetes, Eurotiomycetes)
tax.table <- tax.table[!tax.table$taxonomy %in% 
                         c('Saccharomycetes', 'Eurotiomycetes','Archaeorhizomycetes'),]

#<< Basic >>------------------------------------------------------------------------------
heatmap <- ggplot(data=tax.table,
                       mapping = aes(x = site,
                                     y = taxonomy,
                                     fill = count)) + 
  #fill will add color to heat map
  geom_tile() + # geom_tile creates the heat map using the fill info above
  xlab(label = 'Sample') +
  facet_grid(~ range, switch = 'x', scales = 'free_x', space = 'free_x')
# This will facet on column Depth
# facet_grid creates side by side plots for each of the categories within Depth
# switch = 'x' displays the top labels on the bottom; depth categories moved to bottom
# scales = 'free_x' and space = 'free_x' removes excess space
heatmap

#<< Fancy >>------------------------------------------------------------------------------

heatmap <- ggplot(data=tax.table,
                       mapping = aes(x = site,
                                     y = taxonomy,
                                     fill = count)) + 
  #fill will add color to heat map
  geom_tile() + # geom_tile creates the heat map using the fill info above
  xlab(label = '') +
  facet_grid(~ burn_status, switch = 'x', scales = 'free_x', space = 'free_x') +
  # This will facet on column Depth
  # facet_grid creates side by side plots for each of the categories within Depth
  # switch = 'x' displays the top labels on the bottom; depth categories moved to bottom
  # scales = 'free_x' and space = 'free_x' removes excess space
  scale_fill_gradient(name = 'Frequency',
                      low = '#FFFFFF',
                      high = '#012345') +
  # scale_fill_gradient changes the gradient and color of the fill and the legend title
  theme_bw() + # adds styling to plot; there are a variety of them
  theme(strip.placement = 'outside',
        # moves the depth boxes below the sample names
        plot.title = element_text(hjust = '0.5'),
        plot.subtitle = element_text(hjust = '0.5'),
        axis.title.y = element_blank(),
        # centers title and subtitle and removes y-axis label
        strip.background = element_rect(fill = '#CCCCCC',
                                        color = '#000000')) +
  # changes color of depth category boxes
  #ggtitle(label = 'Microbe Class Diversity', subtitle = 'R workshop') +
  scale_y_discrete(limits = rev(levels(as.factor(c('Agaricomycetes','Dothideomycetes',
                                                   'Leotiomycetes','Pezizomycetes')))))

heatmap

#=========================================================================================
# Hierarchical cluster analysis----
#=========================================================================================
#-----------------------------------------------------------------------------------------
# OTUs grouped by site----
#-----------------------------------------------------------------------------------------
morisita.dist <- vegdist(rangeburn.matrix[4:length(rangeburn.matrix)], method = 'horn',
                         binary = F, na.rm = T)
jaccard.dist <- vegdist(rangeburn.matrix[4:length(rangeburn.matrix)], method = 'jaccard')

morisita.clust <- hclust(morisita.dist)
jaccard.clust <- hclust(jaccard.dist)

plot(morisita.clust,labels = rangeburn.matrix$Site)
plot(jaccard.clust,labels = paste(rangeburn.matrix$Range, rangeburn.matrix$Burn_status))

#-----------------------------------------------------------------------------------------
# Add dendrogram to heat map: shortened-----
#-----------------------------------------------------------------------------------------

#divide your plotting device in 3 plotting regions with different dimensions
layout(matrix(c(1,2), 2, 1, byrow = TRUE), widths = c(4,4), heights = c(8,4)) 

morisita.clust <- hclust(morisita.dist, 'average')

plot(as.dendrogram(morisita.clust), horiz = F, leaflab = rangeburn.matrix$Site)
title(xlab="Dissimilarity (Morisita-Horn)", cex.lab=1)


#=========================================================================================
# Cluster analyses with site descriptions----
#=========================================================================================

#-----------------------------------------------------------------------------------------
# By site-----
#-----------------------------------------------------------------------------------------

pdf(file=paste0(fig.dir,'Cluster_abundance_site.pdf'),
    paper = "special", width = 6, height = 8, family="Times")

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
plot(as.dendrogram(ca.pres), horiz=T, leaflab = 'none')
title(xlab="Dissimilarity (Morisita-Horn)", cex.lab=1)
mtext("Abundance", cex = 1)

# plot(as.dendrogram(ca.pres), horiz=T)
# title(xlab="Dissimilarity (Jaccard)", cex.lab=2)
# mtext("Presence/Absence", cex = 1)

#--data frame for metadata
meta.matrix <- data.frame(range.burn = rownames(matrix.abu),
           range = stsp.matrix$Range,
           burn_status = stsp.matrix$Burn_status)

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
hj <- data.frame(rcol = c(2,2,2,2,1,1,1,1),
                 bcol = c(3,4,4,3,3,4,3,4))
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
text(x=1.3,y=18+0.5, labels="Range", font=2, cex=.8)
text(1.7,16+0.3, labels="Pinaleno Mts.", font=3, cex=.7, pos=4)
text(1.7,15+0.3, labels="Santa Catalina Mts.", font=3, cex=.7, pos=4)

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

#-----------------------------------------------------------------------------------------
# by tree-----
#-----------------------------------------------------------------------------------------

pdf(file=paste0(fig.dir,'Cluster_abundance_tree.pdf'),
    paper = "special", width = 6, height = 8, family="Times")

#divide your plotting device in 3 plotting regions with different dimensions
layout(matrix(c(1,2,3), 1, 3, byrow = TRUE), widths = c(4,2,2)) 
#--abundance
matrix.abu <- as.matrix(stsp.matrix[,6:119])
rownames(matrix.abu) <- stsp.matrix[,5]

#--presence/absence
matrix.pa <- as.matrix(rangeburn.matrix[,4:118])
rownames(matrix.pa) <- rangeburn.matrix[,5]

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
meta.matrix <- data.frame(range.burn = rownames(matrix.abu),
                          range = stsp.matrix$Range,
                          burn_status = stsp.matrix$Burn_status)
clust.labels <- c('pinaleno unburned','pinaleno unburned','santa.catalina unburned',
                  'pinaleno burned', 'santa.catalina unburned', 'pinaleno burned',
                  'pinaleno burned','pinaleno burned','pinaleno unburned','pinaleno burned',
                  'pinaleno unburned','santa.catalina burned','santa.catalina unburned',
                  'santa.catalina burned','santa.catalina burned','pinaleno unburned',
                  'pinaleno burned','santa.catalina burned','santa.catalina burned',
                  'santa.catalina burned','santa.catalina unburned','santa.catalina unburned',
                  'santa.catalina unburned','santa.catalina unburned','santa.catalina burned',
                  'santa.catalina burned','santa.catalina burned','pinaleno unburned',
                  'pinaleno burned','pinaleno burned','santa.catalina burned',
                  'santa.catalina unburned','pinaleno unburned','pinaleno unburned',
                  'pinaleno unburned','pinaleno unburned','pinaleno unburned',
                  'santa.catalina unburned','santa.catalina unburned','pinaleno burned',
                  'pinaleno burned')
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
plot(as.dendrogram(ca.soil), horiz=T, leaflab = 'none')
title(xlab="Dissimilarity (Euclidean)", cex.lab=1)
mtext("Soil", cex = 1)

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
