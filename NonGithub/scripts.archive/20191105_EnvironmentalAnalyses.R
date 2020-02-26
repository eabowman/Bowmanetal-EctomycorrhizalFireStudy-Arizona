## Script created by Liz Bowman July 19, 2017
## for analyzing soil differences between sites, ranges, and burn

#========================================================================================#
# Load data and libraries----
#========================================================================================#

#----------------------------------------------------------------------------------------#
# Load libraries----
#----------------------------------------------------------------------------------------#
# install.packages('ggplot2'); install.packages('tidyr'); install.packages('vegan')
# install.packages('dplyr'); install.packages('gridExtra')
library(ggplot2);library (tidyr);library (vegan);library (dplyr); library(gridExtra)
# install.packages('ggfortify'); installed.packages('plyr'); install.packages('nlme')
library(ggfortify); library(plyr); library(nlme)

#----------------------------------------------------------------------------------------#
# set up paths to directories----
#----------------------------------------------------------------------------------------#
#--path to directory 
dat.dir <- "~/Documents/PhD/2_EM_Fire_effect/data/"
fig.dir <- '~/Documents/PhD/2_EM_Fire_effect/figures_output/'
res.dir <- "~/Documents/PhD/2_EM_Fire_effect/results_output/"

#----------------------------------------------------------------------------------------#
# Soil data----
#----------------------------------------------------------------------------------------#

soil.data <- read.csv('data_output/soil_data.csv')

# #--Soil data
# soil.data <- read.csv(paste0(dat.dir, 'soil_data.csv'),as.is = T)
# 
# #--change < 1.0 in the nitrate column to 0. 
# soil.data[soil.data$no3.n.ppm == '<1.0', 'no3.n.ppm'] <- 0
# soil.data$no3.n.ppm <- as.numeric(soil.data$no3.n.ppm)
# 
# #--add burn status of site
# burn <- c('sc3','sc4','p3','p4')
# soil.data[soil.data$site %in% burn, 'burn_status'] <- 'burned'
# soil.data[!soil.data$site %in% burn, 'burn_status'] <- 'unburned'
# 
# #--add lat and long coordinates
# site.data <- read.csv(paste0(dat.dir,'site_data.csv'), as.is = T)
# for(s in unique(soil.data$site)){
#   soil.data[soil.data$site == s, 'lat'] <- unique(site.data[site.data$site == s, 'lat'])
#   soil.data[soil.data$site == s, 'long'] <- unique(site.data[site.data$site == s, 'long'])
# }
# 
# #--add range information
# soil.data[!soil.data$site %in% c('sc1','sc2','sc3','sc4'), 'range'] <- 'Pinaleno Mts.'
# soil.data[soil.data$site %in% c('sc1','sc2','sc3','sc4'), 'range'] <-
#   'Santa Catalina Mts.'
# 
# #--add forest type
# soil.data$forest <- NA
# soil.data[soil.data$site %in% c('sc1','sc2','sc3','sc4','p2'), 'forest'] <- 'pine oak'
# soil.data[soil.data$site %in% c('p3','p4'), 'forest'] <- 'pine df'
# soil.data[soil.data$site == 'p1', 'forest'] <- 'pine'
# 
# ## PCA of soil data -------
# # only unburned data
# soil.data <- soil.data[soil.data$burn_status == 'unburned',]
# soil.sig <- c('pH.su','po4.p.ppm','mg.ppm','k.ppm','no3.n.ppm')
# soil.perma <- soil.data[colnames(soil.data) %in% soil.sig]
# soil.perma[soil.perma$no3.n.ppm == "<1.0", 'no3.n.ppm'] <- 0.00000001
# soil.perma$no3.n.ppm <- as.numeric(soil.perma$no3.n.ppm)
# soil.pca <- prcomp(soil.perma,
#                    center = TRUE,
#                    scale. = TRUE) 
# #print(soil.pca)
# #plot(soil.pca, type = 'l')
# summary(soil.pca)
# soil.eigenvector <- scores(soil.pca, choices = 1)
# soil.eigenvector <- data.frame(matrix(unlist(soil.eigenvector), nrow=12, byrow=T))
# colnames(soil.eigenvector) <- 'soil.pca'
# soil.data['soil.pca'] <- soil.eigenvector

#----------------------------------------------------------------------------------------#
# Community, climate, and soil data (extrapolated) data----
#----------------------------------------------------------------------------------------#

stsp.matrix <- read.csv('data_output/97%_SitexSpecies_TipAb.csv')
clim.data <- read.csv('data/climate_data.csv')
# #--Load data
# all.data <- read.csv(paste0(dat.dir,'20170806_OTU_data.csv'), as.is = T)
# 
# #--Remove sequences not assigned to an OTU
# otu.data <- all.data[!is.na(all.data$otu.97), ]
# 
# #--Remove sequences not assigned to Ponderosa host
# otu.data <- otu.data[otu.data$Host == 'Ponderosa',]
# 
# #--Add an OTU count column (1 sequence = 1 frequency count)
# otu.data['otu.count'] <- 1
# 
# #--Creates matrix, grouping OTUs by tree number and getting a sum of EM tip abundance 
# # for each OTU
# otu.data %>%
#   select(Sample_name,Burn_status,Range,Site,Tree,otu.97,Tip_count) %>%
#   spread(otu.97,Tip_count) %>%
#   group_by(Tree,Site,Range,Burn_status) %>%
#   select(matches ("otu*")) %>%
#   summarize_each(funs(sum(., na.rm = TRUE))) %>%
#   as.data.frame() -> stsp.matrix
# 
# #--climate data: add to stsp.matrix
# clim.data <- read.csv(paste0(dat.dir, 'climate_data.csv'))
# for(i in unique(stsp.matrix$Site)) {
#   stsp.matrix[stsp.matrix$Site == i, 'Prec.avg'] <- clim.data[clim.data$site == i, 'Prec.avg']
#   stsp.matrix[stsp.matrix$Site == i, 'Temp.avg'] <-
#     clim.data[clim.data$site == i, 'Temp.avg']
# }
# 
# #--reorder matrix
# colnames(stsp.matrix)
# stsp.matrix <- stsp.matrix[c(1:4,121:122,5:120)]

#========================================================================================#
# Assess soil and climate between unburned sites in SCM and PM ranges----
#========================================================================================#

#<<Multiple regression of soil PCA axis one (soil.eigenvector)>>------
## PCA of soil data -------
# only unburned data
soil.sig <- c('pH.su','po4.p.ppm','Log.mg','Log.k','Log.mn','so4.s.ppm')
soil.perma <- soil.data[colnames(soil.data) %in% soil.sig]
soil.pca <- prcomp(soil.perma,
                   center = TRUE,
                   scale. = TRUE)
#print(soil.pca)
#plot(soil.pca, type = 'l')
summary(soil.pca)
soil.eigenvector <- scores(soil.pca, choices = 1)
soil.eigenvector <- data.frame(matrix(unlist(soil.eigenvector), nrow=24, byrow=T))
colnames(soil.eigenvector) <- 'soil.pca'
soil.data['soil.pca'] <- soil.eigenvector

#--Anova to assess differences between range and fire history
soil.lm <- lm(soil.pca ~ range * fire.history, data = soil.data)
summary(soil.lm)

#--plot
range.soil <- ggplot(soil.data, aes(x = fire.history,
                             y = soil.pca)) +
  geom_boxplot() +
  theme_bw() +
  facet_grid(~ range) +
  ylab('PCA of soil variables') +
  xlab(element_blank()) +
  theme(legend.position="none",
        axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size=22, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

ggsave('SoilPCA_UnburnedData_Range.jpeg', plot = range.soil,
       device = 'jpeg', path = 'figures_output/',
       width = 20, height = 15, units = 'cm')

#========================================================================================#
# Climate data----
#========================================================================================#

#<<Assess co-correlation of soil data>>-------
clim.lm <- lm(Prec.avg ~ Temp.avg, data = clim.data)
summary(clim.lm)

#<<T-test and Plot of precipitation>>-------------
#--Range
prec.range.t <- kruskal.test(Prec.avg ~ range, data = clim.data)
prec.range.t

# modify range names
levels(clim.data$range) <- c('Pinaleno Mts','Santa Catalina Mts')

prec.range <- ggplot(clim.data, aes(x = range,
                             y = Prec.avg)) +
  geom_boxplot() +
  theme_bw() +
  ylab('Annual average\nprecipitation (mm)') +
  xlab('Range') +
  theme(legend.position="none",
        axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size=22, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

prec.range

ggsave('Precipitation_Range.jpeg', plot = prec.range,
       device = 'jpeg', path = 'figures_output/',
       width = 20, height = 15, units = 'cm')

#========================================================================================#
# Check for correlation between sig. soil characateristics----
#========================================================================================#

#--pH and Phosphate: correlated
plot(soil.data$pH.su, soil.data$po4.p.ppm)
ph.po4 <- lm(soil.data$po4.p.ppm ~ soil.data$pH.su)
abline(ph.po4)
anova(ph.po4)

#--pH and sulfate:correlated
plot(soil.data$pH.su, soil.data$so4.s.ppm)
ph.so4 <- lm(soil.data$so4.s.ppm ~ soil.data$pH.su)
abline(ph.so4)
anova(ph.so4)

#--pH and boron
plot(soil.data$pH.su, soil.data$b.ppm)
ph.b <- lm(soil.data$b.ppm ~ soil.data$pH.su)
abline(ph.b)
anova(ph.b)

#--Phosphate and sulfate
plot(soil.data$po4.p.ppm, soil.data$so4.s.ppm)
po4.so4 <- lm(soil.data$so4.s.ppm ~ soil.data$po4.p.ppm)
abline(po4.so4)
anova(po4.so4)

#--Phosphate and boron
plot(soil.data$po4.p.ppm, soil.data$b.ppm)
po4.b <- lm(soil.data$b.ppm ~ soil.data$po4.p.ppm)
abline(po4.b)
anova(po4.b)

#--Sulfate and boron: correlated
plot(soil.data$so4.s.ppm, soil.data$b.ppm)
so4.b <- lm(soil.data$b.ppm ~ soil.data$so4.s.ppm)
abline(so4.b)
anova(so4.b)
