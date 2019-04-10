## Script created by Liz Bowman January 22, 2018
## for analyzing soil and climate differences between range and burn status of sites

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
env.matrix <- stsp.matrix[c(1:4,121:length(stsp.matrix))]

source('~/Documents/PhD/2_EM_Fire_effect/scripts/functions.R')

#========================================================================================#
# Load data and libraries----
#========================================================================================#

responses <- c("pH.su","EC.ds.m","ca.ppm","mg.ppm","na.ppm","k.ppm","zn.ppm","fe.ppm",
              "mn.ppm","cu.ppm","ni.ppm","no3.n.ppm","po4.p.ppm","so4.s.ppm","b.ppm",
              "cec.meq.100g")
ln.env.matrix <- env.matrix[1:6]
for(i in responses){
  ln.env.matrix[,i] <- ln(env.matrix[,i])
}

plot(prec)





