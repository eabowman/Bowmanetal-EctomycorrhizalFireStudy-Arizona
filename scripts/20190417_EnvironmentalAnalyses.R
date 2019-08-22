## Script created by Liz Bowman July 19, 2017
## for analyzing soil differences between sites, ranges, and burn

#========================================================================================#
# Load data and libraries----
#========================================================================================#

#----------------------------------------------------------------------------------------#
# Load libraries----
#----------------------------------------------------------------------------------------#
library(ggplot2);library (tidyr);library (vegan);library (dplyr); library(gridExtra)
#install.packages('ggfortify')
library(ggfortify); library(plyr); library(nlme)

#----------------------------------------------------------------------------------------#
# set up paths to directories----
#----------------------------------------------------------------------------------------#
#--path to directory 
dat.dir <- "~/Documents/PhD/2_EM_Fire_effect/data/"
fig.dir <- '~/Documents/PhD/2_EM_Fire_effect/figures_output/'
res.dir <- "~/Documents/PhD/2_EM_Fire_effect/results_output//"

#----------------------------------------------------------------------------------------#
# Soil data----
#----------------------------------------------------------------------------------------#

#--Soil data
soil.data <- read.csv(paste0(dat.dir, 'soil_data.csv'),as.is = T)

#--change < 1.0 in the nitrate column to 0. 
soil.data[soil.data$no3.n.ppm == '<1.0', 'no3.n.ppm'] <- 0
soil.data$no3.n.ppm <- as.numeric(soil.data$no3.n.ppm)

#--add burn status of site
burn <- c('sc3','sc4','p3','p4')
soil.data[soil.data$site %in% burn, 'burn_status'] <- 'burned'
soil.data[!soil.data$site %in% burn, 'burn_status'] <- 'unburned'

#--add lat and long coordinates
site.data <- read.csv(paste0(dat.dir,'site_data.csv'), as.is = T)
for(s in unique(soil.data$site)){
  soil.data[soil.data$site == s, 'lat'] <- unique(site.data[site.data$site == s, 'lat'])
  soil.data[soil.data$site == s, 'long'] <- unique(site.data[site.data$site == s, 'long'])
}

#--add range information
soil.data[!soil.data$site %in% c('sc1','sc2','sc3','sc4'), 'range'] <- 'Pinaleno Mts.'
soil.data[soil.data$site %in% c('sc1','sc2','sc3','sc4'), 'range'] <-
  'Santa Catalina Mts.'

#--add column with range and burn status combined
soil.data$range_burn <- NA
for(i in unique(soil.data$range)){
  for(b in unique(soil.data$burn_status)){
    soil.data[soil.data$range == i & soil.data$burn_status == b, 'range_burn'] <-
      paste(i, b)
  }
}

#--add forest type
soil.data$forest <- NA
soil.data[soil.data$site %in% c('sc1','sc2','sc3','sc4','p2'), 'forest'] <- 'pine oak'
soil.data[soil.data$site %in% c('p3','p4'), 'forest'] <- 'pine df'
soil.data[soil.data$site == 'p1', 'forest'] <- 'pine'

## PCA of soil data -------
# all soil data
# soil.perma <- soil.data[6:22] 
# sig soil data
soil.sig <- c('pH.su','po4.p.ppm','so4.s.ppm','b.ppm','mg.ppm',
              'k.ppm','no3.n.ppm')
soil.perma <- soil.data[colnames(soil.data) %in% soil.sig]
soil.perma[soil.perma$no3.n.ppm == "<1.0", 'no3.n.ppm'] <- 0.00000001
soil.perma$no3.n.ppm <- as.numeric(soil.perma$no3.n.ppm)
soil.pca <- prcomp(soil.perma,
                   center = TRUE,
                   scale. = TRUE) 
#print(soil.pca)
#plot(soil.pca, type = 'l')
summary(soil.pca)
soil.eigenvector <- scores(soil.pca, choices = 1)
soil.eigenvector <- data.frame(matrix(unlist(soil.eigenvector), nrow=24, byrow=T))
colnames(soil.eigenvector) <- 'soil.pca'

#----------------------------------------------------------------------------------------#
# Community, climate, and soil data (extrapolated) data----
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
soil.sig <- c('site','pH.su','po4.p.ppm','so4.s.ppm','b.ppm','mg.ppm',
              'k.ppm','no3.n.ppm')
# all.soil <- c('site','lat','long','pH.su','po4.p.ppm','so4.s.ppm','b.ppm','mg.ppm',
#               'k.ppm','no3.n.ppm','EC.ds.m','ca.ppm','na.ppm','zn.ppm','fe.ppm','mn.ppm',
#               'cu.ppm','ni.ppm','cec.meq.100g')
soil.data.sig <- soil.data[colnames(soil.data)%in% soil.sig]
soil.data.sig[soil.data.sig$no3.n.ppm == '<1.0', 'no3.n.ppm'] <- 0
soil.data.sig$no3.n.ppm <- as.numeric(soil.data.sig$no3.n.ppm)
for(s in unique(soil.data$site)){
  for(f in soil.sig[2:length(soil.sig)]){
    stsp.matrix[stsp.matrix$Site == s, f] <- mean(soil.data.sig[soil.data.sig$site == s, f])
  }
}

#--reorder matrix
colnames(stsp.matrix)
soil.perm <- stsp.matrix[c(1:4,123:129)]
stsp.matrix <- stsp.matrix[c(1:4,121:122,5:120)]

#========================================================================================#
# Assess soil and climate between SCM and PM ranges----
#========================================================================================#

#<<Multiple regression of soil PCA axis one (soil.eigenvector)>>------
# using qualitative variables --> need to create dummy variables
#--Add range and burn_status to soil.eigenvector df
soil.eigenvector$site <- as.factor(soil.data$site)
soil.eigenvector$site <- soil.data$site
soil.eigenvector$tree <- soil.data$tree_number
soil.eigenvector$range <- soil.data$range
soil.eigenvector$burn_status <- soil.data$burn_status
soil.eigenvector$forest <- soil.data$forest
#--Create dummy variables for the predictor variables
soil.eigenvector$range_dum <- as.factor(revalue(soil.eigenvector$range,
                                            c('Pinaleno Mts.' = '1', 'Santa Catalina Mts.' = '2')))
soil.eigenvector$burn_status_dum <- as.factor(revalue(soil.eigenvector$burn_status,
                                                  c('burned' = '1', 'unburned' = '2')))

#--Check to see if a mixed model is needed from "Discovering statistics using R" (p879)
# assess base model
interceptOnly <- gls(soil.pca ~ 1, data = soil.eigenvector, method = 'ML')
summary(interceptOnly)
# add in site as a random variable, assess whether it improves the model
randomInterceptOnly <- lme(soil.pca ~ 1, data = soil.eigenvector,
                           random = ~1|site, method = 'ML')
summary(randomInterceptOnly)
anova(interceptOnly, randomInterceptOnly) # use random effects model
# add range as a fixed variable to the model
randomInterceptRange <- lme(soil.pca ~ range, data = soil.eigenvector,
                            random = ~ 1 | site/tree, method = 'ML')
summary(randomInterceptRange)
# add burn status as a fixed variable
randomInterceptRangeBurn <- lme(soil.pca ~ range * burn_status, data = soil.eigenvector,
                                random = ~ 1 | site/tree, method = 'ML')
summary(randomInterceptRangeBurn)
# assess all models 
anova(randomInterceptOnly, randomInterceptRange, randomInterceptRangeBurn)
### The range only model has the best fit (based on AIC and p-value)

#--Linear multiple regression
lme.soil <- lme(soil.pca ~ range * burn_status,
                data = soil.eigenvector,
                random = ~ 1 | site/tree)
summary(lme.soil)
anova(lme.soil)

#<<Check model>>----
#--Plot random effects, should be random around 0 intercept
plot(ranef(lme.soil))

#--Plot residuals to check for heterstasdicity
plot(lme.soil)

#--plot
soil.eigenvector$range <- factor(soil.eigenvector$range)
levels(soil.eigenvector$range) <- c('Pinaleno Mts.','Santa Catalina Mts')
soil.eigenvector$burn_status <- factor(soil.eigenvector$burn_status)
levels(soil.eigenvector$burn_status) <-  c('FA','FU')

rangeburn.soil <- ggplot(soil.eigenvector, aes(x = burn_status,
                             y = soil.pca)) +
  geom_boxplot() +
  theme_bw() +
  facet_grid(~ range) +
  ylab('PCA of soil variables') +
  xlab('Burn history') +
  theme(legend.position="none",
        axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size=22, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

ggsave('SoilPCA_BurnRange.jpeg', plot = rangeburn.soil,
       device = 'jpeg', path = 'figures_output/',
       width = 20, height = 15, units = 'cm')


#<<PCA of soil data>>-------------
soil.perma <- soil.data[6:22]
soil.pca <- prcomp(soil.perma,
                   center = TRUE,
                   scale. = TRUE)
autoplot(soil.pca,
         data = soil.data, colour = 'range_burn',
         frame = T)
autoplot(soil.pca,
         data = soil.data, colour = 'range',
         loadings = T, loadings.label = T,
         frame = T)
autoplot(soil.pca,
         data = soil.data, colour = 'burn_status',
         frame = T)
autoplot(soil.pca,
         data = soil.data, colour = 'forest',
         frame = T)

#--Cluster analysis of soil----
save.soil <- soil.data
#soil.data <- save.soil
#--distance matrix of soil data
row.names(soil.data) <- paste(row.names(soil.data), soil.data$burn_status)
soil.dist <- dist(soil.data[6:22], method = 'euclidean')

clus <- hclust(soil.dist, 'average')
plot(clus)
range(soil.dist)
cor(soil.dist, cophenetic(clus))

plot(clus)
rect.hclust(clus, 2)
grp <- cutree(clus, 2)

#--Correlation of soil traits (as above) and climate-----
pc <- prcomp(soil.perm[5:length(soil.perm)], scale = T)
pc <- scores(pc, display = 'sites', choices = 1:4)
edis <- vegdist(pc, method = 'euclid')
clim.pc <- dist(stsp.matrix[c('lat','long')], method = 'euclidean')
#clim.pc <- vegdist(wisconsin(sqrt(stsp.matrix[9:length(stsp.matrix)])))
plot(clim.pc, edis)
mantel(clim.pc, edis)

#========================================================================================#
# Climate data----
#========================================================================================#
#--Create matrix of nonduplicated climate data, i.e. by site
stsp.matrix %>%
  select(Site, Burn_status, Range, prec, Tmax) %>%
  distinct(Site, .keep_all = T) -> clim.data

#<<Assess co-correlation of soil data>>-------
clim.lm <- lm(prec ~ Tmax, data = clim.data)
summary(clim.lm)

#<<PCA of climate data: not used, exploratory>>-------------
climate.perma <- stsp.matrix[c('prec','Tmax')]
climate.pca <- prcomp(climate.perma,
                      center = TRUE,
                      scale. = TRUE)
autoplot(climate.pca,
         data = stsp.matrix, colour = 'Range',
         loadings = T, loadings.label = T,
         frame = T)

autoplot(climate.pca,
         data = stsp.matrix, colour = 'Burn_status',
         loadings = T, loadings.label = T,
         frame = T)

#<<T-test and Plot of precipitation>>-------------
#--Range
prec.range.t <- t.test(prec ~ Range, data = clim.data)

prec.range <- ggplot(clim.data, aes(x = Range,
                             y = prec)) +
  geom_boxplot() +
  theme_bw() +
  ylab('Annual average precipitation (mm)') +
  xlab('Range') +
  theme(legend.position="none",
        axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size=22, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

tmax.range.t <- t.test(Tmax ~ Range, data = clim.data)

tmax.range <- ggplot(clim.data, aes(x = Range,
                        y = Tmax)) +
  geom_boxplot() +
  theme_bw() +
  ylab('Max annual temperature (°C)') +
  xlab('Range') +
  theme(legend.position="none",
        axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size=22, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

#--Burn_status
prec.burn.t <- t.test(prec ~ Burn_status, data = clim.data)

prec.burn <- ggplot(clim.data, aes(x = Burn_status,
                                      y = prec)) +
  geom_boxplot() +
  theme_bw() +
  ylab('Annual average precipitation (mm)') +
  xlab('Burn status') +
  theme(legend.position="none",
        axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size=22, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

tmax.burn.t <- t.test(Tmax ~ Burn_status, data = clim.data)

tmax.burn <- ggplot(clim.data, aes(x = Burn_status,
                                      y = Tmax)) +
  geom_boxplot() +
  theme_bw() +
  ylab('Max annual temperature (°C)') +
  xlab('Burn status') +
  theme(legend.position="none",
        axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size=22, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

#<<multiple linear regression of climate>>-------------
#--Check to see if a mixed model is needed from "Discovering statistics using R" (p879)
# assess base model
interceptOnly <- gls(prec ~ 1, data = clim.data, method = 'ML')
summary(interceptOnly)
# add in site as a random variable, assess whether it improves the model
randomInterceptOnly <- lme(prec ~ 1, data = clim.data,
                           random = ~1|Site, method = 'ML')
summary(randomInterceptOnly)
anova(interceptOnly, randomInterceptOnly) # use random effects model
# add range as a fixed variable to the model
randomInterceptRange <- lme(prec ~ Range, data = clim.data,
                            random = ~ 1 | Site, method = 'ML')
summary(randomInterceptRange)
# add burn status as a fixed variable
randomInterceptRangeBurn <- lme(prec ~ Range * Burn_status, data = clim.data,
                                random = ~ 1 | Site, method = 'ML')
summary(randomInterceptRangeBurn)
# assess all models 
anova(randomInterceptOnly, randomInterceptRange, randomInterceptRangeBurn)
### The range only model has the best fit (based on AIC and p-value)

#--Linear multiple regression
lme.clim <- lme(prec ~ Range * Burn_status,
                data = clim.data,
                random = ~ 1 | Site)
summary(lme.clim)
anova(lme.clim)

#<<Check model>>----
#--Plot random effects, should be random around 0 intercept
plot(ranef(lme.clim))

#--Plot residuals to check for heterstasdicity
plot(lme.clim)

clim.data$Range <- factor(clim.data$Range)
levels(clim.data$Range) <- c('Pinaleno Mts.','Santa Catalina Mts')
clim.data$Burn_status <- factor(clim.data$Burn_status)
levels(clim.data$Burn_status) <-  c('FA','FU')

rangeburn.prec <- ggplot(clim.data, aes(x = Burn_status,
                      y = prec)) +
  geom_boxplot() +
  theme_bw() +
  facet_grid(~ Range)+
  ylab('Annual average\nprecipitation (mm)') +
  xlab('Burn history') +
  theme(legend.position="none",
        axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size=22, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

ggsave('Precipitation_BurnRange.jpeg', plot = rangeburn.prec,
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
