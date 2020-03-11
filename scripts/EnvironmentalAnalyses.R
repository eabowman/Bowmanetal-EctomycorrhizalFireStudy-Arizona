## Script created by Liz Bowman July 19, 2017
## for analyzing soil and climate data

#========================================================================================#
# Load data ----
#========================================================================================#

soil.data <- read.csv(paste0(dat.dir,'soil_data.csv'))
stsp.matrix <- read.csv(paste0(dat.dir,'97%_SitexSpecies_TipAb.csv'))
clim.data <- read.csv(paste0(dat.dir,'climate_data.csv'))

#========================================================================================#
# Assess soil and climate between unburned sites in SCM and PM ranges----
#========================================================================================#

#<<Multiple regression of soil PCA axis one (soil.eigenvector)>>------

## PCA of soil data -------
# only unburned data
soil.sig <- c('pH.su','po4.p.ppm','Log.mg','Log.mn','so4.s.ppm','Log.k')
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

#<<Mixed effects model between ranges (random) and fire history >>-------------
soil.model <- aov(soil.pca ~ fire.history * range, data = soil.data)
summary(soil.model)

#--plot
range.soil <- ggplot(soil.data, aes(x = fire.history,
                             y = soil.pca)) +
  geom_boxplot() +
  facet_grid(. ~ range) +
  theme_bw() +
  ylab('PCA of soil variables') +
  xlab(element_blank()) +
  theme(legend.position="none",
        axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size=22, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

ggsave('SoilPCA_AllData_Range.jpeg', plot = range.soil,
       device = 'jpeg', path = 'figures_output/',
       width = 20, height = 15, units = 'cm')

#========================================================================================#
# Climate data----
#========================================================================================#

#<<Assess co-correlation of soil data>>-------
clim.lm <- lm(Prec.avg ~ Temp.avg, data = clim.data)
summary(clim.lm)

#<< linear mixed effects model (prec) ranges (random) and fire history >>-------------
prec.model <- lme(Prec.avg ~ burn_status, data = clim.data, random = ~ 1 | range)
summary(prec.model)
plot(prec.model)

#<< linear mixed effects model (temp) ranges (random) and fire history >>-------------
temp.model <- lme(Temp.avg ~ burn_status, data = clim.data, random = ~ 1 | range)
summary(temp.model)

#<< T-test as function of fire history w/in each range >>-------------
#---Precipitation
# SCM
clim.scm <- clim.data[clim.data$range == 'santa.catalina',]
ttest.scm <- t.test(Prec.avg ~ burn_status, data = clim.scm)
ttest.scm

# PM
clim.pm <- clim.data[clim.data$range == 'pinaleno',]
ttest.pm <- t.test(Prec.avg ~ burn_status, data = clim.pm)
ttest.pm

#---Temperature
# SCM
clim.scm <- clim.data[clim.data$range == 'santa.catalina',]
ttest.scm <- t.test(Temp.avg ~ burn_status, data = clim.scm)
ttest.scm

# PM
clim.pm <- clim.data[clim.data$range == 'pinaleno',]
ttest.pm <- t.test(Temp.avg ~ burn_status, data = clim.pm)
ttest.pm

#<< Plot of precipitation between ranges >>-------------
# modify range names
levels(clim.data$range) <- c('Pinaleno Mts','Santa Catalina Mts')
levels(clim.data$burn_status) <- c('Burned', 'Unburned')

prec.range <- ggplot(clim.data, aes(x = burn_status,
                             y = Prec.avg)) +
  geom_boxplot() +
  facet_grid(. ~ range) +
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
       device = 'jpeg', path = fig.dir,
       width = 20, height = 15, units = 'cm')
