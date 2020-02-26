## Script created by Liz Bowman June 3, 2017
## for analyzing climatic differences

#=========================================================================================
# Load data and libraries----
#=========================================================================================

#-----------------------------------------------------------------------------------------
# Load libraries----
#-----------------------------------------------------------------------------------------
library(ggplot2);library (tidyr);library (vegan);library (dplyr); library(gridExtra)

#-----------------------------------------------------------------------------------------
# set up paths to directories----
#-----------------------------------------------------------------------------------------
#--path to directory where the climate data downloaded from BIOCLIM site is stored
dat.dir <- "~/Documents/PhD/2_EM_Fire_effect/data/"
fig.dir <- '~/Documents/PhD/2_EM_Fire_effect/figures_output/'
res.dir <- "~/Documents/PhD/2_EM_Fire_effect/results_output/"

#----------------------------------------------------------------------------------------#
# Load data and clean up----
#----------------------------------------------------------------------------------------#

clim.data <- read.csv(paste0(dat.dir, 'climate_data.csv'))


#--add range information
clim.data[!clim.data$site %in% c('sc1','sc2','sc3','sc4'), 'range'] <- 'Pinaleno Mts.'
clim.data[clim.data$site %in% c('sc1','sc2','sc3','sc4'), 'range'] <-
  'Santa Catalina Mts.'

#========================================================================================#
# Co-corrrelation of climatic variables----
#========================================================================================#

#--Tmax, Prec: correlated
tmax.prec.lm <- lm(clim.data$Tmax ~ clim.data$prec)
summary(tmax.prec.lm)
plot(clim.data$prec, clim.data$Tmax)

#--Tmin, Prec: correlated
tmin.prec.lm <- lm(clim.data$Tmin ~ clim.data$prec)
summary(tmin.prec.lm)
plot(clim.data$prec, clim.data$Tmin)

#--Tmin, Tmax: correlated (plot)
tmin.tmax.lm <- lm(clim.data$Tmin ~ clim.data$Tmax)
summary(tmin.tmax.lm)
plot(clim.data$Tmax, clim.data$Tmin)

#========================================================================================#
# Differences between range----
#========================================================================================#
source('~/Documents/PhD/2_EM_Fire_effect/scripts/functions.R')

responses <- c('Tmax','Tmin','prec')
ttester(clim.data, 'range', responses)
wilcoxtest(clim.data, 'range', responses)

#<< plots of climatic factors >> ---------------------------------------------------------
# Max temperature---- normal
t.test(clim.data$Tmax ~ clim.data$range)
normtest(clim.data, 'Tmax')
warm <- ggplot(clim.data, aes(x = range, y = Tmax, fill = range)) +
  geom_boxplot() +
  scale_x_discrete(name = "Range") +
  scale_y_continuous(name = "Avg. annual high (°C)") +
  #facet_grid(. ~ range) +
  theme_bw() +
  scale_fill_brewer(palette = "Accent") +
  labs(fill = "Range") +
  theme(legend.position="none",
        axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size=22, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

# avg. temp coldest month---- normal
normtest(clim.data, 'Tmin')
cold <- ggplot(clim.data, aes(x = range, y = Tmin, fill = range)) +
  geom_boxplot() +
  scale_x_discrete(name = "Range") +
  scale_y_continuous(name = "Avg. annual low (°C)") +
  #facet_grid(. ~ range) +
  theme_bw() +
  scale_fill_brewer(palette = "Accent") +
  labs(fill = "Range") +
  theme(legend.position="none",
        axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size=22, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

# Precipitation---- normal
t.test(clim.data$prec ~ clim.data$burn_status)
normtest(clim.data, 'prec')
prec <- ggplot(clim.data, aes(x = range, y = prec, fill = range)) +
  geom_boxplot() +
  scale_x_discrete(name = "Range") +
  scale_y_continuous(name = "Annual avg. \n precipitation (mm)") +
  #facet_grid(. ~ range) +
  theme_bw() +
  scale_fill_manual(values=c("#31A354", "#E5F5E0")) +
  labs(fill = "Range") +
  theme(legend.position="none",
        axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size=20, color = 'black'),
        axis.title = element_text(size = 22),
        strip.text.x = element_text(size = 14))

all.plot <- grid.arrange(warm,cold,prec,
                         nrow = 1, ncol = 3)
ggsave('Prec.pdf', plot = prec, device = 'pdf', path = fig.dir,
       width = 7, height = 5)
