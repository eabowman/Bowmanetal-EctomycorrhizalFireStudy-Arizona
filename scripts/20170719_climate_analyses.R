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
fig.dir <- '~/Documents/PhD/2_EM_Fire_effect/figures/'
res.dir <- "~/Documents/PhD/2_EM_Fire_effect/results/"

#-----------------------------------------------------------------------------------------
# Load data and clean up----
#-----------------------------------------------------------------------------------------

clim.data <- read.csv(paste0(dat.dir, 'climate_data.csv'))


#--add range information
clim.data[!clim.data$site %in% c('sc1','sc2','sc3','sc4'), 'range'] <- 'Pinaleno Mts.'
clim.data[clim.data$site %in% c('sc1','sc2','sc3','sc4'), 'range'] <-
  'Santa Catalina Mts.'

#=========================================================================================
# Differences between range----
#=========================================================================================
source('~/Documents/PhD/2_EM_Fire_effect/scripts/functions.R')

response <- c('BIO10_red','BIO11_red','prec')
ttester(clim.data, 'range', responses)
wilcoxtest(clim.data, 'range', responses)

#<< plots of climatic factors >> ---------------------------------------------------------
# avg. temp warmest month---- normal
normtest(clim.data, 'BIO10_red')
warm <- ggplot(clim.data, aes(x = range, y = BIO10_red, fill = range)) +
  geom_boxplot() +
  scale_x_discrete(name = "Range") +
  scale_y_continuous(name = "Avg. temp. warmest quarter") +
  #facet_grid(. ~ range) +
  theme_bw() +
  scale_fill_brewer(palette = "Accent") +
  labs(fill = "Range") +
  theme(legend.position="none")

# avg. temp coldest month---- normal
normtest(clim.data, 'BIO11_red')
cold <- ggplot(clim.data, aes(x = range, y = BIO11_red, fill = range)) +
  geom_boxplot() +
  scale_x_discrete(name = "Range") +
  scale_y_continuous(name = "Avg. temp. coldest quarter") +
  #facet_grid(. ~ range) +
  theme_bw() +
  scale_fill_brewer(palette = "Accent") +
  labs(fill = "Range") +
  theme(legend.position="none")

# Precipitation---- normal
normtest(clim.data, 'prec')
prec <- ggplot(clim.data, aes(x = range, y = prec, fill = range)) +
  geom_boxplot() +
  scale_x_discrete(name = "Range") +
  scale_y_continuous(name = "Precipitation") +
  #facet_grid(. ~ range) +
  theme_bw() +
  scale_fill_brewer(palette = "Accent") +
  labs(fill = "Range") +
  theme(legend.position="none")

all.plot <- grid.arrange(warm,cold,prec,
                         nrow = 1, ncol = 3)
ggsave('clim.pdf', plot = all.plot, device = 'pdf', path = fig.dir,
       width = 10, height = 5)
