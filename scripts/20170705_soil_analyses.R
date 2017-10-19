## Script created by Liz Bowman July 19, 2017
## for analyzing soil differences between sites, ranges, and burn

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
#--path to directory 
dat.dir <- "~/Documents/PhD/3_EM_Fire_effect/data/"
fig.dir <- '~/Documents/PhD/3_EM_Fire_effect/figures/'
res.dir <- "~/Documents/PhD/3_EM_Fire_effect/results/"

#-----------------------------------------------------------------------------------------
# Load data and clean up----
#-----------------------------------------------------------------------------------------

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

#=========================================================================================
# Anova of soil traits differences based on both burn status and range----
#=========================================================================================

response <- c("pH.su","EC.ds.m","ca.ppm","mg.ppm","na.ppm","k.ppm","zn.ppm","fe.ppm",
               "mn.ppm","cu.ppm","ni.ppm","no3.n.ppm","po4.p.ppm","so4.s.ppm","b.ppm",
               "esp","cec.meq.100g")
anova.rb <- anova.t(soil.data, 'range_burn', response)
anova.rb$p <- round(anova.rb$p, 5)

#<< plots of soil factors >> -------------------------------------------------------------
# pH----
normtest(soil.data, 'pH.su')
ph <- ggplot(soil.data, aes(x = burn_status, y = pH.su, fill = burn_status)) +
  geom_boxplot() +
  scale_x_discrete(name = "Burn status") +
  scale_y_continuous(name = "pH (su)") +
  facet_grid(. ~ range) +
  theme_bw() +
  scale_fill_brewer(palette = "Accent") +
  labs(fill = "Burn status") +
  theme(legend.position="none")

# EC---- not normal
normtest(soil.data, 'EC.ds.m')
ec <- ggplot(soil.data, aes(x = burn_status, y = EC.ds.m, fill = burn_status)) +
  geom_boxplot() +
  scale_x_discrete(name = "Burn status") +
  scale_y_continuous(name = "EC (ds/m)") +
  facet_grid(. ~ range) +
  theme_bw() +
  scale_fill_brewer(palette = "Accent") +
  labs(fill = "Burn status") +
  theme(legend.position="none")

# Ca----not normal
normtest(soil.data, 'ca.ppm')
ca <- ggplot(soil.data, aes(x = burn_status, y = ca.ppm, fill = burn_status)) +
  geom_boxplot() +
  scale_x_discrete(name = "Burn status") +
  scale_y_continuous(name = "Ca (ppm)") +
  facet_grid(. ~ range) +
  theme_bw() +
  scale_fill_brewer(palette = "Accent") +
  labs(fill = "Burn status") +
  theme(legend.position="none")

# Mg----not normal
normtest(soil.data, 'mg.ppm')
mg <- ggplot(soil.data, aes(x = burn_status, y = mg.ppm, fill = burn_status)) +
  geom_boxplot() +
  scale_x_discrete(name = "Burn status") +
  scale_y_continuous(name = "Mg (ppm)") +
  facet_grid(. ~ range) +
  theme_bw() +
  scale_fill_brewer(palette = "Accent") +
  labs(fill = "Burn status") +
  theme(legend.position="none")

# Na---- not normal
normtest(soil.data, 'na.ppm')
na <- ggplot(soil.data, aes(x = burn_status, y = na.ppm, fill = burn_status)) +
  geom_boxplot() +
  scale_x_discrete(name = "Burn status") +
  scale_y_continuous(name = "Na (ppm)") +
  facet_grid(. ~ range) +
  theme_bw() +
  scale_fill_brewer(palette = "Accent") +
  labs(fill = "Burn status") +
  theme(legend.position="none")

# K----not normal
normtest(soil.data, 'k.ppm')
k <- ggplot(soil.data, aes(x = burn_status, y = k.ppm, fill = burn_status)) +
  geom_boxplot() +
  scale_x_discrete(name = "Burn status") +
  scale_y_continuous(name = "K (ppm)") +
  facet_grid(. ~ range) +
  theme_bw() +
  scale_fill_brewer(palette = "Accent") +
  labs(fill = "Burn status") +
  theme(legend.position="none")

# Zn----not normal
normtest(soil.data, 'zn.ppm')
zn <- ggplot(soil.data, aes(x = burn_status, y = zn.ppm, fill = burn_status)) +
  geom_boxplot() +
  scale_x_discrete(name = "Burn status") +
  scale_y_continuous(name = "Zn (ppm)") +
  facet_grid(. ~ range) +
  theme_bw() +
  scale_fill_brewer(palette = "Accent") +
  labs(fill = "Burn status") +
  theme(legend.position="none")

# Fe----not normal
normtest(soil.data, 'fe.ppm')
fe <- ggplot(soil.data, aes(x = burn_status, y = fe.ppm, fill = burn_status)) +
  geom_boxplot() +
  scale_x_discrete(name = "Burn status") +
  scale_y_continuous(name = "Fe (ppm)") +
  facet_grid(. ~ range) +
  theme_bw() +
  scale_fill_brewer(palette = "Accent") +
  labs(fill = "Burn status") +
  theme(legend.position="none")

# Mn----not normal
normtest(soil.data, 'mn.ppm')
mn <- ggplot(soil.data, aes(x = burn_status, y = mn.ppm, fill = burn_status)) +
  geom_boxplot() +
  scale_x_discrete(name = "Burn status") +
  scale_y_continuous(name = "Mn (ppm)") +
  facet_grid(. ~ range) +
  theme_bw() +
  scale_fill_brewer(palette = "Accent") +
  labs(fill = "Burn status") +
  theme(legend.position="none")

# Cu----not normal
normtest(soil.data, 'cu.ppm')
cu <- ggplot(soil.data, aes(x = burn_status, y = cu.ppm, fill = burn_status)) +
  geom_boxplot() +
  scale_x_discrete(name = "Burn status") +
  scale_y_continuous(name = "Cu (ppm)") +
  facet_grid(. ~ range) +
  theme_bw() +
  scale_fill_brewer(palette = "Accent") +
  labs(fill = "Burn status") +
  theme(legend.position="none")

# Ni----not normal
normtest(soil.data, 'ni.ppm')
ni <- ggplot(soil.data, aes(x = burn_status, y = ni.ppm, fill = burn_status)) +
  geom_boxplot() +
  scale_x_discrete(name = "Burn status") +
  scale_y_continuous(name = "Ni (ppm)") +
  facet_grid(. ~ range) +
  theme_bw() +
  scale_fill_brewer(palette = "Accent") +
  labs(fill = "Burn status") +
  theme(legend.position="none")

# Nitrate----not normal
normtest(soil.data, 'no3.n.ppm')
no3 <- ggplot(soil.data, aes(x = burn_status, y = no3.n.ppm, fill = burn_status)) +
  geom_boxplot() +
  scale_x_discrete(name = "Burn status") +
  scale_y_continuous(name = "Nitrate (ppm)") +
  facet_grid(. ~ range) +
  theme_bw() +
  scale_fill_brewer(palette = "Accent") +
  labs(fill = "Burn status") +
  theme(legend.position="none")

# Phosphate----
normtest(soil.data, 'po4.p.ppm')
po4 <- ggplot(soil.data, aes(x = burn_status, y = po4.p.ppm, fill = burn_status)) +
  geom_boxplot() +
  scale_x_discrete(name = "Burn status") +
  scale_y_continuous(name = "Phosphate (ppm)") +
  facet_grid(. ~ range) +
  theme_bw() +
  scale_fill_brewer(palette = "Accent") +
  labs(fill = "Burn status") +
  theme(legend.position="none")

# Sulfate----not normal
normtest(soil.data, 'so4.s.ppm')
so4 <- ggplot(soil.data, aes(x = burn_status, y = so4.s.ppm, fill = burn_status)) +
  geom_boxplot() +
  scale_x_discrete(name = "Burn status") +
  scale_y_continuous(name = "Sulfate (ppm)") +
  facet_grid(. ~ range) +
  theme_bw() +
  scale_fill_brewer(palette = "Accent") +
  labs(fill = "Burn status") +
  theme(legend.position="none")

# B----
normtest(soil.data, 'b.ppm')
b <- ggplot(soil.data, aes(x = burn_status, y = b.ppm, fill = burn_status)) +
  geom_boxplot() +
  scale_x_discrete(name = "Burn status") +
  scale_y_continuous(name = "B (ppm)") +
  facet_grid(. ~ range) +
  theme_bw() +
  scale_fill_brewer(palette = "Accent") +
  labs(fill = "Burn status") +
  theme(legend.position="none")

# esp----not normal
normtest(soil.data, 'esp')
esp <- ggplot(soil.data, aes(x = burn_status, y = esp, fill = burn_status)) +
  geom_boxplot() +
  scale_x_discrete(name = "Burn status") +
  scale_y_continuous(name = "ESP") +
  facet_grid(. ~ range) +
  theme_bw() +
  scale_fill_brewer(palette = "Accent") +
  labs(fill = "Burn status") +
  theme(legend.position="none")

# CEC----not normal
normtest(soil.data, 'cec.meq.100g')
cec <- ggplot(soil.data, aes(x = burn_status, y = cec.meq.100g, fill = burn_status)) +
  geom_boxplot() +
  scale_x_discrete(name = "Burn status") +
  scale_y_continuous(name = "CEC (meq/100g)") +
  facet_grid(. ~ range) +
  theme_bw() +
  scale_fill_brewer(palette = "Accent") +
  labs(fill = "Burn status") +
  theme(legend.position="none")

all.plot <- grid.arrange(ca,mg,k,zn,fe,ni,po4,so4,b,cec,
             nrow = 5, ncol = 2)
ggsave('soil.pdf', plot = all.plot, device = 'pdf', path = fig.dir,
       width = 10, height = 5)

#=========================================================================================
# Differences between ranges----
#=========================================================================================
source('~/Documents/PhD/3_EM_Fire_effect/scripts/functions.R')

responses <- c("pH.su","EC.ds.m","ca.ppm","mg.ppm","na.ppm","k.ppm","zn.ppm","fe.ppm",
               "mn.ppm","cu.ppm","ni.ppm","no3.n.ppm","po4.p.ppm","so4.s.ppm","b.ppm",
               "esp","cec.meq.100g")
ttester(soil.data, 'range', responses)
wilcoxtest(soil.data, 'range', responses)

#=========================================================================================
# Differences between burned and unburned soil----
#=========================================================================================
source('~/Documents/PhD/3_EM_Fire_effect/scripts/functions.R')

responses <- c("pH.su","EC.ds.m","ca.ppm","mg.ppm","na.ppm","k.ppm","zn.ppm","fe.ppm",
               "mn.ppm","cu.ppm","ni.ppm","no3.n.ppm","po4.p.ppm","so4.s.ppm","b.ppm",
               "esp","cec.meq.100g")
ttester(soil.data, 'burn_status', responses)
wilcoxtest(soil.data, 'burn_status', responses)

#=========================================================================================
# soil as a function of distance----
#=========================================================================================

soil.sig <- c('pH.su','po4.p.ppm','so4.s.ppm','b.ppm')

soil.dist <- vegdist(soil.data[colnames(soil.data) %in% soil.sig], method = 'euclidean')
geo.dist <- vegdist(soil.data[c('lat','long')], method = 'euclidean')

mantel(soil.dist,geo.dist)
soil.corr <- mantel.correlog(soil.dist, geo.dist)
plot(soil.corr)

#=========================================================================================
# Check for correlation between sig. soil characateristics----
#=========================================================================================

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
