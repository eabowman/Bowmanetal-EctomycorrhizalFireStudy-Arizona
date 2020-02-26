## Script created by Liz Bowman June 4, 2017
## model fitting of diversity between ranges and burn status of sites

#========================================================================================#
# Load data and libraries----
#========================================================================================#

#----------------------------------------------------------------------------------------#
# Load libraries----
#----------------------------------------------------------------------------------------#
library(ggplot2);library(tidyr);library(vegan);library(dplyr);
library(car); library(rcompanion); library(nlme)

#----------------------------------------------------------------------------------------#
# set up paths to directories----
#----------------------------------------------------------------------------------------#
#--path to directory where the climate data downloaded from BIOCLIM site is stored
dat.dir <- "~/Documents/PhD/2_EM_Fire_effect/data/"
fig.dir <- '~/Documents/PhD/2_EM_Fire_effect/figures_output/'
res.dir <- "~/Documents/PhD/2_EM_Fire_effect/results_output/"


source('~/Documents/PhD/2_EM_Fire_effect/scripts/functions.R')

#----------------------------------------------------------------------------------------#
# Load data and clean up----
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
  select (Sample_name,Burn_status,Range,Site,Tree,otu.97,otu.count) %>%
  spread (otu.97,otu.count) %>%
  group_by(Tree,Site,Range,Burn_status) %>%
  select (matches ("otu.*")) %>%
  summarize_each (funs (sum (., na.rm = TRUE))) %>%
  as.data.frame () -> stsp.matrix

#--make a count table
burn.div <- table(otu.data$Burn_status, otu.data$otu.97)

#--make table for burned data
burned.div <- as.data.frame(burn.div['burned',])
colnames(burned.div) <- 'count'


#--table for unburned data
unburned.div <- as.data.frame(burn.div['unburned',])
colnames(unburned.div) <- 'count'

#--<< Combine range and burn status >> ----------
for(r in unique(stsp.matrix$Range)) {
  for(b in unique(stsp.matrix$Burn_status)) {
    stsp.matrix[stsp.matrix$Burn_status == b & stsp.matrix$Range == r, 'burn_range'] <-
      paste(b, r)
  }
}

#--add forest type
stsp.matrix$forest <- NA
stsp.matrix[stsp.matrix$Site %in% c('sc1','sc2','sc3','sc4','p2'), 'forest'] <- 'pine oak'
stsp.matrix[stsp.matrix$Site %in% c('p3','p4'), 'forest'] <- 'pine df'
stsp.matrix[stsp.matrix$Site == 'p1', 'forest'] <- 'pine'

#--reorder matrix
colnames(stsp.matrix)
stsp.matrix <- stsp.matrix[c(1:4,121,122,5:120)]

#<< w/o singletons: by tree >> ------------------------
nonsingletons <- as.data.frame(which(colSums(stsp.matrix[7:length(stsp.matrix)]) > 1))
nonsingletons <- rownames(nonsingletons)
nonsingletons <- stsp.matrix[nonsingletons]
nonsingletons$Tree <- stsp.matrix$Tree
stsp.wo.singletons <- stsp.matrix[1:6]
stsp.wo.singletons <- left_join(stsp.wo.singletons, nonsingletons, by = NULL)

#<< Creates matrix, grouping OTUs by site, burn status, and range >> -
otu.data %>%
  select(Sample_name,Range,Burn_status,Site,Tree,otu.97,otu.count) %>%
  spread(otu.97,otu.count) %>%
  group_by(Site, Range, Burn_status) %>%
  select(matches("otu*")) %>%
  summarize_each(funs(sum(., na.rm = TRUE))) %>%
  as.data.frame() -> rangeburn.matrix

#--add burn_range column
for(r in unique(rangeburn.matrix$Range)) {
  for(b in unique(rangeburn.matrix$Burn_status)) {
    rangeburn.matrix[rangeburn.matrix$Burn_status == b
                     & rangeburn.matrix$Range == r, 'burn_range'] <-
      paste(b, r)
  }
}

#--add forest type
rangeburn.matrix$forest <- NA
rangeburn.matrix[rangeburn.matrix$site %in% c('sc1','sc2','sc3','sc4','p2'),
                 'forest'] <- 'pine oak'
rangeburn.matrix[rangeburn.matrix$site %in% c('p3','p4'), 'forest'] <- 'pine df'
rangeburn.matrix[rangeburn.matrix$site == 'p1', 'forest'] <- 'pine'

colnames(rangeburn.matrix)
rangeburn.matrix <- rangeburn.matrix[c(1:3, 120, 121, 4:119)]

#--<< w/o singletons: by site >> -
nonsingletons <- as.data.frame(which(colSums(rangeburn.matrix[6:length(rangeburn.matrix)]) > 1))
nonsingletons <- rownames(nonsingletons)
nonsingletons <- rangeburn.matrix[nonsingletons]
nonsingletons$Site <- rangeburn.matrix$Site
rangeburn.wo.singletons <- rangeburn.matrix[1:5]
rangeburn.wo.singletons <- 
  left_join(rangeburn.wo.singletons, nonsingletons, by = NULL)

#--<<table for SR and diversity by range>>----
div.table <- stsp.matrix[7:length(stsp.matrix)]
div.data <- data.frame(tree = stsp.matrix$Tree,
                       site = stsp.matrix$Site,
                       range = stsp.matrix$Range,
                       fire.history = stsp.matrix$Burn_status,
                       forest.type = stsp.matrix$forest,
                       spec.richness = specnumber(div.table),
                       fishers.alpha = fisher.alpha(div.table),
                       shannon = diversity(div.table, index = 'shannon'),
                       Pielous.evenness =  
                       diversity(div.table, index = 'shannon')/log(specnumber(div.table)))

#----------------------------------------------------------------------------------------#
# Create dataframe for results----
#----------------------------------------------------------------------------------------#
DivResults.SequenceBased <- data.frame(div.measure = c(rep('fishers.alpha',2),
                                                       rep('shannon',2)),
                                       range = rep(c('scm','pm'), 2),
                                       df = NA,
                                       t.stat = NA,
                                       p.value = NA)

SpecRichness.SequenceBased <- data.frame(data = c('w/ singletons', 'w/o singletons'),
                                   avg = NA,
                                   sd = NA)
SpecRichnessResults.SequenceBased <- data.frame(data = c('scm', 'pm'),
                                       df = NA,
                                       t.stat = NA,
                                       p.value = NA)

#========================================================================================#
# Species richness----
#========================================================================================#

#--<< Table of species abundance: w/ singletons >> ------------------------
#--Overall
spec <- table(colSums(stsp.matrix[7:length(stsp.matrix)]))
plot(spec, xlab = 'Individuals per species',
     ylab = 'Number of species')

#--<< Stats for species richness >> ------------
#mean
mean.sr.by.tree <- mean(div.data$spec.richness)
SpecRichness.SequenceBased[SpecRichness.SequenceBased$data == 'w/ singletons',
                           'avg'] <- mean.sr.by.tree

#standard deviation
sd.sr.by.tree <- sd(div.data$spec.richness)
SpecRichness.SequenceBased[SpecRichness.SequenceBased$data == 'w/ singletons',
                           'sd'] <- sd.sr.by.tree

#--<< By range and burn >> --------
sr.by.burnrange <- specnumber(stsp.matrix[7:length(stsp.matrix)],
                              groups = stsp.matrix$burn_range)
sr.by.burnrange

#--<< Assess within range differences of species richness based on fire history >>-----
#--T-test Santa Catalina
scm.sr <- div.data[div.data$range == 'santa.catalina',]
scm.sr.t <- t.test(scm.sr$spec.richness ~ scm.sr$fire.history)
SpecRichnessResults.SequenceBased[SpecRichnessResults.SequenceBased$data == 'scm',
                                  'df'] <- scm.sr.t$parameter[1]
SpecRichnessResults.SequenceBased[SpecRichnessResults.SequenceBased$data == 'scm',
                                  't.stat'] <- scm.sr.t$statistic[1]
SpecRichnessResults.SequenceBased[SpecRichnessResults.SequenceBased$data == 'scm',
                                  'p.value'] <- scm.sr.t$p.value[1]
#--T-test Pinaleno
pm.sr <- div.data[div.data$range == 'pinaleno',]
pm.sr.t <- t.test(pm.sr$spec.richness ~ pm.sr$fire.history)
SpecRichnessResults.SequenceBased[SpecRichnessResults.SequenceBased$data == 'pm',
                                  'df'] <- pm.sr.t$parameter[1]
SpecRichnessResults.SequenceBased[SpecRichnessResults.SequenceBased$data == 'pm',
                                  't.stat'] <- pm.sr.t$statistic[1]
SpecRichnessResults.SequenceBased[SpecRichnessResults.SequenceBased$data == 'pm',
                                  'p.value'] <- pm.sr.t$p.value[1]

write.csv(SpecRichnessResults.SequenceBased, 'results_output/species_richness_results.csv',
          row.names = F)

#----------------------------------------------------------------------------------------#
# Plot of Species richness by range and burn_status: tree----
#----------------------------------------------------------------------------------------#
stsp.matrix$Range <- as.factor(stsp.matrix$Range)
levels(stsp.matrix$Range) <- c('Pinaleno Mts.', 'Santa Catalina Mts.')
stsp.matrix$Burn_status <- as.factor(stsp.matrix$Burn_status)
levels(stsp.matrix$Burn_status) <- c('FA', 'FU')
stsp.matrix$spec.richness <- sr.by.tree

sr.plot <- ggplot(stsp.matrix, aes(x = Burn_status, y = spec.richness,
                                                  fill = Burn_status)) +
  geom_boxplot() +
  facet_grid(~ Range) +
  scale_x_discrete(name = "Fire history") +
  scale_y_continuous(name = "Species richness") +
  theme_bw() +
  scale_fill_brewer(palette = "Accent") +
  labs(fill = "Fire history") +
  scale_fill_manual(values=c('#ef8a62','#999999')) +
  theme(legend.position="none",
        axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size=22, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

sr.plot
ggsave('SpecRichness_FireHistory.jpeg', plot = sr.plot,
       device = 'jpeg', path = 'figures_output/',
       width = 20, height = 15, units = 'cm')


#========================================================================================#
# Diversity: Fisher's alpha -------
#========================================================================================#
#--remove outliers (those that have equal species richness and abundance)
fa.out <- c('F11','F12','F13', 'F14','F4','LB021','NF19')
fa.data <- div.data[!div.data$tree %in% fa.out, ]

#--Average FA by fire
mean.FA.burn <- data.frame(fire.history = unique(fa.data$fire.history),
                      mean = NA,
                      sd = NA)
for(i in mean.FA.burn$fire.history){
  mean.i <- mean(fa.data[fa.data$fire.history == i, 'fishers.alpha'])
  sd.i <- sd(fa.data[fa.data$fire.history == i, 'fishers.alpha'])
  mean.FA.burn[mean.FA.burn$fire.history == i, 'mean'] <- mean.i
  mean.FA.burn[mean.FA.burn$fire.history == i, 'sd'] <- sd.i
}

#--Average FA by fire and range
mean.FA.rb <- data.frame(fire.history = rep(unique(fa.data$fire.history),2),
                           range = rep(unique(fa.data$range),each = 2),
                           mean = NA,
                           sd = NA)
for(i in unique(mean.FA.rb$fire.history)){
  for(r in unique(mean.FA.rb$range)){
    mean.i <- mean(fa.data[fa.data$fire.history == i & fa.data$range == r,
                           'fishers.alpha'])
    sd.i <- sd(fa.data[fa.data$fire.history == i & fa.data$range == r,
                       'fishers.alpha'])
    mean.FA.rb[mean.FA.rb$fire.history == i & mean.FA.rb$range == r, 'mean'] <- mean.i
    mean.FA.rb[mean.FA.rb$fire.history == i & mean.FA.rb$range == r, 'sd'] <- sd.i
  }
}

#--Average FA by site
mean.FA.site <- data.frame(site = c('sc1','sc2','sc3','sc4','p1','p2','p3','p4'),
                         fire.history = rep(c('unburned','unburned','burned','burned'),2),
                         range = c(rep('santa.catalina',4),rep('pinaleno',4)),
                         mean = NA,
                         sd = NA)
for(i in unique(mean.FA.site$site)){
  mean.i <- mean(fa.data[fa.data$site == i,'fishers.alpha'])
  sd.i <- sd(fa.data[fa.data$site,'fishers.alpha'])
  mean.FA.site[mean.FA.site$site == i, 'mean'] <- mean.i
  mean.FA.site[mean.FA.site$site == i, 'sd'] <- sd.i
}


#--<< Assess within range differences of fisher's alpha based on fire history >>-----
#--T-test Santa Catalina
scm.fa <- fa.data[fa.data$range == 'santa.catalina',]
scm.fa.t <- t.test(scm.fa$fishers.alpha ~ scm.fa$fire.history)
DivResults.SequenceBased[DivResults.SequenceBased$range == 'scm' & 
                           DivResults.SequenceBased$div.measure == 'fishers.alpha',
                                  'df'] <- scm.fa.t$parameter[1]
DivResults.SequenceBased[DivResults.SequenceBased$range == 'scm' & 
                           DivResults.SequenceBased$div.measure == 'fishers.alpha',
                                  't.stat'] <- scm.fa.t$statistic[1]
DivResults.SequenceBased[DivResults.SequenceBased$range == 'scm' & 
                           DivResults.SequenceBased$div.measure == 'fishers.alpha',
                                  'p.value'] <- scm.fa.t$p.value[1]
#--T-test Pinaleno
pm.fa <- fa.data[fa.data$range == 'pinaleno',]
pm.fa.t <- t.test(pm.fa$fishers.alpha ~ pm.fa$fire.history)
DivResults.SequenceBased[DivResults.SequenceBased$range == 'pm' & 
                           DivResults.SequenceBased$div.measure == 'fishers.alpha',
                         'df'] <- pm.fa.t$parameter[1]
DivResults.SequenceBased[DivResults.SequenceBased$range == 'pm' & 
                           DivResults.SequenceBased$div.measure == 'fishers.alpha',
                         't.stat'] <- pm.fa.t$statistic[1]
DivResults.SequenceBased[DivResults.SequenceBased$range == 'pm' & 
                           DivResults.SequenceBased$div.measure == 'fishers.alpha',
                         'p.value'] <- pm.fa.t$p.value[1]

#----------------------------------------------------------------------------------------#
# Plot of Fisher's alpha by burn status: tree data with singletons----
#----------------------------------------------------------------------------------------#
levels(fa.data$range) <- c('Pinaleno Mts.', 'Santa Catalina Mts.')
levels(fa.data$fire.history) <- c('FA', 'FU')

fa.plot <- ggplot(fa.data, aes(x = fire.history,
                               y = fishers.alpha,
                               fill = fire.history)) +
  geom_boxplot() +
  scale_x_discrete(name = "Fire history") +
  facet_grid(~ range) +
  scale_y_continuous(name = "Fisher's alpha") +
  theme_bw() +
  labs(fill = "Fire history") +
  scale_fill_manual(values=c('#ef8a62','#999999')) +
  theme(legend.position="none",
        axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size=22, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

fa.plot

ggsave('FishersAlpha_FireHistory.jpeg', plot = fa.plot,
       device = 'jpeg', path = 'figures_output/',
       width = 20, height = 15, units = 'cm')

#========================================================================================#
# Diversity: Shannon's diversity -------
#========================================================================================#
#--<< Assess within range differences of Shannon's diversity based on fire history >>-----
#--T-test Santa Catalina
scm.sd <- div.data[div.data$range == 'santa.catalina',]
scm.sd.t <- t.test(scm.sd$shannon ~ scm.sd$fire.history)
DivResults.SequenceBased[DivResults.SequenceBased$range == 'scm' & 
                           DivResults.SequenceBased$div.measure == 'shannon',
                         'df'] <- scm.sd.t$parameter[1]
DivResults.SequenceBased[DivResults.SequenceBased$range == 'scm' & 
                           DivResults.SequenceBased$div.measure == 'shannon',
                         't.stat'] <- scm.sd.t$statistic[1]
DivResults.SequenceBased[DivResults.SequenceBased$range == 'scm' & 
                           DivResults.SequenceBased$div.measure == 'shannon',
                         'p.value'] <- scm.sd.t$p.value[1]
#--T-test Pinaleno
pm.sd <- fa.data[div.data$range == 'pinaleno',]
pm.sd.t <- t.test(pm.sd$spec.richness ~ pm.sd$fire.history)
DivResults.SequenceBased[DivResults.SequenceBased$range == 'pm' & 
                           DivResults.SequenceBased$div.measure == 'shannon',
                         'df'] <- pm.sd.t$parameter[1]
DivResults.SequenceBased[DivResults.SequenceBased$range == 'pm' & 
                           DivResults.SequenceBased$div.measure == 'shannon',
                         't.stat'] <- pm.sd.t$statistic[1]
DivResults.SequenceBased[DivResults.SequenceBased$range == 'pm' & 
                           DivResults.SequenceBased$div.measure == 'shannon',
                         'p.value'] <- pm.sd.t$p.value[1]

write.csv(DivResults.SequenceBased, 'results_output/DiversityResults.csv',
          row.names = F)

#----------------------------------------------------------------------------------------#
# Plot of Shannon's diversity by range and burn_status: tree----
#----------------------------------------------------------------------------------------#
levels(fa.data$range) <- c('Pinaleno Mts.', 'Santa Catalina Mts.')
levels(fa.data$fire.history) <- c('FA', 'FU')

shannon.plot <- ggplot(fa.data, aes(x = fire.history, y = shannon,
                                  fill = fire.history)) +
  geom_boxplot() +
  facet_grid(~ range) +
  scale_x_discrete(name = "Fire history") +
  scale_y_continuous(name = "Shannon diversity index") +
  theme_bw() +
  scale_fill_brewer(palette = "Accent") +
  labs(fill = "Fire history") +
  scale_fill_manual(values=c('#ef8a62','#999999')) +
  theme(legend.position="none",
        axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size=22, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

shannon.plot
ggsave('Shannon_FireHistory.jpeg', plot = shannon.plot,
       device = 'jpeg', path = 'figures_output/',
       width = 20, height = 15, units = 'cm')

