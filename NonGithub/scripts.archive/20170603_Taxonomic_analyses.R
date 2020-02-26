## Script created by Liz Bowman June 3, 2017
## for analyzing taxonomic differences between ranges and burn status of sites

#========================================================================================#
# Load data and libraries----
#========================================================================================#

#----------------------------------------------------------------------------------------#
# Load libraries----
#----------------------------------------------------------------------------------------#
library(ggplot2); library(plyr)

#----------------------------------------------------------------------------------------#
# set up paths to directories----
#----------------------------------------------------------------------------------------#
dat.dir <- "~/Documents/PhD/2_EM_Fire_effect/data/"
fig.dir <- '~/Documents/PhD/2_EM_Fire_effect/figures_output/'
res.dir <- "~/Documents/PhD/2_EM_Fire_effect/results_output/"

source('~/Documents/PhD/2_EM_Fire_effect/scripts/functions.R')

#----------------------------------------------------------------------------------------#
# Load and clean up data----
#----------------------------------------------------------------------------------------#
tax.data <- read.csv(paste0(dat.dir,'20170806_OTU_data.csv'), as.is = T)
tax.data <- tax.data[tax.data$Host == 'Ponderosa',]

#--Add column for range burn combo
for (r in unique(tax.data$Range)){
  for (b in unique(tax.data$Burn_status)){
    tax.data[tax.data$Range == r & tax.data$Burn_status == b, 'rangeburn'] <- paste(r,b)
  }
}

#<< Results dataframe >> -------------------
tax.results <- data.frame(tests = c('burn.range.class','burn.class','scm.class',
                                    'pm.class','burn.range.genus','burn.genus',
                                    'scm.genus', 'pm.genus'),
                          chi.stat = NA,
                          df = NA,
                          p.value = NA)

tax.data %>%
  select(Sample_name, Burn_status, Range, Taxonomy_class, Taxonomy_genus, Tip_count) %>%
  spread(Taxonomy_class, Tip_count, fill = 0) %>%
  group_by(Range) %>%
  summarise(sum = sum(.))-> test

#========================================================================================#
# Class level----
#========================================================================================#

#----------------------------------------------------------------------------------------#
# Chi square test: Significant----
#----------------------------------------------------------------------------------------#
rare <- c('Saccharomycetes','Archaeorhizomycetes','Eurotiomycetes')
#<< Range and Burn >> -----------------------------
# #--Make count table of classes byrange
# rangeburn.tax <- table(tax.data$rangeburn,tax.data$Taxonomy_class)
# rangeburn.tax
# 
# #--Remove rare species  (species with less than 3 occurrences)
# rangeburn.tax <- rangeburn.tax[ , which(!colnames(rangeburn.tax) %in% rare)]
# 
# #--chi square test
# rangeburn.chi <- chisq.test(rangeburn.tax, correct = F)
# tax.results[tax.results$tests == 'burn.range.class', 'chi.stat'] <- rangeburn.chi$statistic[[1]]
# tax.results[tax.results$tests == 'burn.range.class', 'df'] <- rangeburn.chi$parameter[[1]]
# tax.results[tax.results$tests == 'burn.range.class', 'p.value'] <- rangeburn.chi$p.value[[1]]

#<< Range >> -----------------------------
#--Overall between the ranges
range.tab.tax <- table(tax.data$Range, tax.data$Taxonomy_class)

# remove rare classes, Saccharomycetes, Archaeorhizomycetes, Eurotiomycetes
range.tab.tax <- range.tab.tax[,which(!colnames(range.tab.tax) %in% rare)]

range.chi <- chisq.test(range.tab.tax, correct = F)

#--Santa Catalina Mts.
scm.tax <- tax.data[tax.data$Range == "santa.catalina",]

scm.tab.tax <- table(scm.tax$Burn_status, scm.tax$Taxonomy_class)
scm.tab.tax <- scm.tab.tax[, which(!colnames(scm.tab.tax) %in% rare)]

scm.chi <- chisq.test(scm.tab.tax, correct = F)
scm.chi
tax.results[tax.results$tests == 'scm.class', 'chi.stat'] <- range.chi$statistic[[1]]
tax.results[tax.results$tests == 'scm.class', 'df'] <- range.chi$parameter[[1]]
tax.results[tax.results$tests == 'scm.class', 'p.value'] <- range.chi$p.value[[1]]

#--Pinaleno Mts.
pm.tax <- tax.data[tax.data$Range == "pinaleno",]

pm.tab.tax <- table(pm.tax$Burn_status, pm.tax$Taxonomy_class)
pm.tab.tax <- pm.tab.tax[, which(!colnames(pm.tab.tax) %in% rare)]

range.chi <- chisq.test(pm.tab.tax, correct = F)
range.chi
tax.results[tax.results$tests == 'pm.class', 'chi.stat'] <- range.chi$statistic[[1]]
tax.results[tax.results$tests == 'pm.class', 'df'] <- range.chi$parameter[[1]]
tax.results[tax.results$tests == 'pm.class', 'p.value'] <- range.chi$p.value[[1]]

#<< Burn >> -----------------------------
burn.tax <- table(tax.data$Burn_status, tax.data$Taxonomy_class)
burn.tax <- burn.tax[, which(!colnames(burn.tax) %in% rare)]

burn.chi <- chisq.test(burn.tax, correct = F)
burn.chi
tax.results[tax.results$tests == 'burn.class', 'chi.stat'] <- burn.chi$statistic[[1]]
tax.results[tax.results$tests == 'burn.class', 'df'] <- burn.chi$parameter[[1]]
tax.results[tax.results$tests == 'burn.class', 'p.value'] <- burn.chi$p.value[[1]]

#<< Fire history >> -----------------------------
#--burned data
burn.tax <- tax.data[tax.data$Burn_status == 'burned',]
# remove rare classes
burn.tax <- burn.tax[!burn.tax$Taxonomy_class %in% c('Saccharomycetes',
                                                     'Eurotiomycetes'),]
burn.table <- table(burn.tax$Range, burn.tax$Taxonomy_class)
burn.chi <- chisq.test(burn.table, correct = F)
burn.chi

#--unburned data
unburn.tax <- tax.data[tax.data$Burn_status == 'unburned',]
# remove rare classes
unburn.tax <- unburn.tax[!unburn.tax$Taxonomy_class %in% c('Archaeorhizomycetes',
                                                     'Eurotiomycetes'),]
unburn.table <- table(unburn.tax$Range, unburn.tax$Taxonomy_class)
unburn.chi <- chisq.test(unburn.table, correct = F)
unburn.chi

#----------------------------------------------------------------------------------------#
# Plot: burn status and range at class level----
#----------------------------------------------------------------------------------------#

#--remove rare species from the tax data file
# range.tax <- tax.data[which(!tax.data$Taxonomy_class %in% rare), ]
# range.tax <- range.tax[!is.na(range.tax$Taxonomy_class),]
# 
# #Change levels of Burn_status and Range columns
# range.tax$Burn_status <- as.factor(range.tax$Burn_status)
# levels(range.tax$Burn_status) <- c('FA', 'FU')
# range.tax$Range <- as.factor(range.tax$Range)
# levels(range.tax$Range) <- c('Pinaleno Mts.', 'Santa Catalina Mts.')
# 
# #--Bar graph
# rangeburn.class <- ggplot(data = range.tax, 
#                      aes(x = Burn_status,
#                          fill = Taxonomy_class)) + 
#   geom_bar(position = "fill") + 
#   ylab("Proportion of \n sequences per class") +
#   facet_grid(. ~ Range) +
#   theme_bw() +
#   xlab('Fire history') +
#   #ggtitle("Proportion of Classes by Topography") +
#   scale_fill_brewer(palette = "Greens") +
#   guides (fill=guide_legend(title=NULL)) +
#   theme(legend.position='right',
#         axis.title.x = element_text(margin = margin(t = 30)),
#         axis.title.y = element_text(margin = margin(r = 30)),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.text = element_text(size=22, color = 'black'),
#         axis.title = element_text(size = 28),
#         strip.text.x = element_text(size = 14))
# 
# rangeburn.class
# 
# ggsave('TaxonomyClass_SequenceBased.tiff', plot = rangeburn.class,
#        device = 'jpeg', path = fig.dir,
#        width = 20, height = 15, units = 'cm')

#----------------------------------------------------------------------------------------#
# Plot: Insets for figure 3----
#----------------------------------------------------------------------------------------#

#--Isolated unburned data
unburn.tax <- tax.data[tax.data$Burn_status == 'unburned',]
unburn.tax <- unburn.tax[!unburn.tax$Taxonomy_class %in% c('Archaeorhizomycetes',
                                                           'Eurotiomycetes',NA),]
#Change levels of Burn_status and Range columns
unburn.tax[unburn.tax$Range == 'santa.catalina', 'Range'] <- 'Santa Catalina Mts.'
unburn.tax[unburn.tax$Range == 'pinaleno', 'Range'] <- 'Pinaleno Mts.'
unburn.tax$Range <-  factor(unburn.tax$Range, 
       levels = c('Santa Catalina Mts.', 'Pinaleno Mts.'))

#--Bar graph
unburn.class <- ggplot(data = unburn.tax, 
                          aes(x = Range,
                              fill = Taxonomy_class)) + 
  geom_bar(position = "fill") + 
  ylab("Proportion of \n sequences per class") +
  theme_bw() +
  xlab(element_blank()) +
  #ggtitle("Proportion of Classes by Topography") +
  scale_fill_brewer(palette = 2,
                    direction = -1) +
  guides (fill=guide_legend(title=NULL)) +
  theme(legend.position='right',
        # axis.title.x = element_text(margin = margin(t = 30)),
        # axis.title.y = element_text(margin = margin(r = 30)),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        axis.text = element_text(size=20, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

unburn.class

ggsave('TaxonomyClass_SequenceBased_Unburned.tiff', plot = unburn.class,
       device = 'jpeg', path = fig.dir,
       width = 20, height = 15, units = 'cm')

#--Isolated burned data
burn.tax <- tax.data[tax.data$Burn_status == 'burned',]
burn.tax <- burn.tax[!burn.tax$Taxonomy_class %in% c('Saccharomycetes',
                                                     'Eurotiomycetes',NA),]
#Change levels of Burn_status and Range columns
burn.tax[burn.tax$Range == 'santa.catalina', 'Range'] <- 'Santa Catalina Mts.'
burn.tax[burn.tax$Range == 'pinaleno', 'Range'] <- 'Pinaleno Mts.'
burn.tax$Range <-  factor(burn.tax$Range, 
                            levels = c('Santa Catalina Mts.', 'Pinaleno Mts.'))

#--Bar graph
burn.class <- ggplot(data = burn.tax, 
                       aes(x = Range,
                           fill = Taxonomy_class)) + 
  geom_bar(position = "fill") + 
  ylab("Proportion of \n sequences per class") +
  theme_bw() +
  xlab(element_blank()) +
  #ggtitle("Proportion of Classes by Topography") +
  scale_fill_brewer(palette = 6,
                    direction = -1) +
  guides (fill=guide_legend(title=NULL)) +
  theme(legend.position='right',
        # axis.title.x = element_text(margin = margin(t = 30)),
        # axis.title.y = element_text(margin = margin(r = 30)),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        axis.text = element_text(size=20, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

burn.class

ggsave('TaxonomyClass_SequenceBased_Burned.tiff', plot = burn.class,
       device = 'jpeg', path = fig.dir,
       width = 20, height = 15, units = 'cm')

#========================================================================================#
# Difference between ranges: genus level----
#========================================================================================#
#----------------------------------------------------------------------------------------#
# Chi square test: all phyla----
#----------------------------------------------------------------------------------------#
#<< range >> -------------
#--Santa Catalina Mts.
#--Make count table of classes byrange
scm.tab.tax <- table(scm.tax$Burn_status,scm.tax$Taxonomy_genus)
scm.tab.tax

#--Remove rare species  (species with less than 3 occurrences)
rare <- colnames(scm.tab.tax[, colSums(scm.tab.tax) < 4])
scm.tab.tax <- scm.tab.tax[ ,colSums(scm.tab.tax) > 4]

#--chi square test
scm.genus <- chisq.test(scm.tab.tax)
scm.genus
tax.results[tax.results$tests == 'scm.genus', 'chi.stat'] <- range.genus$statistic[[1]]
tax.results[tax.results$tests == 'scm.genus', 'df'] <- range.genus$parameter[[1]]
tax.results[tax.results$tests == 'scm.genus', 'p.value'] <- range.genus$p.value[[1]]

#--Pinaleno Mts.
#--Make count table of classes byrange
pm.tab.tax <- table(pm.tax$Burn_status,pm.tax$Taxonomy_genus)
pm.tab.tax

#--Remove rare species  (species with less than 3 occurrences)
rare <- colnames(pm.tab.tax[, colSums(pm.tab.tax) < 4])
pm.tab.tax <- pm.tab.tax[ ,colSums(pm.tab.tax) > 4]

#--chi square test
pm.genus <- chisq.test(pm.tab.tax)
pm.genus
tax.results[tax.results$tests == 'pm.genus', 'chi.stat'] <- range.genus$statistic[[1]]
tax.results[tax.results$tests == 'pm.genus', 'df'] <- range.genus$parameter[[1]]
tax.results[tax.results$tests == 'pm.genus', 'p.value'] <- range.genus$p.value[[1]]

#<< Fire history >> -------------
#--Make count table of classes byrange
burn.tax <- table(tax.data$Burn_status,tax.data$Taxonomy_genus)
burn.tax

#--Remove rare species  (species with less than 3 occurrences)
rare <- colnames(burn.tax[, colSums(burn.tax) < 4])
burn.tax <- burn.tax[ ,colSums(burn.tax) > 4]

#--chi square test
burn.genus <- chisq.test(burn.tax)
tax.results[tax.results$tests == 'burn.genus', 'chi.stat'] <- burn.genus$statistic[[1]]
tax.results[tax.results$tests == 'burn.genus', 'df'] <- burn.genus$parameter[[1]]
tax.results[tax.results$tests == 'burn.genus', 'p.value'] <- burn.genus$p.value[[1]]

#<< Fire history and Range >> -------------
#--Make count table of classes byrange
burnrange.tax <- table(tax.data$rangeburn,tax.data$Taxonomy_genus)
burnrange.tax

#--Remove rare species  (species with less than 3 occurrences)
rare <- colnames(burnrange.tax[, colSums(burnrange.tax) < 4])
burnrange.tax <- burnrange.tax[ ,colSums(burnrange.tax) > 4]

#--chi square test
burnrange.genus <- chisq.test(range.tax)
burnrange.genus
tax.results[tax.results$tests == 'burn.range.genus', 'chi.stat'] <- 
  burnrange.genus$statistic[[1]]
tax.results[tax.results$tests == 'burn.range.genus', 'df'] <- 
  burnrange.genus$parameter[[1]]
tax.results[tax.results$tests == 'burn.range.genus', 'p.value'] <- 
  burnrange.genus$p.value[[1]]

write.csv(tax.results, paste0(res.dir,'TaxResults_SequenceBased.csv'),
          row.names = F)

#----------------------------------------------------------------------------------------#
# Plot: Range and burn status at genus level----
#----------------------------------------------------------------------------------------#

#<< Ascomyceteous genera >> --------------------------------------------------------------
#--isolate Ascomycete genera
asco <- c('Pezizomycetes', 'Leotiomycetes','Dothideomycetes','Eurotiomycetes')
asco.tax <- tax.data[tax.data$Taxonomy_class %in% asco, ]
asco.tax <- asco.tax[!is.na(asco.tax$Taxonomy_genus), ]

#--Remove rare
asco.tax <- asco.tax[!asco.tax$Taxonomy_genus %in% rare,]

#--Bar graph: Ascomycota
asco <- ggplot(data = asco.tax, 
                      aes(x = Burn_status,
                          fill = Taxonomy_genus)) + 
  geom_bar(position = "fill") + 
  ylab("Proportion of sequences per genus") +
  xlab(NULL) +
  theme_bw() +
  facet_grid(. ~ Range) +
  guides (fill=guide_legend(title=NULL)) +
  theme(legend.position="bottom",
        axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size=22, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

asco

#<< Basid genera >> ----------------------------------------------------------------------
#--isolate Basidiomycete genera
basid <- c('Agaricomycetes')
basid.tax <- tax.data[tax.data$Taxonomy_class %in% basid, ]

#--remove rare species from the tax data file (removed those with 2 occurences or less)
rare <- colnames(burnrange.tax[, colSums(burnrange.tax) < 5])
basid.tax <- basid.tax[!basid.tax$Taxonomy_genus %in% rare,]
basid.tax <- basid.tax[!is.na(basid.tax$Taxonomy_genus),]
#--Bar graph: Ascomycota
basid <- ggplot(data = basid.tax, 
                     aes(x = Burn_status,
                         fill = Taxonomy_genus)) + 
  geom_bar(position = "fill") + 
  ylab(NULL) +
  xlab(NULL) +
  theme_bw() +
  scale_color_gradientn(colours = rainbow(16)) +
  guides (fill=guide_legend(title=NULL)) +
  theme_classic(base_size = 12)  +
  facet_grid(. ~ Range) +
  theme(legend.position="bottom",
        axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size=22, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

basid

#----------------------------------------------------------------------------------------#
# Heat map: Genus level---
#----------------------------------------------------------------------------------------#


#<< By range and fire >> ----
rh.heat <- read.csv('data_output/Taxonomic_GenusRangeFire_HeatMap.csv', as.is = T, header = T)

#--rescale abundance values
for(i in 1:44){
  rh.heat[i, 'ab.rel'] <- rh.heat[i, 'Ab.rf']/sum(rh.heat[1:44,'Ab.rf'])
}

#--reorder phyla so Ascos and Basids group together
rh.heat$Phyla <- as.character(rh.heat$Phyla)
rh.heat$Phyla <- factor(rh.heat$Phyla, levels=unique(rh.heat$Phyla))
rh.heat <- arrange(rh.heat, Phyla)

rh.h <- ggplot(data = rh.heat, aes(x=Fire.history, y=Genus, fill=ab.rel)) + 
  geom_tile(aes(alpha = ab.rel, fill = Phyla), colour = "darkgrey") +
  facet_grid(~ Range) +
  theme_classic(base_size = 12) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size=18, color = 'black'),
        axis.title = element_text(size = 26),
        strip.text.x = element_text(size = 12),
        axis.text.x = element_text(vjust = 0),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

ggsave('HeatMap_RangeFireHistoryGenus.jpeg', plot = rh.h, 
       device = 'jpeg', path = 'figures_output/',
       width = 20, height = 15, units = 'cm')


#<< Range differences: overall >> ----
range.heat <- read.csv('data_output/Taxonomic_GenusRange_HeatMap.csv', as.is = T, header = T)

#--rescale abundance values
for(i in 1:22){
  range.heat[i, 'ab.rel'] <- range.heat[i, 'Ab.range']/sum(range.heat[1:22,'Ab.range'])
}

#--reorder phyla so Ascos and Basids group together
range.heat$Phyla <- as.character(range.heat$Phyla)
range.heat$Phyla <- factor(range.heat$Phyla, levels=unique(range.heat$Phyla))
range.heat <- arrange(range.heat, Phyla)

r.h <- ggplot(data = range.heat, aes(x=Range, y=Genus, fill=ab.rel)) + 
  geom_tile(aes(alpha = ab.rel, fill = Phyla), colour = "darkgrey") +
  theme_classic(base_size = 12) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size=18, color = 'black'),
        axis.title = element_text(size = 26),
        strip.text.x = element_text(size = 12),
        axis.text.x = element_text(vjust = 0),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

ggsave('HeatMap_RangeGenus.jpeg', plot = r.h, 
       device = 'jpeg', path = 'figures_output/',
       width = 20, height = 15, units = 'cm')

#<< Fire activity differences: overall >> ----
fire.heat <- read.csv('data_output/Taxonomic_GenusFireHistory_HeatMap.csv', as.is = T, header = T)

#--rescale abundance values
for(i in 1:22){
  fire.heat[i, 'ab.rel'] <- fire.heat[i, 'Ab.fire']/sum(fire.heat[1:22,'Ab.fire'])
}

#--reorder phyla so Ascos and Basids group together
fire.heat$Phyla <- as.character(fire.heat$Phyla)
fire.heat$Phyla <- factor(fire.heat$Phyla, levels=unique(fire.heat$Phyla))
fire.heat <- arrange(fire.heat, Phyla)

f.h <- ggplot(data = fire.heat, aes(x=Fire.history, y=Genus, fill=ab.rel)) + 
  geom_tile(aes(alpha = ab.rel, fill = Phyla), colour = "darkgrey") +
  theme_classic(base_size = 12) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size=18, color = 'black'),
        axis.title = element_text(size = 26),
        strip.text.x = element_text(size = 12),
        axis.text.x = element_text(vjust = 0),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

ggsave('HeatMap_FireHistoryGenus.jpeg', plot = f.h, 
       device = 'jpeg', path = 'figures_output/',
       width = 20, height = 15, units = 'cm')
