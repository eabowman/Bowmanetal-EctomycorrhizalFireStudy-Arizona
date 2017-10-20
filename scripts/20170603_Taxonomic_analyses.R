## Script created by Liz Bowman June 3, 2017
## for analyzing taxonomic differences between ranges and burn status of sites

#=========================================================================================
# Load data and libraries----
#=========================================================================================

#-----------------------------------------------------------------------------------------
# Load libraries----
#-----------------------------------------------------------------------------------------
library(ggplot2)

#-----------------------------------------------------------------------------------------
# set up paths to directories----
#-----------------------------------------------------------------------------------------
#--path to directory where the climate data downloaded from BIOCLIM site is stored
dat.dir <- "~/Documents/PhD/2_EM_Fire_effect/data/"
fig.dir <- '~/Documents/PhD/2_EM_Fire_effect/figures/'
res.dir <- "~/Documents/PhD/2_EM_Fire_effect/results/"

source('~/Documents/PhD/2_EM_Fire_effect/scripts/functions.R')

#-----------------------------------------------------------------------------------------
# Load data----
#-----------------------------------------------------------------------------------------
tax.data <- read.csv(paste0(dat.dir,'20170806_OTU_data.csv'), as.is = T)

tax.data <- tax.data[tax.data$Host == 'Ponderosa',]

#--Change missing data to NA
tax.data[tax.data$Taxonomy_class == '', 'Taxonomy_class'] <- NA
tax.data[tax.data$Taxonomy_genus == '', 'Taxonomy_genus'] <- NA

#=========================================================================================
# Difference between ranges: class level----
#=========================================================================================
#-----------------------------------------------------------------------------------------
# Chi square test----
#-----------------------------------------------------------------------------------------

#--Make count table of classes byrange
range.tax <- table(tax.data$Range,tax.data$Taxonomy_class)
range.tax

#--Remove rare species
rare <- c('Saccharomycetes','Archaeorhizomycetes')
range.tax <- range.tax[ , which(!colnames(range.tax) %in% rare)]

#--chi square test
chisq.test(range.tax)

#-----------------------------------------------------------------------------------------
# Plot range class level----
#-----------------------------------------------------------------------------------------

# #--remove rare species from the tax data file
# range.tax <- tax.data[which(!tax.data$Taxonomy_class %in% rare), ]
# range.tax <- range.tax[!is.na(range.tax$Taxonomy_class),]
# #--Bar graph
# range.class <- ggplot(data = range.tax, 
#              aes(x = Range,
#                  fill = Taxonomy_class)) + 
#   geom_bar(position = "fill") + 
#   ylab("Proportion of sequences per class") +
#   xlab("Range") +
#   #ggtitle("Proportion of Classes by Topography") +
#   scale_fill_brewer(palette = "Greys") +
#   guides (fill=guide_legend(title=NULL)) +
#   theme_classic(base_size = 12)

#=========================================================================================
# Difference between burn status: class level----
#=========================================================================================
#-----------------------------------------------------------------------------------------
# Chi square test----
#-----------------------------------------------------------------------------------------

#--Make count table of classes byrange
burn.tax <- table(tax.data$Burn_status,tax.data$Taxonomy_class)
burn.tax

#--Remove rare species
rare <- c('Saccharomycetes','Archaeorhizomycetes')
burn.tax <- burn.tax[ , which(!colnames(burn.tax) %in% rare)]

#--chi square test
chisq.test(burn.tax)

#-----------------------------------------------------------------------------------------
# Plot: burn status class level----
#-----------------------------------------------------------------------------------------

# #--remove rare species from the tax data file
# range.tax <- tax.data[which(!tax.data$Taxonomy_class %in% rare), ]
# range.tax <- range.tax[!is.na(range.tax$Taxonomy_class),]
# #--Bar graph
# burn.class <- ggplot(data = range.tax, 
#                       aes(x = Burn_status,
#                           fill = Taxonomy_class)) + 
#   geom_bar(position = "fill") + 
#   ylab("Proportion of sequences per class") +
#   xlab("Burn status") +
#   #ggtitle("Proportion of Classes by Topography") +
#   scale_fill_brewer(palette = "Greys") +
#   guides (fill=guide_legend(title=NULL)) +
#   theme_classic(base_size = 12)

#-----------------------------------------------------------------------------------------
# Plot: burn status and range at class level----
#-----------------------------------------------------------------------------------------

#--remove rare species from the tax data file
range.tax <- tax.data[which(!tax.data$Taxonomy_class %in% rare), ]
range.tax <- range.tax[!is.na(range.tax$Taxonomy_class),]
#--Bar graph
burn.class <- ggplot(data = range.tax, 
                     aes(x = Burn_status,
                         fill = Taxonomy_class)) + 
  geom_bar(position = "fill") + 
  ylab("Proportion of sequences per class") +
  facet_grid(. ~ Range) +
  #ggtitle("Proportion of Classes by Topography") +
  scale_fill_brewer(palette = "Rainbow") +
  guides (fill=guide_legend(title=NULL)) +
  theme_classic(base_size = 12)

#=========================================================================================
# Difference between ranges: genus level----
#=========================================================================================
#-----------------------------------------------------------------------------------------
# Chi square test: all phyla----
#-----------------------------------------------------------------------------------------

#--Make count table of classes byrange
range.tax <- table(tax.data$Range,tax.data$Taxonomy_genus)
range.tax

#--Remove rare species
rare <- c('Ramaria','Amphinema','Boletus','Coltricia','Elaphomyces','Hysterangium',
          'Oidiodendron','Piloderma','Tulasnella','Archaeorhizomyces')
range.tax <- range.tax[ ,which(!colnames(range.tax) %in% rare)]

#--chi square test
chisq.test(range.tax)

#-----------------------------------------------------------------------------------------
# Plot: Range and burn status at genus level----
#-----------------------------------------------------------------------------------------

#<< Ascomyceteous genera >> --------------------------------------------------------------
#--isolate Ascomycete genera
asco <- c('Pezizomycetes', 'Leotiomycetes','Dothideomycetes','Eurotiomycetes')
asco.tax <- tax.data[tax.data$Taxonomy_class %in% asco, ]
asco.tax <- asco.tax[!is.na(asco.tax$Taxonomy_genus), ]

#--remove rare species from the tax data file
#tax.data <- tax.data[which(!tax.data$Taxonomy_class %in% rare), ]
#--Bar graph: Ascomycota
asco <- ggplot(data = asco.tax, 
                      aes(x = Burn_status,
                          fill = Taxonomy_genus)) + 
  geom_bar(position = "fill") + 
  ylab("Proportion of sequences per genus") +
  xlab(NULL) +
  facet_grid(. ~ Range) +
  guides (fill=guide_legend(title=NULL)) +
  theme_classic(base_size = 12) +
  theme(legend.position = 'bottom')

#<< Basid genera >> ----------------------------------------------------------------------
#--isolate Basidiomycete genera
basid <- c('Agaricomycetes')
basid.tax <- tax.data[tax.data$Taxonomy_class %in% basid, ]

#--remove rare species from the tax data file (removed those with 2 occurences or less)
basid.rare <- c('Amphinema','Boletus','Coltricia','Hysterangium','Tulasnella')
basid.tax <- basid.tax[which(!basid.tax$Taxonomy_genus %in% basid.rare), ]
#--Bar graph: Ascomycota
basid <- ggplot(data = basid.tax, 
                     aes(x = Burn_status,
                         fill = Taxonomy_genus)) + 
  geom_bar(position = "fill") + 
  ylab(NULL) +
  xlab(NULL) +
  scale_color_gradientn(colours = rainbow(16)) +
  guides (fill=guide_legend(title=NULL)) +
  theme_classic(base_size = 12)  +
  facet_grid(. ~ Range) +
  theme(legend.position = 'bottom')

#=========================================================================================
# Difference between burn status: genus level----
#=========================================================================================
#-----------------------------------------------------------------------------------------
# Chi square test----
#-----------------------------------------------------------------------------------------

#--Make count table of classes byrange
burn.tax <- table(tax.data$Burn_status,tax.data$Taxonomy_genus)
burn.tax

#--chi square test
chisq.test(burn.tax)

tax.plot <- multiplot(asco, basid, cols=2)

#ggsave('Taxonomic_plots.png', plot = tax.plot, device = 'png', path = fig.dir)

