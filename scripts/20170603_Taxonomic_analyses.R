## Script created by Liz Bowman June 3, 2017
## for analyzing taxonomic differences between ranges and burn status of sites

#========================================================================================#
# Load data and libraries----
#========================================================================================#

#----------------------------------------------------------------------------------------#
# Load libraries----
#----------------------------------------------------------------------------------------#
library(ggplot2)

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

#--Change missing data to NA
tax.data[tax.data$Taxonomy_class == '', 'Taxonomy_class'] <- NA
tax.data[tax.data$Taxonomy_genus == '', 'Taxonomy_genus'] <- NA

#--Add column for range burn combo
for (r in unique(tax.data$Range)){
  for (b in unique(tax.data$Burn_status)){
    tax.data[tax.data$Range == r & tax.data$Burn_status == b, 'rangeburn'] <- paste(r,b)
  }
}

#<< Results dataframe >> -------------------
tax.results <- data.frame(tests = c('burn.range.class','burn.class','range.class',
                                    'burn.range.genus','burn.genus','range.genus'),
                          chi.stat = NA,
                          df = NA,
                          p.value = NA)

#========================================================================================#
# Difference between ranges: class level----
#========================================================================================#
#----------------------------------------------------------------------------------------#
# Chi square test: Not significant----
#----------------------------------------------------------------------------------------#

#--Make count table of classes byrange
range.tax <- table(tax.data$Range,tax.data$Taxonomy_class)
range.tax

#--Remove rare species (species with less than 3 occurrences)
rare <- c('Saccharomycetes','Archaeorhizomycetes')
range.tax <- range.tax[ , which(!colnames(range.tax) %in% rare)]

#--chi square test
range.class <- chisq.test(range.tax)
tax.results[tax.results$tests == 'range.class', 'chi.stat'] <- range.class$statistic[[1]]
tax.results[tax.results$tests == 'range.class', 'df'] <- range.class$parameter[[1]]
tax.results[tax.results$tests == 'range.class', 'p.value'] <- range.class$p.value[[1]]

#----------------------------------------------------------------------------------------#
# Plot range class level----
#----------------------------------------------------------------------------------------#

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

#========================================================================================#
# Difference between burn status: class level----
#========================================================================================#
#----------------------------------------------------------------------------------------#
# Chi square test: Not significant----
#----------------------------------------------------------------------------------------#

#--Make count table of classes byrange
burn.tax <- table(tax.data$Burn_status,tax.data$Taxonomy_class)
burn.tax

#--Remove rare species (species with less than 3 occurrences)
rare <- c('Saccharomycetes','Archaeorhizomycetes')
burn.tax <- burn.tax[ , which(!colnames(burn.tax) %in% rare)]

#--chi square test
burn.class <- chisq.test(burn.tax)
tax.results[tax.results$tests == 'burn.class', 'chi.stat'] <- burn.class$statistic[[1]]
tax.results[tax.results$tests == 'burn.class', 'df'] <- burn.class$parameter[[1]]
tax.results[tax.results$tests == 'burn.class', 'p.value'] <- burn.class$p.value[[1]]

#----------------------------------------------------------------------------------------#
# Plot: burn status class level----
#----------------------------------------------------------------------------------------#

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

#=========================================================================================
# Difference between burn status and range: class level----
#=========================================================================================

#-----------------------------------------------------------------------------------------
# Chi square test: Significant----
#-----------------------------------------------------------------------------------------

#--Make count table of classes byrange
rangeburn.tax <- table(tax.data$rangeburn,tax.data$Taxonomy_class)
rangeburn.tax

#--Remove rare species  (species with less than 3 occurrences)
rare <- c('Saccharomycetes','Archaeorhizomycetes')
rangeburn.tax <- rangeburn.tax[ , which(!colnames(rangeburn.tax) %in% rare)]

#--chi square test
rangeburn.class <- chisq.test(rangeburn.tax)
tax.results[tax.results$tests == 'burn.range.class', 'chi.stat'] <- rangeburn.class$statistic[[1]]
tax.results[tax.results$tests == 'burn.range.class', 'df'] <- rangeburn.class$parameter[[1]]
tax.results[tax.results$tests == 'burn.range.class', 'p.value'] <- rangeburn.class$p.value[[1]]

#=========================================================================================
# Difference between burn status and range: class level----
#=========================================================================================

#-----------------------------------------------------------------------------------------
# Chi square test: Significant----
#-----------------------------------------------------------------------------------------

#--Make count table of classes byrange
rangeburn.tax <- table(tax.data$rangeburn,tax.data$Taxonomy_class)
rangeburn.tax

#--Remove rare species  (species with less than 3 occurrences)
rare <- c('Saccharomycetes','Archaeorhizomycetes')
rangeburn.tax <- rangeburn.tax[ , which(!colnames(rangeburn.tax) %in% rare)]

#--chi square test
rangeburn.class <- chisq.test(rangeburn.tax)
tax.results[tax.results$tests == 'burn.range.class', 'chi.stat'] <- rangeburn.class$statistic[[1]]
tax.results[tax.results$tests == 'burn.range.class', 'df'] <- rangeburn.class$parameter[[1]]
tax.results[tax.results$tests == 'burn.range.class', 'p.value'] <- rangeburn.class$p.value[[1]]

#-----------------------------------------------------------------------------------------
# Plot: burn status and range at class level----
#-----------------------------------------------------------------------------------------

#--remove rare species from the tax data file
range.tax <- tax.data[which(!tax.data$Taxonomy_class %in% rare), ]
range.tax <- range.tax[!is.na(range.tax$Taxonomy_class),]

#Change levels of Burn_status and Range columns
range.tax$Burn_status <- factor(range.tax$Burn_status, levels = c('burned','unburned'),
       labels = c('Burned', 'Unburned'))
range.tax$Range <- factor(range.tax$Range, levels = c('pinaleno','santa.catalina'),
                                labels = c('Pinaleno Mts.', 'Santa Catalina Mts.'))

#--Bar graph
rangeburn.class <- ggplot(data = range.tax, 
                     aes(x = Burn_status,
                         fill = Taxonomy_class)) + 
  geom_bar(position = "fill") + 
  ylab("Proportion of \n sequences per class") +
  facet_grid(. ~ Range) +
  theme_bw() +
  #ggtitle("Proportion of Classes by Topography") +
  scale_fill_brewer(palette = "Greens") +
  guides (fill=guide_legend(title=NULL)) +
  theme(legend.position="none",
        axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size=22, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

rangeburn.class

ggsave('TaxonomyClass_SequenceBased.tiff', plot = rangeburn.class,
       device = 'tiff', path = fig.dir,
       width = 20, height = 20, units = 'cm')

#=========================================================================================
# Difference between ranges: genus level----
#=========================================================================================
#-----------------------------------------------------------------------------------------
# Chi square test: all phyla----
#-----------------------------------------------------------------------------------------

#--Make count table of classes byrange
range.tax <- table(tax.data$Range,tax.data$Taxonomy_genus)
range.tax

#--Remove rare species  (species with less than 3 occurrences)
rare <- colnames(range.tax[, colSums(range.tax) < 4])
range.tax <- range.tax[ ,colSums(range.tax) > 4]

#--chi square test
range.genus <- chisq.test(range.tax)
tax.results[tax.results$tests == 'range.genus', 'chi.stat'] <- range.genus$statistic[[1]]
tax.results[tax.results$tests == 'range.genus', 'df'] <- range.genus$parameter[[1]]
tax.results[tax.results$tests == 'range.genus', 'p.value'] <- range.genus$p.value[[1]]

#=========================================================================================
# Difference between burn status: genus level----
#=========================================================================================
#-----------------------------------------------------------------------------------------
# Chi square test: all phyla----
#-----------------------------------------------------------------------------------------

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

#=========================================================================================
# Difference between burn and range: genus level----
#=========================================================================================
#-----------------------------------------------------------------------------------------
# Chi square test: all phyla----
#-----------------------------------------------------------------------------------------

#--Make count table of classes byrange
burnrange.tax <- table(tax.data$rangeburn,tax.data$Taxonomy_genus)
burnrange.tax

#--Remove rare species  (species with less than 3 occurrences)
rare <- colnames(burnrange.tax[, colSums(burnrange.tax) < 4])
burnrange.tax <- burnrange.tax[ ,colSums(burnrange.tax) > 4]

#--chi square test
burnrange.genus <- chisq.test(range.tax)
tax.results[tax.results$tests == 'burn.range.genus', 'chi.stat'] <- 
  burnrange.genus$statistic[[1]]
tax.results[tax.results$tests == 'burn.range.genus', 'df'] <- 
  burnrange.genus$parameter[[1]]
tax.results[tax.results$tests == 'burn.range.genus', 'p.value'] <- 
  burnrange.genus$p.value[[1]]

write.csv(tax.results, paste0(res.dir,'TaxResults_SequenceBased.csv'),
          row.names = F)

#-----------------------------------------------------------------------------------------
# Plot: Range and burn status at genus level----
#-----------------------------------------------------------------------------------------

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
  theme(legend.position="none",
        axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size=22, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

#<< Basid genera >> ----------------------------------------------------------------------
#--isolate Basidiomycete genera
basid <- c('Agaricomycetes')
basid.tax <- tax.data[tax.data$Taxonomy_class %in% basid, ]

#--remove rare species from the tax data file (removed those with 2 occurences or less)
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
  theme(legend.position="none",
        axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size=22, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))
