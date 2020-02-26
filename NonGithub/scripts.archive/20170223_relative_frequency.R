## Script created by Liz Bowman July 10, 2017
## for analyzing community differences between ranges and burn status of sites

#=========================================================================================
# Load data and libraries----
#=========================================================================================

library(ggplot2); library(tidyr); library(dplyr)

#-----------------------------------------------------------------------------------------
# set up paths to directories----
#-----------------------------------------------------------------------------------------
#--path to directory where the climate data downloaded from BIOCLIM site is stored
dat.dir <- "~/Documents/PhD/3_EM_Fire_effect/data/"
fig.dir <- '~/Documents/PhD/3_EM_Fire_effect/figures/'
res.dir <- "~/Documents/PhD/3_EM_Fire_effect/results/"

#-----------------------------------------------------------------------------------------
# Load data and clean up----
#-----------------------------------------------------------------------------------------
# 
# #--Load data
# all.data <- read.csv(paste0(dat.dir,'20170602_OTU_data.csv'), as.is = T)
# 
# #--Remove sequences not assigned to an OTU
# otu.data <- all.data[!is.na(all.data$OTU.97), ]
# 
# #--Remove sequences not assigned to Ponderosa host
# otu.data <- otu.data[otu.data$Host == 'Ponderosa',]
# 
# #--Add an OTU count column (1 sequence = 1 frequency count)
# otu.data['otu.count'] <- 1
# 
# #--Creates matrix, grouping OTUs by tree number
# otu.data %>%
#   select (Sample_name,Burn_status,Range,Site,Tree,OTU.97,otu.count) %>%
#   spread (OTU.97,otu.count) %>%
#   group_by(Tree,Site,Range,Burn_status) %>%
#   select (matches ("OTU.*")) %>%
#   summarize_each (funs (sum (., na.rm = TRUE))) %>%
#   as.data.frame () -> stsp.matrix
# 
# #--make data frame for relative abundance
# rel.abund <- data.frame(otu = rep(unique(otu.data$OTU.95), each = 4),
#            burn_status = rep(c('burned', 'unburned'), 216),
#            range = rep(c('pinaleno','santa.catalina'), 216),
#            rel.ab = NA)
# rel.abund <- rel.abund[-c(215:216),]
# 
# #--add taxonomic information
# for(o in unique(rel.abund$otu)) {
#   rel.abund[rel.abund$o == o, 'class'] <-
#     unique(otu.data[otu.data$OTU.95 == o, 'Taxonomy_class'])
#   rel.abund[rel.abund$o == o, 'genus'] <-
#     unique(otu.data[otu.data$OTU.95 == o, 'Taxonomy_genus'])
# }
# 
# for (o in unique(rel.abund$otu)) {
#   for (s in c('burned', 'unburned')) {
#     rel.abund[rel.abund$otu == o & rel.abund$burn_status == s, 'rel.ab'] <- 
#     nrow(otu.data[otu.data$OTU.95 == o & otu.data$Burn_status == s, ])/
#       nrow(otu.data[otu.data$Burn_status == s, ])
#   }
# }

#--remove unclassified ones
tax.data <- tax.data[!tax.data$Taxonomy_class == '',]

#--table 
tax.table <- as.data.frame(table(tax.data$Site, tax.data$Taxonomy_class))
colnames(tax.table) <- c('site','taxonomy','count')
for (i in tax.table$site) {
  tax.table[tax.table$site == i, 'range'] <- unique(tax.data[tax.data$Site == i, 'Range'])
  tax.table[tax.table$site == i, 'burn_status'] <-
    unique(tax.data[tax.data$Site == i, 'Burn_status'])
}

tax.table

#=========================================================================================
# Plots----
#=========================================================================================

#--by class
ggplot(tax.table, aes(x = taxonomy, y = count), group = burn_status) +
  geom_col(aes(fill = burn_status, color = burn_status), position = 'dodge') +
  facet_grid(. ~ range) +
  xlab('Class') +
  ylab('Relative frequency') +
  theme(legend.title=element_blank()) +
  scale_y_continuous(limits = c(0,50)) +
  theme(axis.text.x = element_text(angle = 90, hjust =1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text = element_text(size = 16))

ggsave(paste0(fig.dir, 'Relative_frequency.tiff'), plot = last_plot(), device = 'tiff')
