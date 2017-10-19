## Script created by Liz Bowman April 26, 2017 for Emily Burke
## for analyzing abundance of ectomycorrhizal samples in Pinaleno burn study

#=========================================================================================
# Load data
#=========================================================================================

library(ggplot2)

#--load file paths
#--these should match where ever your data is
dat.dir <- '~/Documents/PhD/3_EM_Fire_effect/data/'
res.dir <- '~/Documents/PhD/3_EM_Fire_effect/results/'
fig.dir <- '~/Documents/PhD/3_EM_Fire_effect/figures/'

#-----------------------------------------------------------------------------------------
# Load and clean up abundance data
#-----------------------------------------------------------------------------------------

pin.ab <- read.csv(paste0(dat.dir,'pinaleno_em_data.csv'), as.is = T)

#--create data frame with abundance calculated per tree
ab.tree <- data.frame(tree = unique(pin.ab$tree_number), range = rep('pinaleno',20),
                      abundance = NA)
for(t in ab.tree$tree) {
  ab.tree[ab.tree$tree == t, 'abundance'] <-
    sum(pin.ab[pin.ab$tree_number == t, 'total_col_per_core'])/
    sum(pin.ab[pin.ab$tree_number == t, 'dry_root_weight_per_core_grams'])
  ab.tree[ab.tree$tree == t, 'burn'] <-
    unique(pin.ab[pin.ab$tree_number == t, 'burn_status'])
}

#--remove outlier tree NF13
ab.tree <- ab.tree[!ab.tree$tree == 'NF13',]

#-----------------------------------------------------------------------------------------
# graph based on burn status
#-----------------------------------------------------------------------------------------

ggplot(ab.tree, aes(burn, abundance)) + 
  geom_boxplot() +
  theme_classic(base_size = 14, base_family = '') +
  ylab('EM abundance') +
  xlab('Burn_status') +
  theme(panel.background = element_rect(color = "black"))

t.test(abundance ~ burn, data = ab.tree)
