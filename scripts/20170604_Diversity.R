## Script created by Liz Bowman June 4, 2017
## model fitting of diversity between ranges and burn status of sites

#=========================================================================================
# Load data and libraries----
#=========================================================================================

#-----------------------------------------------------------------------------------------
# Load libraries----
#-----------------------------------------------------------------------------------------
library(ggplot2);library (tidyr);library (vegan);library (dplyr)

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

#--reorder matrix
stsp.matrix <- stsp.matrix[c(1:4,121,5:120)]

#--Creates matrix, grouping OTUs by site and burn status
otu.data %>%
  select(Sample_name,Burn_status,Range,Site,Tree,otu.97,otu.count) %>%
  spread(otu.97,otu.count) %>%
  group_by(Site, Burn_status) %>%
  select(matches("otu*")) %>%
  summarize_each(funs(sum(., na.rm = TRUE))) %>%
  as.data.frame() -> burn.matrix

#--Creates matrix, grouping OTUs by site and range
otu.data %>%
  select(Sample_name,Burn_status,Range,Site,Tree,otu.97,otu.count) %>%
  spread(otu.97,otu.count) %>%
  group_by(Site, Range) %>%
  select(matches("otu*")) %>%
  summarize_each(funs(sum(., na.rm = TRUE))) %>%
  as.data.frame() -> range.matrix

#--Creates matrix, grouping OTUs by site, burn status, and range
otu.data %>%
  select(Sample_name,Range,Burn_status,Site,Tree,otu.97,otu.count) %>%
  spread(otu.97,otu.count) %>%
  group_by(Site, Range, Burn_status) %>%
  select(matches("otu*")) %>%
  summarize_each(funs(sum(., na.rm = TRUE))) %>%
  as.data.frame() -> rangeburn.matrix
# add burn_range column
for(r in unique(rangeburn.matrix$Range)) {
  for(b in unique(rangeburn.matrix$Burn_status)) {
    rangeburn.matrix[rangeburn.matrix$Burn_status == b
                     & rangeburn.matrix$Range == r, 'burn_range'] <-
      paste(b, r)
  }
}
rangeburn.matrix <- rangeburn.matrix[c(1:3, 120, 4:119)]

#=========================================================================================
# Species richness----
#=========================================================================================

# << Table of species abundance >> ------------------------
#--Overall
spec <- table(colSums(stsp.matrix[6:length(stsp.matrix)]))
plot(spec, xlab = 'Individuals per species',
     ylab = 'Number of species')

# #--burned
# spec.burned <- table(colSums(stsp.matrix[stsp.matrix$Burn_status == 'burned' , 
#                                          6:length(stsp.matrix)]))
# spec.burned <- spec.burned[-1]
# plot(spec.burned, xlab = 'Individuals per species',
#      ylab = 'Number of species', main = 'Burned sites')
# 
# #--unburned
# spec.unburned <- table(colSums(stsp.matrix[stsp.matrix$Burn_status == 'unburned' , 
#                                            6:length(stsp.matrix)]))
# spec.unburned <- spec.unburned[-1]
# plot(spec.unburned, xlab = 'Individuals per species',
#      ylab = 'Number of species', main = 'Unburned sites')
# 
# #--Pinaleno
# spec.pin <- table(colSums(stsp.matrix[stsp.matrix$Range == 'pinaleno' , 
#                                          6:length(stsp.matrix)]))
# spec.pin <- spec.pin[-1]
# plot(spec.pin, xlab = 'Individuals per species',
#      ylab = 'Number of species', main = 'Pinaleno Mts.')
# 
# #--Santa Catalina
# spec.scm <- table(colSums(stsp.matrix[stsp.matrix$Range == 'santa.catalina' , 
#                                            6:length(stsp.matrix)]))
# spec.scm <- spec.scm[-1]
# plot(spec.scm, xlab = 'Individuals per species',
#      ylab = 'Number of species', main = 'Santa Catalina Mts.')


# << Burn status >> --------
specnumber(stsp.matrix[6:length(stsp.matrix)], groups = stsp.matrix$Burn_status)

# << Range >> --------
range.s <- specnumber(stsp.matrix[6:length(stsp.matrix)], groups = stsp.matrix$Range)
t.test(range.s) # Significantly different

# << Range and Burn status >> --------
specnumber(stsp.matrix[6:length(stsp.matrix)], groups = stsp.matrix$burn_range)

#<< By tree >> ------------
specnumber(stsp.matrix[6:length(stsp.matrix)], groups = stsp.matrix$Tree)

#=========================================================================================
# Diversity: Fisher's alpha----not normally distributed
#=========================================================================================

div.data <- data.frame(tree = stsp.matrix$Tree,
           range = stsp.matrix$Range,
           burn_status = stsp.matrix$Burn_status,
           burn_range = stsp.matrix$burn_range,
           fa = fisher.alpha(stsp.matrix[6:length(stsp.matrix)]))
#--remove outliers (those that have equal species richness and abundance)
fa.out <- c('F11','F12','F14','F4','LB021','NF19')
div.data <- div.data[!div.data$tree %in% fa.out, ]

#<< Burn status Fisher's alpha >> --------------
#--grouped by site
burn.fa <- fisher.alpha(burn.matrix[3:length(burn.matrix)]) # not sig. diff.
wilcox.test(burn.fa ~ burn.matrix$Burn_status)

#--grouped by tree
wilcox.test(div.data$fa ~ div.data$burn_status) # no sig. diff.
burned.fa <- mean(div.data[div.data$burn_status == 'burned', 'fa'])
unburned.fa <- mean(div.data[div.data$burn_status == 'unburned', 'fa'])

#<< Range Fisher's alpha >> --------------
#--grouped by site
range.fa <- fisher.alpha(range.matrix[3:length(range.matrix)]) # Pinaleno vs SCM order
wilcox.test(range.fa ~ range.matrix$Range) # not sig. different 

#--grouped by tree
wilcox.test(div.data$fa ~ div.data$range) # not sign. diff. w/ burned having higher FA
pm.fa <- mean(div.data[div.data$range == 'pinaleno', 'fa'])
scm.fa <- mean(div.data[div.data$range == 'santa.catalina', 'fa'])

#<< Range and burn status Fisher's alpha >> --------------
#--grouped by site
rangeburn.fa <- fisher.alpha(rangeburn.matrix[5:length(rangeburn.matrix)])
# Pinaleno burned, Pinaleno unburned, SCM burned, SCM unburned order
anova(lm(rangeburn.fa ~ rangeburn.matrix$burn_range))
# not sign. different with Pinaleno being overall higher
# Pinaleno burned higher than unburned; SCM unburned higher than burned

#--grouped by tree
anova(lm(div.data$fa ~ div.data$burn_range)) # not sign. diff. w/ burned having higher FA
pm.burned.fa <- mean(div.data[div.data$burn_range == 'burned pinaleno', 'fa'])
pm.unburned.fa <- mean(div.data[div.data$burn_range == 'unburned pinaleno', 'fa'])
scm.burned.fa <- mean(div.data[div.data$burn_range == 'burned santa.catalina', 'fa'])
scm.unburned.fa <- mean(div.data[div.data$burn_range == 'unburned santa.catalina', 'fa'])

#-----------------------------------------------------------------------------------------
# Plot of Fisher's alpha by range and burn_status----
#-----------------------------------------------------------------------------------------
# By tree
#--remove tree F13
div.data.out <- div.data[!div.data$tree == 'F13',]
 
ggplot(burnrange.fa, aes(x = burn_status, y = fisher.alpha, fill = burn_status)) +
  geom_boxplot() +
  scale_x_discrete(name = "Burn status") +
  scale_y_continuous(name = "Fisher's alpha") +
  facet_grid(. ~ range) +
  theme_bw() +
  scale_fill_brewer(palette = "Accent") +
  labs(fill = "Burn status") +
  theme(legend.position="none")

# By site
burnrange.div <- data.frame(range = rangeburn.matrix$Range,
           burn_status = rangeburn.matrix$Burn_status,
           fisher.alpha = rangeburn.fa)

ggplot(burnrange.div, aes(x = burn_status, y = fisher.alpha, fill = burn_status)) +
  geom_boxplot() +
  scale_x_discrete(name = "Burn status") +
  scale_y_continuous(name = "Fisher's alpha") +
  facet_grid(. ~ range) +
  theme_bw() +
  scale_fill_brewer(palette = "Accent") +
  labs(fill = "Burn status") +
  theme(legend.position="none")

#=========================================================================================
# Diversity: Shannon's diversity---- normally distributed
#=========================================================================================

div.data$shannon <- diversity(stsp.matrix[6:length(stsp.matrix)], index = 'shannon')

#<< Burn status Fisher's alpha >> --------------
burn.shannon <- diversity(burn.matrix[3:length(burn.matrix)], index = 'shannon')
# burned vs. unburned order
t.test(burn.shannon ~ burn.matrix$Burn_status) 
# not sign. different with burned having higher FA

#<< Range Fisher's alpha >> --------------
range.shannon <- diversity(range.matrix[3:length(range.matrix)], index = 'shannon')
# Pinaleno vs SCM order
t.test(range.shannon ~ range.matrix$Range) 
# not sign. different 

#<< Range Fisher's alpha >> --------------
rangeburn.shannon <- diversity(rangeburn.matrix[5:length(rangeburn.matrix)],
                               index = 'shannon')
# Pinaleno burned, Pinaleno unburned, SCM burned, SCM unburned order
anova(lm(rangeburn.shannon ~ rangeburn.matrix$burn_range))
# not sign. different with Pinaleno being overall higher
# Pinaleno burned higher than unburned; SCM unburned higher than burned

#-----------------------------------------------------------------------------------------
# Plot of Shannon's diversity by range and burn_status----
#-----------------------------------------------------------------------------------------

# By tree
ggplot(div.data, aes(x = burn_status, y = shannon, fill = burn_status)) +
  geom_boxplot() +
  scale_x_discrete(name = "Burn status") +
  scale_y_continuous(name = "Shannon's diversity") +
  facet_grid(. ~ range) +
  theme_bw() +
  scale_fill_brewer(palette = "Accent") +
  labs(fill = "Burn status") +
  theme(legend.position="none")

# By site
burnrange.div$shannon <- rangeburn.shannon

ggplot(burnrange.div, aes(x = burn_status, y = shannon, fill = burn_status)) +
  geom_boxplot() +
  scale_x_discrete(name = "Burn status") +
  scale_y_continuous(name = "Shannon's diversity") +
  facet_grid(. ~ range) +
  theme_bw() +
  scale_fill_brewer(palette = "Accent") +
  labs(fill = "Burn status") +
  theme(legend.position="none")

#=========================================================================================
# Fit data to diversity models----
#=========================================================================================

k <-sample(nrow(stsp.matrix), 1)
fish <- fisherfit(stsp.matrix[k,6:length(stsp.matrix)])
fish

prestondistr(stsp.matrix[k,6:length(stsp.matrix)])
radfit(stsp.matrix[2,6:length(stsp.matrix)])

ggplot(data = unburned.div, aes(count)) +
  geom_histogram(binwidth = 2)
