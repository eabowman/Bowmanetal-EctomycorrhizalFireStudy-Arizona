## Script created by Liz Bowman May 9, 2019
## model fitting of diversity between ranges and burn status of sites

#========================================================================================#
# Load data and libraries----
#========================================================================================#

#----------------------------------------------------------------------------------------#
# Load libraries----
#----------------------------------------------------------------------------------------#
library(ggplot2);library(vegan);library(car); library(rcompanion); library(nlme)

#install.packages("tidyverse")
library(tidyverse)

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

#----------------------------------------------------------------------------------------#
# Create dataframe for results----
#----------------------------------------------------------------------------------------#
DivResults.SequenceBased <- data.frame(div.measure = c('specrich.scm', 'specrich.pm',
                                       'fisher.scm','fisher.pm',
                                       'shannon.scm', 'shannon.pm'),
                                       t.stat = NA,
                                       df = NA,
                                       p.value = NA)

SpecRichness.SequenceBased <- data.frame(data = c('w/ singletons', 'w/o singletons'),
                                   avg = NA,
                                   sd = NA)
SpecRichnessResults.SequenceBased <- data.frame(data = c('w/ singletons.range',
                                       'w/ singletons.burn',
                                       'w/o singletons.range',
                                       'w/o singletons.burn'),
                                       df1 = NA,
                                       df2 = NA,
                                       f.stat = NA,
                                       p.value = NA)


#========================================================================================#
# Species richness----
#========================================================================================#
#--Make species richness table
sr.by.tree <- specnumber(stsp.matrix[7:length(stsp.matrix)], 
                         groups = stsp.matrix$Tree)
sp.rich <- stsp.matrix[1:6]
sp.rich['spec.richness'] <- sr.by.tree

#--Create tables based on Pinalneo and Santa Catalina Mts.
scm.sp.rich <- sp.rich[sp.rich$Range == 'santa.catalina',]
pm.sp.rich <- sp.rich[sp.rich$Range == 'pinaleno',]

# << Table of species abundance: w/ singletons >> ------------------------
#--Overall
spec <- table(colSums(stsp.matrix[7:length(stsp.matrix)]))
plot(spec, xlab = 'Individuals per species',
     ylab = 'Number of species')

#<< Mean species richness: By tree >> ------------
#mean
mean.sr.by.tree <- mean(sr.by.tree)
SpecRichness.SequenceBased[SpecRichness.SequenceBased$data == 'w/ singletons',
                           'avg'] <- mean.sr.by.tree

#standard deviation
sd.sr.by.tree <- sd(sr.by.tree)
SpecRichness.SequenceBased[SpecRichness.SequenceBased$data == 'w/ singletons',
                           'sd'] <- sd.sr.by.tree

# << Mean species richness: By range and burn >> --------
sr.by.burnrange <- specnumber(stsp.matrix[7:length(stsp.matrix)],
                              groups = stsp.matrix$burn_range)
sr.by.burnrange


#<< Mean species richness: By tree >> ------------
sr.by.wo.tree <- specnumber(stsp.wo.singletons[7:length(stsp.wo.singletons)], 
                         groups = stsp.wo.singletons$Tree)

mean.sr.by.wo.tree <- mean(sr.by.wo.tree)
SpecRichness.SequenceBased[SpecRichness.SequenceBased$data == 'w/o singletons',
                           'avg'] <- mean.sr.by.wo.tree

sd.sr.by.wo.tree <- sd(sr.by.wo.tree)
SpecRichness.SequenceBased[SpecRichness.SequenceBased$data == 'w/o singletons',
                           'sd'] <- sd.sr.by.wo.tree


#<<Analysis of differences in species richness between FA and FU sites>>-----

#--Santa Catalina Mts.
lme.burn.scm <- lme(spec.richness ~ Burn_status, data = scm.sp.rich,
                random = ~ 1 | Site, method = 'ML')
summary.burn.scm <- summary(lme.burn.scm)
anova.burn.scm <- anova(lme.burn.scm)

sp.scm <- t.test(spec.richness ~ Burn_status, data = scm.sp.rich)

#--Pinaleno Mts.
lme.burn.pm <- lme(spec.richness ~ Burn_status, data = pm.sp.rich,
                    random = ~ 1 | Site, method = 'ML')
summary.burn.pm <- summary(lme.burn.pm)
anova.burn.pm <- anova(lme.burn.pm)

sp.pm <- t.test(spec.richness ~ Burn_status, data = pm.sp.rich)

#--Add results to DivResults.SequencBased dataframe
DivResults.SequenceBased[DivResults.SequenceBased$div.measure == 'specrich.scm',
                         't.stat'] <- sp.scm$statistic
DivResults.SequenceBased[DivResults.SequenceBased$div.measure == 'specrich.scm',
                         'df'] <- sp.scm$parameter
DivResults.SequenceBased[DivResults.SequenceBased$div.measure == 'specrich.scm',
                         'p.value'] <- sp.scm$p.value
DivResults.SequenceBased[DivResults.SequenceBased$div.measure == 'specrich.pm',
                         't.stat'] <- sp.pm$statistic
DivResults.SequenceBased[DivResults.SequenceBased$div.measure == 'specrich.pm',
                         'df'] <- sp.pm$parameter
DivResults.SequenceBased[DivResults.SequenceBased$div.measure == 'specrich.pm',
                         'p.value'] <- sp.pm$p.value

#----------------------------------------------------------------------------------------#
# Plot of Species richness by range and burn_status: tree----
#----------------------------------------------------------------------------------------#
#--Reorder factors SCM, PM; Unburned, Burned
stsp.matrix[stsp.matrix$Range == 'santa.catalina', 'Range'] <- 'Santa Catalina Mts.'
stsp.matrix[stsp.matrix$Range == 'pinaleno', 'Range'] <- 'Pinaleno Mts.'
stsp.matrix$Range <- factor(stsp.matrix$Range, 
                            levels = c('Santa Catalina Mts.', 'Pinaleno Mts.'))

stsp.matrix[stsp.matrix$Burn_status == 'burned', 'Burn_status'] <- 'Burned'
stsp.matrix[stsp.matrix$Burn_status == 'unburned', 'Burn_status'] <- 'Unburned'
stsp.matrix$Burn_status <- factor(stsp.matrix$Burn_status,
                                  levels = c('Unburned','Burned'))

stsp.matrix$spec.richness <- sr.by.tree
stsp.matrix <- stsp.matrix[c(1:6, 124, 7:123)]

# stsp.matrix <- as.tibble(stsp.matrix)
# stsp.matrix %>%
#   ungroup() %>%
#   # 2. Arrange by
#   #   i.  Range
#   #   ii. Fire history
#   arrange(desc(Range),desc(Burn_status)) %>%
#   # 3. Add order column of row numbers
#   mutate(order = row_number())
# 
# stsp.matrix %>%
#   distinct(Range) %>%
#   mutate(Range = fct_relevel(Range, 
#                              c('Santa Catalina Mts', 'Pinaleno Mts.'))) %>%
#   arrange(Range)

sr.plot <- ggplot(stsp.matrix, aes(x = Burn_status, y = spec.richness,
                                                  fill = Burn_status)) +
  geom_boxplot() +
  facet_grid(~ Range) +
  scale_x_discrete(name = "Fire history") +
  scale_y_continuous(name = "Species richness") +
  theme_bw() +
  scale_fill_brewer(palette = "Accent") +
  labs(fill = "Fire history") +
  scale_fill_manual(values=c('#006d2c','#636363')) +
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
# Diversity -------
#========================================================================================#

#<< By Tree: create dataframes >> ------
div.data <- data.frame(tree = stsp.matrix$Tree,
                       site = stsp.matrix$Site,
                       range = stsp.matrix$Range,
                       burn_status = stsp.matrix$Burn_status,
                       burn_range = stsp.matrix$burn_range,
                       forest = stsp.matrix$forest,
                       fa = fisher.alpha(stsp.matrix[8:length(stsp.matrix)]),
                       shannon = diversity(stsp.matrix[8:length(stsp.matrix)],
                                           index = 'shannon'))

#--Check for outliers (i.e. where species richness = species abundance)
out.tree <- div.data[div.data$fa > 50, 'tree']
out.info <- data.frame(tree = out.tree,
             spec.richness = sp.rich[sp.rich$Tree %in% out.tree, 'spec.richness'],
             spec.abundance = rowSums(stsp.matrix[stsp.matrix$Tree %in% out.tree,
                                                8:length(stsp.matrix)]))

div.data.out <- div.data[!div.data$tree %in% out.tree, ]

# << Check for normality: with singletons >> -----------------
#--Fisher's alpha
#Distribution
ggplot(div.data.out, aes(x = fa, colour = range)) +
  geom_density()

#residuals
plot(aov(fa ~ site, data = div.data.out),1)

#qqplot
plot(aov(fa ~ range, data = div.data.out),2)
shapiro.test(log(div.data.out$fa))

#Homogeneity of variance
leveneTest(div.data.out$fa, div.data.out$range,
           center = median)

#--Shannon's diversity index
#Distribution
ggplot(div.data, aes(x = shannon, colour = range)) +
  geom_density()

plotNormalHistogram(div.data$shannon)

#residuals
plot(aov(shannon ~ range, data = div.data),1)

#qqplot
plot(aov(shannon ~ range, data = div.data),2)
shapiro.test(div.data$shannon)

#Homogeneity of variance
leveneTest(div.data$shannon, div.data$range,
           center = median)

# << Diversity analysis >> ------------------

#--Fisher's alpha diversity----
scm.div <- div.data.out[div.data.out$range == 'Santa Catalina Mts.',]
pm.div <- div.data.out[div.data.out$range == 'Pinaleno Mts.',]

#--Santa Catalina Mts.
scm.fa <- t.test(scm.div$fa ~ scm.div$burn_status)

#--Pinaleno Mts.
pm.fa <- t.test(pm.div$fa ~ pm.div$burn_status)

#--Add results to DivResults.SequencBased dataframe
DivResults.SequenceBased[DivResults.SequenceBased$div.measure == 'fisher.scm',
                         't.stat'] <- scm.fa$statistic
DivResults.SequenceBased[DivResults.SequenceBased$div.measure == 'fisher.scm',
                         'df'] <- scm.fa$parameter
DivResults.SequenceBased[DivResults.SequenceBased$div.measure == 'fisher.scm',
                         'p.value'] <- scm.fa$p.value
DivResults.SequenceBased[DivResults.SequenceBased$div.measure == 'fisher.pm',
                         't.stat'] <- pm.fa$statistic
DivResults.SequenceBased[DivResults.SequenceBased$div.measure == 'fisher.pm',
                         'df'] <- pm.fa$parameter
DivResults.SequenceBased[DivResults.SequenceBased$div.measure == 'fisher.pm',
                         'p.value'] <- pm.fa$p.value


#--Shannon's diversity----
scm.shannon <- div.data[div.data$range == 'Santa Catalina Mts.',]
pm.shannon <- div.data[div.data$range == 'Pinaleno Mts.',]

#--Santa Catalina Mts.
scm.sd <- t.test(scm.shannon$shannon ~ scm.shannon$burn_status)

#--Pinaleno Mts.
pm.sd <- t.test(pm.shannon$shannon ~ pm.shannon$burn_status)

#--Add results to DivResults.SequencBased dataframe
DivResults.SequenceBased[DivResults.SequenceBased$div.measure == 'shannon.scm',
                         't.stat'] <- scm.sd$statistic
DivResults.SequenceBased[DivResults.SequenceBased$div.measure == 'shannon.scm',
                         'df'] <- scm.sd$parameter
DivResults.SequenceBased[DivResults.SequenceBased$div.measure == 'shannon.scm',
                         'p.value'] <- scm.sd$p.value
DivResults.SequenceBased[DivResults.SequenceBased$div.measure == 'shannon.pm',
                         't.stat'] <- pm.sd$statistic
DivResults.SequenceBased[DivResults.SequenceBased$div.measure == 'shannon.pm',
                         'df'] <- pm.sd$parameter
DivResults.SequenceBased[DivResults.SequenceBased$div.measure == 'shannon.pm',
                         'p.value'] <- pm.sd$p.value

#----------------------------------------------------------------------------------------#
# Plot of Fisher's alpha by burn status: tree data with singletons----
#----------------------------------------------------------------------------------------#
#--Reorder factors SCM, PM; Unburned, Burned
div.data.out[div.data.out$range == 'santa.catalina', 'range'] <- 'Santa Catalina Mts.'
div.data.out[div.data.out$range == 'pinaleno', 'range'] <- 'Pinaleno Mts.'
div.data.out$range <- factor(div.data.out$range, 
                         levels = c('Santa Catalina Mts.', 'Pinaleno Mts.'))

div.data.out[div.data.out$burn_status == 'burned', 'burn_status'] <- 'Burned'
div.data.out[div.data.out$burn_status == 'unburned', 'burn_status'] <- 'Unburned'
div.data.out$burn_status <- factor(div.data.out$burn_status,
                               levels = c('Unburned','Burned'))

fa.plot <- ggplot(div.data.out, aes(x = burn_status, y = fa,
                                    fill = burn_status)) +
  geom_boxplot() +
  scale_x_discrete(name = "Fire history") +
  facet_grid(~ range) +
  scale_y_continuous(name = "Fisher's alpha") +
  theme_bw() +
  labs(fill = "Fire history") +
  scale_fill_manual(values=c('#006d2c','#636363')) +
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

#----------------------------------------------------------------------------------------#
# Plot of Shannon's diversity by range and burn_status: tree----
#----------------------------------------------------------------------------------------#
#--Reorder factors SCM, PM; Unburned, Burned
div.data[div.data$range == 'santa.catalina', 'range'] <- 'Santa Catalina Mts.'
div.data[div.data$range == 'pinaleno', 'range'] <- 'Pinaleno Mts.'
div.data$range <- factor(div.data$range, 
                            levels = c('Santa Catalina Mts.', 'Pinaleno Mts.'))

div.data[div.data$burn_status == 'burned', 'burn_status'] <- 'Burned'
div.data[div.data$burn_status == 'unburned', 'burn_status'] <- 'Unburned'
div.data$burn_status <- factor(div.data$burn_status,
                                  levels = c('Unburned','Burned'))

shannon.plot <- ggplot(div.data, aes(x = burn_status, y = shannon,
                                     fill = burn_status)) +
  geom_boxplot() +
  facet_grid(~ range) +
  scale_x_discrete(name = "Fire history") +
  scale_y_continuous(name = "Shannon diversity index") +
  theme_bw() +
  scale_fill_brewer(palette = "Accent") +
  labs(fill = "Fire history") +
  scale_fill_manual(values=c('#006d2c','#636363')) +
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

write.csv(DivResults.SequenceBased, 'results_output/DiversityResults.csv',
          row.names = F)



