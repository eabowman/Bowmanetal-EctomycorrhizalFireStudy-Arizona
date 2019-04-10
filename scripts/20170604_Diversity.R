## Script created by Liz Bowman June 4, 2017
## model fitting of diversity between ranges and burn status of sites

#========================================================================================#
# Load data and libraries----
#========================================================================================#

#----------------------------------------------------------------------------------------#
# Load libraries----
#----------------------------------------------------------------------------------------#
library(ggplot2);library(tidyr);library(vegan);library(dplyr); library(car); library(rcompanion)

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

#--reorder matrix
colnames(stsp.matrix)
stsp.matrix <- stsp.matrix[c(1:4,121,5:120)]

#<< w/o singletons: by tree >> ------------------------
nonsingletons <- as.data.frame(which(colSums(stsp.matrix[6:length(stsp.matrix)]) > 1))
nonsingletons <- rownames(nonsingletons)
nonsingletons <- stsp.matrix[nonsingletons]
nonsingletons$Tree <- stsp.matrix$Tree
stsp.wo.singletons <- stsp.matrix[1:5]
stsp.wo.singletons <- left_join(stsp.wo.singletons, nonsingletons, by = NULL)

#<< Creates matrix, grouping OTUs by site, burn status, and range >> -------
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
colnames(rangeburn.matrix)
rangeburn.matrix <- rangeburn.matrix[c(1:3, 120, 4:119)]

#<< w/o singletons: by site >> ------------------------
nonsingletons <- as.data.frame(which(colSums(rangeburn.matrix[5:length(rangeburn.matrix)]) > 1))
nonsingletons <- rownames(nonsingletons)
nonsingletons <- rangeburn.matrix[nonsingletons]
nonsingletons$Site <- rangeburn.matrix$Site
rangeburn.wo.singletons <- rangeburn.matrix[1:4]
rangeburn.wo.singletons <- 
  left_join(rangeburn.wo.singletons, nonsingletons, by = NULL)

#----------------------------------------------------------------------------------------#
# Create dataframe for results----
#----------------------------------------------------------------------------------------#
DivResults.SequenceBased <- data.frame(div.measure = c('fisher.w.range', 'fisher.w.burn',
                                       'fisher.w.interaction','shannon.w.range',
                                       'shannon.w.burn','shannon.w.interaction',
                                       'fisher.wo.range', 'fisher.wo.burn',
                                       'fisher.wo.interaction','shannon.wo.range',
                                       'shannon.wo.burn', 'shannon.wo.interaction'),
                                       df1 = NA,
                                       df2 = NA,
                                       f.stat = NA,
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

# << Table of species abundance: w/ singletons >> ------------------------
#--Overall
spec <- table(colSums(stsp.matrix[6:length(stsp.matrix)]))
plot(spec, xlab = 'Individuals per species',
     ylab = 'Number of species')

#<< Mean species richness: By tree >> ------------
sr.by.tree <- specnumber(stsp.matrix[6:length(stsp.matrix)], 
                         groups = stsp.matrix$Tree)
#mean
mean.sr.by.tree <- mean(sr.by.tree)
SpecRichness.SequenceBased[SpecRichness.SequenceBased$data == 'w/ singletons',
                           'avg'] <- mean.sr.by.tree

#standard deviation
sd.sr.by.tree <- sd(sr.by.tree)
SpecRichness.SequenceBased[SpecRichness.SequenceBased$data == 'w/ singletons',
                           'sd'] <- sd.sr.by.tree

# << Mean species richness: By range and burn >> --------
sr.by.burnrange <- specnumber(stsp.matrix[6:length(stsp.matrix)],
                      groups = stsp.matrix$burn_range)
sr.by.burnrange
anova.burnrange <- 
  aov(sr.by.tree ~ stsp.matrix$Burn_status + stsp.matrix$Range) 
summary.burnrange <- summary(anova.burnrange) # Significantly different

#Add to results dataframe
SpecRichnessResults.SequenceBased[SpecRichnessResults.SequenceBased$data ==
                                    'w/ singletons.burn', 'df1'] <-
  summary.burnrange[[1]]$Df[1]
SpecRichnessResults.SequenceBased[SpecRichnessResults.SequenceBased$data ==
                                    'w/ singletons.burn', 'df2'] <- 
  summary.burnrange[[1]]$Df[3]
SpecRichnessResults.SequenceBased[SpecRichnessResults.SequenceBased$data ==
                                    'w/ singletons.burn', 'f.stat'] <- 
  summary.burnrange[[1]]$`F value`[1]
SpecRichnessResults.SequenceBased[SpecRichnessResults.SequenceBased$data ==
                                    'w/ singletons.burn', 'p.value'] <- 
  summary.burnrange[[1]]$`Pr(>F)`[1]
SpecRichnessResults.SequenceBased[SpecRichnessResults.SequenceBased$data ==
                                    'w/ singletons.range', 'df1'] <-
  summary.burnrange[[1]]$Df[2]
SpecRichnessResults.SequenceBased[SpecRichnessResults.SequenceBased$data ==
                                    'w/ singletons.range', 'df2'] <- 
  summary.burnrange[[1]]$Df[3]
SpecRichnessResults.SequenceBased[SpecRichnessResults.SequenceBased$data ==
                                    'w/ singletons.range', 'f.stat'] <- 
  summary.burnrange[[1]]$`F value`[2]
SpecRichnessResults.SequenceBased[SpecRichnessResults.SequenceBased$data ==
                                    'w/ singletons.range', 'p.value'] <- 
  summary.burnrange[[1]]$`Pr(>F)`[2]


# << Table of species abundance: w/o singletons >> ------------------------
#--Overall
spec <- table(colSums(stsp.wo.singletons[6:length(stsp.wo.singletons)]))
plot(spec, xlab = 'Individuals per species',
     ylab = 'Number of species')

#<< Mean species richness: By tree >> ------------
sr.by.wo.tree <- specnumber(stsp.wo.singletons[6:length(stsp.wo.singletons)], 
                         groups = stsp.wo.singletons$Tree)

mean.sr.by.wo.tree <- mean(sr.by.wo.tree)
SpecRichness.SequenceBased[SpecRichness.SequenceBased$data == 'w/o singletons',
                           'avg'] <- mean.sr.by.wo.tree

sd.sr.by.wo.tree <- sd(sr.by.wo.tree)
SpecRichness.SequenceBased[SpecRichness.SequenceBased$data == 'w/o singletons',
                           'sd'] <- sd.sr.by.wo.tree

# << Mean species richness: By range and burn >> --------
sr.by.wo.burnrange <- specnumber(stsp.wo.singletons[6:length(stsp.wo.singletons)],
                              groups = stsp.wo.singletons$burn_range)
sr.by.wo.burnrange
anova.wo.burnrange <- 
  aov(sr.by.wo.tree ~ stsp.wo.singletons$Burn_status + stsp.wo.singletons$Range) 
summary.wo.burnrange <- summary(anova.wo.burnrange) 

#Add to results dataframe
SpecRichnessResults.SequenceBased[SpecRichnessResults.SequenceBased$data ==
                                    'w/o singletons.burn', 'df1'] <-
  summary.wo.burnrange[[1]]$Df[1]
SpecRichnessResults.SequenceBased[SpecRichnessResults.SequenceBased$data ==
                                    'w/o singletons.burn', 'df2'] <- 
  summary.wo.burnrange[[1]]$Df[3]
SpecRichnessResults.SequenceBased[SpecRichnessResults.SequenceBased$data ==
                                    'w/o singletons.burn', 'f.stat'] <- 
  summary.wo.burnrange[[1]]$`F value`[1]
SpecRichnessResults.SequenceBased[SpecRichnessResults.SequenceBased$data ==
                                    'w/o singletons.burn', 'p.value'] <- 
  summary.wo.burnrange[[1]]$`Pr(>F)`[1]
SpecRichnessResults.SequenceBased[SpecRichnessResults.SequenceBased$data ==
                                    'w/o singletons.range', 'df1'] <-
  summary.wo.burnrange[[1]]$Df[2]
SpecRichnessResults.SequenceBased[SpecRichnessResults.SequenceBased$data ==
                                    'w/o singletons.range', 'df2'] <- 
  summary.wo.burnrange[[1]]$Df[3]
SpecRichnessResults.SequenceBased[SpecRichnessResults.SequenceBased$data ==
                                    'w/o singletons.range', 'f.stat'] <- 
  summary.wo.burnrange[[1]]$`F value`[2]
SpecRichnessResults.SequenceBased[SpecRichnessResults.SequenceBased$data ==
                                    'w/o singletons.range', 'p.value'] <- 
  summary.wo.burnrange[[1]]$`Pr(>F)`[2]


write.csv(SpecRichnessResults.SequenceBased,
          paste0(res.dir, 'SpeciesRichnessResults_SequenceBased.csv'),
          row.names = F)

write.csv(SpecRichness.SequenceBased,
          paste0(res.dir, 'SpeciesRichnessAvg_SequenceBased.csv'),
          row.names = F)
#========================================================================================#
# Diversity: Fisher's alpha -------
#========================================================================================#

#<< By Site >> ------
#--With Singletons included-----
div.site.wSingle.data <- data.frame(site = rangeburn.matrix$Site,
                                    range = rangeburn.matrix$Range,
                                    burn_status = rangeburn.matrix$Burn_status,
                                    burn_range = rangeburn.matrix$burn_range,
                                    fa = fisher.alpha(rangeburn.matrix[5:length(rangeburn.matrix)]),
                                    shannon = diversity(rangeburn.matrix[5:length(rangeburn.matrix)],
                                                        index = 'shannon'))

#--w/o Singletons included-----
div.site.wo.data <- data.frame(site = rangeburn.wo.singletons$Site,
                                    range = rangeburn.wo.singletons$Range,
                                    burn_status = rangeburn.wo.singletons$Burn_status,
                                    burn_range = rep(NA,8),
                                    fa = fisher.alpha(rangeburn.wo.singletons[5:length(rangeburn.wo.singletons)]),
                                    shannon = diversity(rangeburn.wo.singletons[5:length(rangeburn.wo.singletons)],
                                                        index = 'shannon'))
#--Add burn_range column
for(b in c('burned','unburned')){
  for(r in c('pinaleno','santa.catalina')){
    div.site.wo.data[div.site.wo.data$range == r & div.site.wo.data$burn_status == b,
                     'burn_range'] <- paste0(r,'_',b)
  }
}

# << Check for normality: with singletons >> -----------------
#--Fisher's alpha
#Distribution
ggplot(div.site.wSingle.data, aes(x = fa, colour = range)) +
  geom_density()

#residuals
plot(aov(fa ~ site, data = div.site.wSingle.data),1)

#qqplot
plot(aov(fa ~ range, data = div.site.wSingle.data),2)
shapiro.test(log(div.site.wSingle.data$fa))

#Homogeneity of variance
leveneTest(div.site.wSingle.data$fa, div.site.wSingle.data$range,
           center = median)

#--Shannon's diversity index
#Distribution
ggplot(div.site.wSingle.data, aes(x = shannon, colour = range)) +
  geom_density()

plotNormalHistogram(shannon)

#residuals
plot(aov(shannon ~ range, data = div.site.wSingle.data),1)

#qqplot
plot(aov(shannon ~ range, data = div.site.wSingle.data),2)
shapiro.test(log(div.site.wSingle.data$shannon))

#Homogeneity of variance
leveneTest(div.site.wSingle.data$shannon, div.site.wSingle.data$range,
           center = median)

# << Diversity analysis: with singletons >> ------------------
#--Fisher's alpha
#multiple regression
fisher.site.wSingleton.lm <- lm(fa ~ burn_status *range,
                        data = div.site.wSingle.data)
mult.reg.fa.wSingleton <- summary(fisher.site.wSingleton.lm)

#ANOVA
anova.fa.wSingleton <- aov(fisher.site.wSingleton.lm)
fa.wSingleton <- summary(anova.fa.wSingleton)

#Add to results dataframe
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
          'fisher.w.burn', 'df1'] <- fa.wSingleton[[1]]$Df[1]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
          'fisher.w.burn', 'df2'] <- fa.wSingleton[[1]]$Df[4]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
          'fisher.w.burn', 'f.stat'] <- fa.wSingleton[[1]]$`F value`[1]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
          'fisher.w.burn', 'p.value'] <- fa.wSingleton[[1]]$`Pr(>F)`[1]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
          'fisher.w.range', 'df1'] <- fa.wSingleton[[1]]$Df[2]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
          'fisher.w.range', 'df2'] <- fa.wSingleton[[1]]$Df[4]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
          'fisher.w.range', 'f.stat'] <- fa.wSingleton[[1]]$`F value`[2]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
          'fisher.w.range', 'p.value'] <- fa.wSingleton[[1]]$`Pr(>F)`[2]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
          'fisher.w.interaction', 'df1'] <- fa.wSingleton[[1]]$Df[3]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
           'fisher.w.interaction', 'df2'] <- fa.wSingleton[[1]]$Df[4]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
           'fisher.w.interaction', 'f.stat'] <- fa.wSingleton[[1]]$`F value`[3]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
           'fisher.w.interaction', 'p.value'] <- fa.wSingleton[[1]]$`Pr(>F)`[3]

#--Shannon's diversity
#--kruskal wallis test
kruskal.test(shannon ~ burn_range, data = div.site.wSingle.data)

#Add to results dataframe
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'shannon.w.burn', 'df1'] <- fa.wSingleton[[1]]$Df[1]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'shannon.w.burn', 'df2'] <- fa.wSingleton[[1]]$Df[4]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'shannon.w.burn', 'f.stat'] <- fa.wSingleton[[1]]$`F value`[1]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'shannon.w.burn', 'p.value'] <- fa.wSingleton[[1]]$`Pr(>F)`[1]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'shannon.w.range', 'df1'] <- fa.wSingleton[[1]]$Df[2]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'shannon.w.range', 'df2'] <- fa.wSingleton[[1]]$Df[4]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'shannon.w.range', 'f.stat'] <- fa.wSingleton[[1]]$`F value`[2]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'shannon.w.range', 'p.value'] <- fa.wSingleton[[1]]$`Pr(>F)`[2]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'shannon.w.interaction', 'df1'] <- fa.wSingleton[[1]]$Df[3]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'shannon.w.interaction', 'df2'] <- fa.wSingleton[[1]]$Df[4]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'shannon.w.interaction', 'f.stat'] <- fa.wSingleton[[1]]$`F value`[3]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'shannon.w.interaction', 'p.value'] <- fa.wSingleton[[1]]$`Pr(>F)`[3]

# << Check for normality: w/o singletons >> -----------------
#--Fisher's alpha
#Distribution
ggplot(div.site.wo.data, aes(x = fa, colour = range)) +
  geom_density()

#residuals
plot(aov(fa ~ range, data = div.site.wo.data),1)

#qqplot
plot(aov(fa ~ range, data = div.site.wo.data),2)
shapiro.test(log(div.site.wo.data$fa))

#Homogeneity of variance
leveneTest(div.site.wo.data$fa, div.site.wo.data$range,
           center = median)

#--Shannon's diversity index
#Distribution
ggplot(div.site.wo.data, aes(x = shannon, colour = range)) +
  geom_density()

#residuals
plot(aov(shannon ~ range, data = div.site.wo.data),1)

#qqplot
plot(aov(shannon ~ range, data = div.site.wo.data),2)
shapiro.test(log(div.site.wo.data$shannon))

#Homogeneity of variance
leveneTest(div.site.wo.data$shannon, div.site.wo.data$range,
           center = median)

# << Diversity analysis: Site data with singletons >> ------------------
#--Fisher's alpha
#multiple regression
fisher.site.wo.lm <- lm(fa ~ burn_status * range,
                                  data = div.site.wo.data)
mult.reg.fa.wo <- summary(fisher.site.wo.lm)

#ANOVA
anova.fa.wo <- aov(fisher.site.wo.lm)
fa.wo.Singleton <- summary(anova.fa.wo)

#Add to results dataframe
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'fisher.wo.burn', 'df1'] <- fa.wo.Singleton[[1]]$Df[1]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'fisher.wo.burn', 'df2'] <- fa.wo.Singleton[[1]]$Df[4]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'fisher.wo.burn', 'f.stat'] <- fa.wo.Singleton[[1]]$`F value`[1]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'fisher.wo.burn', 'p.value'] <- fa.wo.Singleton[[1]]$`Pr(>F)`[1]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'fisher.wo.range', 'df1'] <- fa.wo.Singleton[[1]]$Df[2]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'fisher.wo.range', 'df2'] <- fa.wo.Singleton[[1]]$Df[4]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'fisher.wo.range', 'f.stat'] <- fa.wo.Singleton[[1]]$`F value`[2]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'fisher.wo.range', 'p.value'] <- fa.wo.Singleton[[1]]$`Pr(>F)`[2]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'fisher.wo.interaction', 'df1'] <- fa.wo.Singleton[[1]]$Df[3]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'fisher.wo.interaction', 'df2'] <- fa.wo.Singleton[[1]]$Df[4]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'fisher.wo.interaction', 'f.stat'] <- fa.wo.Singleton[[1]]$`F value`[3]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'fisher.wo.interaction', 'p.value'] <- fa.wo.Singleton[[1]]$`Pr(>F)`[3]

#--Shannon's diversity
shannon.site.wo.lm <- aov(shannon ~ burn_status * range,
                                   data = div.site.wo.data)
mult.reg.shannon.wo <- summary(shannon.site.wo.lm)

#ANOVA
anova.shannon.wo <- aov(shannon.site.wo.lm)
shannon.wo.Singleton <- summary(anova.shannon.wo)

#Add to results dataframe
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'shannon.wo.burn', 'df1'] <- shannon.wo.Singleton[[1]]$Df[1]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'shannon.wo.burn', 'df2'] <- shannon.wo.Singleton[[1]]$Df[4]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'shannon.wo.burn', 'f.stat'] <- shannon.wo.Singleton[[1]]$`F value`[1]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'shannon.wo.burn', 'p.value'] <- shannon.wo.Singleton[[1]]$`Pr(>F)`[1]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'shannon.wo.range', 'df1'] <- shannon.wo.Singleton[[1]]$Df[2]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'shannon.wo.range', 'df2'] <- shannon.wo.Singleton[[1]]$Df[4]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'shannon.wo.range', 'f.stat'] <- shannon.wo.Singleton[[1]]$`F value`[2]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'shannon.wo.range', 'p.value'] <- shannon.wo.Singleton[[1]]$`Pr(>F)`[2]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'shannon.wo.interaction', 'df1'] <- shannon.wo.Singleton[[1]]$Df[3]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'shannon.wo.interaction', 'df2'] <- shannon.wo.Singleton[[1]]$Df[4]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'shannon.wo.interaction', 'f.stat'] <- shannon.wo.Singleton[[1]]$`F value`[3]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'shannon.wo.interaction', 'p.value'] <- shannon.wo.Singleton[[1]]$`Pr(>F)`[3]

write.csv(DivResults.SequenceBased, paste0(res.dir, 'DivResults_site_SequenceBased.csv'),
          row.names = F)

#<< By Tree >> ------
#--With Singletons included-----
div.tree.wSingle.data <- data.frame(tree = stsp.matrix$Tree,
                                    range = stsp.matrix$Range,
                                    burn_status = stsp.matrix$Burn_status,
                                    burn_range = stsp.matrix$burn_range,
                                    fa = fisher.alpha(stsp.matrix[6:length(stsp.matrix)]),
                                    shannon = diversity(stsp.matrix[6:length(stsp.matrix)], index = 'shannon'))
#--remove outliers (those that have equal species richness and abundance)
fa.out <- c('F11','F12', 'F13', 'F14','F4','LB021','NF19')
div.tree.wSingle.data.o <- div.tree.wSingle.data[!div.tree.wSingle.data$tree %in% fa.out, ]


#--w/o Singletons included-----
div.tree.wo.data <- data.frame(tree = stsp.wo.singletons$Tree,
                               range = stsp.wo.singletons$Range,
                               burn_status = stsp.wo.singletons$Burn_status,
                               burn_range = stsp.wo.singletons$burn_range,
                               fa = fisher.alpha(stsp.wo.singletons[6:length(stsp.wo.singletons)]),
                               shannon = diversity(stsp.wo.singletons[6:length(stsp.wo.singletons)],
                                                   index = 'shannon'))
#--remove outliers (those that have equal species richness and abundance)
fa.out <- c('F11','F12','F13','F14','F4','LB021','LB059','NF19')
div.data.fa <- div.tree.wo.data[!div.tree.wo.data$tree %in% fa.out, ]

# << Diversity analysis: Tree data with singletons >> ------------------
#--Fisher's alpha
#multiple regression
fisher.site.wo.lm <- lm(fa ~ burn_status,
                        data = div.data.fa)
mult.reg.fa.wo <- summary(fisher.site.wo.lm)

#ANOVA
anova.fa.wo <- aov(fisher.site.wo.lm)
fa.wo.Singleton <- summary(anova.fa.wo)
fa.wo.Singleton

#Add to results dataframe
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'fisher.wo.burn', 'df1'] <- fa.wo.Singleton[[1]]$Df[1]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'fisher.wo.burn', 'df2'] <- fa.wo.Singleton[[1]]$Df[4]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'fisher.wo.burn', 'f.stat'] <- fa.wo.Singleton[[1]]$`F value`[1]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'fisher.wo.burn', 'p.value'] <- fa.wo.Singleton[[1]]$`Pr(>F)`[1]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'fisher.wo.range', 'df1'] <- fa.wo.Singleton[[1]]$Df[2]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'fisher.wo.range', 'df2'] <- fa.wo.Singleton[[1]]$Df[4]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'fisher.wo.range', 'f.stat'] <- fa.wo.Singleton[[1]]$`F value`[2]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'fisher.wo.range', 'p.value'] <- fa.wo.Singleton[[1]]$`Pr(>F)`[2]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'fisher.wo.interaction', 'df1'] <- fa.wo.Singleton[[1]]$Df[3]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'fisher.wo.interaction', 'df2'] <- fa.wo.Singleton[[1]]$Df[4]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'fisher.wo.interaction', 'f.stat'] <- fa.wo.Singleton[[1]]$`F value`[3]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'fisher.wo.interaction', 'p.value'] <- fa.wo.Singleton[[1]]$`Pr(>F)`[3]

#--Shannon's diversity
shannon.site.wo.lm <- aov(shannon ~ burn_status * range,
                          data = div.site.wo.data)
mult.reg.shannon.wo <- summary(shannon.site.wo.lm)

#ANOVA
anova.shannon.wo <- aov(shannon.site.wo.lm)
shannon.wo.Singleton <- summary(anova.shannon.wo)

#Add to results dataframe
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'shannon.wo.burn', 'df1'] <- shannon.wo.Singleton[[1]]$Df[1]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'shannon.wo.burn', 'df2'] <- shannon.wo.Singleton[[1]]$Df[4]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'shannon.wo.burn', 'f.stat'] <- shannon.wo.Singleton[[1]]$`F value`[1]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'shannon.wo.burn', 'p.value'] <- shannon.wo.Singleton[[1]]$`Pr(>F)`[1]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'shannon.wo.range', 'df1'] <- shannon.wo.Singleton[[1]]$Df[2]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'shannon.wo.range', 'df2'] <- shannon.wo.Singleton[[1]]$Df[4]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'shannon.wo.range', 'f.stat'] <- shannon.wo.Singleton[[1]]$`F value`[2]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'shannon.wo.range', 'p.value'] <- shannon.wo.Singleton[[1]]$`Pr(>F)`[2]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'shannon.wo.interaction', 'df1'] <- shannon.wo.Singleton[[1]]$Df[3]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'shannon.wo.interaction', 'df2'] <- shannon.wo.Singleton[[1]]$Df[4]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'shannon.wo.interaction', 'f.stat'] <- shannon.wo.Singleton[[1]]$`F value`[3]
DivResults.SequenceBased[DivResults.SequenceBased$div.measure ==
                           'shannon.wo.interaction', 'p.value'] <- shannon.wo.Singleton[[1]]$`Pr(>F)`[3]

write.csv(DivResults.SequenceBased, paste0(res.dir, 'DivResults_tree_SequenceBased.csv'),
          row.names = F)

#----------------------------------------------------------------------------------------#
# Plot of Fisher's alpha by burn status: tree data with singletons----
#----------------------------------------------------------------------------------------#
levels(div.tree.wSingle.data.o$range) <- c('Pinaleno Mts.', 'Santa Catalina Mts.')
levels(div.tree.wSingle.data.o$burn_status) <- c('Burned', 'Unburned')

fa.plot <- ggplot(div.tree.wSingle.data.o, aes(x = burn_status, y = fa,
                                               fill = burn_status)) +
  geom_boxplot() +
  scale_x_discrete(name = "Fire history") +
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
ggsave('Fig2a.jpeg', plot = fa.plot, device = 'jpeg',
       path = '/Volumes/Cenococcum/PhD/Dissertation/Chpt.2/Figures/',
       width = 20, height = 20, units = 'cm')

#----------------------------------------------------------------------------------------#
# Plot of Fisher's alpha by range and burn_status: tree data with singletons----
#----------------------------------------------------------------------------------------#
levels(div.tree.wSingle.data.o$range) <- c('Pinaleno Mts.', 'Santa Catalina Mts.')
levels(div.tree.wSingle.data.o$burn_status) <- c('Burned', 'Unburned')
 
fa.plot <- ggplot(div.tree.wSingle.data.o, aes(x = burn_status, y = fa,
                                  fill = burn_status)) +
  geom_boxplot() +
  scale_x_discrete(name = "Fire history") +
  scale_y_continuous(name = "Fisher's alpha") +
  facet_grid(. ~ range) +
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
ggsave('Fig2b.jpeg', plot = fa.plot, device = 'jpeg',
       path = '/Volumes/Cenococcum/PhD/Dissertation/Chpt.2/Figures/',
       width = 20, height = 20, units = 'cm')

#----------------------------------------------------------------------------------------#
# Plot of Shannon's diversity by range and burn_status: tree----
#----------------------------------------------------------------------------------------#
levels(div.tree.wSingle.data$range) <- c('Pinaleno Mts.', 'Santa Catalina Mts.')
levels(div.tree.wSingle.data$burn_status) <- c('Burned', 'Unburned')

shannon.plot <- ggplot(div.tree.wSingle.data, aes(x = burn_status, y = shannon,
                                  fill = burn_status)) +
  geom_boxplot() +
  scale_x_discrete(name = "Fire history") +
  scale_y_continuous(name = "Shannon diversity index") +
  facet_grid(. ~ range) +
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
ggsave('Fig2c.jpeg', plot = shannon.plot,
       device = 'jpeg', path = '/Volumes/Cenococcum/PhD/Dissertation/Chpt.2/Figures/',
       width = 20, height = 20, units = 'cm')


#========================================================================================#
# Assess normality----
#========================================================================================#
norm.data <- data.frame(tree = stsp.matrix$Tree, 
                        site = stsp.matrix$Site,
                        range = stsp.matrix$Range,
                        burn_status = stsp.matrix$Burn_status,
                        spec.richness = NA)
#--Species richness
norm.data$spec.richness <- specnumber(stsp.matrix[6:length(stsp.matrix)], 
                                      groups = stsp.matrix$Tree)
normtest(norm.data, 'spec.richness')

#--Fisher's alpha


#--Shannon's diversity
