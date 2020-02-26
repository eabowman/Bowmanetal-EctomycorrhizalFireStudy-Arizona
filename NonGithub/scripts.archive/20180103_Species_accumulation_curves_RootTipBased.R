## Script created by Liz Bowman June 5, 2017
## for making species accumulation curves

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
fig.dir <- '~/Documents/PhD/2_EM_Fire_effect/figures_output/'
res.dir <- "~/Documents/PhD/2_EM_Fire_effect/results_output/"

#-----------------------------------------------------------------------------------------
# Load data and clean up
# Make two data frames: one with singletons and one without singletons
#-----------------------------------------------------------------------------------------
#--Load data
all.data <- read.csv(paste0(dat.dir,'20170806_OTU_data.csv'), as.is = T)

#--Remove sequences not assigned to an OTU
otu.data <- all.data[!is.na(all.data$otu.97), ]

#--Remove sequences not assigned to Ponderosa host
otu.data <- otu.data[otu.data$Host == 'Ponderosa',]

#--Creates matrix, grouping OTUs by tree number and getting a sum of EM tip abundance 
# for each OTU
otu.data %>%
  select (Sample_name,Burn_status,Range,Site,Tree,otu.97,Tip_count) %>%
  spread (otu.97,Tip_count) %>%
  group_by(Tree,Site,Range,Burn_status) %>%
  select (matches ("otu.*")) %>%
  summarize_each (funs (sum (., na.rm = TRUE))) %>%
  as.data.frame () -> stsp.matrix

#--<< Combine range and burn status >> ----------
for(r in unique(stsp.matrix$Range)) {
  for(b in unique(stsp.matrix$Burn_status)) {
    stsp.matrix[stsp.matrix$Burn_status == b & stsp.matrix$Range == r, 'burn_range'] <-
      paste(b, r)
  }
}

colnames(stsp.matrix)
stsp.matrix <- stsp.matrix[c(1:4,121,5:120)]

# #=========================================================================================
# # Species accumulation curves by burn status----
# #=========================================================================================
# 
# #<< burned sites >> -----------------
# #--With singletons
# burned.sp <- stsp.matrix[stsp.matrix$Burn_status == 'burned', 6:length(stsp.matrix)]
# #--remove columns with no species occurrences
# burned.sp <- burned.sp[colSums(burned.sp) > 0]
# 
# #--Without singletons
# burned.sp.wo <- stsp.matrix[stsp.matrix$Burn_status == 'burned', 6:length(stsp.matrix)]
# #--remove columns with single species occurrences
# burned.sp.wo <- burned.sp.wo[colSums(burned.sp.wo) > 1]
# #--remove rows with no species occurrences
# burned.sp.wo <- burned.sp.wo[rowSums(burned.sp.wo) > 0]
# 
# #<< unburned sites >> -----------------
# unburned.sp <- stsp.matrix[stsp.matrix$Burn_status == 'unburned', 6:length(stsp.matrix)]
# #--remove columns with no species occurrences
# unburned.sp <- unburned.sp[which(colSums(unburned.sp) > 0)]
# 
# #--Without singletons
# unburned.sp.wo <- stsp.matrix[stsp.matrix$Burn_status == 'unburned',
#                               6:length(stsp.matrix)]
# #--remove columns with single species occurrences
# unburned.sp.wo <- unburned.sp.wo[colSums(unburned.sp.wo) > 1]
# #--remove rows with no species occurrences
# burned.sp.wo <- unburned.sp.wo[rowSums(unburned.sp.wo) > 0]
# 
# #=========================================================================================
# # Species accumulation curves by range----
# #=========================================================================================
# 
# #<< Pinaleno sites >> -----------------
# #--With singletons
# pinaleno.sp <- stsp.matrix[stsp.matrix$Range == 'pinaleno', 6:length(stsp.matrix)]
# #--remove columns with no species occurrences
# pinaleno.sp <- pinaleno.sp[colSums(pinaleno.sp) > 0]
# 
# #--Without singletons
# pinaleno.sp.wo <- stsp.matrix[stsp.matrix$Range == 'pinaleno',
#                               6:length(stsp.matrix)]
# #--remove columns with single species occurrences
# pinaleno.sp.wo <- pinaleno.sp.wo[colSums(pinaleno.sp.wo) > 1]
# #--remove rows with no species occurrences
# pinaleno.sp.wo <- pinaleno.sp.wo[rowSums(pinaleno.sp.wo) > 0]
# 
# #<< Santa Catalina sites >> -----------------
# #--With singletons
# sc.sp <- stsp.matrix[stsp.matrix$Range == 'santa.catalina',
#                            6:length(stsp.matrix)]
# #--remove columns with no species occurrences
# sc.sp <- sc.sp[which(colSums(sc.sp) > 0)]
# 
# #--Without singletons
# sc.sp.wo <- stsp.matrix[stsp.matrix$Range == 'santa.catalina',
#                         6:length(stsp.matrix)]
# #--remove columns with single species occurrences
# sc.sp.wo <- sc.sp.wo[colSums(sc.sp.wo) > 1]
# #--remove rows with no species occurrences
# sc.sp.wo <- sc.sp.wo[rowSums(sc.sp.wo) > 0]

#=========================================================================================
# Species accumulation curves by range and burn_status----
#=========================================================================================

#<< Pinaleno burned sites >> -----------------------------------------
#--With singletons
PMburned.sp <- stsp.matrix[stsp.matrix$burn_range == 'burned pinaleno', 6:length(stsp.matrix)]
#--remove columns with no species occurrences
PMburned.sp <- PMburned.sp[colSums(PMburned.sp) > 0]

#--Without singletons
PMburned.sp.wo <- stsp.matrix[stsp.matrix$burn_range == 'burned pinaleno', 6:length(stsp.matrix)]
#--remove columns with single species occurrences
PMburned.sp.wo <- PMburned.sp.wo[colSums(PMburned.sp.wo) > 1]
#--remove rows with no species occurrences
PMburned.sp.wo <- PMburned.sp.wo[rowSums(PMburned.sp.wo) > 0]

#<< Pinaleno unburned sites >> -----------------------------------------
PMunburned.sp <- stsp.matrix[stsp.matrix$burn_range == 'unburned pinaleno', 6:length(stsp.matrix)]
#--remove columns with no species occurrences
PMunburned.sp <- PMunburned.sp[which(colSums(PMunburned.sp) > 0)]

#--Without singletons
PMunburned.sp.wo <- stsp.matrix[stsp.matrix$burn_range == 'unburned pinaleno',
                                6:length(stsp.matrix)]
#--remove columns with single species occurrences
PMunburned.sp.wo <- PMunburned.sp.wo[colSums(PMunburned.sp.wo) > 1]
#--remove rows with no species occurrences
PMunburned.sp.wo <- PMunburned.sp.wo[rowSums(PMunburned.sp.wo) > 0]

#<< Santa Catalina burned sites >> -----------------------------------------
#--With singletons
SCMburned.sp <- stsp.matrix[stsp.matrix$burn_range == 'burned santa.catalina',
                            6:length(stsp.matrix)]
#--remove columns with no species occurrences
SCMburned.sp <- SCMburned.sp[colSums(SCMburned.sp) > 0]

#--Without singletons
SCMburned.sp.wo <- stsp.matrix[stsp.matrix$burn_range == 'burned santa.catalina',
                               6:length(stsp.matrix)]
#--remove columns with single species occurrences
SCMburned.sp.wo <- SCMburned.sp.wo[colSums(SCMburned.sp.wo) > 1]
#--remove rows with no species occurrences
SCMburned.sp.wo <- SCMburned.sp.wo[rowSums(SCMburned.sp.wo) > 0]

#<< Santa Catalina unburned sites >> -----------------------------------------
SCMunburned.sp <- stsp.matrix[stsp.matrix$burn_range == 'unburned santa.catalina',
                              6:length(stsp.matrix)]
#--remove columns with no species occurrences
SCMunburned.sp <- SCMunburned.sp[which(colSums(SCMunburned.sp) > 0)]

#--Without singletons
SCMunburned.sp.wo <- stsp.matrix[stsp.matrix$burn_range == 'unburned santa.catalina',
                                 6:length(stsp.matrix)]
#--remove columns with single species occurrences
SCMunburned.sp.wo <- SCMunburned.sp.wo[colSums(SCMunburned.sp.wo) > 1]
#--remove rows with no species occurrences
SCMunburned.sp.wo <- SCMunburned.sp.wo[rowSums(SCMunburned.sp.wo) > 0]

#=========================================================================================
# Plot----
#=========================================================================================

#<< Combine Pinaleno burned plots, with singletons and without >>--------------------------
Pin.w.burned <- specaccum(PMburned.sp, sample = min(rowSums(PMburned.sp),  permutations = 999))
Pin.w.burned.df <- data.frame(Sites=Pin.w.burned$sites,
                              Richness=Pin.w.burned$richness,
                              SD=Pin.w.burned$sd)
Pin.wo.burned <- specaccum(PMburned.sp.wo, sample = min(rowSums(PMburned.sp.wo),  permutations = 999))
Pin.wo.burned.df <- data.frame(Sites=Pin.wo.burned$sites,
                               Richness=Pin.wo.burned$richness,
                               SD=Pin.wo.burned$sd)

#<< Combine Pinaleno unburned plots, with singletons and without >>------------------------
Pin.w.unburned <- specaccum(PMunburned.sp, sample = min(rowSums(PMunburned.sp),
                                                        permutations = 999))
Pin.w.unburned.df <- data.frame(Sites=Pin.w.unburned$sites,
                                Richness=Pin.w.unburned$richness,
                                SD=Pin.w.unburned$sd)
Pin.wo.unburned <- specaccum(PMunburned.sp.wo, sample = min(rowSums(PMunburned.sp.wo),
                                                            permutations = 999))
Pin.wo.unburned.df <- data.frame(Sites=Pin.wo.unburned$sites,
                                 Richness=Pin.wo.unburned$richness,
                                 SD=Pin.wo.unburned$sd)

#<< Combine Santa Catalina burned plots, with singletons and without >>---------------------
SCM.w.burned <- specaccum(SCMburned.sp, sample = min(rowSums(SCMburned.sp), permutations = 999))
SCM.w.burned.df <- data.frame(Sites=SCM.w.burned$sites,
                              Richness=SCM.w.burned$richness,
                              SD=SCM.w.burned$sd)
SCM.wo.burned <- specaccum(SCMburned.sp.wo, sample = min(rowSums(SCMburned.sp.wo),permutations = 999))
SCM.wo.burned.df <- data.frame(Sites=SCM.wo.burned$sites,
                               Richness=SCM.wo.burned$richness,
                               SD=SCM.wo.burned$sd)

#<< Combine Santa Catalina unburned plots, with singletons and without >>-------------------
SCM.w.unburned <- specaccum(SCMunburned.sp, sample = min(rowSums(SCMunburned.sp),permutations = 999))
SCM.w.unburned.df <- data.frame(Sites=SCM.w.unburned$sites,
                                Richness=SCM.w.unburned$richness,
                                SD=SCM.w.unburned$sd)
SCM.wo.unburned <- specaccum(SCMunburned.sp.wo, sample = min(rowSums(SCMunburned.sp.wo),
                                                             permutations = 999))
SCM.wo.unburned.df <- data.frame(Sites=SCM.wo.unburned$sites,
                                 Richness=SCM.wo.unburned$richness,
                                 SD=SCM.wo.unburned$sd)

#<< Pinaleno burned plot >> ------------
Pinaleno.burned <- ggplot() +
  geom_point(data=Pin.w.burned.df, aes(x=Sites, y=Richness), size = 3) +
  geom_line(data=Pin.w.burned.df, aes(x=Sites, y=Richness)) +
  geom_ribbon(data=Pin.w.burned.df ,aes(x=Sites,
                                        ymin=(Richness-2*SD),
                                        ymax=(Richness+2*SD)),
              alpha=0.2) +
  geom_point(data=Pin.wo.burned.df, aes(x=Sites, y=Richness), colour = 'darkgrey', size = 3) +
  geom_line(data=Pin.wo.burned.df, aes(x=Sites, y=Richness), colour = 'darkgrey') +
  geom_ribbon(data=Pin.wo.burned.df ,aes(x=Sites,
                                         ymin=(Richness-2*SD),
                                         ymax=(Richness+2*SD)),
              alpha=0.2) +
  theme_bw() +
  expand_limits(y=c(0,50), x = c(0,10)) +
  ylab('OTUs') +
  xlab('Trees sampled') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size=22, color = 'black'),
        axis.title = element_text(size = 28), legend.position = "none") +
  theme(axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)))

ggsave('SpecAccum_PinBurn_RootTipBased.tiff', plot = Pinaleno.burned,
       device = 'tiff', path = fig.dir,
       width = 20, height = 20, units = 'cm')

#<< Pinaleno unburned plot >> ------------
Pinaleno.unburned <- ggplot() +
  geom_point(data=Pin.w.unburned.df, aes(x=Sites, y=Richness), size = 3) +
  geom_line(data=Pin.w.unburned.df, aes(x=Sites, y=Richness)) +
  geom_ribbon(data=Pin.w.unburned.df ,aes(x=Sites,
                                          ymin=(Richness-2*SD),
                                          ymax=(Richness+2*SD)),
              alpha=0.2) +
  geom_point(data=Pin.wo.unburned.df, aes(x=Sites, y=Richness), colour = 'darkgrey', size = 3) +
  geom_line(data=Pin.wo.unburned.df, aes(x=Sites, y=Richness), colour = 'darkgrey') +
  geom_ribbon(data=Pin.wo.unburned.df ,aes(x=Sites,
                                           ymin=(Richness-2*SD),
                                           ymax=(Richness+2*SD)),
              alpha=0.2) +
  theme_bw() +
  expand_limits(y=c(0,50), x = c(0,10)) +
  ylab('OTUs') +
  xlab('Trees sampled') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size=22, color = 'black'),
        axis.title = element_text(size = 28), legend.position = "none") +
  theme(axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)))

ggsave('SpecAccum_PinUnburn_RootTipBased.tiff', plot = Pinaleno.unburned,
       device = 'tiff', path = fig.dir,
       width = 20, height = 20, units = 'cm')

#<< Santa Catalina burned plot >> ------------
Catalina.burned <- ggplot() +
  geom_point(data=SCM.w.burned.df, aes(x=Sites, y=Richness), size = 3) +
  geom_line(data=SCM.w.burned.df, aes(x=Sites, y=Richness)) +
  geom_ribbon(data=SCM.w.burned.df ,aes(x=Sites,
                                        ymin=(Richness-2*SD),
                                        ymax=(Richness+2*SD)),
              alpha=0.2) +
  geom_point(data=SCM.wo.burned.df, aes(x=Sites, y=Richness), colour = 'darkgrey', size = 3) +
  geom_line(data=SCM.wo.burned.df, aes(x=Sites, y=Richness), colour = 'darkgrey') +
  geom_ribbon(data=SCM.wo.burned.df ,aes(x=Sites,
                                         ymin=(Richness-2*SD),
                                         ymax=(Richness+2*SD)),
              alpha=0.2) +
  theme_bw() +
  expand_limits(y=c(0,50), x = c(0,10)) +
  ylab('OTUs') +
  xlab('Trees sampled') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size=22, color = 'black'),
        axis.title = element_text(size = 28), legend.position = "none") +
  theme(axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)))

ggsave('SpecAccum_ScmBurn_RootTipBased.tiff', plot = Catalina.burned,
       device = 'tiff', path = fig.dir,
       width = 20, height = 20, units = 'cm')

#<< Santa Catalina unburned plot >> ------------
Catalina.unburned <- ggplot() +
  geom_point(data=SCM.w.unburned.df, aes(x=Sites, y=Richness), size = 3) +
  geom_line(data=SCM.w.unburned.df, aes(x=Sites, y=Richness)) +
  geom_ribbon(data=SCM.w.unburned.df ,aes(x=Sites,
                                          ymin=(Richness-2*SD),
                                          ymax=(Richness+2*SD)),
              alpha=0.2) +
  geom_point(data=SCM.wo.unburned.df, aes(x=Sites, y=Richness), colour = 'darkgrey', size = 3) +
  geom_line(data=SCM.wo.unburned.df, aes(x=Sites, y=Richness), colour = 'darkgrey') +
  geom_ribbon(data=SCM.wo.unburned.df ,aes(x=Sites,
                                           ymin=(Richness-2*SD),
                                           ymax=(Richness+2*SD)),
              alpha=0.2) +
  theme_bw() +
  expand_limits(y=c(0,50), x = c(0,10)) +
  ylab('OTUs') +
  xlab('Trees sampled') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size=22, color = 'black'),
        axis.title = element_text(size = 28), legend.position = "none") +
  theme(axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)))

ggsave('SpecAccum_ScmUnburn_RootTipBased.tiff', plot = Catalina.unburned,
       device = 'tiff', path = fig.dir,
       width = 20, height = 20, units = 'cm')
