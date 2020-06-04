## Script created by Liz Bowman May 9, 2019
## analyses of species richness and diversity

#========================================================================================#
# Load data ----
#========================================================================================#

stsp.matrix.tip <- read.csv(paste0(dat.dir,'97%_SitexSpecies_TipAb.csv'),
                            as.is = T)

#========================================================================================#
# Species richness----
#========================================================================================#
#--Make species richness table
stsp.matrix.tip <- stsp.matrix.tip[order(stsp.matrix.tip$Tree),]

sr.by.tree <- data.frame(tree = stsp.matrix.tip$Tree,
                         site = stsp.matrix.tip$Site,
                         range = stsp.matrix.tip$Range,
                         fire.history = stsp.matrix.tip$Burn_status,
                         spec.richness = specnumber(stsp.matrix.tip[12:length(stsp.matrix.tip)], 
                                                        groups = stsp.matrix.tip$Tree),
                         shannon = diversity(stsp.matrix.tip[12:length(stsp.matrix.tip)],
                                             index = 'shannon'))
#------------------------------------------------------------------------------------------------------#
# Shannon's diversity-----
#------------------------------------------------------------------------------------------------------#

#<< Multiple linear regression of differences in species richness based on range and fire history >>-----

# Remove outlier
shannon.sr.by.tree <- sr.by.tree[!sr.by.tree$tree %in% c('LB021'),]

# ANOVA analysing diversity as a function of fire history, range, and their interaction
shannon.model <- lm(shannon ~ fire.history * range, data = shannon.sr.by.tree)
anova(shannon.model)

#--T-test of each range separately
# Pinaleno Mts.
pinaleno.shannon <- shannon.sr.by.tree[shannon.sr.by.tree$range == 'pinaleno',]
t.pin.shannon <- t.test(shannon ~ fire.history, data = pinaleno.shannon,
                        alternative = 'greater')

# Santa Catalina Mts.
santacat.shannon <- shannon.sr.by.tree[!shannon.sr.by.tree$range == 'pinaleno',]
t.sc.shannon <- t.test(shannon ~ fire.history, data = santacat.shannon,
                       alternative = 'greater')

#------------------------------------------------------------------------------------------------------#
# Species richness-----
#------------------------------------------------------------------------------------------------------#

#<< Multiple linear regression of differences in species richness based on range and fire history >>-----

# ANOVA analysing diversity as a function of fire history, range, and their interaction
sr.model <- lm(spec.richness ~ fire.history * range, data = sr.by.tree)
anova(sr.model)

#--T-test of each range separately
# Pinaleno Mts.
pinaleno.sr <- sr.by.tree[sr.by.tree$range == 'pinaleno',]
t.pin.sr <- t.test(spec.richness ~ fire.history, data = pinaleno.sr,
                        alternative = 'greater')

# Santa Catalina Mts.
santacat.sr <- sr.by.tree[sr.by.tree$range == 'santa.catalina',]
t.sc.sr <- t.test(spec.richness ~ fire.history, data = santacat.sr,
                       alternative = 'greater')

#----------------------------------------------------------------------------------------#
# Plot of Species richness by fire history: tree----
#----------------------------------------------------------------------------------------#
levels(sr.by.tree$range) <- c('Pinaleno Mts.','Santa Catalina Mts.')
levels(sr.by.tree$fire.history) <- c('Burned','Unburned')

sr.plot <- ggplot(sr.by.tree, 
                  aes(x = fire.history,
                      y = spec.richness)) +
  geom_boxplot() +
  scale_x_discrete(name = "Fire history") +
  scale_y_continuous(name = "Species richness") +
  theme_classic() +
  labs(fill = "Fire history") +
  #scale_fill_manual(values=c('grey','grey')) +
  theme(legend.position="none",
        axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size=22, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14)) +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

sr.plot
ggsave('Fig3A.jpeg', plot = sr.plot,
       device = 'jpeg', path = fig.dir,
       width = 20, height = 15, units = 'cm')

#----------------------------------------------------------------------------------------#
# Plot of Shannon's diversity by fire history ----
#----------------------------------------------------------------------------------------#
levels(sr.by.tree$range) <- c('Pinaleno Mts.','Santa Catalina Mts.')
shannon.sr.by.tree <- sr.by.tree[!sr.by.tree$tree %in% c('F13','LB021'),]

shannon.plot <- ggplot(shannon.sr.by.tree,
                       aes(x = fire.history,
                           y = shannon)) +
  geom_boxplot() +
  scale_x_discrete(name = "Fire history") +
  scale_y_continuous(name = "Shannon's diversity") +
  theme_classic() +
  labs(fill = "Fire history") +
  #scale_fill_manual(values=c('grey','grey')) +
  theme(legend.position="none",
        axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size=22, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14)) +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

shannon.plot

ggsave('Fig3B_Range.jpeg', plot = shannon.plot,
       device = 'jpeg', path = fig.dir,
       width = 20, height = 15, units = 'cm')
