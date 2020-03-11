## Script created by Liz Bowman February 24, 2020
## for analyzing taxonomic differences between ranges and burn status of sites

#========================================================================================#
# Load data ----
#========================================================================================#
tax.data <- read.csv(paste0(dat.dir,'OTUdata.csv'), as.is = T)

#<< Results dataframe >> -------------------
tax.results <- data.frame(tests = c('burn.range.class','burn.class',
                                    'range.class','unburned.class'),
                          chi.stat = NA,
                          df = NA,
                          p.value = NA)

#<< Make data frame with root tip data: genus level >> ------------------------
tax.data %>%
  dplyr::select(Sample_name, Range, Burn_status, Taxonomy_genus, Tip_count) %>%
  filter(!is.na(Taxonomy_genus),
         !Taxonomy_genus %in% c('Piloderma','Hygrophorus','Coltricia','Cantharellus',
                               'Oidiodendron','Clavulina','Cortinarius','Archaerhizomyces','Elaphomyces',
                               'Ramaria','Boletus','Hysterangium','Tulasnella','Laccaria','Archaeorhizomyces')) %>%
  group_by(Range, Burn_status, Taxonomy_genus) %>%
  summarize(total_tip = sum(Tip_count)) %>%
  spread(key = Taxonomy_genus, value = total_tip, fill = 0) -> genus.tax

#========================================================================================#
# Chi-square: class level
#========================================================================================#
#----------------------------------------------------------------------------------------#
# Range and Fire history
#----------------------------------------------------------------------------------------#
#<< Make data frame with root tip data: class level >> ------------------------
tax.data %>%
  dplyr::select(Sample_name, Range, Burn_status, Taxonomy_class, Tip_count) %>%
  filter(!is.na(Taxonomy_class),
         !Taxonomy_class %in% c('Saccharomycetes','Archaeorhizomycetes','Eurotiomycetes')) %>%
  group_by(Range, Burn_status, Taxonomy_class) %>%
  summarize(total_tip = sum(Tip_count)) %>%
  spread(key = Taxonomy_class, value = total_tip, fill = 0) %>%
  data.frame() -> class.tax.br

#--chi square test
class.tax.chi <- class.tax.br[-c(1,2)]
rownames(class.tax.chi) <- c('p.b','p.u','s.b','s.u')
tip.class <- chisq.test(class.tax.chi)
tip.class
tax.results[tax.results$tests == 'burn.range.class', 'chi.stat'] <- tip.class$statistic[[1]]
tax.results[tax.results$tests == 'burn.range.class', 'df'] <- tip.class$parameter[[1]]
tax.results[tax.results$tests == 'burn.range.class', 'p.value'] <- tip.class$p.value[[1]]

#----------------------------------------------------------------------------------------#
# Range 
#----------------------------------------------------------------------------------------#
#<< Make data frame with root tip data: class level >> ------------------------
tax.data %>%
  dplyr::select(Sample_name, Range, Taxonomy_class, Tip_count) %>%
  filter(!is.na(Taxonomy_class),
         !Taxonomy_class %in% c('Saccharomycetes','Archaeorhizomycetes','Eurotiomycetes')) %>%
  group_by(Range, Taxonomy_class) %>%
  summarize(total_tip = sum(Tip_count)) %>%
  spread(key = Taxonomy_class, value = total_tip, fill = 0) %>%
  data.frame() -> class.tax.r

#--chi square test
class.tax.r.chi <- class.tax.r[-1]
tip.class <- chisq.test(class.tax.r.chi)
tip.class
tax.results[tax.results$tests == 'range.class', 'chi.stat'] <- tip.class$statistic[[1]]
tax.results[tax.results$tests == 'range.class', 'df'] <- tip.class$parameter[[1]]
tax.results[tax.results$tests == 'range.class', 'p.value'] <- tip.class$p.value[[1]]

#----------------------------------------------------------------------------------------#
# Fire history
#----------------------------------------------------------------------------------------#
#<< Make data frame with root tip data: class level >> ------------------------
tax.data %>%
  dplyr::select(Sample_name, Burn_status, Taxonomy_class, Tip_count) %>%
  filter(!is.na(Taxonomy_class),
         !Taxonomy_class %in% c('Saccharomycetes','Archaeorhizomycetes','Eurotiomycetes')) %>%
  group_by(Burn_status, Taxonomy_class) %>%
  summarize(total_tip = sum(Tip_count)) %>%
  spread(key = Taxonomy_class, value = total_tip, fill = 0) %>%
  data.frame() -> class.tax.b

#--chi square test
class.tax.b.chi <- class.tax.b[-1]
tip.class <- chisq.test(class.tax.b.chi)
tip.class
tax.results[tax.results$tests == 'burn.class', 'chi.stat'] <- tip.class$statistic[[1]]
tax.results[tax.results$tests == 'burn.class', 'df'] <- tip.class$parameter[[1]]
tax.results[tax.results$tests == 'burn.class', 'p.value'] <- tip.class$p.value[[1]]

#----------------------------------------------------------------------------------------#
# Unburned plots only across both ranges
#----------------------------------------------------------------------------------------#
#<< Make data frame with root tip data: class level >> ------------------------
tax.data %>%
  dplyr::select(Sample_name, Range, Burn_status, Taxonomy_class, Tip_count) %>%
  filter(!is.na(Taxonomy_class),
         !Taxonomy_class %in% c('Saccharomycetes','Archaeorhizomycetes','Eurotiomycetes'),
         Burn_status == 'unburned') %>%
  group_by(Range, Burn_status, Taxonomy_class) %>%
  summarize(total_tip = sum(Tip_count)) %>%
  spread(key = Taxonomy_class, value = total_tip, fill = 0) %>%
  data.frame() -> class.tax.unburned

#--chi square test
class.tax.unburned.chi <- class.tax.unburned[-c(1,2)]
tip.class <- chisq.test(class.tax.unburned.chi)
tip.class
tax.results[tax.results$tests == 'unburned.class', 'chi.stat'] <- tip.class$statistic[[1]]
tax.results[tax.results$tests == 'unburned.class', 'df'] <- tip.class$parameter[[1]]
tax.results[tax.results$tests == 'unburned.class', 'p.value'] <- tip.class$p.value[[1]]

write.csv(tax.results, paste0(res.dir, 'TaxonomicAnalysisResults_ChiSquare.csv'),
          row.names = F)

#========================================================================================#
# Plots
#========================================================================================#
#<< Make data frame with root tip data: class level >> ------------------------
tax.data %>%
  dplyr::select(Sample_name, Range, Burn_status, Taxonomy_class, Tip_count) %>%
  filter(!is.na(Taxonomy_class),!Taxonomy_class %in% c('Saccharomycetes','Archaeorhizomycetes')) %>%
  group_by(Range, Burn_status, Taxonomy_class) %>%
  summarize(total_tip = sum(Tip_count)) -> class.tax.plot

#Change levels of Burn_status and Range columns
class.tax.plot$Burn_status <- factor(class.tax.plot$Burn_status,
                                     levels = c('burned','unburned'),
                                     labels = c('Burned', 'Unburned'))
class.tax.plot$Range <- factor(class.tax.plot$Range,
                               levels = c('pinaleno','santa.catalina'),
                               labels = c('Pinaleno Mts.', 'Santa Catalina Mts.'))

#----------------------------------------------------------------------------------------#
# Range and Fire history
#----------------------------------------------------------------------------------------#
#--Bar graph
Tip.class <- ggplot(data = class.tax.plot, 
                          aes(x = Burn_status,
                              y = total_tip,
                              fill = Taxonomy_class)) + 
  geom_col(position = 'fill') +
  ylab("Proportion of \n sequences per class") +
  theme_classic() +
  xlab('Fire history') +
  facet_grid(. ~ Range) +
  scale_fill_manual(values = c('#7b3294','#c2a5cf','#f7f7f7','#a6dba0','#008837')) +
  theme(legend.position = 'right',
        axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size=22, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14),
        legend.text=element_text(size=16),
        legend.title=element_blank()) +
  theme(strip.text.x = element_blank(),
       strip.text.y = element_blank())

Tip.class

ggsave('Fig5.jpeg', plot = Tip.class,
       device = 'jpeg', path = fig.dir,
       width = 30, height = 20, units = 'cm')

