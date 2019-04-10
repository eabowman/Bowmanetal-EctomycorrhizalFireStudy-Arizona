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

#<< Make data frame with root tip data: class level >> ------------------------
tip_tax <- data.frame(range_burn = c('Santa Catalina burned', 'Santa Catalina unburned',
                          'Pinaleno burned', 'Pinaleno unburned'),
                      range = c('santa.catalina','santa.catalina','pinaleno','pinaleno'),
                      burn = c('burned','unburned','burned','unburned'),
           Agaricomycetes = NA,Dothideomycetes = NA,Pezizomycetes = NA,
           Saccharomycetes = NA,Leotiomycetes = NA,
           Eurotiomycetes = NA,Archaeorhizomycetes = NA)
unique.class <- unique(tax.data$Taxonomy_class)
unique.class <- unique.class[-6]
for(i in unique.class){
  for(c in c('burned', 'unburned')){
    for(d in c('santa.catalina', 'pinaleno')){
      sum.i <- sum(tax.data[tax.data$Range == d & tax.data$Burn_status == c &
                              tax.data$Taxonomy_class == i, 'Tip_count'],
                   na.rm = T)
      tip_tax[tip_tax$range == d & tip_tax$burn == c, i] <- sum.i
    }
  }
}


rownames(tip_tax) <- tip_tax$range_burn
tip.chisq <- tip_tax[4:length(tip_tax)]

#<< Make data frame with root tip data: genus level >> ------------------------
tip.genus.tax <- data.frame(range_burn = c('Santa Catalina burned', 'Santa Catalina unburned',
                                     'Pinaleno burned', 'Pinaleno unburned'),
                      range = c('santa.catalina','santa.catalina','pinaleno','pinaleno'),
                      burn = c('burned','unburned','burned','unburned'),
                      Cenococcum = NA,Russula = NA,Lactarius = NA,Piloderma = NA,
                      Rhizopogon = NA,Sistotrema = NA,Tuber = NA,Sebacina = NA,
                      Tomentella = NA,Inocybe = NA,Tricholoma = NA,Hygrophorus = NA,
                      Hydnobolites = NA,Helvella = NA,Peziza = NA,Phialocephala = NA,
                      Suillus = NA,Amanita = NA,Amphinema = NA,Coltricia = NA,
                      Cantharellus = NA,Oidiodendron = NA,Clavulina = NA,Cortinarius = NA,
                      Elaphomyces = NA,Ramaria = NA, Boletus = NA,
                      Hysterangium = NA, Tulasnella = NA, Laccaria = NA)
unique.gen <- unique(tax.data$Taxonomy_genus)
unique.gen <- unique.gen[-2]
for(i in unique.gen){
  for(c in c('burned', 'unburned')){
    for(d in c('santa.catalina', 'pinaleno')){
      sum.i <- sum(tax.data[tax.data$Range == d & tax.data$Burn_status == c &
                              tax.data$Taxonomy_genus == i, 'Tip_count'],
                   na.rm = T)
      tip.genus.tax[tip.genus.tax$range == d & tip.genus.tax$burn == c, i] <- sum.i
    }
  }
}

rownames(tip.genus.tax) <- tip.genus.tax$range_burn
tip.genus.tax <- tip.genus.tax[4:length(tip.genus.tax)]

#========================================================================================#
# Chi-square: class level
#========================================================================================#
#----------------------------------------------------------------------------------------#
# Range and Fire history
#----------------------------------------------------------------------------------------#
#--Remove rare species (species with less than 3 occurrences)
rare <- c('Saccharomycetes','Archaeorhizomycetes')
tip.chisq <- tip.chisq[ , which(!colnames(tip.chisq) %in% rare)]

#--chi square test
tip.class <- chisq.test(tip.chisq)
tip.class
tax.results[tax.results$tests == 'burn.range.class', 'chi.stat'] <- tip.class$statistic[[1]]
tax.results[tax.results$tests == 'burn.range.class', 'df'] <- tip.class$parameter[[1]]
tax.results[tax.results$tests == 'burn.range.class', 'p.value'] <- tip.class$p.value[[1]]

#----------------------------------------------------------------------------------------#
# Range 
#----------------------------------------------------------------------------------------#
range.class <- data.frame(Agaricomycetes = NA, Dothideomycetes = NA, Pezizomycetes = NA,
           Leotiomycetes = NA, Eurotiomycetes = NA)
for(c in colnames(tip.chisq)){
  sum.c <- sum(tip.chisq[1:2,c])
  range.class['Santa Catalina', c] <- sum.c
  sum.p <- sum(tip.chisq[3:4,c])
  range.class['Pinaleno',c] <- sum.p
}
range.class <- range.class[-1,]

#--chi square test
tip.class <- chisq.test(range.class)
tip.class
tax.results[tax.results$tests == 'range.class', 'chi.stat'] <- tip.class$statistic[[1]]
tax.results[tax.results$tests == 'range.class', 'df'] <- tip.class$parameter[[1]]
tax.results[tax.results$tests == 'range.class', 'p.value'] <- tip.class$p.value[[1]]

#----------------------------------------------------------------------------------------#
# Fire history
#----------------------------------------------------------------------------------------#
fire.class <- data.frame(Agaricomycetes = NA, Dothideomycetes = NA, Pezizomycetes = NA,
                          Leotiomycetes = NA, Eurotiomycetes = NA)
for(c in colnames(tip.chisq)){
  sum.c <- sum(tip.chisq[c(1,3),c])
  fire.class['burned', c] <- sum.c
  sum.p <- sum(tip.chisq[c(2,4),c])
  fire.class['unburned',c] <- sum.p
}
fire.class <- fire.class[-1,]

#--chi square test
tip.class <- chisq.test(fire.class)
tip.class
tax.results[tax.results$tests == 'burn.class', 'chi.stat'] <- tip.class$statistic[[1]]
tax.results[tax.results$tests == 'burn.class', 'df'] <- tip.class$parameter[[1]]
tax.results[tax.results$tests == 'burn.class', 'p.value'] <- tip.class$p.value[[1]]

#========================================================================================#
# Plots
#========================================================================================#

tip.tax <- data.frame(range = rep(c('santa.catalina','pinaleno'), each = 14),
                      burn = rep(c('burned','unburned'),14),
                      taxonomy_class = rep(unique.class,2))
for(r in tip.tax$range){
  for(b in tip.tax$burn){
    for(t in unique.class){
      sum.i <- sum(tax.data[tax.data$Range == r & tax.data$Burn_status == b &
                              tax.data$Taxonomy_class == t, 'Tip_count'],
                   na.rm = T)
      tip.tax[tip.tax$range == r & tip.tax$burn == b &
                tip.tax$taxonomy_class == t, 'count'] <- sum.i
    }
  }
}

tip.tax.ab <- tip.tax[!colnames(tip.tax) %in% rare]

#Change levels of Burn_status and Range columns
tip.tax.ab$burn <- factor(tip.tax.ab$burn, levels = c('burned','unburned'),
                                labels = c('Burned', 'Unburned'))
tip.tax.ab$range <- factor(tip.tax.ab$range, levels = c('pinaleno','santa.catalina'),
                          labels = c('Pinaleno Mts.', 'Santa Catalina Mts.'))

#----------------------------------------------------------------------------------------#
# Range and Fire history
#----------------------------------------------------------------------------------------#
#--Bar graph
Tip.class <- ggplot(data = tip.tax.ab, 
                          aes(x = burn,
                              y = count,
                              fill = taxonomy_class,
                              color = taxonomy_class)) + 
  geom_col(position = 'fill') +
  ylab("Proportion of \n sequences per class") +
  theme_bw() +
  xlab('Fire history') +
  facet_grid(. ~ range) +
  theme(legend.position = 'right',
        axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size=22, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14),
        legend.text=element_text(size=16),
        legend.title=element_blank())

Tip.class

ggsave('TaxonomyClass_RootTip.jpeg', plot = Tip.class,
       device = 'jpeg', path = fig.dir,
       width = 30, height = 20, units = 'cm')

#========================================================================================#
# Chi-square: genus level
#========================================================================================#
#----------------------------------------------------------------------------------------#
# Range and Fire history
#----------------------------------------------------------------------------------------#
rare <- c('Piloderma','Hygrophorus','Coltricia','Cantharellus',
          'Oidiodendron','Clavulina','Cortinarius','Archaerhizomyces','Elaphomyces',
          'Ramaria','Boletus','Hysterangium','Tulasnella','Laccaria','Archaeorhizomyces')
tip.genus.chisq <- tip.genus.tax[ , which(!colnames(tip.genus.tax) %in% rare)]

#--chi square test
tip.genus <- chisq.test(tip.genus.chisq)
tip.genus
tax.results[tax.results$tests == 'burn.range.genus', 'chi.stat'] <- tip.class$statistic[[1]]
tax.results[tax.results$tests == 'burn.range.genus', 'df'] <- tip.class$parameter[[1]]
tax.results[tax.results$tests == 'burn.range.genus', 'p.value'] <- tip.class$p.value[[1]]

#----------------------------------------------------------------------------------------#
# Range 
#----------------------------------------------------------------------------------------#
range.genus <- data.frame(range = c('Santa Catalina', 'Pinaleno'),
                          Cenococcum = NA,Russula = NA,Lactarius = NA,Piloderma = NA,Rhizopogon = NA,
                          Sistotrema = NA,Tuber = NA,Sebacina = NA,Tomentella = NA,Inocybe = NA,
                          Tricholoma = NA,Hygrophorus = NA,Hydnobolites = NA, Helvella = NA, Peziza = NA,
                          Phialocephala = NA,Suillus = NA,Amanita = NA,Amphinema = NA,Coltricia = NA,
                          Cantharellus = NA,Oidiodendron = NA,Clavulina = NA,Cortinarius = NA,
                          Elaphomyces = NA,Ramaria = NA,Boletus = NA,Hysterangium = NA,Tulasnella = NA,
                          Laccaria = NA,Archaeorhizomyces = NA)
for(i in colnames(tip.genus.tax)){
  scm.i <- sum(tip.genus.tax[1:2,i])
  p.i <- sum(tip.genus.tax[3:4,i])
  range.genus[1,i] <- scm.i
  range.genus[2,i] <- p.i
}

row.names(range.genus) <- range.genus[,1]
range.genus <- range.genus[,-1]

rare <- c('Piloderma','Hygrophorus','Coltricia','Cantharellus',
          'Oidiodendron','Clavulina','Cortinarius','Archaerhizomyces','Elaphomyces',
          'Ramaria','Boletus','Hysterangium','Tulasnella','Laccaria','Archaeorhizomyces')
range.genus.chisq <- range.genus[ , which(!colnames(range.genus) %in% rare)]

#--chi square test
tip.genus <- chisq.test(range.genus.chisq)
tip.genus
tax.results[tax.results$tests == 'range.genus', 'chi.stat'] <- tip.class$statistic[[1]]
tax.results[tax.results$tests == 'range.genus', 'df'] <- tip.class$parameter[[1]]
tax.results[tax.results$tests == 'range.genus', 'p.value'] <- tip.class$p.value[[1]]

#----------------------------------------------------------------------------------------#
# Fire history 
#----------------------------------------------------------------------------------------#
burn.genus <- data.frame(range = c('burned', 'unburned'),
                          Cenococcum = NA,Russula = NA,Lactarius = NA,Piloderma = NA,Rhizopogon = NA,
                          Sistotrema = NA,Tuber = NA,Sebacina = NA,Tomentella = NA,Inocybe = NA,
                          Tricholoma = NA,Hygrophorus = NA,Hydnobolites = NA, Helvella = NA, Peziza = NA,
                          Phialocephala = NA,Suillus = NA,Amanita = NA,Amphinema = NA,Coltricia = NA,
                          Cantharellus = NA,Oidiodendron = NA,Clavulina = NA,Cortinarius = NA,
                          Elaphomyces = NA,Ramaria = NA,Boletus = NA,Hysterangium = NA,Tulasnella = NA,
                          Laccaria = NA,Archaeorhizomyces = NA)
for(i in colnames(tip.genus.tax)){
  burn.i <- sum(tip.genus.tax[c(1,3),i])
  unburned.i <- sum(tip.genus.tax[c(2,4),i])
  burn.genus[1,i] <- burn.i
  burn.genus[2,i] <- unburned.i
}

row.names(burn.genus) <- burn.genus[,1]
burn.genus <- burn.genus[,-1]

rare <- c('Piloderma','Hygrophorus','Coltricia','Cantharellus',
          'Oidiodendron','Clavulina','Cortinarius','Archaerhizomyces','Elaphomyces',
          'Ramaria','Boletus','Hysterangium','Tulasnella','Laccaria','Archaeorhizomyces')
burn.genus.chisq <- burn.genus[ , which(!colnames(burn.genus) %in% rare)]

#--chi square test
tip.genus <- chisq.test(burn.genus.chisq)
tip.genus
tax.results[tax.results$tests == 'burn.genus', 'chi.stat'] <- tip.class$statistic[[1]]
tax.results[tax.results$tests == 'burn.genus', 'df'] <- tip.class$parameter[[1]]
tax.results[tax.results$tests == 'burn.genus', 'p.value'] <- tip.class$p.value[[1]]

write.csv(tax.results, './results_output/Taxonomy_roottip_results.csv',
          row.names = F)

#========================================================================================#
# Plots------
#========================================================================================#

tip.genus.tax.plot <- data.frame(range = rep(c('santa.catalina','pinaleno'), each = 34),
                      burn = rep(c('burned','unburned'),34),
                      taxonomy_genus = rep(rep(colnames(burn.genus.chisq), each = 2),2))
for(r in unique(tip.genus.tax.plot$range)){
  for(b in unique(tip.genus.tax.plot$burn)){
    for(t in unique(tip.genus.tax.plot$taxonomy_genus)){
      sum.i <- sum(tax.data[tax.data$Range == r & tax.data$Burn_status == b &
                              tax.data$Taxonomy_genus == t, 'Tip_count'],
                   na.rm = T)
      tip.genus.tax.plot[tip.genus.tax.plot$range == r & tip.genus.tax.plot$burn == b &
                           tip.genus.tax.plot$taxonomy_genus == t, 'count'] <- sum.i
    }
  }
}

#Change levels of Burn_status and Range columns
tip.genus.tax.plot$burn <- factor(tip.genus.tax.plot$burn, levels = c('burned','unburned'),
                                 labels = c('Burned', 'Unburned'))
tip.genus.tax.plot$range <- factor(tip.genus.tax.plot$range, levels = c('pinaleno','santa.catalina'),
                           labels = c('Pinaleno Mts.', 'Santa Catalina Mts.'))

tip.genus.tax.plot <- tip.genus.tax.plot[tip.genus.tax.plot$taxonomy_genus %in%
                     c('Russula','Rhizopogon','Sistotrema','Lactarius','Cenococcum'),]

#----------------------------------------------------------------------------------------#
# Range and Fire history-------
#----------------------------------------------------------------------------------------#
#--Bar graph
Tip.genus <- ggplot(data = tip.genus.tax.plot, 
                    aes(x = burn,
                        y = count,
                        fill = taxonomy_genus,
                        color = taxonomy_genus)) + 
  geom_col(position = 'fill') +
  ylab("Proportion of \n sequences per genus") +
  theme_bw() +
  xlab('Fire history') +
  facet_grid(. ~ range) +
  #ggtitle("Proportion of Classes by Topography") +
  theme(legend.position = 'right',
        axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size=22, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14),
        legend.text=element_text(size=16),
        legend.title = element_blank())

Tip.genus


ggsave('TaxonomyGenus_RootTip.jpeg', plot = Tip.genus,
       device = 'jpeg', path = fig.dir,
       width = 30, height = 20, units = 'cm')

#========================================================================================#
# Heat map: Family------
#========================================================================================#

pdf(file=paste0(fig.dir,'Taxonomyfamily_HeatMap_burn.pdf'),
    paper = "special", width = 6, height = 8, family="Times")

tip.family.tax <- data.frame(range = rep(c('santa.catalina','pinaleno'), each = 52),
                            burn = c(rep('burned',26),rep('unburned',26),
                                     rep('burned',26),rep('unburned',26)),
                            family = rep(c('Amanitaceae','Atheliaceae','Archaeorhizomycetaceae',
                                           'Boletaceae','Cantharellaceae','Gloniaceae',
                                           'Clavulinaceae','Hymenochaetaceae','Cortinariaceae',
                                           'Elaphomycetaceae','Helvellaceae','Pezizaceae',
                                           'Hydnangiaceae','Inocybaceae','Hysterangiaceae',
                                           'Vibrisseaceae','Myxotrichaceae','Russulaceae',
                                           'Gomphaceae',"Rhizopogonaceae","Exidiaceae",
                                          "Hydnaceae","Suillaceae",
                                          "Thelephoraceae","Tuberaceae",
                                          "Tulasnellaceae" ),4),
                            Tip_count = NA,
                            burn_range = NA)
unique.family <- unique(tip.family.tax$family)
for(i in unique.family){
  for(c in c('burned', 'unburned')){
    for(d in c('santa.catalina', 'pinaleno')){
      sum.i <- sum(tax.data[tax.data$Range == d & tax.data$Burn_status == c &
                              tax.data$Taxonomy_family == i, 'Tip_count'],
                   na.rm = T)
      tip.family.tax[tip.family.tax$range == d & tip.family.tax$burn == c &
                       tip.family.tax$family == i, 'Tip_count'] <- sum.i
      burn.range.i <- paste(c,d)
      tip.family.tax[tip.family.tax$range == d & tip.family.tax$burn == c, 'burn_range'] <- burn.range.i
    }
  }
}

rare <- c('Amanitaceae','Archaeorhizomycetaceae','Atheliaceae','Boletaceae','Cantharellaceae',
          'Clavulinaceae','Cortinariaceae','Elaphomycetaceae','Gomphaceae','Hymenochaetaceae',
          'Hysterangiaceae','Myxotrichaceae','Tulasnellaceae')
tip.family.tax <- tip.family.tax[!tip.family.tax$family %in% rare,]

ggplot(tip.family.tax, aes(x = range, y = family)) +
  geom_tile(aes(fill = Tip_count), color = 'white') +
  scale_fill_gradient(low = 'white', high = 'black') +
  xlab('') +
  ylab('') +
  theme_bw() +
  facet_grid(. ~ burn) +
  theme(legend.position="none",
        text = element_text(size = 12),
        axis.text.y = element_text(size = 18, face = 'bold'))
dev.off()  

#========================================================================================#
# Heat map: genus------
#========================================================================================#

pdf(file=paste0(fig.dir,'TaxonomyGenera_HeatMap.pdf'),
    paper = "special", width = 6, height = 8, family="Times")

tip.genus.tax <- data.frame(range = rep(c('santa.catalina','pinaleno'), each = 60),
                            burn = c(rep('burned',30),rep('unburned',30),
                                     rep('burned',30),rep('unburned',30)),
                            genus = rep(c('Cenococcum','Russula','Lactarius','Piloderma','Rhizopogon',
                                      'Sistotrema','Tuber','Sebacina','Tomentella','Inocybe',
                                      'Tricholoma','Hygrophorus','Hydnobolites','Helvella','Peziza',
                                      'Phialocephala','Suillus','Amanita','Amphinema','Coltricia',
                                      'Cantharellus','Oidiodendron','Clavulina','Cortinarius',
                                      'Elaphomyces','Ramaria','Boletus','Hysterangium','Tulasnella',
                                      'Laccaria'),4),
                            Tip_count = NA,
                            burn_range = NA,
                            Rel.ab = NA)
unique.gen <- unique(tip.genus.tax$genus)
for(i in unique.gen){
  for(c in c('burned', 'unburned')){
    for(d in c('santa.catalina', 'pinaleno')){
      sum.i <- sum(tax.data[tax.data$Range == d & tax.data$Burn_status == c &
                              tax.data$Taxonomy_genus == i, 'Tip_count'],
                   na.rm = T)
      tip.genus.tax[tip.genus.tax$range == d & tip.genus.tax$burn == c &
                      tip.genus.tax$genus == i, 'Tip_count'] <- sum.i
      burn.range.i <- paste(c,d)
      tip.genus.tax[tip.genus.tax$range == d & tip.genus.tax$burn == c, 'burn_range'] <- burn.range.i
    }
  }
}

rare <- c('Hysterangium','Coltricia','Clavulina','Piloderma','Peziza',
          'Oidiodendron','Ramaria','Cantharellus','Boletus','Amanita','Amphinema')
tip.genus.tax <- tip.genus.tax[!tip.genus.tax$genus %in% rare,]

total.tips <- sum(tip.genus.tax$Tip_count)
for(i in 1:nrow(tip.genus.tax)){
  rel.abundance.i <- (tip.genus.tax[i, 'Tip_count'])/total.tips
  tip.genus.tax[i,'Rel.ab'] <- rel.abundance.i
}

ggplot(tip.genus.tax, aes(x = burn_range, y = genus)) +
  geom_tile(aes(fill = Rel.ab), color = 'white') +
  scale_fill_gradient(low = 'white', high = 'darkgreen') +
  xlab('') +
  ylab('') +
  theme_bw() +
  theme(legend.position="none",
        text = element_text(size = 12),
        axis.text.y = element_text(size = 14))
dev.off()  
