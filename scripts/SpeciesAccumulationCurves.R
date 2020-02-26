## Script created by Liz Bowman February 13, 2019
## for calculating species accumulation curves

#========================================================================================#
# Load data----
#========================================================================================#

otu.data <- read.csv(paste0(dat.dir,'97%_SitexSpecies_TipAb.csv'),
                     as.is = T)
rownames(otu.data) <- otu.data$Tree

#========================================================================================#
# Species accumulation curves: all data----
#========================================================================================#
#--isolate community data
comm.matrix <- otu.data[12:length(otu.data)]

#--create data frame
all.sa <- data.frame(community = rep(c('Overall','Pinaleno Burned','Pinaleno Unburned',
                             'Santa Catalina Burned','Santa Catalina Unburned'),
                           c(41,10,10,10,11)),
           Sites = NA,
           Richness = NA,
           sd = NA)

#--Overall
overall.sa <- specaccum(comm.matrix)
all.sa[all.sa$community == 'Overall', 'Sites'] <- overall.sa$sites
all.sa[all.sa$community == 'Overall', 'Richness'] <- overall.sa$richness
all.sa[all.sa$community == 'Overall', 'sd'] <- overall.sa$sd

#--Pinaleno burned
# isolate pinaleno burned data
p.b <- rownames(otu.data[otu.data$Range == 'pinaleno' & otu.data$Burn_status == 'burned',])
p.b.data <- comm.matrix[rownames(comm.matrix) %in% p.b,]
pb.sa <- specaccum(p.b.data)
all.sa[all.sa$community == 'Pinaleno Burned', 'Sites'] <- pb.sa$sites
all.sa[all.sa$community == 'Pinaleno Burned', 'Richness'] <- pb.sa$richness
all.sa[all.sa$community == 'Pinaleno Burned', 'sd'] <- pb.sa$sd

#--Pinaleno unburned
# isolate pinaleno unburned data
p.u <- rownames(otu.data[otu.data$Range == 'pinaleno' & otu.data$Burn_status == 'unburned',])
p.u.data <- comm.matrix[rownames(comm.matrix) %in% p.u,]
pu.sa <- specaccum(p.u.data)
all.sa[all.sa$community == 'Pinaleno Unburned', 'Sites'] <- pu.sa$sites
all.sa[all.sa$community == 'Pinaleno Unburned', 'Richness'] <- pu.sa$richness
all.sa[all.sa$community == 'Pinaleno Unburned', 'sd'] <- pu.sa$sd

#--Santa Catalina burned
# isolate Santa Catalina burned data
scm.b <- rownames(otu.data[otu.data$Range == 'santa.catalina' & otu.data$Burn_status == 'burned',])
scm.b.data <- comm.matrix[rownames(comm.matrix) %in% scm.b,]
scmb.sa <- specaccum(scm.b.data)
all.sa[all.sa$community == 'Santa Catalin Burned', 'Sites'] <- scmb.sa$sites
all.sa[all.sa$community == 'Santa Catalina Burned', 'Richness'] <- scmb.sa$richness
all.sa[all.sa$community == 'Santa Catalina Burned', 'sd'] <- scmb.sa$sd

#--Santa Catalina unburned
# isolate pinaleno unburned data
scm.u <- rownames(otu.data[otu.data$Range == 'santa.catalina' & otu.data$Burn_status == 'unburned',])
scm.u.data <- comm.matrix[rownames(comm.matrix) %in% scm.u,]
scmu.sa <- specaccum(scm.u.data)
all.sa[all.sa$community == 'Santa Catalina Unburned', 'Sites'] <- scmu.sa$sites
all.sa[all.sa$community == 'Santa Catalina Unburned', 'Richness'] <- scmu.sa$richness
all.sa[all.sa$community == 'Santa Catalina Unburned', 'sd'] <- scmu.sa$sd

#--plot
ggplot(all.sa, aes(x = Sites,
                   y = Richness,
                   color = community)) +
  geom_smooth(se = F) +
  theme_classic()

jpeg(filename = 'figures_output/Fig2A.jpeg',
     width = 700, height = 600,
     quality = 100)

plot(overall.sa, cex = 4, lwd = 4, xlab = 'Trees')
plot(pb.sa, add = T, col = 2, cex = 4, lwd = 4)
plot(pu.sa, add = T, col = 'darkgreen', cex = 4, lwd = 4)
plot(scmb.sa, add = T, col = 4, cex = 4, lwd = 4)
plot(scmu.sa, add = T, col = 'darkgrey', cex = 4, lwd = 4)
dev.off()

jpeg(filename = 'figures_output/Fig2A_insert.jpeg',
     width = 700, height = 600,
     quality = 100)

plot(pb.sa, col = 2, cex = 4, lwd = 4)
plot(pu.sa, add = T, col = 'darkgreen', cex = 4, lwd = 4)
plot(scmb.sa, add = T, col = 4, cex = 4, lwd = 4)
plot(scmu.sa, add = T, col = 'darkgrey', cex = 4, lwd = 4)
dev.off()

#========================================================================================#
# Species accumulation curves: > 4 occurences----
#========================================================================================#
#--isolate otu with greater than 4 occurrences
comm.matrix.ab <- comm.matrix[colSums(comm.matrix) > 4]

#--create data frame
abund.sa <- data.frame(community = rep(c('Overall','Pinaleno Burned','Pinaleno Unburned',
                                       'Santa Catalina Burned','Santa Catalina Unburned'),
                                     c(41,10,10,10,11)),
                     Sites = NA,
                     Richness = NA,
                     sd = NA)

#--Overall
overall.sa <- specaccum(comm.matrix.ab)
abund.sa[abund.sa$community == 'Overall', 'Sites'] <- overall.sa$sites
abund.sa[abund.sa$community == 'Overall', 'Richness'] <- overall.sa$richness
abund.sa[abund.sa$community == 'Overall', 'sd'] <- overall.sa$sd

#--Pinaleno burned
# isolate pinaleno burned data
p.b <- rownames(otu.data[otu.data$Range == 'pinaleno' & otu.data$Burn_status == 'burned',])
p.b.data <- comm.matrix.ab[rownames(comm.matrix.ab) %in% p.b,]
pb.sa <- specaccum(p.b.data)
abund.sa[abund.sa$community == 'Pinaleno Burned', 'Sites'] <- pb.sa$sites
abund.sa[abund.sa$community == 'Pinaleno Burned', 'Richness'] <- pb.sa$richness
abund.sa[abund.sa$community == 'Pinaleno Burned', 'sd'] <- pb.sa$sd

#--Pinaleno unburned
# isolate pinaleno unburned data
p.u <- rownames(otu.data[otu.data$Range == 'pinaleno' & otu.data$Burn_status == 'unburned',])
p.u.data <- comm.matrix.ab[rownames(comm.matrix.ab) %in% p.u,]
pu.sa <- specaccum(p.u.data)
abund.sa[abund.sa$community == 'Pinaleno Unburned', 'Sites'] <- pu.sa$sites
abund.sa[abund.sa$community == 'Pinaleno Unburned', 'Richness'] <- pu.sa$richness
abund.sa[abund.sa$community == 'Pinaleno Unburned', 'sd'] <- pu.sa$sd

#--Santa Catalina burned
# isolate Santa Catalina burned data
scm.b <- rownames(otu.data[otu.data$Range == 'santa.catalina' & otu.data$Burn_status == 'burned',])
scm.b.data <- comm.matrix.ab[rownames(comm.matrix.ab) %in% scm.b,]
scmb.sa <- specaccum(scm.b.data)
abund.sa[abund.sa$community == 'Santa Catalin Burned', 'Sites'] <- scmb.sa$sites
abund.sa[abund.sa$community == 'Santa Catalina Burned', 'Richness'] <- scmb.sa$richness
abund.sa[abund.sa$community == 'Santa Catalina Burned', 'sd'] <- scmb.sa$sd

#--Santa Catalina unburned
# isolate pinaleno unburned data
scm.u <- rownames(otu.data[otu.data$Range == 'santa.catalina' & otu.data$Burn_status == 'unburned',])
scm.u.data <- comm.matrix.ab[rownames(comm.matrix.ab) %in% scm.u,]
scmu.sa <- specaccum(scm.u.data)
abund.sa[abund.sa$community == 'Santa Catalina Unburned', 'Sites'] <- scmu.sa$sites
abund.sa[abund.sa$community == 'Santa Catalina Unburned', 'Richness'] <- scmu.sa$richness
abund.sa[abund.sa$community == 'Santa Catalina Unburned', 'sd'] <- scmu.sa$sd

#--plot
ggplot(abund.sa, aes(x = Sites,
                   y = Richness,
                   color = community)) +
  geom_smooth() +
  theme_classic()

jpeg(filename = 'figures_output/Fig2B.jpeg',
     width = 700, height = 600,
     quality = 100)

plot(overall.sa, cex = 4, lwd = 4, xlab = 'Trees')
plot(pb.sa, add = T, col = 2, cex = 4, lwd = 4)
plot(pu.sa, add = T, col = 'darkgreen', cex = 4, lwd = 4)
plot(scmb.sa, add = T, col = 4, cex = 4, lwd = 4)
plot(scmu.sa, add = T, col = 'darkgrey', cex = 4, lwd = 4)
dev.off()

jpeg(filename = 'figures_output/Fig2B_insert.jpeg',
     width = 700, height = 600,
     quality = 100)

plot(pb.sa, col = 2, cex = 4, lwd = 4)
plot(pu.sa, add = T, col = 'darkgreen', cex = 4, lwd = 4)
plot(scmb.sa, add = T, col = 4, cex = 4, lwd = 4)
plot(scmu.sa, add = T, col = 'darkgrey', cex = 4, lwd = 4)
dev.off()
