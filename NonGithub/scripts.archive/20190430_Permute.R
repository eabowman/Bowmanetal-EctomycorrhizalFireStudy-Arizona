## Script created by Liz Bowman October 29, 2019
## for analyzing community variation based on
## abiotic and biotic factors

#========================================================================================#
# Load data and libraries----
#========================================================================================#

#----------------------------------------------------------------------------------------#
# Load libraries----
#----------------------------------------------------------------------------------------#
library(ggplot2);library (tidyr);library (vegan);library (dplyr)

#----------------------------------------------------------------------------------------#
# set up paths to directories----
#----------------------------------------------------------------------------------------#
#--path to directory where the climate data downloaded from BIOCLIM site is stored
dat.dir <- "~/Documents/PhD/2_EM_Fire_effect/data/"
fig.dir <- '~/Documents/PhD/2_EM_Fire_effect/figures_output/'
res.dir <- "~/Documents/PhD/2_EM_Fire_effect/results_output/"

#----------------------------------------------------------------------------------------#
# Load data and clean up: Tree level----
#----------------------------------------------------------------------------------------#
stsp.matrix <- read.csv('data_output/97%_SitexSpecies_TipAb.csv', as.is = T)

#--Calculate PCA of soil factors
pca.env <- stsp.matrix[c('pH.su','no3.n.ppm','po4.p.ppm','so4.s.ppm')]
rownames(pca.env) <- stsp.matrix$Tree

pca.env.analysis <- prcomp(pca.env,
                           center = TRUE,
                           scale. = TRUE)
summary(pca.env.analysis)
env.eigenvector <- scores(pca.env.analysis, choices = c(1:4))
env.eigenvector <- data.frame(matrix(unlist(env.eigenvector), nrow=41, byrow=T))
names(env.eigenvector) <- c('PCA1','PCA2','PCA3')
env.eigenvector$tree <- rownames(pca.env)

stsp.matrix['env.pca'] <- NA
for(i in unique(stsp.matrix$Tree)){
  stsp.matrix[stsp.matrix$Tree == i, 'env.pca'] <-
    env.eigenvector[env.eigenvector$tree == i, 'PCA1']
}

# << Create distinct tables for each burn history/range >> --
FA.matrix <- stsp.matrix[stsp.matrix$Burn_status == 'burned',]
FU.matrix <- stsp.matrix[stsp.matrix$Burn_status == 'unburned',]

scm.matrix <- stsp.matrix[stsp.matrix$Range == 'santa.catalina',]
pm.matrix <- stsp.matrix[stsp.matrix$Range == 'pinaleno',]

#========================================================================================#
# Jaccard based dissimilarity index:----
#========================================================================================#

#--isolate otu data
comm.matrix <- stsp.matrix[c(12:127)]

#--comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 2]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 1, ] # remove rows with sums of 0
jaccard.matrix <- stsp.matrix[row.names(comm.matrix),]

#--distance matrix using jaccard index
comm.dist.jaccard <- vegdist(comm.matrix, method = "jaccard", binary = F, na.rm = T)

#--NMDS analysis
jaccard.otu <- metaMDS(comm.dist.jaccard, dist = "bray", permutations = 999,
                       try = 100, trymax = 1000)

#--Permanova
adonis(comm.dist.jaccard ~ Burn_status * prec * forest * env.pca,
       data = jaccard.matrix, permutations = 1000)

#========================================================================================#
# Morisita based dissimilarity index:----
#========================================================================================#

#--isolate otu data
comm.matrix <- stsp.matrix[c(12:127)]

#--comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 2]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 1, ] # remove rows with sums of 0
morisita.matrix <- stsp.matrix[row.names(comm.matrix),]

#--distance matrix using jaccard index
comm.dist.morisita <- vegdist(comm.matrix, method = "morisita", binary = F, na.rm = T)

#--Permanova
adonis(comm.dist.morisita ~ Burn_status * prec * forest * env.pca,
       data = morisita.matrix, permutations = 1000)

#========================================================================================#
# Jaccard based dissimilarity index: Pinaleno----
#========================================================================================#

#--isolate otu data
comm.matrix <- pm.matrix[c(12:127)]

#--comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 2]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 1, ] # remove rows with sums of 0
jaccard.matrix <- pm.matrix[row.names(comm.matrix),]

#--distance matrix using jaccard index
comm.dist.jaccard <- vegdist(comm.matrix, method = "jaccard", binary = F, na.rm = T)

#--NMDS analysis
jaccard.otu <- metaMDS(comm.dist.jaccard, dist = "bray", permutations = 999,
                       try = 100, trymax = 1000)

#--Permanova
adonis(comm.dist.jaccard ~ Burn_status * prec * forest * env.pca,
       data = jaccard.matrix, permutations = 1000)

#========================================================================================#
# Morisita based dissimilarity index: Pinaleno----
#========================================================================================#

#--isolate otu data
comm.matrix <- pm.matrix[c(12:127)]

#--comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 2]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 1, ] # remove rows with sums of 0
morisita.matrix <- pm.matrix[row.names(comm.matrix),]

#--distance matrix using jaccard index
comm.dist.morisita <- vegdist(comm.matrix, method = "morisita", binary = F, na.rm = T)

#--Permanova
pm.morisita <-adonis(comm.dist.morisita ~ Burn_status * prec * env.pca,
       data = morisita.matrix, permutations = 1000)

#========================================================================================#
# Jaccard based dissimilarity index: SCM----
#========================================================================================#

#--isolate otu data
comm.matrix <- scm.matrix[c(12:127)]

#--comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 2]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 1, ] # remove rows with sums of 0
jaccard.matrix <- scm.matrix[row.names(comm.matrix),]

#--distance matrix using jaccard index
comm.dist.jaccard <- vegdist(comm.matrix, method = "jaccard", binary = F, na.rm = T)

#--NMDS analysis
jaccard.otu <- metaMDS(comm.dist.jaccard, dist = "bray", permutations = 999,
                       try = 100, trymax = 1000)

#--Permanova
adonis(comm.dist.jaccard ~ Burn_status * prec * env.pca,
       data = jaccard.matrix, permutations = 1000)

#========================================================================================#
# Morisita based dissimilarity index: SCM----
#========================================================================================#

#--isolate otu data
comm.matrix <- scm.matrix[c(12:127)]

#--comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 2]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 1, ] # remove rows with sums of 0
morisita.matrix <- scm.matrix[row.names(comm.matrix),]

#--distance matrix using jaccard index
comm.dist.morisita <- vegdist(comm.matrix, method = "morisita", binary = F, na.rm = T)

#--Permanova
scm.morisita <- adonis(comm.dist.morisita ~ Burn_status * prec * env.pca,
       data = morisita.matrix, permutations = 1000)
