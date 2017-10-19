otu.97 <- read.csv('~/Documents/PhD/3_EM_Fire_effect/OTU_cluster_files/97%_otu.csv',
                   as.is = T, col.names = T, row.names = F)
otu.97
otu.95 <- read.csv('~/Documents/PhD/3_EM_Fire_effect/OTU_cluster_files/95%_otu.csv',
                   as.is = T)
otu.95

for(i in colnames(otu.97)){
  sample <- as.matrix(otu.97[i])
  for (s in sample) {
    tax.data[tax.data$Sample_name == s, 'OTU.97'] <- colnames(sample)
  }
}

for(i in colnames(otu.95)){
  sample <- as.matrix(otu.95[i])
  for (s in sample) {
    tax.data[tax.data$Sample_name == s, 'OTU.95'] <- colnames(sample)
  }
}
