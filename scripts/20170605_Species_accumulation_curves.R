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
dat.dir <- "~/Documents/PhD/3_EM_Fire_effect/data/"
fig.dir <- '~/Documents/PhD/3_EM_Fire_effect/figures/'
res.dir <- "~/Documents/PhD/3_EM_Fire_effect/results/"

#-----------------------------------------------------------------------------------------
# Load data and clean up----
#-----------------------------------------------------------------------------------------
#--Load data
all.data <- read.csv(paste0(dat.dir,'20170806_OTU_data.csv'), as.is = T)

#--Remove sequences not assigned to an OTU
otu.data <- all.data[!is.na(all.data$otu.97), ]


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

#=========================================================================================
# Species accumulation curves by burn status----
#=========================================================================================

#<< burned sites >> -----------------
#--With singletons
burned.sp <- stsp.matrix[stsp.matrix$Burn_status == 'burned', 5:length(stsp.matrix)]
#--remove columns with no species occurrences
burned.sp <- burned.sp[colSums(burned.sp) > 0]

#--Without singletons
burned.sp.wo <- stsp.matrix[stsp.matrix$Burn_status == 'burned', 5:length(stsp.matrix)]
#--remove columns with single species occurrences
burned.sp.wo <- burned.sp.wo[colSums(burned.sp.wo) > 1]
#--remove rows with no species occurrences
burned.sp.wo <- burned.sp.wo[rowSums(burned.sp.wo) > 0]

#<< unburned sites >> -----------------
unburned.sp <- stsp.matrix[stsp.matrix$Burn_status == 'unburned', 5:length(stsp.matrix)]
#--remove columns with no species occurrences
unburned.sp <- unburned.sp[which(colSums(unburned.sp) > 0)]

#--Without singletons
unburned.sp.wo <- stsp.matrix[stsp.matrix$Burn_status == 'unburned',
                              5:length(stsp.matrix)]
#--remove columns with single species occurrences
unburned.sp.wo <- unburned.sp.wo[colSums(unburned.sp.wo) > 1]
#--remove rows with no species occurrences
burned.sp.wo <- unburned.sp.wo[rowSums(unburned.sp.wo) > 0]

#=========================================================================================
# Species accumulation curves by range----
#=========================================================================================

#<< Pinaleno sites >> -----------------
#--With singletons
pinaleno.sp <- stsp.matrix[stsp.matrix$Range == 'pinaleno', 5:length(stsp.matrix)]
#--remove columns with no species occurrences
pinaleno.sp <- pinaleno.sp[colSums(pinaleno.sp) > 0]

#--Without singletons
pinaleno.sp.wo <- stsp.matrix[stsp.matrix$Range == 'pinaleno',
                              5:length(stsp.matrix)]
#--remove columns with single species occurrences
pinaleno.sp.wo <- pinaleno.sp.wo[colSums(pinaleno.sp.wo) > 1]
#--remove rows with no species occurrences
pinaleno.sp.wo <- pinaleno.sp.wo[rowSums(pinaleno.sp.wo) > 0]

#<< Santa Catalina sites >> -----------------
#--With singletons
sc.sp <- stsp.matrix[stsp.matrix$Range == 'santa.catalina',
                           5:length(stsp.matrix)]
#--remove columns with no species occurrences
sc.sp <- sc.sp[which(colSums(sc.sp) > 0)]

#--Without singletons
sc.sp.wo <- stsp.matrix[stsp.matrix$Range == 'santa.catalina',
                        5:length(stsp.matrix)]
#--remove columns with single species occurrences
sc.sp.wo <- sc.sp.wo[colSums(sc.sp.wo) > 1]
#--remove rows with no species occurrences
sc.sp.wo <- sc.sp.wo[rowSums(sc.sp.wo) > 0]

#=========================================================================================
# Plot----
#=========================================================================================

#--Export as jpeg
png(filename = paste0(fig.dir,"species_accum_curves.jpeg"), width = 1200, height = 1200)
#--Combine into one figure
par(mfrow = c(2,2), "mar"=c(6, 5, 5, 3))

#<< Combine overall burned plots, with singletons and without >>--------------------------
plot(specaccum(burned.sp, sample = min(rowSums(burned.sp),  permutations = 999)),
     xlab = NA, ylab = "OTUs", cex.lab = 2, cex.axis = 2, lwd = 2, yaxt = "n",
     ylim = c(0,100), yaxt = "n")
title(main = 'Burned sites', cex.main = 2)
axis (2, at = seq (0, 100, by = 25), cex.axis = 2, las = 2)
par(new=TRUE)
plot(specaccum(burned.sp.wo, sample = min(rowSums(burned.sp.wo),  permutations = 999)),
     axes = FALSE, xlab = "", ylab = "", lwd = 2, col = "darkblue", ylim = c(0,100))

#<< Combine overall unburned plots, with singletons and without >>------------------------
plot(specaccum(unburned.sp, sample = min(rowSums(unburned.sp),  permutations = 999)),
     xlab = NA, ylab = NA, cex.lab = 2, cex.axis = 2, lwd = 2, 
     ylim = c(0,100), yaxt = "n")
title(main = 'Unburned sites', cex.main = 2)
axis (2, at = seq (0, 100, by = 25), cex.axis = 2, las = 2)
par(new=TRUE)
plot(specaccum(unburned.sp.wo, sample = min(rowSums(unburned.sp.wo),
     permutations = 999)), axes = FALSE, xlab = "", ylab = "", lwd = 2, col = "darkblue",
     ylim = c(0,100))

#<< Combine overall Pinaleno plots, with singletons and without >>------------------------
plot(specaccum(pinaleno.sp, sample = min(rowSums(pinaleno.sp),  permutations = 999)),
     xlab = "Samples", ylab = "OTUs", cex.lab = 2, cex.axis = 2, lwd = 2, 
     ylim = c(0,100), yaxt = "n")
title(main = 'Pinaleno Mts.', cex.main = 2)
axis(2, at = seq(0, 100, by = 25), cex.axis = 2, las = 2)
par(new=TRUE)
plot(specaccum(pinaleno.sp.wo, sample = min(rowSums(pinaleno.sp.wo),
                                              permutations = 999)),
     axes = FALSE, xlab = "", ylab = "", lwd = 2, col = "darkblue", ylim = c(0,100))

#<< Combine overall Santa Catalina plots, with singletons and without >>------------------
plot(specaccum(sc.sp, sample = min(rowSums(sc.sp),  permutations = 999)),
     xlab = "Samples", ylab = NA, cex.lab = 2, cex.axis = 2, lwd = 2,
     ylim = c(0,100), yaxt = "n")
title(main = 'Santa Catalina Mts.', cex.main = 2)
axis(2, at = seq(0, 100, by = 25), cex.axis = 2, las = 2)
par(new=TRUE)
plot(specaccum(sc.sp.wo, sample = min(rowSums(sc.sp.wo),
                                            permutations = 999)),
     axes = FALSE, xlab = "", ylab = "", lwd = 2, col = "darkblue", ylim = c(0,100))

dev.off()
