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
fig.dir <- '~/Documents/PhD/2_EM_Fire_effect/figures/'
res.dir <- "~/Documents/PhD/2_EM_Fire_effect/results/"

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

#--<< Combine range and burn status >> ----------
for(r in unique(stsp.matrix$Range)) {
  for(b in unique(stsp.matrix$Burn_status)) {
    stsp.matrix[stsp.matrix$Burn_status == b & stsp.matrix$Range == r, 'burn_range'] <-
      paste(b, r)
  }
}

stsp.matrix <- stsp.matrix[c(1:4,128,5:127)]

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

#--Export as jpeg
png(filename = paste0(fig.dir,"species_accum_curves.jpeg"), width = 1200, height = 1200)
#--Combine into one figure
par(mfrow = c(2,2), "mar"=c(6, 5, 5, 3))

#<< Combine Pinaleno burned plots, with singletons and without >>--------------------------
plot(specaccum(PMburned.sp, sample = min(rowSums(PMburned.sp),  permutations = 999)),
     xlab = NA, ylab = "OTUs", cex.lab = 2, cex.axis = 2, lwd = 2, yaxt = "n",
     ylim = c(0,100), yaxt = "n")
title(main = 'Pinaleno burned sites', cex.main = 2)
axis (2, at = seq (0, 100, by = 25), cex.axis = 2, las = 2)
par(new=TRUE)
plot(specaccum(PMburned.sp.wo, sample = min(rowSums(PMburned.sp.wo),  permutations = 999)),
     axes = FALSE, xlab = "", ylab = "", lwd = 2, col = "darkblue")

#<< Combine Pinaleno unburned plots, with singletons and without >>------------------------
plot(specaccum(PMunburned.sp, sample = min(rowSums(PMunburned.sp),  permutations = 999)),
     xlab = NA, ylab = NA, cex.lab = 2, cex.axis = 2, lwd = 2, 
     ylim = c(0,100), yaxt = "n")
title(main = 'Pinaleno unburned sites', cex.main = 2)
axis (2, at = seq (0, 100, by = 25), cex.axis = 2, las = 2)
par(new=TRUE)
plot(specaccum(PMunburned.sp.wo, sample = min(rowSums(PMunburned.sp.wo),
     permutations = 999)), axes = FALSE, xlab = "", ylab = "", lwd = 2, col = "darkblue",
     ylim = c(0,100))

#<< Combine Santa Catalina burned plots, with singletons and without >>---------------------
plot(specaccum(SCMburned.sp, sample = min(rowSums(SCMburned.sp),  permutations = 999)),
     xlab = NA, ylab = "OTUs", cex.lab = 2, cex.axis = 2, lwd = 2, yaxt = "n",
     ylim = c(0,100), yaxt = "n")
title(main = 'Santa Catalina burned sites', cex.main = 2)
axis (2, at = seq (0, 100, by = 25), cex.axis = 2, las = 2)
par(new=TRUE)
plot(specaccum(SCMburned.sp.wo, sample = min(rowSums(SCMburned.sp.wo),  permutations = 999)),
     axes = FALSE, xlab = "", ylab = "", lwd = 2, col = "darkblue", ylim = c(0,100))

#<< Combine Santa Catalina unburned plots, with singletons and without >>-------------------
plot(specaccum(SCMunburned.sp, sample = min(rowSums(SCMunburned.sp),  permutations = 999)),
     xlab = NA, ylab = NA, cex.lab = 2, cex.axis = 2, lwd = 2, 
     ylim = c(0,100), yaxt = "n")
title(main = 'Santa Catalina unburned sites', cex.main = 2)
#axis (2, at = seq (0, 100, by = 25), cex.axis = 2, las = 2)
par(new=TRUE)
plot(specaccum(SCMunburned.sp.wo, sample = min(rowSums(SCMunburned.sp.wo),
                                              permutations = 999)), axes = FALSE, xlab = "", ylab = "", lwd = 2, col = "darkblue",
     ylim = c(0,100))

dev.off()
