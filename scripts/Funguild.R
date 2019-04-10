funguild$taxonomy <- NA
for(i in 1:nrow(funguild)){
  otu.i <- funguild[i,'otu']
  p.i <- tax[tax$otu == otu.i, 2]
  c.i <- tax[tax$otu == otu.i, 3]
  o.i <- tax[tax$otu == otu.i, 4]
  g.i <- tax[tax$otu == otu.i, 5]
  s.i <- tax[tax$otu == otu.i, 6]
  tax.i <- paste0(p.i,';',c.i,';',o.i,';',g.i,';',s.i)
  funguild[funguild$otu == otu.i, 'taxonomy'] <- tax.i
}

colnames(funguild)[1] <- 'OTU_ID'

write.csv(funguild, 'FUNGuild/OTU97_withTax.csv')
