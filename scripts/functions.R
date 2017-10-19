## Script created by Liz Bowman July 5, 2017
## functions for analyzing data

library(stats)
#=========================================================================================
# anova.t: function to test differences between samples----
#=========================================================================================

anova.t <- function(data, group, response) {
  anova.output <- data.frame(response = response,
                             F.value =  NA,
                             df.1 = NA,
                             df.2 = NA,
                             p = NA)
  for(r in response){
    simple <- anova(lm(data[,r] ~ data[,group]))
    anova.output[anova.output$response == r, 'F.value'] <- simple$`F value`[1]
    anova.output[anova.output$response == r, 'df.1'] <- simple$Df[1]
    anova.output[anova.output$response == r, 'df.2'] <- simple$Df[2]
    anova.output[anova.output$response == r, 'p'] <- simple$`Pr(>F)`[1]
    # cat('Predictor: ', predictor, '\n')
  }
  return(anova.output)
}

#=========================================================================================
# Ttester: function to test differences for differences between samples----
#=========================================================================================

ttester <- function(data, group, response) {
  ttest.output <- data.frame(response = response,
                             t.stat =  NA,
                             p = NA)
  for(r in response){
    simple <- t.test(data[,r] ~ data[,group])
    ttest.output[ttest.output$response == r, 't.stat'] <- simple$statistic
    ttest.output[ttest.output$response == r, 'p'] <- simple$p.value
    # cat('Predictor: ', predictor, '\n')
  }
  return(ttest.output)
}

wilcoxtest <- ttester <- function(data, group, response) {
  wilcox.output <- data.frame(response = response,
                             t.stat =  NA,
                             p = NA)
  for(r in response){
    simple <- wilcox.test(data[,r] ~ data[,group])
    wilcox.output[wilcox.output$response == r, 't.stat'] <- simple$statistic
    wilcox.output[wilcox.output$response == r, 'p'] <- simple$p.value
    # cat('Predictor: ', predictor, '\n')
  }
  return(wilcox.output)
}


# < plot > -----------------------------
# boxplot.liz <- function(data, x.factor, y.factor, fill.factor,facet.factor) {
#   responses <- NA
#   for(response in y.factor) {
#     plot.r <- paste0(response,'.plot')
#     plot.r <- ggplot(data, aes(x = x.factor, y = response, fill = fill.factor)) +
#       geom_boxplot() +
#       facet_grid(. ~ facet.factor) +
#       theme_bw() +
#       scale_fill_brewer(palette = "Accent") +
#       labs(fill = fill.factor)
#     return(plot.r)
#     responses <- c(responses,paste0(response,'.plot'))
#   }
#   require(gridExtra)
#   grid.arrange(responses,nrow = 5, ncol = 4)
# }

#=========================================================================================
# Normtester: test for normality of variables----
#=========================================================================================

normtest <- function(data, factor){
  par(mfrow=c(1,2))
  ## Have a look at the densities
  plot(density(data[,factor]))
  ## Plot using a qqplot
  qqnorm(data[,factor]); qqline(data[,factor], col = 2) # normality
  test <- shapiro.test(data[,factor]) # normality
  print(factor)
  print(test)
}

## plot residuals
# lm = lm(y ~ x, data=data.set) 
# res = resid(lm)
# plot(x, res) 
# abline(0, 0) 

#=========================================================================================
# Multiplot: to plot multiple graphs in one figure----
#=========================================================================================

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
