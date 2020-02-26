## Script created by Liz Bowman October 20, 2016
## for testing normality of data

#=========================================================================================
# Load packages and data 
#=========================================================================================
script.dir <- '~/Documents/PhD/2_Sky_islands/scripts/'
source(paste0(script.dir, '20170111_Load_data.R'))

#=========================================================================================
# Histogram to assess distribution-------
#=========================================================================================
#--EM diversity:fa
ggplot(div.data.fa, aes(x = fa)) +
  geom_histogram(aes(y = ..density..), color = 'black', fill = 'white') +
  labs(x = 'EM abundance based on dry root weight', y = 'Density') +
  stat_function(fun = dnorm, args = list(mean = mean(div.data.fa$fa)),
                                         sd = sd((div.data.fa$fa)),
                color = 'black', size = 1)

ggplot(div.data.fa, aes(x = log(fa), colour = range)) +
    geom_density()

#--EM diversity: shannon's index
ggplot(div.data.fa, aes(x = shannon)) +
  geom_histogram(aes(y = ..density..), color = 'black', fill = 'white') +
  labs(x = 'EM abundance based on dry root weight', y = 'Density') +
  stat_function(fun = dnorm, args = list(mean = mean(div.data.fa$shannon)),
                sd = sd((div.data.fa$shannon)),
                color = 'black', size = 1)

ggplot(div.data.fa, aes(x = shannon, colour = range)) +
  geom_density()

#=========================================================================================
# Normal distribution-------
#=========================================================================================
#<< EM Diversity: FA >>---------------------------------
#--residuals
plot(aov(log(fa) ~ range, data = div.data.fa),1)
#--qqplot
plot(aov(log(fa) ~ range, data = div.data.fa),2)
shapiro.test(log(div.data.fa$fa))

#<< EM Diversity: Shannon's index >>---------------------------------
#--residuals
plot(aov(shannon ~ range, data = div.data.fa),1)
#--qqplot
plot(aov(shannon ~ range, data = div.data.fa),2)
shapiro.test(div.data.fa$shannon)

#=========================================================================================
# Homogeneity of variance-------
#=========================================================================================

#-----------------------------------------------------------------------------------------
# Ectomycorrhizal diversity
#-----------------------------------------------------------------------------------------
#--FA
leveneTest(div.data.fa$fa, div.data.fa$range, center = median)

#--Shannon's index
leveneTest(div.data.fa$shannon, div.data.fa$range, center = median)


