## Script adapted from script created by Yu-Ling Huang
## modified by Liz Bowman (eabowman on github) Oct. 1, 2019

#========================================================================================#
# Load data ----
#========================================================================================#

sites <- read.csv(paste0(dat.dir,'site_data.csv'), as.is = T)

# isolate site coordinates
sites %>%
  dplyr::select(site, range, burn_status, lat, long) %>%
  distinct(site, .keep_all = T) -> sites
sites$range <- factor(sites$range)

#========================================================================================#
# Download map data ----
#========================================================================================#

# get USA map from GADM, level=0 no state border,
# level=1 state border, level=2 county border...
#US_0<-getData('GADM', country="USA", level=0)
US_1<-getData('GADM', country="USA", level=1)
#US_2<-getData('GADM', country="USA", level=2)
#US_3<-getData('GADM', country="USA", level=3)
#US_4<-getData('GADM', country="USA", level=4)

# get the Arizona state border map from the whole USA map
Arizona<-subset(US_1, NAME_1=="Arizona")
plot(Arizona)

# << make elevation map with sample site points >> ------------
pdf(paste0(fig.dir,'Fig1.pdf'), height = 8, width = 10)

Alt<-getData("worldclim", var="alt", res=0.5, lon=-110, lat=33)
Alt.sub<-crop(Alt, extent(Arizona))
Alt.sub.mask<-mask(Alt.sub, Arizona)
Alt.col<-colorRampPalette(c("white", "gray90","gray80","gray70","gray60","gray50","gray30","black"))

# point color indicate geographic area
color.vec<-c("red", "yellow","blue","green")
plot(Alt.sub.mask$alt_12,legend=FALSE, axes=FALSE, box=FALSE, col=Alt.col(8),
     xlab= "", ylab="")
plot(Arizona,add=TRUE)
axis(1, pos=31, cex.axis=2, tck=-0.02)
axis(2, las=2, cex.axis=2, tck=-0.02)
# mtext(text = 'Longitude', side = 1, line = 4)
# mtext(text = 'Latitude', side = 2, line = 3)
sites.sub <- sites[c(1,5),]
points(sites.sub$long, sites.sub$lat, pch=19, cex=2)

# Legend
par(xpd=TRUE, mar = c(5,4,4,0))
Alt.sub.mask.range<-c(500,maxValue(Alt.sub.mask))
plot(Alt.sub.mask$alt_12, legend.only=TRUE,col=Alt.col(8),
     legend.width=2, legend.shrink=1,
     axis.args=list(at=seq(Alt.sub.mask.range[1], Alt.sub.mask.range[2], 1000),
                    labels=seq(Alt.sub.mask.range[1], Alt.sub.mask.range[2],1000), 
                    cex.axis=1.25),
     legend.args=list(text='Elevation (m)', side=4, font=2, line=3.5, cex=1.25))

dev.off()
