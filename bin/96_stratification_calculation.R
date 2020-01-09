### Stratification indices
require(oce)  # package for handling CTD data

# require(tidyverse)
# take difference near surf and near bottom and divide by water column depth (z) to give density per meter gradient with CTD data
# Simpson and Hunter calculation from Mann & Lazier, Dynamics of Marine Ecosystems: Biological-Physical Interactions in the Ocean

stratify <- function(y) {
  (y$rho[match(max(y$z),y$z)] - y$rho[match(min(y$z),y$z)])/y$bottom[1]
}

# ASGARD
asgard <- read.csv('data/ASGARD2017_SKQ201709S_NiskinBottle_CTD_CHL_NUT.csv')
asgard$rho <- swRho(asgard$Sal00,asgard$T090.degC.,asgard$Pressure.dbar.)

stdf <- data.frame("station"=asgard$X.Station,
                "rho"=asgard$rho,"z"=asgard$Pressure.dbar.,
                "bottom"=asgard$BottomDepth.m.)
stdf = unique(stdf)

asgard.strat <- unclass(by(stdf, stdf$station, stratify))

asgard.strat = data.frame(asgard[match(as.numeric(names(asgard.strat)),asgard$X.Station),], asgard.strat)

ggplot(data=asgard.strat, aes(x=Longitude.W., y=Latitude.S.)) + 
  geom_point(aes(size=asgard.strat, color=asgard.strat)) +
  scale_color_viridis_c()



# AMBON
ambon <- read.csv('AMBON2017_CTD.csv')
ambon$rho <- swRho(ctd$sal00..Salinity..Practical..PSU., ctd$t090C..Temperature..ITS.90..deg.C.., ctd$prDM..Pressure..Digiquartz..db.)

stdf <- data.frame("station"=ambon$StationName,
                   "rho"=ambon$rho,"z"=ambon$prDM..Pressure..Digiquartz..db.,
                   "bottom"=ambon$Bottom.depth..m.)
stdf = unique(stdf)

asgard.strat <- unclass(by(stdf, stdf$station, stratify))
