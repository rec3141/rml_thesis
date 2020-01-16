## map clusters to sampling sites

#require(plyr)
# require(dplyr)
# require(tidyverse)
# require(marmap)
# require(maps)
# require(mapdata)
# require(reshape2)
# require(RColorBrewer)
# require(mapproj)
require(PBSmapping)
require(marmap)
require(scales)

# Bering and Chukchi seas
xmin = 185
xmax = 205 #145
ymin = 62
ymax = 74

# download bathymetry data
noaadat = getNOAA.bathy(xmin-360,xmax-360,ymin,ymax,res=1, keep=TRUE)
noaamat = fortify(noaadat)   # converts b class 'bathy' to data frame
# saveRDS(noaamat,"data/NOAA_Chukchi_bathymetry.RDS")
# noaamat = readRDS("data/NOAA_Chukchi_bathymetry.RDS")

noaadat[noaadat>0] = NA #remove terrestrial terrain


#Shows all available map data
myworld = map_data('world2')
wholemap = ggplot() + geom_polygon(data=myworld,aes(x=long,y=lat,group=factor(group)))

#Zooms in and clips the map
names(myworld) = c("X","Y","PID","POS","region","subregion")
myworld = clipPolys(myworld,xlim=c(xmin,xmax),ylim=c(ymin,ymax),keepExtra=TRUE)

#This little bit of code changes the longitude from positive to negative
#You can't change it earlier, R needs to pull from 'world2' and then you
#can change the sign
myworld$X = myworld$X-360
# head(myworld)
xmin = -175
xmax = -155



#default map
basemap = ggplot() +
  geom_polygon(data=myworld,aes(x=X,y=Y,group=factor(PID))) +
  coord_map(xlim=c(xmin,xmax),ylim = c(ymin,ymax)) +
  theme_bw()

# Plot bathy object using custom ggplot2 functions
bbmap = autoplot(noaadat, geom=c("r"),coast=FALSE) + 
  scale_fill_gradient2(name="Depth",low=muted("blue"),high = muted("blue",l=10), limits=c(-100,0),oob=squish) +
  theme_bw()

#Fix axis labels
basemap = basemap + 
  labs( x = expression("Longitude ("*degree*"W)"), y = expression("Latitude ("*degree*"N)")) +
  scale_x_continuous(breaks=c(-175,-170,-165,-160,-155),labels=c("175","170","165","160","155"))

bbmap = bbmap + 
  labs( x = expression("Longitude ("*degree*"W)"), y = expression("Latitude ("*degree*"N)")) +
  scale_x_continuous(breaks=c(-175,-170,-165,-160,-155),labels=c("175","170","165","160","155"))

# add bathymetry contours to basemap
batmap = basemap + 
  geom_contour(data = noaamat, aes(x=x, y=y, z=z, color=stat(level)),
               breaks=c(-20,-40,-60,-120, -500, -1000, -2000, -3000, -4000),
               size=c(0.5), alpha=0.4) +
  scale_color_gradient(name="Depth",breaks=c(-20, -40, -60), limits=c(-120,0), oob=squish)


#add annotations
batmap = batmap +
  annotate(geom="text", x = -159, y = 68.8, label = "Alaska", fontface = 2, color = "white", size = 6) +
  annotate(geom="label", x = -171.5, y = 71.6, label = 'bold(atop("Chukchi", "Sea"))', fontface = 2, size = 6, parse=TRUE, fill="white", label.size=NA) +
  annotate(geom="label", x = -170, y = 62.5, label = "bold(Bering~Sea)", fontface = 2, size = 6, parse=TRUE, fill="white", label.size=NA) +
  annotate(geom="text", hjust = 0, label = "Utqiaġvik", y = 71.5, x = -157.5, fontface=2, size = 3, color="black", label.size=NA, parse=TRUE) + 
  annotate(geom="text", hjust = 0, label = "Wainright", y = 70.64, x = -159, fontface=2, size = 3, color="white", label.size=NA, parse=TRUE) + 
  annotate(geom="text", hjust = 0, label = "Point~Lay", y = 69.74, x = -162.4, fontface=2, size = 3, color="white", label.size=NA, parse=TRUE) + 
  annotate(geom="text", hjust = 0, label = "Point~Hope", y = 68.35, x = -165.8, fontface=2, size = 3, color="white", label.size=NA, parse=TRUE) + 
  annotate(geom="text", hjust = 0, label = "Kotzebue", y = 67.3, x = -162.5, fontface=2, size = 3, color="white", label.size=NA, parse=TRUE) + 
  annotate(geom="text", hjust = 0, label = "Wales", y = 65.61, x = -167.3, fontface=2, size = 3, color="white", label.size=NA, parse=TRUE) +
  annotate(geom="text", hjust = 0, label = "Nome", y = 64.80, x = -165.7, fontface=2, size = 3, color="white", label.size=NA, parse=TRUE) +
  annotate(geom="text", label = "o", y = 71.387814, x = -156.481068, fontface=2, size = 4, color="blue", label.size=NA) +
  annotate(geom="text", label = "o", y = 70.642511, x = -159.966781, fontface=2, size = 4, color="blue", label.size=NA) +
  annotate(geom="text", label = "o", y = 69.737615, x = -162.967765, fontface=2, size = 4, color="blue", label.size=NA) +
  annotate(geom="text", label = "o", y = 68.348533, x = -166.800202, fontface=2, size = 4, color="blue", label.size=NA) +
  annotate(geom="text", label = "o", y = 66.965460, x = -162.628478, fontface=2, size = 4, color="blue", label.size=NA) +
  annotate(geom="text", label = "o", y = 65.609061, x = -168.079459, fontface=2, size = 4, color="blue", label.size=NA) +
  annotate(geom="text", label = "o", y = 64.503212, x = -165.400909, fontface=2, size = 4, color="blue", label.size=NA)

bbmap = bbmap +
  annotate(geom="text", x = -159, y = 68.8, label = "Alaska", fontface = 2, color = "white", size = 6) +
  annotate(geom="label", x = -171.5, y = 71.6, label = 'bold(atop("Chukchi", "Sea"))', fontface = 2, size = 6, parse=TRUE, fill="white", label.size=NA) +
  annotate(geom="label", x = -170, y = 62.5, label = "bold(Bering~Sea)", fontface = 2, size = 6, parse=TRUE, fill="white", label.size=NA) +
  annotate(geom="text", hjust = 0, label = "Utqiaġvik", y = 71.5, x = -157.5, fontface=2, size = 3, color="black", label.size=NA, parse=TRUE) + 
  annotate(geom="text", hjust = 0, label = "Wainright", y = 70.64, x = -159, fontface=2, size = 3, color="white", label.size=NA, parse=TRUE) + 
  annotate(geom="text", hjust = 0, label = "Point~Lay", y = 69.74, x = -162.4, fontface=2, size = 3, color="white", label.size=NA, parse=TRUE) + 
  annotate(geom="text", hjust = 0, label = "Point~Hope", y = 68.35, x = -165.8, fontface=2, size = 3, color="white", label.size=NA, parse=TRUE) + 
  annotate(geom="text", hjust = 0, label = "Kotzebue", y = 67.3, x = -162.5, fontface=2, size = 3, color="white", label.size=NA, parse=TRUE) + 
  annotate(geom="text", hjust = 0, label = "Wales", y = 65.61, x = -167.3, fontface=2, size = 3, color="white", label.size=NA, parse=TRUE) +
  annotate(geom="text", hjust = 0, label = "Nome", y = 64.80, x = -165.7, fontface=2, size = 3, color="white", label.size=NA, parse=TRUE) +
  annotate(geom="text", label = "o", y = 71.387814, x = -156.481068, fontface=2, size = 4, color="blue", label.size=NA) +
  annotate(geom="text", label = "o", y = 70.642511, x = -159.966781, fontface=2, size = 4, color="blue", label.size=NA) +
  annotate(geom="text", label = "o", y = 69.737615, x = -162.967765, fontface=2, size = 4, color="blue", label.size=NA) +
  annotate(geom="text", label = "o", y = 68.348533, x = -166.800202, fontface=2, size = 4, color="blue", label.size=NA) +
  annotate(geom="text", label = "o", y = 66.965460, x = -162.628478, fontface=2, size = 4, color="blue", label.size=NA) +
  annotate(geom="text", label = "o", y = 65.609061, x = -168.079459, fontface=2, size = 4, color="blue", label.size=NA) +
  annotate(geom="text", label = "o", y = 64.503212, x = -165.400909, fontface=2, size = 4, color="blue", label.size=NA)


basemap = basemap +
  annotate(geom="text", x = -159, y = 68.8, label = "Alaska", fontface = 2, color = "white", size = 6) +
  annotate(geom="label", x = -171.5, y = 71.6, label = 'bold(atop("Chukchi", "Sea"))', fontface = 2, size = 6, parse=TRUE, fill="white", label.size=NA) +
  annotate(geom="label", x = -170, y = 62.5, label = "bold(Bering~Sea)", fontface = 2, size = 6, parse=TRUE, fill="white", label.size=NA) +
  annotate(geom="text", hjust = 0, label = "Utqiaġvik", y = 71.5, x = -157.5, fontface=2, size = 3, color="black", label.size=NA, parse=TRUE) + 
  annotate(geom="text", hjust = 0, label = "Wainright", y = 70.64, x = -159, fontface=2, size = 3, color="white", label.size=NA, parse=TRUE) + 
  annotate(geom="text", hjust = 0, label = "Point~Lay", y = 69.74, x = -162.4, fontface=2, size = 3, color="white", label.size=NA, parse=TRUE) + 
  annotate(geom="text", hjust = 0, label = "Point~Hope", y = 68.35, x = -165.8, fontface=2, size = 3, color="white", label.size=NA, parse=TRUE) + 
  annotate(geom="text", hjust = 0, label = "Kotzebue", y = 67.3, x = -162.5, fontface=2, size = 3, color="white", label.size=NA, parse=TRUE) + 
  annotate(geom="text", hjust = 0, label = "Wales", y = 65.61, x = -167.3, fontface=2, size = 3, color="white", label.size=NA, parse=TRUE) +
  annotate(geom="text", hjust = 0, label = "Nome", y = 64.80, x = -165.7, fontface=2, size = 3, color="white", label.size=NA, parse=TRUE) +
  annotate(geom="text", label = "o", y = 71.387814, x = -156.481068, fontface=2, size = 4, color="blue", label.size=NA) +
  annotate(geom="text", label = "o", y = 70.642511, x = -159.966781, fontface=2, size = 4, color="blue", label.size=NA) +
  annotate(geom="text", label = "o", y = 69.737615, x = -162.967765, fontface=2, size = 4, color="blue", label.size=NA) +
  annotate(geom="text", label = "o", y = 68.348533, x = -166.800202, fontface=2, size = 4, color="blue", label.size=NA) +
  annotate(geom="text", label = "o", y = 66.965460, x = -162.628478, fontface=2, size = 4, color="blue", label.size=NA) +
  annotate(geom="text", label = "o", y = 65.609061, x = -168.079459, fontface=2, size = 4, color="blue", label.size=NA) +
  annotate(geom="text", label = "o", y = 64.503212, x = -165.400909, fontface=2, size = 4, color="blue", label.size=NA)

basemap


# prepare ggplot data frame
batdf = readRDS("data/station_list.RDS")

batdf$ptsize = 1
batdf$ptsize[batdf$station_type=="Hydrographic"] = 0
batdf = batdf[order(batdf$station_type,decreasing=TRUE),] #reorder for plotting purposes

#add station points
stnmap = batmap +
  geom_point(data=batdf, aes(x=lon, y=lat, fill=station_type, size=ptsize, shape=project), stroke=0) +
  scale_shape_manual(name="Project", values=p.shape) +
  guides(shape = guide_legend(override.aes=list(fill="black",color="black",size=5))) +
  guides(fill = guide_legend(override.aes=list(size=5, color=gg_color_hue(3)))) +
  guides(col=guide_legend(ncol = 1,override.aes = list(size=3))) +
  scale_fill_manual(name="Station Type",values=gg_color_hue(3)) +
  scale_size_continuous(range=c(2,3),limits=c(0,1), guide = FALSE) +
  theme(panel.spacing = unit(2, "lines"), strip.text = element_text(size=14)) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.text.align=0, legend.box.just = "left")

#stnmap

ggsave(plot=stnmap, filename="output/05_figure_map_bathy_station.png", device=png(), width=6, height=8)


#add station points to filled bathy map
stnmap2 = bbmap +
  geom_point(data=batdf, aes(x=lon, y=lat, col=station_type, size=ptsize, shape=project)) +
  guides(shape = guide_legend(override.aes=list(shape=p.shape,fill="black",color="black",size=5))) +
  scale_color_manual(name="Station Type",values=gg_color_hue(3)) +
  scale_size_continuous(range=c(2,3),limits=c(0,1), guide = FALSE) +
  theme(panel.spacing = unit(2, "lines"), strip.text = element_text(size=14)) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.text.align=0, legend.box.just = "left")

stnmap2

ggsave(plot=stnmap2, filename="output/05_figure_map_bathy_stations2.png", device=png(), width=6, height=8)



#fix facet labels
depthnames=c(`bottom` = "Bottom",`midwater` = "Midwater",`surface`="Surface",
             `AMBON 2015`="AMBON 215",`AMBON 2017`="AMBON 2017",`ASGARD 2017`="ASGARD 2017",`DBO-NCIS 2017`="DBO-NCIS 2017")

metavar = c("depth_m","lat","lon","temp","salinity","DO",
            "FlECO-AFL(mg/m^3)","chl (ug/l)","phaeo (ug/l)","phaeoratio","PO4(uM)","Sil(uM)",
            "NH4(uM)","N+N (umol/L)","transmission_pct","stratification")
metavarlog = c("FlECO-AFL(mg/m^3)","chl (ug/l)","phaeo (ug/l)","PO4(uM)","Sil(uM)",
               "NO3(uM)","NO2(uM)","NH4(uM)","N+N (umol/L)")
metaplot = c("lat","lon","temp","salinity","stratification","DO","transmission_pct",
             "FlECO-AFL(mg/m^3)","phaeofrac","PO4(uM)","Sil(uM)",
             "N+N (umol/L)","NH4(uM)","depth_type","mass","project")

#AMBON doesn't have NO3, only total
#DBO doesn't have transmission values

z.meta = x.meta
# take some logs?
# z.meta[,metavarlog] = log10(z.meta[,metavarlog])
z.meta = z.meta[,metaplot]
colnames(z.meta) = make.names(colnames(z.meta))
znames = c("Latitude\n(°N)","Longitude\n(°E)","Temperature\n(°C)","Salinity\n(PSU)",
           "Stratification\nIndex","Dissolved Oxygen\n(\u00b5mol/kg)",
           "Beam Transmission\n(%)", "Chl a Fluoresence\n(mg/m\u00B3)",
           "Phaeopigment\nfraction",
           "Phosphate\n(\u00b5M)", "Silicate\n(\u00b5M)","Total Nitrate\n(\u00b5M)",
           "Ammonium\n(\u00b5M)","Depth Bin","Water Mass","Project")


### plot metadata maps
metamp = list()
for(i in 1:(ncol(z.meta)-3)) {
  coltax = colnames(z.meta)[i]
  
  s.shape = c(19,17,18,15)[1:length(unique(x.meta$project))]
  metamp[[coltax]] = basemap + 
    geom_point(data=z.meta, aes_string(x="lon", y="lat", color="mass", size=coltax)) +
    # geom_text(data=z.meta, aes(x= lon, y = lat, label=station), hjust=1, vjust=-.5, size = 1, check_overlap=TRUE) +
    theme(legend.position="right") + theme(panel.spacing = unit(2, "lines"), strip.text = element_text(size=14)) +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.text.align=0, legend.box.just = "left") +
    # guides(col=guide_legend(nrow=1,byrow=TRUE)) + 
    ggtitle(label=znames[i]) + labs(size=znames[i]) + 
    scale_size_continuous(range=c(0,4)) +
    # scale_shape_manual(values=s.shape) +
        guides(fill = guide_legend(order = 2), 
           shape = guide_legend(order = 1),
           size = guide_legend(order=3)) +
    facet_wrap(vars(project,depth_type), labeller = as_labeller(depthnames), ncol = 3)

}

pdf(file="output/05_figure_map_metadata.pdf",width=12,height=18)
invisible(lapply(metamp, print))
dev.off()




### plot cluster and taxa maps
for(tax in names(mat.list)) {

  tx = mat.list[[tax]][["title"]]
  abb = mat.list[[tax]][["short"]]
  x.meta = mat.list[[tax]][["meta"]]
  x.taxa = mat.list[[tax]][["prop"]]
  nclus = clus.list[[tax]]
  clusnames = paste0(abb, 1:nclus)
  x.top.names = mat.list[[tax]][["top.names"]]
  x.top.longnames = mat.list[[tax]][["top.longnames"]] 

  
# cluster map
dmp = basemap +
    geom_point(data=x.meta, aes(x=lon, y=lat, color=cluster), size=2) +
    # geom_text(data=x.meta, aes(x= lon, y = lat, label=station), hjust=1, vjust=-.5, size = 1, check_overlap=TRUE) +
    theme(legend.position="right") + theme(panel.spacing = unit(2, "lines"), strip.text = element_text(size=14)) +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.text.align=0, legend.box.just = "left") +
    guides(col=guide_legend(ncol=1,byrow=TRUE)) + ggtitle(tx) +
    scale_size_continuous(range=c(0,4)) +
    facet_wrap(vars(project, depth_type), labeller = as_labeller(depthnames), ncol = 3)
  
ggsave(plot=dmp, filename=paste0("output/05_figure_map_",tax,"_clusters.png"),width=14,height=18,device=png())


# cluster map with bathy
dmpb = batmap +
  geom_point(data=x.meta, aes(x=lon, y=lat, fill=cluster), size=1.5, shape=21, stroke=0) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.text.align=0, legend.box.just = "left") +
  theme(legend.position="right") +
  guides(shape = guide_legend(override.aes=list(size=3))) +
  guides(color = guide_legend(override.aes=list(size=3))) +
  guides(col=guide_legend(ncol=1,byrow=TRUE)) + ggtitle(tx) +    
  facet_wrap(vars(project, depth_type), labeller = as_labeller(depthnames), ncol=3) +
  theme(panel.spacing = unit(2, "lines"), strip.text = element_text(size=14))

ggsave(plot=dmpb, filename=paste0("output/05_figure_map_bathy_",tax,"_clusters.png"),width=14,height=18,device=png())



# top taxa cluster maps

  df = data.frame(x.meta,x.taxa*1000)
  df[df==0] = NA
  dtopmp = list()

  for(i in 1:length(x.top.names)) {
  coltax = make.names(x.top.names[i])
  ttax = x.top.longnames[i]
  ttax = paste0(long_names[taxa_split$ESV==coltax]," (",sprintf("%.1f",round(colMeans(x.taxa[,x.top.names])[i]*1000,1)), ")")

  dtopmp[[coltax]] = basemap + 
  geom_point(data=df, aes_string(x="lon", y="lat", color="cluster", size=coltax)) +
  geom_text(data=df, aes(x= lon, y = lat, label=station), hjust=1, vjust=-.5, size = 1, check_overlap=TRUE) +
  theme(legend.position="right") + theme(panel.spacing = unit(2, "lines"), strip.text = element_text(size=14)) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.text.align=0, legend.box.just = "left") +
  # guides(col=guide_legend(nrow=1,byrow=TRUE)) + 
  ggtitle(label=ttax) + labs(size="Relative\nAbundance\n(per mille)") + 
    guides(col = guide_legend(order = 2), 
           shape = guide_legend(order = 1),
           size = guide_legend(order = 3)) +
    scale_size_continuous(range=c(0,4)) +
      facet_wrap(vars(project, depth_type), labeller = as_labeller(depthnames), ncol=3) +
    scale_size_continuous(breaks=c(5,50,250,500))

  }

pdf(file=paste0("output/05_figure_map_",tax,"_tops.pdf"),width=12,height=16)
invisible(lapply(dtopmp, print))
dev.off()


graphics.off()


} #end taxa loop
