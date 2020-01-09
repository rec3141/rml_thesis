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
xmax = 210-5
ymin = 62
ymax = 74

# download bathymetry data
noaadat <- getNOAA.bathy(xmin-360,xmax-360,ymin,ymax,res=1, keep=TRUE)
noaamat <- fortify(noaadat)   # converts b class 'bathy' to data frame
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
xmax = -150-5



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
basemap <- basemap + 
  labs( x = expression("Longitude ("*degree*"W)"), y = expression("Latitude ("*degree*"N)")) +
  scale_x_continuous(breaks=c(-175,-170,-165,-160,-155),labels=c("175","170","165","160","155"))

bbmap <- bbmap + 
  labs( x = expression("Longitude ("*degree*"W)"), y = expression("Latitude ("*degree*"N)")) +
  scale_x_continuous(breaks=c(-175,-170,-165,-160,-155),labels=c("175","170","165","160","155"))

# add bathymetry contours to basemap
batmap <- basemap + 
  geom_contour(data = noaamat, aes(x=x, y=y, z=z, color=stat(level)),
               breaks=c(-20,-40,-60,-120, -500, -1000, -2000, -3000, -4000),
               size=c(0.5), alpha=0.4) +
  scale_color_gradient(name="Depth",breaks=c(-20, -40, -60), limits=c(-120,0), oob=squish)


#add annotations
batmap = batmap +
  annotate(geom="text", x = -159, y = 68.3, label = "Alaska", fontface = 2, color = "white", size = 6) +
  annotate(geom="label", x = -171.5, y = 71.6, label = 'bold(atop("Chukchi", "Sea"))', fontface = 2, size = 6, parse=TRUE, fill="white", label.size=NA) +
  annotate(geom="label", x = -170, y = 62.5, label = "bold(Bering~Sea)", fontface = 2, size = 6, parse=TRUE, fill="white", label.size=NA)

bbmap = bbmap +
  annotate(geom="text", x = -159, y = 68.3, label = "Alaska", fontface = 2, color = "white", size = 6) +
  annotate(geom="text", x = -171.5, y = 71.6, label = 'bold(atop("Chukchi", "Sea"))', fontface = 2, size = 6, parse=TRUE, fill="white", label.size=NA) +
  annotate(geom="text", x = -170, y = 62.5, label = "bold(Bering~Sea)", fontface = 2, size = 6, parse=TRUE, fill="white", label.size=NA)

basemap = basemap +
  annotate(geom="text", x = -159, y = 68.3, label = "Alaska", fontface = 2, color = "white", size = 6) +
  annotate(geom="label", x = -171.5, y = 71.6, label = 'bold(atop("Chukchi", "Sea"))', fontface = 2, size = 6, parse=TRUE, fill="white", label.size=NA) +
  annotate(geom="label", x = -170, y = 62.5, label = "bold(Bering~Sea)", fontface = 2, size = 6, parse=TRUE, fill="white", label.size=NA)

# prepare ggplot data frame
batdf = readRDS("data/station_list.RDS")

batdf$ptsize = 1
batdf$ptsize[batdf$station_type=="Hydrographic"] = 0
batdf = batdf[order(batdf$station_type,decreasing=TRUE),] #reorder for plotting purposes

#add station points
stnmap = batmap +
  geom_point(data=batdf, aes(x=lon, y=lat, fill=station_type, size=ptsize, shape=project), stroke=0) +
  scale_shape_manual(name="Project", values=c(21,24,22)) +
  guides(shape = guide_legend(override.aes=list(shape=c(21,24,22),fill="black",color="black",size=5))) +
  guides(fill = guide_legend(override.aes=list(size=5,color=gg_color_hue(3)))) +
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
  guides(shape = guide_legend(override.aes=list(shape=c(21,24,22),fill="black",color="black",size=5))) +
  scale_color_manual(name="Station Type",values=gg_color_hue(3)) +
  scale_size_continuous(range=c(2,3),limits=c(0,1), guide = FALSE) +
  theme(panel.spacing = unit(2, "lines"), strip.text = element_text(size=14)) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.text.align=0, legend.box.just = "left")

stnmap2

ggsave(plot=stnmap2, filename="output/05_figure_map_bathy_stations2.png", device=png(), width=6, height=8)



#fix facet labels
depthnames=c(`bottom` = "Bottom",`midwater` = "Midwater",`surface`="Surface")

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
z.meta[,metavarlog] = log10(z.meta[,metavarlog])
z.meta = z.meta[,metaplot]
colnames(z.meta) = make.names(colnames(z.meta))
znames = c("Latitude\n(°N)","Longitude\n(°E)","Temperature\n(°C)","Salinity\n(PSU)",
           "Stratification Index","Dissolved Oxygen\n(\u00b5mol/kg)",
           "Beam Transmission\n(%)", "Chl a Fluoresence\n(mg/m\u00B3)",
           "Phaeopigment fraction",
           "Phosphate\n(\u00b5M)", "Silicate\n(\u00b5M)","Total Nitrate\n(\u00b5M)",
           "Ammonium\n(\u00b5M)","Depth Bin","Water Mass","Project")


### plot metadata maps
metamp = list()
for(i in 1:(ncol(z.meta)-3)) {
  coltax = colnames(z.meta)[i]
  ttax = znames[i]
  
  metamp[[coltax]] = basemap + 
    geom_point(data=z.meta, aes_string(x="lon", y="lat", color="mass", shape="project", size=coltax)) +
    # geom_text(data=z.meta, aes(x= lon, y = lat, label=station), hjust=1, vjust=-.5, size = 1, check_overlap=TRUE) +
    theme(legend.position="right") + theme(panel.spacing = unit(2, "lines"), strip.text = element_text(size=14)) +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.text.align=0, legend.box.just = "left") +
    # guides(col=guide_legend(nrow=1,byrow=TRUE)) + 
    ggtitle(label=ttax) + labs(size=ttax) + 
    guides(col = guide_legend(order = 2), 
           shape = guide_legend(order = 1),
           size = guide_legend(order=3)) +
    facet_wrap(vars(depth_type), labeller = as_labeller(depthnames))

}

pdf(file="output/05_figure_map_metadata.pdf",width=12,height=6)
invisible(lapply(metamp, print))
dev.off()




### plot cluster and taxa maps


# diatom cluster map
dmp <- basemap +
    geom_point(data=d.meta, aes(x=lon, y=lat, color=cluster, shape=project), size=2) +
    # geom_text(data=d.meta, aes(x= lon, y = lat, label=station), hjust=1, vjust=-.5, size = 1, check_overlap=TRUE) +
    theme(legend.position="right") + theme(panel.spacing = unit(2, "lines"), strip.text = element_text(size=14)) +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.text.align=0, legend.box.just = "left") +
    guides(col=guide_legend(ncol=1,byrow=TRUE)) +
    facet_wrap(vars(depth_type), labeller = as_labeller(depthnames))
  
ggsave(plot=dmp, filename="output/05_figure_map_diatom_clusters.png",width=14,height=8,device=png())


# diatom cluster map with bathy
dmpb <- batmap +
  geom_point(data=d.meta, aes(x=lon, y=lat, fill=cluster, shape=project), size=1.5, stroke=NA) +
  # geom_text(data=d.meta, aes(x= lon, y = lat, label=station), hjust=1, vjust=-.5, size = 1, check_overlap=TRUE) +
  scale_shape_manual(values=c(21,24,22)) +
  guides(shape = guide_legend(override.aes=list(shape=c(21,24,22),fill="black",color="black",size=5))) +
  guides(fill = guide_legend(override.aes=list(size=5,color=gg_color_hue(ndclus)))) +
  theme(legend.position="right") + theme(panel.spacing = unit(2, "lines"), strip.text = element_text(size=14)) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.text.align=0, legend.box.just = "left") +
  guides(col=guide_legend(ncol=1,byrow=TRUE)) +
  facet_wrap(vars(depth_type), labeller = as_labeller(depthnames))

ggsave(plot=dmpb, filename="output/05_figure_map_bathy_diatom_clusters.png",width=14,height=8,device=png())


# picos cluster map
pmp <- basemap + 
  geom_point(data=p.meta, aes(x=lon, y=lat, color=cluster, shape=project), size=2) +
  # geom_text(data=p.meta, aes(x= lon, y = lat, label=station), hjust=1, vjust=-.5, size = 1, check_overlap=TRUE) +
  theme(legend.position="right") + theme(panel.spacing = unit(2, "lines"), strip.text = element_text(size=14)) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.text.align=0, legend.box.just = "left") +
  guides(col=guide_legend(ncol=1,byrow=TRUE)) +
  facet_wrap(vars(depth_type), labeller = as_labeller(depthnames))

ggsave(plot=pmp, filename="output/05_figure_map_picos_clusters.png",width=14,height=8,device=png())



# picos cluster map with bathy
pmpb <- batmap +
  geom_point(data=p.meta, aes(x=lon, y=lat, fill=cluster, shape=project), size=1.5, stroke=NA) +
  # geom_text(data=p.meta, aes(x= lon, y = lat, label=station), hjust=1, vjust=-.5, size = 1, check_overlap=TRUE) +
  scale_shape_manual(values=c(21,24,22)) +
  guides(shape = guide_legend(override.aes=list(shape=c(21,24,22),fill="black",color="black",size=5))) +
  guides(fill = guide_legend(override.aes=list(size=5,color=gg_color_hue(npclus)))) +
  theme(legend.position="right") + theme(panel.spacing = unit(2, "lines"), strip.text = element_text(size=14)) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.text.align=0, legend.box.just = "left") +
  guides(col=guide_legend(ncol=1,byrow=TRUE)) +
  facet_wrap(vars(depth_type), labeller = as_labeller(depthnames))

ggsave(plot=pmpb, filename="output/05_figure_map_bathy_picos_clusters.png",width=14,height=8,device=png())

#####

# top diatom cluster maps
df = data.frame(d.meta,d.prop*1000)
df[df==0] = NA

dtopmp = list()
for(i in 1:length(d.top.names)) {
  coltax = make.names(d.top.names[i])
  ttax = d.top.newnames[i]

  dtopmp[[coltax]] = basemap + 
  geom_point(data=df, aes_string(x="lon", y="lat", color="cluster", shape="project", size=coltax)) +
  geom_text(data=df, aes(x= lon, y = lat, label=station), hjust=1, vjust=-.5, size = 1, check_overlap=TRUE) +
  theme(legend.position="right") + theme(panel.spacing = unit(2, "lines"), strip.text = element_text(size=14)) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.text.align=0, legend.box.just = "left") +
  # guides(col=guide_legend(nrow=1,byrow=TRUE)) + 
  ggtitle(label=ttax) + labs(size="Relative\nAbundance\n(\u2030)") + 
    guides(col = guide_legend(order = 2), 
           shape = guide_legend(order = 1),
           size = guide_legend(order=3)) +
      facet_wrap(vars(depth_type), labeller = as_labeller(depthnames))
}

pdf(file="output/05_figure_map_diatom_tops.pdf",width=12,height=6)
invisible(lapply(dtopmp, print))
dev.off()




# top picos cluster maps
df = data.frame(p.meta,p.prop*1000)
df[df==0] = NA

ptopmp = list()
for(i in 1:length(p.top.names)) {
  coltax = make.names(p.top.names[i])
  ttax = p.top.newnames[i]

  ptopmp[[coltax]] = basemap + 
    geom_point(data=df, aes_string(x="lon", y="lat", color="cluster", shape="project", size=coltax)) +
    geom_text(data=df, aes(x= lon, y = lat, label=station), hjust=1, vjust=-.5, size = 1, check_overlap=TRUE) +
    theme(legend.position="right") + theme(panel.spacing = unit(2, "lines"), strip.text = element_text(size=14)) +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.text.align=0, legend.box.just = "left") +
    # guides(col=guide_legend(nrow=1,byrow=TRUE)) + 
    ggtitle(label=ttax) + labs(size="Relative\nAbundance\n(\u2030)") + 
    guides(col = guide_legend(order = 2),
           shape = guide_legend(order = 1),
           size = guide_legend(order=3)) +
    facet_wrap(vars(depth_type), labeller = as_labeller(depthnames))
}

pdf(file="output/05_figure_map_picos_tops.pdf",width=12,height=6)
invisible(lapply(ptopmp, print))
dev.off()

