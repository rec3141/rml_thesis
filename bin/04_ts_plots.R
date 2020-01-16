## T/S Plots

# require(grid)
# require(gridExtra)
# require(ggplot2)
# library(ggpubr)
# require(stringr)
# require(cowplot)
require(viridis)

# p.shape = c(21,24,25,22)[1:length(unique(sub.meta$project))]
p.shape = c(21,24,22)

# # water mass designations from Pisareva et al. 2015
# wmin = data.frame(x1=c(29.7,30.5, 24, 32,30.5,24,29,30.5,33.64,29,24),
#   x2=c(33.64, 33.64, 29.7,33.64,33.64,30.5,30.5,32,35,30.5,29),
#   y1=c(5, 3, 5,3,0,2,2,3,1,3,2),
#   y2=c(15, 15, 15,0,-1.8,-1.6,3,0,-1.25,5,5))

wmout = matrix(data=NA, ncol=4, nrow=4*nrow(masses))
for(i in 1:nrow(masses)) {
    wmout[4*i-3,] = as.numeric(masses[i, c("xmin","xmax","ymin","ymin")])
    wmout[4*i-2,] = as.numeric(masses[i, c("xmin","xmax","ymax","ymax")])
    wmout[4*i-1,] = as.numeric(masses[i, c("xmin","xmin","ymin","ymax")])
    wmout[4*i-0,] = as.numeric(masses[i, c("xmax","xmax","ymin","ymax")])
}
wmout = data.frame("mass"=rep(masses$mass,rep(4,nrow(masses))),wmout)
colnames(wmout) = colnames(masses)
wmout = data.frame(wmout,stringsAsFactors=FALSE)
wmout = wmout[!(duplicated(wmout) | duplicated(wmout, fromLast = TRUE)), ]

#remove extraneous lines
#wmout = wmout[-c(1,7,21,24,25,28,31,39,44),]

# base T/S plot with segments
baseplot = ggplot(data=sub.meta, aes(x=salinity, y=temp)) +
  geom_segment(data=wmout, aes(x=xmin, xend=xmax, y=ymin, yend=ymax), color='grey') +
  labs(x='Salinity (PSU)', y='Temperature (°C)') +
  theme(legend.position="none") +
  scale_y_continuous(breaks=seq(-2,15,2)) +
  scale_x_continuous(breaks=seq(24,36,2)) + 
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16), legend.text.align=0, legend.box.just = "left")

# water mass plot
mts = baseplot + 
  geom_point(data=sub.meta, aes(x=salinity, y=temp, fill=mass, shape=project), stroke=0) +
  labs(shape='Project', fill='Water Mass') +
  scale_shape_manual(name="Project", na.translate=FALSE,
                     values=p.shape) +
  scale_fill_discrete(name="Water Mass", na.translate=FALSE)

ggsave(plot = mts, filename = "output/04_figure_ts_watermass.png",device = png(),width = 6, height=6)

#ggplot bug when using scale_shape_manual
#https://github.com/tidyverse/ggplot2/issues/2322
mtsl = mts + theme(legend.position="bottom", legend.direction="vertical") + 
  guides( 
    fill=guide_legend(nrow=3,byrow=TRUE, 
                      override.aes = list(shape=21)), 
    shape=guide_legend(override.aes = list(fill="black",shape=p.shape)))
mts_legend <- get_legend(mtsl)
mtsleg = as_ggplot(mts_legend)
ggsave(plot=mtsleg, filename = "output/04_figure_ts_watermass_leg.png",device = png(),width = 8, height=2)


# nutrients plot
nts = baseplot + 
  geom_point(data=sub.meta, aes(x=salinity, y=temp, fill=`N+N (umol/L)`, shape=project, size=`N+N (umol/L)`), stroke=0) + 
  labs(shape='Project', fill="Total Nitrate\n(\u00b5M)", size="Total Nitrate\n(\u00b5M)") +
  scale_shape_manual(name="Project", na.translate=FALSE, guide=FALSE,
                       values=p.shape,
                       breaks=c("ASGARD2017", "AMBON2017", "AMBON2017", "DBONCIS2017"),
                       labels=c("ASGARD 2017", "AMBON 2015", "AMBON 2017", "DBO-NCIS 2017")) +
  scale_fill_viridis_c() +
  # scale_color_discrete(name="Water Mass", na.translate=FALSE, guide=FALSE,
  #                      breaks=c('FACW','ACW','AW','BSW','MWR','RWW','WW'),
  #                      labels=c('FACW','ACW','AW','BSW','MWR','RWW','WW')) +
  # scale_alpha(name="Total Nitrate\n(\u00b5M)",
  #             breaks=c(0,5,10,15),
  #             labels=c('0','5','10','15'),
  #             range = c(0.1, 0.3)) +
  # scale_shape_manual(name="Total Nitrate\n(\u00b5M)", guide=FALSE, values = c(21:23)) +
  scale_radius(breaks=c(0,5,10,15),
             labels=c('0','5','10','15'),
             range=c(1,5))
  
ggsave(plot=nts, filename = "output/04_figure_ts_nutrients.png",device = png(),width = 6, height=6)

ntsl = nts + theme(legend.position="bottom", legend.direction="vertical") + guides(size = guide_legend(reverse = TRUE))
nts_legend <- get_legend(ntsl)
ntsleg = as_ggplot(nts_legend)
ggsave(plot=ntsleg, filename = "output/04_figure_ts_nutrients_leg.png",device = png(),width = 8, height=2)




# DO plot
dots = baseplot + 
  geom_point(data=sub.meta, aes(x=salinity, y=temp, fill=DO, shape=project, size=DO), stroke=0) + 
  labs(shape='Project', fill="Dissolved\nOxygen\n(\u00b5mol/kg)", size="Dissolved\nOxygen\n(\u00b5mol/kg)") +
  scale_shape_manual(name="Project", na.translate=FALSE, guide=FALSE,
                       values=p.shape,
                       breaks=c("ASGARD2017", "AMBON2015", "AMBON2017", "DBONCIS2017"),
                       labels=c("ASGARD 2017", "AMBON 2015", "AMBON 2017", "DBO-NCIS 2017")) +
  scale_fill_viridis_c() +
  # scale_color_discrete(name="Water Mass", na.translate=FALSE, guide=FALSE,
  #                      breaks=c('FACW','ACW','AW','BSW','MWR','RWW','WW'),
  #                      labels=c('FACW','ACW','AW','BSW','MWR','RWW','WW')) +
  scale_radius(breaks=c(200,250, 300, 350, 400), range = c(1,5),
              labels=c('200','250', '300','350','400'))


ggsave(plot=dots, filename = "output/04_figure_ts_oxygen.png",device = png(),width = 6, height=6)

dotsl = dots + theme(legend.position="bottom", legend.direction="vertical") + guides(size = guide_legend(reverse = TRUE))
dots_legend <- get_legend(dotsl)
dotsleg = as_ggplot(dots_legend)
ggsave(plot=dotsleg, filename = "output/04_figure_ts_oxygen_leg.png",device = png(),width = 8, height=2)


# chla plot
cts = baseplot + 
  geom_point(data=sub.meta, aes(x=salinity, y=temp, fill=`FlECO-AFL(mg/m^3)`, shape=project, size=`FlECO-AFL(mg/m^3)`), stroke=0) + 
  labs(shape='Project', fill="Chlorophyll\nFluorescence\n(mg/m\u00B3)", size="Chlorophyll\nFluorescence\n(mg/m\u00B3)") +
  scale_shape_manual(name="Project", na.translate=FALSE, guide=FALSE,
                       values=p.shape,
                       breaks=c("ASGARD2017", "AMBON2015", "AMBON2017", "DBONCIS2017"),
                       labels=c("ASGARD 2017", "AMBON 2017", "AMBON 2017", "DBO-NCIS 2017")) +
  scale_fill_viridis_c(breaks=c(0, 1, 2, 5, 10),
                        labels=c('0','1','2','5','10'),
                        oob = scales::squish, limits=c(0,10)) +
  # scale_color_discrete(name="Water Mass", na.translate=FALSE, guide=FALSE,
  #                      breaks=c('FACW','ACW','AW','BSW','MWR','RWW','WW'),
  #                      labels=c('FACW','ACW','AW','BSW','MWR','RWW','WW')) +
  scale_radius(breaks=c(0, 1, 2, 5, 10),
              labels=c('0','1','2','5','10'),
              range = c(1,5), limits = c(0,10))


ggsave(plot=cts, filename = "output/04_figure_ts_chla.png",device = png(),width = 6, height=6)

ctsl = cts + theme(legend.position="bottom", legend.direction="vertical") + guides(size = guide_legend(reverse = TRUE))
cts_legend <- get_legend(ctsl)
ctsleg = as_ggplot(cts_legend)
ggsave(plot=ctsleg, filename = "output/04_figure_ts_chla_leg.png",device = png(),width = 8, height=2)



# combine and arrange metadata plots
ggar = ggarrange(mts, nts, mtsleg, ntsleg, dots, cts, dotsleg, ctsleg,
   nrow=4,heights=c(4,2.5,4,2.5))
ggar
ggsave(plot=ggar, file="output/04_figure_ts_combined_meta.png",device=png(),width=10,height=8*1.615)





 ### taxonomy clusters T/S plot
plot.list = list()

for(tax in names(mat.list)) {

  x.meta = mat.list[[tax]][["meta"]]

    pts = baseplot + 
    geom_point(data=x.meta, aes(x=salinity, y=temp, fill=cluster, shape=project), stroke=0) + 
    scale_fill_manual(name="Cluster", values=gg_color_hue(length(levels(x.meta$cluster)))) +
    scale_shape_manual(guide=FALSE, values=p.shape)

    plot.list[[tax]][["pts"]] = pts  
  ggsave(plot=pts, filename = paste0("output/04_figure_ts_",tax,".png"),device = png(),width = 6, height=6)
  
  ptsl = pts + theme(legend.position="bottom", legend.direction = "vertical") + guides( 
    fill=guide_legend(nrow=2,byrow=TRUE, 
                      override.aes = list(shape=21)), 
    shape=guide_legend(override.aes = list(fill="black",shape=p.shape)))
  
  pts_legend <- get_legend(ptsl)
  ptsleg = as_ggplot(pts_legend)
  plot.list[[tax]][["leg"]] = ptsleg
  ggsave(plot=ptsleg, filename = paste0("output/04_figure_ts_",tax,"_leg.png"),device = png(), width = 8, height=1)
  
}


ggar2 = ggarrange(plot.list[["diatoms"]][["pts"]], plot.list[["picos"]][["pts"]], plot.list[["diatoms"]][["leg"]], plot.list[["picos"]][["leg"]],
  nrow=2, heights=c(4,1))
ggar2
ggsave(plot=ggar2, file="output/04_figure_ts_combined_clusters.png",device=png(),width=10,height=5*1.25)


#######

graphics.off()

source("bin/05_maps.R")
