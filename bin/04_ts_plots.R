## T/S Plots

# require(grid)
# require(gridExtra)
# require(ggplot2)
# library(ggpubr)
# require(stringr)
# require(cowplot)
require(viridis)


graphics.off()

# water mass designations from Pxx et al.
wmin = data.frame(x1=c(29.7,30.5, 24, 32,30.5,24,29,30.5,33.64,29,24),
  x2=c(33.64, 33.64, 29.7,33.64,33.64,30.5,30.5,32,35,30.5,29),
  y1=c(5, 3, 5,3,0,2,2,3,1,3,2),
  y2=c(15, 15, 15,0,-1.8,-1.6,3,0,-1.25,5,5))

wmout = NULL
for(i in 1:nrow(wmin)) {
  x=unlist(wmin[i,,drop=F])
  wmout = rbind(wmout,
                  c(x[1],x[2],x[3],x[3]),
                  c(x[1],x[2],x[4],x[4]),
                  c(x[1],x[1],x[3],x[4]),
                  c(x[2],x[2],x[3],x[4])
  )
}
colnames(wmout) = colnames(wmin)
wmout = data.frame(wmout)
#remove extraneous lines
wmout = wmout[-c(1,7,21,24,25,28,31,39,44),]

# base T/S plot with segments
baseplot = ggplot(data=meta, aes(x=salinity, y=temp)) +
  geom_segment(data=wmout, aes(x=x1, xend=x2, y=y1, yend=y2), color='grey') +
  labs(x='Salinity (PSU)', y='Temperature (°C)') +
  theme(legend.position="none") +
  scale_y_continuous(breaks=seq(-2,15,2)) +
  scale_x_continuous(breaks=seq(24,36,2)) + 
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16), legend.text.align=0, legend.box.just = "left")

# water mass plot
mts = baseplot + 
  geom_point(data=meta, aes(x=salinity, y=temp, color=mass, shape=project)) +
  labs(shape='Project', col='Water Mass') +
  scale_shape_discrete(name="Project", na.translate=FALSE,
                       breaks=c("ASGARD2017", "AMBON2017", "DBONCIS2017"),
                       labels=c("ASGARD 2017", "AMBON 2017", "DBO-NCIS 2017")) +
  scale_color_discrete(name="Water Mass", na.translate=FALSE,
                       breaks=c('FACW','ACW','AW','BSW','MWR','RWW','WW'),
                       labels=c('FACW','ACW','AW','BSW','MWR','RWW','WW'))

ggsave(plot = mts, filename = "output/04_figure_ts_watermass.png",device = png(),width = 6, height=6)

mtsl = mts + theme(legend.position="bottom", legend.direction="vertical") + 
  guides(col=guide_legend(nrow=3,byrow=TRUE))
mts_legend <- get_legend(mtsl)
mtsleg = as_ggplot(mts_legend)
ggsave(plot=mtsleg, filename = "output/04_figure_ts_watermass_leg.png",device = png(),width = 8, height=2)


# nutrients plot
nts = baseplot + 
  geom_point(data=meta, aes(x=salinity, y=temp, color=`N+N (umol/L)`, shape=project, size=`N+N (umol/L)`)) + 
  labs(shape='Project', col="Total Nitrate\n(\u00b5M)", size="Total Nitrate\n(\u00b5M)") +
  scale_shape_discrete(name="Project", na.translate=FALSE, guide=FALSE,
                       breaks=c("ASGARD2017", "AMBON2017", "DBONCIS2017"),
                       labels=c("ASGARD 2017", "AMBON 2017", "DBO-NCIS 2017")) +
  scale_color_viridis_c() +
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
  geom_point(data=meta, aes(x=salinity, y=temp, color=DO, shape=project, size=DO)) + 
  labs(shape='Project', col="Dissolved\nOxygen\n(\u00b5mol/kg)", size="Dissolved\nOxygen\n(\u00b5mol/kg)") +
  scale_shape_discrete(name="Project", na.translate=FALSE, guide=FALSE,
                       breaks=c("ASGARD2017", "AMBON2017", "DBONCIS2017"),
                       labels=c("ASGARD 2017", "AMBON 2017", "DBO-NCIS 2017")) +
  scale_colour_viridis_c() +
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
  geom_point(data=meta, aes(x=salinity, y=temp, color=`FlECO-AFL(mg/m^3)`, shape=project, size=`FlECO-AFL(mg/m^3)`)) + 
  labs(shape='Project', col="Chlorophyll\nFluorescence\n(mg/m\u00B3)", size="Chlorophyll\nFluorescence\n(mg/m\u00B3)") +
  scale_shape_discrete(name="Project", na.translate=FALSE, guide=FALSE,
                       breaks=c("ASGARD2017", "AMBON2017", "DBONCIS2017"),
                       labels=c("ASGARD 2017", "AMBON 2017", "DBO-NCIS 2017")) +
  scale_color_viridis_c(breaks=c(0, 1, 2, 5, 10),
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




# diatom clusters ts plot
dts = baseplot + 
  geom_point(data=d.meta, aes(x=salinity, y=temp, color=cluster, shape=project))

ggsave(plot=dts, filename = "output/04_figure_ts_diatom.png",device = png(),width = 6, height=6)

dtsl = dts + theme(legend.position="bottom") + 
  scale_color_manual(name="Cluster", values=gg_color_hue(length(unique(d.meta$cluster))),
                       breaks=sort(unique(d.meta$cluster)),
                       labels=sort(unique(d.meta$cluster))) +
  scale_shape_discrete(guide=FALSE)

dts_legend <- get_legend(dtsl)
dtsleg = as_ggplot(dts_legend)
ggsave(plot=dtsleg, filename = "output/04_figure_ts_diatom_leg.png",device = png(),width = 8, height=1)


# picos clusters ts plot
pts = baseplot + 
  geom_point(data=p.meta, aes(x=salinity, y=temp, color=cluster, shape=project)) + 
  scale_color_manual(name="Cluster", values=gg_color_hue(length(unique(p.meta$cluster))),
                     breaks=sort(unique(p.meta$cluster)),
                     labels=sort(unique(p.meta$cluster))) +
  scale_shape_discrete(guide=FALSE)

ggsave(plot=pts, filename = "output/04_figure_ts_picos.png",device = png(),width = 6, height=6)

ptsl = pts + theme(legend.position="bottom") 
pts_legend <- get_legend(ptsl)
ptsleg = as_ggplot(pts_legend)
ggsave(plot=ptsleg, filename = "output/04_figure_ts_picos_leg.png",device = png(), width = 8, height=1)



# combine and arrange plots
ggar = ggarrange(mts, nts, mtsleg, ntsleg, dots, cts, dotsleg, ctsleg,
   nrow=4,heights=c(4,2.5,4,2.5))
ggar
ggsave(plot=ggar, file="output/04_figure_ts_combined_meta.png",device=png(),width=10,height=8*1.615)


ggar2 = ggarrange(dts, pts, dtsleg, ptsleg,
  nrow=2, heights=c(4,1))
ggar2
ggsave(plot=ggar2, file="output/04_figure_ts_combined_clusters.png",device=png(),width=10,height=5*1.25)


#######

source("bin/05_maps.R")
