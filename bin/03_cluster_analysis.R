## Cluster analysis
## dendrogram plus bar plot (big bar plot of all samples and taxonomic make-up with dendrogram)

# require(metafolio)
# library(ggrepel)
require(gplots)
require(vegan)
require(dendextend)    # for cutree
require(phytools)
require(heatmap.plus)
require(ggtree)
require(ggstance)
require(colorspace)
require(ggpubr)
require(compositions)
require(grid)
require(egg) # for ggarrange(), overwrites ggpubr::ggarrange
require(ggbeeswarm)
require(clustsig)

clus.list = list("diatoms" = 7,
                 "picos" = 11,
                 "ciliates" = 13,
                 "fungi" = 21,
                 "dinos" = 15,
                 "other" = 7,
                 "bacteria" = 13,
                 "metazoa" = 7)

cutoff = 0.85

for(tax in names(mat.list)) {

  tx = mat.list[[tax]][["title"]]
  abb = mat.list[[tax]][["short"]]
  x.meta = mat.list[[tax]][["meta"]]
  x.prop = mat.list[[tax]][["prop"]]
  nclus = clus.list[[tax]]
  clusnames = paste0(abb,1:nclus)


#heatmap and clustering

## Using fourth root transform and Bray distance
par(mfrow=c(2,2))
hist(as.matrix(x.taxa[x.taxa>0])^1,breaks=100, main="^1")
hist(as.matrix(x.taxa[x.taxa>0])^.5,breaks=100, main="^1/2")
hist(as.matrix(x.taxa[x.taxa>0])^.25,breaks=100, main="^1/4")
hist(log(as.matrix(x.taxa[x.taxa>0])),breaks=100, main="log")

pdf(file=paste0('output/03_figure_heatmap_',tax,'_manbray.pdf'), width=24, height=24)
hmman <- heatmap.2(as.matrix(x.prop)^.25, trace='none', margins=c(5,5),
                labRow = x.meta$project,
                distfun=function(x) vegdist(x, method="bray", na.rm=TRUE),
                hclustfun=function(x) hclust(x, method="ward.D2"))
dev.off()

## Using CLR transform and Euclidean distance
pdf(file=paste0('output/03_figure_heatmap_',tax,'_clreuc.pdf'), width=24, height=24)
hmclr <- heatmap.2(as.matrix(clr(x.prop)), trace='none', margins=c(5,5),
                labRow = x.meta$project,
                distfun=function(x) vegdist(x, method="euc", na.rm=TRUE), 
                hclustfun=function(x) hclust(x, method="ward.D2"))
dev.off()


## Using CLR transformed dendrograms with Manhattan transformed values
pdf(file=paste0('output/03_figure_heatmap_',tax,'_manbray_with_clr.pdf'), width=24, height=24)
hmmanclr <- heatmap.2(as.matrix(x.prop)^.25, trace='none', margins=c(5,5),
                   labRow = x.meta$project,
                   distfun=function(x) vegdist(x, method="bray", na.rm=TRUE),
                   Rowv = hmclr$rowDendrogram, 
                   Colv = hmclr$colDendrogram)
dev.off()


# the CLR jumbles the data, tending to merge samples that have just one or two major taxa
# the Manhattan transform preserves the structure better, so going forward with it
hm = hmman

hm_cut = dendextend::cutree(hm$rowDendrogram, k=nclus)           # k is number of clusters, can also do by h, height, to cut tree
uniqclus = unique(hm_cut[hm$rowInd]) #need to have it in order

#rename clusters by tree order
for(i in 1:nclus) {
  j = uniqclus[i]
  hm_cut[hm_cut==j] = clusnames[i]
}
x.meta$cluster <- factor(hm_cut,levels = clusnames)

tre = as.dendrogram(hm$rowDendrogram)
hts = heights_per_k.dendrogram(dend = tre)
tre.upper = cut(tre,h = hts[nclus])$upper
tre.upper = force.ultrametric(as.phylo(tre.upper))
tre.upper$tip.label = as.character(levels(x.meta$cluster))
tre.rle = rle(sort(as.character(x.meta$cluster)))
tre.rle = data.frame(tre.rle$values,"N"=tre.rle$lengths)
rownames(tre.rle) = tre.rle[,1]
tre.rle = tre.rle[,-1,drop=F]
tre.rle = tre.rle[order(as.numeric(gsub(abb,"", rownames(tre.rle)))),,drop=F]

# heatmap with rowside color bars
rowcols <- cbind(
  gg_color_hue(length(unique(x.meta$cluster)))[as.numeric(factor(x.meta$cluster))],
  rep("#ffffff",nrow(x.meta)),
  gg_color_hue(length(unique(x.meta$filter)))[as.numeric(factor(x.meta$filter))],
  rep("#ffffff",nrow(x.meta)),
  gg_color_hue(length(unique(x.meta$mass)))[as.numeric(factor(x.meta$mass))],
  rep("#ffffff",nrow(x.meta)),
  gg_color_hue(length(unique(x.meta$depth_type)))[as.numeric(factor(x.meta$depth_type))]
)

#plot heatmap with bars
pdf(file=paste0('output/03_figure_heatmap_',tax,'_bars.pdf'), width=48, height=48)
hm_cuts2 <- heatmap.plus(as.matrix(x.prop)^.25, margins=c(12,12),
                         RowSideColors = rowcols,
                         labRow = x.meta$project,  
                         distfun=function(x) vegdist(x, method="bray", na.rm=TRUE),
                         hclustfun=function(x) hclust(x, method="ward.D2"))
dev.off()

plot(tre.upper)
nodelabels()
tiplabels()

# construct ggtree plot with upper tree and barplots

#rotate to keep tree order
tre.upper = rotateConstr(tre.upper, clusnames)

treeplot <- ggtree(tre.upper, ladderize=F) + geom_tiplab(aes(angle=90,size=24),offset=.5,hjust=0.5) + xlim_tree(3.7)

#have to flip nodes on picos tree to get them in order
#if(tax=="picos") { treeplot = ggtree::flip(treeplot, 11,12) }

tredat <- data.frame(id=tre.upper$tip.label, tre.rle, col=factor(as.character(1:nclus),levels=1:nclus), stringsAsFactors = F)

treebarplot <- facet_plot(data=tredat, treeplot, panel="barplot", geom=geom_barh, mapping = aes(x = N, fill=col), stat='identity') + 
  theme_tree2() + theme(panel.spacing = unit(0.5, "lines")) + theme(legend.position = "none") +
  geom_text(data=tredat,
            aes(x=c(rep(NA,nclus),rep(25,nclus)),y=c(rep(NA,nclus),1:nclus),label=N,angle=90),
            vjust=-0.5,na.rm = T) +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),axis.line=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank())

#plot rotated tree
png(file=paste0("output/03_figure_tree_",tax,".png"),width=6,height=6,units = "in", res=300)
print(treebarplot, vp=viewport(angle=-90, width = unit(.5, "npc"), height = unit(1, "npc")))
dev.off()
# ggsave(treebarplot, filename=paste0("output/03_figure_tree_",tax,".png"),width=6,height=6)





## plot figure 4 -- taxonomic barplots by cluster and depth

# pick the cumulative top % of ESVs from each cluster
x.pick = NULL
for(i in 1:length(unique(x.meta$cluster))) {
    x.here = x.prop[x.meta$cluster==unique(x.meta$cluster)[i],]
    x.here = x.here[,order(colSums(x.here),decreasing=TRUE)]
    
    x.this = which( cumsum(colSums(x.here)/sum(x.here)) < cutoff)
    x.pick = union(x.pick,names(x.this))
}

x.pick = x.pick[order(colSums(x.prop[,x.pick]),decreasing=TRUE)]
x.pick.prop <- x.prop[, x.pick]
x.other.prop <- x.prop[, setdiff(colnames(x.prop),x.pick)]

# for sorting by taxonomy to make nice colors
tmptax = taxa_split[match(x.pick, taxa_split$ESV),]
rownames(tmptax) = x.pick
tmptax$abund = colMeans(x.pick.prop)
tmptax$esvrank = as.numeric(gsub("ESV_","",tmptax$ESV))

tmptax = tmptax[
  order(tmptax$Domain, tmptax$Kingdom, tmptax$Phylum, tmptax$Class, tmptax$Family, tmptax$Genus, tmptax$abund, tmptax$esvrank,
  decreasing=c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE),method="radix")
  ,]

othername = paste0("Other (", sprintf("%.1f",mean(rowSums(x.other.prop)*1000),1),")")
tmptax[othername,] = c(rep("Other",times=ncol(taxa_split)), 0, 0)

#names to plot with, including abundances

x.pick.prop[,othername] = rowSums(x.other.prop)
x.pick.prop = x.pick.prop[,rownames(tmptax)]

colnames(x.pick.prop) = paste0(
  barplot_names[match(colnames(x.pick.prop),taxa_split$ESV)], " (", 
  sprintf("%.1f",round(colMeans(x.pick.prop)*1000,1)), ")")

colnames(x.pick.prop)[ncol(x.pick.prop)] = othername

x.meta$sample = rownames(x.meta)

# reshape data frame for ggplot
df <- data.frame(x.pick.prop,x.meta[,c('depth_type','station','project','cluster')],check.names=FALSE) #taxa

df_long <- reshape2::melt(df, id.vars=c('depth_type','station','project','cluster'), variable.names='taxa', na.rm=TRUE)
df_long$facet1 <- factor(df_long$depth_type, levels=c('bottom', 'midwater', 'surface'), labels=c("bottom","midwater","surface"), ordered=TRUE)
df_long$facet2 <- factor(df_long$cluster, levels=clusnames)

newdat <- df_long[!is.na(df_long$facet2), ]

# construct barplot
bplot <- ggplot(newdat, aes(fill = variable, y = value, x = facet1, label=variable)) +
  geom_bar(stat='identity', position=position_fill(reverse=T)) +
  ylab('Relative Abundance') + xlab('Depth Bin') +
  theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14)) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y = element_text(size=12)) + coord_flip() +
  theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12), strip.text = element_text(size = 12)) + 
  ggtitle(paste0(tx,' Clusters')) + 
  facet_wrap(vars(facet2), scales='fixed', nrow=1) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_blank()) +
  theme(legend.text=element_text(size=8), legend.title=element_blank(), legend.position = 'none')

# this colors the bars and legend by taxonomic level (e.g. 'Genus', 'Family')
bplot = bplot + scale_fill_manual("legend",values = getcolors(tmptax,"Genus"))
#bplot
ggsave(plot = bplot, filename = paste0('output/03_figure_barplot_',tax,'.png'),device = png(),width = 12, height=4)

ppl = bplot + theme(legend.position="bottom") 
ppl_legend <- ggpubr::get_legend(ppl)
ppleg = ggpubr::as_ggplot(ppl_legend)
ppleg
ggsave(ppleg, filename = paste0('output/03_figure_barplot_',tax,'_leg.png'),device = png(),width = 18, height=8,dpi = 300)


#### Generate Table 4
table4 = data.frame(x.meta[,c('depth_type','project','cluster')],value=rep(1,nrow(x.meta)))
table4 = reshape2::dcast(table4,project + depth_type ~ cluster)
table4 = table4[order(table4$depth_type, decreasing=T),]
table4 = table4[order(table4$project),]
colnames(table4)[1:2] = c("Project", "Depth Bin")
write.table(table4,file=paste0("output/03_table_clusters_",tax,".tsv"),sep="\t",quote=F,row.names=F)

#### box plots

metavar = c("depth_m","lat","lon","temp","salinity","DO",
            "FlECO-AFL(mg/m^3)","chl (ug/l)","phaeo (ug/l)","PO4(uM)","Sil(uM)",
            "NH4(uM)","N+N (umol/L)","transmission_pct")
metavarlog = c("FlECO-AFL(mg/m^3)","chl (ug/l)","phaeo (ug/l)","PO4(uM)","Sil(uM)",
               "NO3(uM)","NO2(uM)","NH4(uM)","N+N (umol/L)")
metaplot = c("lat","lon","temp","salinity","DO","transmission_pct",
            "FlECO-AFL(mg/m^3)","PO4(uM)","Sil(uM)",
            "N+N (umol/L)","NH4(uM)","mass","cluster")

#AMBON doesn't have NO3, only total
#DBO doesn't have transmission values

z.meta = x.meta
z.meta[,metavarlog] = log10(z.meta[,metavarlog])
z.meta = z.meta[,metaplot]
colnames(z.meta) = make.names(colnames(z.meta))
znames = c("Latitude\n(°N)","Longitude\n(°E)","Temperature\n(°C)","Salinity\n(PSU)","Dissolved Oxygen\n(\u00b5mol/kg)",
                     "Beam Transmission\n(%)", "Chl a Fluoresence\n(mg/m\u00B3)",
                     "Phosphate\n(\u00b5M)", "Silicate\n(\u00b5M)","Total Nitrate\n(\u00b5M)",
                     "Ammonium\n(\u00b5M)","Water Mass","Cluster")

#plot each metadata variable as a function of the clusters
bptmp = list()
for(i in 1:ncol(z.meta)) {
  mvar = colnames(z.meta)[i]
  if(mvar=="cluster") next
  bptmp[[i]] = ggplot(data=z.meta, aes_string(x="cluster", y=mvar, na.rm=TRUE)) + 
    # geom_violin() + geom_beeswarm()
  geom_boxplot(outlier.shape = NA,) + 
    #geom_jitter(width = 0.3,height = 0,alpha=0.6) +
  geom_quasirandom(width=0.3,bandwidth = 1, aes_string(color="mass"), alpha=1,na.rm = TRUE) +
      theme(axis.text.x=element_text(size = 10),
  axis.text.y=element_text(size = 10),
  axis.title.x=element_text(size = 10),
  axis.title.y=element_text(size = 10)) +
    labs(x=NULL,y=znames[i]) +
    scale_color_discrete(name="Water Mass",drop = FALSE)

  }

pdf(file=paste0("output/03_figure_cluster_boxplots_",tax,".pdf"),width=8,height=3)
invisible(lapply(bptmp, print))
dev.off()

bpar = ggarrange(plots=bptmp,ncol=1)
ggsave(plot=bpar, file=paste0("output/03_figure_cluster_boxplots_combined_",tax,".png"),
       device=png(),width=6,height=30)


### update saved files
mat.list[[tax]][["meta"]] = x.meta
mat.list[[tax]][["top.names"]] = rownames(tmptax)[1:(nrow(tmptax)-1)]
mat.list[[tax]][["top.longnames"]] = colnames(x.pick.prop)



graphics.off()

} #end taxonomy loop


#####


source("bin/04_ts_plots.R")
