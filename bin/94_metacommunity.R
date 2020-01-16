mat.clust = list()
for(i in names(mat.list)) {
  mat.clust[[i]] = mat.list[[i]][["meta"]][,c("cluster","cellid")]
#  print(mat.list[[i]][["meta"]]["BOX_9_77",])
  }

#mat.mrg = do.call(merge, mat.clust,args=list("by"="cellid"))
mat.mrg = do.call(rbind, mat.clust)
mat.agg = aggregate(mat.mrg, by = list(mat.mrg$cellid),FUN=length)
mat.agg = reshape2::dcast(mat.mrg,cellid~cluster)
rownames(mat.agg) = mat.agg$cellid
mat.agg$cellid = NULL
mat.agg[!is.na(mat.agg)] = 1
mat.agg[is.na(mat.agg)] = 0
mat.agg = as.matrix(mat.agg)
class(mat.agg) = "numeric"
mat.agg = mat.agg[rownames(sub.meta),] #reorder to match

# rownames(mat.agg) = sub.meta$station[match(rownames(mat.agg),sub.meta$cellid)]
pdf(file="output/94_figure_heatmap_meta_clusters.pdf", width=12, height=12)
hmmeta = heatmap(mat.agg, scale="none", distfun = function(x) vegdist(x, method="bray", binary=TRUE), hclustfun = function(x) hclust(x, method = "ward.D2"),keep.dendro = TRUE)
dev.off()

heatmap(cor(mat.agg), distfun = function(x) vegdist(x, method="bray", binary=TRUE), hclustfun = function(x) hclust(x, method = "ward.D2"))

nclus=6

hmmeta_cut = dendextend::cutree(hmmeta$Rowv, k=nclus)           # k is number of clusters, can also do by h, height, to cut tree
hmm_uniqclus = unique(hmmeta_cut[hmmeta$rowInd]) #need to have it in order
hmm_clusnames = paste0("M",1:nclus)

#replace levels with cluster names, in tree order
for(i in 1:nclus) {
  j = hmm_uniqclus[i]
  hmmeta_cut[hmmeta_cut==j] = hmm_clusnames[i]
}
sub.meta$metacluster <- factor(hmmeta_cut,levels = hmm_clusnames)

tre = as.dendrogram(hmmeta$Rowv)
hts = heights_per_k.dendrogram(dend = tre)
tre.upper = cut(tre,h = hts[nclus])$upper
tre.upper = force.ultrametric(as.phylo(tre.upper))
tre.upper$tip.label = as.character(levels(x.meta$cluster))
tre.rle = rle(sort(as.character(x.meta$cluster)))
tre.rle = data.frame(tre.rle$values,"N"=tre.rle$lengths)
rownames(tre.rle) = tre.rle[,1]
tre.rle = tre.rle[,-1,drop=F]

# heatmap with rowside color bars
rowcols <- cbind(
  gg_color_hue(length(hmmeta_cut))[order(hmmeta$rowInd)],
  rep("#ffffff",nrow(sub.meta)),
  gg_color_hue(nclus)[as.numeric(factor(sub.meta$metacluster))],
  rep("#ffffff",nrow(sub.meta)),
  gg_color_hue(length(unique(sub.meta$filter)))[as.numeric(factor(sub.meta$filter))],
  rep("#ffffff",nrow(sub.meta)),
  gg_color_hue(length(unique(sub.meta$mass)))[as.numeric(factor(sub.meta$mass))],
  rep("#ffffff",nrow(sub.meta)),
  gg_color_hue(length(unique(sub.meta$depth_type)))[as.numeric(factor(sub.meta$depth_type))]
)

#plot heatmap with bars
pdf(file=paste0('output/94_figure_heatmap_metacluster_bars.pdf'), width=48, height=48)
hmmeta_cuts2 = heatmap.plus(mat.agg, margins=c(12,12),
                         RowSideColors = rowcols,
                         labRow = paste0(sub.meta$depth_type,"_",sub.meta$station),
                         Rowv = hmmeta$Rowv,
                         Colv = hmmeta$Colv)
dev.off()
sub.meta$depth_type = factor(sub.meta$depth_type,levels=c("surface","midwater","bottom"),ordered=T)
dmm = basemap +
  geom_point(data=sub.meta, aes(x=lon, y=lat, color=metacluster, shape=project), size=2) +
  # geom_text(data=x.meta, aes(x= lon, y = lat, label=station), hjust=1, vjust=-.5, size = 1, check_overlap=TRUE) +
  theme(legend.position="right") + theme(panel.spacing = unit(2, "lines"), strip.text = element_text(size=14)) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.text.align=0, legend.box.just = "left") +
  guides(col=guide_legend(ncol=1,byrow=TRUE)) + ggtitle(tx) +
  facet_wrap(vars(project,depth_type))

ggsave(plot=dmm, filename=paste0("output/94_figure_map_metacluster_clusters.png"),width=14,height=18,device=png())



mat.pca = prcomp(mat.agg)
plot(mat.pca$x[,1],mat.pca$x[,2],pch=21,col=NULL,bg=gg_color_hue(length(hmmeta_cut))[order(hmmeta$rowInd)])
#gg_color_hue(nclus)[as.numeric(factor(hmmeta_cut))]

library(Rtsne)
mat.tsne = Rtsne(mat.pca$x,check_duplicates = F)
plot(mat.tsne$Y,pch=21,bg=gg_color_hue(length(hmmeta_cut))[order(hmmeta$rowInd)],col=NULL)

