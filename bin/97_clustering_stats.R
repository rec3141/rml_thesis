

## heirarchical testing
#For reasons outlined in the corresponding paper (Kimes et al. 2017) relating to how the method handles testing when n << p, we recommmend using "euclidean" as the metric, and any of "ward.D2", "single", "average", "complete" as the linkage.
#If metric functions which do not statisfy rotation invariance are desired, e.g. one minus Pearson correlation ("cor") or L1 ("manhattan"), null_alg = "2means" and ci = "2CI" should be specified. The null_alg and ci parameters specify the algorithm for clustering and measure of "cluster strength" used to generate the null distribution for assessing significance. Since the K-means algorithm (2means) optimizes the 2-means CI (2CI), the resulting p-value will be conservative. However, since the hierarchical algorithm is not rotation invariant, using null_alg = "hclust" or ci = "linkage" produces unreliable results. An example for testing using Pearson correlation is given later in this section.

#shc_result <- shc(as.matrix(x.taxa)^.25, matmet=function(x) vegdist(x, method="bray", na.rm=TRUE), linkage="ward.D2", null_alg = "2means", ci = "2CI") #no significant clusters: warnings: results may be meaningless because data have negative entries in method “bray”
#shc_result <- shc(as.matrix(x.taxa), matmet=function(x) vegdist(x, method="bray", na.rm=TRUE), linkage="ward.D2") # gives many significant clusters
#shc_result <- shc(as.matrix(x.taxa), metric="euclidean", linkage="ward.D2") #gives significant clusters based on dominant taxon
#shc_result <- shc(as.matrix(x.taxa), matmet=function(x) vegdist(x+min(x), method="bray", na.rm=TRUE), linkage="ward.D2", null_alg = "2means", ci = "2CI") #no sig clusters, warnings: results may be meaningless because data have negative entries in method “bray”
#shc_result <- shc(as.matrix(x.taxa), matmet=function(x) vegdist(x-min(x), method="bray", na.rm=TRUE), linkage="ward.D2", null_alg = "2means", ci = "2CI")# no significant clusters, no warnings
#shc_result <- shc(as.matrix(x.taxa)^(1/4), matmet=function(x) vegdist(x-min(x), method="bray", na.rm=TRUE), linkage="ward.D2") #gives 6 significant clusters
#shc_result <- shc(as.matrix(x.taxa)^(1/3), matmet=function(x) vegdist(x-min(x), method="bray", na.rm=TRUE), linkage="ward.D2") #gives 5 significant clusters
#shc_result <- shc(as.matrix(x.taxa)^(1/2), matmet=function(x) vegdist(x-min(x), method="bray", na.rm=TRUE), linkage="ward.D2") #gives 4 significant clusters
#shc_result <- shc(as.matrix(x.taxa)^(1/1), matmet=function(x) vegdist(x-min(x), method="bray", na.rm=TRUE), linkage="ward.D2") #gives 4 significant clusters
#shc_result <- shc(0+as.matrix(x.taxa)>0, matmet=function(x) vegdist(x-min(x), method="bray", na.rm=TRUE), linkage="ward.D2") #gives 3 significant clusters
#shc_result <- shc(0+as.matrix(x.taxa)>0, matmet=function(x) vegdist(x-min(x), method="bray", binary=TRUE, na.rm=TRUE), linkage="ward.D2") #gives 40 significant clusters
#shc_result <- shc(log1p(as.matrix(x.taxa)), matmet=function(x) vegdist(x-min(x), method="bray", na.rm=TRUE), linkage="ward.D2") #gives 4 significant clusters
#shc_result <- shc(as.matrix(x.taxa)^(1/4), matmet=function(x) vegdist(x-min(x), method="bray", na.rm=TRUE), linkage="ward.D2", null_alg = "2means", ci = "2CI") #gives 0 significant clusters
#shc_result <- shc(as.matrix(x.taxa)^(1/4), matmet=function(x) vegdist(x-min(x), method="bray", na.rm=TRUE), ci=c("2CI", "linkage", "2CI"), null_alg=c("hclust", "hclust", "2means"), linkage="ward.D2") #gives 6 significant clusters
#plot(shc_result, hang=.1)

#pvc_result <- pvclust(as.matrix(x.taxa)^(1/4), method.dist=function(x) vegan::vegdist(x, method="bray", na.rm=TRUE), method.hclust="ward.D2", nboot=1000, parallel=TRUE) #took >5 minutes, bootstraps all 0...
pvc_result <- pvclust(as.matrix(x.taxa), parallel=TRUE, method.dist=function(x) dist(x), method.hclust="ward.D2", nboot=10, use.cor="all") #took >5 minutes, bootstraps all 0...
pvc_result <- pvclust(as.matrix(x.taxa), parallel=TRUE, method.dist=function(x) dist(x), method.hclust="ward.D2", nboot=10, use.cor="all") #took >5 minutes, bootstraps all 0...
plot(pvc_result)
pvrect(pvc_result, alpha=0.95)

library(pvclust)

mydist = function(x) {
  x[x<0] = 0; x[is.na(x)] = 0; x[is.infinite(x)] = 0; 
  x[rowSums(x)==0,] = mean(x);
  x[,colSums(x)==0] = mean(x);
  as.dist(vegan::vegdist(t(x),method="bray", na.rm=TRUE))
}

# mydist = function(x) dist(x,method="euclidean")
pvc_result <- pvclust((as.matrix(x.taxa[1:50,])^.25), method.dist=mydist, method.hclust="ward.D2", nboot=100, parallel=TRUE)
plot(pvc_result)

pvc_result <- pvclust(t(as.matrix(x.taxa)^.25), method.dist="euclidean", method.hclust="average", nboot=100, parallel=TRUE)
plot(pvc_result)




#significance using SIMPROF
# can't figure out how to get heights out of trees in order to use cutree
hm_sig = simprof(as.matrix(x.prop[1:50,])^.25, num.expected=100, num.simulated=99, silent=FALSE, 
                 increment = 100, method.cluster = "ward.D2", 
                 method.distance = function(x) {x[x<0] = 0; x[is.na(x)] = 0; x[is.infinite(x)] = 0; x[rowSums(x)==0,] = mean(x); vegdist(x, method="bray", na.rm=TRUE)})

hm_col = simprof.plot(hm_sig,leafcolors=gg_color_hue(hm_sig$numgroups)[sample(1:hm_sig$numgroups,hm_sig$numgroups)])

hm_phy = as.phylo(hm_col)

get_branches_heights(hm_col,include_leaves = TRUE)[  as.numeric(unlist(lapply(hm_sig$significantclusters, function(x) try({findMRCA(as.phylo(hm_col),tips=x, type="node")}))))]

maxht = max(
  as.numeric(unlist(lapply(hm_sig$significantclusters, function(x) try({findMRCA(as.phylo(hm_col),tips=x, type="height")}))))
  ,na.rm=TRUE)




###


pv_clust = pvclust(t((as.matrix(x.taxa)^.25)), method.dist=mydist, method.hclust="ward.D2", nboot=100, parallel=TRUE)
plot(pv_clust)

(p <- ggplot(mse, aes(x = n.samp, y = means, group = group)) +
    geom_errorbar(aes(ymax = upper.ci, ymin = lower.ci), width = 0.2) +
    geom_point(aes(shape = group, fill = group), size = 4) + 
    scale_shape_manual(values = c(21, 24:25), name = "") +
    scale_fill_manual(values = c("red", "blue", "chartreuse3"), name = "") + 
    coord_cartesian(ylim = c(0, 0.65)) +
    theme_bw(base_size = 18) +
    labs(x = "Sample size (n)", y = "Multivariate pseudo SE") +
    theme(legend.position = c(0.8, 0.8), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) )

# hm_cut = dendextend::cutree(hm$rowDendrogram, h=maxht)           # k is number of clusters, can also do by h, height, to cut tree

