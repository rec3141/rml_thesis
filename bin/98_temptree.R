library(Biostrings)
library(msa)
library(phangorn)
library(tidyverse)
library(ggtree)
#DECIPHER?
setwd("~/Downloads/thesis-master")

source("prep_data.R")

cor2pvalue = function(r, n) {
  t <- (r*sqrt(n-2))/sqrt(1-r^2)
  p <- 2*(1 - pt(abs(t),(n-2)))
  se <- sqrt((1-r*r)/(n-2))
  out <- list(r, n, t, p, se)
  names(out) <- c("r", "n", "t", "p", "se")
  return(out)
}


#prep_data.r first half

d.class = sort(unique(taxa_split[,4]))
#cl = "Bacillariophyta"

for(cl in d.class) {
  pdf(file=paste0("hm_class_",cl,".pdf"),width=8,height=8)
  try({
    
  
d.pick = grepl(pattern="Mediophyceae",colnames(taxa))
#d.pick = taxa_split[,4] == cl
if(sum(d.pick)<10) next

print(cl)
d.taxa = taxa[,d.pick]
s.pick = rowSums(d.taxa>0)>0 & meta$`station_type`=="S" & meta$filter==0.2
d.taxa = d.taxa[s.pick,]
d.taxa[d.taxa==0] = NA #leaves zeros out of correlations
d.meta = meta[s.pick,]

#to use CCA in heatmap
# taxa.cca <- cca(d.taxa ~ temp + salinity + DO + `FlECO-AFL(mg/m^3)` + `Sil(uM)` + filter + project + mass, data=d.meta, na.action = na.exclude)
# df = taxa.cca$CCA$v
# df = apply(df,2,function(x) rescale(x,to=c(-1,1)))

#do raw correlations
metavars = c("depth_m","lat","lon","temp","salinity","DO",
             "POC (ug/L)","PON (ug/L)","SPM (ug/L)",
             "transmission_pct","poc_pon_ratio")
metavarlog = c("FlECO-AFL(mg/m^3)","chl (ug/l)","phaeo (ug/l)","PO4(uM)","Sil(uM)",
               "NH4(uM)","N+N (umol/L)")

d.meta = cbind(d.meta[,metavars],log1p(d.meta[,metavarlog]))

# find (pairwise complete) correlation matrix between two matrices x and y
d.cor = cor(x=(d.taxa),y=d.meta,use="pairwise",method="pearson")

# df = apply(d.cor,2,function(x) rescale(x,to=c(-1,1)))
df = d.cor
dfc = df[rowSums(is.na(df))<ncol(df),]
hh = hclust(dist(t(dfc)))
df = df[,hh$order]

#get pvalues
d.n <- t(!is.na(d.taxa)) %*% (!is.na(d.meta)) # same as count.pairwise(x,y) from psych package

# get a list with matrices of correlation, pvalues, standard error, etc.
pv = cor2pvalue(d.cor,d.n)$p

# remove low-N pairs
df[d.n < 5] = NA
# to filter by pvalue
df[pv > (0.05/nrow(pv)/ncol(pv))] = NA
#df[pv > 1e-4] = NA
df = df[rowSums(is.na(df))<ncol(df),]

dd.pick = colnames(taxa) %in% rownames(df)
dd.taxa = taxa[,dd.pick]

badseqs = TRUE
i = 1
while(badseqs) {
  seqs.ss <- DNAStringSet(taxa_seqs[dd.pick])
  #Muscle is crash
  sub.aln <- msa(seqs.ss,type="dna","Muscle",order="input",verbose=T)
  #sub.aln <- msa(seqs.ss,type="dna","ClustalOmega",order="input",verbose=T)
  
  sub.aln.c <- as.character(sub.aln@unmasked)
  sub.aln.raw <- t(sapply(strsplit(sub.aln.c,""), tolower))
  rownames(sub.aln.raw) <- colnames(taxa)[dd.pick]
  sub.bin <- as.DNAbin(sub.aln.raw)
  
  # get distances
  sub.dist <- dist.dna(sub.bin,as.matrix=T,model="TN93")
  #sub.max = quantile(sub.dist,prob = seq(0, 1, length = 101))[99]
  sub.max = max(sub.dist,na.rm=T)
  # sub.dist[is.na(sub.dist)] <- sub.max
  # sub.dist[is.infinite(sub.dist)] <- sub.max
  # sub.dist[is.nan(sub.dist)] <- sub.max
  #sub.dist[sub.dist > sub.max] = sub.max
  
  #to remove bad sequences
  if(max(sub.dist,na.rm=T) > 0.6) {
    rem = which.max(colMeans(sub.dist))
    dd.pick[which(colnames(taxa) %in% names(rem))] = FALSE
  } else { badseqs = FALSE }
  
  # make tree
  sub.tree <- bionj(sub.dist)
  sub.tree <- midpoint(sub.tree)
   plot(sub.tree,main=i,show.tip.label=T)

  sub.out = sub.aln.c
  names(sub.out) = rownames(df)
   write.table(sub.out,file="aln.txt",quote=F,row.names=T)
   
    i = i+1
}

#df.tree = prune_taxa(taxa=common_names,sub.tree)
common_names = intersect(rownames(df),sub.tree$tip.label)
df = df[common_names,]

gt = ggtree(sub.tree)

gh = gheatmap(gt, df,width = 1, offset=.15, color = NULL,font.size=1.5,colnames_angle=45, colnames_offset_x = 0, colnames_offset_y = 0, hjust=0, colnames_position = "top") + 
  geom_tiplab(size=1.5,align=TRUE, linesize=.25) +
  scale_fill_viridis_c(limits=range(d.cor,na.rm=T),na.value = "darkgrey")# +
# theme_classic()
#    scale_fill_gradient2(low = "blue",high="green",mid = "black",limits=c(-1,1),na.value = "black")

print(gh)

  
})

dev.off()
  
}

# 
# 
# 
# # Generate a random tree with 30 tips
# #tree <- rtree(30)
# tree = sub.tree
# # Make the original plot
# p <- ggtree(tree)
# 
# # generate some random values for each tip label in the data
# d1 <- data.frame(id=tree$tip.label, val=rnorm(190, sd=3))
# 
# # Make a second plot with the original, naming the new plot "dot",
# # using the data you just created, with a point geom.
# p2 <- facet_plot(p, panel="dot", data=d1, geom=geom_point, aes(x=val), color='red3')
# 
# # Make some more data with another random value.
# d2 <- data.frame(id=tree$tip.label, value = abs(rnorm(190, mean=100, sd=50)))
# 
# # Now add to that second plot, this time using the new d2 data above,
# # This time showing a bar segment, size 3, colored blue.
# p3 <- facet_plot(p2, panel='bar', data=d2, geom=geom_segment,
#                  aes(x=0, xend=value, y=y, yend=y), size=3, color='blue4')
# 
# # Show all three plots with a scale
# p3 + theme_tree2()
# 

