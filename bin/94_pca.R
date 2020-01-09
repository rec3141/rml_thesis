library(scales)
wrap_strings = function(vector_of_strings,width){as.character(sapply(vector_of_strings,FUN=function(x){paste(strwrap(x,width=width), collapse="\n")}))}

metavar = c("depth_m","lat","lon","temp","salinity","DO",
            "FlECO-AFL(mg/m^3)","chl (ug/l)","phaeo (ug/l)","PO4(uM)","Sil(uM)",
            "NH4(uM)","N+N (umol/L)","transmission_pct")
metavarlog = c("FlECO-AFL(mg/m^3)","chl (ug/l)","phaeo (ug/l)","PO4(uM)","Sil(uM)",
               "NO3(uM)","NO2(uM)","NH4(uM)","N+N (umol/L)")
metaplot = c("lat","lon","temp","salinity","DO",
             "FlECO-AFL(mg/m^3)","PO4(uM)","Sil(uM)",
             "N+N (umol/L)","NH4(uM)")

# metaplot = c("temp","salinity","chl (ug/l)","N+N (umol/L)","DO")

metasmall = d.meta
# metasmall$`chl (ug/l)` = log1p(metasmall$`chl (ug/l)`)
# metasmall$`N+N (umol/L)` = log1p(metasmall$`N+N (umol/L)`)
metasmall[,metavarlog] = log1p(metasmall[,metavarlog])
metasmall = metasmall[,metaplot]
metasmall = unique(na.omit(metasmall))
metapca = prcomp(metasmall,scale=T,center=T)$x
metapca[,2] = metapca[,2]*-1
metapca[,1] = metapca[,1]*-1

df = data.frame(d.meta[rownames(metasmall),],metapca)
ggplot(data=df, aes(x=PC1, y=PC2)) + geom_point(aes(col=cluster, shape=project))

ggplot(data=df, aes(x=PC1, y=PC2)) + geom_point(aes(col=mass, shape=project))

plot(metapca[,1],metapca[,2],bg=d.meta[rownames(metasmall),"mass"],pch=21)

#pca_samples = intersect(rownames(metapca),rownames(taxa))
propsmall = taxa[rownames(metapca),]

propsmall = propsmall[rowSums(propsmall)>0,colSums(propsmall)>0]
#propsmall = propsmall[,colMeans(propsmall)>0.0000]
propsmall = propsmall[,colSums(propsmall>0)>=20]

propwtx = apply(propsmall,2,function(y) weighted.mean(x=metapca[,1],w=y))
propsmall = propsmall[,order(propwtx)]
propwtx = propwtx[order(propwtx)]

propwty = apply(propsmall,2,function(y) weighted.mean(x=metapca[,2],w=y))
propwty = propwty[order(propwtx)]

pdf(file="output/figure_pca_plots.pdf")
plot(propwtx,propwty,pch=21,col=rgb(0,0,1,0.5), bg=NULL,cex=scales::rescale(sqrt(colSums(propsmall)),to=c(0.3,3)))

for(i in 1:ncol(metasmall)) {
  plot(metapca[,1],metapca[,2], cex=scales::rescale(metasmall[,i],to=c(0,3)), main=colnames(metasmall)[i], pch=21, bg=rgb(0,0,1,0.5), col=NULL, xlim=c(-3,5),ylim=c(-6,3))
  
}

for(i in 1:ncol(propsmall)) { 
  plot(metapca[,1],metapca[,2], cex=0.3, cex.main=0.5, main=colnames(propsmall)[i], pch=21, bg=rgb(0,0,0,0.5), col=NULL, xlim=c(-3,5),ylim=c(-6,3))
  points(metapca[,1],metapca[,2], cex=scales::rescale(sqrt(propsmall[,i]),to=c(0,3)), pch=21, bg=rgb(0,0,1,0.5), col=NULL, )
  points(propwtx[i],-2,bg="red",pch=21,col=NULL)
  points(-3, propwty[i],bg="red",pch=21,col=NULL)
}
dev.off()
