metad <- readRDS('meta_diatom_esv.RDS')
metap <- readRDS('meta_pico_esv.RDS')
# File options:
## meta_diatom_esv.RDS, meta_diatom_genus.RDS, meta_diatom_fam.RDS
## meta_pico_esv.RDS, meta_pico_genus.RDS, meta_pico_family.RDS

taxad <- readRDS('diatoms.RDS')
# diatoms.RDS, picoeukaryotes.RDS


original_names <- colnames(taxad)
ord_orig =  as.data.frame(read.table(text=colnames(taxad), sep=";", as.is=TRUE, fill=TRUE)$V13)
fam <- as.data.frame(make.unique(read.table(text=colnames(taxad), sep=";", as.is=TRUE, fill=TRUE)$V20))
ge <- as.data.frame(make.unique(read.table(text=colnames(taxad), sep=";", as.is=TRUE, fill=TRUE)$V22))
sp <- as.data.frame(make.unique(read.table(text=colnames(taxad), sep=";", as.is=TRUE, fill=TRUE)$V23))
esv <- as.data.frame(make.unique(read.table(text=colnames(taxad), sep=";", as.is=TRUE, fill=TRUE)$V24))
names <- as.data.frame(cbind(fam,ge,sp,esv))
colnames(names) <- c('family','genus','accession','esv')
names_short <- paste(names$esv)#, names$genus, names$accession, names$esv, sep=';')
colnames(taxad) <- names_short

dot <- gsub("\\..*", "", colnames(taxad))
colnames(taxad) <- dot

#df <- t(rowsum(t(taxad), group = rownames(t(taxad))))
df = taxad
#taxa_transform <- as.data.frame(clr(df))
taxa_d = df[which(rowSums(df>0)>0),which(colSums(df>0)>0)]

common <- intersect(rownames(taxa_d), rownames(metad))
metad <- metad[common, ]
taxa_d <- taxa_d[common, ]

c <- metad[metad$temp<=6, ]  # temperature cutoff 
w <- metad[metad$temp>=6, ]

## for ESV 2
esv <- as.data.frame(taxa_d[colnames(taxa_d)=='ESV_2'])
colnames(esv) <- 'ESV_2'
rownames(esv) <- rownames(taxa_d)
commonc <- intersect(rownames(esv), rownames(c))
commond <- intersect(rownames(esv), rownames(w))
esvc <- esv[commonc, ]
esvw <- esv[commond, ]
t.test(esvc, esvw, paired = FALSE, alternative = "two.sided")


### REC
library(plotrix)
taxad = taxa

pdf(file="mintemp.pdf",width=8,height=16)
mat.save = matrix(data=0,nrow=ncol(taxad),ncol=length(seq(-2,13,1))-1)
for(j in 1:ncol(taxad)) {
esv = as.data.frame(taxad[,j])
colnames(esv)= colnames(taxad)[j]
rownames(esv) = rownames(taxad)
print(j)
tt = NA*(1:length(meta$temp))
temps = sort(meta$temp,na.last = TRUE)
for(i in seq_along(temps)) {
  temp = temps[i]
  if(is.na(temp)) next
    try({
    cm = meta[meta$temp <= temp,"temp",drop=F]
wm = meta[meta$temp >= temp,"temp",drop=F]
commonc <- intersect(rownames(esv), rownames(cm))
commond <- intersect(rownames(esv), rownames(wm))
esvc <- esv[commonc, ]
esvw <- esv[commond, ]
  tt[i] = t.test(esvc, esvw, paired = FALSE, alternative = "two.sided")$p.value
},silent = T)

  }


# 
# par(mfrow=c(2,1))
# try({
#   plot(temps,log10(tt),type='l',xlab="temperature",ylab="log(p.value)",main=colnames(esv),xlim=c(-2,13))
# #  hist(as.numeric(c(cm$temp,wm$temp)),breaks=(temps/10-2))
# #  hist(as.numeric(c(cm$temp,wm$temp)),breaks=(temps/10-2))
#   
#   plot(meta$temp,esv[,1],xlim=c(-2,13))
#     },silent = T)

}

library(viridis)
library(gplots)
library(heatmap.plus)

taxa.ord = order(colSums(taxa),decreasing=T)

taxad = taxa[,taxa.ord]
#take top 95% of relabund
pick.taxa = unname(cumsum((colSums(taxad))/sum(colSums(taxad)))<0.98)
#at least 1% relabund in 2+ samples
pick.taxa = unname(colSums(taxad>1e4)>1)
taxad = taxad[,pick.taxa]

king = taxa_split[taxa.ord,"Class"]
king = king[pick.taxa]
#king = factor(king)
clas = taxa_split[taxa.ord,"Order"]
clas = clas[pick.taxa]

#take all those with at least 5% occurrence in any sample
# taxad = taxad[,which()]
mat.save = matrix(data=0,nrow=ncol(taxad),ncol=length(seq(-2,13,1))-1)
#mat.save = data.frame(stringsAsFactors=F)
for(j in 1:ncol(taxad)) {
  print(j)
    esv = as.data.frame(taxad[,j])
  colnames(esv)= colnames(taxad)[j]
  rownames(esv) = rownames(taxad)

  mat = cbind(meta$temp,esv[,1])
  mat = mat[complete.cases(mat),]
  
  wh = weighted.hist(mat[,1],plotrix::rescale(mat[,2],c(0,1)),breaks=seq(-2,13,1),freq=T,plot = F)  
#  mat.save = rbind(mat.save,cbind(colnames(esv),wh$breaks[1:length(wh$counts)],wh$counts))
  mat.save[j,] = wh$counts
  
    }

rownames(mat.save) = colnames(taxad)
colnames(mat.save) = paste0("Temp_",wh$breaks[1:ncol(mat.save)])
#mat.save = as.data.frame(mat.save,stringsAsFactors=F)

#colnames(mat.save) = c("taxa","temperature","value")
#mat.save[,1] = factor(mat.save[,1])
#mat.save[,2] = as.numeric(as.character(mat.save[,2]))
#mat.save[,3] = as.numeric(as.character(mat.save[,3]))

wm = apply(taxad,2,function(x) weighted.mean(meta$temp,x,na.rm=T))
hist(meta$temp,col="black",breaks=seq(-2,13,1),probability = T,ylim=c(0,0.4))
hist(wm,breaks=seq(-2,13,1),add=T,col=rgb(1,1,1,0.5),probability = T,ylim=c(0,0.4))

# mat.plot = mat.save[rowSums(mat.save)>0,]
mat.plot = mat.save
mat.side = colSums(taxad)
#mat.ord = order(apply(mat.plot,1,which.max))
mat.ord = order(wm)
mat.plot = mat.plot[mat.ord,]
mat.side = mat.side[mat.ord]
mat.king = king[mat.ord]
mat.clas = clas[mat.ord]
colnames(mat.plot) = wh$breaks[1:ncol(mat.save)]




rc1 = viridis(100)[rescale(round(sqrt(sqrt(mat.side))),c(1,100))]
#rc2 = viridis(length(unique(king)))[as.numeric(king)]

rc = data.frame(rc1)

for(i in 1:length(unique(mat.king))) {
  #color by class
  clasnum = (mat.king==unique(mat.king)[i]) * as.numeric(factor(mat.clas)) + 1
  tc = c("white",viridis(length(unique(mat.clas))+20))
  rc2 = tc[clasnum]
rc = cbind(rc,rc2)
}
colnames(rc) = c("relative abundance",unique(mat.king))
#rcplot = rc[,c("relative abundance","Chlorophyta","Haptophyta", "Ochrophyta", "Dinoflagellata","Ciliophora")]
#rcplot = rc[,c("Nitrososphaeria","Bacteroidia","Alphaproteobacteria","Gammaproteobacteria","relative abundance")]
rcplot = rc[,c("Bacillariophyta","Chrysophyceae","Mamiellophyceae","Prymnesiophyceae","Dinophyceae","MAST","Spirotrichea","Syndiniales","relative abundance")]
#heatmap.plus(mat.plot,Colv=NA,Rowv=NA,scale="row",trace="none",col=viridis,RowSideColors = rc, margins=c(12,12),cexRow = 0.1,key = F)
heatmap.plus((mat.plot)^1,Colv=NA,Rowv=NA,col=viridis(12),scale="row",RowSideColors = as.matrix(rcplot), margins=c(12,80),cexRow = 0.1,cexCol=2)

pdf(file="rplot.pdf",width=6,height=16)
heatmap.plus(t(mat.plot)^1,Colv=NA,Rowv=NA,col=viridis(12),scale="col",ColSideColors = as.matrix(rcplot), margins=c(80,24),cexRow = 2,cexCol=0.1,ylab="Temperature (deg. C)")
dev.off()


# library(ggplot2)
# library(ggridges)
# ggplot(mat.save, aes(x = temperature, y = taxa)) + geom_density_ridges(aes(height=..density..,
#                                                                             weight=value),    
#                                                                         scale= 0.95,
#                                                                         stat="density") 
# 
# ggplot(mat.save, aes(x = temperature, y = ..density.., weight = value)) + geom_histogram()
# 


dev.off()




