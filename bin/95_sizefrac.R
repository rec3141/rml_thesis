
#analysis of size fractions
meta.process = meta[meta$`station_type`=="P",]
meta.process$uniqstn = Reduce(paste0,meta.process[,c("project","station","depth_m")])
taxa.process = taxa[rownames(meta.process),]

taxa.coef1 = matrix(data=NA,nrow=ncol(taxa.process),ncol=2)
taxa.coef2 = matrix(data=NA,nrow=ncol(taxa.process),ncol=2)

# zz=data.frame(meta.process$uniqstn,meta.process$filter,taxa.process)
# colnames(zz)[1:2]=c("uniqstn","filter")
# zzz=acast(zz,uniqstn ~ filter,value.var=...,mean)


for(i in 1:ncol(taxa.process)) {
  aa=data.frame(meta.process$uniqstn,meta.process$filter,taxa.process[,i])
  colnames(aa)[1:2]=c("uniqstn","filter")
  bb=acast(aa,uniqstn~filter,mean)
  bb[is.nan(bb)] = NA
  bb = bb[rowSums(bb,na.rm=T)>0,,drop=F]
  bb = rescale(bb,to=c(0,1))
  
  ee = as.data.frame(bb[,c(1,2),drop=F])
  ee$names = rownames(ee)
  gg = as.data.frame(bb[,c(2,3),drop=F])
  gg$names = rownames(gg)
  dd = melt(ee,id.vars='names')
  hh = melt(gg,id.vars='names')
  
  try({
    cc = NULL
    cc = summary(lm(dd$value[dd$variable %in% c("0.2","3")] ~ dd$variable[dd$variable %in% c("0.2","3")]))
    taxa.coef1[i,] = cc$coefficients[2,c(1,4)]
    ff = NULL
    ff = summary(lm(hh$value[hh$variable %in% c("3","20")] ~ hh$variable[hh$variable %in% c("3","20")]))
    taxa.coef2[i,] = ff$coefficients[2,c(1,4)]
  })
  
}

taxa.coef = data.frame(taxa.coef1,taxa.coef2)
rownames(taxa.coef) = colnames(taxa.process)

# plot(sqrt(abs(taxa.coef[,1]))*sign(taxa.coef[,1]),sqrt(abs(taxa.coef[,3]))*sign(taxa.coef[,3]),pch=21,bg="grey")
# plot(taxa.coef[,1],taxa.coef[,3],pch=21,bg="grey",cex=.3)
# plot(taxa.coef[,1],taxa.coef[,3],pch=21,bg=NULL,cex=sqrt(-1*log(taxa.coef[,2])))
# 
# taxa.coef1 = matrix(data=NA,nrow=ncol(taxa.process),ncol=2)
# taxa.coef2 = matrix(data=NA,nrow=ncol(taxa.process),ncol=2)
# taxa.coef3 = matrix(data=NA,nrow=ncol(taxa.process),ncol=2)
# 
# for(i in 1:ncol(taxa.process)) {
#   
#   if(colSums(taxa.process>0)[i] <= 2) next
# 
#   #pick station/depth combos with taxa
#   lm.pick = meta.process$uniqstn %in% unique(meta.process$uniqstn[taxa.process[,i]>0])
#   lm.taxa = taxa.process[,i]
#   lm.taxa1 = as.numeric(scale(lm.taxa[lm.pick & meta.process$filter<5]))
#   lm.taxa2 = as.numeric(scale(lm.taxa[lm.pick & meta.process$filter>2]))
#   lm.taxa3 = as.numeric(scale(lm.taxa[lm.pick & (meta.process$filter<2 | meta.process$filter>5)]))
#   
#         try({
# 
#           # plot(taxa.process[lm.pick,i] ~ meta.process$filter[lm.pick])
# 
#           lm1 = NULL
#           lm1 = lm(lm.taxa1 ~ meta.process[rownames(taxa.process)[lm.pick & meta.process$filter<5],"filter"]);
#           # print(summary(lm1));
#           taxa.coef1[i,] = summary(lm1)$coefficients[2,c(1,4)]
#           
#           lm2 = NULL
#           lm2 = lm(lm.taxa2 ~ meta.process[rownames(taxa.process)[lm.pick & meta.process$filter>2],"filter"]);
#           # print(summary(lm2));
#           taxa.coef2[i,] = summary(lm2)$coefficients[2,c(1,4)]
#           
#           lm3 = NULL
#           lm3 = lm(lm.taxa3 ~ meta.process[rownames(taxa.process)[lm.pick &  (meta.process$filter<2 | meta.process$filter>5)],"filter"]);
#           # print(summary(lm3));
#           taxa.coef3[i,] = summary(lm3)$coefficients[2,c(1,4)]
#           
#         })
#         }

# taxa.coef = cbind(taxa.coef1,taxa.coef2,taxa.coef3)
# rownames(taxa.coef) = colnames(taxa.process)

# par(mfrow=c(1,3))
# plot((sqrt(sqrt(abs(taxa.coef[,1])))*sign(taxa.coef[,1])),log(taxa.coef[,2]),pch=21,bg=rgb(.5,.5,.5,.1),col=NULL)
# plot((sqrt(sqrt(abs(taxa.coef[,3])))*sign(taxa.coef[,3])),log(taxa.coef[,4]),pch=21,bg=rgb(.5,.5,.5,.1),col=NULL)

pval = 1
pval2 = 0.01
#bonferroni / all / 0.0002283105
#pval2 = 0.05/(2*nrow(taxa.coef)-sum(is.na(taxa.coef)))
#bonferroni / FL <-> PA
#pval = 0.05/(1*nrow(taxa.coef)-sum(is.na(taxa.coef[,1])))
#bonferroni / Suspended <-> Sinking
#pval = 0.05/(1*nrow(taxa.coef)-sum(is.na(taxa.coef[,2])))

#much higher confidence in 0.2/3 than 3/20 for bacteria
#slighltly higher confidence for 0.2/3 than 3/2 for protists
#overall lower confidence for protists

plot(  sign(taxa.coef[order(taxa.coef[,2]),1]) * -1*log10(taxa.coef[order(taxa.coef[,2]),2]),cex=0.001*sqrt(colSums(taxa.process)[order(taxa.coef[,2])]),col='blue')
points(sign(taxa.coef[order(taxa.coef[,2]),3]) * -1*log10(taxa.coef[order(taxa.coef[,2]),4]),cex=0.001*sqrt(colSums(taxa.process)[order(taxa.coef[,2])]),col='red')
lines(y=c(log10(0.01),log10(0.01)),x=c(-10,1000))
lines(y=c(log10(0.001),log10(0.001)),x=c(-10,1000))
lines(y=c(log10(pval),log10(pval)),x=c(-10,1000))
lines(y=c(-1*log10(0.01),-1*log10(0.01)),x=c(-10,1000))
lines(y=c(-1*log10(0.001),-1*log10(0.001)),x=c(-10,1000))
lines(y=c(-1*log10(pval),-1*log10(pval)),x=c(-10,1000))

plot(taxa.coef[,3],log10(taxa.coef[,4]),cex=0.001*sqrt(colSums(taxa.process)))
lines(y=c(log10(0.01),log10(0.01)),x=c(-10,1000))
lines(y=c(log10(0.001),log10(0.001)),x=c(-10,1000))
lines(y=c(log10(pval),log10(pval)),x=c(-10,1000))
lines(y=c(1,-20),x=c(0,0))

taxa.coef$new1 = taxa.coef[,1]
taxa.coef$new1[which(taxa.coef[,2] >= pval2)] = 0
taxa.coef$new3 = taxa.coef[,3]
taxa.coef$new3[which(taxa.coef[,4] >= pval2)] = 0

pdf(file="sizefrac-new3.pdf",width=10,height=10)

xlims = c(-0.5,0.5)
ylims = c(-0.5,0.5)
pick.pval = taxa.coef[,2]<pval | taxa.coef[,4]<pval

taxlevel1 = "Order"
taxlevel2 = "Genus"
tax2plot = unique(taxa_split[,taxlevel1])
# tax2plot = unique(taxa_split[taxa_split$Phylum=="Ochrophyta","Family"])
# col=hcl.colors(length(tax2plot), palette = "viridis", alpha = NULL, rev = FALSE, fixup = TRUE)
col=rainbow(length(tax2plot))

plot(taxa.coef[pick.pval,1],taxa.coef[pick.pval,3],bg="grey",col=NULL,pch=21,cex=0.3,xlim=xlims,ylim=ylims,xlab="trend from 0.2um to 3um",ylab="trend from 3um to 20um")
lines(c(-10,10,0,0),c(0,0,-10,10),col='grey')

for(i in 1:length(tax2plot)) {
  pick.taxa = pick.pval & taxa_split[,taxlevel1]==tax2plot[i]
  print(sum(pick.taxa,na.rm=T))
  #  points(taxa.coef[pick.taxa,5],taxa.coef[pick.taxa,6], cex=0.002*sqrt(colSums(taxa.process)[pick.taxa]),col=col[i])
  points(taxa.coef[pick.taxa,5],taxa.coef[pick.taxa,6], cex=0.3,col=col[i])
}
text(x=rep(0.9*xlims[1],length(tax2plot)),y=ylims[2]-1.2*(1:length(tax2plot))/length(tax2plot),labels=tax2plot,col=col)

for(i in 1:length(tax2plot)) {
  pick.taxa = pick.pval & taxa_split[,taxlevel1]==tax2plot[i]
  plot(taxa.coef[pick.taxa,1],taxa.coef[pick.taxa,3],bg="grey",col=NULL,pch=21,cex=0.3,xlim=xlims,ylim=ylims,xlab="trend from 0.2um to 3um",ylab="trend from 3um to 20um",main=tax2plot[i])
  lines(c(-10,10,0,0),c(0,0,-10,10),col='grey')
  
  tax2plot2 = unique(taxa_split[taxa_split[,taxlevel1]==tax2plot[i],taxlevel2])
  col=rainbow(length(tax2plot2))
  for(j in 1:length(tax2plot2)) {
    pick.taxa2 = pick.pval & taxa_split[,taxlevel2]==tax2plot2[j] & taxa_split[,taxlevel1]==tax2plot[i]
    
    points(taxa.coef[pick.taxa2,1],taxa.coef[pick.taxa2,3], cex=0.002*sqrt(colSums(taxa.process)[pick.taxa2]),col=col[j])
    text(x=rep(0.9*xlims[1],length(tax2plot2))[j],y=(0.9*ylims[2]-(1:length(tax2plot2))/length(tax2plot2))[j],labels=tax2plot2[j],col=col[j])
  }
}

dev.off()



taxa.tcoef = cbind(abs(taxa.coef[,1])^(1/2)*sign(taxa.coef[,1]),abs(taxa.coef[,3])^(1/2)*sign(taxa.coef[,3]))
plot(taxa.tcoef[pick.pval,1],taxa.tcoef[pick.pval,2],cex=0.003*sqrt(colSums(taxa.process)))

pul = taxa.tcoef[intersect(which(taxa.tcoef[,1]<0),intersect(which(taxa.tcoef[,2]>0), which(pick.pval))),,drop=F]
pur = taxa.tcoef[intersect(which(taxa.tcoef[,1]>0),intersect(which(taxa.tcoef[,2]>0), which(pick.pval))),,drop=F]
pll = taxa.tcoef[intersect(which(taxa.tcoef[,1]<0),intersect(which(taxa.tcoef[,2]<0), which(pick.pval))),,drop=F]
plr = taxa.tcoef[intersect(which(taxa.tcoef[,1]<0),intersect(which(taxa.tcoef[,2]>0), which(pick.pval))),,drop=F]

taxa_groups = NA*1:ncol(taxa)
taxa_groups[grepl("Bacillariophyta",original_names)] = "Diatom"
taxa_groups[grepl("Haptophyta",original_names) | grepl("Chlorophyta",original_names) | grepl("Cryptophyta",original_names) | grepl("Chrysophy",original_names) | grepl("Raphidophyceae",original_names) | grepl("Eustigmatophyceae",original_names) | grepl("Dictyochophyceae",original_names)] = "Picophytoplankton"

#silva
# taxa_groups[grepl("Diatom",original_names)] = "Diatom"
# taxa_groups[grepl("Haptophyta",original_names) | grepl("Chlorophyta",original_names) | grepl("Cryptomonad",original_names) | grepl("Chrysophy",original_names) | grepl("Raphidophyceae",original_names) | grepl("Eustigmatophyceae",original_names) | grepl("Dictyochophyceae",original_names)] = "Picophytoplankton"
# taxa_groups[grepl("Fungi",original_names) | grepl("Ichthyosporea",original_names)] = "Fungi"
# taxa_groups[grepl("Dinoflagellata",original_names)] = "Dinoflagellates"
# taxa_groups[grepl("Ciliophora",original_names)] = "Ciliates"
# taxa_groups[grepl("Syndiniales",original_names)] = "Syndiniales"
# taxa_groups[grepl("MAST",original_names)] = "MAST"
# taxa_groups[grepl("Amoebozoa",original_names)] = "Amoebozoa"
# taxa_groups[grepl("Phaeophyceae",original_names)] = "Phaeophyceae"
# taxa_groups[grepl("Thecofilosea",original_names)] = "Thecofilosea"
# taxa_groups[grepl("Pelagophyceae",original_names)] = "Pelagophyceae"




# select only survey stations and 0.2um filters -> 392 samples
d.meta = meta[meta$`station_type`=="S" & meta$filter==0.2 ,]
p.meta = d.meta
f.meta = d.meta

}

filtered_samples = intersect(rownames(taxa),rownames(p.meta))



#filter diatoms
d.taxa = taxa[,grepl("Diatom",original_names)]

d.meta = d.meta[filtered_samples,]
d.taxa = d.taxa[filtered_samples,]
d.taxa = d.taxa[,colSums(d.taxa)>0] #remove empty taxa

# limit to taxa present in multiple samples -- doesn't change much
#d.taxa = d.taxa[rowSums(d.taxa)>0,colSums(d.taxa>0)>1]
#d.taxa = d.taxa[,colSums(d.taxa>0)>1]

# do we want proportions within group?
d.prop = as.data.frame(prop.table(as.matrix(d.taxa),margin=1))

## filter picos
#Chlorophyta, Haptophyta, or Chrysophyceae
#p.taxa = taxa[,grepl("Haptophyta",original_names) | grepl("Chlorophyta",original_names) | grepl("Cryptomonad",original_names) | grepl("Chrysophy",original_names)]
p.taxa = taxa[,grepl("Dinoflag",original_names) ]

p.meta = p.meta[filtered_samples,]
p.taxa = p.taxa[filtered_samples,]
p.taxa = p.taxa[,colSums(p.taxa)>0]

# limit to taxa present in multiple samples -- doesn't change much
#p.taxa = p.taxa[,colSums(p.taxa>0)>1]

#proportions within group?
p.prop = as.data.frame(prop.table(as.matrix(p.taxa),margin=1))


## filter fungi
f.taxa = taxa[,grepl("Fungi",original_names) ]

f.meta = f.meta[filtered_samples,]
f.taxa = f.taxa[filtered_samples,]
f.taxa = f.taxa[,colSums(f.taxa)>0]

# limit to taxa present in multiple samples -- doesn't change much
#f.taxa = f.taxa[,colSums(f.taxa>0)>1]

#proportions within group?
f.prop = as.data.frame(prop.table(as.matrix(f.taxa),margin=1))


## filter ciliates
f.taxa = taxa[,grepl("Fungi",original_names) ]

f.meta = f.meta[filtered_samples,]
f.taxa = f.taxa[filtered_samples,]
f.taxa = f.taxa[,colSums(f.taxa)>0]

# limit to taxa present in multiple samples -- doesn't change much
#f.taxa = f.taxa[,colSums(f.taxa>0)>1]

#proportions within group?
f.prop = as.data.frame(prop.table(as.matrix(f.taxa),margin=1))


