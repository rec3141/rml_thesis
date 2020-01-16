## Taxonomy tables

######## Results: Overall Structure (SILVA)

require(ggplot2)

#diatoms = Diatomea
#picos = Chlorophyta, Haptophyta, or Chrysophyceae

sub.lab = rep("Other",ncol(sub.taxa))
sub.lab[grepl("Haptophyta",taxa_split$Kingdom) | grepl("Chlorophyta",taxa_split$Subkingdom) | grepl("Chrysophyceae",taxa_split$Class)] = "Pico"
sub.lab[grepl("Dinoflagellata",taxa_split$Phylum) ] = "Dino"
sub.lab[grepl("Diatomea",taxa_split$Class) ] = "Diatom"
sub.lab[grepl("Ciliophora",taxa_split$Phylum) ] = "Ciliate"

st.agg = aggregate(sub.taxa,by=list(sub.meta$project),FUN=sum)
rownames(st.agg) = st.agg[,1]
st.agg = st.agg[,-1]
st.agg = t(st.agg)

#aggregate at Phylum level
phy.agg = aggregate(st.agg,by = list(taxa_split$Phylum),FUN = sum)

rownames(phy.agg) = phy.agg[,1]
phy.agg = phy.agg[,-1]
phy.agg = 100*phy.agg/colSums(phy.agg)

phy.agg.sd = apply(phy.agg,1,sd)
phy.agg.mean = apply(phy.agg,1,mean)

phy.agg = data.frame(phy.agg,phy.agg.mean,phy.agg.sd)
phy.agg$Phylum = taxa_split$Phylum[match(rownames(phy.agg),taxa_split$Genus)]
phy.agg = phy.agg[order(phy.agg$Phylum, phy.agg$phy.agg.mean,decreasing = T),]
phy.agg = phy.agg[phy.agg$phy.agg.mean>0,]

write.table(format(phy.agg,digits=3,scientific=FALSE),file="output/02_table_taxa_summary.tsv", quote=F, sep="\t",)


######## Results: Diatoms -- Table & Figure - diatom relative abundance (SILVA)


# FAMILY
# Genus
# Relative Abundance (%)
# ESVs

cmesv.list = list()

for(tax in names(mat.list)) {

  x.taxa = mat.list[[tax]][["taxa"]]

# THIS DOES NOT CURRENTLY FILTER OUT BY TAXA_GOOD

#not sure how do deal with Insertae, uncultured, unidentified...
x.pick = match(colnames(x.taxa),taxa_good$ESV)

t2.taxa = sub.taxa[,x.pick]

t.abund = colSums(t2.taxa)

t.agg = data.frame("Kingdom"=as.character(taxa_good$Kingdom[x.pick]), "Family"=as.character(taxa_good$Family[x.pick]),"Genus"=as.character(taxa_good$Genus[x.pick]),"vals"=t.abund,stringsAsFactors=F)
#hack...
t.agg$Genus[grepl("_",t.agg$Genus)] = "Unidentified"
t.agg$Genus[grepl("uncultured",t.agg$Genus)] = "Unidentified"
t.agg$Genus[grepl("Incertae sedis",t.agg$Genus)] = "Unidentified"
t.agg$Genus[grepl("^root$",t.agg$Genus)] = "Unidentified"

#remove empties before counting ESVs
t.agg = t.agg[t.agg$vals>0,]

t.cast = aggregate(. ~ Kingdom + Family + Genus,t.agg, FUN="sum")
t.esvs = aggregate(. ~ Kingdom + Family + Genus,t.agg, FUN="length")
t.cast$esvs = t.esvs$vals

### ESV vs relabund
t.ident = t.cast[t.cast$Genus!="Unidentified",]
cor.test(as.numeric(t.ident$esvs), as.numeric(t.ident$vals), method = "spearman")
nrow(t.ident)

#sort table
t.cast$vals = round(100*t.cast$vals/sum(t.cast$vals),1)
t.cast$order = rep(0,nrow(t.cast))

king.agg = aggregate(. ~ Kingdom,t.agg[,c("Kingdom","vals")], FUN="sum")
king.agg = king.agg[order(king.agg$vals),]
class.agg = aggregate(. ~ Family,t.agg[,c("Family","vals")], FUN="sum")
class.agg = class.agg[order(class.agg$vals),]

for(i in 1:nrow(king.agg)) {
  t.cast$king.order[t.cast$Kingdom==king.agg$Kingdom[i]] = i
}
for(i in 1:nrow(class.agg)) {
  t.cast$class.order[t.cast$Family==class.agg$Family[i]] = i
}

t.cast = t.cast[order(t.cast$king.order, t.cast$class.order, t.cast$Family, t.cast$vals, t.cast$esvs, decreasing=c(TRUE, TRUE, FALSE, TRUE, TRUE)),]
t.cast$vals[t.cast$vals==0] = "<0.1"
t.cast$order = NULL
colnames(t.cast) = c("Kingdom","Family","Genus","Group Relative Abundance","ESVs")

write.table(t.cast,file=paste0("output/02_table_",tax,"_relabund.tsv"),sep="\t",quote=F,row.names=F)

### most abundant at each depth
t.genus = data.frame("Genus"=as.character(taxa_good$Genus[x.pick]),"vals"=t.abund,stringsAsFactors=F)

t.depth = aggregate(sub.taxa,by=list(sub.meta$depth_type), FUN="sum")
rownames(t.depth) = t.depth[,1]
t.depth = t.depth[,-1]
t.depth = t(t.depth)

t.depag = aggregate(t.depth, by=list(1:ncol(taxa) %in% x.pick, taxa_good$Genus), FUN="sum")
t.depag = t.depag[t.depag$Group.1,]
rownames(t.depag) = make.unique(t.depag$Group.2)
t.depag$Group.1 = NULL
t.depag$Group.2 = NULL
t.depag = 100*prop.table(as.matrix(t.depag),2)
t.depag = t.depag[order(rowSums(t.depag), decreasing=TRUE),]

write.table(format(t.depag,digits=2,scientific=FALSE),file=paste0("output/02_table_",tax,"_relabund_depth.tsv"),sep="\t",quote=F,row.names=T)


## taxa cumulative relative abundance plot
t.plot = data.frame("ESVs"=1:nrow(t.agg),"vals"=t.agg$vals)
t.plot = t.plot[order(t.plot$vals,decreasing=TRUE),]
t.plot$Cumulative = cumsum(t.plot$vals)/sum(t.plot$vals)*100
t.plot$ESVs = 1:nrow(t.plot)
t.plot$breakpoint = t.plot$Cumulative-1:nrow(t.plot)
t.plot = t.plot[t.plot$vals>0,]
gd = ggplot(t.plot) + geom_line(aes(x=ESVs,y=Cumulative)) + ylim(0,100) +
  xlab("Number of ESVs") + 
  ylab("Cumulative Percentile") +
  ggtitle(tax)
gd
ggsave(file=paste0("output/02_figure_",tax,"_cumulative.png"),device=png(),width=4,height=4)

cmesv.list[[tax]] = t.plot$Cumulative

} # end taxa loop


# cmesv.list = cmesv.list[c("picos","diatoms")]
t.plot = lapply(names(cmesv.list), function(x) data.frame("taxa"=rep(x,length(cmesv.list[[x]])), "ESVs"=1:length(cmesv.list[[x]]), "Cumulative"=cmesv.list[[x]], stringsAsFactors=F))
t.plot = do.call("rbind", t.plot)

gd = ggplot(t.plot) + geom_line(aes(x=ESVs,y=Cumulative, col=taxa)) + ylim(50,100) + xlim(0,300) +
  xlab("Number of ESVs") + 
  ylab("Cumulative Percentile") +
  ggtitle(tax)
gd

ggsave(file=paste0("output/02_figure_alltax_cumulative.png"),device=png(),width=4,height=4)


graphics.off()
source("bin/03_cluster_analysis.R")




# FAMILY
# Genus
# Relative Abundance (%)
# ESVs

t2.taxa = taxa[rownames(sub.meta),taxa_good$Class=="Diatomea"]

t.taxa = rowSums(t(t2.taxa))
t.match = match(names(t.taxa),taxa_good$ESV)
t.agg = data.frame("Family"=as.character(taxa_good$Family[t.match]),"Genus"=as.character(taxa_good$Genus[t.match]),"vals"=t.taxa,stringsAsFactors=F)
t.agg$Genus[grepl("_",t.agg$Genus)] = "Unidentified"
t.agg$Genus[grepl("uncultured",t.agg$Genus)] = "Unidentified"
t.agg$Genus[grepl("Incertae sedis",t.agg$Genus)] = "Unidentified"
t.agg = t.agg[t.agg$vals>0,]

t.cast = aggregate(. ~ Family + Genus,t.agg, FUN="sum")
t.esvs = aggregate(. ~ Family + Genus,t.agg, FUN="length")
t.cast$esvs = t.esvs$vals
t.cast = t.cast[t.cast$vals>0,]

### ESV vs relabund
t.ident = t.cast[t.cast$Genus!="Unidentified",]
cor.test(as.numeric(t.ident$esvs), as.numeric(t.ident$vals), method = "spearman")
nrow(t.ident)


t.cast$vals = round(100*t.cast$vals/sum(t.cast$vals),1)
t.cast$order = rep(0,nrow(t.cast))
class.agg = sort(colSums(agg.mat(sub.taxa,"Family")))
for(i in 1:length(class.agg)) {
  t.cast$order[t.cast$Family==names(class.agg[i])] = i
}

t.cast = t.cast[order(t.cast$order,t.cast$Family,t.cast$vals,t.cast$esvs,decreasing=TRUE),]
t.cast$vals[t.cast$vals==0] = "<0.1"

write.table(t.cast,file="output/02_table_Diatoms-paper_relabund.tsv",sep="\t",quote=F,row.names=F)


### most abundant at each depth
t.genus = data.frame("Genus"=as.character(taxa_good$Genus[t.match]),"vals"=t.taxa,stringsAsFactors=F)

t.depth = aggregate(sub.taxa,by=list(sub.meta$depth_type), FUN="sum")
rownames(t.depth) = t.depth[,1]
t.depth = t.depth[,-1]
t.depth = t(t.depth)

t.depag = aggregate(t.depth, by=list(taxa_good$Class=="Diatomea", taxa_good$Genus), FUN="sum")
rownames(t.depag) = t.depag$Group.2
t.depag = t.depag[t.depag$Group.1,]
t.depag$Group.1 = NULL
t.depag$Group.2 = NULL
t.depag = 100*prop.table(as.matrix(t.depag),2)
t.depag = t.depag[order(rowSums(t.depag), decreasing=TRUE),]

write.table(format(t.depag,digits=2,scientific=FALSE),file="output/02_table_Diatoms-paper_relabund_depth.tsv",sep="\t",quote=F,row.names=T)


## diatom cumulative relative abundance plot
t.plot = data.frame("ESVs"=1:nrow(t.agg),"vals"=t.agg$vals)
t.plot = t.plot[order(t.plot$vals,decreasing=TRUE),]
t.plot$Cumulative = cumsum(t.plot$vals)/sum(t.plot$vals)*100
t.plot$ESVs = 1:nrow(t.plot)
t.plot$breakpoint = t.plot$Cumulative-1:nrow(t.plot)
t.plot = t.plot[t.plot$vals>0,]
gd = ggplot(t.plot) + geom_line(aes(x=ESVs,y=Cumulative)) + xlab("Number of Diatom ESVs") + ylab("Cumulative Percentile") + ylim(0,100)
gd
ggsave(file="output/02_figure_Diatoms-paper_cumulative.png",device=png(),width=4,height=4)



######## Results: Picoeuks -- Table & Figure - Picoeuk relative abundance (SILVA)

# FAMILY
# Genus
# Relative Abundance (%)
# ESVs

t2.taxa = taxa[rownames(sub.meta),taxa_good$Kingdom=="Haptophyta" | taxa_good$Subkingdom=="Chlorophyta" | taxa_good$Class=="Chrysophyceae"]

t.taxa = rowSums(t(t2.taxa))
t.match = match(names(t.taxa),taxa_split$ESV)
t.agg = data.frame("Kingdom"=as.character(taxa_good$Kingdom[t.match]),"Genus"=as.character(taxa_good$Genus[t.match]),"vals"=t.taxa,stringsAsFactors=F)
t.agg$Genus[grepl("_",t.agg$Genus)] = "Unidentified"
t.agg$Genus[grepl("uncultured",t.agg$Genus)] = "Unidentified"
t.agg$Genus[grepl("Incertae sedis",t.agg$Genus)] = "Unidentified"
t.agg = t.agg[t.agg$vals>0,]

t.cast = aggregate(. ~ Kingdom + Genus,t.agg, FUN="sum")
t.esvs = aggregate(. ~ Kingdom + Genus,t.agg, FUN="length")
t.cast$esvs = t.esvs$vals
t.cast = t.cast[t.cast$vals>0,]

### ESV vs relabund
t.ident = t.cast[t.cast$Genus!="Unidentified",]
cor.test(as.numeric(t.ident$esvs), as.numeric(t.ident$vals), method = "spearman")
nrow(t.ident)
plot(t.ident$esvs ~ t.ident$vals)

t.cast$vals = round(100*t.cast$vals/sum(t.cast$vals),1)
t.cast$order = rep(0,nrow(t.cast))
class.agg = sort(colSums(agg.mat(sub.taxa,"Kingdom")))
for(i in 1:length(class.agg)) {
  t.cast$order[t.cast$Kingdom==names(class.agg[i])] = i
}

t.cast = t.cast[order(t.cast$order,t.cast$Kingdom,t.cast$vals,t.cast$esvs,decreasing=TRUE),]
t.cast$vals[t.cast$vals==0] = "<0.1"

write.table(t.cast,file="output/02_table_Picos-paper_relabund.tsv",sep="\t",quote=F,row.names=F)



### most abundant at each depth
t.genus = data.frame("Genus"=as.character(taxa_good$Genus[t.match]),"vals"=t.taxa,stringsAsFactors=F)

t.depth = aggregate(sub.taxa,by=list(sub.meta$depth_type), FUN="sum")
rownames(t.depth) = t.depth[,1]
t.depth = t.depth[,-1]
t.depth = t(t.depth)

t.depag = aggregate(t.depth, by=list(taxa_good$Kingdom=="Haptophyta" | taxa_good$Subkingdom=="Chlorophyta" | taxa_good$Class=="Chrysophyceae", taxa_good$Genus), FUN="sum")
rownames(t.depag) = make.unique(t.depag$Group.2)
t.depag = t.depag[t.depag$Group.1,]
t.depag$Group.1 = NULL
t.depag$Group.2 = NULL
t.depag = 100*prop.table(as.matrix(t.depag),2)
t.depag = t.depag[order(rowSums(t.depag), decreasing=TRUE),]

write.table(format(t.depag,digits=2,scientific=FALSE),file="output/02_table_Picos-paper_relabund_depth.tsv",sep="\t",quote=F,row.names=T)


## picos cumulative relative abundance plot
t.plot = data.frame("ESVs"=1:nrow(t.agg),"vals"=t.agg$vals)
t.plot = t.plot[order(t.plot$vals,decreasing=TRUE),]
t.plot = t.plot[t.plot$vals>0,]
t.plot$Cumulative = cumsum(t.plot$vals)/sum(t.plot$vals)*100
t.plot$ESVs = 1:nrow(t.plot)
t.plot$breakpoint = t.plot$Cumulative-1:nrow(t.plot)
gp = ggplot(t.plot) + geom_line(aes(x=ESVs,y=Cumulative)) + 
  xlab("Number of Picoeukaryotic ESVs") + ylab("Cumulative Percentile") +
  ylim(0,100)
gp
ggsave(file="output/02_figure_Picos-paper_cumulative.png",device=png(),width=4,height=4)






