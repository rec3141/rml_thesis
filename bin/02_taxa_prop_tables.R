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

write.table(phy.agg,file="output/02_table_taxa_summary.tsv")


######## Results: Diatoms -- Table & Figure - diatom relative abundance (SILVA)

# FAMILY
# Genus
# Relative Abundance (%)
# ESVs

t2.taxa = taxa[rownames(sub.meta),taxa_good$Class=="Diatomea"]

t.taxa = rowSums(t(t2.taxa))
t.match = match(names(t.taxa),short_names)
t.agg = data.frame("Family"=as.character(taxa_good$Family[t.match]),"Genus"=as.character(taxa_good$Genus[t.match]),"vals"=t.taxa,stringsAsFactors=F)
t.agg$Genus[grepl("_",t.agg$Genus)] = "Unidentified"
t.agg = t.agg[t.agg$vals>0,]

t.cast = aggregate(. ~ Family + Genus,t.agg, FUN="sum")
t.esvs = aggregate(. ~ Family + Genus,t.agg, FUN="length")
t.cast$esvs = t.esvs$vals
t.cast = t.cast[t.cast$vals>0,]

t.cast$vals = round(100*t.cast$vals/sum(t.cast$vals),1)
t.cast$order = rep(0,nrow(t.cast))
class.agg = sort(colSums(agg.mat(sub.taxa,"Family")))
for(i in 1:length(class.agg)) {
  t.cast$order[t.cast$Family==names(class.agg[i])] = i
}

t.cast = t.cast[order(t.cast$order,t.cast$Family,t.cast$vals,t.cast$esvs,decreasing=TRUE),]
t.cast$vals[t.cast$vals==0] = "<0.1"

write.table(t.cast,file="output/02_table_diatoms_relabund.tsv",sep="\t",quote=F,row.names=F)


## diatom cumulative relative abundance plot
t.plot = data.frame("ESVs"=1:nrow(t.agg),"vals"=t.agg$vals)
t.plot = t.plot[order(t.plot$vals,decreasing=TRUE),]
t.plot$Cumulative = cumsum(t.plot$vals)/sum(t.plot$vals)*100
t.plot$ESVs = 1:nrow(t.plot)
t.plot$breakpoint = t.plot$Cumulative-1:nrow(t.plot)
t.plot = t.plot[t.plot$vals>0,]
gd = ggplot(t.plot) + geom_line(aes(x=ESVs,y=Cumulative)) + xlab("Number of Diatom ESVs") + ylab("Cumulative Percentile")
gd
ggsave(file="output/02_figure_diatom_cumulative.png",device=png(),width=4,height=4)



######## Results: Picoeuks -- Table & Figure - Picoeuk relative abundance (SILVA)

# FAMILY
# Genus
# Relative Abundance (%)
# ESVs

t2.taxa = taxa[rownames(sub.meta),taxa_good$Kingdom=="Haptophyta" | taxa_good$Subkingdom=="Chlorophyta" | taxa_good$Class=="Chrysophyceae"]

t.taxa = rowSums(t(t2.taxa))
t.match = match(names(t.taxa),short_names)
t.agg = data.frame("Kingdom"=as.character(taxa_good$Kingdom[t.match]),"Genus"=as.character(taxa_good$Genus[t.match]),"vals"=t.taxa,stringsAsFactors=F)
t.agg$Genus[grepl("_",t.agg$Genus)] = "Unidentified"
t.agg = t.agg[t.agg$vals>0,]

t.cast = aggregate(. ~ Kingdom + Genus,t.agg, FUN="sum")
t.esvs = aggregate(. ~ Kingdom + Genus,t.agg, FUN="length")
t.cast$esvs = t.esvs$vals
t.cast = t.cast[t.cast$vals>0,]

t.cast$vals = round(100*t.cast$vals/sum(t.cast$vals),1)
t.cast$order = rep(0,nrow(t.cast))
class.agg = sort(colSums(agg.mat(sub.taxa,"Kingdom")))
for(i in 1:length(class.agg)) {
  t.cast$order[t.cast$Kingdom==names(class.agg[i])] = i
}

t.cast = t.cast[order(t.cast$order,t.cast$Kingdom,t.cast$vals,t.cast$esvs,decreasing=TRUE),]
t.cast$vals[t.cast$vals==0] = "<0.1"

write.table(t.cast,file="output/02_table_picos_relabund.tsv",sep="\t",quote=F,row.names=F)

## picos cumulative relative abundance plot
t.plot = data.frame("ESVs"=1:nrow(t.agg),"vals"=t.agg$vals)
t.plot = t.plot[order(t.plot$vals,decreasing=TRUE),]
t.plot = t.plot[t.plot$vals>0,]
t.plot$Cumulative = cumsum(t.plot$vals)/sum(t.plot$vals)*100
t.plot$ESVs = 1:nrow(t.plot)
t.plot$breakpoint = t.plot$Cumulative-1:nrow(t.plot)
gp = ggplot(t.plot) + geom_line(aes(x=ESVs,y=Cumulative)) + xlab("Number of Picoeukaryotic ESVs") + ylab("Cumulative Percentile")
gp
ggsave(file="output/02_figure_picos_cumulative.png",device=png(),width=4,height=4)

#####

source("bin/03_cluster_analysis.R")


