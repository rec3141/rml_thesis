# prepare metadata and taxonomic tables
# run first

#require(MASS)
#require(plyr)
#require(tidyverse)
# require(reshape2)
# require(gplots)
# require(ggplot2)
# require(compositions)
# require(vegan)
# require(fastcluster)
# require(betareg)
# require(scales)
# require(SDMTools)
# require(e1071)
# require(locClass)

setwd("~/Downloads/thesis-master/201912")
#graphics.off()
#options(device="quartz")

# import some useful functions
source("bin/00_functions.R")

tdb = "SILVA" #or "PR2"

#meta: fixed out of range DO values
#meta: fixed bad AMBON2017 fluorescence calibration by 10*

#run after updating metadata table in excel
#meta <- read.table('metaRDS_denovo.csv',header = T,comment.char = "",sep="\t",row.names=1, check.names = F,stringsAsFactors = F)
#saveRDS(meta,file="meta_denovo.RDS")

meta <- readRDS("data/meta_denovo.RDS")

meta$mass = factor(get_water_mass(meta))
meta$project <- factor(meta$project, levels=c('ASGARD2017', 'AMBON2017', 'DBONCIS2017'))

#calculate stratification index
meta$rho <- swRho(meta$salinity,meta$temp,meta$depth_m)
meta$projcast = paste0(meta$project,meta$CTD_number)

stdf <- data.frame("projcast"=meta$projcast,
                   "rho"=meta$rho,"z"=meta$depth_m,
                   "bottom"=meta$bottom_depth)
strat <- data.frame(stratification=unclass(by(stdf, stdf$projcast, stratify)))
strat$projcast = rownames(strat)
tmp.meta = merge(meta,strat)
rownames(tmp.meta) = rownames(meta)
meta = tmp.meta

#calculate phaeopigment ratio
meta$phaeofrac = meta$`phaeo (ug/l)`/meta$`chl (ug/l)`

# import taxa info from DADA2
seqtab = readRDS("data/seqtab_filt.rds")
tablelist = readRDS("data/table_list.rds") #based on SILVA, not PR2
tablelist[is.na(tablelist)] = ""
nameslist = readRDS("data/names_list.rds")

if(tdb=="PR2") taxadb = "ref_dada2_pr2_version_4.10.0.fasta"
if(tdb=="SILVA") taxadb = "ref_dada2_silva_v132.fasta"

namedlist = nameslist[[taxadb]]
names(namedlist) = colnames(seqtab)

#taxa = seqtab[,tablelist=="16S_prokaryote"]
taxa = seqtab[,tablelist=="18S_protist"]

#remove empty rows and columns
taxa <- taxa[which(rowSums(taxa>0)>0),which(colSums(taxa>0)>0)] # omits columns/rows with 0
taxa <- taxa[which(rowSums(taxa>0)>0),which(colSums(taxa>0)>0)] # omits columns/rows with 0

#get proportional table
taxa = 1e6*as.data.frame(prop.table(taxa,margin=1))

#remove data for which there is missing metadata or taxonomic info
common_samples <- intersect(rownames(meta),rownames(taxa))
meta <- meta[common_samples, ] #916 -> 639
taxa <- taxa[common_samples, ] #734 -> 639

taxa_seqs = colnames(taxa)

bootout = readRDS("data/bootout_edit.rds")
boot = bootout[[taxadb]]
boot = data.frame(boot[taxa_seqs,])

original_names <- unname(namedlist[taxa_seqs])
taxa_split = split_by_character(original_names,";")
#PR2
if(tdb=="PR2") colnames(taxa_split) = c("Domain","Superphylum","Phylum","Class","Order","Family","Genus","Species","ESV")
if(tdb=="SILVA") colnames(taxa_split) = c("root","Domain","Major Clade","Superkingdom","Kingdom","Subkingdom","Infrakingdom","Superphylum","Phylum","Subphylum","Infraphylum","Superclass","Class","Subclass","Infraclass","Superorder","Order","Suborder","Superfamily","Family","Subfamily","Genus","Accession","ESV")

short_names <- paste(taxa_split$Family, taxa_split$Genus, taxa_split$ESV, boot$Genus, sep=' ')
colnames(taxa) <- short_names
rownames(boot) = short_names

#rearrange matrix to taxonomic ordering? not now
#taxa = taxa[,order(original_names)]

### produce taxonomic quality-controlled taxa table
tcutoff = 60 #cutoff for taxonomic id
taxa_good = taxa_split

# remove taxa lacking bootstrap support
for(i in 1:nrow(taxa_good)) {
  curtax=taxa_good[i,1]
  for(j in 1:ncol(boot)) {
    if(boot[i,j] >= tcutoff) {
      curtax = taxa_split[i,j]
      next
    } else {
      taxa_good[i,j] = curtax
    }
  }
}

### Make metadata subsets for diatoms and picos

# select only survey stations and 0.2um filters -> 392 samples
sub.meta = meta[meta$station_type=="S" & meta$filter==0.2,]
sub.meta$depth_type = factor(sub.meta$depth_type, levels=c('bottom', 'mid', 'surf'), labels=c("bottom","midwater","surface"), ordered=TRUE)
sub.meta$project = factor(sub.meta$project, levels=c("ASGARD2017","AMBON2017","DBONCIS2017"),labels=c("ASGARD","AMBON","DBO-NCIS"), ordered=TRUE)

# aggregate by higher taxonomic group for sorting
sub.taxa = taxa[rownames(sub.meta),]
sub.agg = agg.mat(sub.taxa,"Phylum")

#filter diatoms
d.pick = intersect(rownames(taxa),rownames(sub.meta))

d.taxa = taxa[,grepl("Diatom",original_names)]
d.taxa = d.taxa[d.pick,]
d.taxa = d.taxa[rowSums(d.taxa)>0,]

d.meta = sub.meta[rownames(d.taxa),]
d.taxa = d.taxa[,colSums(d.taxa)>0] #remove empty taxa

d.agg = sub.agg[rownames(d.taxa),]

d.meta = data.frame(d.meta,d.agg,check.names=F)

# limit to taxa present in multiple samples -- doesn't change much
#d.taxa = d.taxa[rowSums(d.taxa)>0,colSums(d.taxa>0)>1]
#d.taxa = d.taxa[,colSums(d.taxa>0)>1]

# do we want proportions within group?
d.prop = as.data.frame(prop.table(as.matrix(d.taxa),margin=1))

## filter picos
#Chlorophyta, Haptophyta, or Chrysophyceae

p.pick = intersect(rownames(taxa),rownames(sub.meta))

#p.taxa = taxa[,grepl("Haptophyta",original_names) | grepl("Chlorophyta",original_names) | grepl("Cryptomonad",original_names) | grepl("Chrysophy",original_names)]
p.taxa = taxa[,grepl("Haptophyta",original_names) | grepl("Chlorophyta",original_names) | grepl("Chrysophy",original_names)]
p.taxa = p.taxa[p.pick,]
p.taxa = p.taxa[rowSums(p.taxa)>0,]

p.meta = sub.meta[rownames(p.taxa),]
p.taxa = p.taxa[,colSums(p.taxa)>0] #remove empty taxa
p.agg = sub.agg[rownames(p.taxa),]
p.meta = data.frame(p.meta,p.agg,check.names=F)

# limit to taxa present in multiple samples -- doesn't change much
#p.taxa = p.taxa[,colSums(p.taxa>0)>1]

#proportions within group?
p.prop = as.data.frame(prop.table(as.matrix(p.taxa),margin=1))

######################



