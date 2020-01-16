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
require(oce)

#setwd("~/Downloads/thesis-master/201912")
#graphics.off()
#options(device="quartz")

# import some useful functions
source("bin/00_functions.R")

tdb = "SILVA" #or "PR2"

#meta: fixed out of range DO values
#meta: fixed bad AMBON2015 and AMBON2017 fluorescence calibration by 10*
#meta: added AMBON2015 to metaRDS_denovo_2

#run after updating metadata table in excel ande exporting as tab-delimited csv
#meta <- read.table('data/metaRDS_denovo_2.csv',header = T,comment.char = "",sep="\t",row.names=1, check.names = F,stringsAsFactors = F)
#saveRDS(meta,file="data/meta_denovo_2.RDS")

meta = readRDS("data/meta_denovo_2.RDS")

meta$mass = factor(get_water_mass(meta))
meta$depth_type = factor(meta$depth_type, levels=c('bottom', 'mid', 'surf'), labels=c("bottom","midwater","surface"), ordered=TRUE)
meta$project = factor(meta$project, levels=c("ASGARD2017","AMBON2015","AMBON2017","DBONCIS2017"),labels=c("ASGARD 2017","AMBON 2015","AMBON 2017","DBO-NCIS 2017"), ordered=TRUE)

#calculate stratification index
meta$rho <- swRho(meta$salinity,meta$temp,meta$depth_m)
meta$projcast = paste0(meta$project,meta$CTD_number)

stdf <- data.frame("projcast"=meta$projcast,
                   "rho"=meta$rho,"z"=meta$depth_m,
                   "bottom"=meta$bottom_depth)
strat <- data.frame(stratification=unclass(by(stdf, stdf$projcast, stratify)))
strat$projcast = rownames(strat)

tmp.meta = merge(meta,strat)
rownames(tmp.meta) = tmp.meta$cellid
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

# limit to taxonomic subset
#taxa = seqtab[,tablelist=="16S_prokaryote"]
#taxa = seqtab[,tablelist=="18S_protist"]
taxa = seqtab[,tablelist=="18S_protist" | tablelist=="16S_prokaryote" | tablelist=="18S_metazoa"]

#remove empty rows and columns
taxa = taxa[which(rowSums(taxa>0)>0),which(colSums(taxa>0)>0)] # omits columns/rows with 0 entries

#get proportional table (counts per million)
taxa = 1e6*as.data.frame(prop.table(taxa,margin=1))

#remove data for which there is missing metadata or taxonomic info
common_samples = intersect(rownames(meta),rownames(taxa))
meta = meta[common_samples, ] #916 -> 639
taxa = taxa[common_samples, ] #734 -> 639

# save sequences before renaming taxa
taxa_seqs = colnames(taxa)

# open taxonomic bootstraps
bootout = readRDS("data/bootout_edit.rds")
boot = bootout[[taxadb]]
boot = data.frame(boot[taxa_seqs,])

# some alternative naming schemes for plotting
full_names <- unname(namedlist[taxa_seqs])
taxa_split = split_by_character(full_names,";")
#PR2 taxonomic levels
if(tdb=="PR2") colnames(taxa_split) = c("Domain","Superphylum","Phylum","Class","Order","Family","Genus","Species","ESV")
#SILVA taxonomic levels
if(tdb=="SILVA") colnames(taxa_split) = c("root","Domain","Major Clade","Superkingdom","Kingdom","Subkingdom","Infrakingdom","Superphylum","Phylum","Subphylum","Infraphylum","Superclass","Class","Subclass","Infraclass","Superorder","Order","Suborder","Superfamily","Family","Subfamily","Genus","Accession","ESV")

esvrank = as.numeric(gsub(pattern = "ESV_",replacement="",x=taxa_split$ESV))

short_names = paste0(taxa_split$Genus, " [", boot$Genus, "] ",taxa_split$ESV)
names(short_names) = taxa_split$ESV
long_names = paste0(taxa_split$Phylum, " ", taxa_split$Class, " ", taxa_split$Order, " ", taxa_split$Family, " ", taxa_split$Genus, " [", boot$Genus, "] ",taxa_split$ESV)
names(long_names) = taxa_split$ESV
barplot_names = paste0(taxa_split$Family, " ", taxa_split$Genus, " [", boot$Genus, "] ",taxa_split$ESV)
names(long_names) = taxa_split$ESV

colnames(taxa) = taxa_split$ESV
rownames(boot) = taxa_split$ESV


### produce taxonomic quality-controlled taxa table
tcutoff = 60 #cutoff for taxonomic id
taxa_good = taxa_split

# replace taxa lacking bootstrap support with higher taxonomic level
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
sub.meta = meta[meta$station_type=="S" & meta$filter==0.2 & meta$project != "AMBON 2015",]

# aggregate by higher taxonomic group for sorting
#don't remove any columns from taxa into sub.taxa or things will break downstream
sub.taxa = taxa[rownames(sub.meta),]
sub.agg = agg.mat(sub.taxa,"Phylum")

sub.meta = data.frame(sub.meta,sub.agg,check.names=F)


## filter taxa function
filter_taxa = function(search_terms, x.taxa, x.meta, x.names) {
  
  search_found = sapply(search_terms, function(x) grepl(x,x.names))
  x.cols = rowSums(search_found)>0
  
  x.taxa = x.taxa[, x.cols]
  x.taxa = x.taxa[rowSums(x.taxa)>0,]

  x.meta = x.meta[rownames(x.taxa),]
  x.taxa = x.taxa[,colSums(x.taxa)>0] #remove empty taxa

  #proportions within group
  x.prop = as.data.frame(prop.table(as.matrix(x.taxa),margin=1))
  
  return(list("meta"=x.meta,"taxa"=x.taxa, "prop"=x.prop))
}

mat.list = list()

# #filter diatoms
mat.list[["diatoms"]] = filter_taxa("Diatom",sub.taxa, sub.meta, full_names)
mat.list[["diatoms"]][["title"]] = "Diatom"
mat.list[["diatoms"]][["short"]] = "D"

# ## filter picos
# #Chlorophyta, Haptophyta, or Chrysophyceae
mat.list[["picos"]] = filter_taxa(c("Haptophyta","Chlorophyta","Chrysophy"),sub.taxa, sub.meta, full_names)
mat.list[["picos"]][["title"]] = "Picoeukaryote"
mat.list[["picos"]][["short"]] = "P"

# ## filter flagellates (dinos + MAST + Syndiniales)
mat.list[["dinos"]] = filter_taxa(c("Dinoflagellata","MAST","Protalveolata"),sub.taxa, sub.meta, full_names)
mat.list[["dinos"]][["title"]] = "Dinoflagellate"
mat.list[["dinos"]][["short"]] = "E"

# ## filter fungi (Opisthokonta) and amoeba (Amoebozoa)
mat.list[["fungi"]] = filter_taxa(c("Opisthokonta","Amoebozoa"),sub.taxa, sub.meta, full_names)
mat.list[["fungi"]][["title"]] = "Fungi"
mat.list[["fungi"]][["short"]] = "F"

# ## filter ciliates
mat.list[["ciliates"]] = filter_taxa("Ciliophora",sub.taxa, sub.meta, full_names)
mat.list[["ciliates"]][["title"]] = "Ciliate"
mat.list[["ciliates"]][["short"]] = "C"

# ## filter bacteria and archaea
mat.list[["bacteria"]] = filter_taxa(c("Bacteria;","Archaea;"),sub.taxa, sub.meta, full_names)
mat.list[["bacteria"]][["title"]] = "Bacteria and Archaea"
mat.list[["bacteria"]][["short"]] = "B"

# ## filter metazoa
mat.list[["metazoa"]] = filter_taxa(c("Metazoa"),sub.taxa, sub.meta, full_names)
mat.list[["metazoa"]][["title"]] = "Metazoa"
mat.list[["metazoa"]][["short"]] = "M"


# ## filter others
o.names = setdiff(colnames(taxa), unname(unlist(lapply(mat.list, function(x) colnames(x[["taxa"]])))))
mat.list[["other"]] = filter_taxa("",sub.taxa[,o.names], sub.meta, full_names[match(o.names,taxa_split$ESV)])
mat.list[["other"]][["title"]] = "Other"
mat.list[["other"]][["short"]] = "O"


######################




source("bin/02_taxa_prop_tables.R")

