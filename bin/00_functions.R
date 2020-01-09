#install required packages and dependencies
#install.packages(c("gplots","vegan","dendextend","phytools","heatmap.plus","ggtree","ggstance","colorspace","ggpubr","compositions","grid","egg","ggplot2","oce","viridis","PBSmapping","marmap","scales"),dependencies=TRUE)

# function to calculate stratification
# requires data frame with rho, z, bottom
stratify <- function(y) {
  if(nrow(y)<3) return(NA)
  (y$rho[match(max(y$z),y$z)] - y$rho[match(min(y$z),y$z)])/y$bottom[1]
}

#function to calculate water mass
get_water_mass = function(x) {
  mass = x$salinity + NA
  
  mass[x$salinity>29.7 & x$salinity<33.64 & x$temp>4.5 & x$temp<15] <- 'ACW'
  mass[x$salinity>30.5 & x$salinity<33.64 & x$temp>3 & x$temp<15] <- 'ACW'
  
  mass[x$salinity>24 & x$salinity<29.7 & x$temp>4.5 & x$temp<15] <- 'FACW'
  
  mass[x$salinity>32 & x$salinity<33.64 & x$temp<3 & x$temp>0] <- 'BSW'
  
  mass[x$salinity>30.5 & x$salinity<33.64 & x$temp<0 & x$temp>-1.8] <- 'RWW'
  
  mass[x$salinity>24 & x$salinity<30.5 & x$temp<2 & x$temp>-1.6] <- 'MWR'
  mass[x$salinity>29 & x$salinity<30.5 & x$temp>2 & x$temp<3] <- 'MWR'
  mass[x$salinity>30.5 & x$salinity<32 & x$temp<3 & x$temp>0] <- 'MWR'
  
  mass[x$salinity>33.64 & x$salinity<35 & x$temp > -1.26 & x$temp < 1] <- 'AW'
  
  mass[x$salinity>29 & x$salinity<30.5 & x$temp > 3 & x$temp < 5] <- 'SCW'
  mass[x$salinity>24 & x$salinity<29 & x$temp > 2 & x$temp < 5] <- 'SCW'
  
  return(mass)
}


# get default ggplot color scale
gg_color_hue <- function(n) {
  hues=seq(15, 375, length=n+1)
  hcl(h = hues, l=65, c=100)[1:n]
}

#split a column on a character
split_by_character = function(x,char) {
  as.data.frame(read.table(text=x, sep=char, as.is=TRUE, fill=TRUE))
}

#aggregate a matrix at a specific taxonomic level
agg.mat = function(x,taxlevel) {
  x.agg = aggregate(t(x),by=list(taxa_split[match(colnames(x),short_names),taxlevel]),FUN=sum)
  rownames(x.agg) = x.agg[,1]
  x.agg = x.agg[,-1]
  x.agg = t(x.agg)
  x.agg = data.frame(x.agg)
  return(x.agg)  
}

# get related colors per taxonomic level
getcolors = function(mat,col) {
  uniqfams = rle(as.character(mat[,col]))$values
  famlen = rle(as.character(mat[,col]))$lengths
  
  h = floor(360/length(uniqfams))
  famcols = NULL
  
  for(i in 1:length(uniqfams)) {
    if(uniqfams[i]=="Other") {
      famcols = c(famcols, "black")
    } else {
      famcols = c(famcols, rev(sequential_hcl(famlen[i]+1,h=i*h,power=0.5)[1:(famlen[i])]))
    }
  }
  famcols
}
