#install required packages and dependencies
#install.packages(c("gplots","vegan","dendextend","phytools","heatmap.plus","ggtree","ggstance","colorspace","ggpubr","compositions","grid","egg","ggplot2","oce","viridis","PBSmapping","marmap","scales"),dependencies=TRUE)
#BiocManager::install("pkimes/sigclust2")
#devtools::install_github("nolanlab/Rclusterpp") #clang: error: unsupported option '-fopenmp'
#devtools::install_github("vmikk/metagMisc")

# function to calculate stratification
# requires data frame with rho, z, bottom
stratify <- function(y) {
  if(nrow(y)<3) return(NA)
  (y$rho[match(max(y$z),y$z)] - y$rho[match(min(y$z),y$z)])/y$bottom[1]
}

#water mass definitions based on Pisavera et al. 2015
massmat = 'WACW	WACW	WACW	WACW	WACW	WACW
WACW	WACW	ACW	ACW	ACW	ACW
SCW	SCW	SCW	ACW	ACW	ACW
SCW	MW	MW	MW	BSW	BSW
MW	MW	MW	MW	BSW	AW
MW	MW	MW	RWW	RWW	AW
MW	MW	MW	RWW	RWW	AW
MW	MW	MW	WW	WW	WW
'
massin = read.table(text=massmat,sep="\t",stringsAsFactors=FALSE)
massx = c(24.00, 29.00, 29.70, 30.50, 32.00, 33.64, 35.00)
nx = length(massx)-1
#massy = c(-1.80, -1.60, -1.26, 0.00, 2.00, 3.00, 4.50, 6.00, 14.00) # original
massy = c(-1.80, -1.60, -1.26, 0.00, 1.70, 3.00, 5.00, 6.00, 14.00)
ny = length(massy) - 1
masses = matrix(NA, nrow=nx*ny, ncol=4)
massnames = 1:(nx*ny)

for(i in 1:nx) {
  for(j in 1:ny) {
    masses[(i-1)*ny + j,] = c(massx[i], massx[i+1], massy[j],massy[j+1])
    massnames[(i-1)*ny + j] = massin[ny-j+1,i]
  }
}

masses = data.frame('mass'=massnames, masses, stringsAsFactors=FALSE)
colnames(masses) = c("mass","xmin","xmax","ymin","ymax")

#function to calculate water mass
get_water_mass = function(x) {
  mass = rep(NA,ncol(x))
  
  for(i in 1:nrow(masses)) {
    mass[x$salinity >= masses$xmin[i] & x$salinity < masses$xmax[i] & x$temp >= masses$ymin[i] & x$temp < masses$ymax[i]] = masses$mass[i]
  }
  
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
  x.agg = aggregate(t(x),by=list(taxa_split[match(colnames(x),colnames(taxa)),taxlevel]),FUN=sum)
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
