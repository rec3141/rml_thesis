## Run CCAs

# require(MASS)
# require(plyr)
# require(tidyverse)
# require(reshape2)
# require(gplots)
# require(ggplot2)
# require(compositions)
# require(vegan)


# metadata
# metavar = c("depth_m","lat","lon","temp","salinity","DO",
#             "FlECO-AFL(mg/m^3)","chl (ug/l)","phaeo (ug/l)","PO4(uM)","Sil(uM)",
#             "NH4(uM)","N+N (umol/L)","transmission_pct")
# metavarlog = c("FlECO-AFL(mg/m^3)","chl (ug/l)","phaeo (ug/l)","PO4(uM)","Sil(uM)",
#                "NH4(uM)","N+N (umol/L)")

#################
# do CCAs on ESVs
#################

# DIATOM ESVs

diatoms.cca1 <- cca(d.prop ~ temp + salinity +`N+N (umol/L)` + `FlECO-AFL(mg/m^3)` + DO + Dinoflagellata + Ciliophora + Protalveolata + Chytridiomycota, data=d.meta, na.action = na.exclude)
diatoms.cca1
summary(diatoms.cca1)
plot(diatoms.cca1,display=c("sp","wa","lc","bp","cn"))

diatoms.anova = anova.cca(diatoms.cca1,permutations = 999,by = "margin",parallel = 8)

# Partitioning of scaled Chi-square:
#   Inertia Proportion
# Total           13.00     1.0000
# Constrained      1.48     0.1138
# Unconstrained   11.52     0.8862
# 
# Model: cca(formula = d.taxa ~ temp + salinity + `N+N (umol/L)` + `FlECO-AFL(mg/m^3)` + Dinoflagellata + Ciliophora + Protalveolata + Chytridiomycota, data = d.meta, na.action = na.exclude)
# Df ChiSquare       F Pr(>F)    
# temp                  1    0.2061  8.8562  0.001 ***
#   salinity              1    0.2031  8.7233  0.001 ***
#   `N+N (umol/L)`        1    0.1118  4.8048  0.001 ***
#   `FlECO-AFL(mg/m^3)`   1    0.3483 14.9652  0.001 ***
#   Dinoflagellata        1    0.0821  3.5260  0.001 ***
#   Ciliophora            1    0.0880  3.7813  0.001 ***
#   Protalveolata         1    0.1039  4.4616  0.001 ***
#   Chytridiomycota       1    0.1618  6.9523  0.001 ***
#   Residual            495   11.5220            


diatoms.cca2 <- cca(d.prop ~ mass, data=d.meta, na.action = na.exclude)
diatoms.cca2
plot(diatoms.cca2,display=c("sp","wa","lc","bp","cn"))
diatoms.anova2 = anova.cca(diatoms.cca2,permutations = 999,by = "margin",parallel = 8)

# Inertia Proportion Rank
# Total         13.04067    1.00000     
# Constrained    0.62766    0.04813    6
# Unconstrained 12.41301    0.95187  215
# Model: cca(formula = d.taxa ~ mass, data = d.meta, na.action = na.exclude)
# Df ChiSquare      F Pr(>F)    
# mass       6    0.6277 4.5677  0.001 ***
#   Residual 542   12.4130      

### PICOEUK ESVs

picos.cca1 <- cca(p.prop ~ temp + salinity +`N+N (umol/L)` + `FlECO-AFL(mg/m^3)` + Dinoflagellata + Ciliophora + Protalveolata + Chytridiomycota, data=p.meta, na.action = na.exclude)
picos.cca1
summary(picos.cca1)
plot(picos.cca1,display=c("sp","wa","lc","bp","cn"))

picos.anova = anova.cca(picos.cca1,permutations = 999,by = "margin",parallel = 8)

# Inertia Proportion Rank
# Total         18.21488    1.00000     
# Constrained    1.03241    0.05668    8
# Unconstrained 17.18247    0.94332  234
# Model: cca(formula = p.taxa ~ temp + salinity + `N+N (umol/L)` + `FlECO-AFL(mg/m^3)` + Dinoflagellata + Ciliophora + Protalveolata + Chytridiomycota, data = p.meta, na.action = na.exclude)
# Df ChiSquare      F Pr(>F)    
# temp                  1    0.2259 6.2190  0.001 ***
#   salinity              1    0.1052 2.8955  0.002 ** 
#   `N+N (umol/L)`        1    0.0641 1.7658  0.005 ** 
#   `FlECO-AFL(mg/m^3)`   1    0.0617 1.6986  0.091 .  
# Dinoflagellata        1    0.0646 1.7775  0.008 ** 
#   Ciliophora            1    0.0693 1.9090  0.015 *  
#   Protalveolata         1    0.1015 2.7929  0.005 ** 
#   Chytridiomycota       1    0.1161 3.1970  0.002 ** 
#   Residual            473   17.1825       


picos.cca2 <- cca(p.prop ~ mass, data=p.meta, na.action = na.exclude)
picos.cca2
plot(picos.cca2,display=c("sp","wa","lc","bp","cn"))
picos.anova2 = anova.cca(picos.cca2,permutations = 999,by = "margin",parallel = 8)

# Inertia Proportion Rank
# Total         20.0604     1.0000     
# Constrained    0.5496     0.0274    6
# Unconstrained 19.5109     0.9726  240

# Model: cca(formula = p.taxa ~ mass, data = p.meta, na.action = na.exclude)
# Df ChiSquare      F Pr(>F)  
# mass       6    0.5496 2.4225  0.017 *
#   Residual 516   19.5108              


## DIATOM GENERA
dgenus.taxa = agg.mat(d.taxa,"Genus")

dgenus.cca1 <- cca(dgenus.taxa ~ temp + salinity +`N+N (umol/L)` + `FlECO-AFL(mg/m^3)` + Dinoflagellata + Ciliophora + Protalveolata + Chytridiomycota, data=d.meta, na.action = na.exclude)
dgenus.cca1
summary(dgenus.cca1)
plot(dgenus.cca1,display=c("sp","wa","lc","bp","cn"))

anova.cca(dgenus.cca1,permutations = 999,by = "margin",parallel = 8)

# Inertia Proportion Rank
# Total          3.8677     1.0000     
# Constrained    0.4961     0.1283    8
# Unconstrained  3.3716     0.8717   39

# Model: cca(formula = dgenus.taxa ~ temp + salinity + `N+N (umol/L)` + `FlECO-AFL(mg/m^3)` + Dinoflagellata + Ciliophora + Protalveolata + Chytridiomycota, data = d.meta, na.action = na.exclude)
# Df ChiSquare       F Pr(>F)    
# temp                  1    0.0892 13.0901  0.001 ***
#   salinity              1    0.0378  5.5540  0.001 ***
#   `N+N (umol/L)`        1    0.0274  4.0246  0.001 ***
#   `FlECO-AFL(mg/m^3)`   1    0.1405 20.6349  0.001 ***
#   Dinoflagellata        1    0.0227  3.3332  0.001 ***
#   Ciliophora            1    0.0250  3.6687  0.002 ** 
#   Protalveolata         1    0.0160  2.3492  0.016 *  
#   Chytridiomycota       1    0.0920 13.5001  0.001 ***
#   Residual            495    3.3716        


dgenus.cca2 <- cca(dgenus.taxa ~ mass, data=d.meta, na.action = na.exclude)
dgenus.cca2
plot(dgenus.cca2,display=c("sp","wa","lc","bp","cn"))
anova.cca(dgenus.cca2,permutations = 999,by = "margin",parallel = 8)

# Inertia Proportion Rank
# Total         3.83052    1.00000     
# Constrained   0.14732    0.03846    6
# Unconstrained 3.68321    0.96154   39

# Model: cca(formula = dgenus.taxa ~ mass, data = d.meta, na.action = na.exclude)
# Df ChiSquare      F Pr(>F)    
# mass       6    0.1473 3.6131  0.001 ***
#   Residual 542    3.6832     



## PICOEUK GENERA
pgenus.taxa = agg.mat(p.taxa,"Genus")

pgenus.cca1 <- cca(pgenus.taxa ~ temp + salinity +`N+N (umol/L)` + `FlECO-AFL(mg/m^3)` + Dinoflagellata + Ciliophora + Protalveolata + Chytridiomycota, data=p.meta, na.action = na.exclude)
pgenus.cca1
summary(pgenus.cca1)
plot(pgenus.cca1,display=c("sp","wa","lc","bp","cn"))

anova.cca(pgenus.cca1,permutations = 999,by = "margin",parallel = 8)

# Inertia Proportion Rank
# Total          5.9061     1.0000     
# Constrained    0.6377     0.1080    8
# Unconstrained  5.2684     0.8920   38

# Model: cca(formula = pgenus.taxa ~ temp + salinity + `N+N (umol/L)` + `FlECO-AFL(mg/m^3)` + Dinoflagellata + Ciliophora + Protalveolata + Chytridiomycota, data = p.meta, na.action = na.exclude)
# Df ChiSquare      F Pr(>F)    
# temp                  1    0.1031 9.2569  0.001 ***
#   salinity              1    0.0780 7.0041  0.001 ***
#   `N+N (umol/L)`        1    0.0287 2.5768  0.003 ** 
#   `FlECO-AFL(mg/m^3)`   1    0.0348 3.1289  0.020 *  
#   Dinoflagellata        1    0.0382 3.4337  0.001 ***
#   Ciliophora            1    0.0448 4.0197  0.001 ***
#   Protalveolata         1    0.0410 3.6832  0.001 ***
#   Chytridiomycota       1    0.0852 7.6521  0.001 ***
#   Residual            473    5.2684      
# 

pgenus.cca2 <- cca(pgenus.taxa ~ mass, data=p.meta, na.action = na.exclude)
pgenus.cca2
plot(pgenus.cca2,display=c("sp","wa","lc","bp","cn"))
anova.cca(pgenus.cca2,permutations = 999,by = "margin",parallel = 8)

# Inertia Proportion Rank
# Total         6.02821    1.00000     
# Constrained   0.30162    0.05003    6
# Unconstrained 5.72659    0.94997   38
# 
# Model: cca(formula = pgenus.taxa ~ mass, data = p.meta, na.action = na.exclude)
# Df ChiSquare      F Pr(>F)    
# mass       6    0.3016 4.5296  0.001 ***
#   Residual 516    5.7266            



## DIATOM FAMILY
dfam.taxa = agg.mat(d.taxa,"Family")

dfam.cca1 <- cca(dfam.taxa ~ temp + salinity +`N+N (umol/L)` + `FlECO-AFL(mg/m^3)` + Dinoflagellata + Ciliophora + Protalveolata + Chytridiomycota, data=d.meta, na.action = na.exclude)
dfam.cca1
summary(dfam.cca1)
plot(dfam.cca1,display=c("sp","wa","lc","bp","cn"))

anova.cca(dfam.cca1,permutations = 999,by = "margin",parallel = 8)

# Inertia Proportion Rank
# Total         0.58975    1.00000     
# Constrained   0.07979    0.13530    4
# Unconstrained 0.50996    0.86470    4
# 
# Model: cca(formula = dfam.taxa ~ temp + salinity + `N+N (umol/L)` + `FlECO-AFL(mg/m^3)` + Dinoflagellata + Ciliophora + Protalveolata + Chytridiomycota, data = d.meta, na.action = na.exclude)
# Df ChiSquare       F Pr(>F)    
# temp                  1   0.02252 21.8602  0.001 ***
#   salinity              1   0.00928  9.0098  0.009 ** 
#   `N+N (umol/L)`        1   0.01036 10.0601  0.002 ** 
#   `FlECO-AFL(mg/m^3)`   1   0.00134  1.2994  0.140    
# Dinoflagellata        1   0.00566  5.4902  0.004 ** 
#   Ciliophora            1   0.00220  2.1360  0.113    
# Protalveolata         1   0.00290  2.8130  0.060 .  
# Chytridiomycota       1   0.01947 18.9003  0.001 ***
#   Residual            495   0.50996 


dfam.cca2 <- cca(dfam.taxa ~ mass, data=d.meta, na.action = na.exclude)
dfam.cca2
plot(dfam.cca2,display=c("sp","wa","lc","bp","cn"))
anova.cca(dfam.cca2,permutations = 999,by = "margin",parallel = 8)

# Inertia Proportion Rank
# Total         0.58604    1.00000     
# Constrained   0.02303    0.03930    4
# Unconstrained 0.56301    0.96070    4
# 
# Model: cca(formula = dfam.taxa ~ mass, data = d.meta, na.action = na.exclude)
# Df ChiSquare      F Pr(>F)  
# mass       6   0.02303 3.6953  0.021 *
#   Residual 542   0.56301                
# 



## PICOEUK FAMILY
pfam.taxa = agg.mat(p.taxa,"Family")

pfam.cca1 <- cca(pfam.taxa ~ temp + salinity +`N+N (umol/L)` + `FlECO-AFL(mg/m^3)` + Dinoflagellata + Ciliophora + Protalveolata + Chytridiomycota, data=p.meta, na.action = na.exclude)
pfam.cca1
summary(pfam.cca1)
plot(pfam.cca1,display=c("sp","wa","lc","bp","cn"))

anova.cca(dfam.cca1,permutations = 999,by = "margin",parallel = 8)

# Inertia Proportion Rank
# Total          3.3813     1.0000     
# Constrained    0.4925     0.1457    8
# Unconstrained  2.8888     0.8543   15
# 
# Number of permutations: 999
# 
# Model: cca(formula = dfam.taxa ~ temp + salinity + `N+N (umol/L)` + `FlECO-AFL(mg/m^3)` + Dinoflagellata + Ciliophora + Protalveolata + Chytridiomycota, data = d.meta, na.action = na.exclude)
# Df ChiSquare       F Pr(>F)    
# temp                  1   0.02252 21.8602  0.001 ***
#   salinity              1   0.00928  9.0098  0.013 *  
#   `N+N (umol/L)`        1   0.01036 10.0601  0.001 ***
#   `FlECO-AFL(mg/m^3)`   1   0.00134  1.2994  0.153    
# Dinoflagellata        1   0.00566  5.4902  0.010 ** 
#   Ciliophora            1   0.00220  2.1360  0.113    
# Protalveolata         1   0.00290  2.8130  0.064 .  
# Chytridiomycota       1   0.01947 18.9003  0.001 ***
#   Residual            495   0.50996  









# stuff to use for plotting?
# taxa.pca <- prcomp(as.matrix(unclass(clr(diatoms))))
# biplot(taxa.pca,cex=0.3,)
# 
# 
# 
# # x <- 'NULL'
# xlim=c(-4,4)
# ylim=c(-6,6)
# scl=3
# layout(matrix(1:2, ncol = 2))
# #op <- par(mar = c(5,4,4,1) + 0.1)
# ## site/sample scores
# plot(taxa.cca1, xlim=xlim, ylim=ylim, type="n", scaling=3, main=expression('Taxonomic Groups'), cex=1)
# orditorp(taxa.cca1, display="species", labels=new_name[1:50],scaling=scl,
#          col = "forestgreen", cex = 1, pch = 2, air = 2)
# text(taxa.cca1, scaling = 3, display = "bp")
# ## Species scores
# plot(taxa.cca1, xlim=xlim, ylim=ylim,type = "n", scaling = 3, main = expression('Samples'), cex = 1)
# orditorp(taxa.cca1, display = "sites", labels=x, scaling = scl,
#          col = "darkorchid4", pch = 2, cex = 1, air = 0.5)
# text(taxa.cca1, scaling = 3, display = "bp")
# #par(op)
# layout(1)
# 
# 
# 
# xlim=c(-4,4)
# ylim=c(-6,6)
# scl=3
# layout(matrix(1:2, ncol = 2))
# #op <- par(mar = c(5,4,4,1) + 0.1)
# ## site/sample scores
# plot(taxa.cca1, xlim=xlim, ylim=ylim, type="n", scaling=3, main=expression('Taxonomic Groups'), cex=1)
# orditorp(taxa.cca1, display="species", labels=new_name[1:50],scaling=scl,
#          col = "forestgreen", cex = 1, pch = 2, air = 2)
# text(taxa.cca1, scaling = 3, display = "bp")
# ## Species scores
# plot(taxa.cca1, xlim=xlim, ylim=ylim,type = "n", scaling = 3, main = expression('Samples'), cex = 1)
# orditorp(taxa.cca1, display = "sites", labels=x, scaling = scl,
#          col = "darkorchid4", pch = 2, cex = 1, air = 0.5)
# text(taxa.cca1, scaling = 3, display = "bp")
# #par(op)
# layout(1)
# 
# 
# 
# taxa.plot1 <- plot(taxa.cca1, scaling=1)
