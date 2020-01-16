## plot figure taxonomic barplots by transect and depth

for(tax in names(mat.list)) {
  for(transect in c("DBO3","DBO4","DBO5")) {
  
  tx = mat.list[[tax]][["title"]]
  abb = mat.list[[tax]][["short"]]
  x.meta = mat.list[[tax]][["meta"]]
  x.prop = mat.list[[tax]][["prop"]]
  nclus = clus.list[[tax]]
  clusnames = paste0(abb,1:nclus)

  cutoff = 0.95
  
  y.meta = x.meta[grepl(transect,x.meta$station),]
  y.prop = x.prop[rownames(y.meta),]
  
  # pick the cumulative top % of ESVs from each taxa across transects
  y.here = y.prop[,order(colSums(y.prop),decreasing=TRUE)]
  y.pick = which( cumsum(colSums(y.here)/sum(y.here)) < cutoff)

  y.pick = y.pick[order(colSums(y.prop[,y.pick]),decreasing=TRUE)]
  y.pick = names(y.pick)
  y.pick.prop <- y.prop[, y.pick]
  y.other.prop <- y.prop[, setdiff(colnames(y.prop),y.pick)]
  
  # for sorting by taxonomy to make nice colors
  tmptax = taxa_split[match(y.pick, taxa_split$ESV),]
  rownames(tmptax) = y.pick
  tmptax$abund = colMeans(y.pick.prop)
  tmptax$esvrank = as.numeric(gsub("ESV_","",tmptax$ESV))
  
  tmptax = tmptax[
    order(tmptax$Domain, tmptax$Kingdom, tmptax$Phylum, tmptax$Class, tmptax$Family, tmptax$Genus, tmptax$abund, tmptax$esvrank,
          decreasing=c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE),method="radix")
    ,]

  othername = paste0("Other (", sprintf("%.1f",mean(rowSums(y.other.prop)*1000),1),")")
  tmptax[othername,] = c(rep("Other",times=ncol(taxa_split)), 0, 0)

  #names to plot with, including abundances
  
  y.pick.prop[,othername] = rowSums(y.other.prop)
  y.pick.prop = y.pick.prop[,rownames(tmptax)]
  
  colnames(y.pick.prop) = paste0(
    barplot_names[match(colnames(y.pick.prop),taxa_split$ESV)], " (", 
    sprintf("%.1f",round(colMeans(y.pick.prop)*1000,1)), ")")
  
  colnames(y.pick.prop)[ncol(y.pick.prop)] = othername
  
  y.meta$sample = rownames(y.meta)
  
  # reshape data frame for ggplot
  df = data.frame(y.pick.prop,y.meta[,c('depth_type','station','project')],check.names=FALSE) #taxa
  df_long = reshape2::melt(df, id.vars=c('depth_type','station','project'), variable.names='taxa', na.rm=TRUE)
  
  # construct barplot
  bplot <- ggplot(df_long, aes(fill = variable, y = value, x = depth_type, label=variable)) +
    geom_bar(stat='identity', position=position_fill(reverse=T)) +
    ylab('Relative Abundance') + xlab('Depth Bin') +
    theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14)) +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y = element_text(size=12)) + coord_flip() +
    theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12), strip.text = element_text(size = 12)) + 
    ggtitle(paste0(tx, " along ", transect)) + 
    # facet_wrap(vars(project,station), scales='fixed', nrow=1) +
    facet_grid(project ~ station, scales='fixed') +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_blank()) +
    theme(legend.text=element_text(size=8), legend.title=element_blank(), legend.position = 'none')
  
  # this colors the bars and legend by taxonomic level (e.g. 'Genus', 'Family')
  bplot = bplot + scale_fill_manual("legend",values = getcolors(tmptax,"Genus"))
  #bplot
  ggsave(plot = bplot, filename = paste0('output/92_figure_barplot_',tax,"_",transect,'.png'),device = png(),width = 12, height=4)
  
  ppl = bplot + theme(legend.position="bottom") 
  ppl_legend <- ggpubr::get_legend(ppl)
  ppleg = ggpubr::as_ggplot(ppl_legend)
  ppleg
  ggsave(ppleg, filename = paste0('output/92_figure_barplot_',tax,'_',transect,'_leg.png'),device = png(),width = 18, height=8,dpi = 300)



  } #end transect

  graphics.off()

} #end taxa
