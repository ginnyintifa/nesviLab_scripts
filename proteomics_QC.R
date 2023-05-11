

sumna = function(x)
{
  s = sum(is.na(x))
  return(s)
}



plexes = sort(rep(c(1:13), 17))



bar_number = function(prot_data, prot_annot, tag, pdf_name)
{
   
  
  nplex = length(unique(prot_annot$plex))
  
  prot_nu = rep(0, nplex)
  
  for(i in 1:nplex)
  {
    
    cols = prot_annot$caseID[which(prot_annot$plex == i)]
    
    d = prot_data[, cols]
    sumnad = apply(d, 1, sumna)
    notAllNA = which(sumnad<length(cols))
    
    prot_nu[i] = length(notAllNA)
    
    
  }
  prot_nu_med = median(prot_nu)
  
  
  pdf(pdf_name, useDingbats = F,
      height = 3)
  barplot(prot_nu, xlab = "Plex", ylab = tag,
          names.arg = c(1:nplex),las =1,
          #col = c(rep("steelblue",7), rep("grey", 23)),
          main = paste0("Median: ", prot_nu_med))
  dev.off()
}






get_density = function(prot_data, 
                       tag,
                       pdf_output_filename)
{
  
  
  is.na(prot_data) <- sapply(prot_data, is.infinite)                  
  
  
  pdf(pdf_output_filename, useDingbats = F, height = 3)
  
  ### always need to deal with  -Inf, it might be a result of log2 on 0 
  
  ym = c()
  
  for(i in 1:ncol(prot_data))
  {
    
    d = density(prot_data[, i], na.rm = T)
    
    ym = c(ym , max(d$y))
  }
  
  
  
  
  d = density(prot_data[, 1], na.rm = T)
  plot(density(prot_data[, 1], na.rm = T), ylim = c(0,max(ym)),xlim = c(min(prot_data, na.rm = T), max(prot_data, na.rm = T)), col = "grey",
       main = paste0(tag," density"))
  
  for(i in 1:ncol(prot_data))
  {
    
    this_dens = density(prot_data[, i], na.rm = T)
    lines(this_dens, col = "grey")
  }
  
  
  dev.off()
  
  
  
}




plex_pca_fviz_area  = function(prot_data, prot_annot, pdf_name,
                                   xmina, xmaxa, ymina, ymaxa)
{
  # 
  # # 
  # prot_d =  proteome_data
  # prot_annot = proteomic_annot
  # pdf_name = "/Users/ginny/My Drive/nesviLab_scripts/test_output/proteome_ratio_gene_plex_pca.pdf"
  # xmina = 20
  # xmaxa = 20
  # ymina = 20
  # ymaxa = 10
  # # 
  # 
  
  prot_data = prot_data[, prot_annot$caseID]
  
  
  p_na = apply(prot_data, 1, sumna)
  
  p_d= prot_data[which(p_na ==0),]
  
  nf = nrow(p_d)
  

  pd_pca = prcomp(t(p_d))
  
  
  plexes = prot_annot$plex
  
  
  
  xmin = min(pd_pca$x[,1])-xmina
  xmax = max(pd_pca$x[,1])+xmaxa
  
  ymin = min(pd_pca$x[,2])-ymina
  ymax = max(pd_pca$x[,2])+ymaxa
  
  
  p = fviz(pd_pca,
           element = "ind",
           geom="point", 
           pointsize = 2,
           label = "none",
           # pointshape = type_shape,
           # color = type_col,
           
           habillage = plexes,
           addEllipses = T,
           ellipse.level = 0.9,
           ellipse.alpha = 0,
           invisible = "quali")+
    #  scale_color_manual(name = "subtypes", labels = sts,
    #         values= sts_col$color)+
    xlim(xmin, xmax)+
    ylim(ymin, ymax)+
    
    #scale_shape_manual(name = "subtyp", labels = sts,
    #         values= c(1:length(sts)))+
    labs(title = paste0("# features ", nf))+
    theme_classic()
  
  
  pdf(pdf_name, useDingbats = F, width = 9)
  grid.arrange(p, ncol = 1, nrow = 1)
  dev.off()
  
}






generate_cor_plot = function(prot_data, prot_annot, pdf_name, tag)
{
  
  prot_annot = prot_annot%>%
    dplyr::arrange(caseID)
  
  
  rep_col = prot_annot$caseID[which(prot_annot$isDup == T)]
  rep_col_plex = prot_annot$plex[which(prot_annot$isDup == T)]
  
  
  # 
  # srep_col = prot_annot$caseID[which(prot_annot$specialDup == T)]
  # srep_col_plex = prot_annot$plex[which(prot_annot$specialDup == T)]
  # 
  prot_qc = prot_data[, rep_col]
  
  #prot_ref = prot_data$ReferenceIntensity
  #pqc = sweep(prot_qc, 1, prot_ref)
  
  colnames(prot_qc) = paste(rep_col, rep_col_plex, sep = "_")
  
  
  pqc = prot_qc
  
  
  cor_prot_qc = cor(pqc, use = "pairwise.complete.obs")
  
  
  ggcorr(pqc, label = T, label_round = 2, label_size = 2, size = 2.5, hjust = 1, label_color = "black",
         low = "blue", mid = "white", high = "red") +
    ggtitle(paste0(tag, ": mean(QC r) = ", round(mean(cor_prot_qc, na.rm = T), digits = 2)))+
    theme(plot.margin = unit(c(0.5, 0.5, 1, 2), "pt"))
  
  
  #ggsave( "G:/My Drive/CPTAC_AML/MS_QC/raw_proteome_gene_ratio_all_cor.pdf")
  ggsave(pdf_name)
  
  
}


