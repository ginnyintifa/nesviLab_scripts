### PTM normalized by protein abundance 


ptmSite_proteomeNormalization = function(tidy_phospho_filename, 
                                         tidy_proteome_filename,
                                         output_filename)
{
  
  # tidy_proteome_filename = "/Users/ginny/My Drive/nesviLab_scripts/test_output/ratio_gene_MD_tidy_mapFile_proteogenomics.tsv"
  # tidy_phospho_filename = "/Users/ginny/My Drive/nesviLab_scripts/test_output/ratio_signle-site_MD_tidy_mapFile_proteogenomics.tsv"
  
  
  
  prot_df = fread(tidy_proteome_filename,
                  stringsAsFactors = F, data.table = F)
  psite_df = fread(tidy_phospho_filename,
                   stringsAsFactors = F, data.table = F)
  
  
  pprot = gsub("_.*", "", psite_df$Gene_site)
  inter_prot = intersect(pprot, prot_df$Gene_name)
  psite_df = psite_df%>%
    dplyr::mutate(pprot)%>%
    dplyr::filter(pprot %in% inter_prot)%>%
    dplyr::select(pprot, Gene_site, SequenceWindow, everything())

  inter_sample = intersect(colnames(psite_df), colnames(prot_df))
  
  psite_df = psite_df[,c("pprot", "Gene_site", "SequenceWindow", inter_sample)]
  prot_df = prot_df[,c("Gene_name", inter_sample)]
  
  prot_sort_df = data.frame(Gene_name = inter_prot, stringsAsFactors = F)%>%
    dplyr::left_join(prot_df, by = "Gene_name")
  
  tl = nrow(psite_df)* length(inter_sample)
  all_psite = rep(0, tl)
  all_prot = rep(0, tl)
  record_pos_start = rep(0, length(inter_prot))
  record_pos_end = rep(0, length(inter_prot))
  prot1 = inter_prot[1]
  sel_phos1 = psite_df%>%
    dplyr::filter(pprot == prot1)
  
  phos_d = c(t(as.matrix(sel_phos1[,inter_sample])))
  prot_d = rep(c(as.matrix(prot_sort_df[1, inter_sample])), nrow(sel_phos1))
  
  record_pos_start[1] = 1
  record_pos_end[1] = length(phos_d)
  all_psite[record_pos_start[1]: record_pos_end[1]] = phos_d
  all_prot[record_pos_start[1]: record_pos_end[1]] = prot_d
  
  for(i in 2:length(inter_prot))
  {
    prot = inter_prot[i]
    sel_phos = psite_df%>%
      dplyr::filter(pprot == prot)
    phos_d = c(t(as.matrix(sel_phos[,inter_sample])))
    prot_d = rep(c(as.matrix(prot_sort_df[i, inter_sample])), nrow(sel_phos))
    record_pos_start[i] = record_pos_end[i-1]+1
    record_pos_end[i] = record_pos_start[i]-1+ length(phos_d)
    all_psite[record_pos_start[i]: record_pos_end[i]] = phos_d
    all_prot[record_pos_start[i]: record_pos_end[i]] = prot_d
    if(i%%1000 ==0)
      cat(i, "\n")
  }
  
  
  rm(prot_df)
  
  p_df = data.frame(phospho = all_psite, prot = all_prot,
                    caseID = rep(inter_sample, nrow(psite_df)),
                    geneSite = rep(psite_df$Gene_site, each = length(inter_sample)),
                    sequenceWindow = rep(psite_df$SequenceWindow, each = length(inter_sample)),
                    stringsAsFactors = F)%>%
    na.omit()
  
  rm(psite_df)
  rm(all_psite)
  rm(all_prot)
  
  p = lm(phospho~prot, data = p_df)
  
  # 
  # Call:
  #   lm(formula = phospho ~ prot, data = p_df)
  # 
  # Coefficients:
  #   (Intercept)         prot  
  # -0.02729      0.70707  
  # 

  res = p$residuals
  p_df = p_df%>%
    dplyr::mutate(subPsite = res)
  
  rm(res)
  
  result_cor = cor(p_df$subPsite, p_df$prot)
  cat("correlation with prot after normalization: ",  result_cor, "\n")
  
 # -4.872899e-16
  
  sites_windows = p_df%>%
    dplyr::select(geneSite, sequenceWindow)%>%
    unique()
  
  
    
  cat("get sites_windows", "\n")
  
  rm(p)
  
  dm = matrix(NA, nrow = nrow(sites_windows), ncol = length(inter_sample))
  
  for(i in 1:length(inter_sample))
  {
    this_name = inter_sample[i]
    t = p_df%>%
      dplyr::filter(caseID == this_name)
    
    m_df = sites_windows%>%
      dplyr::left_join(t, by = c("geneSite", "sequenceWindow"))

    dm[,i] = m_df$subPsite
  }
  cat("get dm", "\n")
  rm(p_df)
  
  subpsite = data.frame(sites_windows,  dm, stringsAsFactors = F)
  rm(dm)
  colnames(subpsite) = c("Gene_site", "SequenceWindow", inter_sample)
  cat("get subpsite", "\n")
  write.table(subpsite, output_filename,
              quote = F, row.names = F, sep = "\t")
  
}