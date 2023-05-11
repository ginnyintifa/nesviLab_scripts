

eb.fit_qvalue <- function(dat, design, contrast, lfc){


  
  n <- dim(dat)[1]
  fit <- lmFit(dat, as.matrix(design))
  
  
  
  fit2 <- contrasts.fit(fit, as.matrix(contrast))
  fit.eb <- eBayes(fit2)
  log2FC <- fit.eb$coefficients[, 1]
  log2FC_ori = lfc
  
  df.r <- fit.eb$df.residual
  df.0 <- rep(fit.eb$df.prior, n)
  s2.0 <- rep(fit.eb$s2.prior, n)
  s2 <- (fit.eb$sigma)^2
  s2.post <- fit.eb$s2.post
  t.ord <- fit.eb$coefficients[, 1]/fit.eb$sigma/fit.eb$stdev.unscaled[, 1]
  t.mod <- fit.eb$t[, 1]
  p.ord <- 2*pt(-abs(t.ord), fit.eb$df.residual)
  p.mod <- fit.eb$p.value[, 1]
  
  #q.ord <- qvalue::qvalue(p.ord)$q
  
  q.ord <- qvalue(p.ord,lambda = seq(0.05,max(p.ord, na.rm = T),0.05))$q
  q.mod <- qvalue(p.mod,lambda = seq(0.05,max(p.mod, na.rm = T),0.05))$q
  
  results.eb <- data.frame(
    log2FC,
    log2FC_ori,
    t.ord,
    t.mod,
    p.ord,
    p.mod,
    q.ord,
    q.mod,
    df.r,
    df.0,
    s2.0,
    s2,
    s2.post
  )
  
  # results.eb <- results.eb[order(results.eb$p.mod), ]
  
  return(results.eb)
}





eb.fit <- function(dat, design, contrast, lfc){
  
  
  n <- dim(dat)[1]
  fit <- lmFit(dat, as.matrix(design))
  
  
  
  fit2 <- contrasts.fit(fit, as.matrix(contrast))
  fit.eb <- eBayes(fit2)
  log2FC <- fit.eb$coefficients[, 1]
  log2FC_ori = lfc
  
  df.r <- fit.eb$df.residual
  df.0 <- rep(fit.eb$df.prior, n)
  s2.0 <- rep(fit.eb$s2.prior, n)
  s2 <- (fit.eb$sigma)^2
  s2.post <- fit.eb$s2.post
  t.ord <- fit.eb$coefficients[, 1]/fit.eb$sigma/fit.eb$stdev.unscaled[, 1]
  t.mod <- fit.eb$t[, 1]
  p.ord <- 2*pt(-abs(t.ord), fit.eb$df.residual)
  p.mod <- fit.eb$p.value[, 1]
  
  #q.ord <- qvalue::qvalue(p.ord)$q   ###### there is a bug with qvalue pacakge 
  
  q.ord <- p.adjust(p.ord, method = "BH")
  q.mod <- p.adjust(p.mod, method = "BH")
  

  results.eb <- data.frame(
    log2FC,
    log2FC_ori,
    t.ord,
    t.mod,
    p.ord,
    p.mod,
    q.ord,
    q.mod,
    df.r,
    df.0,
    s2.0,
    s2,
    s2.post
  )
  
  # results.eb <- results.eb[order(results.eb$p.mod), ]
  
  return(results.eb)
}




sumna_not = function(x)
{
  it = is.na(x)
  return(sum(!it))
}


limma_de_table= function(interest_sample1, 
                         interest_sample2,
                         interest_sample1_invert = F, 
                         interest_sample2_invert  = F, 
                         data_table_file, 
                         col_annot_file,
                         adjust_purity = T,
                         table_output_file)
  
{
  # 
  # interest_sample1 = "wgii_high_nccRCC"
  # interest_sample2 = "wgii_low_nccRCC"
  # interest_sample1_invert = F
  # interest_sample2_invert  = F
  # data_table_file = "/Users/ginny/My Drive/nesviLab_scripts/test_output2/ncc_ratio_gene_MD_tidy_mapFile.tsv"
  # col_annot_file = "/Users/ginny/My Drive/nonCCRCC20211020/clinical/col_annot_meta_20220526_wgii.tsv"
  # adjust_purity = T
  # table_output_file = "/Users/ginny/My Drive/nesviLab_scripts/test_output2/ncc_ratio_gene_MD_DE_high_low_wgii.tsv"
  # 
  
  
  # 
  data_table = fread(data_table_file, stringsAsFactors = F, data.table = F)
  
  col_annot = fread(col_annot_file, stringsAsFactors = F, data.table = F)
  
  all_sample = col_annot$caseID
  p_sample = colnames(data_table[,-1])
  
  inter_sample = intersect(all_sample, p_sample)
  
  
  ### trim data to samples in intersection 
  
  col_annot = col_annot%>%
    dplyr::filter(caseID %in% inter_sample)
  
  wi = which(colnames(data_table) %in% inter_sample)
  
  data_table = data_table[, c(1,wi)]
  
  
  ### retrieve samples of interest
  
  s1 = col_annot%>%
    dplyr::filter(tumorClass %in% interest_sample1)
  
  if(interest_sample1_invert == T)
  {
    
    s1 = col_annot%>%
      dplyr::filter(!tumorClass %in% interest_sample1)
    
  }
  s1_ID = intersect(s1$caseID, colnames(data_table))
  
  s1p = data.frame(caseID = s1_ID, stringsAsFactors = F)%>%
    dplyr::left_join(col_annot, by = "caseID")
  s1_purity = s1p$purity
  
  
  
  s2 = col_annot%>%
    dplyr::filter(tumorClass %in% interest_sample2)
  
  if(interest_sample2_invert == T)
  {
    
    s2 = col_annot%>%
      dplyr::filter(!tumorClass %in% interest_sample2)
    
  }
  s2_ID = intersect(s2$caseID, colnames(data_table))
  
  
  
  s2p = data.frame(caseID = s2_ID, stringsAsFactors = F)%>%
    dplyr::left_join(col_annot, by = "caseID")
  s2_purity = s2p$purity
  
  
  
  
  
  if(length(s1_ID)>0 & length(s2_ID)>0)
  {
    
    group = as.factor(c(rep(1, length(s1_ID)),rep(2, length(s2_ID))))
    
    design = model.matrix( ~ group)
    
    if(adjust_purity == T)
    {
      purity = c(s1_purity, s2_purity)
      
      design = model.matrix( ~  group + purity)
      
    }
    
    
    contrast = makeContrasts(0-group2, levels = colnames(design))
    
    
    
    colnames(contrast) = "group2"
    rownames(contrast)[1] = "(Intercept)"
    
    channels <- c(s1_ID, s2_ID)
    
    
    
    td = data_table[, channels]
    gg = data_table$Gene_name
    gg[which(is.na(gg))] = "NA"
    
    rownames(td) = gg

    
    
    m1 = rowMeans(td[, s1_ID], na.rm = T)
    m1[is.nan(m1)] = NA
    m2 = rowMeans(td[, s2_ID], na.rm = T)
    m2[is.nan(m2)] = NA
    lfc = m1-m2
    
    de.r<- eb.fit(td, design, contrast,lfc)
    
    de.r$ID <- rownames(de.r)
    
    ###### so as long as only 1 of the samples has data it can output some signficance, but if t test can not be conducted, 
    de.r = de.r%>%
      dplyr::select(ID, log2FC, everything())
    
    
    colnames(de.r)[1:2] = c("ID","log2fc")
    
    
    ##### add the number of samples with non-missing values in each of the two group
    
    
    
    nmn = rbindlist(lapply(1:nrow(de.r), function(x) {
      
      this_g = de.r$ID[x]
      which_r = which(data_table[,1] == this_g)
      
      d1 = data_table[which_r, s1_ID]
      d2 = data_table[which_r, s2_ID]
      
      n1 = sumna_not(d1)
      n2 = sumna_not(d2)
      
      
      df = data.frame(notNA1 = n1, 
                      notNA2 = n2, 
                      stringsAsFactors = F)
      return(df)
      
      
    }))
    
    result_de = data.frame(de.r, nmn, stringsAsFactors = F)
    
    write.table(result_de, table_output_file,
                quote = F, row.names = F, sep = "\t")
    
    
    return(result_de)
  }else{
    
    
    cat("samples of the class of interest are not present in the data")
  }
}


####possibility to do paired test, by including the tissue ID 





limma_de_table_paired= function(interest_sample1, 
                         interest_sample2,
                         interest_sample1_invert = F, 
                         interest_sample2_invert  = F, 
                         data_table_file, 
                         col_annot_file,
                         adjust_purity = T,
                         table_output_file)
  
{
  # interest_sample1 ="BMA"
  # interest_sample2 = "PB"
  # interest_sample1_invert = F
  # interest_sample2_invert  = F
  # data_table_file = "/Users/ginny/My Drive/CPTAC_AML/pb_bma_not_merge/abundance_gene_MD_geneName.tsv"
  # col_annot_file =  "/Users/ginny/My Drive/CPTAC_AML/pair_pb_bma_annot.tsv"
  # adjust_purity = F
  # table_output_file =  "/Users/ginny/My Drive/CPTAC_AML/DE/pair_pb_bma_prot.tsv"
  # # # 
   
  data_table = fread(data_table_file, stringsAsFactors = F, data.table = F)
  
  col_annot = fread(col_annot_file, stringsAsFactors = F, data.table = F)
  
  all_sample = col_annot$caseID
  p_sample = colnames(data_table[,-1])
  
  inter_sample = intersect(all_sample, p_sample)
  
  
  ### trim data to samples in intersection 
  
  col_annot = col_annot%>%
    dplyr::filter(caseID %in% inter_sample)
  
  wi = which(colnames(data_table) %in% inter_sample)
  
  data_table = data_table[, c(1,wi)]
  
  
  ### retrieve samples of interest
  
  s1 = col_annot%>%
    dplyr::filter(tumorClass %in% interest_sample1)
  
  if(interest_sample1_invert == T)
  {
    
    s1 = col_annot%>%
      dplyr::filter(!tumorClass %in% interest_sample1)
    
  }
  s1_ID = intersect(s1$caseID, colnames(data_table))
  
  s1p = data.frame(caseID = s1_ID, stringsAsFactors = F)%>%
    dplyr::left_join(col_annot, by = "caseID")
  s1_purity = s1p$purity
  s1_tissue = s1p$tissue
  
  
  s2 = col_annot%>%
    dplyr::filter(tumorClass %in% interest_sample2)
  
  if(interest_sample2_invert == T)
  {
    
    s2 = col_annot%>%
      dplyr::filter(!tumorClass %in% interest_sample2)
    
  }
  s2_ID = intersect(s2$caseID, colnames(data_table))
  
  
  
  s2p = data.frame(caseID = s2_ID, stringsAsFactors = F)%>%
    dplyr::left_join(col_annot, by = "caseID")
  s2_purity = s2p$purity
  s2_tissue = s2p$tissue
  
  
  
  
  
  if(length(s1_ID)>0 & length(s2_ID)>0)
  {
    
    group = as.factor(c(rep(1, length(s1_ID)),rep(2, length(s2_ID))))
    tissue = as.factor(c(s1_tissue, s2_tissue))
    
    design = model.matrix( ~ group + tissue)
    
    if(adjust_purity == T)
    {
      purity = c(s1_purity, s2_purity)
      
      design = model.matrix( ~  group + purity + tissue)
      
    }
    
    
    contrast = makeContrasts(0-group2, levels = colnames(design))
    
    
    
    colnames(contrast) = "group2"
    rownames(contrast)[1] = "(Intercept)"
    
    channels <- c(s1_ID, s2_ID)
    
    
    
    td = data_table[, channels]
    
    rownames(td) = data_table[,1]
    
    m1 = rowMeans(td[, s1_ID], na.rm = T)
    m1[is.nan(m1)] = NA
    m2 = rowMeans(td[, s2_ID], na.rm = T)
    m2[is.nan(m2)] = NA
    lfc = m1-m2
    
    de.r<- eb.fit(td, design, contrast,lfc)
    
    de.r$ID <- rownames(de.r)
    
    ###### so as long as only 1 of the samples has data it can output some signficance, but if t test can not be conducted, 
    de.r = de.r%>%
      dplyr::select(ID, log2FC, everything())
    
    
    colnames(de.r)[1:2] = c("ID","log2fc")
    
    
    ##### add the number of samples with non-missing values in each of the two group
    
    
    
    nmn = rbindlist(lapply(1:nrow(de.r), function(x) {
      
      this_g = de.r$ID[x]
      which_r = which(data_table[,1] == this_g)
      
      d1 = data_table[which_r, s1_ID]
      d2 = data_table[which_r, s2_ID]
      
      n1 = sumna_not(d1)
      n2 = sumna_not(d2)
      
      
      df = data.frame(notNA1 = n1, 
                      notNA2 = n2, 
                      stringsAsFactors = F)
      return(df)
      
      
    }))
    
    result_de = data.frame(de.r, nmn, stringsAsFactors = F)
    
    write.table(result_de, table_output_file,
                quote = F, row.names = F, sep = "\t")
    
    
    return(result_de)
  }else{
    
    
    cat("samples of the class of interest are not present in the data")
  }
}





#### the only variable is purity (used in cancer vs normal comparison) 

limma_de_table_purity= function(interest_sample1, 
                                interest_sample2,
                                interest_sample1_invert = F, 
                                interest_sample2_invert  = F, 
                                data_table_file, 
                                col_annot_file,
                                table_output_file)
  
{
  # interest_sample1 = c("AML")
  # interest_sample2 = c("Normal")
  # interest_sample1_invert = F
  # interest_sample2_invert  = F
  # data_table_file  = "G:/My Drive/nonCCRCC20211020/dataFreeze/proteome/protein_MD_allGene_kinases.tsv" ### data table
  # col_annot_file  = "G:/My Drive/nonCCRCC20211020/clinical/col_annot_simp_meta_20211115.tsv"
  # table_output_file = "G:/My Drive/nonCCRCC20211020/DE/phosphosite/kinase/purity_amlvNat.tsv"
  
  data_table = fread(data_table_file, stringsAsFactors = F, data.table = F)
  
  col_annot = fread(col_annot_file, stringsAsFactors = F, data.table = F)
  
  all_sample = col_annot$caseID
  p_sample = colnames(data_table[,-1])
  
  inter_sample = intersect(all_sample, p_sample)
  
  
  ### trim data to samples in intersection 
  
  col_annot = col_annot%>%
    dplyr::filter(caseID %in% inter_sample)
  
  wi = which(colnames(data_table) %in% inter_sample)
  
  data_table = data_table[, c(1,wi)]
  
  
  ### retrieve samples of interest
  
  s1 = col_annot%>%
    dplyr::filter(tumorClass %in% interest_sample1)
  
  if(interest_sample1_invert == T)
  {
    
    s1 = col_annot%>%
      dplyr::filter(!tumorClass %in% interest_sample1)
    
  }
  s1_ID = intersect(s1$caseID, colnames(data_table))
  
  s1p = data.frame(caseID = s1_ID, stringsAsFactors = F)%>%
    dplyr::left_join(col_annot, by = "caseID")
  s1_purity = s1p$purity
  
  
  
  s2 = col_annot%>%
    dplyr::filter(tumorClass %in% interest_sample2)
  
  if(interest_sample2_invert == T)
  {
    
    s2 = col_annot%>%
      dplyr::filter(!tumorClass %in% interest_sample2)
    
  }
  s2_ID = intersect(s2$caseID, colnames(data_table))
  
  
  
  s2p = data.frame(caseID = s2_ID, stringsAsFactors = F)%>%
    dplyr::left_join(col_annot, by = "caseID")
  s2_purity = s2p$purity
  
  purity = c(s1_purity, s2_purity)
  
  
  
  if(length(s1_ID)>0 & length(s2_ID)>0)
  {
    
    # group = as.factor(c(rep(1, length(s1_ID)),rep(2, length(s2_ID))))
    
    design = model.matrix( ~ purity)
    
    
    contrast = makeContrasts(purity, levels = colnames(design))
    
    
    
    colnames(contrast) = "purity"
    rownames(contrast)[1] = "(Intercept)"
    
    channels <- c(s1_ID, s2_ID)
    
    
    
    td = data_table[, channels]
    
    rownames(td) = data_table[,1]
    
    m1 = rowMeans(td[, s1_ID], na.rm = T)
    m1[is.nan(m1)] = NA
    m2 = rowMeans(td[, s2_ID], na.rm = T)
    m2[is.nan(m2)] = NA
    lfc = m1-m2
    
    de.r<- eb.fit(td, design, contrast,lfc)
    
    de.r$ID <- rownames(de.r)
    
    ###### so as long as only 1 of the samples has data it can output some signficance, but if t test can not be conducted, 
    de.r = de.r%>%
      dplyr::select(ID, log2FC, everything())
    
    
    colnames(de.r)[1:2] = c("ID","log2fc")
    
    
    ##### add the number of samples with non-missing values in each of the two group
    
    
    
    nmn = rbindlist(lapply(1:nrow(de.r), function(x) {
      
      this_g = de.r$ID[x]
      which_r = which(data_table[,1] == this_g)
      
      d1 = data_table[which_r, s1_ID]
      d2 = data_table[which_r, s2_ID]
      
      n1 = sumna_not(d1)
      n2 = sumna_not(d2)
      
      
      df = data.frame(notNA1 = n1, 
                      notNA2 = n2, 
                      stringsAsFactors = F)
      return(df)
      
      
    }))
    
    result_de = data.frame(de.r, nmn, stringsAsFactors = F)
    
    write.table(result_de, table_output_file,
                quote = F, row.names = F, sep = "\t")
    
    
    return(result_de)
  }else{
    
    
    cat("samples of the class of interest are not present in the data")
  }
}













volcano_plot_pmod = function(de_table_result, 
                             log2fc_cutoff, 
                             p_mod_cutoff, 
                             text_size, 
                             plot_output_file)
  
{
  #   de_table_result = onco12_limma
  #   log2fc_cutoff = 1
  #   p_mod_cutoff = 0.1
  #   plot_output_file = "G:/My Drive/nonCCRCC20211020/DE/protein/onco1VsOnco2_bc_limma.pdf"
  #   
  col_code  = data.frame(status = c("up","down","none"),
                         color = c("red","blue", "grey"),
                         stringsAsFactors = F)
  
  
  wil_gg = de_table_result%>%
    dplyr::select(ID, log2fc, p.mod)%>%
    na.omit()%>%
    dplyr::mutate(sig = -log10(p.mod))%>%
    dplyr::mutate(mark = case_when(log2fc > log2fc_cutoff & p.mod < p_mod_cutoff  ~ "up",
                                   log2fc < -log2fc_cutoff & p.mod < p_mod_cutoff ~ "down",
                                   T ~ "none"))%>%
    dplyr::mutate(sig = case_when(sig<10 ~ sig, 
                                  T ~ 10))
  
  
  mark_t = as.data.frame(table(wil_gg$mark))%>%
    dplyr::left_join(col_code, by = c("Var1" = "status"))%>%
    dplyr::arrange(Var1)
  
  
  
  ggplot(wil_gg, aes(log2fc, sig))+
    geom_point(aes(color = mark),
               size = 1,
               alpha = 0.6) +
    scale_colour_manual(values=mark_t$color, name = "", labels = mark_t$Var1)+
    geom_hline(yintercept = -log10(p_mod_cutoff), colour = "black", linetype = "dashed", size = 0.3,alpha = 0.5) + 
    geom_vline(xintercept = log2fc_cutoff, colour = "black", linetype="dashed", size = 0.3,alpha = 0.5) + 
    geom_vline(xintercept = -log2fc_cutoff, colour = "black", linetype="dashed", size = 0.3,alpha = 0.5) +
    xlim(min(wil_gg$log2fc), max(wil_gg$log2fc))+
    ylim(0,(max(wil_gg$sig)+0.1*(max(wil_gg$sig)- min(wil_gg$sig)))) +
    geom_text_repel(data = wil_gg[which(wil_gg$mark!="none"),],
                    aes(label = ID,color = mark),
                    size = text_size,
                    segment.size = 0.3,
                    segment.alpha = 0.3,
                    force = 4,
                    max.overlaps = Inf,
                    box.padding = unit(0.35, "lines"),
                    point.padding = unit(0.3, "lines"),
                    seed = 123)+
    theme_bw()+
    theme(text = element_text(size=15)) +
    ggtitle(paste0(" log2fc cutoff: ", log2fc_cutoff, " pvalue cutoff: ", p_mod_cutoff))+
    labs(x="log2 fold change", y = "-log10(pvalue)")
  
  
  ggsave(plot_output_file,
         width = 8, 
         height = 7)
  
}












volcano_plot_qmod = function(de_table_result, 
                             log2fc_cutoff, 
                             q_mod_cutoff, 
                             text_size,
                             ul,# yscale_uppper_limit 
                             plot_output_file)
  
{
  #   de_table_result = onco12_limma
  #   log2fc_cutoff = 1
  #   p_mod_cutoff = 0.1
  #   plot_output_file = "G:/My Drive/nonCCRCC20211020/DE/protein/onco1VsOnco2_bc_limma.pdf"
  #   
  col_code  = data.frame(status = c("up","down","none"),
                         color = c("red","blue", "grey"),
                         stringsAsFactors = F)
  
  
  wil_gg = de_table_result%>%
    dplyr::select(ID, log2fc, q.mod)%>%
    na.omit()%>%
    dplyr::mutate(sig = -log10(q.mod))%>%
    dplyr::mutate(mark = case_when(log2fc > log2fc_cutoff & q.mod < q_mod_cutoff  ~ "up",
                                   log2fc < -log2fc_cutoff & q.mod < q_mod_cutoff ~ "down",
                                   T ~ "none"))%>%
    dplyr::mutate(sig = case_when(sig<ul ~ sig, 
                                  T ~ ul))
  
  
  mark_t = as.data.frame(table(wil_gg$mark))%>%
    dplyr::left_join(col_code, by = c("Var1" = "status"))%>%
    dplyr::arrange(Var1)
  
  
  
  ggplot(wil_gg, aes(log2fc, sig))+
    geom_point(aes(color = mark),
               size = 1,
               alpha = 0.6) +
    scale_colour_manual(values=mark_t$color, name = "", labels = mark_t$Var1)+
    geom_hline(yintercept = -log10(q_mod_cutoff), colour = "black", linetype = "dashed", size = 0.3,alpha = 0.5) + 
    geom_vline(xintercept = log2fc_cutoff, colour = "black", linetype="dashed", size = 0.3,alpha = 0.5) + 
    geom_vline(xintercept = -log2fc_cutoff, colour = "black", linetype="dashed", size = 0.3,alpha = 0.5) +
    xlim(min(wil_gg$log2fc), max(wil_gg$log2fc))+
    ylim(0,(max(wil_gg$sig)+0.1*(max(wil_gg$sig)- min(wil_gg$sig)))) +
    geom_text_repel(data = wil_gg[which(wil_gg$mark!="none"),],
                    aes(label = ID,color = mark),
                    size = text_size,
                    segment.size = 0.3,
                    segment.alpha = 0.3,
                    force = 2,
                    max.overlaps = Inf,
                    box.padding = unit(0.35, "lines"),
                    point.padding = unit(0.3, "lines"),
                    seed = 123)+
    theme_bw()+
    theme(text = element_text(size=15)) +
    ggtitle(paste0(" log2fc cutoff: ", log2fc_cutoff, " qvalue cutoff: ", q_mod_cutoff))+
    labs(x="log2 fold change", y = "-log10(qvalue)")
  
  ggsave(plot_output_file,
         width = 8, 
         height = 7)
  
}




volcano_plot_qmod_positive = function(de_table_result, 
                                      log2fc_cutoff, 
                                      q_mod_cutoff, 
                                      plot_output_file)
  
{
  #   
  col_code  = data.frame(status = c("up","none"),
                         color = c("red", "grey"),
                         stringsAsFactors = F)
  
  
  wil_gg = de_table_result%>%
    dplyr::select(ID, log2fc, q.mod)%>%
    na.omit()%>%
    dplyr::mutate(sig = -log10(q.mod))%>%
    dplyr::mutate(mark = case_when(log2fc > log2fc_cutoff & q.mod < q_mod_cutoff  ~ "up",
                                   #log2fc < -log2fc_cutoff & q.mod < q_mod_cutoff ~ "down",
                                   T ~ "none"))%>%
    dplyr::mutate(sig = case_when(sig<10 ~ sig, 
                                  T ~ 10))
  
  
  mark_t = as.data.frame(table(wil_gg$mark))%>%
    dplyr::left_join(col_code, by = c("Var1" = "status"))%>%
    dplyr::arrange(Var1)
  
  
  
  ggplot(wil_gg, aes(log2fc, sig))+
    geom_point(aes(color = mark),
               size = 1,
               alpha = 0.6) +
    scale_colour_manual(values=mark_t$color, name = "", labels = mark_t$Var1)+
    geom_hline(yintercept = -log10(q_mod_cutoff), colour = "black", linetype = "dashed", size = 0.3,alpha = 0.5) + 
    geom_vline(xintercept = log2fc_cutoff, colour = "black", linetype="dashed", size = 0.3,alpha = 0.5) + 
    geom_vline(xintercept = -log2fc_cutoff, colour = "black", linetype="dashed", size = 0.3,alpha = 0.5) +
    xlim(min(wil_gg$log2fc), max(wil_gg$log2fc))+
    ylim(0,(max(wil_gg$sig)+0.1*(max(wil_gg$sig)- min(wil_gg$sig)))) +
    geom_text_repel(data = wil_gg[which(wil_gg$mark!="none"),],
                    aes(label = ID,color = mark),
                    size = 2,
                    segment.size = 0.3,
                    segment.alpha = 0.3,
                    force = 2,
                    max.overlaps = Inf,
                    box.padding = unit(0.35, "lines"),
                    point.padding = unit(0.3, "lines"),
                    seed = 123)+
    theme_bw()+
    theme(text = element_text(size=15)) +
    ggtitle(paste0(" log2fc cutoff: ", log2fc_cutoff, " qvalue cutoff: ", q_mod_cutoff))+
    labs(x="log2 fold change", y = "-log10(qvalue)")
  
  ggsave(plot_output_file,
         width = 7, 
         height = 7)
  
}




### query names 

volcano_plot_qmod_query= function(de_table_result, 
                                  log2fc_cutoff, 
                                  q_mod_cutoff, 
                                  qnames,
                                  plot_output_file)
  
{
  #   de_table_result = onco12_limma
  #   log2fc_cutoff = 1
  #   p_mod_cutoff = 0.1
  #   plot_output_file = "G:/My Drive/nonCCRCC20211020/DE/protein/onco1VsOnco2_bc_limma.pdf"
  #   
  col_code  = data.frame(status = c("up","none"),
                         color = c("red", "grey"),
                         stringsAsFactors = F)
  
  
  wil_gg = de_table_result%>%
    dplyr::select(ID, log2fc, q.mod)%>%
    na.omit()%>%
    dplyr::mutate(sig = -log10(q.mod))%>%
    dplyr::mutate(mark = case_when(log2fc > log2fc_cutoff & q.mod < q_mod_cutoff  ~ "up",
                                   #log2fc < -log2fc_cutoff & q.mod < q_mod_cutoff ~ "down",
                                   T ~ "none"))%>%
    dplyr::mutate(sig = case_when(sig<10 ~ sig, 
                                  T ~ 10))
  
  
  mark_t = as.data.frame(table(wil_gg$mark))%>%
    dplyr::left_join(col_code, by = c("Var1" = "status"))%>%
    dplyr::arrange(Var1)
  
  
  qwil = wil_gg%>%
    dplyr::filter(ID %in% qnames)
  
  ggplot(wil_gg, aes(log2fc, sig))+
    geom_point(aes(color = mark),
               size = 1,
               alpha = 0.6) +
    scale_colour_manual(values=mark_t$color, name = "", labels = mark_t$Var1)+
    geom_hline(yintercept = -log10(q_mod_cutoff), colour = "black", linetype = "dashed", size = 0.3,alpha = 0.5) + 
    geom_vline(xintercept = log2fc_cutoff, colour = "black", linetype="dashed", size = 0.3,alpha = 0.5) + 
    geom_vline(xintercept = -log2fc_cutoff, colour = "black", linetype="dashed", size = 0.3,alpha = 0.5) +
    xlim(min(wil_gg$log2fc), max(wil_gg$log2fc))+
    ylim(0,(max(wil_gg$sig)+0.1*(max(wil_gg$sig)- min(wil_gg$sig)))) +
    geom_text_repel(data =qwil,
                    aes(label = ID,color = mark),
                    size = 2,
                    segment.size = 0.3,
                    segment.alpha = 0.3,
                    force = 2,
                    max.overlaps = Inf,
                    box.padding = unit(0.35, "lines"),
                    point.padding = unit(0.3, "lines"),
                    seed = 123)+
    theme_bw()+
    theme(text = element_text(size=15)) +
    ggtitle(paste0(" log2fc cutoff: ", log2fc_cutoff, " qvalue cutoff: ", q_mod_cutoff))+
    labs(x="log2 fold change", y = "-log10(qvalue)")
  
  ggsave(plot_output_file,
         width = 7, 
         height = 7)
  
}






volcano_plot_qmod_onlyPositve = function(de_table_result, 
                                         log2fc_cutoff, 
                                         q_mod_cutoff, 
                                         plot_output_file)
  
{
  #   de_table_result = onco12_limma
  #   log2fc_cutoff = 1
  #   p_mod_cutoff = 0.1
  #   plot_output_file = "G:/My Drive/nonCCRCC20211020/DE/protein/onco1VsOnco2_bc_limma.pdf"
  #   
  col_code  = data.frame(status = c("up","down","none"),
                         color = c("red","blue", "grey"),
                         stringsAsFactors = F)
  
  
  wil_gg = de_table_result%>%
    dplyr::filter(log2fc>0)%>%
    dplyr::select(ID, log2fc, q.mod)%>%
    na.omit()%>%
    dplyr::mutate(sig = -log10(q.mod))%>%
    dplyr::mutate(mark = case_when(log2fc > log2fc_cutoff & q.mod < q_mod_cutoff  ~ "up",
                                   log2fc < -log2fc_cutoff & q.mod < q_mod_cutoff ~ "down",
                                   T ~ "none"))%>%
    dplyr::mutate(sig = case_when(sig<10 ~ sig, 
                                  T ~ 10))
  
  
  mark_t = as.data.frame(table(wil_gg$mark))%>%
    dplyr::left_join(col_code, by = c("Var1" = "status"))%>%
    dplyr::arrange(Var1)
  
  
  
  ggplot(wil_gg, aes(log2fc, sig))+
    geom_point(aes(color = mark),
               size = 1,
               alpha = 0.9) +
    scale_colour_manual(values=mark_t$color, name = "", labels = mark_t$Var1)+
    geom_hline(yintercept = -log10(q_mod_cutoff), colour = "black", linetype = "dashed", size = 0.3,alpha = 0.5) + 
    geom_vline(xintercept = log2fc_cutoff, colour = "black", linetype="dashed", size = 0.3,alpha = 0.5) + 
    #geom_vline(xintercept = -log2fc_cutoff, colour = "black", linetype="dashed", size = 0.3,alpha = 0.5) +
    xlim(0, max(wil_gg$log2fc))+
    ylim(0,(max(wil_gg$sig)+0.1*(max(wil_gg$sig)- min(wil_gg$sig)))) +
    geom_text_repel(data = wil_gg[which(wil_gg$mark!="none"),],
                    aes(label = ID,color = mark),
                    size = 3,
                    segment.size = 0.3,
                    segment.alpha = 0.3,
                    force = 2,
                    max.overlaps = Inf,
                    box.padding = unit(0.35, "lines"),
                    point.padding = unit(0.3, "lines"),
                    seed = 123)+
    theme_bw()+
    theme(text = element_text(size=15)) +
    ggtitle(paste0(" log2fc cutoff: ", log2fc_cutoff, " qvalue cutoff: ", q_mod_cutoff))+
    labs(x="log2 fold change", y = "-log10(qvalue)")
  
  ggsave(plot_output_file,
         width = 8, 
         height = 7)
  
}








volcano_plot_pmod_knownSites = function(de_table_result, 
                                        known_sites,
                                        log2fc_cutoff, 
                                        p_mod_cutoff, 
                                        plot_output_file)
  
{
   
  col_code  = data.frame(status = c("up","down","none","know_up","know_down", "know"),
                         color = c("grey60","grey60", "grey70","darkred", "darkblue", "grey44"),
                         stringsAsFactors = F)
  
  
  wil_gg = de_table_result%>%
    dplyr::select(ID, log2fc, p.mod)%>%
    na.omit()%>%
    dplyr::mutate(sig = -log10(p.mod))%>%
    dplyr::mutate(mark = case_when(p.mod <  p_mod_cutoff & log2fc >log2fc_cutoff ~ "up",
                                   p.mod < p_mod_cutoff & log2fc < -log2fc_cutoff ~ "down",
                                   T ~ "none"))%>%
    dplyr::mutate(mark = case_when(ID %in% known_sites & p.mod <  1 & log2fc >0 ~ "know_up",
                                   ID%in% known_sites & p.mod < 1 & log2fc < 0 ~ "know_down",
                                   ID %in% known_sites ~ "know",
                                   T ~ mark))%>%
    dplyr::mutate(sig = case_when(sig<10 ~ sig, 
                                  T ~ 10))
  
  
  t = which(wil_gg$mark == "none" & abs(wil_gg$log2fc)<1)
  
  lt = min(10000, length(t))
  sel_t = sample(t, lt)
  
  wil_gg_m = wil_gg%>%
    dplyr::filter(mark != "none" | abs(log2fc)>1)
  
  wil_gg = rbind(wil_gg_m, wil_gg[sel_t,])
  
  
  wil_gg<- wil_gg[order(as.numeric(factor(wil_gg$mark, levels = c("none","up","down","know","know_up","know_down")))),]
  
  
  mark_t = as.data.frame(table(wil_gg$mark))%>%
    dplyr::left_join(col_code, by = c("Var1" = "status"))%>%
    dplyr::arrange(Var1)
  
  text_df = wil_gg%>%
    dplyr::filter( mark == "know_up"| mark == "know_down"| mark == "know")
  
  ggplot(wil_gg, aes(log2fc, sig))+
    geom_point(aes(color = mark),
               size = 1,
               alpha = 1) +
    scale_colour_manual(values=mark_t$color, name = "", labels = mark_t$Var1)+
    geom_hline(yintercept = -log10(p_mod_cutoff), colour = "black", linetype = "dashed", size = 0.3,alpha = 0.5) + 
    geom_vline(xintercept = log2fc_cutoff, colour = "black", linetype="dashed", size = 0.3,alpha = 0.5) + 
    geom_vline(xintercept = -log2fc_cutoff, colour = "black", linetype="dashed", size = 0.3,alpha = 0.5) +
    xlim(min(wil_gg$log2fc), max(wil_gg$log2fc))+
    ylim(0,(max(wil_gg$sig)+0.01*(max(wil_gg$sig)- min(wil_gg$sig)))) +
    geom_text_repel(data = text_df,
                    aes(label = ID,color = mark),
                    size = 3,
                    segment.size = 0.3,
                    segment.alpha = 0.3,
                    force = 4,
                    max.overlaps = Inf,
                    box.padding = unit(0.35, "lines"),
                    point.padding = unit(0.3, "lines"),
                    seed = 123)+
    theme_bw()+
    theme(text = element_text(size=15)) +
    ggtitle(paste0(" log2fc cutoff: ", log2fc_cutoff, " pvalue cutoff: ", p_mod_cutoff))+
    labs(x="log2 fold change", y = "-log10(pvalue)")
  
  
  ggsave(plot_output_file,
         width = 8, 
         height = 7)
  
}





