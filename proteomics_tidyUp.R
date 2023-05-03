### this script is for proteomics data tidy up, applicable to single site TMT-I report (Gencode data)
### conversion from genecode IDs to geneSymbols can be done with either biomaRt or a prepared file (e.g. the one from Chenwei)

### when people first run FragPipe, in annotation file, please make sure there is no duplicated sample names across plexes, if there are multiple 
### channels for the same biological sample, please give it a suffix, for example, .2, .3 

globalProteome_tidyUp_biomaRt = function(TMTI_input_filename,
                                         replicate_suffix_pattern, 
                                         sel_sample, 
                                         sel_sample_final_name, 
                                output_filename)
{
  
  
#TMTI_input_filename = "/Users/ginny/My Drive/nesviLab_scripts/AML_original_data/proteome/ratio_gene_MD.tsv"
  
  
prot = fread(TMTI_input_filename,
               stringsAsFactors = F, data.table = F)

#replicate_suffix_pattern = c(".2")
if(length(replicate_suffix_pattern)>0)
{
  
sp = paste0("\\",replicate_suffix_pattern, "$")

rep_col = grep(sp, colnames(prot), value = T)
ori_col = gsub(sp, "", rep_col)

cols_to_remove = c()
av_rep = matrix(NA, nrow = nrow(prot), ncol = length(ori_col))

for(i in 1:length(ori_col))
{
  all_col = grep(ori_col[i], colnames(prot), value = T)
  cols_to_remove = c(cols_to_remove, all_col)
  this_mat = prot[, all_col]
  if(length(all_col) >1)
  {
    av_rep[,i] = rowMeans(this_mat, na.rm = T)
    
  }else{
    av_rep[,i] = this_mat
  }
  
}

av_rep[is.na(av_rep)] = NA 

colnames(av_rep) = ori_col

non_rep_col = setdiff(colnames(prot), cols_to_remove)
prot = cbind(prot[, non_rep_col], av_rep)
}

  

human = useMart("ensembl", dataset =  "hsapiens_gene_ensembl")

gid = gsub("\\..*", "", prot$Index)

gn_symbol = getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"),
                  filters = "ensembl_gene_id",
                  values = gid, 
                  mart = human)%>%
  dplyr::group_by(ensembl_gene_id)%>%
  dplyr::filter(row_number() == 1)%>%
  dplyr::mutate(hgnc_symbol = ifelse(nchar(hgnc_symbol)<1, ensembl_gene_id, hgnc_symbol))


prot_tidy = prot%>%
  dplyr::mutate(ensg = gid)%>%
  dplyr::left_join(gn_symbol, by = c("ensg" = "ensembl_gene_id"))%>%
  dplyr::select(hgnc_symbol, ensg,everything())%>%
  dplyr::mutate(hgnc_symbol = ifelse(is.na(hgnc_symbol), Index, hgnc_symbol))%>%
  dplyr::select(-ensg)%>%
  dplyr::select(-Index)%>%
  dplyr::select(-NumberPSM)%>%
  dplyr::select(-ProteinID)%>%
  dplyr::select(-MaxPepProb)%>%
  dplyr::select(-ReferenceIntensity)%>%
  dplyr::group_by(hgnc_symbol)%>%
  dplyr::summarise_if(is.numeric, mean, na.rm = T)%>%
  dplyr::ungroup()%>%
  as.data.frame()%>%
  dplyr::mutate_all(~ifelse(is.nan(.), NA, .))


if(length(sel_sample)>0)
{
  
  prot_tidy = prot_tidy[, c("hgnc_symbol", sel_sample)]
  
}

colnames(prot_tidy) = c("Gene_name", sel_sample_final_name)

write.table(prot_tidy, output_filename, 
            quote = F, row.names = F, sep = "\t")

}


### use an existing mapping file, e.g. the one from Chenwei 


globalProteome_tidyUp_mapFile = function(TMTI_input_filename,
                                         replicate_suffix_pattern, 
                                          map_filename, 
                                         sel_sample, 
                                         sel_sample_final_name, 
                                         output_filename)
{
  
  # TMTI_input_filename = "/Users/ginny/My Drive/nesviLab_scripts/AML_original_data/proteome/ratio_gene_MD.tsv"
  # map_filename ="/Users/ginny/My Drive/CPTAC_AML/GENCODE.V42.basic.CHR.protein.selection.mapping.csv"
  # sel_sample = annot$caseID
  # sel_sample_final_name = annot$Case_ID 
  # output_filename = "/Users/ginny/My Drive/nesviLab_scripts/test_output/ratio_gene_MD_tidy_mapFile_proteogenomics.tsv"
  # 
  
  prot = fread(TMTI_input_filename,
               stringsAsFactors = F, data.table = F)
  if(length(replicate_suffix_pattern)>0)
  {
  sp = paste0("\\",replicate_suffix_pattern, "$")
  
  rep_col = grep(sp, colnames(prot), value = T)
  ori_col = gsub(sp, "", rep_col)
  
  cols_to_remove = c()
  av_rep = matrix(NA, nrow = nrow(prot), ncol = length(ori_col))
  
  for(i in 1:length(ori_col))
  {
    all_col = grep(ori_col[i], colnames(prot), value = T)
    cols_to_remove = c(cols_to_remove, all_col)
    this_mat = prot[, all_col]
    if(length(all_col) >1)
    {
      av_rep[,i] = rowMeans(this_mat, na.rm = T)
      
    }else{
      av_rep[,i] = this_mat
    }
    
  }
  
  av_rep[is.na(av_rep)] = NA 
  
  colnames(av_rep) = ori_col
  
  non_rep_col = setdiff(colnames(prot), cols_to_remove)
  prot = cbind(prot[, non_rep_col], av_rep)
  
  }
  
  
  gmap = fread(map_filename,
               stringsAsFactors = F, data.table = F)%>%
    dplyr::select(Gene_id, Gene_name)%>%
    unique()
  
  prot_tidy = prot%>%
    dplyr::left_join(gmap, by = c("Index" = "Gene_id"))%>%
    dplyr::select(Gene_name, everything())%>%
    dplyr::select(-Index)%>%
    dplyr::select(-NumberPSM)%>%
    dplyr::select(-ProteinID)%>%
    dplyr::select(-MaxPepProb)%>%
    dplyr::select(-ReferenceIntensity)%>%
    dplyr::group_by(Gene_name)%>%
    dplyr::summarise_if(is.numeric, mean, na.rm = T)%>%
    ungroup()%>%
    dplyr::mutate_all(function(x) ifelse(is.nan(x), NA, x))%>%
    as.data.frame()
    
  
  if(length(sel_sample)>0)
  {
    
    prot_tidy = prot_tidy[, c("Gene_name", sel_sample)]
    
  }
  colnames(prot_tidy) = c("Gene_name", sel_sample_final_name)
    
  write.table(prot_tidy, output_filename, 
              quote = F, row.names = F, sep = "\t")
  
}





### the change is to add the peptide flanking sequence 


### with biomart 


phosphoSingleSite_tidyUp_biomaRt = function(TMTI_input_filename,
                                                replicate_suffix_pattern, 
                                                sel_sample, 
                                                sel_sample_final_name, 
                                                output_filename)
{
  
  # TMTI_input_filename = "/Users/ginny/My Drive/nesviLab_scripts/AML_original_data/phosphoproteome/ratio_single-site_MD.tsv"
  # map_filename ="/Users/ginny/My Drive/CPTAC_AML/GENCODE.V42.basic.CHR.protein.selection.mapping.csv"
  # sel_sample = annot$caseID
  # sel_sample_final_name = annot$Case_ID 
  # output_filename = "/Users/ginny/My Drive/nesviLab_scripts/test_output/ratio_signle-site_MD_tidy_mapFile_proteogenomics.tsv"
  # 
  
  phospho = fread(TMTI_input_filename,
                  stringsAsFactors = F, data.table = F)
  if(length(replicate_suffix_pattern)>0)
  {
  sp = paste0("\\",replicate_suffix_pattern, "$")
  
  rep_col = grep(sp, colnames(phospho), value = T)
  ori_col = gsub(sp, "", rep_col)
  
  cols_to_remove = c()
  av_rep = matrix(NA, nrow = nrow(phospho), ncol = length(ori_col))
  
  for(i in 1:length(ori_col))
  {
    all_col = grep(ori_col[i], colnames(phospho), value = T)
    cols_to_remove = c(cols_to_remove, all_col)
    this_mat = phospho[, all_col]
    if(length(all_col) >1)
    {
      av_rep[,i] = rowMeans(this_mat, na.rm = T)
      
    }else{
      av_rep[,i] = this_mat
    }
    
  }
  
  av_rep[is.na(av_rep)] = NA 
  
  colnames(av_rep) = ori_col
  
  non_rep_col = setdiff(colnames(phospho), cols_to_remove)
  phospho = cbind(phospho[, non_rep_col], av_rep)
  
  }
  
  human = useMart("ensembl", dataset =  "hsapiens_gene_ensembl")
  
  gid = gsub("\\..*", "", phospho$Gene)
  
  gn_symbol = getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"),
                    filters = "ensembl_gene_id",
                    values = gid, 
                    mart = human)%>%
    dplyr::group_by(ensembl_gene_id)%>%
    dplyr::filter(row_number() == 1)%>%
    dplyr::mutate(hgnc_symbol = ifelse(nchar(hgnc_symbol)<1, ensembl_gene_id, hgnc_symbol))
  
  
  
  phospho_tidy = phospho%>%
    dplyr::mutate(ensg = gid)%>%
    dplyr::left_join(gn_symbol, by = c("ensg" = "ensembl_gene_id"))%>%
    dplyr::select(hgnc_symbol, ensg,everything())%>%
    dplyr::mutate(hgnc_symbol = ifelse(is.na(hgnc_symbol), Gene, hgnc_symbol))%>%
    dplyr::mutate(site = gsub(".*_","",Index))%>%
    dplyr::mutate(geneSite = paste(hgnc_symbol, site, sep = "_"))%>%
    dplyr::select(-Index)%>%
    dplyr::select(-Gene)%>%
    dplyr::select(-ProteinID)%>%
    dplyr::select(-Peptide)%>%
    dplyr::select(-MaxPepProb)%>%
    dplyr::select(-ReferenceIntensity)%>%
    dplyr::select(-site)%>%
    dplyr::select(geneSite, SequenceWindow, everything())%>%
    dplyr::group_by(geneSite, SequenceWindow)%>%
    dplyr::summarise_if(is.numeric, mean, na.rm = T)%>%
    ungroup()%>%
    dplyr::mutate_all(function(x) ifelse(is.nan(x), NA, x))%>%
    as.data.frame()         ##### somewhere downstream I only take unique row IDs so if i take geneSite as ID it might be problematic
  
  
  if(length(sel_sample)>0)
  {
    phospho_tidy = phospho_tidy[, c("geneSite","SequenceWindow",  sel_sample)]
  }
  
  
  colnames(phospho_tidy) = c("Gene_site","SequenceWindow", sel_sample_final_name)
  
  write.table(phospho_tidy, output_filename, 
              quote = F, row.names = F, sep = "\t")
  
}













phosphoSingleSite_tidyUp_mapFile = function(TMTI_input_filename,
                                                replicate_suffix_pattern, 
                                         map_filename, 
                                         sel_sample, 
                                         sel_sample_final_name, 
                                         output_filename)
{
  
  # TMTI_input_filename = "/Users/ginny/My Drive/nesviLab_scripts/AML_original_data/phosphoproteome/ratio_single-site_MD.tsv"
  # map_filename ="/Users/ginny/My Drive/CPTAC_AML/GENCODE.V42.basic.CHR.protein.selection.mapping.csv"
  # sel_sample = annot$caseID
  # sel_sample_final_name = annot$Case_ID 
  # output_filename = "/Users/ginny/My Drive/nesviLab_scripts/test_output/ratio_signle-site_MD_tidy_mapFile_proteogenomics.tsv"
  # 
  
  phospho = fread(TMTI_input_filename,
               stringsAsFactors = F, data.table = F)
  
  if(length(replicate_suffix_pattern)>0)
  {
  sp = paste0("\\",replicate_suffix_pattern, "$")
  
  rep_col = grep(sp, colnames(phospho), value = T)
  ori_col = gsub(sp, "", rep_col)
  
  cols_to_remove = c()
  av_rep = matrix(NA, nrow = nrow(phospho), ncol = length(ori_col))
  
  for(i in 1:length(ori_col))
  {
    all_col = grep(ori_col[i], colnames(phospho), value = T)
    cols_to_remove = c(cols_to_remove, all_col)
    this_mat = phospho[, all_col]
    if(length(all_col) >1)
    {
      av_rep[,i] = rowMeans(this_mat, na.rm = T)
      
    }else{
      av_rep[,i] = this_mat
    }
    
  }
  
  av_rep[is.na(av_rep)] = NA 
  
  colnames(av_rep) = ori_col
  
  non_rep_col = setdiff(colnames(phospho), cols_to_remove)
  phospho = cbind(phospho[, non_rep_col], av_rep)
  }
  
  ######
  
  gmap = fread(map_filename,
               stringsAsFactors = F, data.table = F)%>%
    dplyr::select(Gene_id, Gene_name)%>%
    unique()
  
  
 phospho_tidy = phospho%>%
    dplyr::left_join(gmap, by = c("Gene" = "Gene_id"))%>%
    dplyr::mutate(site = gsub(".*_","",Index))%>%
    dplyr::mutate(geneSite = paste(Gene_name, site, sep = "_"))%>%
    dplyr::select(-Index)%>%
    dplyr::select(-Gene)%>%
    dplyr::select(-ProteinID)%>%
    dplyr::select(-Peptide)%>%
    dplyr::select(-MaxPepProb)%>%
    dplyr::select(-ReferenceIntensity)%>%
    dplyr::select(-site)%>%
    dplyr::select(geneSite, SequenceWindow, everything())%>%
    dplyr::group_by(geneSite, SequenceWindow)%>%
    dplyr::summarise_if(is.numeric, mean, na.rm = T)%>%
    ungroup()%>%
    dplyr::mutate_all(function(x) ifelse(is.nan(x), NA, x))%>%
    as.data.frame()         ##### somewhere downstream I only take unique row IDs so if i take geneSite as ID it might be problematic
  

  if(length(sel_sample)>0)
  {
    phospho_tidy = phospho_tidy[, c("geneSite","SequenceWindow",  sel_sample)]
  }
 
 
  colnames(phospho_tidy) = c("Gene_site","SequenceWindow", sel_sample_final_name)
  
  write.table(phospho_tidy, output_filename, 
              quote = F, row.names = F, sep = "\t")
  
}










