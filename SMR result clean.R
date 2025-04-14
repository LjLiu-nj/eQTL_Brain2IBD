rm(list=ls())
options(stringsAsFactors = F)
library(stringr)
library(tidyverse)
library(data.table)

fp_output = './result/smr/'
# ao = TwoSampleMR::available_outcomes()
# ao$sample_size = ifelse(
#   is.na(ao$sample_size),ao$ncase+ao$ncontrol,ao$sample_size
# )
# saveRDS(ao,file = './tools/IEU available outcome.rds')

ao = readRDS('./tools/IEU available outcome.rds')

ao %>% filter(consortium=='IIBDGC') %>% 
  filter(population=='European') -> a2
ao %>% filter(consortium=='UKB') %>% 
  filter(population=='European') %>% 
  filter(grepl('inflammatory bowel disease',trait,ignore.case=T))



# 数据说明 ##
# IIBDGC # 
# ieu-a-30: Crohn's disease
# ieu-a-31: Inflammatory bowel disease
# ieu-a-32: Ulcerative colitis

# UKB # 
# 20002_1461: UKB-IBD --
# 20002_1462: UKB-CD --
# 20002_1463: UKB-UC
# phecode-555.2: UKB-UC
# phecode-555.21: UKB-UC(chronic) --
# icd10-K50:UKB-CD(regional enteritis)
# icd10-K51:UKB-UC
# ULCERNAS:UKB-UC --

# FinnGen # 
# CD =c('CHRONSMALL','CHRONLARGE')
# UC='K11_UC_STRICT2'
# IBD = 'K11_KELAIBD'

# SMR results ------------------------------------------------------------------
list.files('I:/smr/',pattern = 'ibdout',full.names = T) %>% 
  map(list.files,pattern='smr$') -> smr_fn_res

names(smr_fn_res) = str_replace(
  list.files('I:/smr/',pattern = 'ibdout'),'\\_ibdout','')

smr_list = lapply(1:length(smr_fn_res),function(i){
  fp=paste0('I:/smr/',names(smr_fn_res)[i],'_ibdout/')
  
  
  sapply(smr_fn_res[[i]],function(x) paste0(fp,x)) %>% 
    map(fread,data.table=F) -> smr_list_tmp
  
  return(smr_list_tmp)
})
names(smr_list) = names(smr_fn_res)
smr_list = unlist(smr_list,recursive = F)

smr_list=lapply(smr_list,function(df){
  
  if(any(grepl('ENSG00',df$Gene))){
    trans_gene = clusterProfiler::bitr(
      df$Gene,fromType = 'ENSEMBL',
      toType = 'SYMBOL',OrgDb = org.Hs.eg.db::org.Hs.eg.db
      )
    df = df %>% inner_join(trans_gene,by=c('Gene'='ENSEMBL')) %>% 
      select(-Gene) %>% 
      rename(Gene=SYMBOL)
  }
  return(df)
  
})

smr_list = lapply(smr_list,function(x){
  # x = smr_list$`Brain_Amygdala.ieu-a-30.smr`
  x = x %>% mutate(fdr = p.adjust(p_SMR,method='fdr')) %>% 
    mutate(Significance = case_when(
      fdr < 0.05  ~ 'FDR',
      p_SMR < 0.05 ~ 'Nominal',
      TRUE ~ 'No'
    )) %>% 
    mutate(HEIDI_test = case_when(p_HEIDI > 0.5 ~ 'Pass',
                                  TRUE ~ 'Fail'))
  # table(x$Significance,x$HEIDI_test)
  return(x)
})

saveRDS(smr_list,file = paste0(fp_output,'Merged SMR result.rds'))





