rm(list=ls())
options(stringsAsFactors = F)
library(stringr)
library(tidyverse)
library(data.table)
library(patchwork)
library(ggsci)
library(RColorBrewer)
library(Seurat)

# fp_sc: scRNA-seq data dir
# fp_output: output dir
fp_output='./result/Fig5-part2/'

sce_obj = readRDS(paste0(fp_sc,'CD colon scRNA.rds'))
col.group = c('#61C0BF','#E9D985','#E76F51')
names(col.group)=c("Heal","NonI","Infl")
table(sce_obj$disease,sce_obj$cell_type)

# add score --------------------------------------------------------------------
gene_list = fgsea::gmtPathways('./tools/h.all.v2023.2.Hs.symbols.gmt')

sce_obj = AddModuleScore(sce_obj,features = gene_list)
colnames(sce_obj@meta.data)[33:82] = names(gene_list)

names(gene_list)

# define subtype ---------------------------------------------------------------
pro_gene = c("NAAA","CCDC88B","NAGLU")
sup_gene = c("CWC15")
markers = intersect(c(pro_gene,sup_gene),rownames(sce_obj))
df = FetchData(sce_obj,vars = markers)
head(df)
df2 = apply(df,2,function(x) ifelse(x>0,'Pos','Neg'))
head(df2)
df2 = as.data.frame(df2)

sce_obj = AddMetaData(sce_obj,metadata = df2)

md = sce_obj@meta.data %>% 
  select(cell_type,all_of(markers),all_of(names(gene_list)))

p.df = read.csv(paste0(fp_output,"DOR summarised result.csv"),row.names = 1)

p.df %>% 
  group_by(gene) %>% 
  top_n(n=1,wt=dor) %>% 
  ungroup() -> cell_gene_pair

save(md,cell_gene_pair,file = paste0(fp_output,'score.rda'))

# plot  -------------------------------------------------------------------
load(paste0(fp_output,'score.rda'))
hallmark = colnames(md)[grepl('^HALLMARK',colnames(md))]

get_stat=function(data){
  # data=md_sub
  res_stat = data %>% 
    pivot_longer(cols = 2:ncol(data),names_to = 'pathway',values_to = 'score') %>% 
    group_by(pathway) %>% 
    rstatix::wilcox_test(score~group)
  res_mean = data %>% 
    pivot_longer(cols = 2:ncol(data),names_to = 'pathway',values_to = 'score') %>% 
    group_by(pathway,group) %>% 
    summarise(mean=mean(score)) %>% 
    pivot_wider(names_from = 'group',values_from = 'mean')
  
  res = res_stat %>% inner_join(res_mean) %>% 
    select(pathway,p,Neg,Pos)
  return(res)
}


res = lapply(1:nrow(cell_gene_pair),function(i){
  gs = cell_gene_pair$gene[i]
  cs = cell_gene_pair$cell_type[i]
  
  md_sub=md %>% filter(cell_type==cs) %>% 
    rename('group'=all_of(gs)) %>% 
    select(group,all_of(hallmark))
  
  res_tmp = get_stat(md_sub)
  return(res_tmp)
})

inf = hallmark[grepl('IL2|IL6|TNF|INTERFERON|INFLAMMATORY',hallmark)]

col.group = c('#D57A66','#A3CBB2')
names(col.group)=c('Pos','Neg')

p_list= lapply(1:nrow(cell_gene_pair),function(i){
  gs = cell_gene_pair$gene[i]
  cs = cell_gene_pair$cell_type[i]
  
  md_sub=md %>% filter(cell_type==cs) %>% 
    rename('group'=all_of(gs)) %>% 
    select(group,all_of(inf))
  
  
  p_tmp = md_sub %>% 
    pivot_longer(cols = 2:ncol(md_sub),names_to = 'pathway',values_to = 'score') %>% 
    mutate(pathway=str_replace(pathway,'HALLMARK\\_','')) %>% 
    ggplot(aes(x=group,y=score,color=group))+
    geom_violin(aes(fill=group),trim = FALSE)+
    geom_boxplot(width=0.2)+
    # theme_base()+
    theme_light(base_size = 12)+
    scale_color_manual(values = col.group)+
    scale_fill_manual(values = col.group)+
    facet_wrap(~pathway,scales = 'free_y')+
    ggpubr::stat_compare_means(aes(group=group),
                               label = 'p.format',label.x.npc = 0.15)+
    labs(x=paste0(cs,'(',gs,' status)'),y='Score')+guides(fill='none',color='none')+
    theme(strip.text.x = element_text(size = 9))
  return(p_tmp)

})

p_list[[1]]
p_list[[2]]
p_list[[3]]
p_list[[4]]
for (i in 1:nrow(cell_gene_pair)) {
  gs = cell_gene_pair$gene[i]
  cs = cell_gene_pair$cell_type[i]
  p = p_list[[i]] 
  ggsave(filename = paste0(fp_output,cs,'_',gs,'_inflammation score.pdf'),
         plot = p,
         width = 17,height = 17/3*2,units = 'cm')
}



