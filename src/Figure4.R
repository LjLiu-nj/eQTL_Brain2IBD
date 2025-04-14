rm(list=ls())
options(stringsAsFactors = F)
library(stringr)
library(tidyverse)
library(data.table)
library(patchwork)
library(ggsci)


# fp_smr: SMR result dir
# fp_output: output dir
fp_output='./result/Fig4/'
# smr result list and FDR adjust -----------------------------------------------
smr_list0 = readRDS(paste0(fp_smr,'Merged SMR result.rds'))
smr_list0 = lapply(smr_list0,function(x){
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

# BrainMeta and sub region co-analysis ----------------------------------------------

# IBD;CD;UC cohorts selection # 
disease_list = list(
  IBD=c('ieu-a-31','categorical-20002-both_sexes-1461','finngen_R8_K11_KELAIBD'),
  CD=c('ieu-a-30','icd10-K50-both_sexes','finngen_R8_CHRONLARGE'),
  UC=c('ieu-a-32','icd10-K51-both_sexes','finngen_R8_K11_UC_STRICT2')
)

# Brain meta and sub region cohorts # 
tissue_list = c('BrainMeta','Colon_Sigmoid','Colon_Transverse')

# select the combination # 
dt = expand.grid(unlist(disease_list),tissue_list)

dt_anno = plyr::ldply(disease_list,data.frame)
colnames(dt_anno)=c('Var3','Var1')
dt = dt %>% inner_join(dt_anno) %>% 
  mutate(traits = paste(dt$Var2,dt$Var1,'smr',sep = '.'))

dt = dt %>% 
  mutate(cohorts = ifelse(grepl('ieu',Var1),'IIBDGC',
                          ifelse(grepl('finngen',Var1),'FinnGen','UKB')))

smr_list = smr_list0[paste(dt$Var2,dt$Var1,'smr',sep = '.')]

# Compare the brain and colon weight among the IBD subtype ---------------------

# * volcano plot ---------------------------------------------------------------
volcano_plot = function(df,bt=0,labels=NULL){
  # df = smr_list$`BrainMeta.ieu-a-31.smr`
  df$labels = ifelse(df$Gene %in% labels, df$Gene ,NA)
  
  df$change = ifelse(df$p_SMR>=0.05,'Not',
                     ifelse(df$b_SMR>bt,'Up','Down'))
  
  
  df$dot.alpha = ifelse(df$change=='Not',0.3,1)
  
  color.value = c("#546de5", "#d2dae2","#ff4757")
  names(color.value) = c('Down','Not','Up')  
  
  p=ggplot(df,aes(x=b_SMR,y=-log10(p_SMR),color=change))+
    geom_point(aes(alpha=dot.alpha),size=1.6)+
    scale_color_manual(values =color.value)+
    theme(text = element_text(size =3))+
    theme_bw()+
    labs(x='b_SMR',y='-log10(p_SMR)')+
    # geom_vline(xintercept=c(-logFC_cutoff,logFC_cutoff),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8)+ 
    guides(alpha='none') +  # 去除透明度的图例
    ggrepel::geom_text_repel(data = df, aes(x = b_SMR,
                                            y =-log10(p_SMR),
                                            label = labels),
                             size = 3,box.padding = unit(0.5, "lines"),
                             point.padding = unit(0.8, "lines"), 
                             colour='black',
                             max.overlaps = 50,
                             segment.color = "black",
                             show.legend = FALSE)
  return(p)
}

volcano_plot(df = smr_list$`BrainMeta.categorical-20002-both_sexes-1461.smr`)

dt2 = dt

vp_list = lapply(1:nrow(dt2),function(x){
  df=smr_list[[dt2[x,'traits']]]
  p_tmp = volcano_plot(df)+
    labs(x='Beta',y='-log10(p value)',color=NULL,
         title = paste0(dt2[x,'Var1'],' eQTL on ',dt2[x,'Var3']))+
    theme(plot.title = element_text(hjust = 0.5))
  return(p_tmp)  
})

wrap_plots(vp_list)+plot_layout(guides = 'collect')
ggsave(filename = paste0(fp_output,'Volcano plot list.pdf'),
       width = 24,height = 20)


# * significant gene number ----------------------------------------------------
get_sig_terms = function(df,adjust=T,HEIDI=T){
  if (adjust){
    df = df %>% filter(fdr<0.05)
  } else(
    df = df %>% filter(p_SMR<0.05)
  )
  if (HEIDI){
    df = df %>% filter(p_HEIDI>0.05)
  } 
  x = df %>% 
    mutate(condition=ifelse(b_SMR>0,'Pro','Sup'))
  x = split(x$Gene,x$condition)
  return(x)
}
get_sig_terms(smr_list$`BrainMeta.ieu-a-30.smr`,adjust = F)
sig_terms = lapply(smr_list,get_sig_terms,adjust=F)

dt2 = dt

sig_terms = sig_terms[dt2$traits]

# venn plot # 
library(ggvenn)
for (i in c('ieu-a-31','ieu-a-30','ieu-a-32')) {
  st_sub = sig_terms[grepl(i,names(sig_terms))]
  
  pro = lapply(st_sub,function(x) x[['Pro']])
  names(pro) = str_split(names(pro),pattern = '\\.',simplify = T)[,1]
  sup = lapply(st_sub,function(x) x[['Sup']])
  names(sup) = str_split(names(sup),pattern = '\\.',simplify = T)[,1]
  
  ggvenn(pro,
         show_percentage = F,
         stroke_color = "white",
         fill_color = c("#ffb2b2","#b2e7cb","#b2d4ec"),
         set_name_color = c("#ff0000","#4a9b83","#1d6295")) -> p1
  ggsave(filename = paste0(fp_output,i,'_','Pro Venn.pdf'),
         width = 5,height = 5,plot=p1)
  
  ggvenn(sup,
         show_percentage = F,
         stroke_color = "white",
         fill_color = c("#ffb2b2","#b2e7cb","#b2d4ec"),
         set_name_color = c("#ff0000","#4a9b83","#1d6295")) -> p2
  ggsave(filename = paste0(fp_output,i,'_','Sup Venn.pdf'),
         width = 5,height = 5,plot=p2)
}

st_shared_list = lapply(c('ieu-a-31','ieu-a-30','ieu-a-32'),function(i){
  st_sub = sig_terms[grepl(i,names(sig_terms))]
  pro = lapply(st_sub,function(x) x[['Pro']])
  names(pro) = str_split(names(pro),pattern = '\\.',simplify = T)[,1]
  sup = lapply(st_sub,function(x) x[['Sup']])
  names(sup) = str_split(names(sup),pattern = '\\.',simplify = T)[,1]
  
  res = list(
    Pro=Reduce(intersect,pro),
    Sup=Reduce(intersect,sup)
  )
  return(res)
})
names(st_shared_list) = c('ieu-a-31','ieu-a-30','ieu-a-32') 
saveRDS(st_shared_list,file = paste0(fp_output,'IIBDGC all shared genes.rds'))

pro = lapply(st_shared_list,function(x) x[['Pro']])
sup = lapply(st_shared_list,function(x) x[['Sup']])

ggvenn(pro,
       show_percentage = F,
       stroke_color = "white",
       fill_color = c("#ffb2b2","#b2e7cb","#b2d4ec"),
       set_name_color = c("#ff0000","#4a9b83","#1d6295")) -> p1
ggsave(filename = paste0(fp_output,'IIBDGC','_','Pro Shared Venn.pdf'),
       width = 5,height = 5,plot=p1)

ggvenn(sup,
       show_percentage = F,
       stroke_color = "white",
       fill_color = c("#ffb2b2","#b2e7cb","#b2d4ec"),
       set_name_color = c("#ff0000","#4a9b83","#1d6295")) -> p2
ggsave(filename = paste0(fp_output,'IIBDGC','_','Sup Shared Venn.pdf'),
       width = 5,height = 5,plot=p2)

# Gene Enrichment --------------------------------------------------------------
shorten_names <- function(x, n_word=4, n_char=200){
  if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > n_char))
  {
    # if (nchar(x) > 40) x <- substr(x, 1, n_char)
    x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],
                     collapse=" "), "...", sep="")
    return(x)
  } 
  else
  {
    return(x)
  }
}
ontology.barplot = function(go_enrich_df,categery_n = 10,n_word=4,title=NULL){
  CPCOLS <- c("#8DA1CB", "#FD8D62", "#66C3A5") #color 
  go_enrich_df <- go_enrich_df[order(go_enrich_df$pvalue),]
  go_enrich_df = go_enrich_df[go_enrich_df$pvalue<0.05,]
  go_enrich_df <- rbind(head(subset(go_enrich_df, ONTOLOGY == "BP"), categery_n),
                        head(subset(go_enrich_df, ONTOLOGY == "CC"), categery_n),
                        head(subset(go_enrich_df, ONTOLOGY == "MF"), categery_n))
  go_enrich_df=na.omit(go_enrich_df)
  
  go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))
  labels=sapply(go_enrich_df$Description,function(x) shorten_names(x,n_word))
  names(labels) = rev(1:nrow(go_enrich_df))
  
  p <- ggplot(data=go_enrich_df, aes(x=number, y=-log10(pvalue), fill=ONTOLOGY)) +
    geom_bar(stat="identity", width=0.8) + coord_flip() + 
    scale_fill_manual(values = CPCOLS) + theme_bw() + 
    scale_x_discrete(labels=labels) +
    labs(x='GO terms',fill='Gene \nOntology') + 
    theme(text=element_text(size=16)) +
    labs(title = title)+
    facet_grid(ONTOLOGY~., scales = "free",space = "free")
  return(p)
}
ontology.dotplot = function(go_enrich_df,categery_n = 10,n_word=4,title=NULL){
  go_enrich_df <- go_enrich_df[order(go_enrich_df$pvalue),]
  go_enrich_df = go_enrich_df[go_enrich_df$pvalue<0.05,]
  go_enrich_df=na.omit(go_enrich_df)
  go_enrich_df = go_enrich_df[1:categery_n,]
  go_enrich_df=na.omit(go_enrich_df)
  go_enrich_df = go_enrich_df[order(go_enrich_df$Count,decreasing = T),]
  go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))
  labels=sapply(go_enrich_df$Description,function(x) shorten_names(x,n_word))
  names(labels) = rev(1:nrow(go_enrich_df))
  
  p=ggplot(data=go_enrich_df, aes(x=number,y=Count, color=pvalue)) +
    geom_point(aes(size=Count)) + coord_flip() + 
    scale_color_continuous(low="red", high="blue",trans='reverse')+
    theme_bw() + 
    scale_x_discrete(labels=labels) +
    labs(x='KEGG terms') + 
    theme(text=element_text(size=16))+
    labs(title = title)
  # p
  return(p)
}

# Seperate enrichment ---------------------------------------------------------- 
for (i in c('pro','sup')) {
  gl=get(i)
  gene_df = plyr::ldply(gl,data.frame)
  colnames(gene_df)=c('cluster','SYMBOL')
  write.csv(gene_df,paste0(fp_output,i,'_IIBDGC significant genes.csv'))
  
  trans_gene = clusterProfiler::bitr(
    unique(gene_df$SYMBOL),fromType = 'SYMBOL',
    toType = 'ENTREZID',OrgDb = org.Hs.eg.db::org.Hs.eg.db)
  gene_df = gene_df %>% inner_join(trans_gene)
  
  gene_list=split(gene_df$ENTREZID,gene_df$cluster)
  gene_list = gene_list[c('ieu-a-31','ieu-a-30','ieu-a-32')]
  
  ego.res = lapply(gene_list, function(gene){
    res_GO=clusterProfiler::enrichGO(
      gene=gene,OrgDb = org.Hs.eg.db::org.Hs.eg.db,
      keyType = "ENTREZID",
      ont='ALL', pAdjustMethod = "BH",
      pvalueCutoff=0.99,
      qvalueCutoff=0.99,
      readable =T)
    return(res_GO)
  })
  
  ekegg.res = lapply(gene_list, function(gene){
    es = clusterProfiler::enrichKEGG(
      gene=gene,organism = 'hsa',
      # universe     = gene_all,
      # use_internal_data =T,
      pvalueCutoff = 0.99,
      qvalueCutoff =0.99)
    es=DOSE::setReadable(es, OrgDb=org.Hs.eg.db::org.Hs.eg.db,keyType='ENTREZID')
    return(es)
  })
  
  p_go = lapply(1:length(ego.res),function(x){
    p_tmp = ontology.barplot(ego.res[[x]]@result,n_word = 10)+
      labs(title = names(ego.res)[x])+
      theme(plot.title = element_text(hjust = 0.5))
    return(p_tmp)
  })
  wrap_plots(p_go)+plot_layout(guides = 'collect')
  ggsave(filename = paste0(fp_output,'IIBDGC','_',i,' Shared Go.pdf'),
         width = 24,height = 6)
  
  p_kegg = lapply(1:length(ekegg.res),function(x){
    p_tmp = ontology.dotplot(ekegg.res[[x]]@result,n_word = 10)+
      labs(title = names(ekegg.res)[x])+
      theme(plot.title = element_text(hjust = 0.5))
    return(p_tmp)
  })
  wrap_plots(p_kegg)
  ggsave(filename = paste0(fp_output,'IIBDGC','_',i,' Shared KEGG.pdf'),
         width = 18,height = 6)
  
}

# Merge enrichment -------------------------------------------------------------
for (i in c('pro','sup')) {
  gl=get(i)
  gene_df = plyr::ldply(gl,data.frame)
  colnames(gene_df)=c('cluster','SYMBOL')
  # write.csv(gene_df,paste0(fp_output,i,'_IIBDGC significant genes.csv'))
  
  trans_gene = clusterProfiler::bitr(
    unique(gene_df$SYMBOL),fromType = 'SYMBOL',
    toType = 'ENTREZID',OrgDb = org.Hs.eg.db::org.Hs.eg.db)
  gene_df = gene_df %>% inner_join(trans_gene)
  
  gene =unique(gene_df$ENTREZID)
  res_GO=clusterProfiler::enrichGO(
    gene=gene,OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    keyType = "ENTREZID",
    ont='ALL', pAdjustMethod = "BH",
    pvalueCutoff=0.99,
    qvalueCutoff=0.99,
    readable =T)
  
  es = clusterProfiler::enrichKEGG(
    gene=gene,organism = 'hsa',
    # universe     = gene_all,
    # use_internal_data =T,
    pvalueCutoff = 0.99,
    qvalueCutoff =0.99)
  es=DOSE::setReadable(es, OrgDb=org.Hs.eg.db::org.Hs.eg.db,keyType='ENTREZID')
  
  
  p_go = ontology.barplot(res_GO@result,n_word = 6)+
    # labs(title = names(ego.res)[x])+
    theme(plot.title = element_text(hjust = 0.5))
  p_go
  ggsave(filename = paste0(fp_output,'IIBDGC','_',i,' Shared Merged GO.pdf'),
         width = 7.5,height = 6)
  
  
  es@result = es@result %>% filter(pvalue<0.05)
  p_kegg = enrichplot::dotplot(es,label_format=100,color='pvalue')
  p_kegg
  ggsave(filename = paste0(fp_output,'IIBDGC','_',i,' Shared Merged KEGG.pdf'),
         width = 6,height = 4.5)
  
}

