rm(list=ls())
options(stringsAsFactors = F)
library(stringr)
library(tidyverse)
library(data.table)
library(patchwork)
library(ggsci)

# fp_smr: SMR result dir
# fp_output: output dir
fp_output='./result/Fig2/'
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

# BrainMeta and Colon co-analysis ----------------------------------------------

# IBD;CD;UC cohorts selection # 
disease_list = list(
  IBD=c('ieu-a-31','categorical-20002-both_sexes-1461','finngen_R8_K11_KELAIBD'),
  CD=c('ieu-a-30','icd10-K50-both_sexes','finngen_R8_CHRONLARGE'),
  UC=c('ieu-a-32','icd10-K51-both_sexes','finngen_R8_K11_UC_STRICT2')
)

# Brain meta and colon cohorts # 
tissue_list = c('BrainMeta','Colon_Sigmoid','Colon_Transverse')

# select the combination # 
dt = expand.grid(unlist(disease_list),tissue_list)

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

dt2 = expand.grid(
  c('BrainMeta','Colon_Sigmoid','Colon_Transverse'),
  c('ieu-a-31','ieu-a-30','ieu-a-32')
)
dt2 = dt2 %>% 
  mutate(traits=paste(Var1,Var2,'smr',sep = '.'),
         Var3 = plyr::mapvalues(
           Var2,from=c('ieu-a-31','ieu-a-30','ieu-a-32'),
           to=c('IBD','CD','UC')))
  
vp_list = lapply(1:nrow(dt2),function(x){
  p_tmp = volcano_plot(df)+
    labs(x='Beta',y='-log10(p value)',color=NULL,
         title = paste0(dt2[x,'Var1'],' eQTL on ',dt2[x,'Var3']))+
    theme(plot.title = element_text(hjust = 0.5))
  return(p_tmp)  
})

wrap_plots(vp_list)+plot_layout(guides = 'collect')
ggsave(filename = paste0(fp_output,'Volcano plot list.pdf'),
       width = 12,height = 10.5)

vp_list = lapply(1:nrow(dt2),function(x){
  df=smr_list[[dt2[x,'traits']]]
  if (dt2[x,'traits']=='BrainMeta.ieu-a-30.smr') {
    df[df$b_SMR< -1,'b_SMR'] = -1
  }
  p_tmp = volcano_plot(df)+
    labs(x='Beta',y='-log10(p value)',color=NULL,
         title = paste0(dt2[x,'Var1'],' eQTL on ',dt2[x,'Var3']))+
    theme(plot.title = element_text(hjust = 0.5))
  return(p_tmp)  
})
wrap_plots(vp_list)+plot_layout(guides = 'collect')
ggsave(filename = paste0(fp_output,'Volcano plot list-2.pdf'),
       width = 12,height = 10.5)

dt3 = dt2 %>% filter(Var1=='BrainMeta')
vp_list = lapply(1:nrow(dt3),function(x){
  df=smr_list[[dt3[x,'traits']]]
  if (dt3[x,'traits']=='BrainMeta.ieu-a-30.smr') {
    df[df$b_SMR< -1,'b_SMR'] = -1
  }
  p_tmp = volcano_plot(df)+
    labs(x='Beta',y='-log10(p value)',color=NULL,
         title = paste0(dt3[x,'Var1'],' eQTL on ',dt3[x,'Var3']))+
    theme(plot.title = element_text(hjust = 0.5))
  return(p_tmp)  
})
wrap_plots(vp_list)+plot_layout(guides = 'collect')
ggsave(filename = paste0(fp_output,'Volcano plot list-3.pdf'),
       width = 12,height = 10.5/3)
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
    mutate(condition=ifelse(df>0,'Pro','Sup'))
  x = split(x$Gene,x$condition)
  return(x)
}
get_sig_terms(smr_list$`BrainMeta.ieu-a-30.smr`,adjust = F)
sig_terms = lapply(smr_list,get_sig_terms,adjust=F)

dt2 = expand.grid(
  c('BrainMeta','Colon_Sigmoid','Colon_Transverse'),
  c('ieu-a-31','ieu-a-30','ieu-a-32')
)
dt2 = dt2 %>% 
  mutate(traits=paste(Var1,Var2,'smr',sep = '.'),
         Var3 = plyr::mapvalues(
           Var2,from=c('ieu-a-31','ieu-a-30','ieu-a-32'),
           to=c('IBD','CD','UC')))

sig_terms = sig_terms[dt2$traits]

df = lapply(sig_terms,unlist,recursive=F)
df = data.frame(lengths(df)) %>% 
  rownames_to_column(var = 'traits') %>% 
  rename(Freq = `lengths.df.`) %>% 
  inner_join(dt2,by='traits')

reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}
scale_x_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_x_discrete(labels = function(x) gsub(reg, "", x), ...)
}
df %>% 
  mutate(Var3 = factor(Var3,levels=c('IBD','CD','UC'))) %>% 
  # mutate(Var3=reorder(Var3, Freq, sum,decreasing=T)) %>% 
  ggplot(aes(x=reorder_within(Var1, Freq,Var3),
             y=Freq))+
  # geom_point(aes(color=Var3))+
  geom_bar(aes(fill=Var3),stat="identity")+
  # scale_color_gradientn(colours = rev(mapal))+
  scale_fill_npg()+
  scale_x_reordered()+
  # coord_fixed(ratio = 1)+
  theme_bw()+
  coord_flip()+
  facet_wrap(~Var3,scales = 'free')+
  labs(x=NULL,y='Significant Genes ',fill=NULL)+
  theme(axis.text.x = element_text(angle = 90,hjust=1,vjust = 0.5))
ggsave(filename = paste0(fp_output,'Significant number Barplot.pdf'),
       width = 7,height = 2)

# * weighted index -------------------------------------------------------------
get_weight = function(df){
  df = df %>% 
    filter(p_SMR<0.05&p_HEIDI>0.05) %>% 
    mutate(ws = b_SMR*(-log10(p_SMR)),
           ws2 = b_SMR*(-log10(p_SMR))*(p_HEIDI))
  return(df)
}

smr_ibd = lapply(smr_list,get_weight)

dt2 = expand.grid(
  c('BrainMeta','Colon_Sigmoid','Colon_Transverse'),
  c('ieu-a-31','ieu-a-30','ieu-a-32')
)
dt2 = dt2 %>% 
  mutate(traits=paste(Var1,Var2,'smr',sep = '.'),
         Var3 = plyr::mapvalues(
           Var2,from=c('ieu-a-31','ieu-a-30','ieu-a-32'),
           to=c('IBD','CD','UC')))

df = lapply(1:nrow(dt2),function(x){
  df_tmp=smr_ibd[[dt2[x,'traits']]]
  df_tmp=colSums(abs(df_tmp[,c('ws','ws2')]))
  df_tmp=data.frame(
    ws = df_tmp[1],ws2=df_tmp[2],traits=dt2[x,'traits'])
  return(df_tmp)
})
df = Reduce(rbind,df)

df = df %>% inner_join(dt2,by='traits')

reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}
scale_x_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_x_discrete(labels = function(x) gsub(reg, "", x), ...)
}
df %>% 
  mutate(Var3 = factor(Var3,levels=c('IBD','CD','UC'))) %>% 
  ggplot(aes(x=reorder_within(Var1, ws,Var3),
             y=ws))+
  # geom_point(aes(color=Var3))+
  geom_bar(aes(fill=Var3),stat="identity")+
  # scale_color_gradientn(colours = rev(mapal))+
  scale_fill_npg()+
  scale_x_reordered()+
  # coord_fixed(ratio = 1)+
  theme_bw()+
  coord_flip()+
  facet_wrap(~Var3,scales = 'free')+
  labs(x=NULL,y='Weighted Significance Score',fill=NULL)+
  theme(axis.text.x = element_text(angle = 90,hjust=1,vjust = 0.5))
ggsave(filename = paste0(fp_output,'Weighted Significance Score-1 Barplot.pdf'),
       width = 7,height = 2)

df %>% 
  mutate(Var3 = factor(Var3,levels=c('IBD','CD','UC'))) %>% 
  ggplot(aes(x=reorder_within(Var1, ws,Var3),
             y=ws))+
  # geom_point(aes(color=Var3))+
  geom_bar(aes(fill=Var3),stat="identity")+
  # scale_color_gradientn(colours = rev(mapal))+
  scale_fill_npg()+
  scale_x_reordered()+
  # coord_fixed(ratio = 1)+
  theme_bw()+
  coord_flip()+
  facet_wrap(~Var3,scales = 'free')+
  labs(x=NULL,y='Weighted Significance Score',fill=NULL)+
  theme(axis.text.x = element_text(angle = 90,hjust=1,vjust = 0.5))
ggsave(filename = paste0(fp_output,'Weighted Significance Score-1 Barplot.pdf'),
       width = 7,height = 2)

df %>% 
  mutate(Var3 = factor(Var3,levels=c('IBD','CD','UC'))) %>% 
  ggplot(aes(x=reorder_within(Var1, ws2,Var3),
             y=ws2))+
  # geom_point(aes(color=Var3))+
  geom_bar(aes(fill=Var3),stat="identity")+
  # scale_color_gradientn(colours = rev(mapal))+
  scale_fill_npg()+
  scale_x_reordered()+
  # coord_fixed(ratio = 1)+
  theme_bw()+
  coord_flip()+
  facet_wrap(~Var3,scales = 'free')+
  labs(x=NULL,y='Weighted Significance Score',fill=NULL)+
  theme(axis.text.x = element_text(angle = 90,hjust=1,vjust = 0.5))
ggsave(filename = paste0(fp_output,'Weighted Significance Score-2 Barplot.pdf'),
       width = 7,height = 2)

