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
fp_output='./result/Fig5-part1/'

# seurat pipeline --------------------------------------------------------------
sce_obj = readRDS(paste0(fp_sc,'CD colon scRNA.rds'))
unique(sce_obj$disease)
sce_obj = subset(sce_obj,subset=disease=="Crohn disease")

# filter Cells #
sce_obj[["percent.mt"]] <- PercentageFeatureSet(sce_obj, pattern = "^MT-")
VlnPlot(sce_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)

min(sce_obj$nFeature_RNA)
min(sce_obj$nCount_RNA)
# Filter cells referring to paper #
# sce_obj <- subset(x=sce_obj, 
#                   subset = nFeature_RNA>100 & nFeature_RNA < 6000& 
#                     nCount_RNA>100&percent.mt < 30 )

Idents(sce_obj)='cell_type'

# DotPlot(sce_obj,features = pro_gene)
# DotPlot(sce_obj,features = sup_gene)
sce_obj = JoinLayers(sce_obj)

# Normalize and scale the data #
sce_obj <- sce_obj %>% 
  NormalizeData() %>% 
  ScaleData() %>% 
  FindVariableFeatures() %>% 
  RunPCA()

ElbowPlot(sce_obj,ndims=40) # Visualize variance along each component

n.pcs <- 30
sce_obj <- sce_obj %>% 
  FindNeighbors(dims = 1:n.pcs, verbose = T) %>%  # Construct Neighbor graph 
  RunTSNE(dims = 1:n.pcs) %>% 
  RunUMAP(dims = 1:n.pcs)

saveRDS(sce_obj,paste0(fp_output,'CD colon scRNA.rds'))

# UMAP by plot1cell ------------------------------------------------------------
library(plot1cell)
sce_obj = readRDS(paste0(fp_output,'CD colon scRNA.rds'))
colnames(sce_obj@meta.data)
Idents(sce_obj)='cell_type'
table(sce_obj$disease,sce_obj$Type)

circ_data <- prepare_circlize_data(sce_obj, scale = 0.8 )
set.seed(1234)

col.celltype = colorRampPalette(pal_npg()(10))(length(levels(sce_obj)))
# col.disease = c('#61C0BF','#F64E60')
col.status = c('#61C0BF','#E76F51','#E9D985')[1:2]


###plot and save figures
pdf(file =  paste0(fp_output,'circlize_plot.pdf'), width =12, height = 12)
plot_circlize(circ_data,do.label = T, pt.size = 0.01, col.use = col.celltype ,
              bg.color = 'white',kde2d.n = 500, repel = T, label.cex = 2)
# add_track(circ_data, group = "disease", colors = col.disease, track_num = 2) ## can change it to one of the columns in the meta data of your seurat object
add_track(circ_data, group = "Type",colors = col.status, track_num = 3) ## can change it to one of the columns in the meta data of your seurat object
dev.off()

par(mfrow=c(1,1))
png(filename =  paste0(fp_output,'circlize_plot.png'),width = 960,height = 960)
plot_circlize(circ_data,do.label = T, pt.size = 0.01, col.use = col.celltype ,
              bg.color = 'white',kde2d.n = 500, repel = T, label.cex = 2)
# add_track(circ_data, group = "disease", colors = col.disease, track_num = 2) ## can change it to one of the columns in the meta data of your seurat object
add_track(circ_data, group = "Type",colors = col.status, track_num = 3) ## can change it to one of the columns in the meta data of your seurat object
dev.off()


names(col.celltype)=levels(sce_obj)
p1 = DimPlot(sce_obj,pt.size = 0.1,reduction = "umap", label = TRUE,
             group.by = 'cell_type',split.by = 'Type')+
  scale_color_manual(values = col.celltype)
p1
ggsave(filename = paste0(fp_output,'umap_plot-2.pdf'), width =11, height = 4,plot = p1)

p2 = DimPlot(sce_obj,pt.size = 0.1,reduction = "umap", label = TRUE,
             group.by = 'cell_type')+
  scale_color_manual(values = col.celltype)
p2
ggsave(filename = paste0(fp_output,'umap_plot.pdf'), width =10, height = 4,plot = p2)

# FET --------------------------------------------------------------------------
# ** Fisher exact test ---------------------------------------------------------
pro_gene = c("NAAA","CCDC88B","NAGLU")
sup_gene = c("CWC15")
gene.can2 = intersect(c(pro_gene,sup_gene),rownames(sce_obj))
colnames(sce_obj@meta.data)
fet.data = FetchData(sce_obj,vars = c('cell_type',gene.can2))
colnames(fet.data)[1] = 'cell_type'

index = apply(fet.data[,-1], 2, function(x) sum(x>0)/length(x))
summary(index)
cell_type = unique(fet.data$cell_type)

library(future.apply)
# plan(multiprocess)
plan(multisession,workers =6)
FET.res = lapply(gene.can2, function(gene){
  future_lapply(cell_type, function(i){
    # gene=gene.can2[1]
    # i = cell_type[1]
    df.tmp = fet.data[,c('cell_type',gene)]
    df.tmp$cell_type = ifelse(df.tmp$cell_type ==i, i,'other')
    df.tmp = df.tmp %>% group_by(cell_type) %>% 
      summarise(epr=sum(get(gene)>0),
                non_exp = sum(get(gene)<= 0))
    df.tmp = df.tmp %>% column_to_rownames(var = 'cell_type')
    
    fet.tmp = fisher.test(df.tmp)
    return(fet.tmp$p.value)
  })
})

names(FET.res) = gene.can2

# ** p.adjust ------------------------------------------------------------------
FET.res2 = lapply(FET.res, function(x){
  p_value = Reduce(c,x)
  return(p_value)
}) #merge the p value

FET.res2 = lapply(FET.res2, function(x){
  p_adj = p.adjust(x,method = 'fdr')
  return(p_adj)
}) # fdr p adjust 

es = lapply(FET.res2, function(x){
  x[x>0.05]=1
  es.tmp = -log10(x)
}) # -log10(p.adj)
es.df = Reduce(rbind.data.frame,es)
rownames(es.df)=names(es);colnames(es.df)=cell_type

# inf problem # 
index = apply(es.df,1,function(x) sum(is.infinite(x)))
table(index)

es2 = lapply(es, function(x){
  # x=es[[1]]
  inf.index = is.infinite(x)
  if (sum(inf.index)>0) {
    # x[!is.infinite(x)]=0 # 先写这个,不然会把1变成0
    # x[is.infinite(x)]=1
    # x[is.infinite(x)]=1.5*max(x[!is.infinite(x)])
    x[is.infinite(x)]=max(unlist(es)[!is.infinite(unlist(es))])*1.5
  }
  return(x)
})

es.df = Reduce(rbind.data.frame,es2)
rownames(es.df)=names(es2);colnames(es.df)=cell_type

saveRDS(es.df,file = paste0(fp_output,'gene cell FET mat.rds'))
save(es.df,gene.can2,file = paste0(fp_output,'gene cell FET mat.rdata'))


# heatmap ----------------------------------------------------------------------
library(ComplexHeatmap)
library(circlize)
library(ggsci)
es.df=readRDS(paste0(fp_output,'gene cell FET mat.rds'))
load(paste0(fp_output,'gene cell FET mat.rdata'))

ee = scale(t(es.df))
ee = ee[,apply(ee,2,function(x) !all(is.na(x)))]

# kmeans.result <- kmeans(t(ee), 6)
# kmeans_df <- data.frame(kmeans_class=kmeans.result$cluster)
# kmeans_df= kmeans_df %>% arrange(kmeans_class) %>% 
#   mutate(kmeans_class=as.character(kmeans_class))
# 
# Top = HeatmapAnnotation(Cluster=kmeans_df$kmeans_class,
#                         annotation_legend_param=list(labels_gp = gpar(fontsize = 10),border = T,
#                                                      title_gp = gpar(fontsize = 10,fontface = "bold"),
#                                                      ncol=1),
#                         border = T,
#                         col=list(Cluster = setNames(pal_npg(alpha=0.7)(length(unique(kmeans_df$kmeans_class))),
#                                                     sort(unique(kmeans_df$kmeans_class)))
#                         ),
#                         show_annotation_name = TRUE,
#                         annotation_name_side="left",
#                         annotation_name_gp = gpar(fontsize = 10))

ht = Heatmap(ee,name='Z-score',
             # top_annotation = Top,
             cluster_rows = T,
             col=colorRamp2(c(-2,0,2),c('#21b6af','white','#eeba4d')),#49b0d9
             color_space = "RGB",
             cluster_columns = T,border = T,
             row_order=NULL,
             row_names_side = 'left',
             column_order=NULL,
             show_column_names = T,
             row_names_gp = gpar(fontsize = 9),
             # column_split = kmeans_df$kmeans_class,
             # column_km = 12,
             gap = unit(1, "mm"),
             column_title = NULL,
             column_title_gp = gpar(fontsize = 10),
             show_heatmap_legend = TRUE,
             heatmap_legend_param=list(labels_gp = gpar(fontsize = 10), border = T,
                                       title_gp = gpar(fontsize = 10, fontface = "bold")),
             column_gap = unit(2,'mm')
) 
ht=draw(ht)
column_order(ht)
gene_cluster <- column_order(ht)
clu_df <- lapply(names(gene_cluster), function(i){
  out <- data.frame(coordinates = colnames(ee)[gene_cluster[[i]]],
                    Cluster = paste0("cluster", i), stringsAsFactors = FALSE)
  return(out)
}) %>%  
  do.call(rbind, .)

clu_df = clu_df %>% column_to_rownames(var = 'coordinates')
getPalette = colorRampPalette(pal_npg(alpha=0.7)(10)) 
Top = HeatmapAnnotation(Cluster=clu_df$Cluster,
                        annotation_legend_param=list(labels_gp = gpar(fontsize = 10),border = T,
                                                     title_gp = gpar(fontsize = 10,fontface = "bold"),
                                                     ncol=1),
                        border = T,
                        col=list(Cluster = setNames(getPalette(length(unique(clu_df$Cluster))),
                                                    sort(unique(clu_df$Cluster)))
                        ),
                        show_annotation_name = TRUE,
                        annotation_name_side="left",
                        annotation_name_gp = gpar(fontsize = 10))
ee = ee[,rownames(clu_df)]
Heatmap(ee,name='Z-score',
        top_annotation = Top,
        cluster_rows = T,
        col=colorRamp2(c(-2,0,2),c('#21b6af','white','#eeba4d')),#49b0d9
        color_space = "RGB",
        cluster_columns = F,border = T,
        row_order=NULL,
        row_names_side = 'left',
        column_order=NULL,
        show_column_names = T,
        row_names_gp = gpar(fontsize = 9),
        column_split = clu_df$Cluster,
        # column_km = 3,
        gap = unit(1, "mm"),
        column_title = NULL,
        column_title_gp = gpar(fontsize = 10),
        show_heatmap_legend = TRUE,
        heatmap_legend_param=list(labels_gp = gpar(fontsize = 10), border = T,
                                  title_gp = gpar(fontsize = 10, fontface = "bold")),
        column_gap = unit(2,'mm')
) 
export::graph2pdf(
  file=paste0(fp_output,'Scaled enrichment score Heatmap.pdf'),
  width=6,height=8
)
save(clu_df,ht,file = paste0(fp_output,'Cluster scRNA.rda'))

# expression-percentage  ------------------------------------------------------- 
ep =future_lapply(cell_type, function(i){
  # i = cell_type[1]
  df.tmp = fet.data %>% dplyr::filter(cell_type==i) %>% 
    dplyr::select(-cell_type)
  ep.tmp = apply(df.tmp, 2, function(x) sum(x>0)/length(x))
  ep.tmp =as.data.frame(ep.tmp)
  ep.tmp = ep.tmp %>% rownames_to_column(var = 'gene') 
  colnames(ep.tmp)[2]=i
  return(ep.tmp)
})

ep.df = Reduce(inner_join,ep)
ep.df = ep.df %>% 
  pivot_longer(cols = 2:ncol(ep.df),
               names_to = 'cell_type',
               values_to = 'Expression_percent')
es.df = es.df %>% t() %>% as.data.frame() %>% 
  rownames_to_column(var ='cell_type') %>% 
  pivot_longer(cols = 2:(length(gene.can2)+1),
               names_to = 'gene',
               values_to = 'Enrichment_score')


p.df = es.df %>% inner_join(ep.df)
ylgn_palette <- brewer.pal(9, "YlGnBu")
p.df %>% 
  mutate(Expression_percent=100*Expression_percent) %>% 
  mutate(gene=reorder(gene, Enrichment_score, sum,decreasing=T)) %>% 
  mutate(cell_type=reorder(cell_type, Expression_percent, sum)) %>% 
  ggplot(aes(x=gene,y=cell_type))+
  geom_point(aes(size=Expression_percent,
                 color=Enrichment_score),
             shape=15)+
  # scale_colour_gradient(low='lightgrey',high ='#B22222' )+ # '#D2413F'
  # scale_colour_gradient2(low='lightgrey',mid = '#FFC0CB',high = '#B22222')+
  scale_color_gradientn(colors =  brewer.pal(9, "YlGnBu"))+
  coord_fixed(ratio = 1)+theme_bw()+
  labs(x=NULL,y='Cell type',color='Enrichment score',size='Expression percent',
       title = 'Core Gene sub-cell Enrichment score')+
  theme(axis.text.x = element_text(angle = 30,hjust=1,vjust = 1),
        plot.title = element_text(hjust=0.5))
ggsave(filename = paste0(fp_output,'FET of core gene sub-cell enr_score.pdf'),
       width=6,height=8)


# ** Average expression dot plot manu ---------------------------------------------

ave_df = AverageExpression(sce_obj,features =gene.can2,
                           assays = 'RNA',group.by = 'cell_type',
                           # slot ='data'
                           slot = 'scale.data'
)
ave_df = as.matrix(ave_df$RNA)
ave_df = ave_df %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene') %>% 
  pivot_longer(cols=2:(length(cell_type)+1),
               names_to = 'cell_type',
               values_to = 'Average')

p.df = es.df %>% inner_join(ep.df)
p.df=ave_df %>% 
  mutate(cell_type=str_replace(cell_type,'\\-','_')) %>% 
  inner_join(p.df)

p.df %>% 
  mutate(gene=reorder(gene, Expression_percent, sum,decreasing=T)) %>% 
  mutate(cell_type=reorder(cell_type, Expression_percent, sum)) %>% 
  ggplot(aes(x=gene,y=cell_type))+
  geom_point(aes(size=Expression_percent*100,
                 color=Average))+
  scale_colour_gradient(low='lightgrey',high = 'royalblue')+
  # scale_color_viridis_c()+
  # scale_color_gradientn(colors =  brewer.pal(9, "YlGnBu"))+
  coord_fixed(ratio = 1)+theme_bw()+
  labs(x=NULL,y='Cell type',color='Average \nExpression',
       size='Expression \npercent',
       title = 'Core Genes sub-cell Average Expression')+
  theme(axis.text.x = element_text(angle = 45,hjust=1,vjust = 1),
        # panel.background = element_rect(fill = "#E6E6E6"),
        plot.title = element_text(hjust=0.5)
  )
ggsave(filename = paste0(fp_output,'FET of core gene sub-cell ave_exp-2.pdf'),
       width=6,height=8)

# dotplot ----------------------------------------------------------------------
mapal = colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(256)
reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}
scale_x_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_x_discrete(labels = function(x) gsub(reg, "", x), ...)
}
p.df %>% 
  mutate(Expression_percent=Expression_percent*100) %>% 
  mutate(gene=reorder(gene, Enrichment_score, sum,decreasing=T)) %>% 
  # mutate(cell_type=reorder(cell_type, Expression_percent, sum)) %>% 
  ggplot(aes(x=reorder_within(cell_type, Enrichment_score,gene),
             y=Enrichment_score))+
  geom_point(aes(size=Expression_percent,
                 color=Enrichment_score))+
  scale_color_gradientn(colours = rev(mapal))+
  scale_x_reordered()+
  # coord_fixed(ratio = 1)+
  theme_bw()+
  coord_flip()+
  facet_wrap(~gene,scales = 'free')+
  labs(x=NULL,y='Enrichment Score',color='Enrichment score',size='Expression percent')+
  theme(plot.title = element_text(hjust=0.5))


d <- p.df %>% 
  mutate(Expression_percent=Expression_percent*100) %>% 
  mutate(gene=reorder(gene, Enrichment_score, sum,decreasing=T)) %>% 
  ungroup() %>%   # As a precaution / handle in a separate .grouped_df method
  arrange(gene, Enrichment_score,Expression_percent) %>%   # arrange by facet variables and continuous values
  mutate(rn = row_number()) # Add a row number variable

d %>% 
  ggplot(aes(x=rn,
             y=Enrichment_score))+
  geom_point(aes(size=Expression_percent,
                 color=Enrichment_score))+
  scale_color_gradientn(colours = rev(mapal))+
  scale_x_continuous(  # This handles replacement of .r for x
    breaks = d$rn,     # notice need to reuse data frame
    labels = d$cell_type
  )+
  theme(axis.text.x = element_text(angle = 90,hjust=1))

hist(d$Enrichment_score)

d %>% 
  group_by(gene) %>% 
  top_n(n=5,wt=Enrichment_score) %>% 
  ggplot(aes(x=rn,
             y=Enrichment_score))+
  geom_point(aes(size=Expression_percent,
                 color=Enrichment_score))+
  scale_color_gradientn(colours = rev(mapal))+
  theme_bw()+
  coord_flip()+
  # facet_wrap(~gene,scales = 'free')+
  facet_grid(gene~.,scales = 'free')+
  scale_x_continuous(  # This handles replacement of .r for x
    breaks = d$rn,     # notice need to reuse data frame
    labels = d$cell_type
  )+
  labs(x=NULL,y='Enrichment Score',color='Enrichment score',size='Expression percent')+
  theme(plot.title = element_text(hjust=0.5))
ggsave(filename = paste0(fp_output,'Top5 FET of core gene sub-cell enr_score dotplot.pdf'),
       width = 4,height = 4)

gene_sub = d %>% 
  filter(Enrichment_score>400) %>% 
  pull(gene) %>% unique()
gene_sub = as.character(gene_sub)

d %>% 
  filter(gene%in%gene_sub) %>% 
  group_by(gene) %>% 
  top_n(n=3,wt=Enrichment_score) %>% 
  ungroup() %>% 
  # filter(Enrichment_score>400) %>%
  ggplot(aes(x=rn,
             y=Enrichment_score))+
  geom_point(aes(size=Expression_percent,
                 color=Enrichment_score))+
  scale_color_gradientn(colours = rev(mapal))+
  theme_bw()+
  coord_flip()+
  # facet_wrap(~gene,scales = 'free')+
  facet_grid(gene~.,scales = 'free')+
  scale_x_continuous(  # This handles replacement of .r for x
    breaks = d$rn,     # notice need to reuse data frame
    labels = d$cell_type
  )+
  labs(x=NULL,y='Enrichment Score',color='Enrichment score',size='Expression percent')+
  theme(plot.title = element_text(hjust=0.5))

ggsave(filename = paste0(fp_output,'FET of core gene sub-cell enr_score dotplot top3.pdf'),
       width = 4,height = 20)

DotPlot(subset(sce_obj,subset=cell_type=='enterocyte'),
        group.by = 'Type',features = "NAAA")
DotPlot(subset(sce_obj,subset=cell_type=='plasma cell'),
        group.by = 'Type',features = "NAAA")
DotPlot(subset(sce_obj,subset=cell_type=='macrophage'),
        group.by = 'Type',features = "NAAA")
DotPlot(subset(sce_obj,subset=cell_type=='CD4-positive, alpha-beta T cell'),
        group.by = 'Type',features = "NAAA")
DotPlot(subset(sce_obj,subset=cell_type=='goblet cell'),
        group.by = 'Type',features = "NAAA")


DotPlot(subset(sce_obj,subset=cell_type=='epithelial cell'),
        group.by = 'Type',features = "CWC15")

DotPlot(subset(sce_obj,subset=cell_type=='plasma cell'),
        group.by = 'Type',features = "NAGLU")

DotPlot(subset(sce_obj,subset=cell_type=='fibroblast'),
        group.by = 'Type',features = "CCDC88B")
