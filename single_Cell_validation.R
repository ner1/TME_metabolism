library(fgsea)
library(KEGGREST)
library(ggplot2)
library(gridExtra)
library(stringr)
library(tidytext)

# my_pathways <- gmtPathways("../data/c2.cp.kegg.v6.2.symbols.gmt")
# KEGG pathway gene set can be downloaded from MSigDB Collection
# http://www.gsea-msigdb.org/gsea/msigdb/collections.jsp

#################################### extract a list of metabolic pathways from KEGG #############
human_pathway_list = keggList("pathway", "hsa")
names(human_pathway_list) = gsub("path:hsa","",names(human_pathway_list))
human_pathway_list = data.frame(path_no = names(human_pathway_list), path_name = human_pathway_list)
human_pathway_list$path_no = as.numeric(as.character(human_pathway_list$path_no))
metabolic_pathways = human_pathway_list[which(human_pathway_list$path_no < 2000),]

metabolic_pathways$path_name = gsub(" ","_", metabolic_pathways$path_name)
metabolic_pathways$path_name = gsub("/","_", metabolic_pathways$path_name)
metabolic_pathways$path_name = gsub(",","_", metabolic_pathways$path_name)
metabolic_pathways$path_name = gsub("-","_", metabolic_pathways$path_name)
metabolic_pathways$path_name = gsub("Homo_sapiens_","_", metabolic_pathways$path_name)
metabolic_pathways$path_name = gsub("human","_", metabolic_pathways$path_name)
metabolic_pathways$path_name = gsub("\\(|\\)","", metabolic_pathways$path_name)
metabolic_pathways$path_name = gsub("___","_", metabolic_pathways$path_name)
metabolic_pathways$path_name = gsub("__","_", metabolic_pathways$path_name)
metabolic_pathways$path_name = gsub("__","_", metabolic_pathways$path_name)
metabolic_pathways$path_name = gsub("__","_", metabolic_pathways$path_name)
metabolic_pathways$path_name = gsub("_$","", metabolic_pathways$path_name)
metabolic_pathways$path_name = toupper(metabolic_pathways$path_name)
metabolic_pathways$path_name = paste("KEGG",metabolic_pathways$path_name, sep = "_")


################# Sigle Cell SKCM expression data (Tirosh et al. 2016)

# skcm <- read.delim(file = "../data/GSE72056_melanoma_single_cell_revised_v2.txt.gz",
#                                sep = "\t",
#                                stringsAsFactors = FALSE)
# cancer_cols = which(skcm[2,] == 2)
# stroma_cols = which(skcm[2,] == 1)
# colnames(skcm)[cancer_cols] = paste(colnames(skcm)[cancer_cols],
#                                                 "_Cancer",sep = "")
# colnames(skcm)[stroma_cols] = paste(colnames(skcm)[stroma_cols],
#                                                 "_Stroma",sep = "")
# skcm = skcm[,c(1,cancer_cols,stroma_cols)]
# skcm = skcm[-c(1:3),]

# ### Difference in median expression 
# result_skcm = lapply(1:length(skcm$Cell), function(x){
#   if (length(which(as.numeric(skcm[x,-1]) <= 0)) < 4514) {
#     print(x)
#     expr_cancer = as.numeric(skcm[x,grep(pattern = "Cancer",colnames(skcm))])
#     expr_stroma = as.numeric(skcm[x,grep(pattern = "Stroma",colnames(skcm))])
#     avg_cancer = median(expr_cancer)
#     sd_cancer = sd(expr_cancer)
#     avg_stroma = median(expr_stroma)
#     sd_stroma = sd(expr_stroma)
#     p.val = t.test(expr_cancer,expr_stroma)$p.value
#     y = cbind(skcm$Cell[x],avg_cancer,sd_cancer,avg_stroma,sd_stroma,p.val)
#     return(y)
#   }
#   })
# result_skcm2 = do.call(rbind, result_skcm)
# save(result_skcm2, file = "../data/validation/result_skcm2.R")

load("../data/validation/result_skcm2.R")
df_skcm = as.data.frame(result_skcm2)
df_skcm$avg_cancer = as.numeric(df_skcm$avg_cancer)
df_skcm$avg_stroma = as.numeric(df_skcm$avg_stroma)
df_skcm$cs = df_skcm$avg_cancer - df_skcm$avg_stroma
df_skcm$padj = p.adjust(df_skcm$p.val,method="BH")

#### GSEA
##############################
gsea_skcm = df_skcm[which(df_skcm$avg_cancer > 1 | df_skcm$avg_stroma > 1),]
gene_rank_skcm = gsea_skcm[order(gsea_skcm$cs, decreasing = TRUE),]
g_rank_skcm = as.numeric(gene_rank_skcm$cs)
names(g_rank_skcm) = gene_rank_skcm$V1
# fgseaRes_skcm <- fgsea(pathways = my_pathways, 
#                   stats = g_rank_skcm,
#                   minSize=15,
#                   maxSize=500,
#                   nperm=100000)
# save(fgseaRes_skcm, file = "../data/validation/fgseaRes_skcm.R")

load("../data/validation/fgseaRes_skcm.R")
fgseaRes_skcm_p = fgseaRes_skcm[fgseaRes_skcm$pathway %in% metabolic_pathways$path_name,]
fgseaRes_skcm_p$logpadj = -log10(fgseaRes_skcm_p$padj) * sign(fgseaRes_skcm_p$NES)
fgseaRes_skcm_p = fgseaRes_skcm_p[order(fgseaRes_skcm_p$logpadj, decreasing = FALSE), ]
fgseaRes_skcm_p$pathway = factor(fgseaRes_skcm_p$pathway, levels = fgseaRes_skcm_p$pathway)
fgseaRes_skcm_p$star = ""
fgseaRes_skcm_p$star[which(fgseaRes_skcm_p$padj <= 0.05 )] = "*"

bars_skcm = ggplot(data = fgseaRes_skcm_p, 
                  aes(x = pathway, y= logpadj)) + 
  geom_bar(stat = "identity", width = 0.6, fill = "#A9A9A9") +
  theme_classic() +
  labs(y = "sign(NES)*log10(p.adjusted)") +
  theme(axis.text.x = element_text(),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text = element_text(colour = "black", size = 12),
        legend.text=element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "bottom",
        plot.title = element_text(size = 12))+
  ggtitle("Cancer vs Stroma (Single-Cell SKCM)") + 
  coord_flip() + 
  geom_text(aes(label=star), colour="black", hjust = 0, vjust=0.7, size=12) + 
  geom_hline(yintercept = c(0,1.30103,-1.30103)) + 
  ylim(c(-4.1,4.1))

######################################
################# Sigle Cell HNSC expression data (Puram et al. 2017)
# 
# hnsc <- read.delim(file = "../data/GSE103322_HNSCC_all_data.txt.gz",
#                    sep = "\t",stringsAsFactors = FALSE)
# cancer_cols = which(hnsc[3,] == 1)
# stroma_cols = which(hnsc[4,] == 1)
# colnames(hnsc)[cancer_cols] = paste(colnames(hnsc)[cancer_cols],
#                                                 "_Cancer",sep = "")
# colnames(hnsc)[stroma_cols] = paste(colnames(hnsc)[stroma_cols],
#                                                 "_Stroma",sep = "")
# hnsc = hnsc[,c(1,cancer_cols,stroma_cols)]
# hnsc = hnsc[-c(1:5),]

# ### Difference in median expression 
# result_hnsc = lapply(1:length(hnsc$X), function(x){
#   if (length(which(as.numeric(hnsc[x,-1]) <= 0)) < 5578) {
#     print(x)
#     expr_cancer = as.numeric(hnsc[x,grep(pattern = "Cancer",colnames(hnsc))])
#     expr_stroma = as.numeric(hnsc[x,grep(pattern = "Stroma",colnames(hnsc))])
#     avg_cancer = median(expr_cancer)
#     sd_cancer = sd(expr_cancer)
#     avg_stroma = median(expr_stroma)
#     sd_stroma = sd(expr_stroma)
#     p.val = t.test(expr_cancer,expr_stroma)$p.value
#     y = cbind(hnsc$X[x],avg_cancer,sd_cancer,avg_stroma,sd_stroma,p.val)
#     return(y) 
#   }
# })
# result_hnsc2 = do.call(rbind, result_hnsc)
# save(result_hnsc2,file = "../data/validation/result_hnsc2.R")

load(file = "../data/validation/result_hnsc2.R")
df_hnsc = as.data.frame(result_hnsc2)
df_hnsc$avg_cancer = as.numeric(df_hnsc$avg_cancer)
df_hnsc$avg_stroma = as.numeric(df_hnsc$avg_stroma)
df_hnsc$cs = df_hnsc$avg_cancer - df_hnsc$avg_stroma
df_hnsc$padj = p.adjust(df_hnsc$p.val,method="BH")

df_hnsc$V1 = gsub(pattern = "[']",x = df_hnsc$V1,replacement = "")

#### GSEA
##############################
gsea_hnsc = df_hnsc[which(df_hnsc$avg_cancer > 1 | df_hnsc$avg_stroma > 1),]
gene_rank_hnsc = gsea_hnsc[order(gsea_hnsc$cs, decreasing = TRUE),]
g_rank_hnsc = as.numeric(gene_rank_hnsc$cs)
names(g_rank_hnsc) = gene_rank_hnsc$V1
# fgseaRes_hnsc <- fgsea(pathways = my_pathways, 
#                   stats = g_rank_hnsc,
#                   minSize=15,
#                   maxSize=500,
#                   nperm=100000)
# save(fgseaRes_hnsc, file = "../data/validation/fgseaRes_hnsc.R")

load(file = "../data/validation/fgseaRes_hnsc.R")
fgseaRes_hnsc_p = fgseaRes_hnsc[fgseaRes_hnsc$pathway %in% metabolic_pathways$path_name,]

fgseaRes_hnsc_p$logpadj = -log10(fgseaRes_hnsc_p$padj) * sign(fgseaRes_hnsc_p$NES)
fgseaRes_hnsc_p = fgseaRes_hnsc_p[order(fgseaRes_hnsc_p$logpadj, decreasing = FALSE), ]
fgseaRes_hnsc_p$pathway = factor(fgseaRes_hnsc_p$pathway, levels = fgseaRes_hnsc_p$pathway)
fgseaRes_hnsc_p$star = ""
fgseaRes_hnsc_p$star[which(fgseaRes_hnsc_p$padj <= 0.05 )] = "*"


bars_hnsc = ggplot(data = fgseaRes_hnsc_p, 
                  aes(x = pathway, y= logpadj)) + 
  geom_bar(stat = "identity", width = 0.6, fill = "#A9A9A9") +
  theme_classic() +
  labs(y = "sign(NES)*log10(p.adjusted)") +
  theme(axis.text.x = element_text(),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text = element_text(colour = "black", size = 12),
        legend.text=element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "bottom",
        plot.title = element_text(size = 12))+
  ggtitle("Cancer vs Stroma (Single-Cell HNSC)") + 
  coord_flip() + 
  geom_text(aes(label=star), colour="black", hjust = 0, vjust=0.7, size=12) + 
  geom_hline(yintercept = c(0,1.30103,-1.30103)) + 
  ylim(c(-4.1,4.1))


############# plot map (Figure 4####################
grid.arrange(bars_skcm,bars_hnsc,nrow = 1)

fgsea_all = merge(fgseaRes_hnsc_p, fgseaRes_skcm_p, all = TRUE)
fgsea_all = fgsea_all[,c("pathway","logpadj_hnsc","logpadj_skcm",
                         "star_skcm","star_hnsc")]

fgsea_all = reshape2::melt(fgsea_all)

bars_gsea = ggplot(data = fgsea_all, 
                     aes(x = reorder_within(pathway, 
                                            cs, tumor_type), 
                         y= cs, 
                         # fill = factor(col)
                     )) + 
  geom_bar(stat = "identity"
           # color = "black"
  ) +
  theme_classic() +
  # scale_fill_manual(values = c("#969696","#737373",
  # "#525252","#252525","#000000")) + 
  facet_wrap(~tumor_type, scales = "free_x") + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_blank(),
        legend.text=element_text(size = 12),
        legend.title = element_text(size = 12),
        # legend.position = "bottom",
        # plot.title = element_text(size = 12)
  ) + 
  ylab("Expr (Cancer - Stroma)") + 
  ylim(c(-6.5,6.5))