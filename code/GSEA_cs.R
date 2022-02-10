library(reshape2)
library(cowplot)
library(ggrepel)
library(ggtern)
library(grid)
library(gridExtra)
library(gtable)
library(KEGGREST)
library(fgsea)

# my_pathways <- gmtPathways("C:/Projects/GSEA/c2.cp.kegg.v6.2.symbols.gmt")
# 
# #################################### extract a list of metabolic pathways from KEGG #############
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

###################

ExprMetabolic = read.table(file = '../data/Expression_Model_All.dat',
                       sep = ',', header = TRUE, stringsAsFactors = FALSE)

ExprMetabolic$ReconGeneName = NULL
ExprMetabolic = ExprMetabolic[!duplicated(ExprMetabolic),]

# Genes in KEGG oxphos with ALIAS
kegg_ab = read.csv(file = "../data/kegg_oxphos_genes_ab.txt",
                   header = FALSE, stringsAsFactors = FALSE)
colnames(kegg_ab) = c("kegg_gene","gene")
ExprMetabolic = merge(ExprMetabolic,kegg_ab, by = "gene", all.x = TRUE)
ExprMetabolic$kegg_gene[which(is.na(ExprMetabolic$kegg_gene))] = ExprMetabolic$gene[which(is.na(ExprMetabolic$kegg_gene))]
ExprMetabolic$gene = ExprMetabolic$kegg_gene

C = ExprMetabolic[,c("gene","cancer_type","cancer")]
C = reshape2::dcast(C, gene ~ cancer_type, value.var = "cancer")
C$Cm = apply(C[,-1], 1, FUN = function(x) median(x))

S = ExprMetabolic[,c("gene","stroma","cancer_type")]
S = reshape2::dcast(S, gene ~ cancer_type, value.var = "stroma")
S$Sm = apply(S[,-1], 1, FUN = function(x) median(x))

N = ExprMetabolic[,c("gene","normal","cancer_type")]
N = reshape2::dcast(N, gene ~ cancer_type, value.var = "normal")
N$Nm = apply(N[,-1], 1, FUN = function(x) median(x))

combinedMedian = merge(N,C,by = "gene")
combinedMedian = combinedMedian[,c("gene","Nm","Cm")]
combinedMedian = merge(combinedMedian, S, by = "gene" )
combinedMedian = combinedMedian[,c("gene","Nm","Cm","Sm")]

expr = combinedMedian
expr$CS = expr$Cm-expr$Sm
gene_rank = expr[order(expr$CS, decreasing = TRUE),]
g_rank = as.numeric(gene_rank$CS)
names(g_rank) = gene_rank$gene
# 
# fgseaRes <- fgsea(pathways = my_pathways,
#                   stats = g_rank,
#                   minSize=15,
#                   maxSize=500,
#                   nperm=100000)
# fgseaRes$logpadj = -log10(fgseaRes$padj) * sign(fgseaRes$NES)
# 
# fgseaRes$leadingEdge = NULL
# write.csv(x = fgseaRes, file = "../data/GSEA_KEGG/fgseaRes_CS_all.csv")

fgseaRes = read.csv(file = "../data/GSEA_KEGG/fgseaRes_CS_all.csv")

fgseaRes = fgseaRes[which(fgseaRes$padj < 0.05),]

fgseaRes = fgseaRes[fgseaRes$pathway %in% metabolic_pathways$path_name,]

for_individual_bars = fgseaRes$pathway[order(fgseaRes$NES, decreasing = TRUE)]

cs_df = fgseaRes[fgseaRes$pathway %in% for_individual_bars,]
cs_df = cs_df[,c("pathway","logpadj","padj")]
rownames(cs_df) = NULL
cs_df = cs_df[order(cs_df$logpadj, decreasing = FALSE), ]
cs_df$pathway = factor(cs_df$pathway, levels = cs_df$pathway)

cols = c("cancer" = "#b2182b", "stroma" = "#2166ac")

cs_df$star = ""
cs_df$star[which(cs_df$padj <= 0.05 )] = "*"

bars_all = ggplot(data = cs_df, aes(x = pathway, y= logpadj)) + 
  geom_bar(stat = "identity", width = 0.6, fill = "#A9A9A9") +
  theme_classic() +
  labs(y = "sign(NES)*log10(p.adjusted)") +
  theme(axis.text.x = element_text(),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_blank(),
        axis.text = element_text(colour = "black", size = 7),
        legend.text=element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.position = "bottom",
        plot.title = element_text(size = 7))+
  ggtitle("Cancer vs Stroma (PANCAN)") + 
  coord_flip() + 
  geom_text(aes(label=star), colour="black", hjust = 0, vjust=0.7, size=6) + 
  geom_hline(yintercept = c(0,1.30103,-1.30103))

# enrich_cs_pan = plotEnrichment(my_pathways[["KEGG_OXIDATIVE_PHOSPHORYLATION"]],
#                g_rank) + labs(title="Oxidative Phosphorylation") +
#   theme_classic() +
#   theme(axis.title = element_text(colour = "black", size = 7),
#         axis.text = element_text(colour = "black", size = 7),
#         title = element_text(colour = "black", size = 7))

######################### cancer - stroma individual #################################

ExprMetabolic$CS = ExprMetabolic$cancer - ExprMetabolic$stroma
cancer_type = unique(ExprMetabolic$cancer_type)

# fgseaRes_ind = vector(length = 20, "list")
# 
# for (ct in 1:length(cancer_type)) {
#   tmp = ExprMetabolic[which(ExprMetabolic$cancer_type == cancer_type[[ct]]),]
#   gene_rank = tmp[order(tmp$CS, decreasing = TRUE),]
#   g_rank = as.numeric(gene_rank$CS)
#   names(g_rank) = gene_rank$gene
# 
#   temp_res <- fgsea(pathways = my_pathways,
#                     stats = g_rank,
#                     minSize=15,
#                     maxSize=500,
#                     nperm=100000)
# 
#   temp_res$cancer_type = cancer_type[[ct]]
#   fgseaRes_ind[[ct]] <- temp_res
#   print(ct)
# }
# 
# fgseaRes_ind = do.call(rbind,fgseaRes_ind)
# fgseaRes_ind$logpval = -log10(fgseaRes_ind$pval)*sign(fgseaRes_ind$NES)
# 
# save(fgseaRes_ind, file = "../data/CCLE/GSEA/KEGG_all_ct/fgseaRes_CS_ind.RData")
# 
# fgseaRes_ind$leadingEdge = NULL
# write.csv(x = fgseaRes_ind, file = "../data/GSEA_KEGG/fgseaRes_CS_ind.csv", row.names = FALSE)

fgseaRes_ind = read.csv(file = "../data/GSEA_KEGG/fgseaRes_CS_ind.csv",
                           stringsAsFactors = FALSE)

cs_df_ind = fgseaRes_ind[which(fgseaRes_ind$pathway %in% for_individual_bars),]
cs_df_ind = cs_df_ind[which(cs_df_ind$pathway %in% for_individual_bars),
                        c("pathway","cancer_type","NES","pval","padj","logpval")]

cb = cs_df_ind[which(cs_df_ind$pathway == "KEGG_OXIDATIVE_PHOSPHORYLATION"),]
cb$cancer_type <- factor( cb$cancer_type, 
                          levels = cb$cancer_type[order(cb$logpval, decreasing = FALSE)])

cb$star = ""
cb$star[which(cb$pval <= 0.005 )] = "*"

cb$leg = "cancer vs stroma"

bars_cs_ind = ggplot(cb,aes(cancer_type, logpval, fill = leg)) + 
    geom_bar(stat = "identity", 
             position = position_dodge(width = 0.5)) +
    theme_classic() +
    scale_fill_manual(values = c("#A9A9A9")) +
    labs(y = "sign(NES)*log10(pval)") +
  theme(legend.position="bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.box.spacing = unit(0.01, "cm"),
        legend.key.size = unit(0.4, "cm"),
        legend.key.width = unit(0.4,"cm"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(colour = "black", size = 7),
        axis.text = element_text(colour = "black", size = 7)) + 
  coord_flip() + 
  geom_text(aes(label=star), colour="black", hjust = 0, vjust=0.7, size=3) +
  geom_hline(yintercept = c(0,1.30103,-1.30103))

####################### individual bars #################################
######################### cancer #################################
ExprMetabolic$CN = ExprMetabolic$cancer - ExprMetabolic$normal

# fgseaRes_cancer = vector(length = 20, "list")
# 
# for (ct in 1:length(cancer_type)) {
#   tmp = ExprMetabolic[which(ExprMetabolic$cancer_type == cancer_type[[ct]]),]
#   gene_rank = tmp[order(tmp$CN, decreasing = TRUE),]
#   g_rank = as.numeric(gene_rank$CN)
#   names(g_rank) = gene_rank$gene
# 
#   temp_res <- fgsea(pathways = my_pathways,
#                     stats = g_rank,
#                     minSize=15,
#                     maxSize=500,
#                     nperm=100000)
# 
#   temp_res$cancer_type = cancer_type[[ct]]
#   fgseaRes_cancer[[ct]] <- temp_res
#   print(cancer_type[[ct]])
# }
# 
# save(fgseaRes_cancer, file = "../data/CCLE/GSEA/KEGG_all_ct/fgseaRes_cancer.RData")
# 
# cancer = do.call(rbind,fgseaRes_cancer)
# cancer$leadingEdge = NULL
# write.csv(x = cancer, file = "../data/GSEA_KEGG/fgseaRes_CN_ind.csv", 
          # row.names = FALSE)

cancer = read.csv(file = "../data/GSEA_KEGG/fgseaRes_CN_ind.csv", stringsAsFactors = FALSE)

cancer$logpval = -log10(cancer$pval)*sign(cancer$NES)
colnames(cancer)[-c(1,8)] <- paste(colnames(cancer)[-c(1,8)],"cancer", sep = "_")

###################### stroma #######################
ExprMetabolic$SN = ExprMetabolic$stroma - ExprMetabolic$normal

# fgseaRes_stroma = vector(length = 20, "list")
# 
# for (ct in 1:length(cancer_type)) {
#   tmp = ExprMetabolic[which(ExprMetabolic$cancer_type == cancer_type[[ct]]),]
#   gene_rank = tmp[order(tmp$SN, decreasing = TRUE),]
#   g_rank = as.numeric(gene_rank$SN)
#   names(g_rank) = gene_rank$gene
# 
#   temp_res <- fgsea(pathways = my_pathways,
#                     stats = g_rank,
#                     minSize=15,
#                     maxSize=500,
#                     nperm=100000)
# 
#   temp_res$cancer_type = cancer_type[[ct]]
#   fgseaRes_stroma[[ct]] <- temp_res
#   print(ct)
# }
# 
# save(fgseaRes_stroma, file = "../data/CCLE/GSEA/KEGG_all_ct/fgseaRes_stroma.RData")
# 
# stroma = do.call(rbind,fgseaRes_stroma)
# stroma$leadingEdge = NULL
# write.csv(x = stroma, file = "../data/GSEA_KEGG/fgseaRes_SN_ind.csv",
# row.names = FALSE)

stroma = read.csv(file = "../data/GSEA_KEGG/fgseaRes_SN_ind.csv",
                  stringsAsFactors = FALSE)

stroma$logpval = -log10(stroma$pval)*sign(stroma$NES)
colnames(stroma)[-c(1,8)] <- paste(colnames(stroma)[-c(1,8)],"stroma", sep = "_")

################ merge both and compare #########################
cs_n = merge(cancer,stroma)
cs_n = cs_n[cs_n$pathway %in% metabolic_pathways$path_name,]
cs_n = cs_n[,c("pathway","cancer_type",
               "NES_stroma","NES_cancer",
               "pval_stroma","pval_cancer",
               "logpval_cancer","logpval_stroma")]
cs_n = reshape2::melt(cs_n)
cs_n$cell = "cell"
cs_n$cell[grep("_cancer", cs_n$variable)] = "cancer vs normal"
cs_n$cell[grep("_stroma", cs_n$variable)] = "stroma vs normal"
cs_n$val = "v"
cs_n$val[grep("NES_", cs_n$variable)] = "NES"
cs_n$val[grep("padj_", cs_n$variable)] = "padj"
cs_n$val[grep("logpval_", cs_n$variable)] = "logpval"
cs_n$variable = NULL
cs_n = reshape2::dcast(cs_n, pathway + cancer_type + cell ~ val, value.var = "value")

### Oxidative Phosphorylation ###
cs_oxphos_can = cs_n[which(cs_n$pathway == "KEGG_OXIDATIVE_PHOSPHORYLATION"),]

cs_oxphos_can$cancer_type = factor(cs_oxphos_can$cancer_type, 
                            levels = cb$cancer_type[order(cb$cancer_type)])

oxphos = ggplot(cs_oxphos_can,aes(cancer_type, logpval, fill = cell)) +
  geom_bar(stat = "identity", width = 0.5,
           position = position_dodge(width=0.5)) +
  # scale_y_continuous(limits = c(-3, 3)) +
  theme_classic() +
  scale_fill_manual(values = c("#b2182b","#2166ac")) +
  labs(y = "sign(NES)*log10(pval)") +
  theme(legend.position="bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.box.spacing = unit(0.01, "cm"),
        legend.key.size = unit(0.4, "cm"),
        legend.key.width = unit(0.4,"cm"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(colour = "black", size = 7),
        axis.text = element_text(colour = "black", size = 7)) +
  coord_flip() + 
  theme(strip.text = element_text(size=7.5)) +
  geom_hline(yintercept = c(0,1.30103,-1.30103))

##############################################################################
###                           COMBINE ALL PLOTS                            ###
##############################################################################

grid.arrange(bars_all + theme(legend.position = "none"),
             arrangeGrob(bars_cs_ind,
                         oxphos + 
                           theme(axis.text.y = element_blank()), 
                         ncol = 2, widths = c(1,1.25),
                         top = textGrob("Oxphos",gp=gpar(fontsize=8))),
             top = textGrob("Cancer vs Stroma",gp=gpar(fontsize=10)))
