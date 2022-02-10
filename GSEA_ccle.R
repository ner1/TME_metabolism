library(reshape2)
library(cowplot)
library(ggrepel)
library(ggtern)
library(grid)
library(gridExtra)
library(gtable)
library(KEGGREST)
library(fgsea)

##################################################################################
###                                    Data                                    ###
##################################################################################

# my_pathways <- gmtPathways("../data/c2.cp.kegg.v6.2.symbols.gmt")

# extract a list of metabolic pathways from KEGG #############
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

#############################################################################################

Expr = read.table(file = '../data/Expression_Model_All.dat',
                  sep = ',', header = TRUE, stringsAsFactors = FALSE)

Expr = Expr[-which(Expr$cancer_type %in% c("KIRP")),]
Expr$ReconGeneName = NULL
Expr = Expr[!duplicated(Expr),]

kegg_ab = read.csv(file = "../data/kegg_oxphos_genes_ab.txt", 
                   header = FALSE, stringsAsFactors = FALSE)

colnames(kegg_ab) = c("kegg_gene","gene")
Expr = merge(Expr,kegg_ab, by = "gene", all.x = TRUE)
Expr$kegg_gene[which(is.na(Expr$kegg_gene))] = Expr$gene[which(is.na(Expr$kegg_gene))]
Expr$gene = Expr$kegg_gene

C = Expr[,c("gene","cancer_type","cancer")]
C = reshape2::dcast(C, gene ~ cancer_type, value.var = "cancer")
C$Cm = apply(C[,-1], 1, FUN = function(x) median(x))

CCLE = Expr[,c("gene","ccle_m","cancer_type")]
CCLE = reshape2::dcast(CCLE, gene ~ cancer_type, value.var = "ccle_m")
CCLE$cclem = apply(CCLE[,-1], 1, FUN = function(x) median(x))

N = Expr[,c("gene","normal","cancer_type")]
N = reshape2::dcast(N, gene ~ cancer_type, value.var = "normal")
N$Nm = apply(N[,-1], 1, FUN = function(x) median(x))

combinedMedian = merge(N,C,by = "gene")
combinedMedian = combinedMedian[,c("gene","Nm","Cm")]
combinedMedian = merge(combinedMedian, CCLE, by = "gene" )
combinedMedian = combinedMedian[,c("gene","Nm","Cm","cclem")]
expr = combinedMedian

## Cancer - CCLE
expr$C_CL = expr$Cm-expr$cclem
gene_rank = expr[order(expr$C_CL, decreasing = TRUE),]
g_rank = as.numeric(gene_rank$C_CL)
names(g_rank) = gene_rank$gene

# fgseaRes_C_CCLE <- fgsea(pathways = my_pathways,
#                          stats = g_rank,
#                          minSize=15,
#                          maxSize=500,
#                          nperm=100000)

# save(fgseaRes_C_CCLE, file = "../data/CCLE/GSEA/KEGG_all_ct/fgsea_C_CCLE_pancan.RData")
# fgseaRes_C_CCLE$leadingEdge = NULL
# write.csv(x = fgseaRes_C_CCLE, file = "../data/CCLE/GSEA/for_web_app/fgseaRes_C_CCLE_all.csv") #for webapp

fgseaRes_C_CCLE = read.csv(file = "../data/GSEA_KEGG/fgseaRes_C_CCLE.csv", stringsAsFactors = FALSE)

fgseaRes_C_CCLE$logpadj = -log10(fgseaRes_C_CCLE$padj)* sign(fgseaRes_C_CCLE$NES)

c_ccle_kegg = fgseaRes_C_CCLE[grep(pattern = "KEGG_", x = fgseaRes_C_CCLE$pathway),]
c_ccle_kegg = c_ccle_kegg[which(c_ccle_kegg$pathway %in% metabolic_pathways$path_name),]
c_ccle_kegg = c_ccle_kegg[which(abs(c_ccle_kegg$NES) > 1.4),]

c_ccle_kegg$pathway = factor(c_ccle_kegg$pathway, 
                             levels = c_ccle_kegg$pathway[order(c_ccle_kegg$logpadj)])

c_ccle_kegg$star = ""
c_ccle_kegg$star[which(c_ccle_kegg$padj <= 0.05 )] = "*"

bars_all = ggplot(data = c_ccle_kegg, aes(x = pathway, y= logpadj)) + 
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
  ggtitle("Cancer vs Cell Lines (PANCAN)") + 
  geom_hline(yintercept = c(-1.30103,0,1.30103)) + 
  coord_flip() + 
  geom_text(aes(label=star), colour="black", hjust = 0, vjust=0.7, size=6)+ ylim(c(-5,5)) 

# enrich_ccle_pan = plotEnrichment(my_pathways[["KEGG_OXIDATIVE_PHOSPHORYLATION"]],
#                                g_rank) + labs(title="Oxidative Phosphorylation") + 
#   theme_classic() + 
#   theme(axis.title = element_text(colour = "black", size = 7),
#         axis.text = element_text(colour = "black", size = 7),
#         title = element_text(colour = "black", size = 7))
# 

# plotGseaTable(pathways = my_pathways["KEGG_OXIDATIVE_PHOSPHORYLATION"],
#               fgseaRes = fgseaRes_C_CCLE,stats = g_rank, gseaParam = 0.8)

##########################################################################
###                          individual heatmaps                       ###
##########################################################################

cancer_type = unique(Expr$cancer_type)

# CANCER - ccle
Expr$C_ccle = Expr$cancer - Expr$ccle

# fgseaRes_c_ccle_ind = vector(length = 19, "list")
# for (ct in 1:length(cancer_type)) {
#   tmp = Expr[which(Expr$cancer_type == cancer_type[[ct]]),]
#   gene_rank = tmp[order(tmp$C_ccle, decreasing = TRUE),]
#   g_rank = as.numeric(gene_rank$C_ccle)
#   names(g_rank) = gene_rank$gene
# 
#   temp_res <- fgsea(pathways = my_pathways,
#                     stats = g_rank,
#                     minSize=15,
#                     maxSize=500,
#                     nperm=100000)
# 
#   temp_res$cancer_type = cancer_type[[ct]]
#   fgseaRes_c_ccle_ind[[ct]] <- temp_res
#   print(ct)
# }
# save(fgseaRes_c_ccle_ind, file = "../data/CCLE/GSEA/KEGG_all_ct/fgseaRes_c_ccle_ind.RData")
# c_ccle_ind = do.call(rbind,fgseaRes_c_ccle_ind)
# c_ccle_ind$leadingEdge = NULL
# write.csv(c_ccle_ind, file = "../data/GSEA_KEGG/fgseaRes_C_CCLE_ind.csv")

c_ccle_ind = read.csv(file = "../data/GSEA_KEGG/fgseaRes_C_CCLE_ind.csv",
                      stringsAsFactors = FALSE)

c_ccle_ind_oxphos = c_ccle_ind[which(c_ccle_ind$pathway == "KEGG_OXIDATIVE_PHOSPHORYLATION"),]

c_ccle_ind_oxphos$logpval = -log10(c_ccle_ind_oxphos$pval) * sign(c_ccle_ind_oxphos$NES)
c_ccle_ind_oxphos$cancer_type = factor(c_ccle_ind_oxphos$cancer_type, 
                                       levels = c_ccle_ind_oxphos$cancer_type[order(c_ccle_ind_oxphos$logpval, decreasing = TRUE)])


c_ccle_ind_oxphos$star = ""
c_ccle_ind_oxphos$star[which(c_ccle_ind_oxphos$pval <= 0.05 )] = "*"
c_ccle_ind_oxphos$leg = "C vs CL"

bars_c_ccle_oxphos = ggplot(data = c_ccle_ind_oxphos, aes(x = cancer_type, y= logpval, fill = leg)) + 
  geom_bar(stat = "identity") + 
  theme_classic() +
  scale_fill_manual(values = c("#A9A9A9")) +
  labs(y = "sign(NES)*log10(pval)") +
  theme(legend.position="right",
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.box.spacing = unit(0.01, "cm"),
        legend.key.size = unit(0.4, "cm"),
        legend.key.width = unit(0.4,"cm"),
        axis.title.y = element_text(colour = "black", size = 7),
        axis.title.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 7),
        axis.text.x = element_blank()) + 
  geom_text(aes(label=star), colour="black", hjust = 0, vjust=0.7, size=6) + 
  ylim(c(-5,5)) +
  geom_hline(yintercept = c(-1.30103,0,1.30103))

# cancer
cancer = read.csv(file = "../data/GSEA_KEGG/fgseaRes_CN_ind.csv", stringsAsFactors = FALSE)
cancer$logpval = -log10(cancer$pval) * sign(cancer$NES)
colnames(cancer)[-c(1,8)] <- paste(colnames(cancer)[-c(1,8)],"cancer", sep = "_")

# ccle 
Expr$CL_N = Expr$ccle - Expr$normal

# fgseaRes_ccle = vector(length = 19, "list")
# for (ct in 1:length(cancer_type)) {
#   tmp = Expr[which(Expr$cancer_type == cancer_type[[ct]]),]
#   gene_rank = tmp[order(tmp$CL_N, decreasing = TRUE),]
#   g_rank = as.numeric(gene_rank$CL_N)
#   names(g_rank) = gene_rank$gene
# 
#   temp_res <- fgsea(pathways = my_pathways,
#                     stats = g_rank,
#                     minSize=15,
#                     maxSize=500,
#                     nperm=100000)
# 
#   temp_res$cancer_type = cancer_type[[ct]]
#   fgseaRes_ccle[[ct]] <- temp_res
#   print(ct)
# }
# save(fgseaRes_ccle, file = "../data/CCLE/GSEA/KEGG_all_ct/fgseaRes_ccle_ind.RData")
# 
# ccle = do.call(rbind,fgseaRes_ccle)
# ccle$leadingEdge = NULL
# write.csv(x = ccle, file = "../data/GSEA_KEGG/fgseaRes_CCLE_N.csv")

ccle = read.csv(file = "../data/GSEA_KEGG/fgseaRes_CCLE_N.csv",
                stringsAsFactors = FALSE)
ccle$X = NULL
ccle$logpval = -log10(ccle$pval) * sign(ccle$NES)
colnames(ccle)[-c(1,8)] <- paste(colnames(ccle)[-c(1,8)],"ccle", sep = "_")
rownames(ccle) = NULL

# merge both and compare 
c_cl_n = merge(cancer,ccle)

c_cl_n = c_cl_n[,c("pathway","cancer_type",
                   "NES_ccle","NES_cancer",
                   "pval_ccle","pval_cancer",
                   "logpval_cancer","logpval_ccle")]
c_cl_n = reshape2::melt(c_cl_n)
c_cl_n$cell = "cell"
c_cl_n$cell[grep("_cancer", c_cl_n$variable)] = "C vs N"
c_cl_n$cell[grep("_ccle", c_cl_n$variable)] = "CL vs N"
c_cl_n$val = "v"
c_cl_n$val[grep("NES_", c_cl_n$variable)] = "NES"
c_cl_n$val[grep("pval_", c_cl_n$variable)] = "pval"
c_cl_n$val[grep("logpval_", c_cl_n$variable)] = "logpval"
c_cl_n$variable = NULL
c_cl_n = reshape2::dcast(c_cl_n, pathway + cancer_type + cell ~ val, value.var = "value")

c_cl_n = c_cl_n[grep(pattern = "KEGG_",x = c_cl_n$pathway),]
c_cl_n = c_cl_n[which(c_cl_n$pathway %in% metabolic_pathways$path_name),]

### Oxidative Phosphorylation ###
c_ccle_oxphos_can = c_cl_n[which(c_cl_n$pathway  == "KEGG_OXIDATIVE_PHOSPHORYLATION"),]

c_ccle_oxphos_can$cancer_type = factor(c_ccle_oxphos_can$cancer_type,
                                       levels = c_ccle_ind_oxphos$cancer_type[order(c_ccle_ind_oxphos$logpval, decreasing = TRUE)])


oxphos_can = ggplot(data = c_ccle_oxphos_can, aes(x = cancer_type, y= logpval, fill = cell)) + 
  geom_bar(stat = "identity",width = 0.5,
           position = position_dodge(width=0.5)) + 
  scale_fill_manual(values = c("#b2182b","#984ea3")) + 
  theme_classic() +
  labs(y = "sign(NES)*log10(pval)") +
  theme(legend.position="right",
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.box.spacing = unit(0.01, "cm"),
        legend.key.size = unit(0.4, "cm"),
        legend.key.width = unit(0.4,"cm"),
        axis.title.y = element_text(colour = "black", size = 7),
        axis.title.x = element_text(colour = "black", size = 7),
        axis.text.y = element_text(colour = "black", size = 7),
        axis.text.x = element_text(colour = "black", size = 7, angle = 45, hjust = 1, vjust = 1)) +
  ylim(c(-5,5)) + 
  geom_hline(yintercept = c(-1.30103,0,1.30103))

##############################################################################
###                           COMBINE ALL PLOTS                            ###
##############################################################################

grid.arrange(bars_all + theme(legend.position = "none"),
             arrangeGrob(bars_c_ccle_oxphos,
                         oxphos_can, ncol = 1, heights = c(1,1.25),
                         top = textGrob("Oxphos",gp=gpar(fontsize=8))),
             ncol = 1, 
             top = textGrob("Cancer Cells vs Cell lines",gp=gpar(fontsize=10)))
