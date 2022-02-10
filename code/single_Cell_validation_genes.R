library(fgsea)
library(KEGGREST)
library(ggplot2)
library(gridExtra)
library(stringr)
library(tidytext)
library(ggpubr)

ExprMetabolic = read.table(file = '../data/Expression_Model_All.dat',
                           sep = ',', header = TRUE, stringsAsFactors = FALSE)

ExprMetabolic$ReconGeneName = NULL
ExprMetabolic = ExprMetabolic[!duplicated(ExprMetabolic),]

################# SKCM
skcm_my = ExprMetabolic[which(ExprMetabolic$cancer_type == "SKCM"),]
skcm_my$cs = skcm_my$cancer - skcm_my$stroma
stroma_skcm = as.character(skcm_my$gene[order(skcm_my$cs, decreasing = FALSE)[1:15]])
cancer_skcm = as.character(skcm_my$gene[order(skcm_my$cs, decreasing = TRUE)[1:15]])

### ################# Sigle Cell SKCM expression data (Tirosh et al. 2016)
# skcm <- read.delim(file = "C:/AndersLab/TME_metabolism/data/GSE72056_melanoma_single_cell_revised_v2.txt.gz",
#                    sep = "\t", 
#                    stringsAsFactors = FALSE)
# cancer_cols = which(skcm[2,] == 2)
# stroma_cols = which(skcm[2,] == 1)
# colnames(skcm)[cancer_cols] = paste(colnames(skcm)[cancer_cols],
#                                     "_Cancer",sep = "")
# colnames(skcm)[stroma_cols] = paste(colnames(skcm)[stroma_cols],
#                                     "_Stroma",sep = "")
# skcm = skcm[,c(1,cancer_cols,stroma_cols)]
# skcm = skcm[-c(1:3),]

########################################
### cancer

# expr_cancer_skcm = skcm[which(skcm$Cell %in% cancer_skcm),]
# save(expr_cancer_skcm, file = "../data/validation/expr_cancer_skcm.R")

load("../data/validation/expr_cancer_skcm.R")
expr_cancer_skcm = reshape2::melt(expr_cancer_skcm, id.vars = "Cell", 
                             variable.name = "sample", value.name = "expr")
expr_cancer_skcm$single_cell_sample = NA
expr_cancer_skcm$single_cell_sample[grep(pattern = "Cancer", 
                                    x = expr_cancer_skcm$sample)] = "Cancer"
expr_cancer_skcm$single_cell_sample[grep(pattern = "Stroma", 
                                    x = expr_cancer_skcm$sample)] = "Stroma"
colnames(expr_cancer_skcm)[1] = "Gene"
result_s <- aggregate(expr ~ Gene, 
                      expr_cancer_skcm[which(expr_cancer_skcm$single_cell_sample == "Cancer"),], 
                      function(x) quantile(x, probs = 0.75))

result_s <- result_s[order(result_s$expr, decreasing = TRUE),]

cancer_level = factor(result_s$Gene,levels=result_s$Gene)
expr_cancer_skcm$Gene = factor(expr_cancer_skcm$Gene, 
                          levels = cancer_level)

box_cancer_skcm = ggplot(expr_cancer_skcm, 
       aes(x = Gene, y = expr, fill = single_cell_sample)) + 
  geom_boxplot(outlier.shape = NA) + 
  theme_classic() + 
  scale_fill_manual(values = c("#b2182b","#2166ac")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  ggtitle("SKCM-CancerGenes") 

#################################################
### stroma
# expr_stroma_skcm = skcm[which(skcm$Cell %in% stroma_skcm),]
# save(expr_stroma_skcm, file = "../data/validation/expr_stroma_skcm.R")

load(file = "../data/validation/expr_stroma_skcm.R")

expr_stroma_skcm = reshape2::melt(expr_stroma_skcm, id.vars = "Cell", 
                             variable.name = "sample", value.name = "expr")
expr_stroma_skcm$single_cell_sample = NA
expr_stroma_skcm$single_cell_sample[grep(pattern = "Cancer", 
                                    x = expr_stroma_skcm$sample)] = "Cancer"
expr_stroma_skcm$single_cell_sample[grep(pattern = "Stroma", 
                                    x = expr_stroma_skcm$sample)] = "Stroma"
colnames(expr_stroma_skcm)[1] = "Gene"
result_s <- aggregate(expr ~ Gene, 
                      expr_stroma_skcm[which(expr_stroma_skcm$single_cell_sample == "Stroma"),], 
                      function(x) quantile(x, probs = 0.75))

result_s <- result_s[order(result_s$expr, decreasing = TRUE),]

stroma_level = factor(result_s$Gene,levels=result_s$Gene)
expr_stroma_skcm$Gene = factor(expr_stroma_skcm$Gene, 
                               levels = stroma_level)

box_stroma_skcm = ggplot(expr_stroma_skcm, 
       aes(x = Gene, y = expr, fill = single_cell_sample)) + 
  geom_boxplot(outlier.shape = NA) + 
  theme_classic() + 
  scale_fill_manual(values = c("#b2182b","#2166ac")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  ggtitle("SKCM-StromaGenes") 

################# HNSC

hnsc_my = ExprMetabolic[which(ExprMetabolic$cancer_type == "HNSC"),]
hnsc_my$cs = hnsc_my$cancer - hnsc_my$stroma
stroma_hnsc = as.character(hnsc_my$gene[order(hnsc_my$cs, decreasing = FALSE)[1:15]])
cancer_hnsc = as.character(hnsc_my$gene[order(hnsc_my$cs, decreasing = TRUE)[1:15]])

################# Sigle Cell HNSC expression data (Puram et al. 2017)
# hnsc <- read.delim(file = "C:/AndersLab/TME_metabolism/data/GSE103322_HNSCC_all_data.txt.gz",
                   # sep = "\t",stringsAsFactors = FALSE)
# cancer_cols = which(hnsc[3,] == 1)
# stroma_cols = which(hnsc[4,] == 1)
# colnames(hnsc)[cancer_cols] = paste(colnames(hnsc)[cancer_cols],
#                                     "_Cancer",sep = "")
# colnames(hnsc)[stroma_cols] = paste(colnames(hnsc)[stroma_cols],
#                                     "_Stroma",sep = "")
# hnsc = hnsc[,c(1,cancer_cols,stroma_cols)]
# hnsc = hnsc[-c(1:5),]
# hnsc$X = gsub(pattern = "[']",x = hnsc$X,replacement = "")

### cancer
# expr_cancer_hnsc = hnsc[which(hnsc$X %in% cancer_hnsc),]
# save(expr_cancer_hnsc, file = "../data/validation/expr_cancer_hnsc.R")

load(file = "../data/validation/expr_cancer_hnsc.R")

expr_cancer_hnsc = reshape2::melt(expr_cancer_hnsc, id.vars = "X", 
                             variable.name = "sample", value.name = "expr")
expr_cancer_hnsc$single_cell_sample = NA
expr_cancer_hnsc$single_cell_sample[grep(pattern = "Cancer", 
                                    x = expr_cancer_hnsc$sample)] = "Cancer"
expr_cancer_hnsc$single_cell_sample[grep(pattern = "Stroma", 
                                    x = expr_cancer_hnsc$sample)] = "Stroma"
colnames(expr_cancer_hnsc)[1] = "Gene"
expr_cancer_hnsc$expr = as.numeric(as.character(expr_cancer_hnsc$expr))

result_h <- aggregate(expr ~ Gene, 
                      expr_cancer_hnsc[which(expr_cancer_hnsc$single_cell_sample == "Cancer"),], 
                      function(x) quantile(x, probs = 0.75))

result_h <- result_h[order(result_h$expr, decreasing = TRUE),]

cancer_level = factor(result_h$Gene,levels=result_h$Gene)
expr_cancer_hnsc$Gene = factor(expr_cancer_hnsc$Gene, 
                               levels = cancer_level)

box_cancer_hnsc = ggplot(expr_cancer_hnsc, 
       aes(x = Gene, y = expr, fill = single_cell_sample)) + 
  geom_boxplot(outlier.shape = NA) +
  theme_classic() + 
  scale_fill_manual(values = c("#b2182b","#2166ac")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  ggtitle("HNSC-CancerGenes") 

### stroma
# expr_stroma_hnsc = hnsc[which(hnsc$X %in% stroma_hnsc),]
# save(expr_stroma_hnsc, file = "../data/validation/expr_stroma_hnsc.R")

load(file = "../data/validation/expr_stroma_hnsc.R")

expr_stroma_hnsc = reshape2::melt(expr_stroma_hnsc, id.vars = "X", 
                             variable.name = "sample", value.name = "expr")
expr_stroma_hnsc$single_cell_sample = NA
expr_stroma_hnsc$single_cell_sample[grep(pattern = "Cancer", 
                                    x = expr_stroma_hnsc$sample)] = "Cancer"
expr_stroma_hnsc$single_cell_sample[grep(pattern = "Stroma", 
                                    x = expr_stroma_hnsc$sample)] = "Stroma"

colnames(expr_stroma_hnsc)[1] = "Gene"
expr_stroma_hnsc$expr = as.numeric(as.character(expr_stroma_hnsc$expr))

result_h <- aggregate(expr ~ Gene, 
                      expr_stroma_hnsc[which(expr_stroma_hnsc$single_cell_sample == "Stroma"),], 
                      function(x) quantile(x, probs = 0.75))

result_h <- result_h[order(result_h$expr, decreasing = TRUE),]

stroma_level = factor(result_h$Gene,levels=result_h$Gene)
expr_stroma_hnsc$Gene = factor(expr_stroma_hnsc$Gene, 
                               levels = stroma_level)

box_stroma_hnsc = ggplot(expr_stroma_hnsc, 
       aes(x = Gene,fill = single_cell_sample, y = expr)) + 
  geom_boxplot(outlier.shape = NA) +
  theme_classic() + 
  scale_fill_manual(values = c("#b2182b","#2166ac")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  ggtitle("HNSC-StromaGenes") 


grid.arrange(box_cancer_skcm, box_stroma_skcm, box_cancer_hnsc, box_stroma_hnsc)