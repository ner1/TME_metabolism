library(preprocessCore)
library(reshape2)
library(data.table)
library(parallel)

##############################################
###           Read in TCGA data            ###
##############################################

# Log2(x+0.001) values from TCGA samples 
tcga_expr <- fread(file = "../data/TCGA/tcga_RSEM_gene_fpkm",
                   data.table = FALSE, 
                   sep = "\t")

###########################################################################
###                        Read GTEX Ovary normal                       ###
###########################################################################
# GTEX data downloaded from Xena

gtex_annotation <- read.delim(file = "../data/GTEX/GTEX_phenotype",
                              sep = "\t", stringsAsFactors = FALSE)

ov = gtex_annotation$Sample[which(gtex_annotation$body_site_detail..SMTSD. == "Ovary")]

gtex_expr = fread(file = "../data/GTEX/gtex_RSEM_gene_fpkm", stringsAsFactors = FALSE,
                  sep = "\t", data.table = FALSE)
gtex_expr = gtex_expr[,c(1,which(colnames(gtex_expr) %in% ov))]

### merge tcga with gtex
expr = merge(tcga_expr, gtex_expr, by = "sample", all.x = TRUE)

rm(gtex_annotation,gtex_expr,tcga_expr,ov)

################################################################
###                       Filter samples                    ###
###############################################################
# separate tumor and normal using link below
# https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes     
# 11 : Solid Tissue Normal  

normal_samples = colnames(expr)[grep(pattern = "-11$",x = colnames(expr))]
normal_samples = c(normal_samples,colnames(expr)[grep(pattern = "GTEX",x = colnames(expr))])

# load CNV data
cnv = fread(file = "../data/TCGA/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes",
            sep = "\t", data.table = FALSE,
            stringsAsFactors = FALSE)

# load and process tumor purity data
filesT = list.files(path = "../data/tumor_purity/",
                    pattern = "_purity.csv")
p=mclapply(filesT, function(f) {
  print(f)
  snv=read.table(paste("../data/from_UMESH/tumor_purity/",f, sep = ""),
                 header=T, sep=",", stringsAsFactors = FALSE)
  colnames(snv) = c("sample","cancer_purity")
  snv$type = gsub(pattern = "_purity.csv",replacement = "",f)
  return(snv)
})

purity = do.call(rbind,p)
rm(p,filesT)

### tumor samples with CNV and Purity data
tumor_samples = intersect(purity$sample, colnames(cnv))
samples_to_include = c(tumor_samples,normal_samples)

expr = expr[,c(1,which(colnames(expr) %in% samples_to_include))]

#############################################
###             Quality Control           ###
#############################################
# convert expression to linear and set negative values to zero
expr[,-1] = (2^expr[,-1]) - 0.001
expr[expr < 0] <- 0

# Exclude genes with zero expression in greater than 10% of individuals (normal and tumor both)
t <- as.data.frame(apply(expr[,-1], 1,
                         function(x) ifelse(sum(x <= 0) > 0.90*(ncol(expr)-1), 1, 0)))
colnames(t)[1] <- "freq"
t$gene <- expr$sample

print(table(t$freq>0))

expr <- expr[expr$sample %in% t$gene[t$freq==0], ]

rm(t)
gc()

# convert ENSG IDs to gene symbols
library(EnsDb.Hsapiens.v79)
expr$ensembl_gene_id = gsub("\\..*","",expr$sample)
EZ = select(EnsDb.Hsapiens.v79, key=expr$ensembl_gene_id,
            columns=c("GENEID", "SYMBOL"),
            keytype="GENEID")

colnames(EZ) = c("ensembl_gene_id","gene_symbol")
expr = merge(expr, EZ, by = "ensembl_gene_id", all.x = TRUE)
rm(EZ)
expr$gene_symbol[which(expr$ensembl_gene_id == "ENSG00000281991")] = "TMEM265"
expr$gene_symbol[which(expr$ensembl_gene_id == "ENSG00000282164")] = "PEG13"
expr$gene_symbol[which(expr$ensembl_gene_id == "ENSG00000282458")] = "WASH5P"
expr$gene_symbol[which(expr$ensembl_gene_id == "ENSG00000282608")] = "ADORA3"

expr = expr[-which(is.na(expr$gene_symbol)),]

expr$ensembl_gene_id = NULL
expr$sample = NULL

# add genes with more than one transcript
df_temp = setDT(expr)
df_sum <- df_temp[, lapply(.SD, sum), by = gene_symbol]

gc()

###########################################################################
###                        Read CCLE and match cancer types            ###
###########################################################################
## data downloaded from:

#https://portals.broadinstitute.org/ccle/data : filename = CCLE_RNAseq_genes_rpkm_20180929.gct.gz

ccle <- fread("../data/CCLE/CCLE_RNAseq_genes_rpkm_20180929.gct",
              stringsAsFactors = FALSE)
ccle$Name = NULL
colnames(ccle)[1] = "gene_symbol"

ccle_sum <- ccle[, lapply(.SD, sum), by = gene_symbol]

# match cancer types
ccle_annotation <- read.delim("../data/CCLE/Cell_lines_annotations_20181226.txt",
                              stringsAsFactors = FALSE)
ccle_annotation = ccle_annotation[,c("CCLE_ID","Site_Primary","Histology","tcga_code")]
ccle_annotation$tcga_code[which(ccle_annotation$tcga_code == "COAD/READ")] = "CRC"

write.csv(x = ccle_annotation, 
          file = "../data/CCLE/annotation_cancer_type.csv")

ccle_annotation = ccle_annotation[which(ccle_annotation$tcga_code %in% purity$type),]

write.csv(x = ccle_annotation, 
          file = "../data/CCLE/annotation_cancer_type.csv")

ccle_sum = ccle_sum[,.SD,.SDcols=c(1,which(colnames(ccle_sum) %in% ccle_annotation$CCLE_ID))]

# merge ccle with tumor
expr = merge(df_sum, ccle_sum, by = "gene_symbol")
rownames(expr) = expr$gene_symbol
expr$gene_symbol = NULL

# upper quartile normalisation
df1 = t(expr)
df2 = apply(df1, 1, function(x){quantile(x, 0.75)})
expr_norm <-as.data.frame(df1 / df2)
colnames(expr_norm) = rownames(expr)
expr_norm = expr_norm * mean(df2)

# log2(x+1) transform
expr_norm = log2(expr_norm+1)	

# clean
rm(df_sum,expr,df2)

# separate tumor, normal and CCLE using link below

ccle_norm = expr_norm[which(rownames(expr_norm) %in% colnames(ccle_sum)[-1]),]
normal_expr = expr_norm[which(rownames(expr_norm) %in% normal_samples),]
tumor_expr = expr_norm[which(rownames(expr_norm) %in% tumor_samples),]

rm(gtex, expr_norm)

fwrite(x = tumor_expr, file = "../data/CCLE/tumor_expr_norm.csv",quote = FALSE, row.names = TRUE)
fwrite(x = normal_expr, file = "../data/CCLE/normal_expr_norm.csv",quote = FALSE,row.names = TRUE)
fwrite(x = ccle_norm, file = "../data/CCLE/ccle_expr_norm.csv",quote = FALSE,row.names = TRUE)

#####################################################################
###                            TUMOR                              ###
#####################################################################
library(data.table)
tumor_expr = fread(file = "../data/CCLE/tumor_expr_norm.csv",
                   data.table = FALSE, sep = ",",
                   stringsAsFactors = FALSE)

tumor_expr = melt(tumor_expr)
colnames(tumor_expr) = c("sample","gene","tumor_expr")
tumor_expr.t = setDT(tumor_expr)

#####################################################################
###                   Read CNV data from Xena                     ###
#####################################################################
cnv = fread(file = "../data/TCGA/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes",
            sep = "\t", data.table = FALSE,
            stringsAsFactors = FALSE)

colnames(cnv)[1] = "gene"
cnv = melt(cnv)
colnames(cnv) = c("gene","sample","cnv")

cnv.t = setDT(cnv)
df_all = merge(tumor_expr.t,cnv.t, by = c("gene","sample"),all.x = TRUE	)

rm(cnv,cnv.t)
gc()

########################################################################
###                             Purity                               ###
########################################################################
# load and process tumor purity data
filesT = list.files(path = "../data/from_UMESH/tumor_purity/",
                    pattern = "_purity.csv")
p=mclapply(filesT, function(f) {
  print(f)
  snv=read.table(paste("../data/from_UMESH/tumor_purity/",f, sep = ""),
                 header=T, sep=",", stringsAsFactors = FALSE)
  colnames(snv) = c("sample","cancer_purity")
  snv$type = gsub(pattern = "_purity.csv",replacement = "",f)
  return(snv)
})

purity = do.call(rbind,p)
rm(p,filesT)

# merge purity with ct_expr
purity.t = setDT(purity)
tumor.table = merge(df_all,purity.t,by = "sample")

rm(df_all,purity,purity.t)

fwrite(x = tumor.table, file = "../data/CCLE/tumor_allD.csv",quote = FALSE)

########################################################################
###                           Deconvolution                          ###
########################################################################

library(parallel)
library(nnls)
library(data.table)
library(pbmcapply)

tumor_df = fread(file = "../data/CCLE/tumor_allD.csv", 
                 sep = ",", 
                 stringsAsFactor = FALSE, 
                 data.table = FALSE)

source("dcv.R")

tumor_df$stroma_purity = 1-tumor_df$cancer_purity

c_t = unique(tumor_df$type)

output <- vector("list",20)

for (i in 1:length(c_t)) {
  
  temp_ct = c_t[i]
  
  temp_df = tumor_df[which(tumor_df$type == temp_ct),]
  
 gene_list = unique(tumor_df$gene)
  
  p=pbmclapply(gene_list, function(g) {
    tmp = temp_df[which(temp_df$gene == g),]

    dcv_estimate = dcv(
      tmp = tmp,
      gene = g,
      expColName = "tumor_expr",
      cnv_data = TRUE
    )
    return(dcv_estimate)
  }
  , mc.cores = 8)
  
  dcv_tmp = do.call(rbind,p)
  dcv_tmp$cancer_type = temp_ct
  
  output[[i]] = dcv_tmp
  print(temp_ct)
fwrite(x = dcv_tmp, file = paste("../data/CCLE/dcv_",temp_ct,".csv"),sep = ",", quote = FALSE)
gc()
}
dcv_all = do.call(rbind,output)

fwrite(x = dcv_all, file = "../data/CCLE/dcv_all.csv",sep = ",", quote = FALSE)
save(dcv_all,file = "../data/CCLE/dcv_all.RData")

#####################################################################
###                            NORMAL                             ###
#####################################################################
library(data.table)
normal_expr = fread(file = "../data/CCLE/normal_expr_norm.csv",
                    data.table = FALSE, sep = ",",
                    stringsAsFactors = FALSE)

## cancer type mapping
colnames(normal_expr)[1] = "sample"
normal_expr = melt(normal_expr)
colnames(normal_expr) = c("sample","gene","normal_expr")

tss = read.csv(file = "../data/from_baruah/tss_tumor_mapping.csv",
               sep = ",",
               stringsAsFactors = FALSE)

tss$type[which(tss$type == "COAD")] = "CRC"
tss$type[which(tss$type == "READ")] = "CRC"

ct = c("BLCA","BRCA","CESC","CRC","ESCA","GBM","HNSC","KIRC","KIRP","LGG",
       "LIHC","LUAD","LUSC","OV","PAAD","PRAD","SKCM","STAD","THCA","UCEC")

tss = tss[which(tss$type %in% ct),]

normal_expr$tss_code = substr(normal_expr$sample, 6, 7)

normal_expr.t = setDT(normal_expr)
tss.t = setDT(tss)
aggregated.table = merge(normal_expr.t, tss.t, by = "tss_code", all.x = TRUE)
aggregated.table$type[grep(pattern="GTEX",x = aggregated.table$sample)] = "OV"
aggregated.table = aggregated.table[-which(is.na(aggregated.table$type)),]

normal_median = aggregated.table[, normal.m:=median(normal_expr), by=list(gene, type)]
normal_median$sample = NULL
normal_median$normal_expr = NULL
normal_median$tss_code = NULL
normal_median = unique(normal_median)

fwrite(x = normal_median, file = "../data/CCLE/normal_median.csv", sep = ",", quote = FALSE)

#####################################################################
###                              CCLE                             ###
#####################################################################
library(data.table)
ccle_expr = fread(file = "../data/CCLE/ccle_expr_norm.csv",
                    data.table = FALSE, sep = ",",
                    stringsAsFactors = FALSE)

## cancer type mapping
colnames(ccle_expr)[1] = "CCLE_ID"
ccle_expr = melt(ccle_expr)
colnames(ccle_expr) = c("CCLE_ID","gene","ccle_expr")

ccle_annotation <- read.delim(file = "../data/CCLE/annotation_cancer_type.csv",
                              sep = ",",
                              stringsAsFactor = FALSE)

ccle_expr = merge(ccle_expr,ccle_annotation)
ccle_expr = ccle_expr[,c("CCLE_ID","tcga_code","gene","ccle_expr")]
colnames(ccle_expr)[2] = "cancer_type"

ccle_df = data.table(ccle_expr)

ccle_median = ccle_df[, ccle.m:=median(ccle_expr), by=list(gene, cancer_type)]
ccle_median$CCLE_ID = NULL
ccle_median$ccle_expr = NULL
ccle_median = unique(ccle_median)

fwrite(x = ccle_median, file = "../data/CCLE/ccle_median.csv", 
       sep = ",", quote = FALSE)
