library(biomaRt)
library(org.Hs.eg.db)

mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

## read in data
# tumor
exprAll = read.delim('../data/CCLE/dcv_all.csv', sep = ",", header = TRUE, stringsAsFactors = FALSE)
colnames(exprAll)[1:2] = c("stroma","cancer")

#normal
exprNormal = read.csv("../data/CCLE/normal_median.csv", 
                      header = TRUE, stringsAsFactors = FALSE)
colnames(exprNormal) = c("gene","cancer_type","normal")

lgg = exprNormal[which(exprNormal$cancer_type == "GBM"),]
lgg$cancer_type = "LGG"
exprNormal = rbind(exprNormal,lgg)

# merge cancer and normal
exprAll = merge(exprAll,exprNormal)

rm(exprNormal,lgg)
gc()

# CCLE
exprAll_CCLE = read.delim("../data/CCLE/ccle_median.csv", 
                          header = TRUE, sep = ",", stringsAsFactors = FALSE)

## merge CCLE, cancer and normal for comparison
exprAll = merge(exprAll,exprAll_CCLE,all.x = TRUE)
exprAll[is.na(exprAll)] = 0

rm(exprAll_CCLE)
gc()

## Convert gene symbols for cancer, tumor, normal and stroma analysis
EZ = select(org.Hs.eg.db, unique(exprAll$gene), c("ENTREZID"), "SYMBOL")
x = EZ$SYMBOL[which(is.na(EZ$ENTREZID))]
colnames(EZ) = c("gene", "entrezgene1")

EZ3 = select(org.Hs.eg.db, unique(x), c("ENTREZID"), "ALIAS")
colnames(EZ3) = c("gene", "entrezgene2")

EZ = merge(EZ,EZ3,all = TRUE)
EZ$entrezgene = NA

EZ$entrezgene = ifelse(is.na(EZ$entrezgene1), EZ$entrezgene2,EZ$entrezgene1)
EZ = EZ[,c("gene","entrezgene")]
x = EZ$gene[which(is.na(EZ$entrezgene))]

EZ2 = getBM(attributes = c("entrezgene_id","hgnc_symbol"),
            filters = "hgnc_symbol",
            values = x,
            mart = mart)

for(i in 1:length(x)){
  if (length(which(EZ2$hgnc_symbol == x[i])) > 0) {
    EZ$entrezgene[which(EZ$gene == x[i])] = EZ2$entrezgene_id[which(EZ2$hgnc_symbol == x[i])]
  }
}

## mart has it wrong
EZ$entrezgene[which(EZ$gene == "MT-CO2")] = 4513
EZ = EZ[which(!is.na(EZ$entrezgene)),]

exprAll_entrez = merge(exprAll,EZ, by = "gene", all = TRUE)
exprAll_entrez = exprAll_entrez[which(!is.na(exprAll_entrez$entrezgene)),]

# write.table(exprAll_entrez, '../data/CCLE/ExprAll_Entrez.txt', sep = '\t', row.names = FALSE)
