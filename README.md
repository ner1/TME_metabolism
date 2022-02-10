# A pan-cancer metabolic atlas of the tumor microenvironment
All data was downloaded from Xena UCSC: 
1. TCGA: Gene expression (tcga_RSEM_gene_fpkm) and copy number variation (Gistic2_CopyNumber_Gistic2_all_thresholded) data was downloaded from:  https://xenabrowser.net/datapages/?cohort=TCGA%20Pan-Cancer%20(PANCAN)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
2. GTEX: Gene expression (gtex_RSEM_gene_fpkm) and phenotype (GTEX_phenotype) information was downloaded from: https://xenabrowser.net/datapages/?cohort=GTEX&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
3. CCLE: Gene expression (CCLE_RNAseq_genes_rpkm_20180929) and phenotype data (Cell_lines_annotations_20181226) was downloaded from: https://xenabrowser.net/datapages/?cohort=Cancer%20Cell%20Line%20Encyclopedia%20(CCLE)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443

Folder “code” contains all the functions and scripts that were used for the analysis. Description of the steps taken in the analysis are detailed below:
1. Processing of the data for analysis was performed using the R script : “dataPreprocessing.R”
2. Gene symbols were converted to Entrez IDs for further analysis using RECON3 (R script : “convert_symbol_to_entrez.R”).
3. Gene expression analysis to highlight cancer and stromal cells specific genes (Figure 2 and 3) was performed using the R script “gene_expression_cs.R”
4. Gene set enrichment analysis to compare cancer and stromal pathways (Figure 4) was performed using the R script “GSEA_cs.R”
5. Gene expression analysis to highlight cancer cell lines specific genes (Figure 5) was performed using the R script “gene_expression_ccle.R”
6. Gene set enrichment analysis was performed to compare metabolic pathways  in cancer cells in vivo vs in vitro (Figure 5) using the R script “GSEA_ccle.R”
7. Reaction activity for cancer cells, stromal cells and cancer cell lines was estimated using constraint based metabolic model RECON3 (Figure 4 and 5). See methods section “Metabolic map generation” for details. Matlab script (“EstimateRxnActivity.m”) and R script (“Escher_map_input.R”) were used for the metabolic map generation. 
8. Validation of top genes predictions using publically available Single Cell RNAseq dataset for SKCM and HNSC (R Script : "single_Cell_validation_genes.R")
9. Validation of gene set enrichment results using publically available Single Cell RNAseq dataset for SKCM and HNSC (R Script : "single_Cell_validation.R")
10. TME_metabolism.Rmd is the R markdown to run the analysis and generate paper figures. 

Folder “data” contains the data used in generation of the figures in the paper:
1. Expression_Model_All.dat: Metabolic gene expression for cancer cells, stromal cells, normal cells and cancer cell lines estimated after deconvolution.
2. tumor_purity: Purity values for tumor samples estimated using TUMERIC (Ghoshdastider et al. 2021).
3. kegg_oxphos_genes_ab.txt: ALIASES for genes in KEGG oxidative phosphorylation pathway
4. GSEA_KEGG: Results from gene set enrichment analysis performed using R package “fgsea”.
 4.1. fgseaRes_C_CCLE.csv: KEGG gene sets enriched in cancer cells as compared to cancer cell lines for all the tumor types.
 4.2. fgseaRes_C_CCLE_ind.csv: KEGG gene sets enriched in cancer cells as compared to cancer cell lines in each tumor types.
 4.3. fgseaRes_CCLE_N.csv: KEGG gene sets enriched in cancer cell lines as compared to normal tissue expression for each tumor types.
 4.4. fgseaRes_CN_ind.csv: KEGG gene sets enriched in cancer cells as compared to normal tissue for each tumor type.
 4.5. fgseaRes_CS_all.csv: KEGG gene sets enriched in cancer cells as compared to stromal cells within tumor for all the tumor types.
 4.6. fgseaRes_CS_ind.csv: KGG gene sets enriched in cancer cells as compared to stromal cells within tumor for each of the tumor types.
 4.7. fgseaRes_SN_ind.csv: KEGG gene sets enriched in stromal cells as compared to the normal tissue for each tumor type. 
