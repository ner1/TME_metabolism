library(reshape2)
library(cowplot)
library(ggrepel)
library(ggtern)
library(grid)
library(gridExtra)
library(gtable)
library(KEGGREST)
library(fgsea)

ccle_norm = read.delim(file = "../data/Expression_Model_All.dat",
                       stringsAsFactors = FALSE, header = TRUE, sep = ",")

ccle_norm$ReconGeneName = NULL
ccle_norm = ccle_norm[!duplicated(ccle_norm),]

# no KIRP in cell line data
ccle_norm = ccle_norm[-which(ccle_norm$cancer_type %in% c("KIRP")),]


######### Gene Level Analysis #########
ccle_norm$cancer.raw = (2^ccle_norm$cancer)-1
ccle_norm$ccle.raw = (2^ccle_norm$ccle_m)-1
ccle_norm$normal.raw = (2^ccle_norm$normal)-1

C = ccle_norm[,c("gene","cancer_type","cancer.raw")]
C = reshape2::dcast(C, gene ~ cancer_type, value.var = "cancer.raw")
C$Cm = apply(C[,-1], 1, FUN = function(x) median(x))

ccle = ccle_norm[,c("gene","cancer_type","ccle.raw")]
ccle = reshape2::dcast(ccle, gene ~ cancer_type, value.var = "ccle.raw")
ccle$ccle_m = apply(ccle[,-1], 1, FUN = function(x) median(x))

N = ccle_norm[,c("gene","normal.raw","cancer_type")]
N = reshape2::dcast(N, gene ~ cancer_type, value.var = "normal.raw")
N$Nm = apply(N[,-1], 1, FUN = function(x) median(x))

combinedMedian = merge(N,C,by = "gene")
combinedMedian = combinedMedian[,c("gene","Nm","Cm")]
combinedMedian = merge(combinedMedian, ccle, by = "gene" )
combinedMedian = combinedMedian[,c("gene","Nm","Cm","ccle_m")]

combinedMedian = combinedMedian[which(combinedMedian$Nm  > 1|
                                        combinedMedian$Cm > 1|
                                        combinedMedian$ccle_m > 1),]

########################################  ternary plot  ###########################################
ternPlot = combinedMedian

tern.colors = list()
tern.colors[""] = "#000000"
tern.colors["cancer"] = "#b2182b"
tern.colors["ccle"] = "#984ea3"
tern.colors = unlist(tern.colors)

combinedMedian$sum = combinedMedian$Nm + combinedMedian$Cm + combinedMedian$ccle_m
combinedMedian$normal.p = (combinedMedian$Nm/combinedMedian$sum)*100
combinedMedian$cancer.p = (combinedMedian$Cm/combinedMedian$sum)*100
combinedMedian$ccle.p = (combinedMedian$ccle_m/combinedMedian$sum)*100

ccle_genes = as.character(combinedMedian$gene[order(combinedMedian$ccle.p, 
                                                    decreasing = TRUE)[1:15]])

ternPlot$cl = ""
ternPlot$cl[which(ternPlot$gene %in% ccle_genes)] = "ccle"

ternPlot$al = 0.9
ternPlot$al[which(ternPlot$gene %in% ccle_genes)] = 1

# For Rendering the Lines, use Segment geometry
lines <- data.frame(x = c(0.5, 0, 0.5), 
                    y = c(0.5, 0.5, 0), 
                    z = c(0, 0.5, 0.5), 
                    xend = c(1, 1, 1)/3, 
                    yend = c(1, 1, 1)/3, 
                    zend = c(1, 1, 1)/3)

a1 = ggtern(data = ternPlot, aes(x = Cm, y = Nm, z = ccle_m)) + 
  geom_point(size=2.5,aes(colour=cl, alpha = al)) + Tlab("N") +
  Llab("C") + Rlab("CCLE") +
  limit_tern(1.1,1.1,1.1) + 
  scale_color_manual(values = tern.colors)  + theme_bw() + 
  theme(axis.title = element_text(colour = "black", size = 10),
        axis.text = element_text(colour = "black", size = 10),
        legend.position = "none") + 
  geom_segment(data = lines, 
               aes(x, y, z, xend = xend, yend = yend, zend = zend), 
               color = 'red', size = 1)

########################################################################################
####              HEAT MAP                                          ###
########################################################################################

ccle_norm$sum = ccle_norm$cancer.raw + ccle_norm$ccle.raw + ccle_norm$normal.raw
ccle_norm$cancer.p = (ccle_norm$cancer.raw/ccle_norm$sum)*100
ccle_norm$ccle.p = (ccle_norm$ccle.raw/ccle_norm$sum)*100
ccle_norm$normal.p = (ccle_norm$normal.raw/ccle_norm$sum)*100

#ccle
ccle_proportion = ccle_norm[,c("gene","ccle.p","cancer_type")]
ccle_proportion = reshape2::dcast(ccle_proportion, gene ~ cancer_type, value.var = "ccle.p")
ccle_proportion[is.na(ccle_proportion)] <- 0
ccle_proportion$ccle_proportion_m = apply(ccle_proportion[,-1], 1, FUN = function(x) median(x))

ccle.p_ind = ccle_proportion[which(ccle_proportion$gene %in% ccle_genes),]
ccle.p_ind$gene = factor(ccle.p_ind$gene, 
                      levels = ccle.p_ind$gene[order(ccle.p_ind$ccle_proportion_m, 
                                                                decreasing = FALSE)])

ccle.p = ccle.p_ind[,c("gene","ccle_proportion_m")]
ccle.p_ind$ccle_proportion_m = NULL

ccle.p_ind = reshape2::melt(ccle.p_ind)
colnames(ccle.p_ind) = c("gene","cancer_type","ccle.p")

ccle.p_ind$ccle.p[which(ccle.p_ind$ccle.p >= 75)] = 75
ccle.p_ind$ccle.p[which(ccle.p_ind$ccle.p <= 25)] = 25

heat_map_ccle = ggplot( ccle.p_ind ) +
  geom_tile( aes( x = cancer_type, y = gene, fill = ccle.p )  ) + 
  scale_fill_gradient2( low = "white", high = "#984ea3", name = "Exp.Prop", limits = c(25,75)) +
  scale_size(breaks=c(0,25,50,75,100),labels=c("0","25","50","75","100"),guide="legend") +
  theme_minimal()+ 
  theme(axis.text.y = element_text(color = "black", size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_text(angle = 45, vjust = 1.2, size = 7,color = "black", hjust = 1),
        axis.title.x = element_blank(),
        legend.position="left",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.box.spacing = unit(0.01, "cm"),
        legend.key.size = unit(0.5, "cm"),
        legend.key.width = unit(0.5,"cm"))

bars_ccle_p = ggplot( ccle.p, aes(x = gene, y = ccle_proportion_m)) + 
  geom_bar(stat = "identity") + 
  theme_classic() + 
  ylab("Exp.prop") +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.title.x = element_text(size = 7)
  )+ 
  ylim(c(0,100)) + 
  coord_flip()

########################### NT5 genes #################################
# lolipop plot

imp_genes = ccle_norm[which(ccle_norm$gene %in% c("NT5M","NT5E")),]

imp_genes = imp_genes[,c("gene","cancer_type","cancer","normal","ccle_m")]
imp_genes = reshape2::melt(imp_genes)

nt5_lp = ggplot(imp_genes, 
       aes(y = value, x = cancer_type, fill = variable, 
           colour = variable)) +
  geom_segment(aes(x = cancer_type, y = 0, xend = cancer_type, 
                   yend = value), color = "grey50", size = 0.5) +
  scale_color_manual(values = c("#b2182b","#33a02c","#984ea3"))+ 
  theme_minimal() + 
  ylab("log2 expr") + 
  theme(axis.text.x = element_text(color = "black",
                                   size = 7, angle = 45, hjust = 1),
        axis.title.x = element_text(size = 7),
        axis.text.y = element_text(color = "black",size = 7),
        axis.title.y = element_text(size = 7),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 7)) + 
  geom_point(size = 1) + 
  facet_grid(gene~. , scales = "free") + 
  theme(strip.text = element_text(size = 7))

grid.arrange(arrangeGrob(ggtern::arrangeGrob(a1,
                                 top = textGrob("Relative expression",gp=gpar(fontsize=8))),
                                 arrangeGrob(heat_map_ccle,bars_ccle_p,
                                 ncol = 2, widths = c(10,2)),
                                 ncol = 2, widths = c(1,2)), 
                     nt5_lp,
                     ncol = 1,
             top = textGrob("Cancer Cell Line Genes",gp=gpar(fontsize=10)))


