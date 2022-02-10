library(reshape2)
library(cowplot)
library(ggrepel)
library(ggtern)
library(grid)
library(gridExtra)
library(gtable)
library(KEGGREST)
library(fgsea)

ExprMetabolic = read.table(file = '../data/Expression_Model_All.dat',
                       sep = ',', header = TRUE, stringsAsFactors = FALSE)

ExprMetabolic$ReconGeneName = NULL
ExprMetabolic = ExprMetabolic[!duplicated(ExprMetabolic),]

ExprMetabolic$cancer.raw = (2^ExprMetabolic$cancer)-1
ExprMetabolic$stroma.raw = (2^ExprMetabolic$stroma)-1
ExprMetabolic$normal.raw = (2^ExprMetabolic$normal)-1

C = ExprMetabolic[,c("gene","cancer_type","cancer.raw")]
C = reshape2::dcast(C, gene ~ cancer_type, value.var = "cancer.raw")
C$Cm = apply(C[,-1], 1, FUN = function(x) median(x))

S = ExprMetabolic[,c("gene","stroma.raw","cancer_type")]
S = reshape2::dcast(S, gene ~ cancer_type, value.var = "stroma.raw")
S$Sm = apply(S[,-1], 1, FUN = function(x) median(x))

N = ExprMetabolic[,c("gene","normal.raw","cancer_type")]
N = reshape2::dcast(N, gene ~ cancer_type, value.var = "normal.raw")
N$Nm = apply(N[,-1], 1, FUN = function(x) median(x))

combinedMedian = merge(N,C,by = "gene")
combinedMedian = combinedMedian[,c("gene","Nm","Cm")]
combinedMedian = merge(combinedMedian, S, by = "gene" )
combinedMedian = combinedMedian[,c("gene","Nm","Cm","Sm")]

combinedMedian = combinedMedian[which(combinedMedian$Nm  > 1|
                                        combinedMedian$Cm > 1|
                                        combinedMedian$Sm > 1),]

########################################  ternary plot  ###########################################
ternPlot = combinedMedian
 
tern.colors = list()
tern.colors[""] = "#000000"
tern.colors["cancer"] = "#b2182b"
tern.colors["stroma"] = "#2166ac"
tern.colors["normal"] = "#33a02c"
tern.colors = unlist(tern.colors)

combinedMedian$sum = combinedMedian$Nm + combinedMedian$Cm + combinedMedian$Sm
combinedMedian$normal.p = (combinedMedian$Nm/combinedMedian$sum)*100
combinedMedian$cancer.p = (combinedMedian$Cm/combinedMedian$sum)*100
combinedMedian$stroma.p = (combinedMedian$Sm/combinedMedian$sum)*100

stroma_genes = as.character(combinedMedian$gene[order(combinedMedian$stroma.p, 
                                                      decreasing = TRUE)[1:15]])
cancer_genes = as.character(combinedMedian$gene[order(combinedMedian$cancer.p, 
                                                      decreasing = TRUE)[1:15]])

ternPlot$cl = ""
ternPlot$cl[which(ternPlot$gene %in% stroma_genes)] = "stroma"
ternPlot$cl[which(ternPlot$gene %in% cancer_genes)] = "cancer"

ternPlot$al = 0.9
ternPlot$al[which(ternPlot$gene %in% stroma_genes)] = 1
ternPlot$al[which(ternPlot$gene %in% cancer_genes)] = 1

# For Rendering the Lines, use Segment geometry
lines <- data.frame(x = c(0.5, 0, 0.5), 
                    y = c(0.5, 0.5, 0), 
                    z = c(0, 0.5, 0.5), 
                    xend = c(1, 1, 1)/3, 
                    yend = c(1, 1, 1)/3, 
                    zend = c(1, 1, 1)/3)

ggtern(data = ternPlot, aes(x = Cm, y = Nm, z = Sm)) + 
  geom_point(size=2.5,aes(colour=cl, alpha = al)) + Tlab("N") + Llab("C") + Rlab("S") +
  limit_tern(1.1,1.1,1.1) + 
  scale_color_manual(values = tern.colors)  + theme_bw() + 
  theme(axis.title = element_text(colour = "black", size = 7),
        axis.text = element_text(colour = "black", size = 7),
        legend.position = "none") + 
  geom_segment(data = lines, 
               aes(x, y, z, xend = xend, yend = yend, zend = zend), 
               color = 'red', size = 1) + 
  ggtitle("Relative Expression")

########################################### HEAT MAP #############################################
ExprMetabolic$sum = ExprMetabolic$cancer.raw + ExprMetabolic$stroma.raw + ExprMetabolic$normal.raw
ExprMetabolic$cancer.p = (ExprMetabolic$cancer.raw/ExprMetabolic$sum)*100
ExprMetabolic$stroma.p = (ExprMetabolic$stroma.raw/ExprMetabolic$sum)*100
ExprMetabolic$normal.p = (ExprMetabolic$normal.raw/ExprMetabolic$sum)*100

#cancer
cancer_proportion = ExprMetabolic[,c("gene","cancer_type","cancer.p")]
cancer_proportion = reshape2::dcast(cancer_proportion, gene ~ cancer_type, value.var = "cancer.p")
cancer_proportion[is.na(cancer_proportion)] <- 0
cancer_proportion$cancer_proportion_m = apply(cancer_proportion[,-1], 1, FUN = function(x) median(x))
cancer_proportion = cancer_proportion[which(cancer_proportion$gene %in% combinedMedian$gene),]

c.p_ind = cancer_proportion[which(cancer_proportion$gene %in% cancer_genes),]
c.p_ind$gene = factor(c.p_ind$gene, levels = c.p_ind$gene[order(c.p_ind$cancer_proportion_m, decreasing = FALSE)])

c.p = c.p_ind[,c("gene","cancer_proportion_m")]
c.p_ind$cancer_proportion_m = NULL

c.p_ind = reshape2::melt(c.p_ind)
colnames(c.p_ind) = c("gene","cancer_type","cancer.p")

c.p_ind$cancer.p[which(c.p_ind$cancer.p <= 25)] = 25
c.p_ind$cancer.p[which(c.p_ind$cancer.p >= 75)] = 75

heat_map_cancer = ggplot( c.p_ind ) +
  geom_tile( aes( x = cancer_type, y = gene, fill = cancer.p )  ) + 
  scale_fill_gradient2( low = "white", high = "#b2182b", na.value="", name = "Exp.Prop", limits = c(25,75)) +
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

bars_cancer_p = ggplot( c.p, aes(x = gene, y = cancer_proportion_m)) + 
  geom_bar(stat = "identity") + 
  theme_classic() + 
  ylab("Exp.prop") + 
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.title.x = element_text(size = 7))+ 
  ylim(c(0,100)) + 
  coord_flip()

#stroma
stroma_proportion = ExprMetabolic[,c("gene","stroma.p","cancer_type")]
stroma_proportion = reshape2::dcast(stroma_proportion, gene ~ cancer_type, value.var = "stroma.p")
stroma_proportion[is.na(stroma_proportion)] <- 0
stroma_proportion$stroma_proportion_m = apply(stroma_proportion[,-1], 1, FUN = function(x) median(x))

s.p_ind = stroma_proportion[which(stroma_proportion$gene %in% stroma_genes),]
s.p_ind$gene = factor(s.p_ind$gene, levels = s.p_ind$gene[order(s.p_ind$stroma_proportion_m, decreasing = FALSE)])

s.p = s.p_ind[,c("gene","stroma_proportion_m")]
s.p_ind$stroma_proportion_m = NULL

s.p_ind = reshape2::melt(s.p_ind)
colnames(s.p_ind) = c("gene","cancer_type","stroma.p")

s.p_ind$stroma.p[which(s.p_ind$stroma.p <= 25)] = 25
s.p_ind$stroma.p[which(s.p_ind$stroma.p >= 75)] = 75

heat_map_stroma = ggplot( s.p_ind ) +
  geom_tile( aes( x = cancer_type, y = gene, fill = stroma.p )  ) + 
  scale_fill_gradient2( low = "white", high = "#2166ac", name = "Exp.Prop", limits = c(25,75)) +
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

bars_stroma_p = ggplot( s.p, aes(x = gene, y = stroma_proportion_m)) + 
  geom_bar(stat = "identity") + 
  theme_classic() + 
  ylab("Exp.prop") + 
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.title.x = element_text(size = 7))+
  ylim(c(0,100)) + 
  # scale_y_continuous(breaks=seq(0,50,100)) + 
  coord_flip()

##############################################################################
###                          Lollipop Plots                                ###
##############################################################################

stroma_genes = ExprMetabolic[which(ExprMetabolic$gene %in% c("IDO1","IL4I1","TDO2")),]

stroma_genes = stroma_genes[,c("gene","cancer_type","stroma","cancer","normal")]
stroma_genes = reshape2::melt(stroma_genes)

lollypop_stroma = ggplot(stroma_genes, aes(y = value, x = cancer_type, fill = variable, colour = variable)) +
  geom_segment(aes(x = cancer_type, y = 0, xend = cancer_type, yend = value), color = "grey50", size = 0.5) +
  scale_color_manual(values = c("#2166ac","#b2182b","#33a02c"))+ theme_minimal() + 
  ylab("log2 expr") + 
  theme(axis.text.x = element_text(color = "black",size = 7, angle = 45, hjust = 1),
        axis.title.x = element_text(size = 7),
        axis.text.y = element_text(color = "black",size = 7),
        axis.title.y = element_text(size = 7),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 7)) + 
  geom_point(size = 1) + 
  facet_grid(gene~. , scales = "free") + 
  theme(strip.text = element_text(size = 7))


cancer_genes = ExprMetabolic[which(ExprMetabolic$gene %in% c("RRM2","SLC25A10","NT5M")),]

cancer_genes = cancer_genes[,c("gene","cancer_type","stroma","cancer","normal")]
cancer_genes = reshape2::melt(cancer_genes)

lollypop_cancer = ggplot(cancer_genes, aes(y = value, x = cancer_type, fill = variable, colour = variable)) +
  geom_segment(aes(x = cancer_type, y = 0, xend = cancer_type, yend = value), color = "grey50", size = 0.5) +
  scale_color_manual(values = c("#2166ac","#b2182b","#33a02c"))+ theme_minimal() + 
  ylab("log2 expr") + 
  theme(axis.text.x = element_text(color = "black",size = 7, angle = 45, hjust = 1),
        axis.title.x = element_text(size = 7),
        axis.text.y = element_text(color = "black",size = 7),
        axis.title.y = element_text(size = 7),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 7)) + 
  geom_point(size = 1) + 
  facet_grid(gene~. , scales = "free") + 
  theme(strip.text = element_text(size = 7))


############################################################################################
###               Ternary plot Mitchondrial Carrier (MC, SLC25) Family                   ###
############################################################################################
# List of mitochondrial carriers:
#        http://www.tcdb.org/search/result.php?tc=2.A.29#ref34546069

slc25_tern = ternPlot[grep(pattern = "SLC25",x=ternPlot$gene),]
slc_tern = ggtern(data = slc25_tern, aes(x = Cm, y = Nm, z = Sm)) + 
  geom_point(size=2.5, aes(color = cl)) + Tlab("N") + Llab("C") + Rlab("S") +
  scale_color_manual(values = c("black","#b2182b")) + 
  limit_tern(1.1,1.1,1.1) + 
  theme_bw() + 
  theme(axis.title = element_text(colour = "black", size = 7),
        axis.text = element_text(colour = "black", size = 7),
        legend.position = "none") + 
  geom_segment(data = lines, 
               aes(x, y, z, xend = xend, yend = yend, zend = zend), 
               color = 'red', size = 1) + 
  geom_text(data = slc25_tern[which(slc25_tern$cl != ""),], 
            aes(label = slc25_tern$gene[which(slc25_tern$cl != "")]),
            hjust = 0.5, vjust = 0.5, size = 3)

############################################################################################
###                       Ternary plot Tryptophan metabolism genes                       ###
############################################################################################
# my_pathways <- gmtPathways("C:/Projects/GSEA/c2.cp.kegg.v6.2.symbols.gmt")
# tryp_genes = my_pathways$KEGG_TRYPTOPHAN_METABOLISM
tryp_genes <- c("MAOB","MAOA","IDO2","AOX1","ALDH1B1","AANAT","ALDH2",
                "WARS","IDO1","CAT", "ACAT2","ACAT1","IL4I1","HADH",
                "OGDH","TPH1","HADHA","DDC","AFMID","CYP1A1","CYP1A2",
                "CYP1B1","ASMT","ECHS1","ALDH9A1","KMO","ALDH3A2","WARS2",
                "EHHADH","OGDHL","GCDH","ALDH7A1","ABP1","INMT","TDO2",
                "HAAO","KYNU","AADAT","TPH2","ACMSD")

tryp_tern = ternPlot[which(ternPlot$gene %in% tryp_genes),]
tryp_TernPlot = ggtern(data = tryp_tern, aes(x = Cm, y = Nm, z = Sm)) + 
  geom_point(size=2.5, aes(color = cl)) + Tlab("N") + Llab("C") + Rlab("S") +
  scale_color_manual(values = tern.colors) + 
  limit_tern(1.1,1.1,1.1) + 
  theme_bw() + 
  theme(axis.title = element_text(colour = "black", size = 7),
        axis.text = element_text(colour = "black", size = 7),
        legend.position = "none") + 
  geom_segment(data = lines, 
               aes(x, y, z, xend = xend, yend = yend, zend = zend), 
               color = 'red', size = 1) + 
  geom_text(data = tryp_tern[which(tryp_tern$cl != ""),], aes(label = tryp_tern$gene[which(tryp_tern$cl != "")]),
            hjust = 0.5, vjust = 0.5, size = 3)


#####################################################################################
###                                   Combine plots                               ###
#####################################################################################
lay = rbind(c(1,1),
            c(2,3))

grid.arrange(arrangeGrob(heat_map_cancer,bars_cancer_p,
                                 ncol = 2, widths = c(10,1.5),
                                 top = textGrob("Cancer Genes",
                                                gp=gpar(fontsize=8))),
                     lollypop_cancer,
             ggtern::arrangeGrob(slc_tern,
                                 top = textGrob("Mitonchondrial Carrier Genes",
                                                gp=gpar(fontsize=8))),
                     layout_matrix = lay,
                     top = textGrob("Cancer specific genes",gp=gpar(fontsize=10)))


grid.arrange(arrangeGrob(heat_map_stroma,bars_stroma_p,
                                 ncol = 2, widths = c(10,1.5),
                                 top = textGrob("Stromal Genes",
                                                gp=gpar(fontsize=8))),
                     lollypop_stroma,
             ggtern::arrangeGrob(tryp_TernPlot,
                                 top = textGrob("Tryptophan Metabolism Genes",
                                                gp=gpar(fontsize=8))),
                     layout_matrix = lay,
                     top = textGrob("Stroma specific genes",gp=gpar(fontsize=10)))

