
ExprReactions = read.table(file = '../data/CCLE/metabolic_data/min_max_new_model/ExprReactions_All.dat', 
                                sep = ',', header = TRUE, stringsAsFactors = FALSE)

ExprReactions$cancer.score = ExprReactions$cancer - ExprReactions$normal
ExprReactions$stroma.score = ExprReactions$stroma - ExprReactions$normal
ExprReactions$cs = ExprReactions$cancer - ExprReactions$stroma
ExprReactions$cancer_ccle.score = ExprReactions$cancer - ExprReactions$CCLE
ExprReactions$ccle.score = ExprReactions$CCLE - ExprReactions$normal

### Pan-Cancer

# Cancer
cancer = ExprReactions[-which(ExprReactions$cancer_type %in% c("LGG","GBM")),
                       c("Reaction.Abbreviation","cancer.score","cancer_type")]
cancer = dcast(cancer, Reaction.Abbreviation ~ cancer_type, value.var = "cancer.score")
cancer$m = apply(cancer[,-c(1)], 1, FUN = function(x) median(x, na.rm = TRUE))
cancer = cancer[,c("Reaction.Abbreviation","m")]

# write.csv(cancer,file = "../data/CCLE/metabolic_data/escher_maps/rxns_input/cancer.txt",
          # quote = FALSE,
          # row.names = FALSE)

# stroma
stroma = ExprReactions[-which(ExprReactions$cancer_type %in% c("LGG","GBM")),
                       c("Reaction.Abbreviation","stroma.score","cancer_type")]
stroma = dcast(stroma, Reaction.Abbreviation ~ cancer_type, value.var = "stroma.score")
stroma$m = apply(stroma[,-c(1)], 1, FUN = function(x) median(x, na.rm = TRUE))
stroma = stroma[,c("Reaction.Abbreviation","m")]

# write.csv(stroma,file = "../data/CCLE/metabolic_data/escher_maps/rxns_input/stroma.txt",
          # quote = FALSE,
          # row.names = FALSE)

# CS
cs.df = ExprReactions[-which(ExprReactions$cancer_type %in% c("LGG","GBM")),
                      c("Reaction.Abbreviation","cs","cancer_type")]
cs.df = dcast(cs.df, Reaction.Abbreviation ~ cancer_type, value.var = "cs")
cs.df$m = apply(cs.df[,-c(1)], 1, FUN = function(x) median(x, na.rm = TRUE))
cs.df = cs.df[,c("Reaction.Abbreviation","m")]

# write.csv(cs.df,file = "../data/CCLE/metabolic_data/escher_maps/rxns_input/cs.txt",
          # quote = FALSE,
          # row.names = FALSE)

# Cancer vs CCLE
# No KIRP in CCLE
cancer_ccle = ExprReactions[-which(ExprReactions$cancer_type %in% c("KIRP","LGG","GBM")),
                            c("Reaction.Abbreviation","cancer_ccle.score","cancer_type")]
cancer_ccle = dcast(cancer_ccle, 
                    Reaction.Abbreviation ~ cancer_type, 
                    value.var = "cancer_ccle.score")
cancer_ccle$m = apply(cancer_ccle[,-c(1)], 1, 
                      FUN = function(x) median(x, na.rm = TRUE))
cancer_ccle = cancer_ccle[,c("Reaction.Abbreviation","m")]

# write.csv(cancer_ccle,file = "../data/CCLE/metabolic_data/escher_maps/rxns_input/cancer_ccle.txt",
          # quote = FALSE,
          # row.names = FALSE)

# CCLE
ccle = ExprReactions[-which(ExprReactions$cancer_type %in% c("KIRP","LGG","GBM","CRC")),
                     c("Reaction.Abbreviation","ccle.score","cancer_type")]
ccle = dcast(ccle, Reaction.Abbreviation ~ cancer_type, value.var = "ccle.score")
ccle$m = apply(ccle[,-c(1)], 1, 
               FUN = function(x) median(x, na.rm = TRUE))
ccle = ccle[,c("Reaction.Abbreviation","m")]

# write.csv(ccle,file = "../data/CCLE/metabolic_data/escher_maps/rxns_input/ccle.txt",
          # quote = FALSE,
          # row.names = FALSE)



#### for transporters in the map

ExprMetabolic = read.table(file = '../data/CCLE/metabolic_data/Expression_Model_All.dat',
                           sep = ',', header = TRUE, stringsAsFactors = FALSE)

ExprMetabolic = ExprMetabolic[-which(ExprMetabolic$cancer_type %in% c("LGG","GBM")),]

ExprMetabolic$ReconGeneName = NULL
ExprMetabolic = ExprMetabolic[!duplicated(ExprMetabolic),]
ExprMetabolic$cn = ExprMetabolic$cancer - ExprMetabolic$normal
ExprMetabolic$sn = ExprMetabolic$stroma - ExprMetabolic$normal