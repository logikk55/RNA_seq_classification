library(TCGA2STAT)

# NOPE
#rnaseq.laml <- getTCGA(disease="LAML", data.type="RNASeq", type="count", clinical=TRUE)

# rnaseq.acc <- getTCGA(disease="ACC", data.type="RNASeq2", type="count", clinical=TRUE)
# summary(rnaseq.acc$clinical)
# # very unbalanced histological data
# 
# rnaseq.blca <- getTCGA(disease="BLCA", data.type="RNASeq2", type="count", clinical=TRUE)
# summary(rnaseq.blca$clinical)
# # very unbalanced hist

rnaseq.brca <- getTCGA(disease="BRCA", data.type="RNASeq2", type="count", clinical=TRUE)
summary(rnaseq.brca$clinical)
# good hist

rnaseq.lgg <- getTCGA(disease="LGG", data.type="RNASeq2", type="count", clinical=TRUE)
summary(rnaseq.lgg$clinical)
# very good hist, 194, 130, 191

# rnaseq.cesc <- getTCGA(disease="CESC", data.type="RNASeq2", type="count", clinical=TRUE)
# summary(rnaseq.cesc$clinical)
# # poor hist, 6, 254, 6, 21, 3, 17

# rnaseq.chol <- getTCGA(disease="CHOL", data.type="RNASeq2", type="count", clinical=TRUE)
# summary(rnaseq.chol$clinical)
# # poor hist, less < 50 samples. 36 to one type, 7 to one, 2 to one

# rnaseq.coad <- getTCGA(disease="COAD", data.type="RNASeq2", type="count", clinical=TRUE)
# summary(rnaseq.coad$clinical)
# # 391 to one type, 62 to one, Descent

rnaseq.coadread <- getTCGA(disease="COADREAD", data.type="RNASeq2", type="count", clinical=TRUE)
summary(rnaseq.coadread$clinical)
# good hist: 391, 62, 152, 13

# rnaseq.dlbc <- getTCGA(disease="DLBC", data.type="RNASeq2", type="count", clinical=TRUE)
# summary(rnaseq.dlbc$clinical)
# # poor hist: 41, 3, 4

# rnaseq.gbm <- getTCGA(disease="GBM", data.type="RNASeq2", type="count", clinical=TRUE)
# summary(rnaseq.gbm$clinical)
# # poor hist: 31, 20, 544

# rnaseq.gbmlgg <- getTCGA(disease="GBMLGG", data.type="RNASeq2", type="count", clinical=TRUE)
# summary(rnaseq.gbmlgg$clinical)
# GMBLGG# good hist: 194, 31, 130, 191 <<<----- same set as the LGG

# rnaseq.hnsc <- getTCGA(disease="HNSC", data.type="RNASeq2", type="count", clinical=TRUE)
# summary(rnaseq.hnsc$clinical)
# # poorhist: 517, 10, 1

rnaseq.kipan <- getTCGA(disease="KIPAN", data.type="RNASeq2", type="count", clinical=TRUE)
summary(rnaseq.kipan$clinical)
# good hist: 113, 537, 291

# rnaseq.lihc <- getTCGA(disease="LIHC", data.type="RNASeq2", type="count", clinical=TRUE)
# summary(rnaseq.lihc$clinical)
# # poor hist: 3, 367, 7

# rnaseq.luad <- getTCGA(disease="LUAD", data.type="RNASeq2", type="count", clinical=TRUE)
# summary(rnaseq.luad$clinical)
# # poor hist types and numbers

# rnaseq.lusc <- getTCGA(disease="LUSC", data.type="RNASeq2", type="count", clinical=TRUE)
# summary(rnaseq.lusc$clinical)
# # poor hist

rnaseq.meso <- getTCGA(disease="MESO", data.type="RNASeq2", type="count", clinical=TRUE)
summary(rnaseq.meso$clinical)
# descent hist: 23, 5, 57, 2

# rnaseq.ov <- getTCGA(disease="OV", data.type="RNASeq2", type="count", clinical=TRUE)
# summary(rnaseq.ov$clinical)
# # poor hist

# rnaseq.paad <- getTCGA(disease="PAAD", data.type="RNASeq2", type="count", clinical=TRUE)
# summary(rnaseq.paad$clinical)
# # poor histtypes and numbers

# rnaseq.pcpg <- getTCGA(disease="PCPG", data.type="RNASeq2", type="count", clinical=TRUE)
# summary(rnaseq.pcpg$clinical)
# # poor histtypes and numbers: 18, 13, 148

# rnaseq.prad <- getTCGA(disease="PRAD", data.type="RNASeq2", type="count", clinical=TRUE)
# summary(rnaseq.prad$clinical)
# # poor histtypes and numbers

# rnaseq.read <- getTCGA(disease="READ", data.type="RNASeq2", type="count", clinical=TRUE)
# summary(rnaseq.read$clinical)
# # poor histtypes and numbers

# rnaseq.read <- getTCGA(disease="READ", data.type="RNASeq2", type="count", clinical=TRUE)
# summary(rnaseq.read$clinical)
# # poor histtypes and numbers

rnaseq.sarc <- getTCGA(disease="SARC", data.type="RNASeq2", type="count", clinical=TRUE)
summary(rnaseq.sarc$clinical)
# descent hist types: 105, 59, 29, 25, 21, 9

# rnaseq.skcm <- getTCGA(disease="SKCM", data.type="RNASeq2", type="count", clinical=TRUE)
# summary(rnaseq.skcm$clinical)
# # poor hist types

# rnaseq.stad <- getTCGA(disease="STAD", data.type="RNASeq2", type="count", clinical=TRUE)
# summary(rnaseq.stad$clinical)
# # poor hist types

# rnaseq.tgct <- getTCGA(disease="TGCT", data.type="RNASeq2", type="count", clinical=TRUE)
# summary(rnaseq.tgct$clinical)
# # poor hist types

rnaseq.thca <- getTCGA(disease="THCA", data.type="RNASeq2", type="count", clinical=TRUE)
summary(rnaseq.thca$clinical)
# descent hist: 358, 102, 36

rnaseq.thym <- getTCGA(disease="THYM", data.type="RNASeq2", type="count", clinical=TRUE)
summary(rnaseq.thym$clinical)
# descent hist: 17, 38, 15, 31, 12, 11 (balanced)

rnaseq.ucec <- getTCGA(disease="UCEC", data.type="RNASeq2", type="count", clinical=TRUE)
summary(rnaseq.ucec$clinical)
# descent hist: 411, 22, 115

rnaseq.ucs <- getTCGA(disease="UCS", data.type="RNASeq2", type="count", clinical=TRUE)
summary(rnaseq.ucs$clinical)
# descent hist: 24, 20, 13

# rnaseq.uvm <- getTCGA(disease="UVM", data.type="RNASeq2", type="count", clinical=TRUE)
# summary(rnaseq.uvm$clinical)
# # poor hist

















