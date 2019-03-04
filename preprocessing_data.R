process <- function(df_counts, df_clinic){
  counts <- df_counts
  clinic <- df_clinic
  names_clin <- rownames(clinic)
  colnames(counts) <- substring(colnames(counts),1,12)
  matches <- match(names_clin, colnames(counts))
  matches <- matches[!is.na(matches)]
  counts <- counts[, matches]
  clinic <- clinic[which(names_clin %in% colnames(counts)), ]
  return(list(counts, clinic))
}

remove_rowzeros <- function(counts){
  trh <- counts > 0
  keep <- rowSums(trh) > 0
  counts <- counts[keep, ]
  counts
}

rowCPM_above_one <- function(counts){
  library(edgeR, quietly=T)
  y <- DGEList(counts)
  CPM <- cpm(y, keep.lib.sizes=FALSE)
  trh <- CPM > 1
  keep <- rowSums(trh) >= 1
  y <- DGEList(y[keep, ])
  as.data.frame(y$counts)
}

log2_pseudo_counts <- function(counts){
  counts <- log2(counts+1)
  counts
}

# Acute myeloid leukemia
# load(file='~/LST/data rdata/LAML_clinical_RNASeq.RData')
# aml_counts <- as.data.frame(rnaseq.aml$dat)
# aml_clinical <- rnaseq.aml$clinical
# processed_aml <- process(aml_counts, aml_clinical)

# Breast invasive carcinoma
load(file='~/LST/data rdata/BCRA_clinical_RNASeq.RData')
brca_counts <- as.data.frame(rnaseq.brca$dat)
brca_clinical <- rnaseq.brca$clinical
processed_brca <- process(brca_counts, brca_clinical)

# Lung adenocarcinoma
# load(file='~/LST/data rdata/LUAD_clinical_RNASeq.RData')
# luad_counts <- as.data.frame(rnaseq.luad$dat)
# luad_clinical <- rnaseq.luad$clinical
# processed_luad <- process(luad_counts, luad_clinical)
# 
# # Lungs squamous cell carcinoma
# load(file='~/LST/data rdata/LUSC_clinical_RNASeq.RData')
# lusc_counts <- as.data.frame(rnaseq.lusc$dat)
# lusc_clinical <- rnaseq.lusc$clinical
# processed_lusc <- process(lusc_counts, lusc_clinical)
# 
# # Ovarian serous cystadenocarcinoma
# load(file='~/LST/data rdata/OV_clinical_RNASeq.RData')
# ov_counts <- as.data.frame(rnaseq.ovarian$dat)
# ov_clinical <- rnaseq.ovarian$clinical
# processed_ov <- process(ov_counts, ov_clinical)
# 
# # Thyroid carcinoma
# load(file='~/LST/data rdata/THCA_clinical_RNASeq.RData')
# thca_counts <- as.data.frame(rnaseq.thca$dat)
# thca_clinical <- rnaseq.thca$clinical
# processed_thca <- process(thca_counts, thca_clinical)


#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
# Remove non-expressing genes
brca_counts <- processed_brca[[1]]
brca_counts <- remove_rowzeros(brca_counts)
brca_counts <- rowCPM_above_one(brca_counts)

# Transform to log2-pseudo-counts
# brca_counts <- log2_pseudo_counts(brca_counts)

#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
# FILTER OUT THE SUBTYPES
library(dplyr, quietly = T)
library(tibble)

# summary(processed_brca[[2]])

clin_count_bind <- rbind(brca_counts, t(processed_brca[[2]]))
data <- as.data.frame(t(clin_count_bind))

ductal_brca <- data %>% 
  rownames_to_column('Sample') %>% 
  filter(histologicaltype=='infiltrating ductal carcinoma') %>% 
  select(-seq(ncol(data)-16, ncol(data)+1)) %>%
  column_to_rownames('Sample')
ductal_brca <- as.data.frame(t(ductal_brca))
write.table(ductal_brca, file="~/LST/ductal_brca_raw.csv", sep = "\t")

lobular_brca <- data %>% 
  rownames_to_column('Sample') %>% 
  filter(histologicaltype=='infiltrating lobular carcinoma') %>% 
  select(-seq(ncol(data)-16, ncol(data)+1)) %>%
  column_to_rownames('Sample')
lobular_brca <- as.data.frame(t(lobular_brca))
write.table(lobular_brca, file="~/LST/lobular_brca_raw.csv", sep = "\t")

mucinous_brca <- data %>% 
  rownames_to_column('Sample') %>% 
  filter(histologicaltype=='mucinous carcinoma') %>% 
  select(-seq(ncol(data)-16, ncol(data)+1)) %>%
  column_to_rownames('Sample')
mucinous_brca <- as.data.frame(t(mucinous_brca))
write.table(mucinous_brca, file="~/LST/mucinous_brca_raw.csv", sep = "\t")

mixed_brca <- data %>% 
  rownames_to_column('Sample') %>% 
  filter(histologicaltype=='mixed histology (please specify)') %>% 
  select(-seq(ncol(data)-16, ncol(data)+1)) %>%
  column_to_rownames('Sample')
mixed_brca <- as.data.frame(t(mixed_brca))
write.table(mixed_brca, file="~/LST/mixed_brca_raw.csv", sep = "\t")

others_brca <- data %>% 
  rownames_to_column('Sample') %>% 
  filter(histologicaltype=='other, specify') %>% 
  select(-seq(ncol(data)-16, ncol(data)+1)) %>%
  column_to_rownames('Sample')
others_brca <- as.data.frame(t(others_brca))
write.table(others_brca, file="~/LST/others_brca_raw.csv", sep = "\t")


#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
# SUBTYPES
summary(processed_aml[[2]])
summary(processed_brca[[2]])
summary(processed_luad[[2]])
summary(processed_lusc[[2]])
summary(processed_ov[[2]])
summary(processed_thca[[2]])

#################################################################
#################################################################
#################################################################
#################################################################
#################################################################











