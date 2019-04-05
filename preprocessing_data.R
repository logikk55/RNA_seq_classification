library(dplyr, quietly = T)
library(tibble)
library(edgeR, quietly=T)


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

filter_classes <- function(counts, clinical, clinical_filter){
  clin_count_bind <- rbind(counts, t(clinical))
  data <- as.data.frame(t(clin_count_bind))
  
  filtered_counts <- data %>% 
    rownames_to_column('Sample') %>% 
    filter(pathologicstage==clinical_filter) %>% #switch the filter value based on what  youre filtering
    select(-seq(ncol(data)-ncol(clinical), ncol(data)+1)) %>%
    column_to_rownames('Sample')
  filtered_counts <- as.data.frame(t(filtered_counts))
  filtered_counts
}

#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
# Acute myeloid leukemia
load(file='~/LST/data rdata/LAML_clinical_RNASeq.RData')
aml_counts <- as.data.frame(rnaseq.aml$dat)
aml_clinical <- rnaseq.aml$clinical
processed_aml <- process(aml_counts, aml_clinical)
summary(processed_aml[[2]])

aml_counts <- processed_aml[[1]]
aml_counts <- remove_rowzeros(aml_counts)
aml_counts <- rowCPM_above_one(aml_counts)
aml_clinical <- processed_aml[[2]]
aml_male <- filter_classes(counts=aml_counts, clinical=aml_clinical, clinical_filter='male')
aml_female <- filter_classes(counts=aml_counts, clinical=aml_clinical, clinical_filter='female')
write.table(aml_male, file="~/LST/aml_raw_male.csv", sep = "\t")
write.table(aml_female, file="~/LST/aml_raw_female.csv", sep = "\t")


#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
# Breast invasive carcinoma
# load(file='~/LST/data rdata/BCRA_clinical_RNASeq.RData')
# brca_counts <- as.data.frame(rnaseq.brca$dat)
# brca_clinical <- rnaseq.brca$clinical
# processed_brca <- process(brca_counts, brca_clinical)

#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
# Lung adenocarcinoma
load(file='~/LST/data rdata/LUAD_clinical_RNASeq.RData')
luad_counts <- as.data.frame(rnaseq.luad$dat)
luad_clinical <- rnaseq.luad$clinical
processed_luad <- process(luad_counts, luad_clinical)
summary(processed_luad[[2]])

luad_counts <- processed_luad[[1]]
luad_counts <- remove_rowzeros(luad_counts)
luad_counts <- rowCPM_above_one(luad_counts)
luad_clinical <- processed_luad[[2]]

luad_adenocarcinoma_NS <- filter_classes(counts=luad_counts, clinical=luad_clinical, clinical_filter='lung adenocarcinoma- not otherwise specified (nos)')
luad_adenocarcinoma_mix <- filter_classes(counts=luad_counts, clinical=luad_clinical, clinical_filter='lung adenocarcinoma mixed subtype')
luad_bronchio_mix <- filter_classes(counts=luad_counts, clinical=luad_clinical, clinical_filter='lung bronchioloalveolar carcinoma nonmucinous')
luad_papillary_adeno <- filter_classes(counts=luad_counts, clinical=luad_clinical, clinical_filter='lung papillary adenocarcinoma')
luad_acinar_adeno <- filter_classes(counts=luad_counts, clinical=luad_clinical, clinical_filter='lung acinar adenocarcinoma')
luad_clear_adeno <- filter_classes(counts=luad_counts, clinical=luad_clinical, clinical_filter='lung clear cell adenocarcinoma')
write.table(luad_adenocarcinoma_NS, file="~/LST/luad_raw_adeno_not_specified.csv", sep = "\t")
write.table(luad_adenocarcinoma_mix, file="~/LST/luad_raw_adeno_mixed.csv", sep = "\t")
write.table(luad_bronchio_mix, file="~/LST/luad_raw_brochioloalveolar.csv", sep = "\t")
write.table(luad_papillary_adeno, file="~/LST/luad_raw_papillary_adeno.csv", sep = "\t")
write.table(luad_acinar_adeno, file="~/LST/luad_raw_acinar_adeno.csv", sep = "\t")
write.table(luad_clear_adeno, file="~/LST/luad_raw_clear_adeno.csv", sep = "\t")

#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
# Lungs squamous cell carcinoma
load(file='~/LST/data rdata/LUSC_clinical_RNASeq.RData')
lusc_counts <- as.data.frame(rnaseq.lusc$dat)
lusc_clinical <- rnaseq.lusc$clinical
processed_lusc <- process(lusc_counts, lusc_clinical)
summary(processed_lusc[[2]])

lusc_counts <- processed_lusc[[1]]
lusc_counts <- remove_rowzeros(lusc_counts)
lusc_counts <- rowCPM_above_one(lusc_counts)
lusc_clinical <- processed_lusc[[2]]

lusc_stage_ib <- filter_classes(counts=lusc_counts, clinical=lusc_clinical, clinical_filter='stage ib')
lusc_stage_ia <- filter_classes(counts=lusc_counts, clinical=lusc_clinical, clinical_filter='stage ia')
lusc_stage_iib <- filter_classes(counts=lusc_counts, clinical=lusc_clinical, clinical_filter='stage iib')
lusc_stage_iiia <- filter_classes(counts=lusc_counts, clinical=lusc_clinical, clinical_filter='stage iiia')
lusc_stage_iia <- filter_classes(counts=lusc_counts, clinical=lusc_clinical, clinical_filter='stage iia')
write.table(lusc_stage_ib, file="~/LST/lusc_raw_stage_ib.csv", sep = "\t")
write.table(lusc_stage_ia, file="~/LST/lusc_raw_stage_ia.csv", sep = "\t")
write.table(lusc_stage_iib, file="~/LST/lusc_raw_stage_iib.csv", sep = "\t")
write.table(lusc_stage_iiia, file="~/LST/lusc_raw_stage_iiia.csv", sep = "\t")
write.table(lusc_stage_iia, file="~/LST/lusc_raw_stage_iia.csv", sep = "\t")

#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
# # Ovarian serous cystadenocarcinoma
# load(file='~/LST/data rdata/OV_clinical_RNASeq.RData')
# ov_counts <- as.data.frame(rnaseq.ovarian$dat)
# ov_clinical <- rnaseq.ovarian$clinical
# processed_ov <- process(ov_counts, ov_clinical)

#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
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










