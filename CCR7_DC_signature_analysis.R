# =============================================================================
# CCR7⁺ DC Signature Analysis - Cutaneous + Acral Melanoma
# FULL REPRODUCIBLE SCRIPT (GitHub-ready)
# =============================================================================
# Author: Luis Castillo
# Date: April 2026
# Manuscript term: CCR7⁺ DC signature (technical name: mRegDC)
# =============================================================================

setwd("/Users/lcastillo6/Desktop/Melanoma human data analysis")

library(readxl)
library(data.table)
library(GSVA)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(survival)
library(survminer)
library(broom)

out_dir <- "mRegDC_Plots_Final_v4"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

signature_name <- "CCR7+ DC signature"

# 1. Load data
expr_full <- fread("RNA-CancerCell-MORRISON1-combat_batch_corrected-logcpm-all_samples.tsv", data.table = FALSE)
rownames(expr_full) <- expr_full[[1]]
expr_full <- expr_full[,-1]

meta_rna <- fread("RNA-CancerCell-MORRISON1-metadata.tsv", data.table = FALSE)
meta_subjects <- fread("Subjects-CancerCell-MORRISON1-metadata.tsv", data.table = FALSE)

mregdc_list <- read_excel("mRegDC gene list (Yang et al. 2025, Nat Commun).xlsx", sheet = 1)
mregdc_genes <- unique(na.omit(mregdc_list[[1]]))

# 2. Compute CCR7⁺ DC signature scores
expr_mat <- as.matrix(expr_full)
storage.mode(expr_mat) <- "numeric"

gsva_param <- gsvaParam(expr_mat, list(mRegDC = mregdc_genes), minSize = 10, maxSize = 1000)
ssgsea_scores <- gsva(gsva_param, verbose = TRUE)

mregdc_df <- data.frame(sample_id = colnames(ssgsea_scores),
                        mRegDC_score = as.numeric(ssgsea_scores["mRegDC", ]))

# 3. Merge + Filter
meta_merged <- meta_rna %>%
  left_join(mregdc_df, by = c("sample.id" = "sample_id")) %>%
  left_join(meta_subjects, by = "subject.id", suffix = c("", "_subj")) %>%
  filter(!is.na(mRegDC_score)) %>%
  filter(grepl("PD1|Combo", treatment.regimen.name, ignore.case = TRUE)) %>%
  filter(grepl("cutaneous|acral", sample.tumor.type, ignore.case = TRUE)) %>%
  mutate(
    Timepoint_Group = case_when(
      grepl("pre|baseline", timepoint.id, ignore.case = TRUE) ~ "Pre-treatment",
      grepl("on|post", timepoint.id, ignore.case = TRUE)      ~ "On_Post-treatment",
      TRUE ~ "Other"
    ),
    response_binary = ifelse(bor %in% c("CR", "PR"), "Responder", "Non-responder"),
    response_binary = factor(response_binary, levels = c("Non-responder", "Responder")),
    mRegDC_high_low = factor(ifelse(mRegDC_score >= median(mRegDC_score, na.rm = TRUE), "High", "Low"), 
                             levels = c("High", "Low")),
    Tumor_Type = sample.tumor.type
  )

cat("FINAL DATASET:", nrow(meta_merged), "samples\n")
print(table(meta_merged$Timepoint_Group, meta_merged$response_binary))

# 4. Mandatory DC Abundance Estimation (simple robust marker average)
cat("\n=== Computing DC Abundance (Marker Gene Average) ===\n")
dc_markers <- c("HLA-DRA", "CD1C", "CLEC9A", "BATF3", "IRF8", "ITGAX", "CD83")
dc_markers <- intersect(dc_markers, rownames(expr_mat))

if (length(dc_markers) > 0) {
  meta_merged$DC_abundance <- colMeans(expr_mat[dc_markers, meta_merged$sample.id, drop = FALSE], na.rm = TRUE)
  cat("Used", length(dc_markers), "DC marker genes for abundance estimation.\n")
} else {
  meta_merged$DC_abundance <- NA
  cat("Warning: No DC markers found.\n")
}

# 5. Styled Analysis Function
run_analysis <- function(data, title) {
  if(nrow(data) < 10) return(NULL)
  
  cat("\n=== ", title, " (n =", nrow(data), ") ===\n")
  wilcox_p <- wilcox.test(mRegDC_score ~ response_binary, data = data)$p.value
  cat("Wilcoxon p =", format.pval(wilcox_p, digits = 3), "\n")
  
  p_box <- ggplot(data, aes(x = response_binary, y = mRegDC_score, fill = response_binary)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.95, color = "#000000", linewidth = 0.95, fatten = 2.2) +
    geom_jitter(width = 0.3, size = 4.5, shape = 21, color = "black", fill = "white", stroke = 1.2) +
    stat_compare_means(method = "wilcox.test", label = "p.format", size = 5.5,
                       label.y = max(data$mRegDC_score, na.rm = TRUE) * 1.12) +
    scale_fill_manual(values = c("Responder" = "#e31a1c", "Non-responder" = "#000000")) +
    theme_minimal(base_size = 15) +
    labs(title = paste(signature_name, "by RECIST Response -", title),
         subtitle = paste("n =", nrow(data), "samples | Cutaneous + Acral only"),
         y = paste(signature_name, "ssGSEA Score"), x = "") +
    theme(legend.position = "none")
  
  print(p_box)
  ggsave(file.path(out_dir, paste0("CCR7_DC_Response_", gsub("/", "_", title), ".pdf")),
         p_box, width = 7, height = 6.5, dpi = 300)
  
  # KM curves
  data_surv <- data %>% filter(!is.na(os) | !is.na(pfs)) %>% mutate(OS_event = 1, PFS_event = 1)
  if(nrow(data_surv) > 20) {
    fit_os  <- survfit(Surv(os, OS_event) ~ mRegDC_high_low, data = data_surv)
    fit_pfs <- survfit(Surv(pfs, PFS_event) ~ mRegDC_high_low, data = data_surv)
    
    p_os <- ggsurvplot(fit_os, data = data_surv, pval = TRUE, risk.table = TRUE,
                       title = paste("Overall Survival by", signature_name, "-", title),
                       palette = c("#e31a1c", "#000000"),
                       pval.coord = c(0.05, 0.05), pval.size = 5.5)
    
    p_pfs <- ggsurvplot(fit_pfs, data = data_surv, pval = TRUE, risk.table = TRUE,
                        title = paste("Progression-Free Survival by", signature_name, "-", title),
                        palette = c("#e31a1c", "#000000"),
                        pval.coord = c(0.05, 0.05), pval.size = 5.5)
    
    print(p_os)
    print(p_pfs)
    ggsave(file.path(out_dir, paste0("CCR7_DC_OS_KM_", gsub("/", "_", title), ".pdf")), 
           plot = p_os$plot, width = 8, height = 6, dpi = 300)
    ggsave(file.path(out_dir, paste0("CCR7_DC_PFS_KM_", gsub("/", "_", title), ".pdf")), 
           plot = p_pfs$plot, width = 8, height = 6, dpi = 300)
  }
}

# Run plots
run_analysis(filter(meta_merged, Timepoint_Group == "Pre-treatment"), "Pre-treatment")
run_analysis(filter(meta_merged, Timepoint_Group == "On_Post-treatment"), "On_Post-treatment")

# 6. Export final table
summary_table_final <- meta_merged %>%
  select(sample.id, Timepoint_Group, treatment.regimen.name, bor, response,
         cohort, os, pfs, mRegDC_score, previous.treatment, Tumor_Type, DC_abundance) %>%
  arrange(Timepoint_Group, bor) %>%
  rename(Sample_ID = sample.id, Timepoint = Timepoint_Group,
         Treatment = treatment.regimen.name, BOR = bor, Cohort = cohort,
         OS_days = os, PFS_days = pfs, Prior_treatment = previous.treatment)

write.csv(summary_table_final, "mRegDC_PD1_Focused_Summary_Table_CUTANEOUS_ACRAL.csv", row.names = FALSE)

cat("\n✅ Full analysis complete! All outputs use 'CCR7+ DC signature' nomenclature.\n")
cat("Files saved in:", out_dir, "\n")
sessionInfo()
# =============================================================================
# Create Supplementary Table: Missing mRegDC Genes
# =============================================================================

cat("\n=== Creating Supplementary Table for Missing Genes ===\n")

mregdc_list <- read_excel("mRegDC gene list (Yang et al. 2025, Nat Commun).xlsx", sheet = 1)
mregdc_genes <- unique(na.omit(mregdc_list[[1]]))

expr_genes <- rownames(expr_full)   # already loaded earlier in the script

overlap_genes <- intersect(mregdc_genes, expr_genes)
missing_genes <- setdiff(mregdc_genes, overlap_genes)

supp_missing <- data.frame(
  Missing_Gene = missing_genes,
  stringsAsFactors = FALSE
)

write.csv(supp_missing, "Supplementary_Table_Missing_mRegDC_genes.csv", row.names = FALSE)

cat("Supplementary table created!\n")
cat("File saved as: Supplementary_Table_Missing_mRegDC_genes.csv\n")
cat("Total missing genes:", nrow(supp_missing), "\n")