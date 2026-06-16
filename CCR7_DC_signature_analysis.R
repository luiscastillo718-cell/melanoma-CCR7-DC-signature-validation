# =============================================================================
# CCR7⁺ DC Signature Analysis - ALL Subtypes (Pre & On/Post-treatment)
# Final clean version - Dual mRegDC signatures only
# =============================================================================
# Ready to run in a fresh RStudio session
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
expr_full <- fread("RNA-CancerCell-MORRISON1-combat_batch_corrected-logcpm-all_samples.tsv",
                   data.table = FALSE)
rownames(expr_full) <- expr_full[[1]]
expr_full <- expr_full[,-1]

meta_rna <- fread("RNA-CancerCell-MORRISON1-metadata.tsv", data.table = FALSE)
meta_subjects <- fread("Subjects-CancerCell-MORRISON1-metadata.tsv", data.table = FALSE)

# Primary mRegDC list
mregdc_list <- read_excel("mRegDC gene list (Yang et al. 2025, Nat Commun).xlsx", sheet = 1)
mregdc_genes <- unique(na.omit(mregdc_list[[1]]))

# Maier et al. 2020 mRegDC signature
maier_mregdc_list <- read_excel("/Users/lcastillo6/Desktop/Melanoma human data analysis/Human DC genes validated by Maier et al. Nature 2020.xlsx", sheet = 1)
maier_mregdc_genes <- unique(na.omit(maier_mregdc_list[[1]]))

cat("Primary mRegDC genes:", length(mregdc_genes), "\n")
cat("Maier mRegDC genes:", length(maier_mregdc_genes), "\n")

# 2. Compute both signatures
expr_mat <- as.matrix(expr_full)
storage.mode(expr_mat) <- "numeric"

gsva_param1 <- gsvaParam(expr_mat, list(mRegDC = mregdc_genes), minSize = 10, maxSize = 1000)
ssgsea1 <- gsva(gsva_param1, verbose = TRUE)

gsva_param2 <- gsvaParam(expr_mat, list(Maier_mRegDC = maier_mregdc_genes), minSize = 10, maxSize = 1000)
ssgsea2 <- gsva(gsva_param2, verbose = TRUE)

mregdc_df <- data.frame(sample_id = colnames(ssgsea1),
                        mRegDC_score = as.numeric(ssgsea1["mRegDC", ]))

maier_df <- data.frame(sample_id = colnames(ssgsea2),
                       Maier_mRegDC_score = as.numeric(ssgsea2["Maier_mRegDC", ]))

# 3. Merge and prepare data
meta_merged <- meta_rna %>%
  left_join(mregdc_df, by = c("sample.id" = "sample_id")) %>%
  left_join(maier_df, by = c("sample.id" = "sample_id")) %>%
  left_join(meta_subjects, by = "subject.id", suffix = c("", "_subj")) %>%
  filter(!is.na(mRegDC_score)) %>%
  filter(grepl("PD1|Combo", treatment.regimen.name, ignore.case = TRUE)) %>%
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
    CCR7_group = factor(mRegDC_high_low,
                        levels = c("High", "Low"),
                        labels = c("CCR7+DC High", "CCR7+DC Low")),
    Tumor_Type = sample.tumor.type,
    OS_time = os,
    PFS_time = pfs,
    OS_event = ifelse(!is.na(os), 1, NA),
    PFS_event = ifelse(!is.na(pfs), 1, NA)
  )

cat("FINAL DATASET:", nrow(meta_merged), "samples\n\n")
cat("Tumor subtype distribution:\n")
print(table(meta_merged$Tumor_Type))

# Refresh subsets
pre_data   <- meta_merged %>% filter(Timepoint_Group == "Pre-treatment")
onpost_data <- meta_merged %>% filter(Timepoint_Group == "On_Post-treatment")

cat("Sample sizes:\n")
cat("  Pre-treatment:     ", nrow(pre_data), "\n")
cat("  On/Post-treatment: ", nrow(onpost_data), "\n\n")

# 4. Analysis function (no number-at-risk CSVs)
run_analysis <- function(data, title) {
  if (nrow(data) < 10) return(NULL)
  
  cat("=== ", title, " (n =", nrow(data), ") ===\n")
  wilcox_p <- wilcox.test(mRegDC_score ~ response_binary, data = data)$p.value
  cat("Wilcoxon p (mRegDC) =", format.pval(wilcox_p, digits = 3), "\n")
  
  p_box <- ggplot(data, aes(x = response_binary, y = mRegDC_score, fill = response_binary)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.95, color = "#000000", linewidth = 0.95, fatten = 2.2) +
    geom_jitter(width = 0.3, size = 4.5, shape = 21, color = "black", fill = "white", stroke = 1.2) +
    stat_compare_means(method = "wilcox.test", label = "p.format", size = 5.5,
                       label.y = max(data$mRegDC_score, na.rm = TRUE) * 1.12) +
    scale_fill_manual(values = c("Responder" = "#e31a1c", "Non-responder" = "#000000")) +
    theme_minimal(base_size = 15) +
    labs(title = paste(signature_name, "by RECIST Response -", title),
         subtitle = paste("n =", nrow(data), "samples | All subtypes"),
         y = paste(signature_name, "ssGSEA Score"), x = "") +
    theme(legend.position = "none")
  
  print(p_box)
  ggsave(file.path(out_dir, paste0("CCR7_DC_Response_", gsub(" ", "_", title), ".pdf")),
         p_box, width = 7, height = 6.5, dpi = 300)
  
  data_surv <- data %>% filter(!is.na(OS_time) | !is.na(PFS_time)) %>%
    mutate(OS_event = ifelse(!is.na(OS_time), OS_event, 0),
           PFS_event = ifelse(!is.na(PFS_time), PFS_event, 0))
  
  if (nrow(data_surv) > 20) {
    fit_os  <- survfit(Surv(OS_time, OS_event) ~ CCR7_group, data = data_surv)
    fit_pfs <- survfit(Surv(PFS_time, PFS_event) ~ CCR7_group, data = data_surv)
    
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
    ggsave(file.path(out_dir, paste0("CCR7_DC_OS_KM_", gsub(" ", "_", title), ".pdf")),
           plot = p_os$plot, width = 8, height = 6, dpi = 300)
    ggsave(file.path(out_dir, paste0("CCR7_DC_PFS_KM_", gsub(" ", "_", title), ".pdf")),
           plot = p_pfs$plot, width = 8, height = 6, dpi = 300)
  }
}

# 5. Run the two analyses
run_analysis(pre_data,   "All subtypes - Pre-treatment")
run_analysis(onpost_data, "All subtypes - On_Post-treatment")

# 6. Correlation between the two mRegDC signatures
cat("\n=== Spearman correlation between mRegDC_score and Maier_mRegDC_score ===\n")
cor_pre <- cor.test(pre_data$mRegDC_score, pre_data$Maier_mRegDC_score, method = "spearman")
cor_onpost <- cor.test(onpost_data$mRegDC_score, onpost_data$Maier_mRegDC_score, method = "spearman")

cat("Pre-treatment:\n"); print(cor_pre)
cat("\nOn/Post-treatment:\n"); print(cor_onpost)

# 7. Export main summary table
summary_table_final <- meta_merged %>%
  select(sample.id, Timepoint_Group, treatment.regimen.name, bor, response,
         cohort, OS_time, PFS_time, OS_event, PFS_event,
         mRegDC_score, Maier_mRegDC_score, previous.treatment, Tumor_Type) %>%
  arrange(Timepoint_Group, bor) %>%
  rename(Sample_ID = sample.id, Timepoint = Timepoint_Group,
         Treatment = treatment.regimen.name, BOR = bor, Cohort = cohort,
         Prior_treatment = previous.treatment)

write.csv(summary_table_final, "mRegDC_Dual_mRegDC_Signatures_Summary_Table.csv", row.names = FALSE)
cat("\nMain summary table exported.\n")

# 8. Missing genes table
cat("\n=== Creating Supplementary Table for Missing mRegDC Genes ===\n")
missing_genes <- setdiff(mregdc_genes, rownames(expr_full))
write.csv(data.frame(Missing_Gene = missing_genes), 
          "Supplementary_Table_Missing_mRegDC_genes.csv", row.names = FALSE)
cat("Created with", length(missing_genes), "missing genes.\n")

cat("\n✅ Final clean script completed successfully!\n")
cat("All plots saved in:", out_dir, "\n")
sessionInfo()