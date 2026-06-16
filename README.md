# CCR7⁺ DC Signature Validation in Melanoma (Bulk RNA-seq)

**Dual independent mRegDC signatures – Final validated analysis**

This repository contains the complete, reproducible analysis of **two independent mRegDC signatures** in an independent harmonized bulk RNA-seq cohort of melanoma patients treated with anti-PD-1 ± anti-CTLA-4 (Campbell et al., Cancer Cell 2023).

### Key Results
- **Primary mRegDC signature** (886 genes) and **Maier et al. (Nature 2020) mRegDC signature** (294 genes)
- Strong concordance between signatures:
  - Pre-treatment: Spearman ρ = **0.87** (p < 2.2e-16)
  - On/Post-treatment: Spearman ρ = **0.92** (p < 2.2e-16)
- Significant association with RECIST response in both pre-treatment (p = 0.000619) and on/post-treatment (p = 0.000159) samples
- Full stratification by treatment timing with Kaplan–Meier OS and PFS analyses

### Repository Contents
- `CCR7_DC_signature_analysis.R` — Main reproducible R script
- `mRegDC_Dual_mRegDC_Signatures_Summary_Table.csv` — Patient-level results
- `mRegDC_Plots_Final_v4/` — Publication-ready PDF figures
- `Supplementary_Table_Missing_mRegDC_genes.csv` — List of genes not detected

### How to Reproduce
```bash
Rscript CCR7_DC_signature_analysis.R
