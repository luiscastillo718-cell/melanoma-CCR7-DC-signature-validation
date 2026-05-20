# CCR7⁺ DC Signature Validation in Cutaneous + Acral Melanoma

**Orthogonal bulk RNA-seq validation for manuscript revision**

This repository contains the complete, reproducible analysis of the **CCR7⁺ DC signature** (multi-gene mRegDC module from Yang et al., Nat Commun 2025) in an independent harmonized bulk RNA-seq cohort of cutaneous and acral melanoma patients (Campbell et al., Cancer Cell 2023).

### Key Results
- Pre-treatment: Wilcoxon p = 0.0259 (n = 244)
- On/post-treatment: Wilcoxon p = 0.00337 (n = 59)
- Multivariable logistic regression: OR = 3.33 (p = 0.0239) pre-treatment; OR = 31.6 (p = 0.0063) on/post-treatment
- Strong correlation with estimated DC abundance (Spearman ρ = 0.858, p < 2.2e-16)
- 752 / 886 genes (84.9 %) of the original mRegDC module were detected

### Repository Contents
- `CCR7_DC_signature_analysis.R` – Main reproducible script
- `mRegDC_PD1_Focused_Summary_Table_CUTANEOUS_ACRAL.csv` – Patient-level results
- `mRegDC_Plots_Final_v4/` – All publication-ready PDF figures
- `Supplementary_Table_Missing_mRegDC_genes.csv` – List of 134 missing genes

### How to Reproduce
```bash
Rscript CCR7_DC_signature_analysis.R