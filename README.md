# CCR7⁺ DC Signature Validation in Cutaneous + Acral Melanoma

**Orthogonal bulk RNA-seq validation for manuscript revision**

This repository contains the complete, reproducible analysis of the **CCR7⁺ DC signature** (multi-gene mRegDC module from Yang et al., *Nat Commun* 2025) in an independent harmonized bulk RNA-seq cohort of cutaneous and acral melanoma patients (Campbell et al., *Cancer Cell* 2023).

## Biological Context of the mRegDC / CCR7⁺ DC Signature

The **mRegDC** (mature DCs enriched in immunoregulatory molecules; also referred to as the CCR7⁺ DC or mature/migratory DC state) represents a conserved transcriptional program acquired by conventional dendritic cell (cDC) subsets in the tumor microenvironment.

Multiple peer-reviewed single-cell RNA-seq studies have established that this mature regulatory state is **expressed by / derived from both cDC1 and cDC2 lineages** upon tumor antigen uptake or exposure to TME signals:

- **Maier et al. (2020)** — Foundational study: "We find that the mregDC program is expressed by canonical DC1s and DC2s upon uptake of tumour antigens." (*Nature* 580, 257–262). DOI: [10.1038/s41586-020-2134-y](https://doi.org/10.1038/s41586-020-2134-y)

- **Gerhard et al. (2021)** — Demonstrated conservation of tumor-infiltrating DC states (including mregDC/DC3) across solid human cancers (*J Exp Med* 218(1): e20200264). DOI: [10.1084/jem.20200264](https://doi.org/10.1084/jem.20200264)

- **Yang et al. (2025)** — Trajectory and pathway analyses in metastatic melanoma: both cDC1 and cDC2 possess the potential to transition into mregDCs; mregDC likely derives from cDC1 with contributions from cDC2 (*Nat Commun* 16, 8151). DOI: [10.1038/s41467-025-62878-5](https://doi.org/10.1038/s41467-025-62878-5)

- Additional supporting evidence: mregDCs represent a matured state of both cDC1 and cDC2 subsets (Mazzoccoli et al., *Cancers* 2024; reviews in *FEBS Lett* 2025 and *Sci Immunol* 2022).

This state is characterized by high expression of maturation/migration markers (e.g., *CCR7*, *CD40*, *IL12B*) together with immunoregulatory molecules (e.g., *PD-L1*, *PD-L2*). The high *CCR7* expression confers migratory capacity, consistent with the "migratory DC" terminology used in some contexts. The present validation analysis tests the clinical relevance of this cDC1/cDC2-derived program in an independent melanoma cohort treated with PD-1-based regimens.

## Key Results
- Pre-treatment: Wilcoxon p = 0.0259 (n = 244)
- On/post-treatment: Wilcoxon p = 0.00337 (n = 59)
- Multivariable logistic regression: OR = 3.33 (p = 0.0239) pre-treatment; OR = 31.6 (p = 0.0063) on/post-treatment
- Strong correlation with estimated DC abundance (Spearman ρ = 0.858, p < 2.2e-16)
- 752 / 886 genes (84.9 %) of the original mRegDC module were detected

## Repository Contents
- `CCR7_DC_signature_analysis.R` – Main reproducible script
- `mRegDC_PD1_Focused_Summary_Table_CUTANEOUS_ACRAL.csv` – Patient-level results
- `mRegDC_Plots_Final_v4/` – All publication-ready PDF figures
- `Supplementary_Table_Missing_mRegDC_genes.csv` – List of 134 missing genes

## How to Reproduce
```bash
Rscript CCR7_DC_signature_analysis.R
```

> **Note**: The script requires the original expression matrix, metadata files, and the mRegDC gene list from Yang et al. (2025). These input files are not included in the repository due to data access restrictions but are described in the script comments.