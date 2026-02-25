# Differential Abundance Analysis (Low-N Confirmatory + Exploratory)

## Methods
A low-sample confirmatory framework was applied to species-level and genus-level abundance profiles.
First, candidate features were prespecified using biologically relevant taxa and turmeric-associated candidates from MaAsLin3 screening results.
For each feature, within-participant baseline-to-6-month change was calculated and compared between turmeric treatment and control groups.
The primary endpoint used an exact Wilcoxon rank-sum test on change scores (Treatment delta vs Control delta), with permutation testing as sensitivity analysis.
Effect sizes were summarized as Cohen's d on change scores. Multiple testing correction across each prespecified panel used Benjamini-Hochberg FDR.
Confirmatory significance was defined as q < 0.05, and exploratory evidence as q < 0.25.

## Results
Species panel: 12 features analyzed; 0 with q < 0.05 and 2 with q < 0.25.
Genus panel: 12 features analyzed; 0 with q < 0.05 and 3 with q < 0.25.
Top-ranked species-level features (by q then p): Blautia_wexlerae, Dorea_formicigenerans, Eubacterium_rectale.
Top-ranked genus-level features (by q then p): Blautia, Alistipes, Dorea.
Across both levels, nominal signals were present but no robust FDR-confirmed findings (q < 0.05) were observed, consistent with limited power in a small repeated-measures cohort.

## Figure
See `lowN_confirmatory_species_genus_summary.png` for a combined species/genus effect-size summary plot.
