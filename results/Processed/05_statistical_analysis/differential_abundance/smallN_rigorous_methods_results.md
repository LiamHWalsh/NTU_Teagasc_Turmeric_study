# Differential Abundance Analysis (Rigorous Small-Sample Framework)

## Methods
A pre-specified small-sample framework was used to limit multiplicity, preserve paired structure, and provide exact inference where feasible.
Species-level and genus-level abundance profiles were analysed separately after prevalence filtering (non-zero abundance in >=20% of retained samples).
Primary feature-level inference used within-participant baseline-to-6-month change (delta = month6 - baseline) and compared treatment versus control deltas by exact Wilcoxon rank-sum tests.
Permutation p-values were computed by exhaustive label permutation of group assignments under fixed group sizes.
Effect size was summarized using Cohen's d and rank-biserial correlation. Uncertainty in mean change difference was estimated by bootstrap 95% confidence intervals (10,000 resamples).
False discovery control used Benjamini-Hochberg correction within each pre-specified panel. Confirmatory evidence was defined as q < 0.05 and exploratory evidence as q < 0.25.
To test overall compositional shifts beyond individual features, clr-transformed profiles were analysed using an exact label-permutation test on participant-level clr change vectors (pseudo-F and R2 reported).
Stability of feature-level findings was assessed by leave-one-participant-out re-analysis across all panel features.

## Results
Species panel: 12 features analysed; 0 with q < 0.05 and 2 with q < 0.25.
Genus panel: 12 features analysed; 0 with q < 0.05 and 5 with q < 0.25.
Top species-level features (q then p): Blautia_wexlerae, Roseburia_intestinalis, Eubacterium_rectale.
Top genus-level features (q then p): Blautia, Alistipes, Ruminococcus.
Global clr exact permutation test: Species p = 0.259, R2 = 0.112; Genus p = 0.627, R2 = 0.0923.
Overall, directional shifts were observed in biologically plausible taxa, but no feature reached strict confirmatory FDR significance (q < 0.05).
Given sample size constraints, findings should be interpreted as prioritized hypotheses with quantified uncertainty and stability metrics rather than definitive treatment effects.

## Figure
See `smallN_rigorous_species_genus_summary.png` for combined species/genus effect-size and FDR evidence summary.
