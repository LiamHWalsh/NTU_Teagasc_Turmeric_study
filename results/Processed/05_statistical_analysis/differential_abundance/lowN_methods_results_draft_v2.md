# Differential Abundance Analysis (Low-N Confirmatory and Exploratory Framework)

## Methods
Given the limited repeated-measures sample size, a pre-specified low-N analysis framework was used to prioritize estimation stability and control of false discoveries. Species-level and genus-level analyses were run separately on filtered relative-abundance tables.

Candidate features were selected a priori from two sources: (i) biologically relevant taxa linked to short-chain fatty acid production or prior turmeric/microbiome hypotheses, and (ii) turmeric-associated features from MaAsLin3 screening output. For each selected feature, within-participant change was computed as:

\[
\Delta = \text{abundance}_{6\text{ months}} - \text{abundance}_{baseline}
\]

The primary endpoint compared change scores between treatment and control groups (Treatment delta vs Control delta) using the exact Wilcoxon rank-sum test. A permutation Wilcoxon test (10,000 resamples) was used as a sensitivity analysis. Effect size was summarized as Cohen’s d on change scores.

To reduce sparse-feature instability, features were required to meet a minimum prevalence threshold (non-zero abundance in at least 20% of retained samples). Multiple testing correction within each pre-specified panel used the Benjamini-Hochberg false discovery rate (FDR). Confirmatory significance was defined as q < 0.05; q < 0.25 was treated as exploratory only.

## Results
The paired analysis included 11 participants with both baseline and 6-month samples (treatment n = 8, control n = 3).

At species level (12-feature confirmatory panel), 2 features reached exploratory FDR support (q < 0.25), but none met confirmatory FDR significance (q < 0.05). The strongest species-level signals were:

- Blautia wexlerae: p = 0.012, q = 0.145, d = 2.45
- Dorea formicigenerans: p = 0.024, q = 0.145, d = 1.89
- Eubacterium rectale: p = 0.085, q = 0.339, d = 1.47

At genus level (12-feature confirmatory panel), 3 genera met exploratory FDR support (q < 0.25), and none met q < 0.05:

- Blautia: p = 0.012, q = 0.073, d = 1.93 (positive treatment-associated change)
- Alistipes: p = 0.012, q = 0.073, d = -1.77 (negative treatment-associated change)
- Dorea: p = 0.048, q = 0.194, d = 1.38 (positive treatment-associated change)

Overall, the directionality of exploratory signals is biologically coherent (increases in selected Lachnospiraceae-related taxa and decreases in Alistipes/Bacteroidetes-linked signal), but inferential certainty remains limited. No association survived strict multiple-testing correction at q < 0.05, and all feature-level claims should therefore be interpreted as hypothesis-generating.

## Figure
The file `lowN_confirmatory_species_genus_summary.png` provides a combined species/genus effect-size summary. Points are ranked by FDR-adjusted evidence and show effect magnitude (Cohen’s d), nominal significance (-log10 p), and FDR evidence bands (q < 0.05, q < 0.25, q >= 0.25).
