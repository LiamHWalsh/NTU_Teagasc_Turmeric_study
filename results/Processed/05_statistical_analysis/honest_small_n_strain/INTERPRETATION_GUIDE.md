# Honest Small-Sample Microbiome and Functional Analysis

## Study Design
- Retained participants: 16
- Paired baseline->6 month participants: 11
- Treatment pairs: 8, Control pairs: 3

## Statistical Position
- This is an effect-estimation analysis for a low-N study.
- No feature-level significance threshold claims are made.
- Uncertainty and sensitivity diagnostics are reported explicitly.

## Transformation and Compositional Strategy
- Primary effect scale: raw (for interpretability).
- Secondary transform sensitivity: log1p (directional concordance reported).
- Taxonomy CLR sensitivity enabled: TRUE.

## Methods (Concise)
- Taxonomy filtering: prevalence >= 0.2, mean abundance >= 0.001
- Pathway filtering: prevalence >= 0.2, mean abundance >= 1
- Paired delta per participant: value_6m - value_baseline.
- Effect metrics: Cohen's d, probability of superiority, 95% bootstrap CI (2000 resamples).
- LOO summary: sign agreement with full-sample effect (participant-sensitivity, not validation).
- Feature ranking: transparent multi-rank (effect size, CI width, LOO sign agreement).

### Strain Panel
- Features analyzed: 29
- Median |d|: 0.71
- Median CI width: 0.786
- Median LOO sign agreement: 1
- Raw-vs-log1p sign agreement: 1
- Top features (effect size + precision + LOO sensitivity):
  - SGB1934: d=-1.4, CI [-0.59, -0.08], width=0.5, LOO sign agreement=1, d_raw=-1.4, d_log1p=-1.41, d_clr=-0.86, primary-vs-CLR sign agreement=TRUE
  - SGB2290: d=-1.04, CI [-0.4, -0.01], width=0.39, LOO sign agreement=1, d_raw=-1.04, d_log1p=-1.03, d_clr=-0.68, primary-vs-CLR sign agreement=TRUE
  - SGB2328: d=-0.85, CI [-0.31, -0.01], width=0.3, LOO sign agreement=1, d_raw=-0.85, d_log1p=-0.87, d_clr=-0.9, primary-vs-CLR sign agreement=TRUE
  - SGB4804: d=1.07, CI [0.01, 0.5], width=0.49, LOO sign agreement=1, d_raw=1.07, d_log1p=1.08, d_clr=1.07, primary-vs-CLR sign agreement=TRUE
  - SGB4828: d=1.12, CI [0.07, 0.87], width=0.8, LOO sign agreement=1, d_raw=1.12, d_log1p=1.28, d_clr=1.29, primary-vs-CLR sign agreement=TRUE
  - SGB4816: d=-0.82, CI [-0.45, -0.01], width=0.45, LOO sign agreement=1, d_raw=-0.82, d_log1p=-0.84, d_clr=-0.41, primary-vs-CLR sign agreement=TRUE

### Genus Panel
- Features analyzed: 5
- Median |d|: 0.83
- Median CI width: 9.261
- Median LOO sign agreement: 1
- Raw-vs-log1p sign agreement: 1
- Top features (effect size + precision + LOO sensitivity):
  - Blautia: d=1.92, CI [5.34, 14.6], width=9.26, LOO sign agreement=1, d_raw=1.92, d_log1p=1.64, d_clr=1.3, primary-vs-CLR sign agreement=TRUE
  - Alistipes: d=-1.76, CI [-8.06, -2.57], width=5.49, LOO sign agreement=1, d_raw=-1.76, d_log1p=-1.96, d_clr=-1.35, primary-vs-CLR sign agreement=TRUE
  - Parabacteroides: d=-0.75, CI [-0.88, 0.19], width=1.07, LOO sign agreement=1, d_raw=-0.75, d_log1p=-0.76, d_clr=-0.61, primary-vs-CLR sign agreement=TRUE
  - Bacteroides: d=-0.83, CI [-8.37, 1.07], width=9.44, LOO sign agreement=1, d_raw=-0.83, d_log1p=-0.78, d_clr=-0.42, primary-vs-CLR sign agreement=TRUE
  - Bifidobacterium: d=-0.83, CI [-17.01, 0.38], width=17.39, LOO sign agreement=1, d_raw=-0.83, d_log1p=-0.95, d_clr=0.03, primary-vs-CLR sign agreement=FALSE

### Pathway Panel
- Features analyzed: 339
- Median |d|: 0.67
- Median CI width: 106.452
- Median LOO sign agreement: 1
- Raw-vs-log1p sign agreement: 0.94
- Top features (effect size + precision + LOO sensitivity):
  - PWY66-389: phytol degradation: d=-1.34, CI [-5.11, -0.43], width=4.68, LOO sign agreement=1, d_raw=-1.34, d_log1p=-1.59
  - REDCITCYC: TCA cycle VI (Helicobacter): d=-1.08, CI [-6.4, -0.42], width=5.98, LOO sign agreement=1, d_raw=-1.08, d_log1p=-1.34
  - PWY-7254: TCA cycle VII (acetate-producers): d=-1.11, CI [-8.43, -1.11], width=7.31, LOO sign agreement=1, d_raw=-1.11, d_log1p=-1.23
  - SO4ASSIM-PWY: assimilatory sulfate reduction I: d=-1.23, CI [-13.43, -1.64], width=11.79, LOO sign agreement=1, d_raw=-1.23, d_log1p=-1.85
  - PWY-7209: superpathway of pyrimidine ribonucleosides degradation: d=2.13, CI [20.9, 58.56], width=37.66, LOO sign agreement=1, d_raw=2.13, d_log1p=0.71
  - GOLPDLCAT-PWY: superpathway of glycerol degradation to 1,3-propanediol: d=2.35, CI [36.42, 85.98], width=49.55, LOO sign agreement=1, d_raw=2.35, d_log1p=1.05

## Taxonomy-Function Integration (Exploratory)
- Delta-based Spearman correlations with bootstrap CI and permutation p-values are provided.
- BH-adjusted q-values are included as multiplicity context across tested pairs.
- Low N means these are hypothesis-prioritization signals, not confirmatory findings.
- Top candidate pairs:
  - Blautia -> superpathway of glycerol degradation to 1,3-propanediol: rho=0.92, CI [0.65, 1], p_perm=5e-04, q_perm=0.025, n=11, rho_tx=0.81, rho_ct=NA
  - Alistipes -> guanosine nucleotides degradation III: rho=-0.83, CI [-0.97, -0.41], p_perm=0.004, q_perm=0.08, n=11, rho_tx=-0.57, rho_ct=NA
  - Alistipes -> 2-oxobutanoate degradation I: rho=0.89, CI [0.56, 1], p_perm=0.001, q_perm=0.0333, n=11, rho_tx=0.93, rho_ct=NA
  - Alistipes -> TCA cycle V (2-oxoglutarate synthase): rho=0.82, CI [0.3, 0.99], p_perm=0.003, q_perm=0.075, n=11, rho_tx=0.83, rho_ct=NA
  - Blautia -> superpathway of pyrimidine ribonucleosides degradation: rho=0.76, CI [0.29, 0.97], p_perm=0.01, q_perm=0.0923, n=11, rho_tx=0.45, rho_ct=NA
  - Alistipes -> superpathway of glycerol degradation to 1,3-propanediol: rho=-0.74, CI [-0.92, -0.24], p_perm=0.012, q_perm=0.0923, n=11, rho_tx=-0.38, rho_ct=NA
  - Bacteroides -> chondroitin sulfate degradation I (bacterial): rho=0.91, CI [0.54, 1], p_perm=5e-04, q_perm=0.025, n=11, rho_tx=0.88, rho_ct=NA
  - Alistipes -> 1,4-dihydroxy-6-naphthoate biosynthesis II: rho=0.78, CI [0.39, 0.92], p_perm=0.008, q_perm=0.0888, n=11, rho_tx=0.76, rho_ct=NA
  - Blautia -> guanosine nucleotides degradation III: rho=0.72, CI [0.12, 0.93], p_perm=0.0175, q_perm=0.125, n=11, rho_tx=0.33, rho_ct=NA
  - Bifidobacterium -> D-glucarate degradation I: rho=-0.8, CI [-0.96, -0.4], p_perm=0.0075, q_perm=0.0888, n=11, rho_tx=-0.84, rho_ct=NA

## What This Supports
- Ranked hypotheses for follow-up based on effect size, precision, and participant-level sensitivity.
- Directional signal stability checks across transform and compositional sensitivity analyses.
- Supplementary Figures S1-S4 map to complete panels: S1 (all strain-resolved features), S2 (all genera), S3 (all pathways), S4 (all tested genus-pathway pairs).

## What This Does Not Support
- Confirmatory causal claims from this dataset alone.
- Population-level generalization without independent replication.

## Recommended Next Step
- Replicate top-ranked taxonomy and pathway candidates in a larger independent cohort.

Report generated: 2026-03-05
