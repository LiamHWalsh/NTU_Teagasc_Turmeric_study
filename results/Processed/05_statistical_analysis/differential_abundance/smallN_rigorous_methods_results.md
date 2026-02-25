# Differential Abundance Analysis (Rigorous Small-Sample Framework)

## Methods
Given the limited sample size and repeated-measures design, we used a pre-specified, small-sample analytical framework designed to maximise inferential validity while minimizing overfitting and multiplicity burden.

Species-level and genus-level profiles were analysed separately. Features were filtered to retain taxa present at non-zero abundance in at least 20% of retained samples. For feature-level inference, we analysed within-participant change from baseline to 6 months (delta = month6 - baseline) and compared delta distributions between turmeric-exposed and non-exposed participants using exact Wilcoxon rank-sum tests. This approach preserves paired structure within participants and avoids large-sample distributional assumptions that are unstable at low n.

In parallel, exact permutation p-values were computed by exhaustive label permutation under fixed group sizes. Effect magnitude was reported as Cohen's d and rank-biserial correlation. Uncertainty in mean between-group delta differences was quantified with bootstrap 95% confidence intervals (10,000 resamples).

Multiple testing was controlled within each pre-specified panel (species, genus) using Benjamini-Hochberg false discovery rate (FDR) adjustment. We defined confirmatory evidence as q < 0.05 and exploratory evidence as q < 0.25. To assess whether treatment-associated structure was detectable at the multivariate compositional level, we additionally applied an exact label-permutation test to participant-level clr change vectors (reporting pseudo-F and R2).

Robustness of feature-level findings was evaluated using leave-one-participant-out (LOO) re-analysis across all panel features, summarised by median p-value, p-value range, fraction of LOO fits with p < 0.05, and sign consistency of direction estimates.

## Results
The paired analysis included 11 baseline-to-6-month participant pairs (8 turmeric-exposed, 3 non-exposed) for both species and genus panels.

Species panel: 12 features were tested. No species met confirmatory FDR significance (q < 0.05). Two species met exploratory evidence criteria (q < 0.25): *Blautia wexlerae* and *Roseburia intestinalis* (both q = 0.0727).

Genus panel: 12 features were tested. No genus met confirmatory FDR significance (q < 0.05). Five genera met exploratory evidence criteria (q < 0.25): *Blautia*, *Alistipes* (both q = 0.0727), *Ruminococcus*, *Dorea* (both q = 0.1455), and unclassified *Lachnospiraceae* (q = 0.2036).

At the global compositional level, exact permutation testing on clr change vectors did not show statistically significant treatment-associated structure: species pseudo-F = 1.133, R2 = 0.1118, p = 0.259 (165 permutations); genus pseudo-F = 0.915, R2 = 0.0923, p = 0.627 (165 permutations).

LOO diagnostics indicated that the strongest exploratory signals (e.g., *Blautia wexlerae*, *Blautia*, *Alistipes*) were directionally stable, but this stability should be interpreted cautiously because small group sizes make individual observations highly influential.

Overall, the data support directional, biologically plausible hypotheses rather than definitive differential-abundance claims. Under this low-n setting, effect sizes and stability metrics are more informative than dichotomous significance calls, and findings should be treated as prioritized candidates for independent replication.

## Figure
`smallN_rigorous_species_genus_summary.png` summarises feature-level effect sizes and adjusted FDR evidence across species and genus panels.
