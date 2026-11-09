# Mutengo et al. (2026) — *Schistosoma haematobium* in school-age children, Lusaka

Analysis script for:

> Mutengo MM, Mwansa J, Chisanga K, Zulu E, Mulunda NR, Muchinga J, Rodríguez E, Sánchez S, Castro L, Perteguer MJ, Carmena D, Sotillo J (2026). *Schistosoma haematobium* infection among school age children in two peri-urban areas of Lusaka, Zambia: A cross-sectional study. *PLoS Neglected Tropical Diseases* 20(10): e0014013. https://doi.org/10.1371/journal.pntd.0014013

## Contents

| Script | Description |
|---|---|
| `prevalence_ci.py` | Prevalence and 95% confidence interval using the normal (Wald) approximation. Run without arguments, it reproduces Table 2 of the article. |

The article links to an identical copy at the root of this repository, [`Prevalence_CI_Calculator.py`](../../Prevalence_CI_Calculator.py), which is kept at that path so the link remains valid.

## Usage

Requires Python 3 only, with no additional packages.

```bash
python prevalence_ci.py 48 208   # positives, total examined
python prevalence_ci.py          # reproduces Table 2
```

Output for a single estimate:

```
  48/208   23.08% (95% CI: 17.35-28.80)
```

## Method

For *x* positives out of *n* children examined, prevalence is *p* = *x*/*n* and the 95% confidence interval is *p* ± 1.96 √(*p*(1 − *p*)/*n*). The Wald interval is adequate for the sample sizes and prevalences reported here but performs poorly for very small samples or prevalences close to 0% or 100%; in those cases a Wilson or exact (Clopper–Pearson) interval is preferable.

## Data

Participant-level demographic and diagnostic data are provided as S1 Table of the article.

## Citation

If you use this script, please cite the article above.
