# BIOSCI 738 — Significance article project

iNaturalist citizen-science analysis on the **Biota of University of Auckland City Campus** project. Two assessment hooks:

- **IOA Project Discussion** — 10 min oral, 2026-05-01, 12.5%
- **Significance article (final)** — written piece + GitHub repo, 2026-06-05

## Working research question

> *How does observer effort heterogeneity bias estimates of urban biodiversity composition on the University of Auckland City Campus?*

Hits all four rubric criteria — and stephen_thorpe's 36% share *is* the story rather than a problem to apologise for.

## Layout

```
significance/
├── data/
│   ├── raw/observations_raw.json      # 310MB full API payloads
│   └── clean/observations_flat.csv    # 5,689 obs × 44 cols
├── scripts/
│   ├── 01_pull_inat.py                # cursor-paginated pull (regenerable)
│   └── 02_eda_snapshot.py             # EDA → notes/eda_snapshot.md
├── notes/
│   └── eda_snapshot.md                # dataset summary; observer Pareto, March bioblitz
├── figures/                           # ggplot output target
└── STATUS.md                          # IOA prep timeline
```

## Dataset shape (anchor numbers for the oral)

- 5,689 observations · 594 observers · 1,982 taxa
- Top observer (`stephen_thorpe`) = 36% of all obs; top 10 = 52%; 47% of observers are one-and-done
- 41% research-grade · 31% casual · 28% needs-id
- 53% introduced · 34% native · 18% endemic · 150 threatened
- Spatial: ~1.1 km × 1.1 km bbox; 8 m median positional accuracy
- Temporal: 2015 onward dominant; March bioblitz spike (1,351 vs ~400 in adjacent months)

## Charlotte-vocabulary checklist for the IOA

- Frame distribution choice from the response form first (counts → Poisson, then NB if dispersion fails)
- Open with the simulation: "each observer has a baseline rate λᵢ from a Gamma; obs are Poisson(λᵢ × effort)"
- Random effect for observer or location is non-negotiable here — nesting is extreme
- Limitations to name unprompted: observer-effort heterogeneity, path/building spatial bias, charismatic-species skew, research-grade vs casual ID quality, bioblitz temporal clumping
- Never say "significant"
