# STATUS — Significance project

**Last updated:** 2026-04-28

## Timeline to IOA — DONE

**IOA delivered 2026-04-28.** Final article + repo still due 2026-06-05.

## Findings from 04-28 EDA (talking points for the IOA)

**Lifecycle / funnel** (notes/lifecycle_eda.md):
- After upload: 100% → 99.6% has any ID → 67% has any external ID → 41% research grade.
- Lag distributions are **heavy-tailed at every stage** (median in hours, p99 in years). Not exponential — assume lognormal/Weibull if modelling time-to-event.
- Funnel-row caveat: 2 *external* supports ≠ RG threshold because iNat counts the observer's own ID toward the 2/3 rule. ~940 obs reach RG via observer + 1 external.

**Identifier roles** (notes/group_trips_hypothesis.md, notes/lloyd_and_bioblitz.md):
- Three distinct roles: observer / leading-improver / supporter.
- **lloyd_esler** is the top external IDer (438 IDs, 196 distinct observers helped) but **0pp RG effect** — he's 78% supporting, mostly redundant 3rd/4th agreement on already-converging obs. He's the Invercargill mass-IDer, not local.
- **stephen_thorpe** has only 128 external IDs but **+18pp RG effect** (39%→57%) — 30% improving / 22% leading / 46% supporting. He establishes/changes community taxon, attracts agreement.
- Failed hypothesis: I expected Lloyd's no-effect was IDing at family/genus. Wrong — he IDs at species 87% of the time. Real answer is the supporting-vs-improving role split.

**Bioblitzes are course assignments** (notes/lloyd_and_bioblitz.md):
- All 10 high-observer-count days (≥10 distinct observers) cluster in **Feb/March on Tu/Wed/Thu** — 1 weekend day across the whole dataset.
- **86.9% of bioblitz observers (152/175) are one-and-done** — never observed before or since.
- Cohort flips ~100% each year (only 2/63 returning in 2026).
- Bioblitz obs are 33% RG vs 43% non-bioblitz; **Fungi take a -33pp hit**, Insecta -14pp, Aves +6pp.
- 2026-03-24 alone is 798 obs (14% of the dataset!) with 53% needs_id.

## Open decisions

- Confirm research question vs. two alternatives (not yet pitched)
- Decide observer effort model: random effect `(1|user_login)` vs offset `offset(log(n_obs_by_user))`
- **NEW:** consider adding `(1|bioblitz_event_date)` as a third level — it absorbs the course-cohort effect that `(1|user_login)` can't (since blitz observers each have only 1 obs)
- Whether to add spatial term (Matérn / CAR) — likely defer, mention as alternative

## Known limitations to surface in IOA

1. Observer-effort heterogeneity (stephen_thorpe = 36%)
2. Spatial bias toward paths/buildings
3. Research-grade (41%) vs casual (31%) ID quality — and **identifier-role asymmetry**: who IDs you matters more than how many ID you
4. Charismatic-species taxonomic skew (Aves 91% RG, Plantae 22%, Fungi worst-affected by bioblitzes)
5. March bioblitz temporal clumping — **specifically a recurring course assignment with single-event tourist observers**, not opportunistic citizen science
6. Right-censoring: 2026 cohort is 47% still needs_id and may yet flip
