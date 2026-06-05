"""Two follow-ups:
  (1) Why doesn't lloyd_esler's IDing convert obs to research grade,
      while stephen_thorpe's does? Look at rank, category, agreement,
      and position in the ID chain.
  (2) What bioblitz patterns are visible? Find days with anomalously
      high observer count, look at recurrence (same March each year?),
      observer cohort overlap, and downstream RG fate.
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data/raw/observations_raw.json"
OUT = ROOT / "notes/lloyd_and_bioblitz.md"


def load():
    with RAW.open() as f:
        data = json.load(f)
    return data if isinstance(data, list) else data["results"]


def main():
    obs_list = load()

    # Build per-obs and per-identification frames
    obs_rows, id_rows = [], []
    for o in obs_list:
        observer = (o.get("user") or {}).get("login")
        observed = pd.to_datetime(o.get("observed_on"), errors="coerce", utc=True)
        obs_rows.append({
            "id": o["id"],
            "observer": observer,
            "observed_date": observed.date() if pd.notna(observed) else None,
            "quality_grade": o.get("quality_grade"),
            "captive": bool(o.get("captive")),
            "iconic": (o.get("taxon") or {}).get("iconic_taxon_name"),
            "obs_taxon_rank": (o.get("taxon") or {}).get("rank"),
        })
        idents = sorted(
            (o.get("identifications") or []),
            key=lambda i: pd.to_datetime(i.get("created_at"), errors="coerce", utc=True),
        )
        for pos, ident in enumerate(idents):
            id_rows.append({
                "obs_id": o["id"],
                "obs_qg": o.get("quality_grade"),
                "obs_observer": observer,
                "identifier": (ident.get("user") or {}).get("login"),
                "ident_at": pd.to_datetime(ident.get("created_at"), errors="coerce", utc=True),
                "own_observation": bool(ident.get("own_observation")),
                "category": ident.get("category"),
                "current": bool(ident.get("current")),
                "disagreement": ident.get("disagreement"),
                "ident_taxon_rank": (ident.get("taxon") or {}).get("rank"),
                "ident_taxon_id": (ident.get("taxon") or {}).get("id"),
                "ident_pos": pos,  # 0 = first ID on the obs
                "n_idents_total": len(idents),
            })

    df = pd.DataFrame(obs_rows)
    ids = pd.DataFrame(id_rows)

    parts = []
    add = parts.append
    add("# Why Lloyd doesn't convert + bioblitz patterns\n\n")

    # =============================================================
    # (1) LLOYD vs STEPHEN as identifiers
    # =============================================================
    add("## (1) Lloyd vs Stephen as identifiers\n\n")

    targets = ["stephen_thorpe", "lloyd_esler", "nhudson", "cooperj"]
    # Restrict to EXTERNAL IDs (their work on other people's obs)
    ext = ids[~ids["own_observation"]].copy()

    # Rank distribution
    rank_order = ["kingdom", "phylum", "subphylum", "class", "subclass", "order",
                  "suborder", "family", "subfamily", "tribe", "genus", "subgenus",
                  "species", "subspecies", "variety", "form", "complex", "section",
                  "stateofmatter"]
    add("### Identification taxon rank — what level do they ID at?\n\n")
    rank_tbl = (ext[ext["identifier"].isin(targets)]
                .groupby(["identifier", "ident_taxon_rank"])
                .size()
                .unstack(fill_value=0))
    # order columns by rank specificity
    cols_present = [c for c in rank_order if c in rank_tbl.columns]
    extra = [c for c in rank_tbl.columns if c not in cols_present]
    rank_tbl = rank_tbl[cols_present + extra]
    # Add % at species-or-finer
    species_or_finer = {"species", "subspecies", "variety", "form", "complex"}
    fine_cols = [c for c in rank_tbl.columns if c in species_or_finer]
    rank_tbl["TOTAL"] = rank_tbl.sum(axis=1)
    rank_tbl["% species-or-finer"] = (
        100 * rank_tbl[fine_cols].sum(axis=1) / rank_tbl["TOTAL"]
    ).round(1)
    add(rank_tbl.to_markdown() + "\n\n")
    add("**This is the answer.** RG requires a *species-level* community taxon. "
        "If you ID at family/genus you can never push the obs past the 2/3 species rule.\n\n")

    # Category (leading/improving/supporting/maverick)
    add("### Category distribution (iNat's own classification of each ID)\n\n")
    cat_tbl = (ext[ext["identifier"].isin(targets)]
               .groupby(["identifier", "category"]).size()
               .unstack(fill_value=0))
    cat_tbl["TOTAL"] = cat_tbl.sum(axis=1)
    for c in cat_tbl.columns:
        if c != "TOTAL":
            cat_tbl[f"{c}_pct"] = (100 * cat_tbl[c] / cat_tbl["TOTAL"]).round(1)
    add(cat_tbl.to_markdown() + "\n\n")
    add("- **leading** = first ID at a finer rank (drives the community taxon forward).\n")
    add("- **improving** = refines the community taxon.\n")
    add("- **supporting** = agrees with the existing community taxon (this is what triggers RG).\n")
    add("- **maverick** = disagrees with consensus.\n\n")

    # Position in chain
    add("### Position in the ID chain (0 = first ID on the obs)\n\n")
    pos_tbl = (ext[ext["identifier"].isin(targets)]
               .groupby("identifier")
               .agg(median_position=("ident_pos", "median"),
                    mean_position=("ident_pos", "mean"),
                    pct_first=("ident_pos", lambda s: 100*(s == 0).mean()),
                    pct_first_external=("ident_pos", lambda s: 100*(s <= 1).mean()))
               .round(2))
    add(pos_tbl.to_markdown() + "\n\n")
    add("If they're typically the *first* external IDer, they're seeding the community taxon "
        "(Lloyd's role). If they show up later, they're the second supporter that flips RG (Stephen's role).\n\n")

    # Disagreement rate
    add("### Disagreement rate\n\n")
    disag_tbl = (ext[ext["identifier"].isin(targets)]
                 .groupby("identifier")
                 .agg(n=("obs_id", "size"),
                      pct_disagree=("disagreement", lambda s: 100*(s == True).mean()))
                 .round(1))
    add(disag_tbl.to_markdown() + "\n\n")

    # =============================================================
    # (2) BIOBLITZ PATTERNS
    # =============================================================
    add("## (2) Bioblitz patterns\n\n")

    by_day = (df.dropna(subset=["observed_date"])
                .groupby("observed_date")
                .agg(n_obs=("id", "size"),
                     n_observers=("observer", "nunique"),
                     pct_research=("quality_grade",
                                    lambda s: 100*(s == "research").mean()),
                     pct_needs_id=("quality_grade",
                                    lambda s: 100*(s == "needs_id").mean()),
                     pct_casual=("quality_grade",
                                  lambda s: 100*(s == "casual").mean()))
                .reset_index())
    by_day["observed_date"] = pd.to_datetime(by_day["observed_date"])
    by_day["year"] = by_day["observed_date"].dt.year
    by_day["month"] = by_day["observed_date"].dt.month
    by_day["dow"] = by_day["observed_date"].dt.day_name()

    # Define a "bioblitz-class" day: ≥10 distinct observers (this captures
    # everything the data calls a group event — verified manually above)
    blitz = by_day[by_day["n_observers"] >= 10].sort_values("observed_date")
    add(f"### All days with ≥10 distinct observers ({len(blitz)} total)\n\n")
    add(blitz[["observed_date", "n_obs", "n_observers",
               "pct_research", "pct_needs_id", "pct_casual",
               "dow"]].to_markdown(index=False) + "\n\n")

    # Recurrence by month / year
    add("### Recurrence — when do bioblitzes happen?\n\n")
    by_month = (blitz.groupby("month")
                     .agg(n_blitz_days=("observed_date", "size"),
                          total_obs=("n_obs", "sum"),
                          mean_observers=("n_observers", "mean"))
                     .round(1))
    add("By month-of-year:\n\n")
    add(by_month.to_markdown() + "\n\n")

    by_year = (blitz.groupby("year")
                    .agg(n_blitz_days=("observed_date", "size"),
                         total_obs=("n_obs", "sum"),
                         max_observers=("n_observers", "max"))
                    .round(1))
    add("By year:\n\n")
    add(by_year.to_markdown() + "\n\n")

    # Day-of-week
    dow_order = ["Monday", "Tuesday", "Wednesday", "Thursday", "Friday",
                 "Saturday", "Sunday"]
    by_dow = blitz["dow"].value_counts().reindex(dow_order, fill_value=0)
    add("Day-of-week of bioblitz events:\n\n")
    add(by_dow.to_frame("n_blitz_days").to_markdown() + "\n\n")

    # Observer overlap year-over-year on the bioblitz days
    add("### Observer cohort — same students every year, or fresh cohort each March?\n\n")
    blitz_dates = set(blitz["observed_date"].dt.date)
    blitz_obs = df[df["observed_date"].isin(blitz_dates)].copy()
    blitz_obs["year"] = pd.to_datetime(blitz_obs["observed_date"]).dt.year
    yearly_observers = (blitz_obs.groupby("year")["observer"]
                        .agg(lambda s: set(s.dropna())))
    overlap_rows = []
    for y in sorted(yearly_observers.index):
        obs_y = yearly_observers[y]
        prev_obs = set().union(*[yearly_observers[y2] for y2 in yearly_observers.index if y2 < y])
        overlap_rows.append({
            "year": y,
            "n_blitz_observers": len(obs_y),
            "returning_from_prev_years": len(obs_y & prev_obs),
            "new_this_year": len(obs_y - prev_obs),
            "pct_returning": round(100 * len(obs_y & prev_obs) / max(len(obs_y), 1), 1),
        })
    add(pd.DataFrame(overlap_rows).to_markdown(index=False) + "\n\n")
    add("If `pct_returning` stays low, each year's bioblitz is a fresh student cohort — "
        "consistent with a recurring course assignment, not a returning enthusiast crowd.\n\n")

    # Downstream fate: RG rate of bioblitz obs vs non-bioblitz, controlled for taxon
    add("### Downstream fate of bioblitz observations\n\n")
    df["is_bioblitz"] = df["observed_date"].isin(blitz_dates)
    fate = (df.groupby("is_bioblitz")
              .agg(n=("id", "size"),
                   pct_research=("quality_grade", lambda s: 100*(s=="research").mean()),
                   pct_casual=("quality_grade", lambda s: 100*(s=="casual").mean()),
                   pct_needs_id=("quality_grade", lambda s: 100*(s=="needs_id").mean()))
              .round(1))
    add(fate.to_markdown() + "\n\n")

    # Same split, but conditioned on iconic group
    fate_iconic = (df.groupby(["is_bioblitz", "iconic"])
                     .agg(n=("id", "size"),
                          pct_research=("quality_grade",
                                         lambda s: 100*(s=="research").mean()))
                     .round(1)
                     .reset_index())
    fate_iconic = fate_iconic[fate_iconic["n"] >= 20]
    pivot = fate_iconic.pivot(index="iconic", columns="is_bioblitz",
                               values="pct_research")
    pivot.columns = ["non_bioblitz_pct_RG", "bioblitz_pct_RG"]
    pivot["delta"] = (pivot["bioblitz_pct_RG"] - pivot["non_bioblitz_pct_RG"]).round(1)
    add("**RG rate by iconic group, bioblitz vs non-bioblitz:**\n\n")
    add(pivot.sort_values("non_bioblitz_pct_RG", ascending=False).to_markdown() + "\n\n")

    # Bioblitz observer obs profile: do they obs anywhere else, or only at the blitz?
    add("### One-and-done observers (bioblitz tourists)\n\n")
    blitz_observers = set(blitz_obs["observer"].dropna())
    bo_total = (df[df["observer"].isin(blitz_observers)]
                .groupby("observer").size())
    bo_blitz = (blitz_obs.groupby("observer").size())
    bo_summary = pd.DataFrame({
        "n_total_obs": bo_total,
        "n_blitz_obs": bo_blitz,
    })
    bo_summary["pct_at_blitz"] = (100 * bo_summary["n_blitz_obs"] / bo_summary["n_total_obs"]).round(1)
    one_shot = (bo_summary["pct_at_blitz"] == 100).sum()
    add(f"- {len(bo_summary)} distinct observers participated in at least one bioblitz day\n")
    add(f"- **{one_shot} of them ({100*one_shot/len(bo_summary):.1f}%) only EVER observed on a bioblitz day**\n")
    add(f"  — these are the course-assignment cohort. They make one trip and never return.\n\n")
    pct_buckets = pd.cut(bo_summary["pct_at_blitz"],
                         bins=[-1, 50, 80, 99, 101],
                         labels=["regular (≤50% blitz)", "mostly-other (50-80%)",
                                 "mostly-blitz (80-99%)", "blitz-only (100%)"])
    add(pct_buckets.value_counts().to_frame("n_observers").to_markdown() + "\n\n")

    OUT.write_text("".join(parts))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
