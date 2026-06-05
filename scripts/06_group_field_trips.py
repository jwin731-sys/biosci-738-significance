"""Test the hypothesis: stephen_thorpe (or another power-IDer) is leading
field trips for students, and those co-observation events show unusually
high RG rates because the leader IDs everyone's uploads on the spot.

Three angles:
  (A) "Group day" detection: observed_on dates where many distinct users
      were active in the same small window. Compare RG rates.
  (B) Power-IDer reach: how many OTHER users' observations does each
      heavy IDer appear on? Does being IDed by stephen_thorpe predict RG?
  (C) Co-observation clusters: same observed_on date + spatial proximity
      + multiple users → likely guided trip.
"""
from __future__ import annotations

import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data/raw/observations_raw.json"
OUT = ROOT / "notes/group_trips_hypothesis.md"


def load():
    with RAW.open() as f:
        data = json.load(f)
    return data if isinstance(data, list) else data["results"]


def main():
    obs_list = load()
    print(f"loaded {len(obs_list)} obs")

    # Per-obs row + per-identification long table
    rows = []
    id_rows = []
    for o in obs_list:
        observer = (o.get("user") or {}).get("login")
        observed = pd.to_datetime(o.get("observed_on"), errors="coerce", utc=True)
        rows.append({
            "id": o["id"],
            "observer": observer,
            "observed_date": observed.date() if pd.notna(observed) else None,
            "lat": o.get("geojson", {}).get("coordinates", [None, None])[1] if o.get("geojson") else None,
            "lon": o.get("geojson", {}).get("coordinates", [None, None])[0] if o.get("geojson") else None,
            "quality_grade": o.get("quality_grade"),
            "captive": bool(o.get("captive")),
            "iconic": (o.get("taxon") or {}).get("iconic_taxon_name"),
        })
        for ident in (o.get("identifications") or []):
            id_rows.append({
                "obs_id": o["id"],
                "observer": observer,
                "identifier": (ident.get("user") or {}).get("login"),
                "ident_at": pd.to_datetime(ident.get("created_at"), errors="coerce", utc=True),
                "own_observation": bool(ident.get("own_observation")),
                "category": ident.get("category"),
                "current": bool(ident.get("current")),
                "obs_quality_grade": o.get("quality_grade"),
            })

    df = pd.DataFrame(rows)
    ids = pd.DataFrame(id_rows)

    # -------------------------------------------------------------------
    # (A) Group-day detection
    # A "group day" = an observed_date with many distinct observers active
    # -------------------------------------------------------------------
    by_day = (df.dropna(subset=["observed_date"])
                .groupby("observed_date")
                .agg(n_obs=("id", "size"),
                     n_observers=("observer", "nunique"),
                     pct_research=("quality_grade",
                                    lambda s: 100*(s=="research").mean()),
                     pct_casual=("quality_grade",
                                  lambda s: 100*(s=="casual").mean()),
                     pct_needs_id=("quality_grade",
                                    lambda s: 100*(s=="needs_id").mean()),
                     )
                .reset_index())

    # Bucket days by how many observers were active
    by_day["observer_bucket"] = pd.cut(
        by_day["n_observers"],
        bins=[0, 1, 2, 4, 9, 1000],
        labels=["1 (solo)", "2", "3-4", "5-9", "10+"],
    )

    bucket_summary = (by_day.groupby("observer_bucket", observed=True)
                            .agg(n_days=("observed_date", "size"),
                                 total_obs=("n_obs", "sum"),
                                 mean_pct_research=("pct_research", "mean"),
                                 mean_pct_casual=("pct_casual", "mean"),
                                 mean_pct_needs_id=("pct_needs_id", "mean"))
                            .round(1))

    # Now do it OBS-LEVEL (weighted by obs, not by day):
    # tag each obs with the n_observers on its observed_date
    day_observers = by_day.set_index("observed_date")["n_observers"]
    df["n_observers_that_day"] = df["observed_date"].map(day_observers)
    df["day_bucket"] = pd.cut(
        df["n_observers_that_day"],
        bins=[0, 1, 2, 4, 9, 1000],
        labels=["1 (solo)", "2", "3-4", "5-9", "10+"],
    )
    obs_bucket = (df.groupby("day_bucket", observed=True)
                    .agg(n_obs=("id", "size"),
                         pct_research=("quality_grade",
                                        lambda s: 100*(s=="research").mean()),
                         pct_casual=("quality_grade",
                                      lambda s: 100*(s=="casual").mean()),
                         pct_needs_id=("quality_grade",
                                        lambda s: 100*(s=="needs_id").mean()))
                    .round(1))

    # Top 15 highest-observer-count days (these are the smoking gun for trips)
    top_days = (by_day.sort_values("n_observers", ascending=False).head(15))

    # -------------------------------------------------------------------
    # (B) Power-IDer reach
    # For each identifier, count distinct OTHER observers they've IDed for
    # -------------------------------------------------------------------
    ext_ids = ids[~ids["own_observation"]].copy()  # external IDs only
    reach = (ext_ids.groupby("identifier")
                    .agg(n_idents=("obs_id", "size"),
                         n_obs=("obs_id", "nunique"),
                         n_observers_helped=("observer", "nunique"))
                    .sort_values("n_idents", ascending=False)
                    .head(15))

    # Does being externally-IDed by stephen_thorpe predict RG?
    # Restrict to obs where Stephen is NOT the observer
    not_stephen = df[df["observer"] != "stephen_thorpe"].copy()
    stephen_ided_obs = set(
        ext_ids[ext_ids["identifier"] == "stephen_thorpe"]["obs_id"]
    )
    not_stephen["stephen_ided"] = not_stephen["id"].isin(stephen_ided_obs)
    stephen_id_effect = (not_stephen.groupby("stephen_ided")
                                    .agg(n=("id", "size"),
                                         pct_research=("quality_grade",
                                                        lambda s: 100*(s=="research").mean()),
                                         pct_casual=("quality_grade",
                                                      lambda s: 100*(s=="casual").mean()),
                                         pct_needs_id=("quality_grade",
                                                        lambda s: 100*(s=="needs_id").mean()))
                                    .round(1))

    # Same test for the next-tier power IDer (whoever it is)
    runner_up = reach.index[1] if "stephen_thorpe" in reach.index[:3] else reach.index[0]
    if reach.index[0] != "stephen_thorpe":
        runner_up = reach.index[0]
    runner_obs = set(ext_ids[ext_ids["identifier"] == runner_up]["obs_id"])
    not_runner = df[df["observer"] != runner_up].copy()
    not_runner[f"{runner_up}_ided"] = not_runner["id"].isin(runner_obs)
    runner_effect = (not_runner.groupby(f"{runner_up}_ided")
                               .agg(n=("id", "size"),
                                    pct_research=("quality_grade",
                                                   lambda s: 100*(s=="research").mean()))
                               .round(1))

    # -------------------------------------------------------------------
    # (C) Co-observation clusters: same date + Stephen present + others
    # -------------------------------------------------------------------
    stephen_dates = set(df[df["observer"] == "stephen_thorpe"]["observed_date"].dropna())
    df["stephen_was_out"] = df["observed_date"].isin(stephen_dates)
    stephen_day_effect = (df[df["observer"] != "stephen_thorpe"]
                          .groupby("stephen_was_out")
                          .agg(n=("id", "size"),
                               pct_research=("quality_grade",
                                              lambda s: 100*(s=="research").mean()),
                               pct_casual=("quality_grade",
                                            lambda s: 100*(s=="casual").mean()),
                               pct_needs_id=("quality_grade",
                                              lambda s: 100*(s=="needs_id").mean()))
                          .round(1))

    # Stronger test: dates where Stephen AND >=3 other distinct observers were active
    # AND Stephen IDed several of the other-observer obs that day
    stephen_obs_by_date = df[df["observer"] == "stephen_thorpe"].groupby("observed_date").size()
    others_by_date = (df[df["observer"] != "stephen_thorpe"]
                      .groupby("observed_date")
                      .agg(n_others=("id", "size"),
                           n_other_observers=("observer", "nunique")))
    trip_candidates = others_by_date.join(stephen_obs_by_date.rename("n_stephen"), how="inner")
    trip_candidates = trip_candidates[
        (trip_candidates["n_stephen"] >= 1) &
        (trip_candidates["n_other_observers"] >= 3)
    ].sort_values("n_other_observers", ascending=False).head(20)

    # For each candidate trip date, what frac of other-observer obs did Stephen ID?
    trip_dates = trip_candidates.index.tolist()
    extra = []
    for d in trip_dates:
        day_obs = df[(df["observed_date"] == d) & (df["observer"] != "stephen_thorpe")]
        day_obs_ids = set(day_obs["id"])
        stephen_ided_on_day = len(day_obs_ids & stephen_ided_obs)
        rg_rate = 100 * day_obs["quality_grade"].eq("research").mean()
        extra.append({
            "date": d,
            "stephen_ided_on_day": stephen_ided_on_day,
            "frac_stephen_ided": round(100*stephen_ided_on_day/max(len(day_obs_ids),1), 1),
            "pct_research_others": round(rg_rate, 1),
        })
    trip_table = (trip_candidates.reset_index()
                                  .merge(pd.DataFrame(extra), left_on="observed_date", right_on="date")
                                  .drop(columns=["date"]))

    # Baseline: RG rate for non-stephen obs on dates where Stephen wasn't out
    baseline_rg = 100*df[(df["observer"] != "stephen_thorpe") &
                         (~df["observed_date"].isin(stephen_dates))]["quality_grade"].eq("research").mean()

    # -------------------------------------------------------------------
    # Write report
    # -------------------------------------------------------------------
    parts = []
    add = parts.append
    add("# Hypothesis: power-IDers leading field trips boost RG rates\n\n")
    add("Specifically: does stephen_thorpe (or another heavy IDer) take students/groups out, "
        "ID their obs on the spot, and produce co-observation events with anomalously high RG rates?\n\n")

    add("## (A) RG rate vs. how many distinct observers were active that day\n\n")
    add("**Day-level (each row = one calendar date, equal weight):**\n\n")
    add(bucket_summary.to_markdown() + "\n\n")
    add("**Obs-level (each obs weighted; reflects what fraction of the dataset comes from group days):**\n\n")
    add(obs_bucket.to_markdown() + "\n\n")

    add("**Top 15 highest-observer-count days (likely guided events):**\n\n")
    add(top_days.to_markdown(index=False) + "\n\n")

    add("## (B) Power-IDer reach\n\n")
    add("Top 15 external identifiers — how many distinct observers they help:\n\n")
    add(reach.to_markdown() + "\n\n")
    add("### Does being externally-IDed by stephen_thorpe predict RG (excluding his own obs)?\n\n")
    add(stephen_id_effect.to_markdown() + "\n\n")
    add(f"### Same test for the runner-up IDer (`{runner_up}`):\n\n")
    add(runner_effect.to_markdown() + "\n\n")

    add("## (C) Stephen-was-out-that-day effect\n\n")
    add("For obs by everyone EXCEPT stephen_thorpe, split by whether Stephen also "
        "made an observation on the same date:\n\n")
    add(stephen_day_effect.to_markdown() + "\n\n")
    add(f"Baseline (obs by non-Stephen, on dates Stephen wasn't out): RG = {baseline_rg:.1f}%\n\n")

    add("**Likely guided-trip dates (Stephen + ≥3 other observers active):**\n\n")
    add(trip_table.to_markdown(index=False) + "\n\n")
    add("Read this table as: on these days, what fraction of NON-Stephen observations "
        "did Stephen end up IDing externally? And what was the RG rate for those non-Stephen obs?\n\n")

    add("## How to read this in the oral\n\n")
    add("- If the 10+ observer day-bucket has dramatically higher RG than solo days → group events drive RG.\n")
    add("- If `stephen_ided=True` obs have much higher RG → he's a one-man identification engine.\n")
    add("- If both are true and overlap → it's the *combination*: he runs trips and also IDs other people's lone uploads from his armchair.\n")
    add("- The trip-candidates table will show specific dates worth eyeballing — match against UoA academic calendar for ENVSCI/BIOSCI fieldwork.\n")

    OUT.write_text("".join(parts))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
