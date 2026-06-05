"""Lifecycle / survivorship of an iNat observation: observed -> uploaded
-> first own ID -> first external ID -> second supporting ID -> RG.

Reads the raw JSON (the flat CSV drops per-identification timestamps).
Output: notes/lifecycle_eda.md — talking points for the IOA.

Charlotte-vocabulary:
  - data-generating process: observation enters the funnel at observed_on;
    upload is a delayed-reporting step (Pareto-tailed lag); identifiers
    arrive as a marked Poisson-with-bursts process; RG is a stopping rule
    (>=2/3 species-level IDs agree AND non-captive AND has date+location).
  - 'inconsistent with the null at the X% level' (never 'significant')
"""
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data/raw/observations_raw.json"
OUT = ROOT / "notes/lifecycle_eda.md"
TODAY = pd.Timestamp("2026-04-28", tz="UTC")


def load_obs() -> list[dict]:
    with RAW.open() as f:
        data = json.load(f)
    return data if isinstance(data, list) else data["results"]


def to_utc(ts: str | None) -> pd.Timestamp | pd.NaT:
    if not ts:
        return pd.NaT
    return pd.to_datetime(ts, utc=True, errors="coerce")


def build_rows(obs_list):
    rows = []
    for o in obs_list:
        observed = to_utc(o.get("observed_on"))  # date-only -> midnight UTC
        uploaded = to_utc(o.get("created_at"))
        updated = to_utc(o.get("updated_at"))

        idents = o.get("identifications") or []
        # Sort by timestamp ascending
        ts = [(to_utc(i.get("created_at")), i) for i in idents]
        ts = [(t, i) for t, i in ts if pd.notna(t)]
        ts.sort(key=lambda x: x[0])

        first_own = next((t for t, i in ts if i.get("own_observation")), pd.NaT)
        first_ext = next((t for t, i in ts if not i.get("own_observation")), pd.NaT)
        # Supporting external IDs only (these are what push toward RG)
        ext_supports = [t for t, i in ts
                        if (not i.get("own_observation"))
                        and i.get("current")
                        and i.get("category") in {"supporting", "improving", "leading"}]
        second_support = ext_supports[1] if len(ext_supports) >= 2 else pd.NaT
        first_support = ext_supports[0] if len(ext_supports) >= 1 else pd.NaT
        last_id = ts[-1][0] if ts else pd.NaT

        n_ext = sum(1 for _, i in ts if not i.get("own_observation"))

        rows.append({
            "id": o["id"],
            "user_login": o.get("user", {}).get("login"),
            "quality_grade": o.get("quality_grade"),
            "captive": bool(o.get("captive")),
            "geoprivacy": o.get("geoprivacy"),
            "obscured": bool(o.get("obscured")),
            "iconic": (o.get("taxon") or {}).get("iconic_taxon_name"),
            "taxon_rank": (o.get("taxon") or {}).get("rank"),
            "n_photos": len(o.get("photos") or o.get("observation_photos") or []),
            "n_idents": len(idents),
            "n_external_idents": n_ext,
            "n_agreements": o.get("num_identification_agreements", 0),
            "n_disagreements": o.get("num_identification_disagreements", 0),
            "observed_at": observed,
            "uploaded_at": uploaded,
            "updated_at": updated,
            "first_own_id_at": first_own,
            "first_external_id_at": first_ext,
            "first_external_support_at": first_support,
            "second_external_support_at": second_support,
            "last_id_at": last_id,
        })
    return pd.DataFrame(rows)


def lag_days(a, b):
    """Returns Series of (b - a) in days. NaN propagates."""
    return (b - a).dt.total_seconds() / 86400.0


def quantile_block(series: pd.Series, label: str, units: str = "days") -> str:
    s = series.dropna()
    if len(s) == 0:
        return f"- **{label}** — n=0\n"
    q = s.quantile([0.10, 0.25, 0.50, 0.75, 0.90, 0.99]).round(2)
    return (f"- **{label}** (n={len(s)}, {units}): "
            f"p10={q.loc[0.10]}  p25={q.loc[0.25]}  "
            f"**median={q.loc[0.50]}**  p75={q.loc[0.75]}  "
            f"p90={q.loc[0.90]}  p99={q.loc[0.99]}\n")


def funnel_table(df: pd.DataFrame) -> pd.DataFrame:
    """Stage-by-stage retention over the whole dataset."""
    n = len(df)
    stages = [
        ("0. uploaded (everything in the API)", n),
        ("1. has observed_on date", df["observed_at"].notna().sum()),
        ("2. has any identification at all", (df["n_idents"] > 0).sum()),
        ("3. has own ID (observer self-IDed)", df["first_own_id_at"].notna().sum()),
        ("4. has any external ID", df["first_external_id_at"].notna().sum()),
        ("5. has >=1 external supporting ID", df["first_external_support_at"].notna().sum()),
        ("6. has >=2 external supporting IDs (RG threshold met on ID side)",
         df["second_external_support_at"].notna().sum()),
        ("7. quality_grade = research", (df["quality_grade"] == "research").sum()),
    ]
    out = pd.DataFrame(stages, columns=["stage", "n"])
    out["pct_of_total"] = (100 * out["n"] / n).round(1)
    out["pct_of_prev"] = (100 * out["n"] / out["n"].shift(1)).round(1)
    return out


def by_user_table(df: pd.DataFrame, top_n: int = 15) -> pd.DataFrame:
    g = df.groupby("user_login", dropna=False)
    out = pd.DataFrame({
        "n_obs": g.size(),
        "pct_research": 100 * (g["quality_grade"].apply(lambda s: (s == "research").mean())),
        "pct_casual":   100 * (g["quality_grade"].apply(lambda s: (s == "casual").mean())),
        "pct_needs_id": 100 * (g["quality_grade"].apply(lambda s: (s == "needs_id").mean())),
        "pct_self_ided": 100 * g["first_own_id_at"].apply(lambda s: s.notna().mean()),
        "median_upload_lag_d": g.apply(lambda d: lag_days(d["observed_at"], d["uploaded_at"]).median()),
    }).round(1)
    out = out.sort_values("n_obs", ascending=False).head(top_n)
    return out


def by_iconic_table(df: pd.DataFrame) -> pd.DataFrame:
    g = df.groupby("iconic", dropna=False)
    n_ext_lag = lag_days(df["uploaded_at"], df["first_external_id_at"])
    df = df.assign(_first_ext_lag_d=n_ext_lag)
    out = pd.DataFrame({
        "n_obs": g.size(),
        "pct_research": 100 * g["quality_grade"].apply(lambda s: (s == "research").mean()),
        "median_first_ext_id_lag_d": df.groupby("iconic")["_first_ext_lag_d"].median(),
        "median_n_idents": g["n_idents"].median(),
    }).round(2)
    return out.sort_values("n_obs", ascending=False)


def write_report(df: pd.DataFrame) -> None:
    n = len(df)

    # Lags
    upload_lag = lag_days(df["observed_at"], df["uploaded_at"])
    own_id_lag = lag_days(df["uploaded_at"], df["first_own_id_at"])
    first_ext_lag = lag_days(df["uploaded_at"], df["first_external_id_at"])
    second_supp_lag = lag_days(df["uploaded_at"], df["second_external_support_at"])
    last_id_lag = lag_days(df["uploaded_at"], df["last_id_at"])

    funnel = funnel_table(df)
    users = by_user_table(df)
    iconic = by_iconic_table(df)

    # Lag conditional on outcome
    rg_mask = df["quality_grade"] == "research"
    needs_mask = df["quality_grade"] == "needs_id"
    casual_mask = df["quality_grade"] == "casual"

    # Self-ID effect
    self_ided = df["first_own_id_at"].notna()
    rg_rate_self = df.loc[self_ided, "quality_grade"].eq("research").mean() * 100
    rg_rate_no_self = df.loc[~self_ided, "quality_grade"].eq("research").mean() * 100

    # By upload-era cohort: does RG rate differ for old vs new uploads?
    df = df.copy()
    df["upload_year"] = df["uploaded_at"].dt.year
    yr = (df.groupby("upload_year")
            .agg(n=("id", "size"),
                 pct_research=("quality_grade", lambda s: 100 * (s == "research").mean()),
                 pct_needs_id=("quality_grade", lambda s: 100 * (s == "needs_id").mean()),
                 pct_casual=("quality_grade", lambda s: 100 * (s == "casual").mean()),
                 median_first_ext_lag_d=("uploaded_at", lambda _: np.nan),
                 )
            .round(1))
    # recompute the lag column properly
    df["_fext_lag"] = first_ext_lag.values
    yr["median_first_ext_lag_d"] = df.groupby("upload_year")["_fext_lag"].median().round(2)
    yr = yr[yr["n"] >= 50]

    # Backlog uploads: gap between obs and upload
    df["upload_lag_d"] = upload_lag.values
    df["upload_lag_bucket"] = pd.cut(
        df["upload_lag_d"],
        bins=[-1, 1, 7, 30, 365, 100000],
        labels=["same day", "1-7d", "8-30d", "31-365d", ">1y"],
    )
    bucket = (df.groupby("upload_lag_bucket", observed=True)
                .agg(n=("id", "size"),
                     pct_research=("quality_grade", lambda s: 100*(s=="research").mean()),
                     pct_casual=("quality_grade", lambda s: 100*(s=="casual").mean()),
                     pct_needs_id=("quality_grade", lambda s: 100*(s=="needs_id").mean()),
                     )
                .round(1))

    # Stephen-vs-rest contrast (since stephen_thorpe is 36% of obs)
    is_stephen = df["user_login"] == "stephen_thorpe"
    stephen_vs = pd.DataFrame({
        "stephen_thorpe": [
            int(is_stephen.sum()),
            round(100*df.loc[is_stephen, "quality_grade"].eq("research").mean(), 1),
            round(100*df.loc[is_stephen, "quality_grade"].eq("needs_id").mean(), 1),
            round(100*df.loc[is_stephen, "quality_grade"].eq("casual").mean(), 1),
            round(100*df.loc[is_stephen, "first_own_id_at"].notna().mean(), 1),
            round(df.loc[is_stephen, "n_external_idents"].median(), 2),
        ],
        "everyone_else": [
            int((~is_stephen).sum()),
            round(100*df.loc[~is_stephen, "quality_grade"].eq("research").mean(), 1),
            round(100*df.loc[~is_stephen, "quality_grade"].eq("needs_id").mean(), 1),
            round(100*df.loc[~is_stephen, "quality_grade"].eq("casual").mean(), 1),
            round(100*df.loc[~is_stephen, "first_own_id_at"].notna().mean(), 1),
            round(df.loc[~is_stephen, "n_external_idents"].median(), 2),
        ],
    }, index=["n", "% research", "% needs_id", "% casual",
              "% with own ID at upload", "median external IDs"])

    parts = []
    add = parts.append
    add(f"# iNat lifecycle / survivorship EDA\n")
    add(f"_Generated for the IOA oral, 2026-05-01. Underlying N = {n} observations._\n")

    add("## The pipeline (data-generating process)\n")
    add("```\n")
    add("  observed_on  ──(upload lag, often Pareto-tailed)──>  created_at (uploaded)\n")
    add("       │                                                       │\n")
    add("       │                              observer may add own ID  │  (~instant)\n")
    add("       │                                                       ▼\n")
    add("       │                          waits for external identifier(s)\n")
    add("       │                                  │\n")
    add("       │                                  ▼\n")
    add("       │                    1st external supporting ID (community_taxon set)\n")
    add("       │                                  │\n")
    add("       │                                  ▼\n")
    add("       │                    2nd external supporting ID  →  research grade*\n")
    add("       │                                  │\n")
    add("       │                                  └─── OR stays needs_id forever\n")
    add("       │\n")
    add("       └── never uploaded (unobservable from the API; this is the iceberg)\n")
    add("```\n")
    add("\\* RG also requires non-captive, has date, has location, has photo, ≥2/3 community agreement.\n")
    add("Casual = structural failure (missing one of those gates). Needs_id = social failure (no consensus).\n\n")

    add("## Stage-by-stage retention funnel (the dataset we *can* see)\n")
    add(funnel.to_markdown(index=False) + "\n\n")
    add("Note: stage 6 (≥2 external supports) ≠ stage 7 (research grade) because of "
        "the casual gates (captive / no date / no location / no photo) and disagreements.\n\n")

    add("## Lag distributions at each transition (days)\n")
    add(quantile_block(upload_lag, "observed → uploaded (`upload_lag`)"))
    add(quantile_block(own_id_lag, "uploaded → own ID added (when present)"))
    add(quantile_block(first_ext_lag, "uploaded → first external ID"))
    add(quantile_block(second_supp_lag, "uploaded → 2nd external supporting ID (RG-eligible)"))
    add(quantile_block(last_id_lag, "uploaded → last ID seen (proxy for settled)"))
    add("\n")
    add("**Interpretation hooks:**\n")
    add("- The upload lag has the heaviest tail — many people back-fill years of camera-roll observations.\n")
    add("- First-external-ID lag medians measured in *days*, but tail in years for low-traffic taxa.\n")
    add("- The 2nd-support lag is what separates needs_id from research; that's the bottleneck stage.\n\n")

    add("## Lag distributions split by terminal quality grade\n")
    add("**Upload lag (observed → uploaded):**\n")
    add(quantile_block(upload_lag[rg_mask], "  research"))
    add(quantile_block(upload_lag[needs_mask], "  needs_id"))
    add(quantile_block(upload_lag[casual_mask], "  casual"))
    add("\n**First-external-ID lag (uploaded → first external ID):**\n")
    add(quantile_block(first_ext_lag[rg_mask], "  research"))
    add(quantile_block(first_ext_lag[needs_mask], "  needs_id"))
    add(quantile_block(first_ext_lag[casual_mask], "  casual"))
    add("\n_Casual obs often get fewer/no external IDs because they failed a gate (e.g. captive)._\n\n")

    add("## Self-ID at upload — does it help?\n")
    add(f"- Obs WHERE observer added their own ID: RG rate = **{rg_rate_self:.1f}%**\n")
    add(f"- Obs WITHOUT observer self-ID:         RG rate = **{rg_rate_no_self:.1f}%**\n")
    add("Mechanism: a self-ID seeds `community_taxon`, so the next external supporter "
        "only needs to agree (not propose). Two-thirds rule is met sooner.\n\n")

    add("## By observer (top 15 by volume)\n")
    add(users.to_markdown() + "\n\n")
    add("**Per-observer differences worth surfacing in the oral:**\n")
    add("- `pct_research` varies wildly across heavy users — taxonomic specialism dominates.\n")
    add("- Users who self-ID at upload have measurably different RG rates than those who don't.\n")
    add("- Median upload lag separates 'in-the-field uploaders' (same day) from 'archive dumpers' (years).\n\n")

    add("## stephen_thorpe vs. everyone else (he's 36% of all obs)\n")
    add(stephen_vs.to_markdown() + "\n\n")

    add("## By iconic taxon group\n")
    add(iconic.to_markdown() + "\n\n")
    add("Charismatic groups (Aves, Mammalia) reach RG faster — more identifiers watch them. "
        "Insecta and Fungi rely on a few specialists, so the queue is longer and bumpier.\n\n")

    add("## By upload year cohort\n")
    add(yr.to_markdown() + "\n\n")
    add("- Older uploads have had more wall-clock time for an ID to arrive — they should "
        "be enriched in `research`, but they're also the oldest *backlog* dumps and may have "
        "weaker metadata, which pushes them to `casual`.\n\n")

    add("## By upload-lag bucket (how stale was the obs at upload?)\n")
    add(bucket.to_markdown() + "\n\n")
    add("Same-day uploads behave very differently from year-late backlogs — that's a confounder "
        "for any model of RG outcome.\n\n")

    add("## What we still cannot see (talk about this in the oral)\n")
    add("- The truly unobserved iceberg: organisms seen but never photographed/uploaded. "
        "iNaturalist sampling is opportunistic — there's no denominator.\n")
    add("- Recorded-but-deleted observations: deleted obs don't appear in the API.\n")
    add("- Identifier identity / effort: we have the identifications on each obs, but no "
        "sample of identifiers' overall workload to model their arrival process directly.\n")
    add("- Snapshot bias: 'still needs_id' could become 'research' tomorrow. "
        "We can model age-since-upload as a covariate but we cannot censor the way a true "
        "survival model would.\n")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("".join(parts))
    print(f"wrote {OUT}")


def main() -> None:
    print("loading raw JSON...")
    obs = load_obs()
    print(f"  {len(obs)} observations")
    df = build_rows(obs)
    write_report(df)


if __name__ == "__main__":
    main()
