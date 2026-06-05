"""Quick excitement-cluster scan for the iNat City Campus dataset.

Tests two questions:
1. Which species cluster *more tightly* than expected given their own seasonal
   envelope? ("Excitable species")
2. Which observers are disproportionately *first-in-cluster*? ("Trigger-prone
   observers")

Detrending: the per-species null preserves each species' own month-of-year
distribution and randomizes day-within-month + year. This preserves seasonality
(March bioblitz, winter species, etc.) but destroys within-species fine-grained
temporal clustering — which is exactly the signal we want.

The cross-observer pair count uses Δt ≤ 7 days AND different observers, which
rules out one-person batch dumps (e.g., stephen_thorpe submitting 50 obs in a
single afternoon).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data/clean/observations_flat.csv"

WINDOW_DAYS = 7
N_PERMS = 500
MIN_OBS_SPECIES = 8     # minimum obs/species to test
MIN_OBS_OBSERVER = 10   # minimum obs/observer to test
RNG = np.random.default_rng(42)


def load() -> pd.DataFrame:
    df = pd.read_csv(
        DATA,
        usecols=["observed_on", "user_login", "taxon_name",
                 "taxon_iconic_taxon_name", "taxon_rank", "quality_grade"],
    )
    df["date"] = pd.to_datetime(df["observed_on"], errors="coerce")
    df = df.dropna(subset=["date", "user_login", "taxon_name"]).copy()
    df["day"] = (df["date"].astype("int64") // 10**9 // 86400).astype(np.int64)
    df["month"] = df["date"].dt.month.astype(np.int8)
    return df.reset_index(drop=True)


def cross_observer_pairs(days: np.ndarray, obs: np.ndarray, window: int) -> int:
    """Count ordered pairs (i, j) with i<j, days[j]-days[i]<=window, obs[i]!=obs[j].
    Inputs need not be sorted; sorts internally."""
    n = len(days)
    if n < 2:
        return 0
    order = np.argsort(days, kind="stable")
    d = days[order]
    o = obs[order]
    cnt = 0
    j = 0
    for i in range(n):
        if j <= i:
            j = i + 1
        while j < n and d[j] - d[i] <= window:
            j += 1
        for k in range(i + 1, j):
            if o[i] != o[k]:
                cnt += 1
    return cnt


def first_in_cluster_mask(days: np.ndarray, obs: np.ndarray, window: int) -> np.ndarray:
    """For each obs, True if at least one DIFFERENT observer observed the same
    species within (t, t+window] days. (Forward-looking only — i.e. this obs is
    'before' the cluster, candidate trigger.)"""
    n = len(days)
    mask = np.zeros(n, dtype=bool)
    if n < 2:
        return mask
    order = np.argsort(days, kind="stable")
    d = days[order]
    o = obs[order]
    out = np.zeros(n, dtype=bool)
    j = 0
    for i in range(n):
        if j <= i:
            j = i + 1
        while j < n and d[j] - d[i] <= window:
            j += 1
        for k in range(i + 1, j):
            if o[k] != o[i]:
                out[i] = True
                break
    mask[order] = out
    return mask


def shuffle_dates_preserving_phenology(orig_days: np.ndarray, orig_months: np.ndarray,
                                        global_days_by_month: dict[int, np.ndarray],
                                        rng: np.random.Generator) -> np.ndarray:
    """For each obs, draw a new day-of-dataset from the global pool of days
    that fall in the same month-of-year as the original. Preserves per-obs
    month, randomizes year + day-within-month."""
    new = np.empty_like(orig_days)
    for i, m in enumerate(orig_months):
        pool = global_days_by_month[int(m)]
        new[i] = rng.choice(pool)
    return new


def main() -> None:
    df = load()
    print(f"loaded {len(df):,} obs with valid dates", file=sys.stderr)
    print(f"date range: {df['date'].min().date()} → {df['date'].max().date()}", file=sys.stderr)
    print(f"unique species: {df['taxon_name'].nunique():,}", file=sys.stderr)
    print(f"unique observers: {df['user_login'].nunique():,}", file=sys.stderr)

    # global day-of-dataset pool, indexed by month-of-year (preserves seasonality)
    global_days_by_month: dict[int, np.ndarray] = {
        m: df.loc[df["month"] == m, "day"].to_numpy(np.int64)
        for m in range(1, 13)
    }
    for m in range(1, 13):
        global_days_by_month.setdefault(m, np.array([], dtype=np.int64))

    # ---- per-species excitement test ----
    species_results = []
    species_groups = df.groupby("taxon_name")
    species_to_test = [s for s, g in species_groups if len(g) >= MIN_OBS_SPECIES]
    print(f"testing {len(species_to_test)} species (≥{MIN_OBS_SPECIES} obs)", file=sys.stderr)

    for s in species_to_test:
        g = species_groups.get_group(s)
        days = g["day"].to_numpy(np.int64)
        observers = g["user_login"].to_numpy()
        months = g["month"].to_numpy(np.int8)
        n = len(g)

        observed = cross_observer_pairs(days, observers, WINDOW_DAYS)

        # null: keep per-obs month, redraw day from global month pool
        null_counts = np.empty(N_PERMS, dtype=np.int64)
        for k in range(N_PERMS):
            null_days = shuffle_dates_preserving_phenology(
                days, months, global_days_by_month, RNG
            )
            null_counts[k] = cross_observer_pairs(null_days, observers, WINDOW_DAYS)

        mu, sd = null_counts.mean(), null_counts.std()
        z = (observed - mu) / sd if sd > 0 else 0.0
        p_one_sided = float((null_counts >= observed).mean())

        iconic = g["taxon_iconic_taxon_name"].mode()
        iconic = iconic.iloc[0] if len(iconic) else "?"

        species_results.append({
            "taxon": s,
            "iconic": iconic,
            "n_obs": n,
            "n_observers": g["user_login"].nunique(),
            "obs_pairs": observed,
            "null_mean": mu,
            "null_std": sd,
            "z": z,
            "p_perm": p_one_sided,
        })

    sp = pd.DataFrame(species_results).sort_values("z", ascending=False)
    print("\n=== Top 25 'excitable' species (high z = clusters more than seasonal null) ===")
    print(sp.head(25).to_string(index=False, float_format=lambda x: f"{x:.2f}"))

    print("\n=== Bottom 10 'avoidant' species (negative z = clusters less than null) ===")
    print(sp.tail(10).to_string(index=False, float_format=lambda x: f"{x:.2f}"))

    sp.to_csv(ROOT / "notes/excitement_species.csv", index=False)
    print(f"\nwrote {ROOT / 'notes/excitement_species.csv'}", file=sys.stderr)

    # ---- per-observer first-in-cluster test ----
    # For the OBSERVED data, mark each obs as first-in-cluster within its species
    df["fic"] = False
    for s in species_to_test:
        g = species_groups.get_group(s)
        idx = g.index.to_numpy()
        days = g["day"].to_numpy(np.int64)
        observers = g["user_login"].to_numpy()
        mask = first_in_cluster_mask(days, observers, WINDOW_DAYS)
        df.loc[idx, "fic"] = mask

    obs_summary = (df.groupby("user_login")
                   .agg(n_obs=("day", "size"),
                        fic_obs=("fic", "sum"),
                        n_species=("taxon_name", "nunique"))
                   .reset_index())
    obs_summary["fic_rate"] = obs_summary["fic_obs"] / obs_summary["n_obs"]
    obs_summary = obs_summary[obs_summary["n_obs"] >= MIN_OBS_OBSERVER].copy()
    print(f"\ntesting {len(obs_summary)} observers (≥{MIN_OBS_OBSERVER} obs)", file=sys.stderr)

    # observer null: permute the *observer* labels across all observations within
    # each species block, recompute per-observer fic count. Tests whether each
    # observer is disproportionately first-in-cluster vs. random observer assignment.
    null_fic_per_observer = np.zeros((N_PERMS, len(obs_summary)), dtype=np.int64)
    obs_idx = {u: i for i, u in enumerate(obs_summary["user_login"])}

    for k in range(N_PERMS):
        df_perm_fic = np.zeros(len(df), dtype=bool)
        # for each species, shuffle observers among that species' obs
        for s in species_to_test:
            g = species_groups.get_group(s)
            idx = g.index.to_numpy()
            days = g["day"].to_numpy(np.int64)
            observers_perm = RNG.permutation(g["user_login"].to_numpy())
            mask = first_in_cluster_mask(days, observers_perm, WINDOW_DAYS)
            # assign mask back: which OBSERVER (in permuted assignment) gets the fic credit
            for off, fic in enumerate(mask):
                if fic:
                    o = observers_perm[off]
                    if o in obs_idx:
                        null_fic_per_observer[k, obs_idx[o]] += 1

    obs_summary["null_mean"] = null_fic_per_observer.mean(axis=0)
    obs_summary["null_std"] = null_fic_per_observer.std(axis=0)
    obs_summary["z"] = np.where(
        obs_summary["null_std"] > 0,
        (obs_summary["fic_obs"] - obs_summary["null_mean"]) / obs_summary["null_std"],
        0.0,
    )

    obs_summary = obs_summary.sort_values("z", ascending=False)
    print("\n=== Top 20 'trigger-prone' observers (high z = first-in-cluster more than chance) ===")
    cols = ["user_login", "n_obs", "n_species", "fic_obs", "null_mean", "z"]
    print(obs_summary.head(20)[cols].to_string(index=False, float_format=lambda x: f"{x:.2f}"))

    obs_summary.to_csv(ROOT / "notes/excitement_observers.csv", index=False)
    print(f"\nwrote {ROOT / 'notes/excitement_observers.csv'}", file=sys.stderr)


if __name__ == "__main__":
    main()
