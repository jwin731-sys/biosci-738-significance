"""Contagion probes — answers two questions:

  Q1: If stephen_thorpe is removed, does another observer take his place
      as 'trigger-prone'? (Pareto-redistribution test.)
  Q2: Are there *bursts* — same species, multiple different observers within
      hours, suggesting physical co-presence (class trip, bioblitz)?
      That's the user's 'people going out together' hypothesis.

Q2 is the cleaner test for *virality*: if obs cluster within Δt ≤ 6 hours by
multiple different observers more than seasonal+observer null predicts, that's
evidence of either (a) shared physical attention or (b) social-media propagation
within a few hours. Either way, more than mere correlation.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data/clean/observations_flat.csv"
RNG = np.random.default_rng(42)
N_PERMS = 300


def load() -> pd.DataFrame:
    df = pd.read_csv(DATA, usecols=[
        "observed_on", "time_observed_at", "user_login", "taxon_name",
        "taxon_iconic_taxon_name", "longitude", "latitude",
    ])
    df["t"] = pd.to_datetime(df["time_observed_at"], errors="coerce", utc=True)
    df.loc[df["t"].isna(), "t"] = pd.to_datetime(
        df.loc[df["t"].isna(), "observed_on"], errors="coerce", utc=True)
    df = df.dropna(subset=["t", "user_login", "taxon_name"]).copy()
    df["sec"] = df["t"].astype("int64") // 10**9
    df["day"] = df["sec"] // 86400
    df["month"] = df["t"].dt.month.astype(np.int8)
    df["date"] = df["t"].dt.date
    return df.reset_index(drop=True)


def cross_observer_pairs(seconds: np.ndarray, obs: np.ndarray, window_sec: int) -> int:
    n = len(seconds)
    if n < 2: return 0
    order = np.argsort(seconds, kind="stable")
    s = seconds[order]; o = obs[order]
    cnt = 0; j = 0
    for i in range(n):
        if j <= i: j = i + 1
        while j < n and s[j] - s[i] <= window_sec: j += 1
        for k in range(i + 1, j):
            if o[i] != o[k]: cnt += 1
    return cnt


def first_in_cluster_mask(seconds: np.ndarray, obs: np.ndarray, window_sec: int) -> np.ndarray:
    n = len(seconds)
    out = np.zeros(n, dtype=bool)
    if n < 2: return out
    order = np.argsort(seconds, kind="stable")
    s = seconds[order]; o = obs[order]
    flag = np.zeros(n, dtype=bool)
    j = 0
    for i in range(n):
        if j <= i: j = i + 1
        while j < n and s[j] - s[i] <= window_sec: j += 1
        for k in range(i + 1, j):
            if o[k] != o[i]:
                flag[i] = True
                break
    out[order] = flag
    return out


def trigger_prone_scan(df: pd.DataFrame, window_days: int = 7,
                       min_obs_species: int = 8, min_obs_observer: int = 10,
                       label: str = "") -> pd.DataFrame:
    window_sec = window_days * 86400
    species_groups = df.groupby("taxon_name")
    species_to_test = [s for s, g in species_groups if len(g) >= min_obs_species]

    df = df.copy()
    df["fic"] = False
    for s in species_to_test:
        g = species_groups.get_group(s)
        idx = g.index.to_numpy()
        secs = g["sec"].to_numpy(np.int64)
        observers = g["user_login"].to_numpy()
        df.loc[idx, "fic"] = first_in_cluster_mask(secs, observers, window_sec)

    summ = (df.groupby("user_login")
            .agg(n_obs=("sec", "size"),
                 fic_obs=("fic", "sum"),
                 n_species=("taxon_name", "nunique"))
            .reset_index())
    summ["fic_rate"] = summ["fic_obs"] / summ["n_obs"]
    summ = summ[summ["n_obs"] >= min_obs_observer].copy()
    obs_idx = {u: i for i, u in enumerate(summ["user_login"])}

    null_fic = np.zeros((N_PERMS, len(summ)), dtype=np.int64)
    for k in range(N_PERMS):
        for s in species_to_test:
            g = species_groups.get_group(s)
            secs = g["sec"].to_numpy(np.int64)
            obs_perm = RNG.permutation(g["user_login"].to_numpy())
            mask = first_in_cluster_mask(secs, obs_perm, window_sec)
            for off, fic in enumerate(mask):
                if fic and obs_perm[off] in obs_idx:
                    null_fic[k, obs_idx[obs_perm[off]]] += 1

    summ["null_mean"] = null_fic.mean(axis=0)
    summ["null_std"] = null_fic.std(axis=0)
    summ["z"] = np.where(summ["null_std"] > 0,
                          (summ["fic_obs"] - summ["null_mean"]) / summ["null_std"], 0.0)
    summ = summ.sort_values("z", ascending=False)
    print(f"\n=== Top trigger-prone {label} (n_observers={len(summ)}) ===")
    print(summ.head(15)[["user_login","n_obs","n_species","fic_obs","null_mean","z"]]
          .to_string(index=False, float_format=lambda x: f"{x:.2f}"))
    return summ


def burst_test(df: pd.DataFrame, window_hours: float = 6.0,
               min_obs_species: int = 5) -> None:
    """For each species, count *tight bursts*: pairs within window_hours by
    different observers. Compare to null where dates of obs are drawn from
    same month-of-year pool (preserves bioblitz/seasonality at the daily
    level but destroys hour-of-day clustering)."""
    window_sec = int(window_hours * 3600)
    species_groups = df.groupby("taxon_name")
    species_to_test = [s for s, g in species_groups if len(g) >= min_obs_species]

    by_month_secs: dict[int, np.ndarray] = {
        m: df.loc[df["month"] == m, "sec"].to_numpy(np.int64) for m in range(1, 13)}
    for m in range(1, 13):
        by_month_secs.setdefault(m, np.array([], dtype=np.int64))

    obs_total, null_means, null_stds, n_obs_list = [], [], [], []
    species_kept = []
    for s in species_to_test:
        g = species_groups.get_group(s)
        secs = g["sec"].to_numpy(np.int64)
        observers = g["user_login"].to_numpy()
        months = g["month"].to_numpy(np.int8)
        observed = cross_observer_pairs(secs, observers, window_sec)
        if observed == 0 and len(g) < 30:
            continue
        nulls = np.empty(N_PERMS, dtype=np.int64)
        for k in range(N_PERMS):
            new_secs = np.empty_like(secs)
            for i, m in enumerate(months):
                pool = by_month_secs[int(m)]
                new_secs[i] = RNG.choice(pool)
            nulls[k] = cross_observer_pairs(new_secs, observers, window_sec)
        obs_total.append(observed)
        null_means.append(nulls.mean())
        null_stds.append(nulls.std())
        n_obs_list.append(len(g))
        species_kept.append(s)

    res = pd.DataFrame({
        "taxon": species_kept, "n_obs": n_obs_list,
        "burst_pairs": obs_total, "null_mean": null_means, "null_std": null_stds,
    })
    res["z"] = np.where(res["null_std"] > 0,
                         (res["burst_pairs"] - res["null_mean"]) / res["null_std"], 0.0)
    res = res.sort_values("z", ascending=False)
    total_obs = res["burst_pairs"].sum()
    total_null = res["null_mean"].sum()
    print(f"\n=== Burst test: same species, Δt ≤ {window_hours}h, different observers ===")
    print(f"Total bursts observed: {total_obs:,.0f}")
    print(f"Total bursts under seasonal null (mean): {total_null:,.0f}")
    print(f"Ratio observed/null: {total_obs/total_null:.2f}x" if total_null else "null=0")
    print("\nTop 20 species by excess-burst z:")
    print(res.head(20).to_string(index=False, float_format=lambda x: f"{x:.2f}"))
    res.to_csv(ROOT / "notes/contagion_bursts.csv", index=False)


def bioblitz_pattern(df: pd.DataFrame) -> None:
    """How concentrated are bursts on specific dates?"""
    daily = df.groupby("date").size().sort_values(ascending=False)
    print(f"\n=== Top 15 most-observed days ===")
    print(daily.head(15).to_string())
    print(f"\nTotal obs: {len(df):,}")
    print(f"Top day: {daily.iloc[0]} obs ({daily.iloc[0]/len(df)*100:.1f}% of dataset)")
    print(f"Top 5 days combined: {daily.head(5).sum()} obs ({daily.head(5).sum()/len(df)*100:.1f}%)")


def main() -> None:
    df = load()
    print(f"loaded {len(df):,} obs", file=sys.stderr)

    bioblitz_pattern(df)

    print("\n" + "="*70)
    print("Q1: Pareto-redistribution test — strip stephen_thorpe")
    print("="*70)
    df_no_stephen = df[df["user_login"] != "stephen_thorpe"].copy()
    print(f"After removing stephen_thorpe: {len(df_no_stephen):,} obs "
          f"({len(df_no_stephen)/len(df)*100:.1f}% of original)")
    summ_orig = trigger_prone_scan(df, label="(full dataset)")
    summ_no = trigger_prone_scan(df_no_stephen, label="(stephen excluded)")
    summ_orig.to_csv(ROOT / "notes/contagion_observers_full.csv", index=False)
    summ_no.to_csv(ROOT / "notes/contagion_observers_no_stephen.csv", index=False)

    print("\n" + "="*70)
    print("Q2: Burst test — physical co-presence proxy (Δt ≤ 6h)")
    print("="*70)
    burst_test(df, window_hours=6.0, min_obs_species=5)

    print("\n" + "="*70)
    print("Q2b: Tight burst — Δt ≤ 1h (very strong co-presence signal)")
    print("="*70)
    burst_test(df, window_hours=1.0, min_obs_species=5)


if __name__ == "__main__":
    main()
