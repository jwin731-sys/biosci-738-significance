"""Quality-grade 'survival' snapshot for the iNat City Campus dataset.

Backup talking-track for the IOA. Frames the question as: given that an obs
exists at scrape time, what is the probability it has reached *research grade*?

This is NOT a true survival analysis — we have only snapshot quality_grade,
not the transition timestamps. The cleanest defensible framing is therefore
discrete-time logistic regression with `obs_age_days = today - created_at` as
one covariate. Older obs that haven't crossed the threshold are 'stuck'; young
obs may still be maturing. We can't disentangle those two without ID-history.

Two views:
  (A) Full: research_grade vs (casual ∪ needs_id) — includes the structural
      casual obs (captive, missing date/location, geoprivacy obscured).
  (B) Conditional: research_grade vs needs_id only — the obs that *can* reach
      research grade. This is the more interpretable 'time-to-consensus' view.

Charlotte-vocabulary notes for the oral:
  - data-generating process: obs enters at created_at; identifiers arrive as
    a marked Poisson process; transitions to research grade once 2/3 species-
    level IDs agree AND obs has date/location/non-captive flags. We model the
    *outcome* of that process at scrape time.
  - 'the data are inconsistent with the null at the X% level' (never 'sig')
  - varying effect would be on observer / iconic taxon — we report cluster-
    robust SEs by user_login since obs are deeply nested in observers.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data/clean/observations_flat.csv"
TODAY = pd.Timestamp("2026-04-28", tz="Pacific/Auckland")

USE_COLS = [
    "observed_on", "created_at", "quality_grade", "captive",
    "n_photos", "identifications_count",
    "num_identification_agreements", "num_identification_disagreements",
    "comments_count", "user_login", "user_observations_count",
    "taxon_iconic_taxon_name", "taxon_rank", "geoprivacy", "obscured",
]


def load() -> pd.DataFrame:
    df = pd.read_csv(DATA, usecols=USE_COLS)
    df["created_at"] = pd.to_datetime(df["created_at"], errors="coerce", utc=True)
    df["observed_on_dt"] = pd.to_datetime(df["observed_on"], errors="coerce")
    df = df.dropna(subset=["created_at", "quality_grade", "user_login"]).copy()
    df["age_days"] = (TODAY - df["created_at"]).dt.total_seconds() / 86400.0
    df["age_years"] = df["age_days"] / 365.25
    df["log_age"] = np.log1p(df["age_days"].clip(lower=0))
    df["is_research"] = (df["quality_grade"] == "research").astype(int)
    df["is_casual"] = (df["quality_grade"] == "casual").astype(int)
    df["is_needs_id"] = (df["quality_grade"] == "needs_id").astype(int)
    df["captive_bin"] = df["captive"].fillna(False).astype(int)
    df["species_rank"] = (df["taxon_rank"] == "species").astype(int)
    df["log_user_obs"] = np.log1p(df["user_observations_count"].fillna(0))
    df["log_n_photos"] = np.log1p(df["n_photos"].fillna(0))
    df["iconic"] = df["taxon_iconic_taxon_name"].fillna("Unknown")
    return df.reset_index(drop=True)


def km_style_curve(age_days: np.ndarray, is_research: np.ndarray,
                   bin_edges: np.ndarray) -> pd.DataFrame:
    """Empirical P(research | age) by age bin. NOT a real KM curve — we don't
    have transition times — but it makes the same picture: cumulative share
    of obs that have 'survived' to research grade as a function of age."""
    rows = []
    for lo, hi in zip(bin_edges[:-1], bin_edges[1:]):
        m = (age_days >= lo) & (age_days < hi)
        n = int(m.sum())
        if n == 0:
            continue
        p = float(is_research[m].mean())
        se = float(np.sqrt(p * (1 - p) / n)) if n > 0 else np.nan
        rows.append({
            "age_lo_days": lo, "age_hi_days": hi,
            "n": n, "p_research": p, "se": se,
            "ci_lo": max(0, p - 1.96 * se), "ci_hi": min(1, p + 1.96 * se),
        })
    return pd.DataFrame(rows)


def fit_logistic_iterative(X: np.ndarray, y: np.ndarray,
                           ridge: float = 1e-4, n_iter: int = 50) -> np.ndarray:
    """IRLS for logistic regression. Adds a tiny ridge for numeric stability.
    Avoids dragging in statsmodels for a one-shot script."""
    n, p = X.shape
    beta = np.zeros(p)
    for _ in range(n_iter):
        eta = X @ beta
        eta = np.clip(eta, -30, 30)
        mu = 1.0 / (1.0 + np.exp(-eta))
        W = mu * (1 - mu) + 1e-9
        z = eta + (y - mu) / W
        XtW = X.T * W
        H = XtW @ X + ridge * np.eye(p)
        g = XtW @ z
        new_beta = np.linalg.solve(H, g)
        if np.max(np.abs(new_beta - beta)) < 1e-7:
            beta = new_beta
            break
        beta = new_beta
    return beta


def cluster_robust_se(X: np.ndarray, y: np.ndarray, beta: np.ndarray,
                      cluster: np.ndarray, ridge: float = 1e-4) -> np.ndarray:
    """CR1-style cluster-robust SEs by `cluster` (observer). Standard sandwich:
    V = (X'WX)^-1 [Σ_g X_g' u_g u_g' X_g] (X'WX)^-1 with u_g = Σ_i (y - mu)*x.
    """
    eta = np.clip(X @ beta, -30, 30)
    mu = 1.0 / (1.0 + np.exp(-eta))
    W = mu * (1 - mu) + 1e-9
    bread = np.linalg.inv(X.T * W @ X + ridge * np.eye(X.shape[1]))
    resid = (y - mu)
    meat = np.zeros((X.shape[1], X.shape[1]))
    df = pd.DataFrame({"r": resid, "g": cluster})
    for g, idx in df.groupby("g").groups.items():
        idx = np.asarray(idx)
        u = (X[idx].T * resid[idx]).sum(axis=1)
        meat += np.outer(u, u)
    G = len(np.unique(cluster))
    n, k = X.shape
    adj = (G / (G - 1)) * ((n - 1) / (n - k)) if G > 1 else 1.0
    V = adj * bread @ meat @ bread
    return np.sqrt(np.diag(V))


def design(df: pd.DataFrame, drop_iconic: str = "Plantae") -> tuple[np.ndarray, list[str]]:
    """Numeric design + iconic dummies (drop one to avoid collinearity)."""
    base = pd.DataFrame({
        "Intercept": 1.0,
        "log_age_days": df["log_age"],
        "log_n_photos": df["log_n_photos"],
        "identifications_count": df["identifications_count"],
        "log_user_obs": df["log_user_obs"],
        "captive": df["captive_bin"],
        "species_rank": df["species_rank"],
    })
    iconic_dummies = pd.get_dummies(df["iconic"], prefix="iconic", drop_first=False)
    if f"iconic_{drop_iconic}" in iconic_dummies.columns:
        iconic_dummies = iconic_dummies.drop(columns=[f"iconic_{drop_iconic}"])
    keep = [c for c in iconic_dummies.columns
            if iconic_dummies[c].sum() >= 30]  # rare iconic groups dropped
    iconic_dummies = iconic_dummies[keep]
    X = pd.concat([base, iconic_dummies.astype(float)], axis=1)
    return X.to_numpy(dtype=float), list(X.columns)


def fit_and_report(df: pd.DataFrame, label: str) -> pd.DataFrame:
    X, names = design(df)
    y = df["is_research"].to_numpy(float)
    beta = fit_logistic_iterative(X, y)
    se = cluster_robust_se(X, y, beta, df["user_login"].to_numpy())
    z = beta / se
    or_ = np.exp(beta)
    ci_lo = np.exp(beta - 1.96 * se)
    ci_hi = np.exp(beta + 1.96 * se)
    out = pd.DataFrame({
        "term": names, "beta": beta, "se_cluster": se, "z": z,
        "OR": or_, "OR_ci_lo": ci_lo, "OR_ci_hi": ci_hi,
    })
    print(f"\n=== {label} ===  (n={len(df):,}, P(research)={y.mean():.3f}, "
          f"clusters={df['user_login'].nunique():,})")
    print(out.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    return out


def main() -> None:
    df = load()
    print(f"loaded {len(df):,} obs", file=sys.stderr)

    # ---- view A: research vs everything ----
    print("\n=== quality_grade composition ===")
    print(df["quality_grade"].value_counts())

    print("\n=== P(research) by iconic taxon (≥30 obs) ===")
    g = (df.groupby("iconic")
         .agg(n=("is_research", "size"), p_research=("is_research", "mean"))
         .query("n >= 30")
         .sort_values("p_research", ascending=False))
    print(g.to_string(float_format=lambda x: f"{x:.3f}"))

    # 'survival' curve: P(research) as a function of obs age
    edges = np.array([0, 30, 90, 180, 365, 730, 1825, 3650, 100000])
    curve = km_style_curve(df["age_days"].to_numpy(), df["is_research"].to_numpy(), edges)
    print("\n=== P(research-grade) by obs age bucket ===")
    print(curve.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    # ---- regressions ----
    out_full = fit_and_report(df, "View A: research vs (casual ∪ needs_id)")
    out_cond = fit_and_report(df[df["is_casual"] == 0].copy(),
                              "View B: research vs needs_id only (excl. casual)")

    # ---- write outputs ----
    notes = ROOT / "notes"
    notes.mkdir(exist_ok=True)
    curve.to_csv(notes / "survival_age_curve.csv", index=False)
    out_full.to_csv(notes / "survival_logit_full.csv", index=False)
    out_cond.to_csv(notes / "survival_logit_conditional.csv", index=False)
    print(f"\nwrote {notes / 'survival_age_curve.csv'}", file=sys.stderr)
    print(f"wrote {notes / 'survival_logit_full.csv'}", file=sys.stderr)
    print(f"wrote {notes / 'survival_logit_conditional.csv'}", file=sys.stderr)


if __name__ == "__main__":
    main()
