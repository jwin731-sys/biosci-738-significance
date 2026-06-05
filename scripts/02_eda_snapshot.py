"""Lightweight EDA snapshot for the iNat City Campus pull.
Writes a markdown report to notes/eda_snapshot.md.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DF_PATH = ROOT / "data" / "clean" / "observations_flat.csv"
OUT = ROOT / "notes" / "eda_snapshot.md"

df = pd.read_csv(DF_PATH, low_memory=False)
df["observed_on"] = pd.to_datetime(df["observed_on"], errors="coerce")
df["created_at"] = pd.to_datetime(df["created_at"], errors="coerce", utc=True)
df["year"] = df["observed_on"].dt.year
df["month"] = df["observed_on"].dt.month
df["yyyy_mm"] = df["observed_on"].dt.to_period("M").astype(str)


def section(title: str) -> str:
    return f"\n## {title}\n\n"


lines: list[str] = ["# iNat City Campus — EDA snapshot", ""]
lines.append(f"Total observations: **{len(df):,}**  ")
lines.append(f"Unique observers: **{df['user_login'].nunique():,}**  ")
lines.append(f"Unique taxa (any rank): **{df['taxon_id'].nunique():,}**  ")
lines.append(f"Date range observed_on: **{df['observed_on'].min().date()} → {df['observed_on'].max().date()}**  ")

lines.append(section("Quality grade"))
qg = df["quality_grade"].value_counts(dropna=False)
lines.append(qg.to_markdown())

lines.append(section("Iconic taxon (broad group)"))
ic = df["taxon_iconic_taxon_name"].value_counts(dropna=False)
lines.append(ic.to_markdown())

lines.append(section("Top 20 taxa by observation count"))
top_taxa = (
    df.groupby(["taxon_name", "taxon_iconic_taxon_name", "taxon_rank"])
    .size()
    .sort_values(ascending=False)
    .head(20)
    .reset_index(name="n")
)
lines.append(top_taxa.to_markdown(index=False))

lines.append(section("Observer effort distribution (Pareto check)"))
obs_per_user = df.groupby("user_login").size().sort_values(ascending=False)
lines.append(f"- Top observer: **{obs_per_user.index[0]}** with **{int(obs_per_user.iloc[0])}** observations")
lines.append(f"- Top 10 observers contribute **{int(obs_per_user.head(10).sum()):,}** of {len(df):,} ({obs_per_user.head(10).sum()/len(df)*100:.1f}%)")
lines.append(f"- Median obs/user: **{int(obs_per_user.median())}**, mean: **{obs_per_user.mean():.1f}**")
lines.append(f"- Single-observation observers: **{int((obs_per_user == 1).sum())}** ({(obs_per_user == 1).mean()*100:.1f}% of users)")

lines.append("\n**Top 15 observers:**\n")
lines.append(obs_per_user.head(15).reset_index(name="n").to_markdown(index=False))

lines.append(section("Observations by year"))
yr = df.groupby("year").size().reset_index(name="n").sort_values("year")
lines.append(yr.to_markdown(index=False))

lines.append(section("Observations by month-of-year (seasonality)"))
mo = df.groupby("month").size().reset_index(name="n")
lines.append(mo.to_markdown(index=False))

lines.append(section("Spatial — bounding box & precision"))
lines.append(f"- Latitude range: **{df['latitude'].min():.5f} to {df['latitude'].max():.5f}**")
lines.append(f"- Longitude range: **{df['longitude'].min():.5f} to {df['longitude'].max():.5f}**")
lines.append(f"- Records with lat/lon: **{df[['latitude','longitude']].dropna().shape[0]:,}**")
lines.append(f"- Obscured: **{int(df['obscured'].sum())}** ({df['obscured'].mean()*100:.1f}%)")
acc = df["public_positional_accuracy"].dropna()
lines.append(f"- Public positional accuracy median: **{acc.median():.0f} m**, 90th pct: **{acc.quantile(0.9):.0f} m**, max: **{acc.max():.0f} m**")

lines.append(section("Native / introduced flagging"))
for col in ("taxon_native", "taxon_introduced", "taxon_endemic", "taxon_threatened"):
    s = df[col]
    lines.append(f"- {col}: True={int((s == True).sum())}, False={int((s == False).sum())}, NA={int(s.isna().sum())}")

lines.append(section("Identification structure (potential nesting)"))
lines.append(f"- Mean identifications per obs: **{df['identifications_count'].mean():.2f}**")
lines.append(f"- Obs with 0 identifications: **{int((df['identifications_count'] == 0).sum())}**")
lines.append(f"- Obs with disagreements: **{int((df['num_identification_disagreements'] > 0).sum())}**")

OUT.write_text("\n".join(str(x) for x in lines))
print(f"Wrote {OUT}")
print(f"  ({OUT.stat().st_size/1024:.1f} KB)")
