"""Pull all observations from the iNaturalist project
'biota-of-university-of-auckland-city-campus' and save raw JSON + flat CSV.

Run: python3 scripts/01_pull_inat.py
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import pandas as pd
import requests

PROJECT_SLUG = "biota-of-university-of-auckland-city-campus"
API = "https://api.inaturalist.org/v1/observations"
PER_PAGE = 200
USER_AGENT = "biosci738-significance-prep/0.1 (jwin74@gmail.com)"

ROOT = Path(__file__).resolve().parent.parent
RAW_DIR = ROOT / "data" / "raw"
CLEAN_DIR = ROOT / "data" / "clean"
RAW_JSON = RAW_DIR / "observations_raw.json"
CLEAN_CSV = CLEAN_DIR / "observations_flat.csv"


def fetch_all() -> list[dict]:
    """Cursor-paginate via id_above to avoid the 10k page*per_page cap."""
    session = requests.Session()
    session.headers["User-Agent"] = USER_AGENT

    all_results: list[dict] = []
    id_above = 0
    page_num = 0
    while True:
        page_num += 1
        params = {
            "project_id": PROJECT_SLUG,
            "per_page": PER_PAGE,
            "order_by": "id",
            "order": "asc",
            "id_above": id_above,
        }
        r = session.get(API, params=params, timeout=60)
        r.raise_for_status()
        payload = r.json()
        results = payload.get("results", [])
        if not results:
            break
        all_results.extend(results)
        id_above = results[-1]["id"]
        total = payload.get("total_results", "?")
        print(f"  page {page_num}: +{len(results)}  (cumulative {len(all_results)}/{total})")
        if len(results) < PER_PAGE:
            break
        time.sleep(1.0)
    return all_results


def flatten(obs: dict) -> dict:
    taxon = obs.get("taxon") or {}
    user = obs.get("user") or {}
    geo = obs.get("geojson") or {}
    coords = geo.get("coordinates") or [None, None]
    photos = obs.get("photos") or []
    place_ids = obs.get("place_ids") or []
    return {
        "id": obs.get("id"),
        "uuid": obs.get("uuid"),
        "observed_on": obs.get("observed_on"),
        "observed_on_details_year": (obs.get("observed_on_details") or {}).get("year"),
        "observed_on_details_month": (obs.get("observed_on_details") or {}).get("month"),
        "observed_on_details_day": (obs.get("observed_on_details") or {}).get("day"),
        "observed_on_details_hour": (obs.get("observed_on_details") or {}).get("hour"),
        "observed_on_details_week": (obs.get("observed_on_details") or {}).get("week"),
        "observed_on_details_day_of_week": (obs.get("observed_on_details") or {}).get("day"),
        "created_at": obs.get("created_at"),
        "time_observed_at": obs.get("time_observed_at"),
        "time_zone": obs.get("time_zone"),
        "quality_grade": obs.get("quality_grade"),
        "captive": obs.get("captive"),
        "license_code": obs.get("license_code"),
        "geoprivacy": obs.get("geoprivacy"),
        "obscured": obs.get("obscured"),
        "positional_accuracy": obs.get("positional_accuracy"),
        "public_positional_accuracy": obs.get("public_positional_accuracy"),
        "longitude": coords[0],
        "latitude": coords[1],
        "place_guess": obs.get("place_guess"),
        "user_id": user.get("id"),
        "user_login": user.get("login"),
        "user_name": user.get("name"),
        "user_orcid": user.get("orcid"),
        "user_observations_count": user.get("observations_count"),
        "identifications_count": obs.get("identifications_count"),
        "num_identification_agreements": obs.get("num_identification_agreements"),
        "num_identification_disagreements": obs.get("num_identification_disagreements"),
        "comments_count": obs.get("comments_count"),
        "faves_count": obs.get("faves_count"),
        "taxon_id": taxon.get("id"),
        "taxon_name": taxon.get("name"),
        "taxon_rank": taxon.get("rank"),
        "taxon_iconic_taxon_name": taxon.get("iconic_taxon_name"),
        "taxon_preferred_common_name": taxon.get("preferred_common_name"),
        "taxon_threatened": taxon.get("threatened"),
        "taxon_introduced": taxon.get("introduced"),
        "taxon_native": taxon.get("native"),
        "taxon_endemic": taxon.get("endemic"),
        "n_photos": len(photos),
        "n_place_ids": len(place_ids),
        "uri": obs.get("uri"),
    }


def main() -> None:
    print(f"Pulling project '{PROJECT_SLUG}' from {API}")
    obs = fetch_all()
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    CLEAN_DIR.mkdir(parents=True, exist_ok=True)
    RAW_JSON.write_text(json.dumps(obs))
    print(f"Wrote {len(obs)} raw observations -> {RAW_JSON} ({RAW_JSON.stat().st_size/1e6:.1f} MB)")

    df = pd.DataFrame([flatten(o) for o in obs])
    df.to_csv(CLEAN_CSV, index=False)
    print(f"Wrote flat CSV -> {CLEAN_CSV} ({len(df)} rows, {len(df.columns)} cols)")


if __name__ == "__main__":
    main()
