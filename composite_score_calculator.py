#!/usr/bin/env python3
"""
Composite Evolutionary Rate Score Calculator
=============================================

Computes a composite evolutionary-risk score per genome by combining:
  - published molecular-clock substitution rate (per site per year)
  - BV-BRC specialty-gene AMR count
  - BV-BRC AMR phenotypic resistance fraction

Three normalization modes are supported via --normalization-mode:

  absolute  (default)
      Fixed log-space bounds for rate (1e-9 to 1e-4) and a fixed AMR gene
      cap. Missing resistance -> 0.5. Scores are comparable across runs.

  batch
      Log-space min-max for rate, linear min-max for AMR, min-max for
      resistance fraction; all computed within the current run. Missing
      values default to 0.5. Not comparable across runs.

  legacy
      Reproduces the earlier prototype: linear min-max on raw
      rate_per_site_per_year, linear min-max on raw amr_genes, and raw
      resistance_fraction with missing -> 0.0. Only meaningful within a
      single run; unstable for very small batches (n < ~5).
"""

import sys
import argparse
import numpy as np
import pandas as pd
from collections import Counter, OrderedDict
import requests

from parser import parse_tsv, parse_organism
from rate_database import DynamicRateDatabase
from bvbrc_api import fetch_genome_metadata, fetch_amr_phenotypes, fetch_specialty_genes
from utils import save_output

# =============================================================================
# CONFIGURATION
# =============================================================================

W_RATE = 0.50
W_AMR = 0.25
W_RES = 0.25

BV_BRC_BASE = "https://www.bv-brc.org/api"

# Absolute-mode bounds
RATE_SCORE_MIN = 1e-9
RATE_SCORE_MAX = 1e-4
AMR_GENES_SCORE_MAX = 100.0


# =============================================================================
# BV-BRC GENOME ID RESOLUTION
# =============================================================================

def _query_first_genome_id(rql_filter):
    try:
        r = requests.get(
            f"{BV_BRC_BASE}/genome/?{rql_filter}&select(genome_id)&limit(1)",
            headers={"Accept": "application/json"},
            timeout=20,
        )
        r.raise_for_status()
        d = r.json()
        if d:
            return d[0].get("genome_id")
    except Exception:
        return None
    return None


def resolve_genome_ids(items):
    print("=" * 70)
    print("RESOLVING BV-BRC GENOME IDs")
    print("=" * 70)
    for it in items:
        gid = it.get("genome_id")
        it["resolution_source"] = ""
        it["proxy_used"] = False
        it["skip_reason"] = ""
        it["bvbrc_available"] = False

        if gid:
            it["resolution_source"] = "pre_supplied_genome_id"
            it["bvbrc_available"] = True
            continue

        asm = (it.get("assembly_accession") or "").strip()
        taxon_id = (it.get("taxon_id") or "").strip()
        species_name = (it.get("species_name") or "").strip()
        genus = (it.get("genus") or "").strip()

        if asm:
            q = requests.utils.quote(asm, safe="._-")
            gid = _query_first_genome_id(f"eq(assembly_accession,{q})")
            if gid:
                it["resolution_source"] = "assembly_accession"

        if not gid and taxon_id:
            tid = taxon_id.strip()
            if tid.isdigit():
                gid = _query_first_genome_id(f"eq(taxon_id,{tid})")
                if gid:
                    it["resolution_source"] = "ncbi_taxon_id"
                    it["proxy_used"] = True

        if not gid and species_name:
            q = requests.utils.quote(species_name, safe="._-")
            gid = _query_first_genome_id(f"keyword({q})")
            if gid:
                it["resolution_source"] = "species_keyword"
                it["proxy_used"] = True

        if not gid and genus:
            q = requests.utils.quote(genus, safe="._-")
            gid = _query_first_genome_id(f"keyword({q})")
            if gid:
                it["resolution_source"] = "genus_keyword"
                it["proxy_used"] = True

        if gid:
            it["genome_id"] = gid
            it["bvbrc_available"] = True
        else:
            it["skip_reason"] = "no_bvbrc_genome_id_match"
            it["resolution_source"] = "unresolved"
            print(
                f"  ⚠ Unresolved in BV-BRC: "
                f"{it.get('input_label') or it.get('species_name') or 'unknown'}"
            )


# =============================================================================
# NORMALIZATION FUNCTIONS
# =============================================================================

# ---- absolute mode --------------------------------------------------------

def _normalize_rate_absolute(values):
    """Log-space scaling to fixed bounds [RATE_SCORE_MIN, RATE_SCORE_MAX]."""
    s = pd.to_numeric(values, errors="coerce")
    out = pd.Series(0.5, index=s.index, dtype=float)
    valid = s[(s.notna()) & (s > 0)]
    if valid.empty:
        return out
    log_lo = np.log10(RATE_SCORE_MIN)
    log_hi = np.log10(RATE_SCORE_MAX)
    clipped = valid.clip(lower=RATE_SCORE_MIN, upper=RATE_SCORE_MAX)
    out.loc[valid.index] = ((np.log10(clipped) - log_lo) / (log_hi - log_lo)).clip(0, 1)
    return out


def _normalize_amr_absolute(values):
    """Linear scaling against fixed AMR_GENES_SCORE_MAX."""
    s = pd.to_numeric(values, errors="coerce")
    out = pd.Series(0.5, index=s.index, dtype=float)
    valid = s.dropna()
    if valid.empty:
        return out
    out.loc[valid.index] = (valid / AMR_GENES_SCORE_MAX).clip(0, 1)
    return out


# ---- batch mode -----------------------------------------------------------

def _normalize_rate_batch(values):
    """Log-space min-max within the current batch (clipped to absolute bounds first)."""
    s = pd.to_numeric(values, errors="coerce")
    out = pd.Series(0.5, index=s.index, dtype=float)
    valid = s[(s.notna()) & (s > 0)]
    if valid.empty:
        return out
    clipped = valid.clip(lower=RATE_SCORE_MIN, upper=RATE_SCORE_MAX)
    log_valid = np.log10(clipped)
    lo, hi = log_valid.min(), log_valid.max()
    if lo == hi:
        out.loc[valid.index] = 0.5
        return out
    out.loc[valid.index] = ((log_valid - lo) / (hi - lo)).clip(0, 1)
    return out


def _normalize_amr_batch(values):
    """Linear min-max within the current batch."""
    s = pd.to_numeric(values, errors="coerce")
    out = pd.Series(0.5, index=s.index, dtype=float)
    valid = s[(s.notna()) & (s >= 0)]
    if valid.empty:
        return out
    lo, hi = valid.min(), valid.max()
    if lo == hi:
        out.loc[valid.index] = 0.5
        return out
    out.loc[valid.index] = ((valid - lo) / (hi - lo)).clip(0, 1)
    return out


def _normalize_resistance_batch(values):
    """Min-max on resistance fraction within the current batch; missing -> 0.5."""
    s = pd.to_numeric(values, errors="coerce").clip(0, 1)
    out = pd.Series(0.5, index=s.index, dtype=float)
    valid = s.dropna()
    if valid.empty:
        return out
    lo, hi = valid.min(), valid.max()
    if lo == hi:
        out.loc[valid.index] = 0.5
        return out
    out.loc[valid.index] = ((valid - lo) / (hi - lo)).clip(0, 1)
    return out


# ---- legacy mode (reproduces prototype behavior) --------------------------

def _normalize_rate_legacy(values):
    """Linear min-max on raw rate_per_site_per_year. No clipping, no log transform."""
    s = pd.to_numeric(values, errors="coerce")
    out = pd.Series(np.nan, index=s.index, dtype=float)
    valid = s.dropna()
    if valid.empty:
        return out
    lo, hi = valid.min(), valid.max()
    if lo == hi:
        out.loc[valid.index] = 0.0
        return out
    out.loc[valid.index] = (valid - lo) / (hi - lo)
    return out


def _normalize_amr_legacy(values):
    """Linear min-max on raw amr_genes counts."""
    s = pd.to_numeric(values, errors="coerce")
    out = pd.Series(np.nan, index=s.index, dtype=float)
    valid = s.dropna()
    if valid.empty:
        return out
    lo, hi = valid.min(), valid.max()
    if lo == hi:
        out.loc[valid.index] = 0.0
        return out
    out.loc[valid.index] = (valid - lo) / (hi - lo)
    return out


# =============================================================================
# COMPOSITE SCORE
# =============================================================================

def compute_scores(items, gdata, adata, sdata, rate_db,
                   strict_bvbrc=False, normalization_mode="absolute"):
    print("=" * 70)
    print(f"STEP 4: Computing composite score (mode={normalization_mode})")
    print("=" * 70)

    rows = []
    for it in items:
        gid = it.get("genome_id")
        has_bvbrc = bool(gid and gid in gdata)
        if strict_bvbrc and not has_bvbrc:
            continue

        if has_bvbrc:
            meta = gdata[gid]
            amr = adata.get(gid, [])
            spec = sdata.get(gid, [])
            gname = meta.get("genome_name", "")
            glen = meta.get("genome_length", np.nan)
            genus, species = parse_organism(gname)
            if not genus:
                genus, species = it.get("genus"), it.get("species")

            sc = Counter(r.get("property", "") for r in spec)
            amr_g = sc.get("Antibiotic Resistance", 0)
            rc = sum(1 for r in amr if r.get("resistant_phenotype") == "Resistant")
            tot = sum(
                1 for r in amr
                if r.get("resistant_phenotype") in ["Resistant", "Susceptible"]
            )
            rfrac = rc / tot if tot else np.nan
        else:
            gname = ""
            glen = np.nan
            genus, species = it.get("genus"), it.get("species")
            amr_g = np.nan
            rfrac = np.nan

        rates, mlevel = rate_db.lookup(genus, species)
        rps = rates["rate_mid"]
        label = (
            f"{genus} {species or ''}".strip()
            if genus
            else (it.get("species_name") or it.get("input_label") or "")
        )

        rows.append({
            "genome_id": gid or "",
            "label": label,
            "input_label": it.get("input_label", ""),
            "species_name_input": it.get("species_name", ""),
            "mapping_level": it.get("mapping_level", ""),
            "mapping_description": it.get("mapping_description", ""),
            "strain": it.get("strain", ""),
            "selection_method": it.get("selection_method", ""),
            "taxon_id": it.get("taxon_id", ""),
            "resolution_source": it.get("resolution_source", ""),
            "bvbrc_available": has_bvbrc,
            "proxy_used": bool(it.get("proxy_used")),
            "skip_reason": it.get("skip_reason", ""),
            "amr_genes": amr_g,
            "resistance_fraction": rfrac,
            "rate_per_site_per_year": rps,
            "rate_match_level": mlevel,
            "genome_length_bp": glen,
            "reference": rates["reference"],
            "bvbrc_genome_name": gname,
        })

    df = pd.DataFrame(rows)
    if df.empty:
        return df

    # Warn on tiny batches when using within-batch modes
    if normalization_mode in ("batch", "legacy") and len(df) < 5:
        print(
            f"  ⚠ {normalization_mode!r} normalization requested with only "
            f"{len(df)} record(s); within-batch scores are unstable for small n. "
            f"Consider --normalization-mode absolute."
        )

    if normalization_mode == "batch":
        rate_score = _normalize_rate_batch(df["rate_per_site_per_year"])
        amr_score = _normalize_amr_batch(df["amr_genes"])
        resistance_score = _normalize_resistance_batch(df["resistance_fraction"])
    elif normalization_mode == "absolute":
        rate_score = _normalize_rate_absolute(df["rate_per_site_per_year"])
        amr_score = _normalize_amr_absolute(df["amr_genes"])
        resistance_score = pd.to_numeric(
            df["resistance_fraction"], errors="coerce"
        ).fillna(0.5).clip(0, 1)
    elif normalization_mode == "legacy":
        # Reproduces the earlier prototype exactly:
        #   linear min-max on raw rate and amr,
        #   raw resistance_fraction with missing -> 0.0.
        rate_score = _normalize_rate_legacy(df["rate_per_site_per_year"])
        amr_score = _normalize_amr_legacy(df["amr_genes"])
        resistance_score = pd.to_numeric(
            df["resistance_fraction"], errors="coerce"
        ).fillna(0.0).clip(0, 1)
    else:
        raise ValueError(
            f"Unsupported normalization_mode: {normalization_mode!r}. "
            "Expected 'absolute', 'batch', or 'legacy'."
        )

    # Persist per-component scores for downstream inspection
    df["rate_score"] = rate_score
    df["amr_score"] = amr_score
    df["resistance_score"] = resistance_score
    df["normalization_mode"] = normalization_mode
    df["composite_score"] = (
        W_RATE * rate_score
        + W_AMR * amr_score
        + W_RES * resistance_score
    )

    return df


# =============================================================================
# MAIN
# =============================================================================

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--tsv", help="Input manifest file")
    p.add_argument("--output", default="evolutionary_rates_quantified.csv")
    p.add_argument(
        "--strict-bvbrc", action="store_true",
        help="Only include genomes with BV-BRC metadata/feature support",
    )
    p.add_argument(
        "--normalization-mode",
        choices=["absolute", "batch", "legacy"],
        default="absolute",
        help=(
            "Composite score normalization mode. "
            "'absolute' uses fixed bounds (cross-run comparable, default). "
            "'batch' uses log-space min-max within the current run "
            "(within-run only; missing -> 0.5). "
            "'legacy' uses linear min-max within the current run "
            "(reproduces the earlier prototype; missing resistance -> 0.0)."
        ),
    )
    args = p.parse_args()

    if not args.tsv:
        print("Error: Please provide a structural dataset file input mapping using --tsv")
        sys.exit(1)

    items = parse_tsv(args.tsv)
    resolve_genome_ids(items)

    gids = list(OrderedDict.fromkeys(i["genome_id"] for i in items if i.get("genome_id")))
    if not gids and args.strict_bvbrc:
        print("No BV-BRC genome IDs resolved in strict mode.")
        sys.exit(1)

    if gids:
        gdata = fetch_genome_metadata(gids)
        adata = fetch_amr_phenotypes(gids)
        sdata = fetch_specialty_genes(gids)
    else:
        print("No BV-BRC genome IDs resolved; continuing with non-BV-BRC-compatible mode.")
        gdata, adata, sdata = {}, {}, {}

    rate_db = DynamicRateDatabase()

    df = compute_scores(
        items, gdata, adata, sdata, rate_db,
        strict_bvbrc=args.strict_bvbrc,
        normalization_mode=args.normalization_mode,
    )
    save_output(df, items, args.output)

    print("\nProcessing Pipeline Complete.\n")


if __name__ == "__main__":
    main()
