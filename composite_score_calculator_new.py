#!/usr/bin/env python3
"""
Composite Evolutionary Rate Score Calculator
=============================================

Computes a composite evolutionary-risk score per genome by combining:
  - published molecular-clock substitution rate (per site per year)
  - BV-BRC specialty-gene AMR count
  - BV-BRC AMR phenotypic resistance fraction

Normalization modes (via --normalization-mode):

  absolute  (default)   Fixed log-space bounds for rate, fixed AMR cap.
                        Missing resistance -> 0.5. Cross-run comparable.

  batch                 Log-space min-max for rate, linear min-max for AMR
                        and resistance, all within the current run.
                        Missing values default to 0.5.

  legacy                Linear min-max on raw rate and amr_genes, raw
                        resistance_fraction with missing -> 0.0. Reproduces
                        the earlier prototype exactly. Output columns and
                        their order match the prototype CSV schema:
                          genome_id, label, genome_name, genome_length_bp,
                          gc_content, cds_count, amr_genes, virulence_genes,
                          resistant_count, susceptible_count,
                          resistance_fraction, rate_per_site_per_year,
                          rate_low, rate_high, snps_per_genome_per_year,
                          norm_rate_minmax, norm_abs_rate_minmax,
                          composite_score, reference

Phenotype-aware rate lookup:
    Manifest phenotype qualifiers (MRSA / MSSA / MRSE / MSSE / VRE /
    VRSA / ESBL / CRE) are inferred from an explicit 'phenotype' column
    or from 'strain' / 'selection_method' / 'input_label' / 'species_name',
    and passed to DynamicRateDatabase.lookup(genus, species, qualifier).
"""

import sys
import re
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

RATE_SCORE_MIN = 1e-9
RATE_SCORE_MAX = 1e-4
AMR_GENES_SCORE_MAX = 100.0

KNOWN_PHENOTYPE_TOKENS = ["MRSA", "MSSA", "MRSE", "MSSE", "VRE", "VRSA", "ESBL", "CRE"]
_PHENOTYPE_REGEX = re.compile(
    r"\b(" + "|".join(KNOWN_PHENOTYPE_TOKENS) + r")\b",
    re.IGNORECASE,
)

# Prototype output schema - used only in legacy mode
LEGACY_OUTPUT_COLUMNS = [
    "genome_id", "label", "genome_name", "genome_length_bp", "gc_content",
    "cds_count", "amr_genes", "virulence_genes", "resistant_count",
    "susceptible_count", "resistance_fraction", "rate_per_site_per_year",
    "rate_low", "rate_high", "snps_per_genome_per_year",
    "norm_rate_minmax", "norm_abs_rate_minmax", "composite_score", "reference",
]


# =============================================================================
# PHENOTYPE INFERENCE
# =============================================================================

def infer_phenotype(item):
    for key in ("phenotype", "qualifier"):
        val = (item.get(key) or "").strip()
        if val:
            return val.upper()
    for key in ("strain", "selection_method", "input_label", "species_name",
                "acdb_strain_designation", "ncbi_strain"):
        val = item.get(key) or ""
        m = _PHENOTYPE_REGEX.search(str(val))
        if m:
            return m.group(1).upper()
    return None


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

def _normalize_rate_absolute(values):
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
    s = pd.to_numeric(values, errors="coerce")
    out = pd.Series(0.5, index=s.index, dtype=float)
    valid = s.dropna()
    if valid.empty:
        return out
    out.loc[valid.index] = (valid / AMR_GENES_SCORE_MAX).clip(0, 1)
    return out


def _normalize_rate_batch(values):
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


def _minmax_linear(values, fill=np.nan):
    """Linear min-max used by the prototype (legacy mode)."""
    s = pd.to_numeric(values, errors="coerce")
    out = pd.Series(fill, index=s.index, dtype=float)
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
# ROW ASSEMBLY
# =============================================================================

def _build_rows(items, gdata, adata, sdata, rate_db, strict_bvbrc, mode):
    """
    Assemble per-genome feature rows.

    In legacy mode the rows carry the extra prototype columns
    (gc_content, cds_count, virulence_genes, resistant_count,
    susceptible_count, rate_low, rate_high, snps_per_genome_per_year,
    genome_name). In other modes only the standard columns are populated.
    """
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
            gc = meta.get("gc_content", np.nan)
            cds = meta.get("patric_cds", np.nan)
            genus, species = parse_organism(gname)
            if not genus:
                genus, species = it.get("genus"), it.get("species")

            sc = Counter(r.get("property", "") for r in spec)
            amr_g = sc.get("Antibiotic Resistance", 0)
            vir_g = sc.get("Virulence Factor", 0) + sc.get("Virulance factor", 0)
            rc = sum(1 for r in amr if r.get("resistant_phenotype") == "Resistant")
            sc_count = sum(1 for r in amr if r.get("resistant_phenotype") == "Susceptible")
            tot = rc + sc_count
            if mode == "legacy":
                rfrac = rc / tot if tot else 0.0
            else:
                rfrac = rc / tot if tot else np.nan
        else:
            meta = {}
            gname = ""
            glen = np.nan
            gc = np.nan
            cds = np.nan
            genus, species = it.get("genus"), it.get("species")
            amr_g = np.nan
            vir_g = np.nan
            rc = 0
            sc_count = 0
            rfrac = np.nan if mode != "legacy" else 0.0

        phenotype = infer_phenotype(it)
        rates, mlevel = rate_db.lookup(genus, species, phenotype)
        rps = rates["rate_mid"]
        rlo = rates.get("rate_low", np.nan)
        rhi = rates.get("rate_high", np.nan)

        label_species = f"{genus} {species or ''}".strip() if genus else ""
        if phenotype and label_species:
            label = f"{label_species} {phenotype}"
        elif label_species:
            label = label_species
        else:
            label = it.get("species_name") or it.get("input_label") or ""

        row = {
            "genome_id": gid or "",
            "label": label,
            "input_label": it.get("input_label", ""),
            "species_name_input": it.get("species_name", ""),
            "phenotype": phenotype or "",
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
        }

        if mode == "legacy":
            row.update({
                "genome_name": gname,
                "gc_content": gc,
                "cds_count": cds,
                "virulence_genes": vir_g,
                "resistant_count": rc,
                "susceptible_count": sc_count,
                "rate_low": rlo,
                "rate_high": rhi,
                "snps_per_genome_per_year": (
                    rps * glen
                    if pd.notna(rps) and pd.notna(glen) else np.nan
                ),
            })

        rows.append(row)

    return pd.DataFrame(rows)


# =============================================================================
# COMPOSITE SCORE
# =============================================================================

def compute_scores(items, gdata, adata, sdata, rate_db,
                   strict_bvbrc=False, normalization_mode="absolute"):
    print("=" * 70)
    print(f"STEP 4: Computing composite score (mode={normalization_mode})")
    print("=" * 70)

    df = _build_rows(items, gdata, adata, sdata, rate_db,
                     strict_bvbrc, normalization_mode)
    if df.empty:
        return df

    if normalization_mode in ("batch", "legacy") and len(df) < 5:
        print(
            f"  ⚠ {normalization_mode!r} normalization with only {len(df)} "
            f"record(s) is unstable; consider --normalization-mode absolute."
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
        # Prototype behavior: linear min-max on rate_per_site_per_year
        # (norm_rate_minmax) and on snps_per_genome_per_year
        # (norm_abs_rate_minmax). Composite uses norm_rate_minmax for the
        # rate component, linear min-max on amr_genes for the AMR component,
        # and raw resistance_fraction (missing -> 0.0) for the resistance
        # component.
        rate_score = _minmax_linear(df["rate_per_site_per_year"])
        abs_rate_score = _minmax_linear(df["snps_per_genome_per_year"])
        amr_score = _minmax_linear(df["amr_genes"])
        resistance_score = pd.to_numeric(
            df["resistance_fraction"], errors="coerce"
        ).fillna(0.0).clip(0, 1)

        df["norm_rate_minmax"] = rate_score
        df["norm_abs_rate_minmax"] = abs_rate_score
    else:
        raise ValueError(
            f"Unsupported normalization_mode: {normalization_mode!r}. "
            "Expected 'absolute', 'batch', or 'legacy'."
        )

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
# OUTPUT
# =============================================================================

def _save_legacy(df, filename):
    """Save the prototype's fixed 19-column schema in fixed order."""
    print("=" * 70)
    print(f"Saving legacy-schema output to {filename}")
    print("=" * 70)
    # Ensure every expected column exists (NaN if unavailable)
    for col in LEGACY_OUTPUT_COLUMNS:
        if col not in df.columns:
            df[col] = np.nan
    df[LEGACY_OUTPUT_COLUMNS].to_csv(filename, index=False)
    print(f"  ✓ Saved {len(df)} rows × {len(LEGACY_OUTPUT_COLUMNS)} columns")


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
            "'absolute' (default): fixed bounds, cross-run comparable. "
            "'batch': log-space min-max within run (missing -> 0.5). "
            "'legacy': linear min-max within run, prototype 19-column output "
            "schema (missing resistance -> 0.0)."
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

    if args.normalization_mode == "legacy":
        _save_legacy(df, args.output)
    else:
        save_output(df, items, args.output)

    print("\nProcessing Pipeline Complete.\n")


if __name__ == "__main__":
    main()
