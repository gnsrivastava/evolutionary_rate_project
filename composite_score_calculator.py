#!/usr/bin/env python3
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

W_RATE = 0.50; W_AMR = 0.25; W_RES = 0.25
BV_BRC_BASE = "https://www.bv-brc.org/api"

def resolve_genome_ids(items):
    print("="*70); print("RESOLVING BV-BRC GENOME IDs"); print("="*70)
    for it in items:
        if it.get("genome_id"): continue
        gid = None
        if it["assembly_accession"]:
            try:
                r = requests.get(f"{BV_BRC_BASE}/genome/?eq(assembly_accession,{it['assembly_accession']})&select(genome_id)&limit(1)", 
                                 headers={"Accept":"application/json"}, timeout=20)
                d = r.json()
                if d: gid = d[0].get("genome_id")
            except: pass
        if not gid and it["species_name"]:
            try:
                r = requests.get(f"{BV_BRC_BASE}/genome/?keyword({it['species_name']})&select(genome_id)&limit(1)", 
                                 headers={"Accept":"application/json"}, timeout=20)
                d = r.json()
                if d: gid = d[0].get("genome_id")
            except: pass
        if gid:
            it["genome_id"] = gid

def compute_scores(items, gdata, adata, sdata, rate_db):
    print("="*70); print("STEP 4: Computing composite score"); print("="*70)
    rows = []
    for it in items:
        gid = it.get("genome_id")
        if not gid or gid not in gdata: continue
        
        meta = gdata[gid]
        amr = adata.get(gid, [])
        spec = sdata.get(gid, [])
        
        gname = meta.get("genome_name", "")
        glen = meta.get("genome_length", 0)
        genus, species = parse_organism(gname)
        
        # Call the dynamic database engine logic block here
        rates, mlevel = rate_db.lookup(genus, species)
        
        sc = Counter(r.get("property","") for r in spec)
        amr_g = sc.get("Antibiotic Resistance", 0)
        
        rc = sum(1 for r in amr if r.get("resistant_phenotype") == "Resistant")
        tot = sum(1 for r in amr if r.get("resistant_phenotype") in ["Resistant", "Susceptible"])
        rfrac = rc / tot if tot else 0.0
        
        rps = rates["rate_mid"]
        
        rows.append({
            "genome_id": gid, "label": f"{genus} {species or ''}".strip(),
            "amr_genes": amr_g, "resistance_fraction": rfrac,
            "rate_per_site_per_year": rps, "rate_match_level": mlevel,
            "genome_length_bp": glen, "reference": rates["reference"]
        })
        
    df = pd.DataFrame(rows)
    if df.empty: return df
    
    def mm(a):
        return (a - a.min()) / (a.max() - a.min()) if a.max() != a.min() else np.full_like(a, 0.5)

    if len(df) >= 2:
        df["composite_score"] = (W_RATE * mm(df["rate_per_site_per_year"].values) +
                                 W_AMR * mm(df["amr_genes"].values.astype(float)) +
                                 W_RES * df["resistance_fraction"].values)
    else:
        df["composite_score"] = 0.5
        
    return df

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--tsv", help="Input manifest file")
    p.add_argument("--output", default="evolutionary_rates_quantified.csv")
    args = p.parse_args()
    
    if not args.tsv:
        print("Error: Please provide a structural dataset file input mapping using --tsv"); sys.exit(1)
        
    items = parse_tsv(args.tsv)
    resolve_genome_ids(items)
    
    gids = list(OrderedDict.fromkeys(i["genome_id"] for i in items if i.get("genome_id")))
    if not gids:
        print("No IDs resolved contextually."); sys.exit(1)
        
    gdata = fetch_genome_metadata(gids)
    adata = fetch_amr_phenotypes(gids)
    sdata = fetch_specialty_genes(gids)
    
    # Initialize the dynamic database orchestrator
    rate_db = DynamicRateDatabase()
    
    df = compute_scores(items, gdata, adata, sdata, rate_db)
    save_output(df, items, args.output)
    
    print("\nProcessing Pipeline Complete.\n")

if __name__ == "__main__":
    main()
