#!/usr/bin/env python3
import requests

BV_BRC_BASE = "https://www.bv-brc.org/api"

def fetch_genome_metadata(gids):
    print("="*70); print("STEP 1: Fetching genome metadata"); print("="*70)
    out = {}
    fields = "genome_id,genome_name,genome_length,contigs,gc_content,patric_cds,genome_status"
    for g in gids:
        try:
            r = requests.get(f"{BV_BRC_BASE}/genome/?eq(genome_id,{g})&select({fields})",
                             headers={"Accept":"application/json"}, timeout=30)
            r.raise_for_status(); d = r.json()
            if d: 
                out[g] = d[0]
                print(f"  ✓ {g}: {d[0].get('genome_name','N/A')} ({d[0].get('genome_length',0):,} bp)")
        except Exception as e: 
            print(f"  ✗ {g}: {e}")
    return out

def fetch_amr_phenotypes(gids):
    print("="*70); print("STEP 2: Fetching AMR phenotypes"); print("="*70)
    out = {}
    fields = "genome_id,antibiotic,resistant_phenotype"
    for g in gids:
        try:
            r = requests.get(f"{BV_BRC_BASE}/genome_amr/?eq(genome_id,{g})&select({fields})",
                             headers={"Accept":"application/json"}, timeout=30)
            r.raise_for_status(); d = r.json(); out[g] = d
            print(f"  ✓ {g}: {len(d)} AMR records fetched")
        except Exception as e: 
            print(f"  ✗ {g}: {e}"); out[g] = []
    return out

def fetch_specialty_genes(gids):
    print("="*70); print("STEP 3: Fetching specialty genes"); print("="*70)
    out = {}
    fields = "genome_id,property,gene"
    for g in gids:
        try:
            r = requests.get(f"{BV_BRC_BASE}/sp_gene/?eq(genome_id,{g})&select({fields})&limit(500)",
                             headers={"Accept":"application/json"}, timeout=30)
            r.raise_for_status(); d = r.json(); out[g] = d
            print(f"{g}: {len(d)} specialty genes resolved")
        except Exception as e: 
            print(f"  ✗ {g}: {e}"); out[g] = []
    return out
