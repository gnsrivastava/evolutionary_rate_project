#!/usr/bin/env python3
import re
import pandas as pd

_COL_ALIASES = {
    "organism": ["acdb_organism","organism","organism_name","strain","species","name"],
    "species":  ["species_query","species_name","scientific_name"],
    "assembly": ["assembly_accession","assembly","accession"],
    "taxon_id": ["ncbi_taxon_id","taxon_id","tax_id"],
    "genome_id":["genome_id","bvbrc_id","patric_id"],
}

def _find_col(cols, key):
    lc = {c.lower().strip():c for c in cols}
    for a in _COL_ALIASES.get(key,[]):
        if a.lower() in lc: return lc[a.lower()]
    return None

def parse_organism(name):
    if not name: return None, None
    parts = name.strip().split()
    genus = parts[0] if parts else None
    species = parts[1] if len(parts)>1 else None
    if species and species.lower() in ("sp.","spp.","cf.","subsp."):
        species = None
    return genus, species

def parse_tsv(path):
    sep = "\t" if path.lower().endswith(".tsv") else ","
    df = pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False)
    c_org = _find_col(df.columns,"organism")
    c_sp = _find_col(df.columns,"species")
    c_asm = _find_col(df.columns,"assembly")
    c_gid = _find_col(df.columns,"genome_id")
    
    if not c_org and not c_sp and not c_gid:
        c_org = df.columns[0]
        
    items = []
    for _, row in df.iterrows():
        org = (row.get(c_org,"") if c_org else "").strip()
        sp = (row.get(c_sp,"") if c_sp else "").strip()
        asm = (row.get(c_asm,"") if c_asm else "").strip()
        gid = (row.get(c_gid,"") if c_gid else "").strip()
        
        best = sp or org
        genus, species = parse_organism(best)
        items.append({
            "input_label": org or sp, "species_name": best,
            "genus": genus, "species": species, "assembly_accession": asm,
            "qualifiers": [], "genome_id": gid, "taxon_id": ""
        })
    return items
