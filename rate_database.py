#!/usr/bin/env python3
from literature_search import search_evolutionary_rates
from rate_extractor import synthesize_rate_record
from rate_cache import load_cache, save_cache

# Your hardcoded published baselines serve as primary lookups
STATIC_RATE_DB = {
    ("Staphylococcus","aureus","MRSA"): {"rate_low":1.2e-6,"rate_mid":1.4e-6, "rate_high":1.57e-6,"reference":"Harris 2010 Science"},
    ("Staphylococcus","aureus","MSSA"): {"rate_low":1.1e-6,"rate_mid":1.3e-6, "rate_high":1.5e-6, "reference":"Holden 2013 Genome Res"},
    ("Staphylococcus","aureus",None):   {"rate_low":1.1e-6,"rate_mid":1.3e-6, "rate_high":1.57e-6,"reference":"Harris 2010; Holden 2013 combined"},
    ("Klebsiella","pneumoniae",None):   {"rate_low":3.0e-7,"rate_mid":3.65e-7,"rate_high":4.4e-7, "reference":"Wyres 2019 Nat Commun"},
    ("Escherichia","coli",None):        {"rate_low":2.0e-7,"rate_mid":4.4e-7, "rate_high":8.3e-7, "reference":"Duchene 2016 MBE"},
}

DEFAULT_RATE = {
    "rate_low":3.0e-7, "rate_mid":8.0e-7, "rate_high":2.0e-6,
    "reference":"Pan-bacterial default; Duchene 2016 MBE Table 1 median",
}

class DynamicRateDatabase:
    def __init__(self):
        self.static_db = STATIC_RATE_DB
        self.dynamic_cache = load_cache()

    def get_cache_key(self, genus, species, qualifier):
        return f"{genus or ''}||{species or ''}||{qualifier or ''}"

    def lookup(self, genus, species, qualifier=None):
        # 1. Attempt Static Hardcoded Matrix Lookups
        for key, level in [((genus, species, qualifier), "species+qualifier"),
                           ((genus, species, None), "species"),
                           ((genus, None, None), "genus")]:
            if key in self.static_db:
                return self.static_db[key], level

        # 2. Attempt Dynamic JSON Cache Lookups
        for g, s, q, lvl in [(genus, species, qualifier, "dynamic_species+qualifier"),
                             (genus, species, None, "dynamic_species"),
                             (genus, None, None, "dynamic_genus")]:
            ckey = self.get_cache_key(g, s, q)
            if ckey in self.dynamic_cache:
                return self.dynamic_cache[ckey], lvl

        # 3. Dynamic Live Literature Synthesis Execution Call
        if genus:
            print(f"Species missing from baseline profile. Querying dynamic literature APIs for: {genus} {species or ''}")
            articles = search_evolutionary_rates(genus, species)
            rate_rec = synthesize_rate_record(genus, species, articles)
            
            if rate_rec:
                ckey = self.get_cache_key(genus, species, qualifier)
                self.dynamic_cache[ckey] = rate_rec
                save_cache(self.dynamic_cache)
                return rate_rec, "dynamic_extracted"
                
        return DEFAULT_RATE, "default"
