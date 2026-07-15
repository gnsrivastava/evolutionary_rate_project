#!/usr/bin/env python3
from literature_search import search_evolutionary_rates
from rate_extractor import synthesize_rate_record
from rate_cache import load_cache, save_cache

# Your hardcoded published baselines serve as primary lookups
STATIC_RATE_DB = {
     ("Staphylococcus","aureus","MRSA"):      {"rate_low":1.2e-6,"rate_mid":1.4e-6, "rate_high":1.57e-6,"reference":"Harris 2010 Science; Nubel 2010 PNAS"},
    ("Staphylococcus","aureus","MSSA"):      {"rate_low":1.1e-6,"rate_mid":1.3e-6, "rate_high":1.5e-6, "reference":"Holden 2013 Genome Res"},
    ("Staphylococcus","aureus",None):         {"rate_low":1.1e-6,"rate_mid":1.3e-6, "rate_high":1.57e-6,"reference":"Harris 2010; Holden 2013 combined"},
    ("Staphylococcus","epidermidis","MRSE"):  {"rate_low":1.2e-6,"rate_mid":1.4e-6, "rate_high":1.6e-6, "reference":"Meric 2018 Nat Commun; Espadinha 2019"},
    ("Staphylococcus","epidermidis","MSSE"):  {"rate_low":1.0e-6,"rate_mid":1.2e-6, "rate_high":1.4e-6, "reference":"Meric 2018 Nat Commun"},
    ("Staphylococcus","epidermidis",None):    {"rate_low":1.0e-6,"rate_mid":1.3e-6, "rate_high":1.6e-6, "reference":"Meric 2018 Nat Commun"},
    ("Staphylococcus",None,None):             {"rate_low":1.0e-6,"rate_mid":1.3e-6, "rate_high":1.6e-6, "reference":"Staphylococcus genus avg"},
    ("Klebsiella","pneumoniae",None):         {"rate_low":3.0e-7,"rate_mid":3.65e-7,"rate_high":4.4e-7, "reference":"Wyres 2019 Nat Commun; Duchene 2016 MBE"},
    ("Klebsiella","oxytoca",None):            {"rate_low":3.0e-7,"rate_mid":4.0e-7, "rate_high":5.0e-7, "reference":"Enterobacterales range approx"},
    ("Klebsiella","aerogenes",None):          {"rate_low":3.0e-7,"rate_mid":4.0e-7, "rate_high":5.0e-7, "reference":"Enterobacterales range approx"},
    ("Klebsiella","ornithinolytica",None):    {"rate_low":3.0e-7,"rate_mid":4.0e-7, "rate_high":5.0e-7, "reference":"Enterobacterales range approx"},
    ("Klebsiella",None,None):                 {"rate_low":3.0e-7,"rate_mid":3.8e-7, "rate_high":5.0e-7, "reference":"Klebsiella genus (Wyres 2019)"},
    ("Pseudomonas","aeruginosa",None):        {"rate_low":5.0e-7,"rate_mid":6.5e-7, "rate_high":2.6e-6, "reference":"Dettman 2016; Marvig 2014"},
    ("Pseudomonas","testosteroni",None):      {"rate_low":5.0e-7,"rate_mid":7.0e-7, "rate_high":2.0e-6, "reference":"Pseudomonas genus approx"},
    ("Pseudomonas",None,None):                {"rate_low":5.0e-7,"rate_mid":6.5e-7, "rate_high":2.6e-6, "reference":"Pseudomonas genus (Dettman 2016)"},
    ("Escherichia","coli",None):              {"rate_low":2.0e-7,"rate_mid":4.4e-7, "rate_high":8.3e-7, "reference":"Duchene 2016 MBE; Reeves 2011"},
    ("Escherichia",None,None):                {"rate_low":2.0e-7,"rate_mid":4.4e-7, "rate_high":8.3e-7, "reference":"Escherichia genus (Duchene 2016)"},
    ("Enterococcus","faecium",None):          {"rate_low":3.0e-6,"rate_mid":5.0e-6, "rate_high":9.0e-6, "reference":"Raven 2016; Pinholt 2019"},
    ("Enterococcus","faecalis",None):         {"rate_low":1.0e-6,"rate_mid":2.5e-6, "rate_high":4.0e-6, "reference":"Pinholt 2019; Raven 2016"},
    ("Enterococcus",None,None):               {"rate_low":1.0e-6,"rate_mid":3.5e-6, "rate_high":9.0e-6, "reference":"Enterococcus genus avg"},
    ("Acinetobacter","baumannii",None):       {"rate_low":5.0e-7,"rate_mid":9.5e-7, "rate_high":2.5e-6, "reference":"Snitkin 2011 Sci Transl Med; Wright 2016"},
    ("Acinetobacter","calcoaceticus",None):   {"rate_low":5.0e-7,"rate_mid":9.0e-7, "rate_high":2.0e-6, "reference":"A. baumannii complex approx"},
    ("Acinetobacter",None,None):              {"rate_low":5.0e-7,"rate_mid":9.5e-7, "rate_high":2.5e-6, "reference":"Acinetobacter genus (Snitkin 2011)"},
    ("Mycobacterium","tuberculosis",None):    {"rate_low":1.0e-8,"rate_mid":5.0e-8, "rate_high":1.0e-7, "reference":"Ford 2011 PLoS Pathog; Walker 2013 Lancet"},
    ("Mycobacterium","abscessus",None):       {"rate_low":5.0e-8,"rate_mid":2.0e-7, "rate_high":5.0e-7, "reference":"Bryant 2016 Science approx"},
    ("Mycobacterium",None,None):              {"rate_low":1.0e-8,"rate_mid":1.0e-7, "rate_high":5.0e-7, "reference":"Mycobacterium genus avg"},
    ("Salmonella","enterica",None):           {"rate_low":1.5e-7,"rate_mid":3.4e-7, "rate_high":6.0e-7, "reference":"Okoro 2012 Nat Genet; Duchene 2016"},
    ("Salmonella","typhi",None):              {"rate_low":1.0e-7,"rate_mid":2.0e-7, "rate_high":4.0e-7, "reference":"Wong 2015 Nat Genet; Holt 2008"},
    ("Salmonella","typhimurium",None):        {"rate_low":2.0e-7,"rate_mid":3.5e-7, "rate_high":6.5e-7, "reference":"Okoro 2012; Mather 2013"},
    ("Salmonella",None,None):                 {"rate_low":1.5e-7,"rate_mid":3.4e-7, "rate_high":6.0e-7, "reference":"Salmonella genus (Okoro 2012)"},
    ("Enterobacter","cloacae",None):          {"rate_low":3.0e-7,"rate_mid":5.0e-7, "rate_high":1.0e-6, "reference":"Harada 2021; Peirano 2018"},
    ("Enterobacter","hormaechei",None):       {"rate_low":3.0e-7,"rate_mid":5.0e-7, "rate_high":1.0e-6, "reference":"E. cloacae complex range approx"},
    ("Enterobacter",None,None):               {"rate_low":3.0e-7,"rate_mid":5.0e-7, "rate_high":1.0e-6, "reference":"Enterobacter genus approx"},
    ("Neisseria","gonorrhoeae",None):         {"rate_low":5.0e-7,"rate_mid":1.0e-6, "rate_high":2.0e-6, "reference":"Grad 2014; Sánchez-Busó 2019"},
    ("Neisseria",None,None):                  {"rate_low":5.0e-7,"rate_mid":9.0e-7, "rate_high":2.0e-6, "reference":"Neisseria genus avg"},
    ("Streptococcus","pneumoniae",None):      {"rate_low":5.0e-7,"rate_mid":1.6e-6, "rate_high":3.0e-6, "reference":"Croucher 2013 Science; Chewapreecha 2014"},
    ("Streptococcus","pyogenes",None):        {"rate_low":5.0e-7,"rate_mid":1.2e-6, "rate_high":2.0e-6, "reference":"Nasser 2014 PNAS"},
    ("Streptococcus",None,None):              {"rate_low":5.0e-7,"rate_mid":1.4e-6, "rate_high":3.0e-6, "reference":"Streptococcus genus avg"},
    ("Haemophilus","influenzae",None):        {"rate_low":5.0e-7,"rate_mid":1.0e-6, "rate_high":2.0e-6, "reference":"De Chiara 2014; Connor 2015"},
    ("Haemophilus","ducreyi",None):           {"rate_low":3.0e-7,"rate_mid":7.0e-7, "rate_high":1.5e-6, "reference":"Approx; limited data"},
    ("Haemophilus",None,None):                {"rate_low":3.0e-7,"rate_mid":8.5e-7, "rate_high":2.0e-6, "reference":"Haemophilus genus avg"},
    ("Moraxella","catarrhalis",None):         {"rate_low":5.0e-7,"rate_mid":9.0e-7, "rate_high":1.8e-6, "reference":"Approx; limited data"},
    ("Moraxella",None,None):                  {"rate_low":5.0e-7,"rate_mid":9.0e-7, "rate_high":1.8e-6, "reference":"Moraxella genus avg"},
    ("Listeria","monocytogenes",None):       {"rate_low":2.0e-7,"rate_mid":4.5e-7, "rate_high":7.0e-7, "reference":"Moura 2016 Nat Microbiol"},
    ("Listeria",None,None):                   {"rate_low":2.0e-7,"rate_mid":4.5e-7, "rate_high":7.0e-7, "reference":"Listeria genus avg"},
    ("Shigella","sonnei",None):              {"rate_low":3.0e-7,"rate_mid":6.0e-7, "rate_high":1.0e-6, "reference":"Holt 2012 Nat Genet"},
    ("Shigella",None,None):                   {"rate_low":3.0e-7,"rate_mid":6.0e-7, "rate_high":1.0e-6, "reference":"Shigella genus (Holt 2012)"},
    ("Bacillus","anthracis",None):            {"rate_low":5.0e-8,"rate_mid":1.5e-7, "rate_high":4.0e-7, "reference":"Kenefic 2009; Pearson 2004"},
    ("Bacillus","cereus",None):              {"rate_low":1.0e-7,"rate_mid":3.0e-7, "rate_high":6.0e-7, "reference":"B. cereus group approx"},
    ("Bacillus",None,None):                   {"rate_low":5.0e-8,"rate_mid":2.0e-7, "rate_high":6.0e-7, "reference":"Bacillus genus avg"},
    ("Brucella","melitensis",None):           {"rate_low":5.0e-8,"rate_mid":2.0e-7, "rate_high":5.0e-7, "reference":"Foster 2009; Tan 2015"},
    ("Brucella",None,None):                   {"rate_low":5.0e-8,"rate_mid":2.0e-7, "rate_high":5.0e-7, "reference":"Brucella genus avg"},
    ("Pasteurella","multocida",None):        {"rate_low":3.0e-7,"rate_mid":6.0e-7, "rate_high":1.2e-6, "reference":"Pasteurellaceae approx"},
    ("Pasteurella",None,None):                {"rate_low":3.0e-7,"rate_mid":6.0e-7, "rate_high":1.2e-6, "reference":"Pasteurella genus avg"},
    ("Proteus","vulgaris",None):             {"rate_low":2.0e-7,"rate_mid":4.0e-7, "rate_high":8.0e-7, "reference":"Enterobacterales range approx"},
    ("Proteus",None,None):                    {"rate_low":2.0e-7,"rate_mid":4.0e-7, "rate_high":8.0e-7, "reference":"Proteus genus avg"},
    ("Stenotrophomonas","maltophilia",None):  {"rate_low":5.0e-7,"rate_mid":1.0e-6, "rate_high":2.0e-6, "reference":"Esposito 2017 approx"},
    ("Stenotrophomonas",None,None):           {"rate_low":5.0e-7,"rate_mid":1.0e-6, "rate_high":2.0e-6, "reference":"Stenotrophomonas genus avg"},
    ("Clostridioides","difficile",None):     {"rate_low":8.0e-7,"rate_mid":1.4e-6, "rate_high":2.5e-6, "reference":"He 2013 Nat Genet"},
    ("Clostridium","difficile",None):        {"rate_low":8.0e-7,"rate_mid":1.4e-6, "rate_high":2.5e-6, "reference":"He 2013 Nat Genet"},
    ("Clostridioides",None,None): {"rate_low":8.0e-7,"rate_mid":1.4e-6,"rate_high":2.5e-6,"reference":"C. difficile (He 2013)"},
    ("Clostridium",None,None):    {"rate_low":8.0e-7,"rate_mid":1.4e-6,"rate_high":2.5e-6,"reference":"C. difficile (He 2013)"},
    ("Campylobacter","jejuni",None):         {"rate_low":2.0e-6,"rate_mid":3.2e-6, "rate_high":5.0e-6, "reference":"Wilson 2009; Mourkas 2020"},
    ("Campylobacter",None,None):              {"rate_low":2.0e-6,"rate_mid":3.2e-6, "rate_high":5.0e-6, "reference":"Campylobacter genus avg"},
    ("Helicobacter","pylori",None):          {"rate_low":1.0e-5,"rate_mid":2.0e-5, "rate_high":6.0e-5, "reference":"Kennemann 2011 PNAS; Didelot 2013"},
    ("Helicobacter",None,None):               {"rate_low":1.0e-5,"rate_mid":2.0e-5, "rate_high":6.0e-5, "reference":"Helicobacter genus avg"},
    ("Vibrio","cholerae",None):              {"rate_low":5.0e-7,"rate_mid":8.0e-7, "rate_high":1.5e-6, "reference":"Mutreja 2011 Nature"},
    ("Vibrio",None,None):                     {"rate_low":5.0e-7,"rate_mid":8.0e-7, "rate_high":1.5e-6, "reference":"Vibrio genus avg"},
    ("Legionella","pneumophila",None):       {"rate_low":2.0e-7,"rate_mid":5.0e-7, "rate_high":9.0e-7, "reference":"David 2016; McAdam 2014"},
    ("Legionella",None,None):                 {"rate_low":2.0e-7,"rate_mid":5.0e-7, "rate_high":9.0e-7, "reference":"Legionella genus avg"},
    ("Serratia","marcescens",None):          {"rate_low":3.0e-7,"rate_mid":5.5e-7, "rate_high":1.0e-6, "reference":"Moradigaravand 2016"},
    ("Serratia",None,None):                   {"rate_low":3.0e-7,"rate_mid":5.5e-7, "rate_high":1.0e-6, "reference":"Serratia genus avg"},
    ("Bordetella","pertussis",None):          {"rate_low":2.0e-7,"rate_mid":5.0e-7, "rate_high":8.0e-7, "reference":"Bart 2014 PNAS"},
    ("Bordetella",None,None):                 {"rate_low":2.0e-7,"rate_mid":5.0e-7, "rate_high":8.0e-7, "reference":"Bordetella genus avg"},
    ("Citrobacter",None,None):                {"rate_low":3.0e-7,"rate_mid":5.0e-7, "rate_high":9.0e-7, "reference":"Citrobacter genus (Enterobacterales range)"},
    ("Morganella",None,None):                 {"rate_low":3.0e-7,"rate_mid":5.0e-7, "rate_high":8.0e-7, "reference":"Morganella genus avg"},
    ("Burkholderia",None,None):               {"rate_low":1.0e-7,"rate_mid":7.0e-7, "rate_high":2.0e-6, "reference":"Lieberman 2011; Pearson 2020"},
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
