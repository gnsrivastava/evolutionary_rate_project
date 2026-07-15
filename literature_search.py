#!/usr/bin/env python3
import requests
import time

def search_evolutionary_rates(genus, species, limit=5):
    """
    Queries Europe PMC for abstracts matching the species and evolutionary rate keywords.
    """
    if not genus:
        return []
    
    species_query = f'"{genus} {species}"' if species else f'"{genus}"'
    # Target molecular clock, substitution rate, or SNPs per year literature
    query = f'{species_query} AND ("molecular clock" OR "substitution rate" OR "mutation rate" OR "evolutionary rate" OR "SNPs per year")'
    
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    params = {
        "query": query,
        "format": "json",
        "pageSize": limit,
        "resultType": "core"
    }
    
    try:
        response = requests.get(url, params=params, timeout=20)
        response.raise_for_status()
        data = response.json()
        results = data.get("resultList", {}).get("result", [])
        
        articles = []
        for r in results:
            articles.append({
                "id": r.get("id"),
                "title": r.get("title", ""),
                "abstract": r.get("abstractText", ""),
                "authorString": r.get("authorString", ""),
                "journal": r.get("journalInfo", {}).get("journal", {}).get("title", ""),
                "year": r.get("pubYear", "")
            })
        return articles
    except Exception as e:
        print(f"Literature search failed for {genus} {species or ''}: {e}")
        return []
