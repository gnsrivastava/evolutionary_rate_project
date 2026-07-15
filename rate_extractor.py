#!/usr/bin/env python3
import re

# Regex matching patterns like 1.2 x 10^-6, 1.2e-6, 3.4 × 10−7
RATE_PATTERN = re.compile(
    r'\b(\d+(?:\.\d+)?)\s*(?:[x××]\\s*10\s*[-\u2212\u2013]\s*(\d+)|\s*e\s*[-\u2212\u2013]\s*(\d+))\b', 
    re.IGNORECASE
)

def extract_rates_from_text(text):
    """
    Parses text to extract numerical values that look like substitution rates.
    Returns a list of float numbers.
    """
    if not text:
        return []
    
    matches = RATE_PATTERN.findall(text)
    extracted = []
    
    for base, exp1, exp2 in matches:
        exponent = exp1 if exp1 else exp2
        try:
            val = float(base) * (10 ** -int(exponent))
            # Evolutionary rates for bacteria typically fall between 1e-9 and 1e-4
            if 1e-9 <= val <= 1e-4:
                extracted.append(val)
        except ValueError:
            continue
            
    return extracted

def synthesize_rate_record(genus, species, articles):
    """
    Analyzes literature search matches, extracts numbers, and computes low/mid/high boundaries.
    """
    all_rates = []
    best_ref = "Dynamic Literature Search"
    
    for art in articles:
        text_to_scan = f"{art['title']} {art['abstract']}"
        rates = extract_rates_from_text(text_to_scan)
        if rates:
            all_rates.extend(rates)
            if best_ref == "Dynamic Literature Search":
                auth = art['authorString'].split(',')[0] if art['authorString'] else "Unknown"
                best_ref = f"{auth} {art['year']} {art['journal']}".strip()

    if not all_rates:
        return None

    # De-duplicate and sort
    all_rates = sorted(list(set(all_rates)))
    
    if len(all_rates) == 1:
        mid = all_rates[0]
        low = mid * 0.8
        high = mid * 1.2
    elif len(all_rates) == 2:
        low = all_rates[0]
        high = all_rates[1]
        mid = (low + high) / 2
    else:
        low = all_rates[0]
        mid = all_rates[len(all_rates)//2]
        high = all_rates[-1]
        
    return {
        "rate_low": low,
        "rate_mid": mid,
        "rate_high": high,
        "reference": f"{best_ref} (Extracted dynamically)"
    }
