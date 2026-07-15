#!/usr/bin/env python3
import os
import json

CACHE_DIR = "cache"
CACHE_FILE = os.path.join(CACHE_DIR, "rate_cache.json")

def load_cache():
    if not os.path.exists(CACHE_DIR):
        os.makedirs(CACHE_DIR, exist_ok=True)
    if not os.path.exists(CACHE_FILE):
        return {}
    try:
        with open(CACHE_FILE, 'r') as f:
            data = json.load(f)
            # Reconvert stringified dictionary keys back to internal operational formats if necessary
            return data
    except Exception:
        return {}

def save_cache(cache_data):
    if not os.path.exists(CACHE_DIR):
        os.makedirs(CACHE_DIR, exist_ok=True)
    try:
        with open(CACHE_FILE, 'w') as f:
            json.dump(cache_data, f, indent=4)
    except Exception as e:
        print(f"Failed to save local cache file: {e}")
