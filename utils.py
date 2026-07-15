#!/usr/bin/env python3
import pandas as pd

def save_output(df, items, fn):
    print("="*70); print(f"STEP 5: Saving → {fn}"); print("="*70)
    if df.empty: return
    df.to_csv(fn, index=False)
    print(f"  ✓ Summary saved to → {fn}")
