# P17_filter_absorption.py
"""
P17: Filter objects with absorption in Hβ or [OIII] 5007.
"""
import numpy as np
import pandas as pd
from astropy.io import fits
from pathlib import Path
import sys, gc
sys.path.insert(0, str(Path(__file__).parent))
from config import *


def check_absorption(f, z):
    """Checks whether Hβ or [OIII] 5007 are in absorption."""
    try:
        with fits.open(SPECTRA / f, memmap=False) as hdul:
            d = hdul[1].data
            cn = hdul[1].columns.names
            wc = next((c for c in cn if 'wave' in c.lower()), cn[0])
            fc = next((c for c in cn if 'flux' in c.lower() and 'err' not in c.lower()), cn[1])
            w = d[wc].copy(); fl = d[fc].copy()
        if np.nanmax(w) < 100: w = w * 1e4
        w_rest = w / (1 + z)
        good = np.isfinite(fl)
        if good.sum() < 50: return None
        w_rest = w_rest[good]; fl = fl[good]
    except Exception:
        return None

    # Window around each line
    lam_hb = LAMBDA_VAC['Hβ']
    lam_oiii = LAMBDA_VAC['[OIII]5007']

    # Hβ
    mask = (w_rest > lam_hb - 30) & (w_rest < lam_hb + 30)
    if mask.sum() < 5: return None
    fl_hb = fl[mask]
    cont_hb = np.percentile(fl_hb, 50)
    # Maximum and minimum peak
    peak_hb = np.max(fl_hb - cont_hb)
    min_hb = np.min(fl_hb - cont_hb)

    # [OIII] 5007
    mask = (w_rest > lam_oiii - 30) & (w_rest < lam_oiii + 30)
    if mask.sum() < 5: return None
    fl_oiii = fl[mask]
    cont_oiii = np.percentile(fl_oiii, 50)
    peak_oiii = np.max(fl_oiii - cont_oiii)
    min_oiii = np.min(fl_oiii - cont_oiii)

    return {
        'peak_hb': peak_hb,
        'min_hb': min_hb,
        'peak_oiii': peak_oiii,
        'min_oiii': min_oiii,
        'absorption_hb': min_hb < -3 * np.std(fl_hb),
        'absorption_oiii': min_oiii < -3 * np.std(fl_oiii),
    }


def main():
    print("=" * 70)
    print("P17: FILTER ABSORPTION")
    print("=" * 70)

    mz = pd.read_csv(ANALYSIS / 'z_multiline.csv')
    mz['ID'] = mz['ID'].astype(str).str.zfill(8)
    z_map = dict(zip(mz['ID'], mz['z_multiline']))

    tier_a = pd.read_csv(ANALYSIS / 'tier_A_z5_v2.csv')
    tier_a['ID'] = tier_a['ID'].astype(str).str.zfill(8)

    results = []
    for _, row in tier_a.iterrows():
        sid = row['ID']
        if sid not in z_map: continue
        z = z_map[sid]
        if z <= 0: continue
        r = check_absorption(row['file'], z)
        if r is None: continue
        results.append({
            'ID': sid, 'file': row['file'], 'z': z,
            'field': row['field'], 'grating': row['grating'],
            **r,
        })
        gc.collect()

    df = pd.DataFrame(results)
    print(f"\nObjects analyzed: {len(df)}")
    print(f"With absorption in Hβ: {df['absorption_hb'].sum()}")
    print(f"With absorption in [OIII] 5007: {df['absorption_oiii'].sum()}")
    print(f"With absorption in both: {(df['absorption_hb'] & df['absorption_oiii']).sum()}")
    print()

    # Save list of clean objects
    clean = df[~df['absorption_hb'] & ~df['absorption_oiii']]
    print(f"CLEAN objects (no absorption): {len(clean)}")
    clean.to_csv(ANALYSIS / 'tier_A_z5_clean.csv', index=False)
    print(f"Saved: tier_A_z5_clean.csv")

    # Save all
    df.to_csv(ANALYSIS / 'absorption_check.csv', index=False)


if __name__ == '__main__':
    main()
