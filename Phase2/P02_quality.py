"""
P02: Quality control.
- Peak S/N per spectrum
- Tiers: A (S/N>10), B (5<S/N<10), C (S/N<5)
- z > 5 filter
- Output: quality_metrics_v2.csv, tier_A_z5_v2.csv        """

import numpy as np
import pandas as pd
from astropy.io import fits
from pathlib import Path
import re
import sys
sys.path.insert(0, str(Path(__file__).parent))
from config import *

try:
    from scipy.ndimage import median_filter
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


def peak_snr(w, fl, er):
    if len(fl) < 10:
        return 0.0
    if HAS_SCIPY:
        fl_s = fl.astype(np.float64)
    else:
        fl_s = fl
    snr = (fl_s - np.median(fl)) / er
    return float(np.max(snr))


def main():
    print("=" * 70)
    print("P02: QUALITY CONTROL")
    print("=" * 70)

    p01 = pd.read_csv(ANALYSIS / 'p01_clean_spectra.csv')
    print(f"Spectra from P01: {len(p01)}")

    # Load z from catalog
    cat = pd.read_csv(ANALYSIS / 'master_features_v4.csv')
    cat['NIRSpec_ID_str'] = cat['NIRSpec_ID'].astype(str).str.zfill(8)
    z_map = dict(zip(cat['NIRSpec_ID_str'], cat['z_Spec']))
    print(f"z_map: {len(z_map)} entries")

    results = []
    for _, row in p01.iterrows():
        f = SPECTRA / row['file']
        m = re.search(r'(\d{8})', row['file'])
        if not m:
            continue
        sid = m.group(1)
        z = z_map.get(sid)

        try:
            with fits.open(f) as hdul:
                d = hdul[1].data
                cn = hdul[1].columns.names
                wc = next((c for c in cn if 'wave' in c.lower()), cn[0])
                fc = next((c for c in cn if 'flux' in c.lower() and 'err' not in c.lower()), cn[1])
                ec = next((c for c in cn if 'err' in c.lower()), cn[2])
                w = d[wc].copy(); fl = d[fc].copy(); er = d[ec].copy()
            if np.nanmax(w) < 100:
                w = w * 1e4
            good = np.isfinite(fl) & np.isfinite(er) & (er > 0)
            if good.sum() < 50:
                continue
            w = w[good]; fl = fl[good]; er = er[good]
        except Exception:
            continue

        snr = peak_snr(w, fl, er)

        if snr > SNR_MIN_TIER_A:
            tier = 'A'
        elif snr > SNR_MIN_TIER_B:
            tier = 'B'
        else:
            tier = 'C'

        results.append({
            'file': row['file'],
            'ID': sid,
            'z': z,
            'snr_peak': snr,
            'tier': tier,
            'grating': row['grating'],
            'field': 'goods-s' if 'goods-s' in row['file'] else 'goods-n',
        })

    df = pd.DataFrame(results)
    print(f"\nTotal processed: {len(df)}")
    print(f"Tier A (S/N>{SNR_MIN_TIER_A}): {(df['tier']=='A').sum()}")
    print(f"Tier B (S/N>{SNR_MIN_TIER_B}): {(df['tier']=='B').sum()}")
    print(f"Tier C: {(df['tier']=='C').sum()}")
    print()

    # Objects with valid z
    df_z = df.dropna(subset=['z'])
    print(f"With valid z: {len(df_z)}")
    print(f"With z > 5:   {(df_z['z'] > Z_MIN).sum()}")
    print(f"With z < 5:   {(df_z['z'] < Z_MIN).sum()}")
    print()

    # Tier A with z ∈ [5, 8]
    tier_a_z5 = df[(df['tier'] == 'A') & (df['z'] > Z_MIN) & (df['z'] < Z_MAX)].copy()
    print(f"Tier A with z ∈ [{Z_MIN}, {Z_MAX}]: {len(tier_a_z5)}")
    print()

    # Deduplicate by ID (keep the highest S/N)
    tier_a_unique = tier_a_z5.sort_values('snr_peak', ascending=False).drop_duplicates('ID')
    print(f"Unique Tier A z>5: {len(tier_a_unique)}")
    print()
    print("Field distribution:")
    print(tier_a_unique['field'].value_counts())
    print()

    df.to_csv(OUT_QUALITY, index=False)
    tier_a_unique.to_csv(ANALYSIS / 'tier_A_z5_v2.csv', index=False)
    print(f"Saved: {OUT_QUALITY}")
    print(f"Saved: tier_A_z5_v2.csv")


if __name__ == '__main__':
    main()
