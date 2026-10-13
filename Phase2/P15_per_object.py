"""
P15: Per-object offset.
"""
import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.optimize import curve_fit
from pathlib import Path
import sys, gc
sys.path.insert(0, str(Path(__file__).parent))
from config import *


def measure_offset(f, z):
    """Measures Hβ vs [OIII] 5007 offset in a single object."""
    try:
        with fits.open(SPECTRA / f, memmap=False) as hdul:
            d = hdul[1].data
            cn = hdul[1].columns.names
            wc = next((c for c in cn if 'wave' in c.lower()), cn[0])
            fc = next((c for c in cn if 'flux' in c.lower() and 'err' not in c.lower()), cn[1])
            ec = next((c for c in cn if 'err' in c.lower()), cn[2])
            w = d[wc].copy(); fl = d[fc].copy(); er = d[ec].copy()
        if np.nanmax(w) < 100: w = w * 1e4
        w_rest = w / (1 + z)
        good = np.isfinite(fl) & np.isfinite(er) & (er > 0)
        if good.sum() < 50: return None
        w_rest = w_rest[good]; fl = fl[good]; er = er[good]
    except Exception:
        return None

    # Measure Hβ
    lam_hb = LAMBDA_VAC['Hβ']
    mask = (w_rest > lam_hb - 30) & (w_rest < lam_hb + 30)
    if mask.sum() < 5: return None
    w_hb = w_rest[mask]; fl_hb = fl[mask]
    cont = np.percentile(fl_hb, 10)
    fl_hb_sub = np.clip(fl_hb - cont, 0, None)
    if fl_hb_sub.sum() <= 0: return None
    mu_hb = np.sum(w_hb * fl_hb_sub) / np.sum(fl_hb_sub)

    # Measure [OIII] 5007
    lam_oiii = LAMBDA_VAC['[OIII]5007']
    mask = (w_rest > lam_oiii - 30) & (w_rest < lam_oiii + 30)
    if mask.sum() < 5: return None
    w_oiii = w_rest[mask]; fl_oiii = fl[mask]
    cont = np.percentile(fl_oiii, 10)
    fl_oiii_sub = np.clip(fl_oiii - cont, 0, None)
    if fl_oiii_sub.sum() <= 0: return None
    mu_oiii = np.sum(w_oiii * fl_oiii_sub) / np.sum(fl_oiii_sub)

    d_hb = mu_hb - lam_hb
    d_oiii = mu_oiii - lam_oiii
    diff = d_hb - d_oiii
    dalpha = -diff / (2 * lam_hb)
    return d_hb, d_oiii, diff, dalpha


def main():
    print("=" * 70)
    print("P15: PER-OBJECT OFFSET")
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
        r = measure_offset(row['file'], z)
        if r is None: continue
        d_hb, d_oiii, diff, dalpha = r
        results.append({
            'ID': sid, 'file': row['file'], 'z': z, 'field': row['field'],
            'grating': row['grating'],
            'd_hb': d_hb, 'd_oiii': d_oiii, 'diff': diff, 'dalpha': dalpha,
        })
        gc.collect()

    df = pd.DataFrame(results)
    print(f"\nObjects with measurement: {len(df)}")
    print()
    print("Per-object Δα/α distribution:")
    print(df['dalpha'].describe())
    print()
    print("Distribution by field:")
    for field in ['goods-s', 'goods-n']:
        sub = df[df['field'] == field]
        print(f"  {field}: N={len(sub)}, mean={sub['dalpha'].mean()*1e6:+.2f} ppm, "
              f"std={sub['dalpha'].std()*1e6:.2f} ppm")
    print()
    print("Distribution by grating:")
    for grat in ['G395M', 'G395H']:
        sub = df[df['grating'] == grat]
        if len(sub) > 0:
            print(f"  {grat}: N={len(sub)}, mean={sub['dalpha'].mean()*1e6:+.2f} ppm, "
                  f"std={sub['dalpha'].std()*1e6:.2f} ppm")

    df.to_csv(ANALYSIS / 'per_object_offset.csv', index=False)
    print(f"\nSaved: per_object_offset.csv")


if __name__ == '__main__':
    main()
