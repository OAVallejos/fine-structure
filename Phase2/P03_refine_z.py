"""
P03: z_Spec bias correction.
- Measures Δz(Hβ) and Δz([OIII]) in the sample
- Does NOT refine object by object
- Corrects the systematic bias
- Output: refined_z.csv (with corrected z)
"""

import numpy as np
import pandas as pd
from astropy.io import fits
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from config import *


def measure_dz(f, z_spec):
    """Measures Δz using the Hβ and [OIII] 5007 peaks."""
    try:
        with fits.open(f) as hdul:
            d = hdul[1].data
            cn = hdul[1].columns.names
            wc = next((c for c in cn if 'wave' in c.lower()), cn[0])
            fc = next((c for c in cn if 'flux' in c.lower() and 'err' not in c.lower()), cn[1])
            w = d[wc].copy(); fl = d[fc].copy()
        if np.nanmax(w) < 100:
            w = w * 1e4
        good = np.isfinite(fl)
        if good.sum() < 50:
            return None
        w = w[good]; fl = fl[good]
    except Exception:
        return None

    dz = {}
    for name, lam_lab in [('Hb', LAMBDA_VAC['Hβ']), ('OIII', LAMBDA_VAC['[OIII]5007'])]:
        lam_obs = lam_lab * (1 + z_spec)
        mask = (w > lam_obs - 50) & (w < lam_obs + 50)
        if mask.sum() < 5:
            continue
        w_win = w[mask]; fl_win = fl[mask]
        # Subtract local continuum
        cont = np.percentile(fl_win, 10)
        fl_sub = fl_win - cont
        fl_sub = np.clip(fl_sub, 0, None)
        if fl_sub.sum() <= 0:
            continue
        cen = np.sum(w_win * fl_sub) / np.sum(fl_sub)
        z_peak = cen / lam_lab - 1
        dz[name] = z_peak - z_spec

    if 'Hb' in dz and 'OIII' in dz:
        return {'dz_hb': dz['Hb'], 'dz_oiii': dz['OIII'],
                'dz_avg': (dz['Hb'] + dz['OIII']) / 2}
    return None


def main():
    print("=" * 70)
    print("P03: Z_SPEC BIAS CORRECTION")
    print("=" * 70)

    tier_a = pd.read_csv(ANALYSIS / 'tier_A_z5_v2.csv')
    print(f"Tier A z>5: {len(tier_a)}")

    # Measure Δz in the sample
    rows = []
    for _, row in tier_a.iterrows():
        f = SPECTRA / row['file']
        if not f.exists():
            continue
        r = measure_dz(f, row['z'])
        if r is None:
            continue
        rows.append({
            'ID': str(row['ID']).zfill(8),
            'file': row['file'],
            'z_Spec': row['z'],
            'dz_hb': r['dz_hb'],
            'dz_oiii': r['dz_oiii'],
            'dz_avg': r['dz_avg'],
        })

    df = pd.DataFrame(rows)
    print(f"\nValid measurements: {len(df)}/{len(tier_a)}")

    if len(df) > 0:
        dz_mean = df['dz_avg'].mean()
        dz_std = df['dz_avg'].std()
        dz_sem = dz_std / np.sqrt(len(df))

        print()
        print("=== Z_SPEC BIAS ===")
        print(f"  Δz(Hβ):      {df['dz_hb'].mean():+.6f} ± {df['dz_hb'].std():.6f}")
        print(f"  Δz([OIII]):  {df['dz_oiii'].mean():+.6f} ± {df['dz_oiii'].std():.6f}")
        print(f"  Mean Δz:     {dz_mean:+.6f} ± {dz_std:.6f}")
        print(f"  SEM:         {dz_sem:.6f}")
        print(f"  Significance: {dz_mean/dz_sem:+.2f}σ")
        print()
        print(f"=== Hβ vs [OIII] COMPARISON ===")
        diff = df['dz_hb'] - df['dz_oiii']
        print(f"  Δz(Hβ) − Δz([OIII]): {diff.mean():+.6f} ± {diff.std():.6f}")
        print(f"  t = {diff.mean() / (diff.std() / np.sqrt(len(df))):+.2f}σ")
        print()

        # Correct z_Spec by subtracting the bias
        df['z_Corrected'] = df['z_Spec'] + dz_mean
        print(f"=== CORRECTION ===")
        print(f"  z_Corrected = z_Spec − ({dz_mean:+.6f})")
        print(f"  Example: z_Spec=5.500000 → z_Corrected={5.5 - dz_mean:.6f}")

    df.to_csv(OUT_REFINED_Z, index=False)
    print(f"\nSaved: {OUT_REFINED_Z}")


if __name__ == '__main__':
    main()
