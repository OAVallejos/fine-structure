# P13_multiline_z.py
import numpy as np
import pandas as pd
from astropy.io import fits
from pathlib import Path
import sys, gc
sys.path.insert(0, str(Path(__file__).parent))
from config import *

# Lines to use
LINES = {
    'Hβ': LAMBDA_VAC['Hβ'],
    'Hα': LAMBDA_VAC['Hα'],
    '[OIII]5007': LAMBDA_VAC['[OIII]5007'],
    '[NII]6584': LAMBDA_VAC['[NII]6584'],
}

def measure_z_multiline(f, z_spec):
    """Measures z using multiple lines."""
    try:
        with fits.open(f, memmap=False) as hdul:
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

    z_per_line = {}
    for name, lam_teo in LINES.items():
        lam_obs = lam_teo * (1 + z_spec)
        # Narrow window (±30 Å observed)
        mask = (w > lam_obs - 50) & (w < lam_obs + 50)
        if mask.sum() < 5:
            continue
        w_win = w[mask]; fl_win = fl[mask]
        # Subtract continuum
        cont = np.percentile(fl_win, 10)
        fl_sub = fl_win - cont
        fl_sub = np.clip(fl_sub, 0, None)
        if fl_sub.sum() <= 0:
            continue
        # Centroid
        cen = np.sum(w_win * fl_sub) / np.sum(fl_sub)
        z_line = cen / lam_teo - 1
        z_per_line[name] = z_line

    return z_per_line


def main():
    print("=" * 70)
    print("P13: MULTI-LINE Z")
    print("=" * 70)

    tier_a = pd.read_csv(ANALYSIS / 'tier_A_z5_v2.csv')
    tier_a['ID'] = tier_a['ID'].astype(str).str.zfill(8)
    print(f"Tier A: {len(tier_a)}")

    results = []
    for _, row in tier_a.iterrows():
        f = SPECTRA / row['file']
        if not f.exists():
            continue
        z_per_line = measure_z_multiline(f, row['z'])
        if z_per_line is None or len(z_per_line) < 2:
            continue
        # Weighted average
        z_vals = list(z_per_line.values())
        z_avg = np.mean(z_vals)
        z_std = np.std(z_vals) if len(z_vals) > 1 else 0
        results.append({
            'ID': row['ID'],
            'z_Spec': row['z'],
            'z_Hβ': z_per_line.get('Hβ', np.nan),
            'z_Hα': z_per_line.get('Hα', np.nan),
            'z_OIII5007': z_per_line.get('[OIII]5007', np.nan),
            'z_NII6584': z_per_line.get('[NII]6584', np.nan),
            'z_multiline': z_avg,
            'z_std': z_std,
            'n_lines': len(z_per_line),
        })
        gc.collect()

    df = pd.DataFrame(results)
    print(f"\nObjects with multi-line z: {len(df)}")
    if len(df) > 0:
        print(f"Δz(Hβ)        = {(df['z_Hβ'] - df['z_Spec']).mean():+.6f} ± {(df['z_Hβ'] - df['z_Spec']).std():.6f}")
        print(f"Δz(Hα)        = {(df['z_Hα'] - df['z_Spec']).mean():+.6f} ± {(df['z_Hα'] - df['z_Spec']).std():.6f}")
        print(f"Δz([OIII]5007)= {(df['z_OIII5007'] - df['z_Spec']).mean():+.6f} ± {(df['z_OIII5007'] - df['z_Spec']).std():.6f}")
        print(f"Δz([NII]6584) = {(df['z_NII6584'] - df['z_Spec']).mean():+.6f} ± {(df['z_NII6584'] - df['z_Spec']).std():.6f}")
        print(f"Δz(multiline) = {(df['z_multiline'] - df['z_Spec']).mean():+.6f} ± {(df['z_multiline'] - df['z_Spec']).std():.6f}")

    df.to_csv(ANALYSIS / 'z_multiline.csv', index=False)
    print(f"\nSaved: z_multiline.csv")


if __name__ == '__main__':
    main()
