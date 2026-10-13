"""
P03b: Refine z by cross-correlation.
- Template: Hβ + [OIII] 4959 + [OIII] 5007 (Gaussians)
- Finds the z that maximizes the correlation
"""
import numpy as np
import pandas as pd
from astropy.io import fits
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from config import *


def make_template(wave, z, sigma=3.0):
    """Template of Hβ + [OIII] in the observed frame."""
    tmpl = np.zeros_like(wave)
    for lam, amp in [(LAMBDA_VAC['Hβ'], 1.0),
                     (LAMBDA_VAC['[OIII]4959'], 0.5),
                     (LAMBDA_VAC['[OIII]5007'], 1.5)]:
        lam_obs = lam * (1 + z)
        tmpl += amp * np.exp(-0.5 * ((wave - lam_obs) / sigma)**2)
    return tmpl


def refine_z_crosscorr(w, fl, z_spec, z_range=0.01, n_steps=200):
    """Finds the z that maximizes the correlation with the template."""
    # Subtract continuum
    cont = np.median(fl)
    fl_sub = fl - cont

    # z grid
    z_grid = np.linspace(z_spec - z_range, z_spec + z_range, n_steps)
    corr = np.zeros_like(z_grid)

    for i, z in enumerate(z_grid):
        tmpl = make_template(w, z)
        # Normalized correlation
        num = np.sum(fl_sub * tmpl)
        den = np.sqrt(np.sum(fl_sub**2) * np.sum(tmpl**2))
        corr[i] = num / den if den > 0 else 0

    # Maximum
    idx = np.argmax(corr)
    z_best = z_grid[idx]
    return z_best, corr[idx]


def main():
    print("=" * 70)
    print("P03b: REFINE Z BY CROSS-CORRELATION")
    print("=" * 70)

    tier_a = pd.read_csv(ANALYSIS / 'tier_A_z5_v2.csv')
    tier_a['ID'] = tier_a['ID'].astype(str).str.zfill(8)
    print(f"Tier A: {len(tier_a)}")

    results = []
    for _, row in tier_a.iterrows():
        f = SPECTRA / row['file']
        if not f.exists():
            continue
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
                continue
            w = w[good]; fl = fl[good]
        except Exception:
            continue

        try:
            z_best, corr_max = refine_z_crosscorr(w, fl, row['z'])
            results.append({
                'ID': row['ID'],
                'file': row['file'],
                'z_Spec': row['z'],
                'z_Xcorr': z_best,
                'dz': z_best - row['z'],
                'corr_max': corr_max,
            })
        except Exception:
            continue

    df = pd.DataFrame(results)
    print(f"\nRefined: {len(df)}/{len(tier_a)}")
    if len(df) > 0:
        print(f"Mean Δz:  {df['dz'].mean():+.6f} ± {df['dz'].std():.6f}")
        print(f"Median Δz: {df['dz'].median():+.6f}")
        print(f"corr_max:   {df['corr_max'].mean():.4f} ± {df['corr_max'].std():.4f}")
        print()
        print("First 10:")
        print(df[['ID', 'z_Spec', 'z_Xcorr', 'dz', 'corr_max']].head(10).to_string())

    df.to_csv(ANALYSIS / 'refined_z_xcorr.csv', index=False)
    print(f"\nSaved: refined_z_xcorr.csv")


if __name__ == '__main__':
    main()
