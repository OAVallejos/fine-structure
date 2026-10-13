"""
P04: Rest-frame stack with bias-corrected z.
- Wide grid (4800–6700 Å)
- Hα filter within the range
- Output: stack_v2.fits
"""

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from config import *


def main():
    print("=" * 70)
    print("P04: BUILD STACK")
    print("=" * 70)

    grid = np.linspace(GRID_MIN, GRID_MAX, GRID_N)
    print(f"Grid: {GRID_MIN}–{GRID_MAX} Å, {GRID_N} points")

    # Load corrected z
    rz = pd.read_csv(OUT_REFINED_Z)
    rz['ID'] = rz['ID'].astype(str).str.zfill(8)
    z_map = dict(zip(rz['ID'], rz['z_Corrected']))
    print(f"Objects with corrected z: {len(z_map)}")

    # Load tier A
    tier_a = pd.read_csv(ANALYSIS / 'tier_A_z5_v2.csv')
    tier_a['ID'] = tier_a['ID'].astype(str).str.zfill(8)

    num = np.zeros_like(grid)
    den = np.zeros_like(grid)
    n_ok = 0
    n_skip = 0
    n_ha_out = 0

    for _, row in tier_a.iterrows():
        sid = row['ID']
        if sid not in z_map:
            n_skip += 1
            continue
        z = z_map[sid]
        if z <= 0:
            n_skip += 1
            continue

        # Check Hα is within range
        lam_ha_obs = LAMBDA_VAC['Hα'] * (1 + z)
        if not (X1D_WAVE_MIN < lam_ha_obs < X1D_WAVE_MAX):
            n_ha_out += 1
            continue

        try:
            with fits.open(SPECTRA / row['file']) as hdul:
                d = hdul[1].data
                cn = hdul[1].columns.names
                wc = next((c for c in cn if 'wave' in c.lower()), cn[0])
                fc = next((c for c in cn if 'flux' in c.lower() and 'err' not in c.lower()), cn[1])
                ec = next((c for c in cn if 'err' in c.lower()), cn[2])
                w = d[wc].copy(); fl = d[fc].copy(); er = d[ec].copy()
            if np.nanmax(w) < 100:
                w = w * 1e4
            w_rest = w / (1 + z)
            good = np.isfinite(fl) & np.isfinite(er) & (er > 0)
            if good.sum() < 50:
                n_skip += 1
                continue
            f_int = np.interp(grid, w_rest[good], fl[good], left=0, right=0)
            e_int = np.interp(grid, w_rest[good], er[good], left=np.inf, right=np.inf)
            wgt = 1.0 / e_int**2
            wgt[~np.isfinite(wgt)] = 0
            num += f_int * wgt
            den += wgt
            n_ok += 1
        except Exception:
            n_skip += 1
            continue

    sf = np.where(den > 0, num / den, 0)
    se = np.where(den > 0, 1.0 / np.sqrt(den), np.inf)

    t = Table()
    t['WAVE_REST'] = grid
    t['FLUX'] = sf
    t['FLUX_ERR'] = se
    t.write(OUT_STACK, overwrite=True)

    print(f"\nStack: {n_ok} objects")
    print(f"Skip: {n_skip}")
    print(f"Hα out of range: {n_ha_out}")
    print(f"Saved: {OUT_STACK}")

    # Check lines
    print("\nDetected lines:")
    for name, lam in [('Hβ', LAMBDA_VAC['Hβ']),
                      ('[OIII]4959', LAMBDA_VAC['[OIII]4959']),
                      ('[OIII]5007', LAMBDA_VAC['[OIII]5007']),
                      ('Hα', LAMBDA_VAC['Hα']),
                      ('[NII]6584', LAMBDA_VAC['[NII]6584'])]:
        mask = (grid > lam - 8) & (grid < lam + 8)
        if mask.sum() > 0:
            peak = np.max(sf[mask])
            cont = np.median(sf[(grid > 5200) & (grid < 5400)])
            print(f"  {name:12s}  λ={lam:.2f}  peak={peak:.3e}  cont={cont:.3e}")


if __name__ == '__main__':
    main()
