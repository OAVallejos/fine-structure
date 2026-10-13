# P16_hb_profiles.py
"""
P16: Visualize Hβ profiles in individual objects.
"""
import numpy as np
import pandas as pd
from astropy.io import fits
from pathlib import Path
import sys, gc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).parent))
from config import *


def main():
    mz = pd.read_csv(ANALYSIS / 'z_multiline.csv')
    mz['ID'] = mz['ID'].astype(str).str.zfill(8)
    z_map = dict(zip(mz['ID'], mz['z_multiline']))

    tier_a = pd.read_csv(ANALYSIS / 'tier_A_z5_v2.csv')
    tier_a['ID'] = tier_a['ID'].astype(str).str.zfill(8)

    # Select 12 random objects
    np.random.seed(42)
    sample = tier_a.sample(12)

    fig, axes = plt.subplots(4, 3, figsize=(15, 12))
    axes = axes.flatten()

    for i, (_, row) in enumerate(sample.iterrows()):
        sid = row['ID']
        if sid not in z_map: continue
        z = z_map[sid]
        if z <= 0: continue

        try:
            with fits.open(SPECTRA / row['file'], memmap=False) as hdul:
                d = hdul[1].data
                cn = hdul[1].columns.names
                wc = next((c for c in cn if 'wave' in c.lower()), cn[0])
                fc = next((c for c in cn if 'flux' in c.lower() and 'err' not in c.lower()), cn[1])
                w = d[wc].copy(); fl = d[fc].copy()
            if np.nanmax(w) < 100: w = w * 1e4
            w_rest = w / (1 + z)
            good = np.isfinite(fl)
            if good.sum() < 50: continue
            w_rest = w_rest[good]; fl = fl[good]
        except Exception:
            continue

        # Window around Hβ and [OIII]
        mask = (w_rest > 4840) & (w_rest > 4840) & (w_rest < 5030)
        if mask.sum() < 20: continue

        ax = axes[i]
        ax.plot(w_rest[mask], fl[mask] / np.median(fl[mask]), 'k-', lw=0.5)
        ax.axvline(LAMBDA_VAC['Hβ'], color='b', ls='--', lw=0.8, alpha=0.7)
        ax.axvline(LAMBDA_VAC['[OIII]4959'], color='g', ls='--', lw=0.8, alpha=0.7)
        ax.axvline(LAMBDA_VAC['[OIII]5007'], color='r', ls='--', lw=0.8, alpha=0.7)
        ax.set_title(f'ID {sid}  z={z:.3f}', fontsize=9)
        ax.set_xlim(4840, 5030)
        ax.grid(alpha=0.2)
        ax.tick_params(labelsize=7)

    plt.tight_layout()
    out = ANALYSIS / 'figures' / 'fig5_hb_profiles.png'
    plt.savefig(out, dpi=100, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out}")


if __name__ == '__main__':
    main()
