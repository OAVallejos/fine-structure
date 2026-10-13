# P20_lrd_profiles.py
"""
P20: Visualize Hα in individual LRDs.
"""
import numpy as np
import pandas as pd
from astropy.io import fits
from pathlib import Path
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).parent))
from config import *

SPECTRA_LRD = ROOT / 'spectra_lrd_candidates'

def main():
    df = pd.read_csv(ANALYSIS / 'lrd_z4.csv')
    df['ID'] = df['ID'].astype(str).str.zfill(8)
    df = df[df['snr_peak'] > 10]

    np.random.seed(42)
    sample = df.sample(min(12, len(df)))

    fig, axes = plt.subplots(4, 3, figsize=(15, 12))
    axes = axes.flatten()

    for i, (_, row) in enumerate(sample.iterrows()):
        z = row['z']
        f = SPECTRA_LRD / row['file']
        if not f.exists(): continue

        try:
            with fits.open(f, memmap=False) as hdul:
                d = hdul[1].data
                cn = hdul[1].columns.names
                wc = next((c for c in cn if 'wave' in c.lower()), cn[0])
                fc = next((c for c in cn if 'flux' in c.lower() and 'err' not in c.lower()), cn[1])
                w = d[wc].copy(); fl = d[fc].copy()
            if np.nanmax(w) < 100: w = w * 1e4
            w_rest = w / (1 + z)
            good = np.isfinite(fl)
            w_rest = w_rest[good]; fl = fl[good]
        except Exception:
            continue

        # Hα + [NII] window
        mask = (w_rest > 6400) & (w_rest < 6800)
        if mask.sum() < 20: continue

        ax = axes[i]
        ax.plot(w_rest[mask], fl[mask] / np.median(fl[mask]), 'k-', lw=0.5)
        ax.axvline(LAMBDA_VAC['Hα'], color='r', ls='--', lw=0.8, alpha=0.7)
        ax.axvline(LAMBDA_VAC['[NII]6584'], color='b', ls='--', lw=0.8, alpha=0.7)
        ax.set_title(f'ID {row["ID"]}  z={z:.3f}', fontsize=9)
        ax.set_xlim(6400, 6800)
        ax.grid(alpha=0.2)
        ax.tick_params(labelsize=7)

    plt.tight_layout()
    out = ANALYSIS / 'figures' / 'fig7_lrd_ha_profiles.png'
    plt.savefig(out, dpi=100, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out}")


if __name__ == '__main__':
    main()
