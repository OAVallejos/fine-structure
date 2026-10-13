# P26_plot_lrd_examples.py
"""
P26: Plot of LRD examples with Hα + [NII].
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from config import *


def main():
    SPECTRA_DIR = ROOT / 'spectra_lrd_candidates'
    OUT_DIR = ANALYSIS / 'figures'
    OUT_DIR.mkdir(exist_ok=True)

    TARGETS = {
        '00109798': {'z': 4.524403, 'snr': 84.7},
        '00289178': {'z': 4.237908, 'snr': 155.1}
    }

    LAM_HA = 6564.63
    LAM_NII_6584 = 6585.27

    fig, axes = plt.subplots(2, 1, figsize=(10, 8))
    fig.suptitle('LRD Candidates: Hα + [NII] (G395M)', fontsize=14, fontweight='bold')

    for idx, (obj_id, props) in enumerate(TARGETS.items()):
        ax = axes[idx]
        z = props['z']

        files = list(SPECTRA_DIR.glob(f"*{obj_id}*.fits*"))
        if not files:
            print(f"⚠️  Not found: {obj_id}")
            continue

        print(f"File: {files[0].name}")
        with fits.open(files[0], memmap=False) as hdul:
            d = hdul[1].data
            cn = hdul[1].columns.names
            wc = next((c for c in cn if 'wave' in c.lower()), cn[0])
            fc = next((c for c in cn if 'flux' in c.lower() and 'err' not in c.lower()), cn[1])
            ec = next((c for c in cn if 'err' in c.lower()), cn[2])
            wave_obs = d[wc].copy()
            flux = d[fc].copy()
            err = d[ec].copy()

        if np.nanmax(wave_obs) < 100:
            wave_obs = wave_obs * 1e4

        ha_obs = LAM_HA * (1 + z)
        nii_obs = LAM_NII_6584 * (1 + z)

        ax.plot(wave_obs, flux, color='black', lw=0.8, label=f'z = {z:.3f}')
        ax.fill_between(wave_obs, flux - err, flux + err, color='gray', alpha=0.3, label='1σ')
        ax.axvline(ha_obs, color='red', ls='--', lw=1.5, label='Hα')
        ax.axvline(nii_obs, color='blue', ls=':', lw=1.5, label='[NII] 6584')

        ax.set_ylabel('Flux (W/m²)')
        ax.set_title(f'ID: {obj_id} (S/N ~ {props["snr"]:.0f})')
        ax.legend(loc='upper right', fontsize=9)
        ax.grid(alpha=0.3)

        delta = 1500
        ax.set_xlim(ha_obs - delta, ha_obs + delta)

    axes[-1].set_xlabel('Observed Wavelength (Å)')
    plt.tight_layout()
    out = OUT_DIR / 'fig_lrd_examples_halpha.png'
    plt.savefig(out, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out}")


if __name__ == '__main__':
    main()
