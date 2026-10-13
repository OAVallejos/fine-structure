# P01_load.py
"""
P01: Load and clean x1d.
- Check units (µm → Å)
- Clean NaN/Inf
- Sort by wavelength
- Output: p01_clean_spectra.csv
"""

import numpy as np
import pandas as pd
from astropy.io import fits
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from config import *


def load_x1d(f):
    """Load an x1d. Return (w, fl, er, grating, filt) or None."""
    try:
        with fits.open(f) as hdul:
            hdr0 = hdul[0].header
            d = hdul[1].data
            cn = hdul[1].columns.names
            wc = next((c for c in cn if 'wave' in c.lower()), cn[0])
            fc = next((c for c in cn if 'flux' in c.lower() and 'err' not in c.lower()), cn[1])
            ec = next((c for c in cn if 'err' in c.lower()), cn[2] if len(cn) > 2 else cn[1])
            w = d[wc].copy()
            fl = d[fc].copy()
            er = d[ec].copy()
            grating = hdr0.get('GRATING', 'UNKNOWN')
            filt = hdr0.get('FILTER', 'UNKNOWN')
    except Exception:
        return None

    # Units: µm → Å
    if np.nanmax(w) < 100:
        w = w * 1e4

    # Clean
    good = np.isfinite(fl) & np.isfinite(er) & (er > 0) & np.isfinite(w)
    if good.sum() < 50:
        return None

    w = w[good]
    fl = fl[good]
    er = er[good]

    # Sort by wavelength
    idx = np.argsort(w)
    w = w[idx]
    fl = fl[idx]
    er = er[idx]

    return w, fl, er, grating, filt


def main():
    print("=" * 70)
    print("P01: LOAD AND CLEAN")
    print("=" * 70)

    files = sorted(SPECTRA.glob('*.fits'))
    print(f"Files in {SPECTRA}: {len(files)}")

    results = []
    n_fail = 0
    for f in files:
        r = load_x1d(f)
        if r is None:
            n_fail += 1
            continue
        w, fl, er, grating, filt = r
        results.append({
            'file': f.name,
            'wave_min': float(w.min()),
            'wave_max': float(w.max()),
            'n_pixels': int(len(w)),
            'grating': grating,
            'filter': filt,
        })

    df = pd.DataFrame(results)
    print(f"\nValid files: {len(df)}")
    print(f"Discarded files: {n_fail}")
    print()
    print("Grating:")
    print(df['grating'].value_counts())
    print()
    print("Wavelength range (observed):")
    print(f"  global min:   {df['wave_min'].min():.0f} Å")
    print(f"  global max:   {df['wave_max'].max():.0f} Å")
    print(f"  median min:   {df['wave_min'].median():.0f} Å")
    print(f"  median max:   {df['wave_max'].median():.0f} Å")
    print()
    print("N pixels:")
    print(f"  median: {df['n_pixels'].median():.0f}")
    print(f"  min:    {df['n_pixels'].min()}")
    print(f"  max:    {df['n_pixels'].max()}")
    print()

    out = ANALYSIS / 'p01_clean_spectra.csv'
    df.to_csv(out, index=False)
    print(f"Saved: {out}")


if __name__ == '__main__':
    main()
