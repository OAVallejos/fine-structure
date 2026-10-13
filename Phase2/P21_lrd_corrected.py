# P21_lrd_corrected.py
"""
P21: LRD stack with corrected z.
"""
import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.optimize import curve_fit
from astropy.table import Table
from pathlib import Path
import sys, gc
sys.path.insert(0, str(Path(__file__).parent))
from config import *

SPECTRA_LRD = ROOT / 'spectra_lrd_candidates'
GRID = np.linspace(6400, 6800, 500)
DZ_CORR = 0.0014  # measured offset

def model_ha_nii(x, A_ha, mu_ha, sig, A_nii, c):
    D_NII = LAMBDA_VAC['[NII]6584'] - LAMBDA_VAC['Hα']
    g_ha = A_ha * np.exp(-0.5 * ((x - mu_ha) / sig)**2)
    g_nii = A_nii * np.exp(-0.5 * ((x - (mu_ha + D_NII)) / sig)**2)
    return g_ha + g_nii + c

def main():
    print("=" * 70)
    print("P21: LRD STACK WITH CORRECTED Z")
    print("=" * 70)
    print(f"Applied correction: Δz = +{DZ_CORR}")

    df = pd.read_csv(ANALYSIS / 'lrd_z4.csv')
    df['ID'] = df['ID'].astype(str).str.zfill(8)
    df = df[df['snr_peak'] > 10]
    print(f"LRDs: {len(df)}")

    num = np.zeros_like(GRID); den = np.zeros_like(GRID)
    n_ok = 0

    for _, row in df.iterrows():
        z = row['z'] + DZ_CORR
        f = SPECTRA_LRD / row['file']
        if not f.exists(): continue
        try:
            with fits.open(f, memmap=False) as hdul:
                d = hdul[1].data
                cn = hdul[1].columns.names
                wc = next((c for c in cn if 'wave' in c.lower()), cn[0])
                fc = next((c for c in cn if 'flux' in c.lower() and 'err' not in c.lower()), cn[1])
                ec = next((c for c in cn if 'err' in c.lower()), cn[2])
                w = d[wc].copy(); fl = d[fc].copy(); er = d[ec].copy()
            if np.nanmax(w) < 100: w = w * 1e4
            w_rest = w / (1 + z)
            good = np.isfinite(fl) & np.isfinite(er) & (er > 0)
            if good.sum() < 50: continue
            f_int = np.interp(GRID, w_rest[good], fl[good], left=0, right=0)
            e_int = np.interp(GRID, w_rest[good], er[good], left=np.inf, right=np.inf)
            wgt = 1.0 / e_int**2
            wgt[~np.isfinite(wgt)] = 0
            num += f_int * wgt; den += wgt
            n_ok += 1
        except Exception:
            continue
        gc.collect()

    sf = np.where(den > 0, num / den, 0)
    se = np.where(den > 0, 1.0 / np.sqrt(den), np.inf)
    print(f"Stack: {n_ok} LRDs")

    # Save
    t = Table()
    t['WAVE_REST'] = GRID
    t['FLUX'] = sf
    t['FLUX_ERR'] = se
    t.write(ANALYSIS / 'stack_lrd_corrected.fits', overwrite=True)
    print(f"Saved: stack_lrd_corrected.fits")

    # Fit
    print()
    print("=" * 70)
    print("FIT Hα + [NII] (corrected z)")
    print("=" * 70)

    cont = np.median(sf[(GRID > 6700) & (GRID < 6800)])
    def amp(lam, w=5):
        m = (GRID > lam - w) & (GRID < lam + w)
        return max(np.max(sf[m]) - cont, abs(cont) * 0.5)

    p0 = [amp(LAMBDA_VAC['Hα'], 10), LAMBDA_VAC['Hα'], 5.0,
          amp(LAMBDA_VAC['[NII]6584'], 10) * 0.3, cont]

    bounds = (
        [0, LAMBDA_VAC['Hα'] - 10, 1.5, 0, cont * 0.5],
        [np.inf, LAMBDA_VAC['Hα'] + 10, 20, np.inf, cont * 1.5]
    )

    try:
        popt, pcov = curve_fit(model_ha_nii, GRID, sf, p0=p0, sigma=se,
                               absolute_sigma=True, bounds=bounds, maxfev=100000)
        perr = np.sqrt(np.diag(pcov))
        chi2 = np.sum(((sf - model_ha_nii(GRID, *popt)) / se)**2)
        dof = len(GRID) - len(popt)

        print(f"χ²_red = {chi2/dof:.3f}")
        print(f"Hα:       A={popt[0]:.3e}, μ={popt[1]:.3f}, σ={popt[2]:.3f}")
        print(f"[NII]:    A={popt[3]:.3e}")
        print()
        d_ha = popt[1] - LAMBDA_VAC['Hα']
        print(f"Hα:  Δλ = {d_ha:+.4f} ± {perr[1]:.4f} Å  ({d_ha/perr[1]:+.2f}σ)")
        print(f"Ratio NII/Hα = {popt[3]/popt[0]:.3f}")

        # Δα/α (assuming NII at theoretical position)
        q_ha = 0.5
        dalpha = d_ha / (q_ha * LAMBDA_VAC['Hα'])
        e_dalpha = perr[1] / (q_ha * LAMBDA_VAC['Hα'])
        print(f"\nΔα/α = {dalpha*1e6:+.2f} ± {e_dalpha*1e6:.2f} ppm")
        if abs(dalpha/e_dalpha) < 2:
            print(f"→ CONSISTENT WITH ZERO")
            print(f"→ 2σ limit: |Δα/α| < {2*e_dalpha*1e6:.2f} ppm")
        else:
            print(f"→ {dalpha/e_dalpha:+.1f}σ")
    except Exception as ex:
        print(f"FIT FAILED: {ex}")


if __name__ == '__main__':
    main()
