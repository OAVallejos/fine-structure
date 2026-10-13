"""
P07: α test in subsamples.
- GS vs GN
- z bins
- S/N bins
- Without outliers
- Memory-safe (Termux)
"""

import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.optimize import curve_fit
from pathlib import Path
import sys, gc
sys.path.insert(0, str(Path(__file__).parent))
from config import *


def model(w, a0, a1,
          A_hb, mu_hb, sig_hb,
          A_4959, mu_4959, sig_4959,
          A_5007, mu_5007, sig_5007,
          A_ha, mu_ha, sig_ha_nii, A_nii):
    cont = a0 + a1 * (w - 5700)
    D_NII = LAMBDA_VAC['[NII]6584'] - LAMBDA_VAC['Hα']
    g_hb   = A_hb   * np.exp(-0.5 * ((w - mu_hb)   / sig_hb)**2)
    g_4959 = A_4959 * np.exp(-0.5 * ((w - mu_4959) / sig_4959)**2)
    g_5007 = A_5007 * np.exp(-0.5 * ((w - mu_5007) / sig_5007)**2)
    g_ha   = A_ha   * np.exp(-0.5 * ((w - mu_ha)   / sig_ha_nii)**2)
    g_nii  = A_nii  * np.exp(-0.5 * ((w - (mu_ha + D_NII)) / sig_ha_nii)**2)
    return cont + g_hb + g_4959 + g_5007 + g_ha + g_nii


def build_stack(df, z_map, grid):
    num = np.zeros_like(grid)
    den = np.zeros_like(grid)
    n_ok = 0
    for _, row in df.iterrows():
        sid = str(row['ID']).zfill(8)
        if sid not in z_map:
            continue
        z = z_map[sid]
        if z <= 0:
            continue
        lam_ha_obs = LAMBDA_VAC['Hα'] * (1 + z)
        if not (X1D_WAVE_MIN < lam_ha_obs < X1D_WAVE_MAX):
            continue
        try:
            with fits.open(SPECTRA / row['file'], memmap=False) as hdul:
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
                continue
            f_int = np.interp(grid, w_rest[good], fl[good], left=0, right=0)
            e_int = np.interp(grid, w_rest[good], er[good], left=np.inf, right=np.inf)
            wgt = 1.0 / e_int**2
            wgt[~np.isfinite(wgt)] = 0
            num += f_int * wgt
            den += wgt
            n_ok += 1
            del w, fl, er, w_rest, good, f_int, e_int, wgt
        except Exception:
            continue
    sf = np.where(den > 0, num / den, 0)
    se = np.where(den > 0, 1.0 / np.sqrt(den), np.inf)
    return sf, se, n_ok


def fit_alpha(df, label, z_map, grid):
    if len(df) < 20:
        print(f"  {label:25s}  N={len(df):4d}  — too few")
        return None
    sf, se, n_ok = build_stack(df, z_map, grid)
    if not np.any(np.isfinite(se)):
        print(f"  {label:25s}  N={len(df):4d}  — empty")
        return None
    cont = np.median(sf[(grid > 5200) & (grid < 5400)])
    def amp(lam, w=5):
        m = (grid > lam - w) & (grid < lam + w)
        return max(np.max(sf[m]) - cont, abs(cont) * 0.5)
    p0 = [cont, 0.0,
          amp(LAMBDA_VAC['Hβ']), LAMBDA_VAC['Hβ'], 2.5,
          amp(LAMBDA_VAC['[OIII]4959']), LAMBDA_VAC['[OIII]4959'], 2.7,
          amp(LAMBDA_VAC['[OIII]5007']), LAMBDA_VAC['[OIII]5007'], 2.7,
          amp(LAMBDA_VAC['Hα'], 10), LAMBDA_VAC['Hα'], 3.0,
          amp(LAMBDA_VAC['[NII]6584'], 10) * 0.3]
    try:
        popt, pcov = curve_fit(model, grid, sf, p0=p0, sigma=se,
                               absolute_sigma=True, maxfev=100000)
        perr = np.sqrt(np.diag(pcov))
        d_hb = popt[3] - LAMBDA_VAC['Hβ']
        d_5007 = popt[9] - LAMBDA_VAC['[OIII]5007']
        e_hb = perr[3]; e_5007 = perr[9]
        diff = d_hb - d_5007
        e_diff = np.sqrt(e_hb**2 + e_5007**2)
        dalpha = -diff / (2 * LAMBDA_VAC['Hβ'])
        e_dalpha = e_diff / (2 * LAMBDA_VAC['Hβ'])
        sig = diff / e_diff if e_diff > 0 else 0
        print(f"  {label:25s}  N={n_ok:4d}  "
              f"Δα/α = {dalpha*1e6:+8.4f} ± {e_dalpha*1e6:6.4f} ppm  "
              f"({sig:+.3f}σ)")
        return dalpha, e_dalpha
    except Exception as ex:
        print(f"  {label:25s}  N={len(df):4d}  — fit failed: {ex}")
        return None
    finally:
        del sf, se
        gc.collect()


def main():
    print("=" * 90)
    print("P07: SUBSAMPLES")
    print("=" * 90)

    grid = np.linspace(GRID_MIN, GRID_MAX, GRID_N)

    rz = pd.read_csv(OUT_REFINED_Z)
    rz['ID'] = rz['ID'].astype(str).str.zfill(8)
    z_map = dict(zip(rz['ID'], rz['z_Corrected']))

    tier_a = pd.read_csv(ANALYSIS / 'tier_A_z5_v2.csv')
    tier_a['ID'] = tier_a['ID'].astype(str).str.zfill(8)

    print(f"\n[1] By field:")
    for field in ['goods-s', 'goods-n']:
        sub = tier_a[tier_a['field'] == field]
        fit_alpha(sub, field.upper(), z_map, grid)
        gc.collect()

    print(f"\n[2] By z bins:")
    for zlo, zhi in [(5.0, 5.5), (5.5, 6.0), (6.0, 6.5), (6.5, 8.0)]:
        sub = tier_a[(tier_a['z'] >= zlo) & (tier_a['z'] < zhi)]
        fit_alpha(sub, f"z∈[{zlo},{zhi}]", z_map, grid)
        gc.collect()

    print(f"\n[3] By S/N quartiles:")
    qs = tier_a['snr_peak'].quantile([0, 0.25, 0.5, 0.75, 1.0]).values
    for i in range(4):
        sub = tier_a[(tier_a['snr_peak'] >= qs[i]) & (tier_a['snr_peak'] < qs[i+1])]
        fit_alpha(sub, f"S/N Q{i+1} [{qs[i]:.1f},{qs[i+1]:.1f}]", z_map, grid)
        gc.collect()

    print(f"\n[4] Without outliers (S/N > 100):")
    sub = tier_a[tier_a['snr_peak'] <= 100]
    fit_alpha(sub, "S/N ≤ 100", z_map, grid)
    gc.collect()

    print(f"\n[5] Full sample:")
    fit_alpha(tier_a, "ALL", z_map, grid)

    # Save
    print()
    print("Note: results shown on stdout, save to CSV if needed")


if __name__ == '__main__':
    main()
