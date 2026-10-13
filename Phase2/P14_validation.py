"""
P14: Validation with independent samples.
- Field (GS vs GN)
- Grating (G395M vs G395H)
- Class (AGN vs starburst)
- Bootstrap split (half vs half)
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


def load_spectrum(f, z):
    try:
        with fits.open(SPECTRA / f, memmap=False) as hdul:
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
            return None
        return w_rest[good], fl[good], er[good]
    except Exception:
        return None


def build_stack(df, z_map, grid):
    num = np.zeros_like(grid); den = np.zeros_like(grid)
    n_ok = 0
    for _, row in df.iterrows():
        sid = str(row['ID']).zfill(8)
        if sid not in z_map: continue
        z = z_map[sid]
        if z <= 0: continue
        lam_ha_obs = LAMBDA_VAC['Hα'] * (1 + z)
        if not (X1D_WAVE_MIN < lam_ha_obs < X1D_WAVE_MAX): continue
        spec = load_spectrum(row['file'], z)
        if spec is None: continue
        w_rest, fl, er = spec
        f_int = np.interp(grid, w_rest, fl, left=0, right=0)
        e_int = np.interp(grid, w_rest, er, left=np.inf, right=np.inf)
        wgt = 1.0 / e_int**2
        wgt[~np.isfinite(wgt)] = 0
        num += f_int * wgt; den += wgt
        n_ok += 1
    sf = np.where(den > 0, num / den, 0)
    se = np.where(den > 0, 1.0 / np.sqrt(den), np.inf)
    return sf, se, n_ok


def fit_alpha(df, z_map, grid, label):
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
        print(f"  {label:25s}  — fit failed: {ex}")
        return None
    finally:
        del sf, se; gc.collect()


def main():
    print("=" * 90)
    print("P14: VALIDATION WITH INDEPENDENT SAMPLES")
    print("=" * 90)

    grid = np.linspace(GRID_MIN, GRID_MAX, GRID_N)

    # Load multi-line z
    mz = pd.read_csv(ANALYSIS / 'z_multiline.csv')
    mz['ID'] = mz['ID'].astype(str).str.zfill(8)
    z_map = dict(zip(mz['ID'], mz['z_multiline']))

    # Load Tier A
    tier_a = pd.read_csv(ANALYSIS / 'tier_A_z5_v2.csv')
    tier_a['ID'] = tier_a['ID'].astype(str).str.zfill(8)

    # Load classification
    print("\n[1] Sample classification...")
    mf = pd.read_csv(ANALYSIS / 'master_features_v4.csv')
    mf['NIRSpec_ID_str'] = mf['NIRSpec_ID'].astype(str).str.zfill(8)
    mf['ratio_OIII_HB'] = mf['O3_5007_flux'] / mf['HB_4861_flux']
    def classify(r):
        if pd.isna(r): return 'unknown'
        if r > 5: return 'AGN'
        elif r > 3: return 'composite'
        else: return 'starburst'
    mf['class'] = mf['ratio_OIII_HB'].apply(classify)
    tier_a = tier_a.merge(mf[['NIRSpec_ID_str', 'class']],
                          left_on='ID', right_on='NIRSpec_ID_str', how='left')

    # ============================================================
    # [1] Field: GS vs GN
    # ============================================================
    print("\n[1] Validation by field (independent):")
    for field in ['goods-s', 'goods-n']:
        sub = tier_a[tier_a['field'] == field]
        fit_alpha(sub, z_map, grid, field.upper())

    # ============================================================
    # [2] Grating: G395M vs G395H
    # ============================================================
    print("\n[2] Validation by grating (independent):")
    for grat in ['G395M', 'G395H']:
        sub = tier_a[tier_a['grating'] == grat]
        fit_alpha(sub, z_map, grid, grat)

    # ============================================================
    # [3] Class: AGN vs starburst
    # ============================================================
    print("\n[3] Validation by class (independent):")
    for cls in ['AGN', 'composite', 'starburst']:
        sub = tier_a[tier_a['class'] == cls]
        fit_alpha(sub, z_map, grid, f'class={cls}')

    # ============================================================
    # [4] Bootstrap split (random)
    # ============================================================
    print("\n[4] Bootstrap split (half vs half, 3 splits):")
    for i in range(3):
        np.random.seed(i)
        idx = np.random.permutation(len(tier_a))
        half = len(tier_a) // 2
        sub1 = tier_a.iloc[idx[:half]]
        sub2 = tier_a.iloc[idx[half:]]
        fit_alpha(sub1, z_map, grid, f'split{i+1}a')
        fit_alpha(sub2, z_map, grid, f'split{i+1}b')

    print("\n" + "=" * 90)
    print("END OF VALIDATION")
    print("=" * 90)


if __name__ == '__main__':
    main()
