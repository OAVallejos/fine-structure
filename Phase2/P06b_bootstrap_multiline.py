"""
P06: Bootstrap of the α test (memory-safe).
- Resamples in batches
- Frees memory after each batch
- Saves results incrementally
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
            # Free
            del w, fl, er, w_rest, good, f_int, e_int, wgt
        except Exception:
            continue
    sf = np.where(den > 0, num / den, 0)
    se = np.where(den > 0, 1.0 / np.sqrt(den), np.inf)
    return sf, se, n_ok


def fit_stack(sf, se, grid):
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
        if popt[8] / popt[5] < 1.5 or popt[8] / popt[5] > 4:
            return None
        return popt, perr
    except Exception:
        return None


def main():
    print("=" * 70)
    print("P06: BOOTSTRAP (memory-safe)")
    print("=" * 70)

    grid = np.linspace(GRID_MIN, GRID_MAX, GRID_N)

    # Load corrected z
    rz = pd.read_csv(ANALYSIS / "z_multiline.csv")
    rz['ID'] = rz['ID'].astype(str).str.zfill(8)
    z_map = dict(zip(rz['ID'], rz['z_multiline']))
    print(f"Objects with corrected z: {len(z_map)}")

    # Load tier A
    tier_a = pd.read_csv(ANALYSIS / 'tier_A_z5_v2.csv')
    tier_a['ID'] = tier_a['ID'].astype(str).str.zfill(8)
    print(f"Tier A: {len(tier_a)}")

    # Save the df indexed to avoid recomputing
    files_list = tier_a[['ID', 'file']].values.tolist()

    # Original stack
    print("\nOriginal stack...")
    sf, se, n_ok = build_stack(tier_a, z_map, grid)
    print(f"Stack: {n_ok} objects")
    gc.collect()

    res = fit_stack(sf, se, grid)
    if res is None:
        print("FIT FAILED")
        return
    popt0, perr0 = res

    d_hb_0 = popt0[3] - LAMBDA_VAC['Hβ']
    d_5007_0 = popt0[9] - LAMBDA_VAC['[OIII]5007']
    diff_0 = d_hb_0 - d_5007_0
    dalpha_0 = -diff_0 / (2 * LAMBDA_VAC['Hβ'])

    print(f"\nORIGINAL STACK:")
    print(f"  Hβ       Δλ = {d_hb_0:+.4f} ± {perr0[3]:.4f} Å")
    print(f"  OIII5007 Δλ = {d_5007_0:+.4f} ± {perr0[9]:.4f} Å")
    print(f"  Difference = {diff_0:+.4f} Å")
    print(f"  Δα/α = {dalpha_0:+.3e}  ({dalpha_0*1e6:+.2f} ppm)")

    # Bootstrap in batches
    print(f"\nBootstrap ({N_BOOTSTRAP} iter, batches of {N_BATCH})...")
    dalpha_boot = []
    n_total = len(tier_a)
    files_arr = np.array([r[1] for r in files_list])
    ids_arr = np.array([r[0] for r in files_list])

    for batch_start in range(0, N_BOOTSTRAP, N_BATCH):
        batch_end = min(batch_start + N_BATCH, N_BOOTSTRAP)
        for i in range(batch_start, batch_end):
            idx = np.random.RandomState(i).choice(n_total, size=n_total, replace=True)
            sample = pd.DataFrame({'ID': ids_arr[idx], 'file': files_arr[idx]})
            sf_b, se_b, n_b = build_stack(sample, z_map, grid)
            r = fit_stack(sf_b, se_b, grid)
            if r is None:
                continue
            popt, _ = r
            d_hb = popt[3] - LAMBDA_VAC['Hβ']
            d_5007 = popt[9] - LAMBDA_VAC['[OIII]5007']
            diff = d_hb - d_5007
            dalpha = -diff / (2 * LAMBDA_VAC['Hβ'])
            dalpha_boot.append(dalpha)
            del sf_b, se_b, sample
        gc.collect()
        print(f"  Batch {batch_start//N_BATCH + 1}: {len(dalpha_boot)} accumulated")

    dalpha_boot = np.array(dalpha_boot)
    print(f"\nBOOTSTRAP ({len(dalpha_boot)} valid):")
    print(f"  Mean     = {dalpha_boot.mean():+.3e}  ({dalpha_boot.mean()*1e6:+.2f} ppm)")
    print(f"  Median   = {np.median(dalpha_boot):+.3e}  ({np.median(dalpha_boot)*1e6:+.2f} ppm)")
    print(f"  Std      = {dalpha_boot.std():.3e}  ({dalpha_boot.std()*1e6:.2f} ppm)")
    print(f"  2σ       = [{np.percentile(dalpha_boot, 2.5):+.3e}, {np.percentile(dalpha_boot, 97.5):+.3e}]")

    pd.DataFrame({'dalpha_over_alpha': dalpha_boot}).to_csv(OUT_BOOTSTRAP, index=False)
    print(f"\nSaved: {OUT_BOOTSTRAP}")

    mean = dalpha_boot.mean()
    std = dalpha_boot.std()
    print()
    print("=" * 70)
    print("VERDICT")
    print("=" * 70)
    if abs(mean / std) < 2:
        print(f"  Δα/α = {mean*1e6:+.2f} ± {std*1e6:.2f} ppm")
        print(f"  CONSISTENT WITH ZERO")
        print(f"  2σ limit: |Δα/α| < {2*std:.3e}")
    else:
        print(f"  Δα/α = {mean*1e6:+.2f} ± {std*1e6:.2f} ppm")
        print(f"  Residual: {mean/std:+.2f}σ")


if __name__ == '__main__':
    main()
