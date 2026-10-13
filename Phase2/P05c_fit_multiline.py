"""
P05: Multi-line fit of the stack v2.
- Hβ + [OIII] 4959 + [OIII] 5007 + Hα + [NII]
- Hα in VACUUM (6564.63)
- [NII] tied to Hα
- α test with Hβ vs [OIII] 5007
- Output: fit_results_v2.csv
"""

import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.optimize import curve_fit
from pathlib import Path
import sys
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


def main():
    print("=" * 70)
    print("P05: FIT MULTI-LINE")
    print("=" * 70)

    stack = fits.open(ANALYSIS / "stack_v2_multiline.fits")
    d = stack[1].data
    wave = d['WAVE_REST']
    flux = d['FLUX']
    err = d['FLUX_ERR']

    cont = np.median(flux[(wave > 5200) & (wave < 5400)])
    print(f"Estimated continuum: {cont:.3e}")
    print()

    def amp(lam, w=5):
        m = (wave > lam - w) & (wave < lam + w)
        return max(np.max(flux[m]) - cont, abs(cont) * 0.5)

    p0 = [cont, 0.0,
          amp(LAMBDA_VAC['Hβ']), LAMBDA_VAC['Hβ'], 2.5,
          amp(LAMBDA_VAC['[OIII]4959']), LAMBDA_VAC['[OIII]4959'], 2.7,
          amp(LAMBDA_VAC['[OIII]5007']), LAMBDA_VAC['[OIII]5007'], 2.7,
          amp(LAMBDA_VAC['Hα'], 10), LAMBDA_VAC['Hα'], 3.0,
          amp(LAMBDA_VAC['[NII]6584'], 10) * 0.3]

    popt, pcov = curve_fit(model, wave, flux, p0=p0, sigma=err,
                           absolute_sigma=True, maxfev=200000)
    perr = np.sqrt(np.diag(pcov))
    chi2 = np.sum(((flux - model(wave, *popt)) / err)**2)
    dof = len(wave) - len(popt)

    print(f"χ² = {chi2:.1f}  dof = {dof}  χ²_red = {chi2/dof:.3f}")
    print()

    print("=== OFFSETS ===")
    results = []
    for name, mu_idx, lam_lab in [
        ('Hβ', 3, LAMBDA_VAC['Hβ']),
        ('[OIII]4959', 6, LAMBDA_VAC['[OIII]4959']),
        ('[OIII]5007', 9, LAMBDA_VAC['[OIII]5007']),
        ('Hα', 12, LAMBDA_VAC['Hα']),
    ]:
        d_off = popt[mu_idx] - lam_lab
        sig = d_off / perr[mu_idx]
        results.append({
            'line': name, 'lambda': lam_lab,
            'lambda_fit': popt[mu_idx],
            'delta_lambda': d_off, 'err': perr[mu_idx],
            'significance': sig,
        })
        print(f"  {name:12s}  λ_fit={popt[mu_idx]:.3f}  Δλ={d_off:+.4f} ± {perr[mu_idx]:.4f} Å  ({sig:+.2f}σ)")

    print()
    print("=== AMPLITUDES ===")
    for name, A_idx in [('Hβ', 2), ('[OIII]4959', 5), ('[OIII]5007', 8), ('Hα', 11), ('[NII]6584', 14)]:
        print(f"  {name:12s}  A = {popt[A_idx]:.3e}")

    print()
    print("=== RATIOS ===")
    print(f"  [OIII]5007/[OIII]4959 = {popt[8]/popt[5]:.3f}")
    print(f"  [OIII]5007/Hβ         = {popt[8]/popt[2]:.3f}")
    print(f"  Hα/Hβ                 = {popt[11]/popt[2]:.3f}")
    print(f"  [NII]6584/Hα          = {popt[14]/popt[11]:.3f}")

    print()
    print("=== α TEST (Hβ vs [OIII] 5007) ===")
    d_hb = popt[3] - LAMBDA_VAC['Hβ']
    d_5007 = popt[9] - LAMBDA_VAC['[OIII]5007']
    e_hb = perr[3]; e_5007 = perr[9]
    diff = d_hb - d_5007
    e_diff = np.sqrt(e_hb**2 + e_5007**2)
    dalpha = -diff / (2 * LAMBDA_VAC['Hβ'])
    e_dalpha = e_diff / (2 * LAMBDA_VAC['Hβ'])
    print(f"  Δλ(Hβ)        = {d_hb:+.4f} ± {e_hb:.4f} Å")
    print(f"  Δλ([OIII]5007) = {d_5007:+.4f} ± {e_5007:.4f} Å")
    print(f"  Difference    = {diff:+.4f} ± {e_diff:.4f} Å  ({diff/e_diff:+.2f}σ)")
    print(f"  Δα/α = {dalpha:+.3e} ± {e_dalpha:.3e}  ({dalpha*1e6:+.2f} ± {e_dalpha*1e6:.2f} ppm)")

    # Save
    df = pd.DataFrame(results)
    df.to_csv(OUT_FIT, index=False)
    print(f"\nSaved: {OUT_FIT}")


if __name__ == '__main__':
    main()
