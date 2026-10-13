"""
P08: Quantification of systematics.
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


def main():
    print("=" * 70)
    print("P08: SYSTEMATICS")
    print("=" * 70)

    # 1. z_Spec bias
    print("\n[1] z_Spec bias:")
    rz = pd.read_csv(OUT_REFINED_Z)
    print(f"  Mean Δz: {rz['dz_avg'].mean():+.6f} ± {rz['dz_avg'].std():.6f}")
    print(f"  Δz(Hβ):   {rz['dz_hb'].mean():+.6f} ± {rz['dz_hb'].std():.6f}")
    print(f"  Δz([OIII]): {rz['dz_oiii'].mean():+.6f} ± {rz['dz_oiii'].std():.6f}")
    print(f"  Δz(Hβ) − Δz([OIII]): {(rz['dz_hb'] - rz['dz_oiii']).mean():+.6f} "
          f"± {(rz['dz_hb'] - rz['dz_oiii']).std():.6f}")

    # 2. Calibration (residual offset from the fit)
    print("\n[2] NIRSpec calibration:")
    fit = pd.read_csv(OUT_FIT)
    print(fit.to_string(index=False))

    # 3. LSF (instrumental width)
    print("\n[3] LSF (line width):")
    stack = fits.open(OUT_STACK)
    wave = stack[1].data['WAVE_REST']
    flux = stack[1].data['FLUX']
    err = stack[1].data['FLUX_ERR']
    cont = np.median(flux[(wave > 5200) & (wave < 5400)])
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
                           absolute_sigma=True, maxfev=100000)
    perr = np.sqrt(np.diag(pcov))
    for name, sig_idx, lam_lab in [
        ('Hβ', 4, LAMBDA_VAC['Hβ']),
        ('[OIII]4959', 7, LAMBDA_VAC['[OIII]4959']),
        ('[OIII]5007', 10, LAMBDA_VAC['[OIII]5007']),
        ('Hα', 13, LAMBDA_VAC['Hα']),
    ]:
        sig = popt[sig_idx]
        # σ_obs² = σ_intr² + σ_LSF²  → σ_LSF = sqrt(σ_obs² - σ_intr²)
        # Approximation: σ_intr ~ 1.5 Å (AGN lines)
        sig_lsf = np.sqrt(max(sig**2 - 1.5**2, 0))
        print(f"  {name:12s}  σ_obs = {sig:.3f} Å  →  σ_LSF ~ {sig_lsf:.3f} Å")

    # 4. Blending ([OIII] 4959)
    print("\n[4] Blending [OIII] 4959:")
    ratio_5007_4959 = popt[8] / popt[5]
    print(f"  Ratio 5007/4959 = {ratio_5007_4959:.3f}")
    print(f"  Theoretical (low density) = 2.98")
    print(f"  Deviation: {(ratio_5007_4959 - 2.98) / 2.98 * 100:+.1f}%")

    # 5. Cosmic variance
    print("\n[5] Cosmic variance (between objects):")
    # Use the dispersion of Δz as a proxy
    print(f"  σ(Δz) = {rz['dz_avg'].std():.6f}")
    print(f"  Contribution to σ(Δα/α) from z: ~ {rz['dz_avg'].std():.2e}")
    print(f"  (without stacking correction)")

    # Save
    results = {
        'dz_mean': rz['dz_avg'].mean(),
        'dz_std': rz['dz_avg'].std(),
        'dz_hb_mean': rz['dz_hb'].mean(),
        'dz_oiii_mean': rz['dz_oiii'].mean(),
        'dz_diff_mean': (rz['dz_hb'] - rz['dz_oiii']).mean(),
        'dz_diff_std': (rz['dz_hb'] - rz['dz_oiii']).std(),
        'ratio_5007_4959': float(ratio_5007_4959),
        'sigma_hb': float(popt[4]),
        'sigma_4959': float(popt[7]),
        'sigma_5007': float(popt[10]),
        'sigma_ha': float(popt[13]),
    }
    pd.DataFrame([results]).to_csv(OUT_SYSTEMATICS, index=False)
    print(f"\nSaved: {OUT_SYSTEMATICS}")


if __name__ == '__main__':
    main()
