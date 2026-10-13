# P23_lrd_debug.py
"""
P23: Debug of the LRD stack.
"""
import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.optimize import curve_fit
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from config import *


def model(w, A_ha, mu_ha, sig, A_nii6584, c):
    """Hα + [NII] 6548,6584 (ratio 2.94)."""
    D_6548 = 6549.86 - 6564.63
    D_6584 = 6585.27 - 6564.63
    A_nii6548 = A_nii6584 / 2.94

    g_ha = A_ha * np.exp(-0.5 * ((w - mu_ha) / sig)**2)
    g_nii6548 = A_nii6548 * np.exp(-0.5 * ((w - (mu_ha + D_6548)) / sig)**2)
    g_nii6584 = A_nii6584 * np.exp(-0.5 * ((w - (mu_ha + D_6584)) / sig)**2)
    return g_ha + g_nii6548 + g_nii6584 + c


def main():
    print("=" * 70)
    print("P23: DEBUG STACK LRDs")
    print("=" * 70)

    stack = fits.open(ANALYSIS / 'stack_lrd_corrected.fits')
    d = stack[1].data
    wave = d['WAVE_REST']
    flux = d['FLUX']
    err = d['FLUX_ERR']

    # Line range
    GRID = np.linspace(6500, 6650, 200)
    f_int = np.interp(GRID, wave, flux, left=0, right=0)
    e_int = np.interp(GRID, wave, err, left=np.inf, right=np.inf)

    cont = np.median(f_int[(GRID > 6620) & (GRID < 6650)])

    # Find the peak
    idx_peak = np.argmax(f_int)
    lam_peak = GRID[idx_peak]
    print(f"Observed peak at: {lam_peak:.2f} Å")
    print(f"Continuum: {cont:.3e}")
    print()

    # Test 1: peak is Hα
    z_if_ha = lam_peak / 6564.63 - 1
    print(f"If the peak is Hα: z = {z_if_ha:.5f}")
    print(f"  (original z_Spec: 4.170790, corrected z: 4.172190)")

    # Test 2: peak is [NII] 6548
    z_if_nii = lam_peak / 6549.86 - 1
    print(f"If the peak is [NII]6548: z = {z_if_nii:.5f}")

    # Test 3: peak is [NII] 6584
    z_if_nii2 = lam_peak / 6585.27 - 1
    print(f"If the peak is [NII]6584: z = {z_if_nii2:.5f}")
    print()

    # Fit with Hα as reference
    print("=" * 70)
    print("FIT with Hα + [NII] 6548,6584")
    print("=" * 70)

    def amp(lam, w=10):
        m = (GRID > lam - w) & (GRID < lam + w)
        return max(np.max(f_int[m]) - cont, abs(cont) * 0.5)

    p0 = [amp(6564.63, 20), 6564.63, 5.0, amp(6585.27, 10) * 0.3, cont]

    bounds = (
        [0, 6555, 2.0, 0, cont * 0.5],
        [np.inf, 6575, 20, np.inf, cont * 1.5]
    )

    try:
        popt, pcov = curve_fit(model, GRID, f_int, p0=p0, sigma=e_int,
                               absolute_sigma=True, bounds=bounds, maxfev=100000)
        perr = np.sqrt(np.diag(pcov))
        chi2 = np.sum(((f_int - model(GRID, *popt)) / e_int)**2)
        dof = len(GRID) - len(popt)

        print(f"χ²_red = {chi2/dof:.3f}")
        print()
        print(f"Hα:       A={popt[0]:.3e}, μ={popt[1]:.4f}, σ={popt[2]:.3f}")
        print(f"[NII]6584: A={popt[3]:.3e}")
        print(f"Continuum: c={popt[4]:.3e}")
        print()
        print(f"Hα offset: {popt[1] - 6564.63:+.4f} Å")
        print(f"Ratio [NII]6584/Hα: {popt[3]/popt[0]:.3f}")
    except Exception as ex:
        print(f"FIT FAILED: {ex}")


if __name__ == '__main__':
    main()
