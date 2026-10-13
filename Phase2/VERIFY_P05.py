# VERIFY_P05.py
import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.optimize import curve_fit

stack = fits.open('/root/lrd_workspace/analysis/stack_v2.fits')
d = stack[1].data
wave = d['WAVE_REST']; flux = d['FLUX']; err = d['FLUX_ERR']

LAB = {'Hβ': 4862.68, 'OIII4959': 4960.30, 'OIII5007': 5008.24, 'Hα': 6564.63}

def model(w, a0, a1, A_hb, mu_hb, sig_hb, A_4959, mu_4959, sig_4959,
          A_5007, mu_5007, sig_5007, A_ha, mu_ha, sig_ha_nii, A_nii):
    cont = a0 + a1 * (w - 5700)
    D_NII = 6585.27 - 6564.63
    g_hb   = A_hb   * np.exp(-0.5 * ((w - mu_hb)   / sig_hb)**2)
    g_4959 = A_4959 * np.exp(-0.5 * ((w - mu_4959) / sig_4959)**2)
    g_5007 = A_5007 * np.exp(-0.5 * ((w - mu_5007) / sig_5007)**2)
    g_ha   = A_ha   * np.exp(-0.5 * ((w - mu_ha)   / sig_ha_nii)**2)
    g_nii  = A_nii  * np.exp(-0.5 * ((w - (mu_ha + D_NII)) / sig_ha_nii)**2)
    return cont + g_hb + g_4959 + g_5007 + g_ha + g_nii

cont = np.median(flux[(wave > 5200) & (wave < 5400)])
def amp(lam, w=5):
    m = (wave > lam - w) & (wave < lam + w)
    return max(np.max(flux[m]) - cont, abs(cont) * 0.5)

p0 = [cont, 0.0,
      amp(4862.68), 4862.68, 2.5,
      amp(4960.30), 4960.30, 2.7,
      amp(5008.24), 5008.24, 2.7,
      amp(6564.63, 10), 6564.63, 3.0, amp(6585.27, 10) * 0.3]

popt, pcov = curve_fit(model, wave, flux, p0=p0, sigma=err,
                       absolute_sigma=True, maxfev=200000)
perr = np.sqrt(np.diag(pcov))

# Linearity test
lam = np.array([popt[3], popt[6], popt[9], popt[12]])
dlam = lam - np.array([LAB['Hβ'], LAB['OIII4959'], LAB['OIII5007'], LAB['Hα']])
elam = np.array([perr[3], perr[6], perr[9], perr[12]])

A = np.vstack([lam, np.ones_like(lam)]).T
W = np.diag(1.0/elam**2)
coef = np.linalg.solve(A.T @ W @ A, A.T @ W @ dlam)
resid = dlam - A @ coef
chi2 = np.sum((resid/elam)**2)

print("=" * 70)
print("P05 VERIFICATION")
print("=" * 70)
print(f"\nOffset vs λ:")
print(f"  Slope:     {coef[0]:.3e} Å/Å")
print(f"  Intercept: {coef[1]:+.4f} Å")
print(f"  χ²_red = {chi2/2:.3f}")
print()
for i, name in enumerate(['Hβ', 'OIII4959', 'OIII5007', 'Hα']):
    print(f"  {name:12s}  λ={lam[i]:.2f}  Δλ={dlam[i]:+.4f}  resid={resid[i]:+.4f}  ({resid[i]/elam[i]:+.2f}σ)")

print()
if chi2/2 < 2:
    print("→ The 4 lines follow a linear trend")
    print("→ NON-LINEAR NIRSpec calibration")
else:
    print("→ The 4 lines do NOT follow a linear trend")
    print("→ [OIII] 5007 is an outlier → something else")
