python3 << 'EOF'
import numpy as np
from astropy.io import fits
from scipy.optimize import curve_fit
import os

ANALYSIS = os.path.expanduser('~/lrd_workspace/analysis')

with fits.open(os.path.join(ANALYSIS, 'stacked_tierA_final.fits')) as h:
    t = h[1].data
    wave, flux, err = t['WAVE_REST'], t['FLUX'], t['FLUX_ERR']

SCALE = 1e-20
flux, err = flux/SCALE, err/SCALE

mask = np.isfinite(flux) & np.isfinite(err) & (err > 0)
wave, flux, err = wave[mask], flux[mask], err[mask]

fit_mask = (wave > 4830) & (wave < 5040)
w_fit, f_fit, e_fit = wave[fit_mask], flux[fit_mask], err[fit_mask]

LAM_HB, LAM_4959, LAM_5007 = 4862.68, 4960.30, 5008.24

def model(w, a0, a1, A_hb, mu_hb, sig_hb, A_4959, mu_4959, sig_4959, A_5007, mu_5007, sig_5007):
    return (a0 + a1*(w-4900)
            + A_hb   * np.exp(-0.5*((w-mu_hb)/sig_hb)**2)
            + A_4959 * np.exp(-0.5*((w-mu_4959)/sig_4959)**2)
            + A_5007 * np.exp(-0.5*((w-mu_5007)/sig_5007)**2))

# Fit of the real stack
p0 = [0.06, 0.0, 1.0, LAM_HB, 2.5, 1.6, LAM_4959, 2.7, 3.7, LAM_5007, 2.7]
popt, pcov = curve_fit(model, w_fit, f_fit, p0=p0, sigma=e_fit, maxfev=50000, absolute_sigma=True)

print(f"Fit of the real stack:")
print(f"  mu_Hb = {popt[3]:.4f}, mu_4959 = {popt[6]:.4f}, mu_5007 = {popt[9]:.4f}")

# Inject shift equivalent to Δα/α
injected_values = [1e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3]
results = []

for da_inj in injected_values:
    # Shift in all lines: Δλ/λ = -2 Δα/α
    shift_frac = -2.0 * da_inj
    mu_hb_inj   = popt[3] * (1 + shift_frac)
    mu_4959_inj = popt[6] * (1 + shift_frac)
    mu_5007_inj = popt[9] * (1 + shift_frac)
    
    # Mock with shifted lines
    f_mock = model(w_fit, popt[0], popt[1],
                   popt[2],  mu_hb_inj,   popt[4],
                   popt[5],  mu_4959_inj, popt[7],
                   popt[8],  mu_5007_inj, popt[10])
    
    # Fit of the mock
    p0_mock = [popt[0], popt[1], popt[2], mu_hb_inj, popt[4],
               popt[5], mu_4959_inj, popt[7], popt[8], mu_5007_inj, popt[10]]
    try:
        popt_mock, pcov_mock = curve_fit(model, w_fit, f_mock, p0=p0_mock,
                                          sigma=e_fit, maxfev=50000, absolute_sigma=True)
        perr_mock = np.sqrt(np.diag(pcov_mock))
        
        # Recover Δα/α from each line
        da_hb   = -0.5 * (popt_mock[3] - popt[3]) / popt[3]
        da_4959 = -0.5 * (popt_mock[6] - popt[6]) / popt[6]
        da_5007 = -0.5 * (popt_mock[9] - popt[9]) / popt[9]
        
        # Weighted average
        da_arr = np.array([da_hb, da_4959, da_5007])
        err_arr = np.array([0.5*perr_mock[3]/popt[3],
                            0.5*perr_mock[6]/popt[6],
                            0.5*perr_mock[9]/popt[9]])
        w = 1/err_arr**2
        da_mean = np.sum(w*da_arr)/np.sum(w)
        da_err  = 1/np.sqrt(np.sum(w))
        
        results.append((da_inj, da_mean, da_err))
        print(f"\nΔα/α injected = {da_inj:.1e}")
        print(f"  Recovered: {da_mean:+.3e} ± {da_err:.3e}")
        print(f"  Bias: {(da_mean - da_inj):.2e}  ({abs(da_mean-da_inj)/da_err:.2f}σ)")
    except Exception as e:
        print(f"Δα/α = {da_inj:.1e}: failed ({e})")

# Save CSV
import pandas as pd
df = pd.DataFrame(results, columns=['da_injected', 'da_recovered', 'da_error'])
df.to_csv(os.path.join(ANALYSIS, 'forecast_sensitivity_curve.csv'), index=False)
print(f"\n💾 {os.path.join(ANALYSIS, 'forecast_sensitivity_curve.csv')}")
EOF