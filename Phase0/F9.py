python3 << 'EOF'
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.optimize import curve_fit, minimize
import os

BASE = os.path.expanduser('~/lrd_workspace')
stack_fits = os.path.join(BASE, 'analysis', 'stacked_tierA_final.fits')

# 1. Load the REAL stack you just generated
with fits.open(stack_fits) as h:
    data = h[1].data
    wave = data['WAVE_REST']
    flux_real = data['FLUX']
    err_real = data['FLUX_ERR']

print(f"📊 Loaded real stack. Typical continuum S/N: {1/np.median(err_real):.1f}")

# 2. Gaussian model for [O III] 5007
def gaussian_model(wave, amp, mu, sigma, continuum):
    return continuum + amp * np.exp(-0.5 * ((wave - mu) / sigma)**2)

# Fit the real line to get baseline parameters
mask_line = (wave > 4980) & (wave < 5030)
p0 = [np.max(flux_real[mask_line]), 5008.0, 3.0, np.median(flux_real)]
popt, pcov = curve_fit(gaussian_model, wave[mask_line], flux_real[mask_line], sigma=err_real[mask_line], p0=p0)
amp_true, mu_true, sigma_true, cont_true = popt
print(f"✅ Real [O III] 5007 line fitted: Center = {mu_true:.2f} Å, Sigma = {sigma_true:.2f} Å")

# 3. SIGNAL INJECTION: Simulate Delta_alpha/alpha = 5e-5
# Relation: Delta_lambda / lambda ≈ -2 * (Delta_alpha / alpha)
delta_alpha_injected = 5.0e-5
delta_lambda_frac = -2.0 * delta_alpha_injected
mu_injected = mu_true * (1.0 + delta_lambda_frac)

print(f"💉 Injecting shift equivalent to Δα/α = {delta_alpha_injected:.1e}")
print(f"   Expected shift: {(mu_injected - mu_true):.4f} Å")

# Create mock spectrum with the shift + realistic noise from the stack
flux_mock = gaussian_model(wave, amp_true, mu_injected, sigma_true, cont_true)
np.random.seed(42)
flux_noisy = flux_mock + np.random.normal(0, err_real, len(wave))

# 4. RECOVERY: Fit the mock to see if we recover the injected value
def chi2(params, w, f, e):
    amp, mu, sigma, cont = params
    model = gaussian_model(w, amp, mu, sigma, cont)
    return np.sum(((f - model) / e)**2)

result = minimize(chi2, popt, args=(wave[mask_line], flux_noisy[mask_line], err_real[mask_line]), 
                  method='L-BFGS-B', bounds=[(0, None), (4980, 5030), (1, 10), (None, None)])

amp_fit, mu_fit, sigma_fit, cont_fit = result.x
delta_alpha_recovered = -0.5 * (mu_fit - mu_true) / mu_true

# Compute uncertainty (sigma) from the fit covariance matrix
sigma_mu = np.sqrt(pcov[1, 1]) # Uncertainty in mu
sigma_delta_alpha = 0.5 * sigma_mu / mu_true

print(f"🎯 Recovered: Center = {mu_fit:.4f} ± {sigma_mu:.4f} Å")
print(f"📉 Recovered Δα/α = {delta_alpha_recovered:.1e} ± {sigma_delta_alpha:.1e}")

# 5. Visualization (Panel 2 of the Killer Figure)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7), gridspec_kw={'height_ratios': [2, 1]})

ax1.plot(wave[mask_line], flux_real[mask_line], color='black', alpha=0.5, label='Real Stack')
ax1.plot(wave[mask_line], flux_noisy[mask_line], color='dodgerblue', label=f'Mock (injected Δα/α = {delta_alpha_injected:.1e})')
ax1.plot(wave[mask_line], gaussian_model(wave[mask_line], amp_fit, mu_fit, sigma_fit, cont_fit), 
         color='red', linestyle='--', label=f'Fit (recovered Δα/α = {delta_alpha_recovered:.1e})')
ax1.axvline(mu_true, color='black', linestyle=':', label='True Center')
ax1.set_ylabel('Flux (W/m³)')
ax1.legend(fontsize=9)
ax1.set_title(f'Sensitivity Forecast: JADES DR3 Stack (N=330, S/N≈30)')

residuals = (flux_noisy[mask_line] - gaussian_model(wave[mask_line], amp_fit, mu_fit, sigma_fit, cont_fit)) / err_real[mask_line]
ax2.plot(wave[mask_line], residuals, color='gray', marker='.', linestyle='none')
ax2.axhline(0, color='red', linestyle='--')
ax2.set_xlabel('Rest-Frame Wavelength (Å)')
ax2.set_ylabel('Residuals ($\sigma$)')
ax2.set_ylim(-3, 3)

plt.tight_layout()
out_forecast = os.path.join(BASE, 'analysis', 'forecast_sensitivity.png')
plt.savefig(out_forecast, dpi=300)
print(f"🖼️ Forecast plot saved to: {out_forecast}")
EOF