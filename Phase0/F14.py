cd ~/lrd_workspace/analysis
python3 << 'EOF'
import numpy as np
from astropy.io import fits
from scipy.optimize import curve_fit
import os

ANALYSIS = os.path.expanduser('~/lrd_workspace/analysis')

# --- Load stack ---
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

def chi2_of(popt, model, w, f, e):
    resid = (f - model(w, *popt)) / e
    return np.sum(resid**2), len(w) - len(popt)

print("="*70)
print("TEST A: FREE MODEL (11 params)")
print("="*70)

def model_free(w, a0, a1, A_hb, mu_hb, sig_hb, A_4959, mu_4959, sig_4959, A_5007, mu_5007, sig_5007):
    return (a0 + a1*(w-4900)
            + A_hb   * np.exp(-0.5*((w-mu_hb)/sig_hb)**2)
            + A_4959 * np.exp(-0.5*((w-mu_4959)/sig_4959)**2)
            + A_5007 * np.exp(-0.5*((w-mu_5007)/sig_5007)**2))

p0 = [0.06, 0.0, 1.0, LAM_HB, 2.5, 1.6, LAM_4959, 2.7, 3.7, LAM_5007, 2.7]
popt_free, pcov_free = curve_fit(model_free, w_fit, f_fit, p0=p0, sigma=e_fit,
                                  maxfev=50000, absolute_sigma=True)
perr_free = np.sqrt(np.diag(pcov_free))
chi2_free, dof_free = chi2_of(popt_free, model_free, w_fit, f_fit, e_fit)
print(f"  χ² = {chi2_free:.2f}  dof = {dof_free}  χ²_red = {chi2_free/dof_free:.2f}")
print()
print(f"  mu_Hb   = {popt_free[3]:+.4f} ± {perr_free[3]:.4f}")
print(f"  mu_4959 = {popt_free[6]:+.4f} ± {perr_free[6]:.4f}")
print(f"  mu_5007 = {popt_free[9]:+.4f} ± {perr_free[9]:.4f}")
print(f"  A_5007/A_4959 = {popt_free[8]/popt_free[5]:.3f}")
print(f"  A_5007/A_Hb   = {popt_free[8]/popt_free[2]:.3f}")

print("\n" + "="*70)
print("TEST B: COUPLED WIDTHS (common sig)")
print("="*70)

def model_same_width(w, a0, a1, A_hb, mu_hb, A_4959, mu_4959, A_5007, mu_5007, sig):
    return (a0 + a1*(w-4900)
            + A_hb   * np.exp(-0.5*((w-mu_hb)/sig)**2)
            + A_4959 * np.exp(-0.5*((w-mu_4959)/sig)**2)
            + A_5007 * np.exp(-0.5*((w-mu_5007)/sig)**2))

p0_b = [0.06, 0.0, 1.0, LAM_HB, 1.6, LAM_4959, 3.7, LAM_5007, 2.7]
popt_b, pcov_b = curve_fit(model_same_width, w_fit, f_fit, p0=p0_b, sigma=e_fit,
                           maxfev=50000, absolute_sigma=True)
chi2_b, dof_b = chi2_of(popt_b, model_same_width, w_fit, f_fit, e_fit)
print(f"  χ² = {chi2_b:.2f}  dof = {dof_b}  χ²_red = {chi2_b/dof_b:.2f}")
print(f"  Δχ² vs free = {chi2_b - chi2_free:+.2f}  (2 fewer params)")

print("\n" + "="*70)
print("TEST C: COUPLED SHIFTS (common shift to all 3 lines)")
print("="*70)

def model_same_shift(w, a0, a1, A_hb, sig_hb, A_4959, sig_4959, A_5007, sig_5007, shift):
    return (a0 + a1*(w-4900)
            + A_hb   * np.exp(-0.5*((w-(LAM_HB+shift))/sig_hb)**2)
            + A_4959 * np.exp(-0.5*((w-(LAM_4959+shift))/sig_4959)**2)
            + A_5007 * np.exp(-0.5*((w-(LAM_5007+shift))/sig_5007)**2))

p0_c = [0.06, 0.0, 1.0, 2.5, 1.6, 2.7, 3.7, 2.7, -0.08]
popt_c, pcov_c = curve_fit(model_same_shift, w_fit, f_fit, p0=p0_c, sigma=e_fit,
                           maxfev=50000, absolute_sigma=True)
chi2_c, dof_c = chi2_of(popt_c, model_same_shift, w_fit, f_fit, e_fit)
perr_c = np.sqrt(np.diag(pcov_c))
print(f"  χ² = {chi2_c:.2f}  dof = {dof_c}  χ²_red = {chi2_c/dof_c:.2f}")
print(f"  Δχ² vs free = {chi2_c - chi2_free:+.2f}  (2 fewer params)")
print(f"  common shift = {popt_c[8]:+.4f} ± {perr_c[8]:.4f} Å")

print("\n" + "="*70)
print("TEST D: 5007/4959 RATIO FIXED AT 2.98")
print("="*70)

def model_fixed_ratio(w, a0, a1, A_hb, mu_hb, sig_hb, A_4959, mu_4959, sig_4959, mu_5007, sig_5007):
    A_5007 = A_4959 * 2.98
    return (a0 + a1*(w-4900)
            + A_hb   * np.exp(-0.5*((w-mu_hb)/sig_hb)**2)
            + A_4959 * np.exp(-0.5*((w-mu_4959)/sig_4959)**2)
            + A_5007 * np.exp(-0.5*((w-mu_5007)/sig_5007)**2))

p0_d = [0.06, 0.0, 1.0, LAM_HB, 2.5, 1.6, LAM_4959, 2.7, LAM_5007, 2.7]
popt_d, pcov_d = curve_fit(model_fixed_ratio, w_fit, f_fit, p0=p0_d, sigma=e_fit,
                           maxfev=50000, absolute_sigma=True)
chi2_d, dof_d = chi2_of(popt_d, model_fixed_ratio, w_fit, f_fit, e_fit)
print(f"  χ² = {chi2_d:.2f}  dof = {dof_d}  χ²_red = {chi2_d/dof_d:.2f}")
print(f"  Δχ² vs free = {chi2_d - chi2_free:+.2f}  (1 fewer param)")

print("\n" + "="*70)
print("SUMMARY: WHICH MODEL WINS?")
print("="*70)
print(f"  A) Free (11p):           χ² = {chi2_free:.1f}  χ²_red = {chi2_free/dof_free:.2f}")
print(f"  B) Common widths (9p):   χ² = {chi2_b:.1f}  χ²_red = {chi2_b/dof_b:.2f}  Δχ² = {chi2_b-chi2_free:+.1f}")
print(f"  C) Common shift (9p):    χ² = {chi2_c:.1f}  χ²_red = {chi2_c/dof_c:.2f}  Δχ² = {chi2_c-chi2_free:+.1f}")
print(f"  D) Fixed ratio (10p):    χ² = {chi2_d:.1f}  χ²_red = {chi2_d/dof_d:.2f}  Δχ² = {chi2_d-chi2_free:+.1f}")

print()
print("Interpretation:")
print("  - If χ²_red > 10 in ALL models → the Gaussian model does not capture the shape")
print("  - If χ²_red < 2 in some model → that model is the correct one")
print("  - Δχ² > 25 (3σ, 2 dof) → the restricted model is worse, the restriction does NOT hold")
print("  - Δχ² < 6 (1σ, 2 dof)  → the restricted model is equally good, the restriction DOES hold")
EOF