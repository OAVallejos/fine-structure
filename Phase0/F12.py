cd ~/lrd_workspace/analysis
python3 << 'EOF'
import numpy as np
from astropy.io import fits
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter1d
import os
import pandas as pd

ANALYSIS = os.path.expanduser('~/lrd_workspace/analysis')

# --- 1. Load the real stack ---
with fits.open(os.path.join(ANALYSIS, 'stacked_tierA_final.fits')) as h:
    t = h[1].data
    wave, flux, err = t['WAVE_REST'], t['FLUX'], t['FLUX_ERR']

# Normalize for numerical stability
SCALE = 1e-20
flux_n = flux / SCALE
err_n  = err  / SCALE

mask = np.isfinite(flux_n) & np.isfinite(err_n) & (err_n > 0)
wave, flux_n, err_n = wave[mask], flux_n[mask], err_n[mask]

fit_mask = (wave > 4830) & (wave < 5040)
w_fit, f_fit, e_fit = wave[fit_mask], flux_n[fit_mask], err_n[fit_mask]

LAM_HB, LAM_4959, LAM_5007 = 4862.68, 4960.30, 5008.24

# --- 2. FIT model (the one always used) ---
def model_fit(w, a0, a1, A_hb, mu_hb, sig_hb, A_4959, mu_4959, sig_4959, A_5007, mu_5007, sig_5007):
    cont = a0 + a1 * (w - 4900)
    return (cont
            + A_hb   * np.exp(-0.5 * ((w - mu_hb)   / sig_hb)**2)
            + A_4959 * np.exp(-0.5 * ((w - mu_4959) / sig_4959)**2)
            + A_5007 * np.exp(-0.5 * ((w - mu_5007) / sig_5007)**2))

# --- 3. Fit of the REAL stack as reference ---
p0 = [0.06, 0.0, 1.0, LAM_HB, 2.5, 1.6, LAM_4959, 2.7, 3.7, LAM_5007, 2.7]
popt_real, pcov_real = curve_fit(model_fit, w_fit, f_fit, p0=p0, sigma=e_fit,
                                 maxfev=50000, absolute_sigma=True)
perr_real = np.sqrt(np.diag(pcov_real))

print("="*70)
print("FIT OF THE REAL STACK")
print("="*70)
for n, v, e in zip(['a0','a1','A_Hb','mu_Hb','sig_Hb','A_4959','mu_4959','sig_4959','A_5007','mu_5007','sig_5007'],
                    popt_real, perr_real):
    print(f"  {n:10s} = {v:10.4f} ± {e:.4f}")

# --- 4. INJECTION MODEL (different from the fit model) ---
# Key differences:
#   a) Lines are convolved with an additional Gaussian LSF (sigma_lsf)
#   b) There is a z_Spec error that shifts the ENTIRE spectrum
#   c) There is a continuum "curvature" term (a2 * (w-4900)^2) that the fit does not model

def model_inject(w, a0, a1, a2, A_hb, mu_hb, sig_hb,
                 A_4959, mu_4959, sig_4959, A_5007, mu_5007, sig_5007,
                 sigma_lsf):
    """Injection model with extra LSF and continuum curvature."""
    cont = a0 + a1 * (w - 4900) + a2 * (w - 4900)**2
    g_hb   = A_hb   * np.exp(-0.5 * ((w - mu_hb)   / np.sqrt(sig_hb**2   + sigma_lsf**2))**2)
    g_4959 = A_4959 * np.exp(-0.5 * ((w - mu_4959) / np.sqrt(sig_4959**2 + sigma_lsf**2))**2)
    g_5007 = A_5007 * np.exp(-0.5 * ((w - mu_5007) / np.sqrt(sig_5007**2 + sigma_lsf**2))**2)
    return cont + g_hb + g_4959 + g_5007

# Injection model parameters
p_inj_base = [
    popt_real[0],   # a0
    popt_real[1],   # a1
    0.0,            # a2 (curvature, the fit does NOT model it)
    popt_real[2],   # A_hb
    popt_real[3],   # mu_hb
    popt_real[4],   # sig_hb
    popt_real[5],   # A_4959
    popt_real[6],   # mu_4959
    popt_real[7],   # sig_4959
    popt_real[8],   # A_5007
    popt_real[9],   # mu_5007
    popt_real[10],  # sig_5007
    1.0,            # sigma_lsf (Å) — extra convolution
]

# --- 5. Realistic errors ---
Z_SPEC_ERR = 0.001       # typical JADES z_Spec error
SIGMA_LSF  = 1.0         # Å of extra LSF (unmodeled resolution)
CURVATURE  = 0.0001      # continuum curvature (a2)

# --- 6. CONSERVATIVE injection-recovery loop ---
injected_values = [1e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3]
n_trials = 20  # multiple noise realizations

results = []

print("\n" + "="*70)
print("CONSERVATIVE FORECAST (systematics included)")
print("="*70)
print(f"  z_Spec error: ±{Z_SPEC_ERR}")
print(f"  Extra LSF:    {SIGMA_LSF} Å")
print(f"  Curvature:    {CURVATURE}")
print(f"  Trials:       {n_trials}")
print()

for da_inj in injected_values:
    recovered = []

    for trial in range(n_trials):
        np.random.seed(trial)

        # a) z_Spec error → shifts ALL lines
        z_err = np.random.normal(0, Z_SPEC_ERR)
        shift_frac = z_err / (1 + 5.5)  # z~5.5 typical
        delta_lambda_z = shift_frac  # relative shift in λ

        # b) Physical shift from Δα/α
        da_shift_frac = -2.0 * da_inj

        # c) Total shift
        total_shift_frac = delta_lambda_z + da_shift_frac

        mu_hb_inj   = popt_real[3] * (1 + total_shift_frac)
        mu_4959_inj = popt_real[6] * (1 + total_shift_frac)
        mu_5007_inj = popt_real[9] * (1 + total_shift_frac)

        # d) Build mock injection (with LSF and curvature)
        p_inj = list(p_inj_base)
        p_inj[4]  = mu_hb_inj
        p_inj[7]  = mu_4959_inj
        p_inj[10] = mu_5007_inj
        p_inj[2]  = CURVATURE
        p_inj[12] = SIGMA_LSF

        f_mock = model_inject(w_fit, *p_inj)

        # e) Add realistic noise
        f_mock_noisy = f_mock + np.random.normal(0, e_fit, len(w_fit))

        # f) FIT with the simple model (no LSF, no curvature)
        p0_mock = [popt_real[0], popt_real[1],
                   popt_real[2], mu_hb_inj,   popt_real[4],
                   popt_real[5], mu_4959_inj, popt_real[7],
                   popt_real[8], mu_5007_inj, popt_real[10]]

        try:
            popt_mock, pcov_mock = curve_fit(model_fit, w_fit, f_mock_noisy, p0=p0_mock,
                                              sigma=e_fit, maxfev=50000, absolute_sigma=True)
            perr_mock = np.sqrt(np.diag(pcov_mock))

            # g) Recover Δα/α from each line
            da_hb   = -0.5 * (popt_mock[3] - popt_real[3]) / popt_real[3]
            da_4959 = -0.5 * (popt_mock[6] - popt_real[6]) / popt_real[6]
            da_5007 = -0.5 * (popt_mock[9] - popt_real[9]) / popt_real[9]

            # h) Weighted average
            errs = np.array([0.5*perr_mock[3]/popt_real[3],
                             0.5*perr_mock[6]/popt_real[6],
                             0.5*perr_mock[9]/popt_real[9]])
            w = 1/errs**2
            da_mean = np.sum(w * np.array([da_hb, da_4959, da_5007])) / np.sum(w)
            recovered.append(da_mean)
        except Exception:
            continue

    if len(recovered) > 0:
        recovered = np.array(recovered)
        bias = np.mean(recovered) - da_inj
        scatter = np.std(recovered)
        results.append({
            'da_injected': da_inj,
            'da_mean': np.mean(recovered),
            'bias': bias,
            'scatter': scatter,
            'n_ok': len(recovered),
        })
        scatter_str = f"{scatter:.3e}"
        print(f"Δα/α = {da_inj:.1e}: recovered = {np.mean(recovered):+.3e} ± {scatter_str}  "
              f"(bias = {bias:+.2e}, {abs(bias)/scatter if scatter>0 else 0:.1f}σ)")

# --- 7. Save CSV ---
df = pd.DataFrame(results)
out_csv = os.path.join(ANALYSIS, 'forecast_conservative.csv')
df.to_csv(out_csv, index=False)
print(f"\n💾 {out_csv}")

# --- 8. Real sensitivity ---
print("\n" + "="*70)
print("REAL SENSITIVITY (conservative)")
print("="*70)
print(f"  Statistical σ (ideal pipeline): ~2.8 × 10⁻⁷")
print(f"  Conservative σ (with systematics): ~{np.mean([r['scatter'] for r in results]):.1e}")
print()
print("Interpretation:")
print("  - Statistical error is the theoretical floor of the pipeline")
print("  - Conservative error includes: z_Spec, LSF, curvature, noise")
print("  - The REAL error lies between both, closer to the conservative one")
EOF