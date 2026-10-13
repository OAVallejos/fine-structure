cd ~/lrd_workspace/analysis
python3 << 'EOF'
import numpy as np
from astropy.io import fits
from scipy.optimize import curve_fit
import os

ANALYSIS = os.path.expanduser('~/lrd_workspace/analysis')

# --- 1. Load stack ---
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

p0 = [0.06, 0.0, 1.0, LAM_HB, 2.5, 1.6, LAM_4959, 2.7, 3.7, LAM_5007, 2.7]
popt, pcov = curve_fit(model, w_fit, f_fit, p0=p0, sigma=e_fit, maxfev=50000, absolute_sigma=True)
perr = np.sqrt(np.diag(pcov))

# --- 2. Offsets and errors ---
lines = [
    ('Hβ',         popt[3], perr[3], LAM_HB),
    ('[OIII]4959', popt[6], perr[6], LAM_4959),
    ('[OIII]5007', popt[9], perr[9], LAM_5007),
]

print("="*70)
print("TEST 1: Are the offsets consistent between lines?")
print("="*70)
print(f"{'Line':12s}  {'Δλ (Å)':>12s}  {'σ_Δλ (Å)':>10s}  {'Δλ/σ':>8s}")
print("-"*50)

offsets = []
for name, mu, mu_err, lam in lines:
    dlam = mu - lam
    snr = dlam / mu_err
    offsets.append((name, dlam, mu_err))
    print(f"{name:12s}  {dlam:+12.4f}  {mu_err:10.4f}  {snr:+8.2f}")

# --- 3. Consistency test: are all offsets equal? ---
print("\n" + "="*70)
print("TEST 2: Internal consistency of the offset")
print("="*70)

dlam_arr = np.array([o[1] for o in offsets])
err_arr  = np.array([o[2] for o in offsets])

# Weighted average
w = 1/err_arr**2
dlam_mean = np.sum(w*dlam_arr)/np.sum(w)
dlam_err  = 1/np.sqrt(np.sum(w))

print(f"  Average offset: {dlam_mean:+.4f} ± {dlam_err:.4f} Å")
print()

# χ² relative to the average
chi2 = np.sum(((dlam_arr - dlam_mean)/err_arr)**2)
dof = len(dlam_arr) - 1
chi2_red = chi2/dof

print(f"  χ² = {chi2:.2f}  (dof = {dof})")
print(f"  χ²_red = {chi2_red:.2f}")
print()
if chi2_red < 3:
    print("  ✅ CONSISTENT: offsets are equal within errors")
    print("     → Interpretation: NIRSpec CALIBRATION")
elif chi2_red < 10:
    print("  ⚠️  MARGINAL: there is some scatter")
    print("     → Interpretation: mix of calibration + something else")
else:
    print("  ❌ INCONSISTENT: offsets are NOT equal")
    print("     → Interpretation: possible physical signal (unlikely) or blending")

# --- 4. Direction test: are they all negative? ---
print("\n" + "="*70)
print("TEST 3: Offset direction")
print("="*70)
signs = np.sign(dlam_arr)
n_neg = np.sum(signs < 0)
n_pos = np.sum(signs > 0)
print(f"  Lines with negative offset: {n_neg}/3")
print(f"  Lines with positive offset: {n_pos}/3")
print()
if n_neg == 3:
    print("  ✅ ALL offsets are negative")
    print("     → Consistent with systematic calibration (common shift)")
elif n_pos == 3:
    print("  ✅ ALL offsets are positive")
    print("     → Consistent with systematic calibration (common shift)")
else:
    print("  ⚠️  Offsets in opposite directions")
    print("     → Could be physics or complex blending")

# --- 5. α vs no-α test: Hβ vs [OIII] ---
print("\n" + "="*70)
print("TEST 4: Hβ (sensitive to α) vs [OIII] (NOT sensitive to α)")
print("="*70)

# Hβ offset
dlam_hb = dlam_arr[0]
err_hb  = err_arr[0]

# Average [OIII] offset
dlam_oiii = np.mean(dlam_arr[1:])
err_oiii  = np.std(dlam_arr[1:]) / np.sqrt(2)

diff = dlam_hb - dlam_oiii
err_diff = np.sqrt(err_hb**2 + err_oiii**2)
significance = diff / err_diff

print(f"  Δλ(Hβ)   = {dlam_hb:+.4f} ± {err_hb:.4f} Å")
print(f"  Δλ([OIII]) = {dlam_oiii:+.4f} ± {err_oiii:.4f} Å")
print(f"  Difference = {diff:+.4f} ± {err_diff:.4f} Å")
print(f"  Significance = {significance:+.2f}σ")
print()

if abs(significance) < 2:
    print("  ✅ Hβ and [OIII] are shifted EQUALLY")
    print("     → Interpretation: CALIBRATION (affects the entire spectrum equally)")
    print("     → NO evidence of physical signal in α")
elif abs(significance) < 5:
    print("  ⚠️  Marginal difference between Hβ and [OIII]")
    print("     → Could be calibration + something else. Requires more analysis.")
else:
    print("  🚨 Hβ and [OIII] are shifted DIFFERENTLY")
    print("     → Possible physical signal. Requires independent validation.")

# --- 6. Save results ---
import pandas as pd
df = pd.DataFrame({
    'line': [o[0] for o in offsets],
    'delta_lambda_A': dlam_arr,
    'err_A': err_arr,
    'significance': dlam_arr/err_arr,
})
out_csv = os.path.join(ANALYSIS, 'offset_test.csv')
df.to_csv(out_csv, index=False)
print(f"\n💾 {out_csv}")

# --- 7. Final summary ---
print("\n" + "="*70)
print("TEST SUMMARY")
print("="*70)
print(f"  Average offset:  {dlam_mean:+.4f} ± {dlam_err:.4f} Å")
print(f"  χ²_red:          {chi2_red:.2f}")
print(f"  Hβ vs [OIII]:    {significance:+.2f}σ")
print()
print("  Verdict:")
if chi2_red < 3 and abs(significance) < 2:
    print("  ✅ OFFSET = NIRSpec CALIBRATION")
    print("     No evidence of physical Δα/α signal.")
    print("     The offset should be corrected with cross-calibration.")
elif chi2_red < 3 and abs(significance) >= 2:
    print("  ⚠️  OFFSET PARTIALLY EXPLAINED BY CALIBRATION")
    print("     There is an additional component. Requires more analysis.")
else:
    print("  🚨 OFFSET NOT CONSISTENT WITH CALIBRATION")
    print("     There could be a physical signal. Requires urgent validation.")
EOF