import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.optimize import curve_fit
import warnings
warnings.filterwarnings('ignore')

N_BOOTSTRAP = 100
GRID = np.linspace(4800, 6700, 1900)
LAB = {'Hb': 4862.68, 'OIII4959': 4960.30, 'OIII5007': 5008.24, 'Ha': 6564.63}

tier_A = pd.read_csv('lrd_workspace/analysis/tier_A_ha_valid.csv')
tier_A['ID'] = tier_A['ID'].astype(str).str.zfill(8)

cat = pd.read_csv('lrd_workspace/analysis/master_features_v4.csv')
cat['NIRSpec_ID_str'] = cat['NIRSpec_ID'].astype(str).str.zfill(8)
z_map = dict(zip(cat['NIRSpec_ID_str'], cat['z_Spec']))

print(f"Tier A: {len(tier_A)} objects")
print(f"z_map: {len(z_map)} entries")
print(f"Example ID: {repr(tier_A['ID'].iloc[0])}")
print()

def load_spectrum(f, z):
    try:
        with fits.open(f'lrd_workspace/spectra/{f}') as hdul:
            d = hdul[1].data
            cn = hdul[1].columns.names
            wc = next((c for c in cn if 'wave' in c.lower()), cn[0])
            fc = next((c for c in cn if 'flux' in c.lower() and 'err' not in c.lower()), cn[1])
            ec = next((c for c in cn if 'err' in c.lower()), cn[2] if len(cn) > 2 else cn[1])
            w = d[wc].copy(); fl = d[fc].copy(); er = d[ec].copy()
        if np.nanmax(w) < 100: w = w * 1e4
        w_rest = w / (1 + z)
        good = np.isfinite(fl) & np.isfinite(er) & (er > 0)
        if good.sum() < 50: return None
        return w_rest[good], fl[good], er[good]
    except Exception:
        return None

def build_stack(sample_df):
    num = np.zeros_like(GRID); den = np.zeros_like(GRID)
    n_ok = 0; n_skip = 0
    for _, row in sample_df.iterrows():
        z = z_map.get(row['ID'])
        if z is None or z <= 0:
            n_skip += 1; continue
        spec = load_spectrum(row['file'], z)
        if spec is None:
            n_skip += 1; continue
        w_rest, fl, er = spec
        f_int = np.interp(GRID, w_rest, fl, left=0, right=0)
        e_int = np.interp(GRID, w_rest, er, left=np.inf, right=np.inf)
        wgt = 1.0 / e_int**2
        wgt[~np.isfinite(wgt)] = 0
        num += f_int * wgt; den += wgt
        n_ok += 1
    sf = np.where(den > 0, num / den, 0)
    se = np.where(den > 0, 1.0 / np.sqrt(den), np.inf)
    return sf, se, n_ok, n_skip

def model(w, a0, a1,
          A_hb, mu_hb, sig_hb,
          A_4959, mu_4959, sig_4959,
          A_5007, mu_5007, sig_5007,
          A_ha, mu_ha, sig_ha_nii, A_nii):
    cont = a0 + a1 * (w - 5700)
    D_NII = 6585.27 - 6564.63
    g_hb   = A_hb   * np.exp(-0.5 * ((w - mu_hb)   / sig_hb)**2)
    g_4959 = A_4959 * np.exp(-0.5 * ((w - mu_4959) / sig_4959)**2)
    g_5007 = A_5007 * np.exp(-0.5 * ((w - mu_5007) / sig_5007)**2)
    g_ha   = A_ha   * np.exp(-0.5 * ((w - mu_ha)   / sig_ha_nii)**2)
    g_nii  = A_nii  * np.exp(-0.5 * ((w - (mu_ha + D_NII)) / sig_ha_nii)**2)
    return cont + g_hb + g_4959 + g_5007 + g_ha + g_nii

def fit_stack(sf, se):
    cont = np.median(sf[(GRID > 5200) & (GRID < 5400)])
    def amp(lam, w=5):
        m = (GRID > lam - w) & (GRID < lam + w)
        return max(np.max(sf[m]) - cont, abs(cont) * 0.5)
    p0 = [cont, 0.0,
          amp(4862.68), 4862.68, 2.5,
          amp(4960.30), 4960.30, 2.7,
          amp(5008.24), 5008.24, 2.7,
          amp(6564.63, 10), 6564.63, 3.0, amp(6585.27, 10) * 0.3]
    try:
        popt, pcov = curve_fit(model, GRID, sf, p0=p0, sigma=se,
                               absolute_sigma=True, maxfev=200000)
        perr = np.sqrt(np.diag(pcov))
        if popt[8] / popt[5] < 1.5 or popt[8] / popt[5] > 4:
            return None
        return popt, perr
    except Exception:
        return None

print("Building original stack...")
sf, se, n_ok, n_skip = build_stack(tier_A)
print(f"Stack: {n_ok} objects stacked, {n_skip} skipped")

if n_ok == 0:
    print("ERROR: 0 objects. Fix paths/IDs first.")
    exit(1)

res = fit_stack(sf, se)
if res is None:
    print("ORIGINAL FIT FAILED"); exit(1)
popt0, perr0 = res

d_hb_0 = popt0[3] - LAB['Hb']
d_5007_0 = popt0[9] - LAB['OIII5007']
diff_0 = d_hb_0 - d_5007_0
dalpha_0 = -diff_0 / (2 * LAB['Hb'])

print(f"\nORIGINAL STACK:")
print(f"  Hβ       Δλ = {d_hb_0:+.4f} ± {perr0[3]:.4f} Å")
print(f"  OIII5007 Δλ = {d_5007_0:+.4f} ± {perr0[9]:.4f} Å")
print(f"  Difference = {diff_0:+.4f} Å")
print(f"  Δα/α = {dalpha_0:+.3e}  ({dalpha_0*1e6:+.2f} ppm)")

print(f"\nRunning bootstrap ({N_BOOTSTRAP} iterations)...")
dalpha_boot = []
for i in range(N_BOOTSTRAP):
    sample = tier_A.sample(n=len(tier_A), replace=True, random_state=i)
    sf_b, se_b, n_b, _ = build_stack(sample)
    r = fit_stack(sf_b, se_b)
    if r is None: continue
    popt, _ = r
    d_hb = popt[3] - LAB['Hb']
    d_5007 = popt[9] - LAB['OIII5007']
    diff = d_hb - d_5007
    dalpha = -diff / (2 * LAB['Hb'])
    dalpha_boot.append(dalpha)
    if (i+1) % 10 == 0:
        print(f"  {i+1}/{N_BOOTSTRAP}")

dalpha_boot = np.array(dalpha_boot)
print(f"\nBOOTSTRAP RESULTS ({len(dalpha_boot)} valid):")
print(f"  Mean    = {dalpha_boot.mean():+.3e}  ({dalpha_boot.mean()*1e6:+.2f} ppm)")
print(f"  Median  = {np.median(dalpha_boot):+.3e}")
print(f"  Std     = {dalpha_boot.std():.3e}  ({dalpha_boot.std()*1e6:.2f} ppm)")
print(f"  2σ      = [{np.percentile(dalpha_boot, 2.5):+.3e}, {np.percentile(dalpha_boot, 97.5):+.3e}]")

pd.DataFrame({'dalpha_over_alpha': dalpha_boot}).to_csv(
    'lrd_workspace/analysis/bootstrap_alpha.csv', index=False)

print("\n" + "="*70)
print("VERDICT")
print("="*70)
if abs(dalpha_boot.mean() / dalpha_boot.std()) < 2:
    print(f"  Δα/α = {dalpha_boot.mean()*1e6:+.2f} ± {dalpha_boot.std()*1e6:.2f} ppm")
    print(f"  Consistent with ZERO")
    print(f"  2σ limit: |Δα/α| < {2*dalpha_boot.std():.3e}")
else:
    print(f"  Δα/α = {dalpha_boot.mean()*1e6:+.2f} ± {dalpha_boot.std()*1e6:.2f} ppm")
    print(f"  Significant: {dalpha_boot.mean()/dalpha_boot.std():+.2f}σ")