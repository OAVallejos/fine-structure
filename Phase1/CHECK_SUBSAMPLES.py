# CHECK_SUBSAMPLES.py
import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.optimize import curve_fit
import warnings
warnings.filterwarnings('ignore')

GRID = np.linspace(4800, 6700, 1900)
LAB = {'Hb': 4862.68, 'OIII4959': 4960.30, 'OIII5007': 5008.24, 'Ha': 6564.63}

tier_A = pd.read_csv('lrd_workspace/analysis/tier_A_ha_valid.csv')
tier_A['ID'] = tier_A['ID'].astype(str).str.zfill(8)

cat = pd.read_csv('lrd_workspace/analysis/master_features_v4.csv')
cat['NIRSpec_ID_str'] = cat['NIRSpec_ID'].astype(str).str.zfill(8)
z_map = dict(zip(cat['NIRSpec_ID_str'], cat['z_Spec']))

qm = pd.read_csv('lrd_workspace/analysis/quality_metrics.csv')
snr_map = dict(zip(qm['file'], qm['snr_peak']))
tier_A['snr'] = tier_A['file'].map(snr_map)

def load_spectrum(f, z):
    try:
        with fits.open(f'lrd_workspace/spectra/{f}') as hdul:
            d = hdul[1].data
            cn = hdul[1].columns.names
            wc = next((c for c in cn if 'wave' in c.lower()), cn[0])
            fc = next((c for c in cn if 'flux' in c.lower() and 'err' not in c.lower()), cn[1])
            ec = next((c for c in cn if 'err' in c.lower()), cn[2])
            w = d[wc].copy(); fl = d[fc].copy(); er = d[ec].copy()
        if np.nanmax(w) < 100: w = w * 1e4
        w_rest = w / (1 + z)
        good = np.isfinite(fl) & np.isfinite(er) & (er > 0)
        if good.sum() < 50: return None
        return w_rest[good], fl[good], er[good]
    except Exception:
        return None

def build_stack(df):
    num = np.zeros_like(GRID); den = np.zeros_like(GRID)
    n_ok = 0
    for _, row in df.iterrows():
        z = z_map.get(row['ID'])
        if z is None or z <= 0: continue
        spec = load_spectrum(row['file'], z)
        if spec is None: continue
        w_rest, fl, er = spec
        f_int = np.interp(GRID, w_rest, fl, left=0, right=0)
        e_int = np.interp(GRID, w_rest, er, left=np.inf, right=np.inf)
        wgt = 1.0 / e_int**2
        wgt[~np.isfinite(wgt)] = 0
        num += f_int * wgt; den += wgt
        n_ok += 1
    sf = np.where(den > 0, num / den, 0)
    se = np.where(den > 0, 1.0 / np.sqrt(den), np.inf)
    return sf, se, n_ok

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

def fit_alpha(df, label):
    if len(df) < 20:
        print(f"  {label:22s}  N={len(df):4d}  — too few")
        return None
    sf, se, n_ok = build_stack(df)   # ← FIX: 3 values
    if not np.any(np.isfinite(se)):
        print(f"  {label:22s}  N={len(df):4d}  — empty stack")
        return None
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
        d_hb = popt[3] - LAB['Hb']
        d_5007 = popt[9] - LAB['OIII5007']
        e_hb = perr[3]
        e_5007 = perr[9]
        diff = d_hb - d_5007
        e_diff = np.sqrt(e_hb**2 + e_5007**2)
        dalpha = -diff / (2 * LAB['Hb'])
        e_dalpha = e_diff / (2 * LAB['Hb'])
        sig = diff / e_diff if e_diff > 0 else 0
        print(f"  {label:22s}  N={len(df):4d}  "
              f"Δα/α = {dalpha:+.6e} ± {e_dalpha:.6e}  "
              f"({dalpha*1e6:+8.4f} ± {e_dalpha*1e6:6.4f} ppm)  "
              f"({sig:+.3f}σ)")
        return dalpha, e_dalpha
    except Exception as ex:
        print(f"  {label:22s}  N={len(df):4d}  — fit failed: {ex}")
        return None

print("="*90)
print("α TEST ON SUBSAMPLES")
print("="*90)

print("\n[1] By field:")
for field in ['goods-s', 'goods-n']:
    sub = tier_A[tier_A['field'] == field]
    fit_alpha(sub, field.upper())

print("\n[2] By z bins:")
for zlo, zhi in [(5.0, 5.5), (5.5, 6.0), (6.0, 6.5), (6.5, 8.0)]:
    sub = tier_A[(tier_A['z'] >= zlo) & (tier_A['z'] < zhi)]
    fit_alpha(sub, f"z∈[{zlo},{zhi}]")

print("\n[3] By S/N quartiles:")
qs = tier_A['snr'].quantile([0, 0.25, 0.5, 0.75, 1.0]).values
for i in range(4):
    sub = tier_A[(tier_A['snr'] >= qs[i]) & (tier_A['snr'] < qs[i+1])]
    fit_alpha(sub, f"S/N Q{i+1} [{qs[i]:.1f},{qs[i+1]:.1f}]")

print("\n[4] Without outliers (S/N > 100):")
sub = tier_A[tier_A['snr'] <= 100]
fit_alpha(sub, "S/N ≤ 100")

print("\n[5] Full sample:")
fit_alpha(tier_A, "ALL")