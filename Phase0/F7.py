python3 << 'EOF'
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d
import os, re, pandas as pd
from astropy.table import Table

BASE = os.path.expanduser('~/lrd_workspace')
SPECTRA_DIR = os.path.join(BASE, 'spectra')
CATALOGS_DIR = os.path.expanduser('~/catalogs_jades')
ANALYSIS_DIR = os.path.join(BASE, 'analysis')

with open(os.path.join(ANALYSIS_DIR, 'tier_A_files.txt')) as f:
    tier_a_files = [l.strip() for l in f if l.strip()]
print(f"📂 Tier A: {len(tier_a_files)} spectra")

cat_s = Table.read(os.path.join(CATALOGS_DIR, 'hlsp_jades_jwst_nirspec_goods-s_gratings-line-fluxes_v1.1_catalog.fits'), hdu=1).to_pandas()
cat_n = Table.read(os.path.join(CATALOGS_DIR, 'hlsp_jades_jwst_nirspec_goods-n_gratings-line-fluxes_v1.1_catalog.fits'), hdu=1).to_pandas()
cats = pd.concat([cat_s, cat_n], ignore_index=True)
cats['NIRSpec_ID_str'] = cats['NIRSpec_ID'].astype(str).str.zfill(8)
z_map = {}
for _, row in cats.iterrows():
    z = row.get('z_Spec', np.nan)
    if pd.notna(z) and z > 0:
        z_map[row['NIRSpec_ID_str']] = float(z)
print(f"🔑 IDs with z: {len(z_map)}")

def extract_id(fname):
    m = re.search(r'-(\d{8})_', fname)
    return m.group(1) if m else None

wave_rest_grid = np.linspace(4700, 5200, 1500)
TARGET = 5008.24  # [OIII] 5007

sum_flux_w = np.zeros_like(wave_rest_grid)
sum_inv_var = np.zeros_like(wave_rest_grid)
n_ok = n_no_z = n_fail = n_shift = 0

for i, fname in enumerate(tier_a_files):
    if i % 50 == 0:
        print(f"  {i}/{len(tier_a_files)}...")
    spec_id = extract_id(fname)
    if spec_id is None or spec_id not in z_map:
        n_no_z += 1; continue
    z_spec = z_map[spec_id]
    try:
        with fits.open(os.path.join(SPECTRA_DIR, fname)) as h:
            wave_obs = h[1].data['WAVELENGTH'] * 1e10
            flux = h[1].data['FLUX']
            err  = h[1].data['FLUX_ERR']
        
        # FIX 1: mask BEFORE interpolating
        m = np.isfinite(wave_obs) & np.isfinite(flux) & np.isfinite(err) & (err > 0)
        wave_obs, flux, err = wave_obs[m], flux[m], err[m]
        if len(wave_obs) < 100:
            n_fail += 1; continue
        
        wave_rest = wave_obs / (1.0 + z_spec)
        
        # FIX 2: safe fill_value
        fg = interp1d(wave_rest, flux, bounds_error=False, fill_value=0.0)(wave_rest_grid)
        eg = interp1d(wave_rest, err,  bounds_error=False, fill_value=np.inf)(wave_rest_grid)
        
        # FIX 3: peak detection ONLY if significant
        win = (wave_rest_grid > 4990) & (wave_rest_grid < 5025)
        fg_smooth = gaussian_filter1d(fg, sigma=2)
        if win.sum() > 5:
            peak_idx = np.argmax(np.where(win, fg_smooth, -np.inf))
            noise = np.median(np.abs(fg[~win]))
            if fg_smooth[peak_idx] > 3 * noise and fg_smooth[peak_idx] > 0:
                shift = TARGET - wave_rest_grid[peak_idx]
                # Apply shift ONLY if reasonable (|shift| < 30 Å)
                if abs(shift) < 30:
                    fg = interp1d(wave_rest + shift, flux, bounds_error=False, fill_value=0.0)(wave_rest_grid)
                    eg = interp1d(wave_rest + shift, err,  bounds_error=False, fill_value=np.inf)(wave_rest_grid)
                    n_shift += 1
        
        # FIX 4: has_data = only where NOT inf
        has_data = np.isfinite(fg) & np.isfinite(eg) & (eg > 0) & (eg != np.inf)
        w = np.zeros_like(fg)
        w[has_data] = 1.0 / (eg[has_data]**2)
        sum_flux_w += np.where(has_data, fg * w, 0.0)
        sum_inv_var += w
        n_ok += 1
    except Exception:
        n_fail += 1

print(f"\n✅ Stacked: {n_ok}  (with shift: {n_shift})")
print(f"⚠️  No z: {n_no_z}")
print(f"❌ Failed: {n_fail}")

final_flux = np.divide(sum_flux_w, sum_inv_var, out=np.zeros_like(sum_flux_w), where=sum_inv_var>0)
final_err  = np.divide(1.0, np.sqrt(sum_inv_var), out=np.full_like(sum_inv_var, np.inf), where=sum_inv_var>0)

# Save
cols = [
    fits.Column(name='WAVE_REST', format='D', array=wave_rest_grid),
    fits.Column(name='FLUX',      format='D', array=final_flux),
    fits.Column(name='FLUX_ERR',  format='D', array=final_err),
]
fits.BinTableHDU.from_columns(cols).writeto(os.path.join(ANALYSIS_DIR, 'stacked_tierA_aligned.fits'), overwrite=True)

# Plot
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(wave_rest_grid, final_flux, 'k-', lw=1.2, label=f'Stack N={n_ok}')
ax.fill_between(wave_rest_grid, final_flux-final_err, final_flux+final_err, color='dodgerblue', alpha=0.3, label=r'$1\sigma$')
for name, lam in [('Hβ',4862.68),('[OIII]4959',4960.30),('[OIII]5007',5008.24)]:
    ax.axvline(lam, color='red', ls='--', alpha=0.7, lw=1.5)
    idx = np.argmin(np.abs(wave_rest_grid - lam))
    ax.text(lam+8, final_flux[idx], name, color='darkred', rotation=90, va='center', fontweight='bold')
ax.set_title(f'JADES DR3 Stacked (Aligned, N={n_ok})', fontsize=13)
ax.set_xlabel('Rest-Frame Wavelength (Å)')
ax.set_ylabel('Flux (W/m³)')
ax.legend(loc='upper left'); ax.grid(True, alpha=0.3); ax.set_xlim(4800, 5100)
plt.savefig(os.path.join(ANALYSIS_DIR, 'stacked_lrd_tierA_aligned.png'), dpi=300, bbox_inches='tight')
print(f"🖼️  PNG: {os.path.join(ANALYSIS_DIR, 'stacked_lrd_tierA_aligned.png')}")
print(f"💾 FITS: {os.path.join(ANALYSIS_DIR, 'stacked_tierA_aligned.fits')}")
EOF