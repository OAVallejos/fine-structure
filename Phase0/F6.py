cd ~/lrd_workspace/analysis
python3 << 'EOF'
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import interp1d
import os
import re
import pandas as pd
from astropy.table import Table

BASE = os.path.expanduser('~/lrd_workspace')
SPECTRA_DIR = os.path.join(BASE, 'spectra')
CATALOGS_DIR = os.path.expanduser('~/catalogs_jades')
ANALYSIS_DIR = os.path.join(BASE, 'analysis')

# 1. Load Tier A
with open(os.path.join(ANALYSIS_DIR, 'tier_A_files.txt')) as f:
    tier_a_files = [l.strip() for l in f if l.strip()]
print(f"📂 Tier A: {len(tier_a_files)} spectra")

# 2. Load catalogs and build z index by ID
cat_s = Table.read(os.path.join(CATALOGS_DIR, 'hlsp_jades_jwst_nirspec_goods-s_gratings-line-fluxes_v1.1_catalog.fits'), hdu=1).to_pandas()
cat_n = Table.read(os.path.join(CATALOGS_DIR, 'hlsp_jades_jwst_nirspec_goods-n_gratings-line-fluxes_v1.1_catalog.fits'), hdu=1).to_pandas()
cats = pd.concat([cat_s, cat_n], ignore_index=True)
cats['NIRSpec_ID_str'] = cats['NIRSpec_ID'].astype(str).str.zfill(8)

z_map = {}
for _, row in cats.iterrows():
    z = row.get('z_Spec', np.nan)
    if pd.notna(z) and z > 0:
        z_map[row['NIRSpec_ID_str']] = float(z)
print(f"🔑 IDs with valid z_spec: {len(z_map)}")

def extract_id(fname):
    m = re.search(r'-(\d{8})_', fname)
    return m.group(1) if m else None

wave_rest_grid = np.linspace(4750, 5150, 1000)
sum_flux_w = np.zeros_like(wave_rest_grid)
sum_inv_var = np.zeros_like(wave_rest_grid)
n_ok = n_no_z = n_fail = 0

for i, fname in enumerate(tier_a_files):
    if i % 50 == 0:
        print(f"  {i}/{len(tier_a_files)}...")
    spec_id = extract_id(fname)
    if spec_id is None or spec_id not in z_map:
        n_no_z += 1
        continue
    z_spec = z_map[spec_id]
    try:
        with fits.open(os.path.join(SPECTRA_DIR, fname)) as h:
            wave_obs = h[1].data['WAVELENGTH'] * 1e10
            flux = h[1].data['FLUX']
            err = h[1].data['FLUX_ERR']
        mask = np.isfinite(flux) & np.isfinite(err) & (err > 0)
        if mask.sum() < 100:
            n_fail += 1
            continue
        wave_obs, flux, err = wave_obs[mask], flux[mask], err[mask]
        wave_rest = wave_obs / (1.0 + z_spec)
        fg = interp1d(wave_rest, flux, bounds_error=False, fill_value=np.nan)(wave_rest_grid)
        eg = interp1d(wave_rest, err,  bounds_error=False, fill_value=np.nan)(wave_rest_grid)
        valid = np.isfinite(fg) & np.isfinite(eg) & (eg > 0)
        w = np.zeros_like(fg); w[valid] = 1.0/eg[valid]**2
        sum_flux_w += np.where(valid, fg*w, 0.0)
        sum_inv_var += w
        n_ok += 1
    except Exception:
        n_fail += 1

print(f"\n✅ Stacked: {n_ok}")
print(f"⚠️  No z_spec: {n_no_z}")
print(f"❌ Failed: {n_fail}")

final_flux = np.divide(sum_flux_w, sum_inv_var, out=np.zeros_like(sum_flux_w), where=sum_inv_var>0)
final_err = np.divide(1.0, np.sqrt(sum_inv_var), out=np.full_like(sum_inv_var, np.inf), where=sum_inv_var>0)

# Save FITS
cols = [
    fits.Column(name='WAVE_REST', format='D', array=wave_rest_grid),
    fits.Column(name='FLUX',      format='D', array=final_flux),
    fits.Column(name='FLUX_ERR',  format='D', array=final_err),
]
fits.BinTableHDU.from_columns(cols).writeto(os.path.join(ANALYSIS_DIR, 'stacked_tierA_restframe.fits'), overwrite=True)

# Killer Figure
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(wave_rest_grid, final_flux, color='black', lw=1.2, label=f'Stack N={n_ok}')
ax.fill_between(wave_rest_grid, final_flux-final_err, final_flux+final_err, color='dodgerblue', alpha=0.3, label=r'$1\sigma$')
for name, lam in [('Hβ',4862.68),('[OIII]4959',4960.30),('[OIII]5007',5008.24)]:
    ax.axvline(lam, color='red', ls='--', alpha=0.7, lw=1.5)
    idx = np.argmin(np.abs(wave_rest_grid - lam))
    ax.text(lam+8, final_flux[idx], name, color='darkred', rotation=90, va='center', fontweight='bold')
ax.set_title(f'JADES DR3 Stacked NIRSpec Rest-Frame (Tier A, N={n_ok})', fontsize=13)
ax.set_xlabel('Rest-Frame Wavelength (Å)')
ax.set_ylabel('Flux (W/m³)')
ax.legend(loc='upper left'); ax.grid(True, alpha=0.3); ax.set_xlim(4800, 5100)
plt.savefig(os.path.join(ANALYSIS_DIR, 'stacked_lrd_tierA.png'), dpi=300, bbox_inches='tight')
print(f"\n🖼️  PNG: {os.path.join(ANALYSIS_DIR, 'stacked_lrd_tierA.png')}")
print(f"💾 FITS: {os.path.join(ANALYSIS_DIR, 'stacked_tierA_restframe.fits')}")
EOF