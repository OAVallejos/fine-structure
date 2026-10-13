cd ~/lrd_workspace
python3 << 'EOF'
import numpy as np
import pandas as pd
from astropy.table import Table
import os
import re

# --- 1. Load gratings catalog with z_Spec ---
cat_s = Table.read(os.path.expanduser('~/catalogs_jades/hlsp_jades_jwst_nirspec_goods-s_gratings-line-fluxes_v1.1_catalog.fits'), hdu=1).to_pandas()
cat_n = Table.read(os.path.expanduser('~/catalogs_jades/hlsp_jades_jwst_nirspec_goods-n_gratings-line-fluxes_v1.1_catalog.fits'), hdu=1).to_pandas()
cat_s['Field'] = 'GOODS-S'
cat_n['Field'] = 'GOODS-N'
grat = pd.concat([cat_s, cat_n], ignore_index=True)

# Filter z >= 4
grat_z4 = grat[grat['z_Spec'] >= 4.0].copy()
print(f"Spectra with z >= 4 in catalog: {len(grat_z4)}")

# --- 2. Load photometry (both fields) ---
def load_phot(cat_path):
    size = Table.read(cat_path, hdu=3).to_pandas()
    kron = Table.read(cat_path, hdu=8).to_pandas()
    return size[['ID', 'RA', 'DEC', 'F444W_RHALF']].merge(
        kron[['ID', 'F444W_KRON', 'F277W_KRON']], on='ID', how='inner'
    )

df_n = load_phot(os.path.expanduser('~/catalogs_jades/hlsp_jades_jwst_nircam_goods-n_photometry_v1.0_catalog.fits'))
df_s = load_phot(os.path.expanduser('~/catalogs_jades/hlsp_jades_jwst_nircam_goods-s-deep_photometry_v2.0_catalog.fits'))
df_phot = pd.concat([df_n, df_s], ignore_index=True)

# Convert fluxes to magnitudes
def flux_to_mag(f):
    return -2.5 * np.log10(f / 3631e9) if f > 0 else np.nan
df_phot['mag_F444W'] = df_phot['F444W_KRON'].apply(flux_to_mag)
df_phot['mag_F277W'] = df_phot['F277W_KRON'].apply(flux_to_mag)
df_phot['color_F444W_F277W'] = df_phot['mag_F444W'] - df_phot['mag_F277W']

# --- 3. Cross-match by NIRCam_ID ---
grat_z4['NIRSpec_ID_str'] = grat_z4['NIRSpec_ID'].astype(str).str.zfill(8)
grat_z4['NIRCam_ID_int'] = pd.to_numeric(grat_z4['NIRCam_ID'], errors='coerce')

df_merge = grat_z4.merge(
    df_phot[['ID', 'mag_F444W', 'mag_F277W', 'color_F444W_F277W', 'F444W_RHALF']],
    left_on='NIRCam_ID_int', right_on='ID', how='left'
)

# --- 4. Apply LRD cuts ---
mask_valid = df_merge['mag_F444W'].notna()
mask_color = df_merge['color_F444W_F277W'] > 0.5
mask_mag = df_merge['mag_F444W'].between(20, 29)
mask_compact = df_merge['F444W_RHALF'] < 0.15  # arcsec (compact)

mask_lrd = mask_valid & mask_color & mask_mag & mask_compact

print(f"\n{'='*60}")
print("LRD CANDIDATES z >= 4")
print('='*60)
print(f"  Total z >= 4:              {len(df_merge)}")
print(f"  With photometric match:    {mask_valid.sum()}")
print(f"  Color > 0.5:               {mask_color.sum()}")
print(f"  Magnitude 20-29:           {mask_mag.sum()}")
print(f"  Compact (RHALF < 0.15\"):  {mask_compact.sum()}")
print(f"  🎯 LRD candidates:         {mask_lrd.sum()}")

if mask_lrd.sum() > 0:
    print(f"\n  Breakdown by field:")
    print(df_merge[mask_lrd]['Field'].value_counts().to_string())
    
    print(f"\n  z distribution:")
    print(df_merge[mask_lrd]['z_Spec'].describe())

# Save
df_merge['is_LRD_candidate'] = mask_lrd
out_csv = os.path.expanduser('~/lrd_workspace/catalogs/lrd_candidates_z4.csv')
df_merge.to_csv(out_csv, index=False)
print(f"\n💾 {out_csv}")
EOF