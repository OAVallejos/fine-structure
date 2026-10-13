cd ~/lrd_workspace/analysis
python3 << 'EOF'
import numpy as np
import pandas as pd
from astropy.table import Table
import os
import re

# --- 1. Load photometry ---
def load_photometry(cat_path):
    size = Table.read(cat_path, hdu=3).to_pandas()
    kron = Table.read(cat_path, hdu=8).to_pandas()
    return size[['ID', 'RA', 'DEC', 'F444W_RHALF', 'F277W_RHALF']].merge(
        kron[['ID', 'F444W_KRON', 'F277W_KRON']],
        on='ID', how='inner'
    )

print("Loading catalogs...")
df_n = load_photometry(os.path.expanduser('~/catalogs_jades/hlsp_jades_jwst_nircam_goods-n_photometry_v1.0_catalog.fits'))
df_s = load_photometry(os.path.expanduser('~/catalogs_jades/hlsp_jades_jwst_nircam_goods-s-deep_photometry_v2.0_catalog.fits'))
df_phot = pd.concat([df_n, df_s], ignore_index=True)
print(f"Total sources: {len(df_phot)}")

# --- 2. Detect if they are fluxes or magnitudes (by median) ---
median_f444 = df_phot['F444W_KRON'].dropna().median()
print(f"\nMedian F444W_KRON: {median_f444:.4f}")

if 0 < median_f444 < 1000:
    print("→ They are FLUXES (nJy). Converting to magnitudes...")
    def flux_to_mag(f):
        return -2.5 * np.log10(f / 3631e9) if f > 0 else np.nan
    df_phot['mag_F444W'] = df_phot['F444W_KRON'].apply(flux_to_mag)
    df_phot['mag_F277W'] = df_phot['F277W_KRON'].apply(flux_to_mag)
else:
    print("→ They are MAGNITUDES")
    df_phot['mag_F444W'] = df_phot['F444W_KRON']
    df_phot['mag_F277W'] = df_phot['F277W_KRON']

# Verify that the magnitudes are reasonable
print(f"\nmag_F444W statistics:")
print(f"  min: {df_phot['mag_F444W'].min():.2f}")
print(f"  median: {df_phot['mag_F444W'].median():.2f}")
print(f"  max: {df_phot['mag_F444W'].max():.2f}")

df_phot['color_F444W_F277W'] = df_phot['mag_F444W'] - df_phot['mag_F277W']
print(f"\ncolor F444W-F277W statistics:")
print(f"  median: {df_phot['color_F444W_F277W'].median():.2f}")
print(f"  90th percentile: {df_phot['color_F444W_F277W'].quantile(0.9):.2f}")
print(f"  99th percentile: {df_phot['color_F444W_F277W'].quantile(0.99):.2f}")

# --- 3. Load gratings ---
cat_s_grat = Table.read(os.path.expanduser('~/catalogs_jades/hlsp_jades_jwst_nirspec_goods-s_gratings-line-fluxes_v1.1_catalog.fits'), hdu=1).to_pandas()
cat_n_grat = Table.read(os.path.expanduser('~/catalogs_jades/hlsp_jades_jwst_nirspec_goods-n_gratings-line-fluxes_v1.1_catalog.fits'), hdu=1).to_pandas()
cat_s_grat['Field'] = 'GOODS-S'
cat_n_grat['Field'] = 'GOODS-N'
grat = pd.concat([cat_s_grat, cat_n_grat], ignore_index=True)
grat['NIRSpec_ID_str'] = grat['NIRSpec_ID'].astype(str).str.zfill(8)
grat['NIRCam_ID_int'] = pd.to_numeric(grat['NIRCam_ID'], errors='coerce')

# Tier A
spectra_ids = set()
with open(os.path.expanduser('~/lrd_workspace/analysis/tier_A_files.txt')) as f:
    for line in f:
        m = re.search(r'-(\d{8})_', line)
        if m:
            spectra_ids.add(m.group(1))
grat_tierA = grat[grat['NIRSpec_ID_str'].isin(spectra_ids)].copy()
print(f"\nTier A in gratings: {len(grat_tierA)}")

# --- 4. Merge ---
df_merge = grat_tierA.merge(
    df_phot[['ID', 'mag_F444W', 'mag_F277W', 'color_F444W_F277W', 'F444W_RHALF']],
    left_on='NIRCam_ID_int',
    right_on='ID',
    how='left'
)

# --- 5. LRD cuts ---
mask_valid = df_merge['mag_F444W'].notna() & df_merge['mag_F277W'].notna()
mask_color = df_merge['color_F444W_F277W'] > 1.5
mask_mag = df_merge['mag_F444W'].between(20, 30)
mask_lrd = mask_color & mask_mag & mask_valid

print(f"\n{'='*60}")
print("LRD CUT RESULTS (CORRECTED)")
print('='*60)
print(f"  Total Tier A:                   {len(df_merge)}")
print(f"  With photometric match:         {mask_valid.sum()}")
print(f"  Red color (F444W-F277W > 1.5):  {mask_color.sum()}")
print(f"  Valid magnitude (20-30):        {mask_mag.sum()}")
print(f"  🎯 Confirmed LRDs:              {mask_lrd.sum()}")
print(f"\n  Breakdown by field:")
print(df_merge[mask_lrd]['Field'].value_counts().to_string())

df_merge['is_LRD'] = mask_lrd
out_csv = os.path.expanduser('~/lrd_workspace/analysis/lrd_classification.csv')
df_merge.to_csv(out_csv, index=False)
print(f"\n💾 {out_csv}")

print(f"\nTop 10 LRDs by color:")
top = df_merge[mask_lrd].nlargest(10, 'color_F444W_F277W')
print(top[['NIRSpec_ID_str', 'Field', 'mag_F444W', 'mag_F277W',
           'color_F444W_F277W', 'F444W_RHALF']].to_string(index=False))
EOF