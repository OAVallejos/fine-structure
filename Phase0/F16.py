cd ~/lrd_workspace/analysis
python3 << 'EOF'
import numpy as np
import pandas as pd
from astropy.table import Table
import os

# --- 1. Load GOODS-N DR3 photometric catalog ---
CAT = os.path.expanduser('~/catalogs_jades/hlsp_jades_jwst_nircam_goods-n_photometry_v1.0_catalog.fits')

# Key extensions
ext_size = 3   # SIZE (RHALF)
ext_kron = 8   # KRON_CONV (PSF-matched magnitudes)
ext_photoz = 9 # PHOTOZ

print("Loading catalog...")
size = Table.read(CAT, hdu=ext_size).to_pandas()
kron = Table.read(CAT, hdu=ext_kron).to_pandas()
photoz = Table.read(CAT, hdu=ext_photoz).to_pandas()

print(f"Total sources: {len(size)}")
print(f"SIZE columns: {size.columns.tolist()[:10]}...")
print(f"KRON_CONV columns: {kron.columns.tolist()[:10]}...")

# --- 2. Extract the columns we need ---
# Check exactly which columns exist
print("\n--- Columns with F444W ---")
cols_f444 = [c for c in kron.columns if 'F444W' in c.upper()]
print(cols_f444[:10])

print("\n--- Columns with F277W ---")
cols_f277 = [c for c in kron.columns if 'F277W' in c.upper()]
print(cols_f277[:10])

print("\n--- Columns with RHALF ---")
cols_rhal = [c for c in size.columns if 'RHALF' in c.upper()]
print(cols_rhal[:10])

# --- 3. Merge by ID ---
df = size[['ID', 'RA', 'DEC'] + cols_rhal].merge(
    kron[['ID'] + cols_f444 + cols_f277],
    on='ID', how='inner'
)
print(f"\nMerge: {len(df)} sources")

# --- 4. Apply LRD cuts ---
# We need magnitudes, not fluxes. Check if they are mag or flux
print("\n--- Example F444W values ---")
print(df[cols_f444[:3]].describe())

# If they are magnitudes (values ~20-30), apply color:
# F444W - F277W > 1.5
col_f444 = [c for c in df.columns if 'F444W_KRON' in c.upper() or 'F444W' in c.upper()][0]
col_f277 = [c for c in df.columns if 'F277W_KRON' in c.upper() or 'F277W' in c.upper()][0]

# Check if they are magnitudes or fluxes
sample_val = df[col_f444].dropna().iloc[0]
print(f"\nExample value {col_f444}: {sample_val}")
if sample_val > 100:
    print("→ They are FLUXES (nJy), converting to magnitudes")
    def flux_to_mag(f):
        return -2.5 * np.log10(f / 3631e9) if f > 0 else np.nan
    df['mag_F444W'] = df[col_f444].apply(flux_to_mag)
    df['mag_F277W'] = df[col_f277].apply(flux_to_mag)
else:
    print("→ They are MAGNITUDES")
    df['mag_F444W'] = df[col_f444]
    df['mag_F277W'] = df[col_f277]

# Color
df['color_F444W_F277W'] = df['mag_F444W'] - df['mag_F277W']

# --- 5. Load ID list of the 330 grating spectra ---
spectra_ids = set()
with open(os.path.expanduser('~/lrd_workspace/analysis/tier_A_files.txt')) as f:
    for line in f:
        import re
        m = re.search(r'-(\d{8})_', line)
        if m:
            spectra_ids.add(m.group(1))

print(f"\nTier A spectra IDs: {len(spectra_ids)}")

# --- 6. Cross-match with the photometric catalog ---
# The NIRCam_ID of the gratings catalog is the same as the ID of the photometric catalog
# But we need to match by file ID or by NIRCam_ID
# Load the gratings catalog to get the NIRCam_ID
cat_grat = Table.read(os.path.expanduser('~/catalogs_jades/hlsp_jades_jwst_nirspec_goods-n_gratings-line-fluxes_v1.1_catalog.fits'), hdu=1).to_pandas()
cat_grat['NIRSpec_ID_str'] = cat_grat['NIRSpec_ID'].astype(str).str.zfill(8)
cat_grat['NIRCam_ID_str'] = cat_grat['NIRCam_ID'].astype(str)

# Filter only Tier A
tier_a_grat = cat_grat[cat_grat['NIRSpec_ID_str'].isin(spectra_ids)]
print(f"Tier A in gratings catalog: {len(tier_a_grat)}")

# Merge with photometry by NIRCam_ID
df_merge = tier_a_grat.merge(
    df[['ID', 'mag_F444W', 'mag_F277W', 'color_F444W_F277W'] + cols_rhal],
    left_on='NIRCam_ID',
    right_on='ID',
    how='left'
)
print(f"Merge with photometry: {df_merge['mag_F444W'].notna().sum()}/{len(df_merge)} with data")

# --- 7. Apply LRD cuts ---
mask_lrd = (
    (df_merge['color_F444W_F277W'] > 1.5) &
    (df_merge['mag_F444W'].between(22, 29)) &
    (df_merge['mag_F444W'].notna())
)
n_lrd = mask_lrd.sum()
print(f"\n🎯 Confirmed LRDs: {n_lrd}/{len(df_merge)}")

# Save
df_merge['is_LRD'] = mask_lrd
out_csv = os.path.expanduser('~/lrd_workspace/analysis/lrd_classification.csv')
df_merge.to_csv(out_csv, index=False)
print(f"💾 {out_csv}")
EOF