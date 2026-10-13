#!/usr/bin/env python3
"""
build_master_dataset_v4.py (DEFINITIVE)
Builds the Master Table fixing the 8-digit regex and removing the non-existent radius column.
"""
import os
import re
import pandas as pd
import numpy as np
from astropy.table import Table

WORKSPACE = os.path.expanduser('~/lrd_workspace')
ANALYSIS = os.path.join(WORKSPACE, 'analysis')

print("🔨 Building Master Table v4 (Definitive)...")

# 1. Load NIRSpec catalogs
print("   [1/5] Loading NIRSpec catalogs...")
cat_s = Table.read(os.path.expanduser('~/catalogs_jades/hlsp_jades_jwst_nirspec_goods-s_gratings-line-fluxes_v1.1_catalog.fits'), hdu=1).to_pandas()
cat_n = Table.read(os.path.expanduser('~/catalogs_jades/hlsp_jades_jwst_nirspec_goods-n_gratings-line-fluxes_v1.1_catalog.fits'), hdu=1).to_pandas()

cat_s['Field'] = 'GOODS-S'
cat_n['Field'] = 'GOODS-N'
df_nirspec = pd.concat([cat_s, cat_n], ignore_index=True)

# Ensure NIRSpec_ID is an 8-digit string for the merge
df_nirspec['NIRSpec_ID'] = df_nirspec['NIRSpec_ID'].astype(str).str.zfill(8)
print(f"   Total NIRSpec entries: {len(df_nirspec)}")

# 2. Load quality_metrics and EXTRACT ID with the CORRECT regex (8 digits)
print("   [2/5] Loading quality_metrics.csv and extracting IDs...")
df_metrics = pd.read_csv(os.path.join(ANALYSIS, 'quality_metrics.csv'))

def extract_id_8digits(filename):
    # Looks for exactly 8 digits between dashes: e.g. "-00000607_"
    match = re.search(r'-(\d{8})_', str(filename))
    if match:
        return match.group(1)
    return None

df_metrics['NIRSpec_ID'] = df_metrics['file'].apply(extract_id_8digits)
df_metrics['NIRSpec_ID'] = df_metrics['NIRSpec_ID'].astype(str).str.zfill(8)

valid_ids = df_metrics['NIRSpec_ID'].notna().sum()
print(f"   ✅ IDs extracted successfully: {valid_ids} / {len(df_metrics)}")

# 3. Load photometric catalogs
print("   [3/5] Loading photometric catalogs...")
phot_s = Table.read(os.path.expanduser('~/catalogs_jades/hlsp_jades_jwst_nircam_goods-s-deep_photometry_v2.0_catalog.fits'), hdu=8).to_pandas()
phot_n = Table.read(os.path.expanduser('~/catalogs_jades/hlsp_jades_jwst_nircam_goods-n_photometry_v1.0_catalog.fits'), hdu=8).to_pandas()

phot_s['Field'] = 'GOODS-S'
phot_n['Field'] = 'GOODS-N'
df_phot = pd.concat([phot_s, phot_n], ignore_index=True)

# Compute magnitudes (avoiding fragmentation)
if df_phot['F444W_KRON'].median() < 1000:
    mag_f444w = 23.9 - 2.5 * np.log10(df_phot['F444W_KRON'].replace(0, np.nan))
    mag_f277w = 23.9 - 2.5 * np.log10(df_phot['F277W_KRON'].replace(0, np.nan))
else:
    mag_f444w = df_phot['F444W_KRON']
    mag_f277w = df_phot['F277W_KRON']

phot_additions = pd.DataFrame({
    'MAG_F444W': mag_f444w,
    'MAG_F277W': mag_f277w,
    'COLOR_F444W_F277W': mag_f444w - mag_f277w
})
df_phot = pd.concat([df_phot, phot_additions], axis=1)

# Ensure NIRCam_ID is an 8-digit string
df_phot['ID'] = df_phot['ID'].astype(str).str.zfill(8)

# 4. Master merge
print("   [4/5] Merging datasets...")

# Merge 1: metrics + nirspec (by NIRSpec_ID)
df_master = df_metrics.merge(df_nirspec, on='NIRSpec_ID', how='left')
match_nirspec = df_master['z_Spec'].notna().sum()
print(f"   ✅ Match with NIRSpec catalog: {match_nirspec} / {len(df_master)}")

# Merge 2: with photometry (by NIRCam_ID and Field)
if 'NIRCam_ID' in df_master.columns:
    df_master['NIRCam_ID'] = df_master['NIRCam_ID'].astype(str).str.zfill(8)

    # NOTE: 'F444W_RHALF' was removed because the diagnostic confirmed it does not exist in this extension.
    phot_cols = ['ID', 'Field', 'MAG_F444W', 'COLOR_F444W_F277W']

    df_master = df_master.merge(
        df_phot[phot_cols].rename(columns={'ID': 'NIRCam_ID_phot'}),
        left_on=['NIRCam_ID', 'Field'],
        right_on=['NIRCam_ID_phot', 'Field'],
        how='left'
    )
    df_master.drop(columns=['NIRCam_ID_phot'], inplace=True, errors='ignore')
    match_phot = df_master['MAG_F444W'].notna().sum()
    print(f"   ✅ Match with Photometric catalog: {match_phot} / {len(df_master)}")
else:
    print("   ⚠️  WARNING: NIRCam_ID not found for merge with photometry")

# 5. Cleanup for CatBoost
print("   [5/5] Preparing for CatBoost...")
cat_cols = ['GRATING', 'FILTER', 'Field', 'TIER', 'PID', 'SRCTYPE']
for col in cat_cols:
    if col in df_master.columns:
        df_master[col] = df_master[col].fillna('Missing').astype(str)

# Save
out_path = os.path.join(ANALYSIS, 'master_features_v4.csv')
df_master.to_csv(out_path, index=False)

print(f"\n✅ MASTER TABLE V4 SUCCESSFULLY CREATED: {out_path}")
print(f"   Total rows: {len(df_master)}")
print(f"   Total columns: {len(df_master.columns)}")

# Final diagnostic
print("\n📊 COLUMNS WITH >50% VALID DATA (Ready for CatBoost):")
for col in df_master.columns:
    n_valid = df_master[col].notna().sum()
    pct = 100 * n_valid / len(df_master)
    if pct > 50:
        print(f"   {col:30s}: {n_valid:5d} / {len(df_master)}  ({pct:.0f}%)")