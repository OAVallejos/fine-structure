cd ~/lrd_workspace
python3 << 'EOF'
import numpy as np
import pandas as pd
from astropy.table import Table
import os

cat_s_grat = Table.read(os.path.expanduser('~/catalogs_jades/hlsp_jades_jwst_nirspec_goods-s_gratings-line-fluxes_v1.1_catalog.fits'), hdu=1).to_pandas()
cat_n_grat = Table.read(os.path.expanduser('~/catalogs_jades/hlsp_jades_jwst_nirspec_goods-n_gratings-line-fluxes_v1.1_catalog.fits'), hdu=1).to_pandas()
cat_s_grat['Field'] = 'GOODS-S'
cat_n_grat['Field'] = 'GOODS-N'

def load_phot(cat_path):
    size = Table.read(cat_path, hdu=3).to_pandas()
    kron = Table.read(cat_path, hdu=8).to_pandas()
    df = size[['ID', 'RA', 'DEC', 'F444W_RHALF']].merge(
        kron[['ID', 'F444W_KRON', 'F277W_KRON']], on='ID', how='inner'
    )
    def flux_to_mag(f):
        return -2.5 * np.log10(f / 3631e9) if f > 0 else np.nan
    df['mag_F444W'] = df['F444W_KRON'].apply(flux_to_mag)
    df['mag_F277W'] = df['F277W_KRON'].apply(flux_to_mag)
    df['color_F444W_F277W'] = df['mag_F444W'] - df['mag_F277W']
    return df

df_phot_s = load_phot(os.path.expanduser('~/catalogs_jades/hlsp_jades_jwst_nircam_goods-s-deep_photometry_v2.0_catalog.fits'))
df_phot_n = load_phot(os.path.expanduser('~/catalogs_jades/hlsp_jades_jwst_nircam_goods-n_photometry_v1.0_catalog.fits'))

cat_s_grat['NIRSpec_ID_str'] = cat_s_grat['NIRSpec_ID'].astype(str).str.zfill(8)
cat_s_grat['NIRCam_ID_int'] = pd.to_numeric(cat_s_grat['NIRCam_ID'], errors='coerce')
grat_s = cat_s_grat[cat_s_grat['z_Spec'] >= 4.0].merge(
    df_phot_s[['ID', 'mag_F444W', 'mag_F277W', 'color_F444W_F277W', 'F444W_RHALF']],
    left_on='NIRCam_ID_int', right_on='ID', how='left'
)

cat_n_grat['NIRSpec_ID_str'] = cat_n_grat['NIRSpec_ID'].astype(str).str.zfill(8)
cat_n_grat['NIRCam_ID_int'] = pd.to_numeric(cat_n_grat['NIRCam_ID'], errors='coerce')
grat_n = cat_n_grat[cat_n_grat['z_Spec'] >= 4.0].merge(
    df_phot_n[['ID', 'mag_F444W', 'mag_F277W', 'color_F444W_F277W', 'F444W_RHALF']],
    left_on='NIRCam_ID_int', right_on='ID', how='left'
)

df_merge = pd.concat([grat_s, grat_n], ignore_index=True)

mask_valid = df_merge['mag_F444W'].notna()
mask_color = df_merge['color_F444W_F277W'] > 0.5
mask_mag = df_merge['mag_F444W'].between(20, 29)
mask_compact = df_merge['F444W_RHALF'] < 0.15
mask_lrd = mask_valid & mask_color & mask_mag & mask_compact

print(f"{'='*60}")
print("LRD CANDIDATES z >= 4 (CORRECTED MERGE)")
print('='*60)
print(f"  Total z >= 4:              {len(df_merge)}")
print(f"  With photometric match:    {mask_valid.sum()}")
print(f"  Color > 0.5:               {mask_color.sum()}")
print(f"  Magnitude 20-29:           {mask_mag.sum()}")
print(f"  Compact (RHALF < 0.15\"):  {mask_compact.sum()}")
print(f"  🎯 LRD candidates:         {mask_lrd.sum()}")

print(f"\n  Breakdown by field:")
print(df_merge[mask_lrd]['Field'].value_counts().to_string())

df_merge['is_LRD_candidate'] = mask_lrd
out_csv = os.path.expanduser('~/lrd_workspace/catalogs/lrd_candidates_z4_v2.csv')
df_merge.to_csv(out_csv, index=False)
print(f"\n💾 {out_csv}")
EOF