# FILTER_HA.py
import numpy as np
import pandas as pd
from astropy.io import fits
import warnings
warnings.filterwarnings('ignore')

tier_A = pd.read_csv('lrd_workspace/analysis/tier_A_unique_z5.csv')
tier_A['ID'] = tier_A['ID'].astype(str).str.zfill(8)

LAB_HA = 6564.63

# Filter: Hα must be within range
valid = []
for _, row in tier_A.iterrows():
    lam_ha = LAB_HA * (1 + row['z'])
    try:
        with fits.open(f'lrd_workspace/spectra/{row["file"]}') as hdul:
            w = hdul[1].data['WAVELENGTH']
        if np.nanmax(w) < 100: w = w * 1e4
        if w.min() < lam_ha < w.max():
            valid.append(row)
    except Exception:
        pass

df_valid = pd.DataFrame(valid)
print(f"Total objects: {len(tier_A)}")
print(f"Objects with Hα within range: {len(df_valid)}")
print(f"Objects discarded: {len(tier_A) - len(df_valid)}")
print()
print(f"z distribution (valid Hα):")
print(df_valid['z'].describe())
print()
print(f"field distribution:")
print(df_valid['field'].value_counts())
print()
print(f"grating distribution:")
print(df_valid['grating'].value_counts())

df_valid.to_csv('lrd_workspace/analysis/tier_A_ha_valid.csv', index=False)
print(f"\nSaved: tier_A_ha_valid.csv")
