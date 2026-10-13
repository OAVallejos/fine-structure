# CHECK_Z65.py
import numpy as np
import pandas as pd
from astropy.io import fits
import warnings
warnings.filterwarnings('ignore')

tier_A = pd.read_csv('lrd_workspace/analysis/tier_A_unique_z5.csv')
tier_A['ID'] = tier_A['ID'].astype(str).str.zfill(8)

# Bin z > 6.5
sub = tier_A[tier_A['z'] > 6.5]
print(f"Objects with z > 6.5: {len(sub)}")
print()

LAB_HA = 6564.63
LAB_HB = 4862.68
LAB_OIII = 5008.24

print(f"{'ID':10s} {'z':>7s} {'Hα_obs':>10s} {'Hβ_obs':>10s} {'[OIII]':>10s} {'x1d_range':>22s} {'Hα?':>5s}")
n_ha_in = 0
for _, row in sub.iterrows():
    z = row['z']
    lam_ha = LAB_HA * (1 + z)
    lam_hb = LAB_HB * (1 + z)
    lam_oiii = LAB_OIII * (1 + z)
    try:
        with fits.open(f'lrd_workspace/spectra/{row["file"]}') as hdul:
            w = hdul[1].data['WAVELENGTH']
        if np.nanmax(w) < 100: w = w * 1e4
        wmin, wmax = w.min(), w.max()
        ha_in = 'YES' if (wmin < lam_ha < wmax) else 'NO'
        if ha_in == 'YES': n_ha_in += 1
        print(f"{row['ID']:10s} {z:7.3f} {lam_ha:10.0f} {lam_hb:10.0f} {lam_oiii:10.0f} "
              f"[{wmin:8.0f},{wmax:8.0f}] {ha_in:>5s}")
    except Exception as e:
        print(f"{row['ID']:10s} {z:7.3f} ERROR: {e}")

print()
print(f"Hα within range: {n_ha_in}/{len(sub)}")