# P24_eda_lrd.py
"""
P24: Exhaustive EDA of LRDs.
- LRD criteria (color, magnitude, size, z)
- Line behavior vs z
- Outliers
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits
from pathlib import Path
import sys, re, gc
sys.path.insert(0, str(Path(__file__).parent))
from config import *

SPECTRA_LRD = ROOT / 'spectra_lrd_candidates'
OUT_DIR = ANALYSIS / 'figures'
OUT_DIR.mkdir(exist_ok=True)

# ============================================================
# 1. LOAD CATALOG
# ============================================================
print("=" * 70)
print("P24: EDA OF LRDs")
print("=" * 70)

df = pd.read_csv(ANALYSIS / 'lrd_z4.csv')
df['ID'] = df['ID'].astype(str).str.zfill(8)
print(f"LRDs loaded: {len(df)}")

# Load photometric features
cat = pd.read_csv(ROOT / 'catalogs' / 'lrd_candidates_z4_v2.csv')
cat['ID'] = cat['NIRSpec_ID'].astype(str).str.zfill(8)

# Merge
df = df.merge(cat[['ID', 'mag_F444W', 'mag_F277W',
                    'color_F444W_F277W', 'F444W_RHALF',
                    'is_LRD_candidate']],
              on='ID', how='left')

# ============================================================
# 2. LRD CRITERIA
# ============================================================
print("\n" + "=" * 70)
print("LRD CRITERIA")
print("=" * 70)

c1 = df['color_F444W_F277W'] > 0.5
c2 = (df['mag_F444W'] > 20) & (df['mag_F444W'] < 29)
c3 = df['F444W_RHALF'] < 0.15
c4 = df['z'] > 4.0

print(f"\nColor (F444W-F277W > 0.5):     {c1.sum()}/{len(df)}")
print(f"Magnitude (20 < F444W < 29):   {c2.sum()}/{len(df)}")
print(f"Size (RHALF < 0.15\"):          {c3.sum()}/{len(df)}")
print(f"Redshift (z > 4):              {c4.sum()}/{len(df)}")
print(f"ALL criteria:                  {(c1 & c2 & c3 & c4).sum()}/{len(df)}")

df['LRD_confirmado'] = c1 & c2 & c3 & c4

# ============================================================
# 3. DASHBOARD
# ============================================================
fig, axes = plt.subplots(2, 3, figsize=(18, 10))
fig.suptitle('EDA: LRDs in JADES DR3', fontsize=16, fontweight='bold')

# (0,0) z distribution
ax = axes[0, 0]
ax.hist(df['z'], bins=20, color='purple', alpha=0.7, edgecolor='black')
ax.axvline(4.0, color='r', ls='--', label='z = 4')
ax.axvline(df['z'].median(), color='g', ls='-', label=f'Median = {df["z"].median():.2f}')
ax.set_xlabel('z')
ax.set_ylabel('N')
ax.set_title('z distribution')
ax.legend()

# (0,1) Color vs magnitude
ax = axes[0, 1]
sc = ax.scatter(df['mag_F444W'], df['color_F444W_F277W'],
                c=df['z'], cmap='viridis', s=30, alpha=0.7, edgecolor='k')
ax.axhline(0.5, color='r', ls='--', label='color = 0.5')
ax.axvline(20, color='orange', ls=':', label='mag = 20')
ax.set_xlabel('mag_F444W')
ax.set_ylabel('F444W − F277W')
ax.set_title('Color vs Magnitude')
ax.invert_xaxis()
plt.colorbar(sc, ax=ax, label='z')
ax.legend(fontsize=8)

# (0,2) Size (RHALF)
ax = axes[0, 2]
ax.hist(df['F444W_RHALF'].dropna(), bins=20, color='teal', alpha=0.7, edgecolor='black')
ax.axvline(0.15, color='r', ls='--', label='0.15"')
ax.set_xlabel('F444W RHALF (arcsec)')
ax.set_ylabel('N')
ax.set_title('Effective size')
ax.legend()

# (1,0) S/N
ax = axes[1, 0]
ax.hist(df['snr_peak'].dropna(), bins=30, color='coral', alpha=0.7, edgecolor='black')
ax.axvline(10, color='r', ls='--', label='S/N = 10')
ax.set_xlabel('Peak S/N')
ax.set_ylabel('N')
ax.set_title('S/N distribution')
ax.set_yscale('log')
ax.legend()

# (1,1) Lines vs z (observed Hβ)
ax = axes[1, 1]
lam_hb = 4862.68
lam_ha = 6564.63
df['HB_obs'] = lam_hb * (1 + df['z'])
df['HA_obs'] = lam_ha * (1 + df['z'])
ax.scatter(df['z'], df['HB_obs'] / 1e4, c='blue', s=20, alpha=0.7, label='Hβ')
ax.scatter(df['z'], df['HA_obs'] / 1e4, c='red', s=20, alpha=0.7, label='Hα')
ax.axhline(2.87, color='gray', ls='--', lw=0.8, label='G395M min')
ax.axhline(5.27, color='gray', ls='--', lw=0.8, label='G395M max')
ax.set_xlabel('z')
ax.set_ylabel('λ_obs (µm)')
ax.set_title('Lines vs z')
ax.legend(fontsize=8)
ax.grid(alpha=0.2)

# (1,2) Hβ in range
ax = axes[1, 2]
df['HB_in_range'] = (df['HB_obs'] > 28700) & (df['HB_obs'] < 52687)
df['HA_in_range'] = (df['HA_obs'] > 28700) & (df['HA_obs'] < 52687)
ax.scatter(df['z'], df['HB_in_range'].astype(int) + 0.02 * np.random.randn(len(df)),
           c='blue', s=30, alpha=0.7, label='Hβ in range')
ax.scatter(df['z'], df['HA_in_range'].astype(int) - 0.02 * np.random.randn(len(df)),
           c='red', s=30, alpha=0.7, label='Hα in range')
ax.set_xlabel('z')
ax.set_ylabel('In G395M range')
ax.set_title('Spectral coverage vs z')
ax.legend(fontsize=8)
ax.grid(alpha=0.2)

plt.tight_layout()
out = OUT_DIR / 'fig9_eda_lrd.png'
plt.savefig(out, dpi=150, bbox_inches='tight')
plt.close()
print(f"\nSaved: {out}")

# ============================================================
# 4. OUTLIERS
# ============================================================
print("\n" + "=" * 70)
print("OUTLIERS")
print("=" * 70)

# z outliers
z_q1 = df['z'].quantile(0.25)
z_q3 = df['z'].quantile(0.75)
z_iqr = z_q3 - z_q1
z_outliers = df[(df['z'] < z_q1 - 1.5*z_iqr) | (df['z'] > z_q3 + 1.5*z_iqr)]
print(f"\nz outliers: {len(z_outliers)}")
if len(z_outliers) > 0:
    print(z_outliers[['ID', 'z', 'snr_peak']].to_string())

# S/N outliers
snr_q1 = df['snr_peak'].quantile(0.25)
snr_q3 = df['snr_peak'].quantile(0.75)
snr_iqr = snr_q3 - snr_q1
snr_outliers = df[(df['snr_peak'] < snr_q1 - 1.5*snr_iqr) | (df['snr_peak'] > snr_q3 + 1.5*snr_iqr)]
print(f"\nS/N outliers: {len(snr_outliers)}")
if len(snr_outliers) > 0:
    print(snr_outliers[['ID', 'z', 'snr_peak']].head(20).to_string())

# color outliers
color_q1 = df['color_F444W_F277W'].quantile(0.25)
color_q3 = df['color_F444W_F277W'].quantile(0.75)
color_iqr = color_q3 - color_q1
color_outliers = df[(df['color_F444W_F277W'] < color_q1 - 1.5*color_iqr) |
                    (df['color_F444W_F277W'] > color_q3 + 1.5*color_iqr)]
print(f"\nColor outliers: {len(color_outliers)}")
if len(color_outliers) > 0:
    print(color_outliers[['ID', 'z', 'color_F444W_F277W', 'snr_peak']].head(20).to_string())

# ============================================================
# 5. LINE BEHAVIOR VS Z
# ============================================================
print("\n" + "=" * 70)
print("LINE BEHAVIOR VS Z")
print("=" * 70)

# Load gratings catalog with fluxes
gratings_s = pd.read_csv(ROOT / 'catalogs' / 'lrd_candidates_z4_v2.csv')
gratings_s['ID'] = gratings_s['NIRSpec_ID'].astype(str).str.zfill(8)

# Ratios per z bin
bins = [(4.0, 4.3), (4.3, 4.6), (4.6, 4.9), (4.9, 5.5), (5.5, 7.0)]
print("\nLine ratios per z bin:")
print(f"{'Bin':12s} {'N':>4s} {'OIII/Hβ':>10s} {'Hα/Hβ':>10s} {'NII/Hα':>10s}")
for zlo, zhi in bins:
    sub = df[(df['z'] >= zlo) & (df['z'] < zhi)]
    if len(sub) == 0: continue
    # Merge with fluxes
    m = sub.merge(gratings_s[['ID', 'HB_4861_flux', 'O3_5007_flux',
                              'HA_6563_flux', 'N2_6584_flux']],
                  on='ID', how='left')
    oiii_hb = (m['O3_5007_flux'] / m['HB_4861_flux']).median()
    ha_hb = (m['HA_6563_flux'] / m['HB_4861_flux']).median()
    nii_ha = (m['N2_6584_flux'] / m['HA_6563_flux']).median()
    print(f"[{zlo:.1f},{zhi:.1f}) {len(sub):4d} {oiii_hb:10.3f} {ha_hb:10.3f} {nii_ha:10.3f}")

# Save
df.to_csv(ANALYSIS / 'eda_lrd.csv', index=False)
print(f"\nSaved: eda_lrd.csv")


if __name__ == '__main__':
    main()
