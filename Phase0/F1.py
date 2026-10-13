#!/usr/bin/env python3
"""
Phase 1: Quality assessment of x1d spectra for LRD candidates.
Runs from any directory using absolute paths.
"""
import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits

# --- Absolute paths ---
BASE = os.path.expanduser('~/lrd_workspace')
SPECTRA_DIR = os.path.join(BASE, 'spectra')
OUT_DIR = os.path.join(BASE, 'analysis')
os.makedirs(OUT_DIR, exist_ok=True)

files = sorted(glob.glob(os.path.join(SPECTRA_DIR, '*.fits')))
print(f"Spectra directory: {SPECTRA_DIR}")
print(f"Total files to evaluate: {len(files)}\n")

if len(files) == 0:
    print("❌ No files found. Check that ~/lrd_workspace/spectra/ contains .fits files")
    raise SystemExit(1)

snr_peaks = []
valid_pixels = []
filenames = []

for f in files:
    with fits.open(f) as h:
        wave = h[1].data['WAVELENGTH']
        flux = h[1].data['FLUX']
        flux_err = h[1].data['FLUX_ERR']
        
        mask = np.isfinite(flux) & np.isfinite(flux_err) & (flux_err > 0)
        wave_c = wave[mask]
        flux_c = flux[mask]
        err_c = flux_err[mask]
        
        # Trim edges (instrumental garbage)
        if len(wave_c) > 200:
            wave_c = wave_c[50:-50]
            flux_c = flux_c[50:-50]
            err_c = err_c[50:-50]
        
        valid_pixels.append(len(wave_c))
        filenames.append(os.path.basename(f))
        
        if len(wave_c) > 100:
            snr = np.abs(flux_c) / err_c
            snr_peaks.append(float(np.max(snr)))
        else:
            snr_peaks.append(0.0)

# Save metrics
df = pd.DataFrame({
    'file': filenames,
    'n_valid_pixels': valid_pixels,
    'snr_peak': snr_peaks
})
csv_path = os.path.join(OUT_DIR, 'quality_metrics.csv')
df.to_csv(csv_path, index=False)
print(f"💾 Metrics saved to: {csv_path}\n")

# Visualization
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

axes[0].hist(valid_pixels, bins=50, color='steelblue', edgecolor='black')
axes[0].set_title('Valid Pixels per Spectrum')
axes[0].set_xlabel('Number of pixels')
axes[0].set_ylabel('Frequency')

axes[1].hist(snr_peaks, bins=50, color='darkorange', edgecolor='black')
axes[1].set_title('Peak S/N per Spectrum')
axes[1].set_xlabel('Max S/N')
axes[1].set_ylabel('Frequency')
axes[1].axvline(x=5, color='red', linestyle='--', label='S/N = 5')
axes[1].legend()

plt.tight_layout()
png_path = os.path.join(OUT_DIR, 'data_quality_check.png')
plt.savefig(png_path, dpi=300)
print(f"✅ Plot saved to: {png_path}\n")

# Summary
snr_arr = np.array(snr_peaks)
print("📊 Summary:")
print(f"  Total spectra: {len(files)}")
print(f"  S/N > 3:  {np.sum(snr_arr > 3)}")
print(f"  S/N > 5:  {np.sum(snr_arr > 5)}")
print(f"  S/N > 10: {np.sum(snr_arr > 10)}")
print(f"  Mean S/N:  {np.mean(snr_arr):.2f}")
print(f"  Median S/N: {np.median(snr_arr):.2f}")
print()
print("🏆 Top 20 spectra by peak S/N:")
top = df.sort_values('snr_peak', ascending=False).head(20)
print(top.to_string(index=False))