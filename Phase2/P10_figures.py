# P10_figures.py
"""
P10: Diagnostic figures for the paper.
- Stack with fit and residuals
- Bootstrap histogram
- Forest plot of subsamples
"""
import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.optimize import curve_fit
from pathlib import Path
import sys, gc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).parent))
from config import *

STACK_FILE = ANALYSIS / 'stack_v2_multiline.fits'
OUT_DIR = ANALYSIS / 'figures'
OUT_DIR.mkdir(exist_ok=True)


def model(w, a0, a1,
          A_hb, mu_hb, sig_hb,
          A_4959, mu_4959, sig_4959,
          A_5007, mu_5007, sig_5007,
          A_ha, mu_ha, sig_ha_nii, A_nii):
    cont = a0 + a1 * (w - 5700)
    D_NII = LAMBDA_VAC['[NII]6584'] - LAMBDA_VAC['Hα']
    g_hb   = A_hb   * np.exp(-0.5 * ((w - mu_hb)   / sig_hb)**2)
    g_4959 = A_4959 * np.exp(-0.5 * ((w - mu_4959) / sig_4959)**2)
    g_5007 = A_5007 * np.exp(-0.5 * ((w - mu_5007) / sig_5007)**2)
    g_ha   = A_ha   * np.exp(-0.5 * ((w - mu_ha)   / sig_ha_nii)**2)
    g_nii  = A_nii  * np.exp(-0.5 * ((w - (mu_ha + D_NII)) / sig_ha_nii)**2)
    return cont + g_hb + g_4959 + g_5007 + g_ha + g_nii


def gaussian(w, A, mu, sig):
    return A * np.exp(-0.5 * ((w - mu) / sig)**2)


def load_and_fit():
    stack = fits.open(STACK_FILE)
    d = stack[1].data
    wave = d['WAVE_REST']
    flux = d['FLUX']
    err = d['FLUX_ERR']

    cont = np.median(flux[(wave > 5200) & (wave < 5400)])
    def amp(lam, w=5):
        m = (wave > lam - w) & (wave < lam + w)
        return max(np.max(flux[m]) - cont, abs(cont) * 0.5)

    p0 = [cont, 0.0,
          amp(LAMBDA_VAC['Hβ']), LAMBDA_VAC['Hβ'], 2.5,
          amp(LAMBDA_VAC['[OIII]4959']), LAMBDA_VAC['[OIII]4959'], 2.7,
          amp(LAMBDA_VAC['[OIII]5007']), LAMBDA_VAC['[OIII]5007'], 2.7,
          amp(LAMBDA_VAC['Hα'], 10), LAMBDA_VAC['Hα'], 3.0,
          amp(LAMBDA_VAC['[NII]6584'], 10) * 0.3]

    popt, pcov = curve_fit(model, wave, flux, p0=p0, sigma=err,
                           absolute_sigma=True, maxfev=200000)
    perr = np.sqrt(np.diag(pcov))
    return wave, flux, err, popt, perr


# ============================================================
# FIGURE 1: Stack + fit + residuals
# ============================================================
def fig1_stack_fit():
    print("Figure 1: stack + fit + residuals...")
    wave, flux, err, popt, perr = load_and_fit()

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True,
                                    gridspec_kw={'height_ratios': [3, 1]})

    # Top panel: data + fit
    ax1.errorbar(wave, flux, yerr=err, fmt='.', ms=1, color='gray',
                 alpha=0.5, label='Stack (212 galaxies)', zorder=1)
    ax1.plot(wave, model(wave, *popt), 'r-', lw=1.5, label='Multi-line fit', zorder=3)

    # Individual components
    cont = popt[0] + popt[1] * (wave - 5700)
    ax1.plot(wave, cont, 'b--', lw=0.8, alpha=0.5, label='Continuum')

    D_NII = LAMBDA_VAC['[NII]6584'] - LAMBDA_VAC['Hα']
    ax1.plot(wave, cont + gaussian(wave, popt[2], popt[3], popt[4]),
             'g--', lw=0.7, alpha=0.6, label='Hβ')
    ax1.plot(wave, cont + gaussian(wave, popt[5], popt[6], popt[7]),
             'm--', lw=0.7, alpha=0.6, label='[OIII] 4959')
    ax1.plot(wave, cont + gaussian(wave, popt[8], popt[9], popt[10]),
             'c--', lw=0.7, alpha=0.6, label='[OIII] 5007')
    ax1.plot(wave, cont + gaussian(wave, popt[11], popt[12], popt[13]),
             'y--', lw=0.7, alpha=0.6, label='Hα')
    ax1.plot(wave, cont + gaussian(wave, popt[14], popt[12] + D_NII, popt[13]),
             'orange', ls='--', lw=0.7, alpha=0.6, label='[NII] 6584')

    # Theoretical λ markers
    for name, lam in [('Hβ', LAMBDA_VAC['Hβ']),
                      ('[OIII]', LAMBDA_VAC['[OIII]4959']),
                      ('[OIII]', LAMBDA_VAC['[OIII]5007']),
                      ('Hα', LAMBDA_VAC['Hα']),
                      ('[NII]', LAMBDA_VAC['[NII]6584'])]:
        ax1.axvline(lam, color='k', ls=':', lw=0.5, alpha=0.3)

    ax1.set_ylabel('Flux (W/m²)', fontsize=12)
    ax1.set_title('Rest-frame stack of 212 galaxies at z ∈ [5, 8] (JADES DR3)',
                  fontsize=13)
    ax1.legend(fontsize=8, ncol=3, loc='upper right')
    ax1.set_xlim(4800, 6700)
    ax1.grid(alpha=0.2)

    # Bottom panel: residuals
    resid = (flux - model(wave, *popt)) / err
    ax2.axhline(0, color='k', lw=0.5)
    ax2.axhline(1, color='r', ls='--', lw=0.5, alpha=0.5)
    ax2.axhline(-1, color='r', ls='--', lw=0.5, alpha=0.5)
    ax2.plot(wave, resid, '.', ms=1, color='gray', alpha=0.5)
    ax2.set_xlabel('Rest-frame wavelength (Å)', fontsize=12)
    ax2.set_ylabel('Residual (σ)', fontsize=12)
    ax2.set_ylim(-5, 5)
    ax2.grid(alpha=0.2)

    plt.tight_layout()
    out = OUT_DIR / 'fig1_stack_fit.png'
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out}")


# ============================================================
# FIGURE 2: Bootstrap histogram
# ============================================================
def fig2_bootstrap():
    print("Figure 2: bootstrap histogram...")
    boot = pd.read_csv(OUT_BOOTSTRAP)
    vals = boot['dalpha_over_alpha'].values * 1e6  # ppm

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(vals, bins=25, color='steelblue', alpha=0.7, edgecolor='black')
    ax.axvline(0, color='k', ls='-', lw=1, label='Δα/α = 0')
    ax.axvline(vals.mean(), color='r', ls='--', lw=1.5,
               label=f'Mean = {vals.mean():+.2f} ppm')
    ax.axvspan(vals.mean() - vals.std(), vals.mean() + vals.std(),
               alpha=0.2, color='r', label=f'1σ = {vals.std():.2f} ppm')

    ax.set_xlabel('Δα/α (ppm)', fontsize=12)
    ax.set_ylabel('Frequency', fontsize=12)
    ax.set_title('Bootstrap of the α test (100 iterations, 212 galaxies)',
                 fontsize=13)
    ax.legend(fontsize=10)
    ax.grid(alpha=0.2)

    plt.tight_layout()
    out = OUT_DIR / 'fig2_bootstrap.png'
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out}")


# ============================================================
# FIGURE 3: Forest plot of subsamples
# ============================================================
def fig3_forest():
    print("Figure 3: forest plot of subsamples...")

    # P07 data (from output)
    labels = [
        'GOODS-S', 'GOODS-N',
        'z ∈ [5.0, 5.5]', 'z ∈ [5.5, 6.0]',
        'z ∈ [6.0, 6.5]', 'z ∈ [6.5, 8.0]',
        'S/N Q1', 'S/N Q2', 'S/N Q3', 'S/N Q4',
        'Full sample'
    ]
    means = [2.68, 1.75, -0.20, 0.83, -14.14, 16.67,
             2.65, -6.14, 17.77, 3.18, 3.02]
    errs = [2.35, 2.88, 2.45, 2.86, 9.48, 5.72,
            13.74, 7.37, 4.72, 2.13, 5.73]

    fig, ax = plt.subplots(figsize=(8, 6))

    y_pos = np.arange(len(labels))
    colors = ['steelblue'] * 10 + ['red']

    ax.errorbar(means, y_pos, xerr=errs, fmt='o', color='k',
                capsize=4, ms=5, elinewidth=1.5)
    ax.axvline(0, color='k', ls='-', lw=0.8)
    ax.axvspan(-5.73, 5.73, alpha=0.1, color='gray',
               label='1σ of the full sample')

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel('Δα/α (ppm)', fontsize=12)
    ax.set_title('α test in subsamples', fontsize=13)
    ax.legend(fontsize=9, loc='lower right')
    ax.grid(alpha=0.2, axis='x')

    # Mark outliers
    for i, (m, e) in enumerate(zip(means, errs)):
        if abs(m/e) > 2.5:
            ax.annotate('', xy=(m, i), xytext=(m + 1.5*np.sign(m), i),
                        arrowprops=dict(arrowstyle='->', color='red', lw=1))

    plt.tight_layout()
    out = OUT_DIR / 'fig3_forest.png'
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out}")


# ============================================================
# FIGURE 4: Z distribution
# ============================================================
def fig4_z_dist():
    print("Figure 4: z distribution...")
    tier_a = pd.read_csv(ANALYSIS / 'tier_A_z5_v2.csv')

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # z histogram
    ax1.hist(tier_a['z'], bins=20, color='steelblue',
             alpha=0.7, edgecolor='black')
    ax1.axvline(5.0, color='r', ls='--', lw=1, label='z = 5')
    ax1.axvline(7.0, color='orange', ls='--', lw=1, label='z = 7')
    ax1.set_xlabel('z', fontsize=12)
    ax1.set_ylabel('N galaxies', fontsize=12)
    ax1.set_title(f'z distribution (N = {len(tier_a)})', fontsize=13)
    ax1.legend()
    ax1.grid(alpha=0.2)

    # S/N histogram
    ax2.hist(tier_a['snr_peak'], bins=30, color='darkorange',
             alpha=0.7, edgecolor='black')
    ax2.axvline(10, color='g', ls='--', lw=1, label='Tier A (S/N > 10)')
    ax2.set_xlabel('Peak S/N', fontsize=12)
    ax2.set_ylabel('N galaxies', fontsize=12)
    ax2.set_title('Peak S/N distribution', fontsize=13)
    ax2.set_yscale('log')
    ax2.legend()
    ax2.grid(alpha=0.2)

    plt.tight_layout()
    out = OUT_DIR / 'fig4_z_snr.png'
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out}")


def main():
    print("=" * 70)
    print("P10: DIAGNOSTIC FIGURES")
    print("=" * 70)
    print(f"Stack: {STACK_FILE}")
    print(f"Output: {OUT_DIR}")
    print()

    fig1_stack_fit()
    fig2_bootstrap()
    fig3_forest()
    fig4_z_dist()

    print()
    print("=" * 70)
    print("FIGURES COMPLETE")
    print("=" * 70)
    print(f"Directory: {OUT_DIR}")
    print()
    for f in sorted(OUT_DIR.glob('*.png')):
        size = f.stat().st_size / 1024
        print(f"  {f.name}  ({size:.1f} KB)")


if __name__ == '__main__':
    main()
