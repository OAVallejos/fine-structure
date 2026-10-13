# config.py
"""
Central configuration of pipeline v2.
All constants in a single place.
"""

from pathlib import Path

# ============================================================
# PATHS
# ============================================================
ROOT = Path.home() / 'lrd_workspace'
ANALYSIS = ROOT / 'analysis'
CATALOGS = ROOT / 'catalogs'
SPECTRA = ROOT / 'spectra'
SPECTRA_LRD = ROOT / 'spectra_lrd_candidates'
ALL_FITS = Path.home() / 'todos_los_fits'
CATALOGS_JADES = Path.home() / 'catalogs_jades'

# ============================================================
# LABORATORY λ — ALL IN VACUUM
# ============================================================
LAMBDA_VAC = {
    'Hδ':        4102.89,
    'Hγ':        4341.69,
    'Hβ':        4862.68,
    'Hα':        6564.63,
    '[OII]3727': 3727.09,
    '[OIII]4959': 4960.30,
    '[OIII]5007': 5008.24,
    '[OI]6300':  6301.50,
    '[NII]6584': 6585.27,
    '[SII]6718': 6718.29,
    '[SII]6733': 6732.67,
    '[NeIII]3869': 3869.86,
}

# ============================================================
# STACK GRID
# ============================================================
GRID_MIN = 4800.0
GRID_MAX = 6700.0
GRID_N = 1900

# ============================================================
# CRITERIA
# ============================================================
SNR_MIN_TIER_A = 10.0
SNR_MIN_TIER_B = 5.0
Z_MIN = 5.0
Z_MAX = 8.0

# Observed range of the x1d (JADES DR3 G395M/H)
X1D_WAVE_MIN = 28700.0
X1D_WAVE_MAX = 52687.0

# ============================================================
# BOOTSTRAP
# ============================================================
N_BOOTSTRAP = 100
N_BOOTSTRAP_QUICK = 20

# ============================================================
# OUTPUTS
# ============================================================
OUT_QUALITY = ANALYSIS / 'quality_metrics_v2.csv'
OUT_REFINED_Z = ANALYSIS / 'refined_z.csv'
OUT_STACK = ANALYSIS / 'stack_v2.fits'
OUT_FIT = ANALYSIS / 'fit_results_v2.csv'
OUT_BOOTSTRAP = ANALYSIS / 'bootstrap_alpha_v2.csv'
OUT_SUBSAMPLES = ANALYSIS / 'subsamples_v2.csv'
OUT_SYSTEMATICS = ANALYSIS / 'systematics_v2.csv'
OUT_FORECAST = ANALYSIS / 'forecast_v2.csv'
OUT_REPORT = ANALYSIS / 'report_v2.md'

# Bootstrap in batches (Termux memory-safe)
N_BATCH = 10
