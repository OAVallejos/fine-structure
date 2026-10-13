"""
P09: Forecast with systematics.
"""

import numpy as np
import pandas as pd
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from config import *


def main():
    print("=" * 70)
    print("P09: FORECAST WITH SYSTEMATICS")
    print("=" * 70)

    # Load results
    boot = pd.read_csv(OUT_BOOTSTRAP)
    sys_df = pd.read_csv(OUT_SYSTEMATICS)

    dalpha_boot = boot['dalpha_over_alpha'].values
    boot_std = dalpha_boot.std()
    boot_mean = dalpha_boot.mean()

    print("\n[1] Statistical sensitivity (bootstrap):")
    print(f"  σ(Δα/α) = {boot_std:.3e}  ({boot_std*1e6:.2f} ppm)")
    print(f"  Δα/α = {boot_mean*1e6:+.2f} ± {boot_std*1e6:.2f} ppm")

    # z_Spec contribution
    dz_std = sys_df['dz_std'].iloc[0]
    n_objects = 206
    # In the stack, the z error averages out
    dz_effective = dz_std / np.sqrt(n_objects)
    # Conversion to Δα/α: Δλ/λ = -2 Δα/α → Δα/α = -Δλ/(2λ)
    # Δλ/λ = Δz (by definition of redshift)
    # So Δα/α ~ Δz / 2
    dalpha_from_z = dz_effective / 2
    print(f"\n[2] z_Spec contribution:")
    print(f"  σ(Δz) per object: {dz_std:.3e}")
    print(f"  Effective σ(Δz) in stack: {dz_effective:.3e}")
    print(f"  Contribution to σ(Δα/α): {dalpha_from_z:.3e}  ({dalpha_from_z*1e6:.4f} ppm)")

    # Calibration contribution
    # Residual offset of Hβ vs [OIII] 5007
    fit = pd.read_csv(OUT_FIT)
    d_hb = fit[fit['line'] == 'Hβ']['delta_lambda'].values[0]
    e_hb = fit[fit['line'] == 'Hβ']['err'].values[0]
    d_5007 = fit[fit['line'] == '[OIII]5007']['delta_lambda'].values[0]
    e_5007 = fit[fit['line'] == '[OIII]5007']['err'].values[0]
    diff = d_hb - d_5007
    e_diff = np.sqrt(e_hb**2 + e_5007**2)
    dalpha_cal = e_diff / (2 * LAMBDA_VAC['Hβ'])
    print(f"\n[3] Calibration contribution (formal):")
    print(f"  Hβ vs [OIII] 5007 difference: {diff:+.4f} ± {e_diff:.4f} Å")
    print(f"  Formal σ(Δα/α): {dalpha_cal:.3e}  ({dalpha_cal*1e6:.4f} ppm)")

    # Blending contribution
    print(f"\n[4] Blending contribution ([OIII] 4959):")
    ratio = sys_df['ratio_5007_4959'].iloc[0]
    print(f"  Ratio 5007/4959 = {ratio:.3f} (theoretical 2.98)")
    print(f"  Blending: 4959 is not used in the α test")
    print(f"  Impact on σ(Δα/α): negligible (4959 excluded)")

    # Total sensitivity
    print(f"\n[5] Total sensitivity:")
    print(f"  Bootstrap (real):         σ = {boot_std*1e6:.2f} ppm")
    print(f"  Formal fit (underestimates): σ = {dalpha_cal*1e6:.2f} ppm")
    print(f"  Bootstrap/formal ratio:   {boot_std/dalpha_cal:.2f}")

    # 2σ limit
    limit_2sigma = 2 * boot_std
    print(f"\n[6] 2σ limit:")
    print(f"  |Δα/α| < {limit_2sigma:.3e}")
    print(f"  |Δα/α| < {limit_2sigma*1e6:.2f} ppm")

    # Save
    results = {
        'boot_mean': boot_mean,
        'boot_std': boot_std,
        'boot_sigma_ppm': boot_std * 1e6,
        'dalpha_from_z': dalpha_from_z,
        'dalpha_cal': dalpha_cal,
        'limit_2sigma': limit_2sigma,
        'n_objects': n_objects,
    }
    pd.DataFrame([results]).to_csv(OUT_FORECAST, index=False)
    print(f"\nSaved: {OUT_FORECAST}")


if __name__ == '__main__':
    main()
