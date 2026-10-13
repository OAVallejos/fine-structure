python3 << 'EOF'
from astropy.table import Table
import os
import re

BASE = os.path.expanduser('~/lrd_workspace')
cat_s = os.path.join(os.path.expanduser('~'), 'catalogs_jades', 'hlsp_jades_jwst_nirspec_goods-s_gratings-line-fluxes_v1.1_catalog.fits')
t = Table.read(cat_s, hdu=1)

# IDs of the 3
ids = ['201249', '202208', '210003']
peaks = [30190, 32291, 33943]  # Å

for id_, lam_obs in zip(ids, peaks):
    # Search in the catalog
    mask = [str(x).zfill(8) == id_.zfill(8) or str(x) == id_ for x in t['NIRSpec_ID']]
    if any(mask):
        row = t[mask][0]
        z = row['z_Spec']
        print(f"\n📄 ID {id_}:")
        print(f"   z_Spec = {z:.4f}")
        print(f"   λ_obs = {lam_obs} Å")
        # Try candidate lines
        for line, lam_rest in [('[OIII]5007', 5007), ('Hβ 4861', 4861), ('Hα 6563', 6563), ('[OII]3727', 3727)]:
            z_implied = lam_obs / lam_rest - 1
            if abs(z_implied - z) < 0.05:
                print(f"   ✅ {line} (λ_rest={lam_rest}) → z_implied={z_implied:.3f}")
EOF