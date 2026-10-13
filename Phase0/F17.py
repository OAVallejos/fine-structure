python3 << 'EOF'
from astropy.io import fits
import os

for name, path in [
    ('GOODS-N', '~/catalogs_jades/hlsp_jades_jwst_nircam_goods-n_photometry_v1.0_catalog.fits'),
    ('GOODS-S', '~/catalogs_jades/hlsp_jades_jwst_nircam_goods-s-deep_photometry_v2.0_catalog.fits'),
]:
    full = os.path.expanduser(path)
    print(f"\n{'='*70}")
    print(f"{name}: {os.path.basename(full)}")
    print('='*70)
    with fits.open(full) as h:
        for i, ext in enumerate(h):
            if ext.data is not None and hasattr(ext.data, 'names'):
                print(f"\nExt {i} ({ext.name}): {len(ext.data)} rows")
                print(f"  Total columns: {len(ext.data.names)}")
                
                # Search for F444W, F277W, RHALF
                f444 = [c for c in ext.data.names if 'F444W' in c.upper()]
                f277 = [c for c in ext.data.names if 'F277W' in c.upper()]
                rhal = [c for c in ext.data.names if 'RHALF' in c.upper()]
                
                if f444:
                    print(f"  → F444W: {f444[:8]}")
                if f277:
                    print(f"  → F277W: {f277[:8]}")
                if rhal:
                    print(f"  → RHALF: {rhal[:8]}")
EOF