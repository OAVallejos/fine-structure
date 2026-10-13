python3 << 'EOF'
from astroquery.mast import Observations
import numpy as np

# Missing IDs
missing_ids = ['00005321', '00006519', '00009325', '00009901', 
               '00107028', '00000721', '00002845', '00008239',
               '00009429', '00009569', '00009883']  # add more if they appear

obs = Observations.query_criteria(provenance_name='jades', instrument_name='NIRSpec')
products = Observations.get_product_list(obs)

fnames = np.array(products['productFilename'])
mask = np.array([
    any(mid in f for mid in missing_ids) and ('x1d' in f.lower())
    for f in fnames
])
filtered = products[mask]
print(f"Products found: {len(filtered)}")
for f in filtered['productFilename']:
    print(f"  {f}")

if len(filtered) > 0:
    Observations.download_products(
        filtered,
        download_dir='./lrd_workspace/spectra_lrd_candidates',
        mrp_only=False
    )
    print("✅ Download complete")
EOF