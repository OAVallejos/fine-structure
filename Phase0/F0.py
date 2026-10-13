from astroquery.mast import Observations
import numpy as np
import os
import shutil

print('Consultando MAST...')
obs = Observations.query_criteria(
    provenance_name='jades',
    instrument_name='NIRSpec'
)
products = Observations.get_product_list(obs)

fnames = np.array(products['productFilename'])
mask = np.array([
    ('g395m' in str(f).lower() or 'g395h' in str(f).lower()) and 'x1d' in str(f).lower()
    for f in fnames
])

filtered = products[mask]
print(f'Descargando {len(filtered)} espectros x1d...')

manifest = Observations.download_products(
    filtered,
    download_dir='./temp_download',
    mrp_only=False
)

# Aplanar
target_dir = './todos_los_fits'
os.makedirs(target_dir, exist_ok=True)

count = 0
for root, dirs, files in os.walk('./temp_download'):
    for file in files:
        if file.endswith('.fits'):
            src_path = os.path.join(root, file)
            dst_path = os.path.join(target_dir, file)
            if os.path.exists(dst_path):
                base, ext = os.path.splitext(file)
                dst_path = os.path.join(target_dir, f'{base}_dup{count}{ext}')
            shutil.copy2(src_path, dst_path)
            os.remove(src_path)
            count += 1
            if count % 100 == 0:
                print(f'  Movidos: {count}')

shutil.rmtree('./temp_download', ignore_errors=True)
print(f'¡Listo! {count} archivos en ./todos_los_fits')
