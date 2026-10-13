import numpy as np
import pandas as pd
import re
from collections import Counter

tier_A = pd.read_csv('lrd_workspace/analysis/tier_A_files.txt', header=None, names=['file'])
print(f"Tier A: {len(tier_A)} archivos")

def get_id(f):
    m = re.search(r'(\d{8})', f)
    return m.group(1) if m else None

tier_A['ID'] = tier_A['file'].apply(get_id)
tier_A['field'] = tier_A['file'].apply(lambda x: 'goods-s' if 'goods-s' in x else 'goods-n')
tier_A['grating'] = tier_A['file'].apply(
    lambda x: 'G395M' if 'g395m' in x.lower() else ('G395H' if 'g395h' in x.lower() else 'other'))

# ID ya es string con 8 dígitos (viene del regex)
id_counts = Counter(tier_A['ID'])
duplicated = {k: v for k, v in id_counts.items() if v > 1}

print(f"\nIDs únicos: {len(id_counts)}")
print(f"IDs duplicados: {len(duplicated)}")

print(f"\nField:")
print(tier_A['field'].value_counts())

cat = pd.read_csv('lrd_workspace/analysis/master_features_v4.csv')
cat['NIRSpec_ID_str'] = cat['NIRSpec_ID'].astype(str).str.zfill(8)
z_map = dict(zip(cat['NIRSpec_ID_str'], cat['z_Spec']))

tier_A['z'] = tier_A['ID'].map(z_map)
print(f"\nCon z válido: {tier_A['z'].notna().sum()}/{len(tier_A)}")
print(f"z > 5: {(tier_A['z'] > 5).sum()}")
print(f"z < 5: {(tier_A['z'] < 5).sum()}")

unique_z5 = tier_A[tier_A['z'] > 5].drop_duplicates('ID').copy()
# ¡CLAVE! Guardar ID como string con zfill
unique_z5['ID'] = unique_z5['ID'].astype(str).str.zfill(8)

print(f"\nObjetos únicos con z > 5: {len(unique_z5)}")
print(f"\nID ejemplo: {repr(unique_z5['ID'].iloc[0])}")

unique_z5.to_csv('lrd_workspace/analysis/tier_A_unique_z5.csv', index=False)
print(f"\nGuardado: tier_A_unique_z5.csv")
