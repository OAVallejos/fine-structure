python3 << 'EOF'
import os
import re
import pandas as pd
from collections import Counter

target = os.path.expanduser('~/lrd_workspace/spectra_lrd_candidates')
files = [f for f in os.listdir(target) if f.endswith('.fits')]

ids = []
for f in files:
    m = re.search(r'-(\d{8})_', f)
    if m:
        ids.append(m.group(1))

counter = Counter(ids)
print(f"Files: {len(files)}")
print(f"Unique IDs: {len(set(ids))}")
print(f"IDs with multiple files: {len([k for k,v in counter.items() if v>1])}")

# Candidates from CSV
df = pd.read_csv(os.path.expanduser('~/lrd_workspace/catalogs/lrd_candidates_z4_v2.csv'))
candidates = set(df[df['is_LRD_candidate'] == True]['NIRSpec_ID'].astype(str).str.zfill(8))
present = set(ids)
missing = candidates - present

print(f"\nTotal candidates: {len(candidates)}")
print(f"Present: {len(present)}")
print(f"Missing: {len(missing)}")
if missing:
    print(f"\nMissing IDs:")
    for id_ in sorted(missing):
        print(f"  {id_}")

# Search those IDs in ALL folders
print(f"\n{'='*60}")
print("SEARCHING FOR MISSING IDS IN ALL FOLDERS")
print('='*60)
dirs = [
    os.path.expanduser('~/todos_los_fits'),
    os.path.expanduser('~/lrd_workspace/spectra'),
    os.path.expanduser('~/lrd_workspace/spectra_lrd_candidates'),
]
import glob
for id_ in sorted(missing):
    found = False
    for d in dirs:
        files = glob.glob(f'{d}/*{id_}*')
        if files:
            print(f"\n{id_}:")
            for f in files:
                print(f"  {os.path.basename(f)}")
            found = True
    if not found:
        print(f"\n{id_}: ❌ Not in any folder")
EOF