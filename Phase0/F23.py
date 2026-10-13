python3 << 'EOF'
import os
import re
import glob
from collections import Counter

# Directories to inspect
dirs = [
    os.path.expanduser('~/todos_los_fits'),
    os.path.expanduser('~/lrd_workspace/spectra'),
    os.path.expanduser('~/lrd_workspace/spectra_lrd_candidates'),
]

all_ids = []

for d in dirs:
    if not os.path.exists(d):
        continue
    files = glob.glob(f'{d}/*.fits')
    print(f"\n{d}: {len(files)} files")
    
    ids_in_dir = []
    for f in files:
        # Extract 8-digit ID
        m = re.search(r'-(\d{8})_', os.path.basename(f))
        if m:
            ids_in_dir.append(m.group(1))
    
    # Count duplicates within the directory
    counter = Counter(ids_in_dir)
    dups = {k: v for k, v in counter.items() if v > 1}
    print(f"  Unique IDs: {len(set(ids_in_dir))}")
    print(f"  Duplicate IDs: {len(dups)}")
    if dups:
        for k, v in list(dups.items())[:5]:
            print(f"    {k}: {v} files")
    
    all_ids.extend(ids_in_dir)

# Duplicates across directories
print(f"\n{'='*60}")
print("GLOBAL SUMMARY")
print('='*60)
print(f"Total files: {len(all_ids)}")
print(f"Global unique IDs: {len(set(all_ids))}")

# IDs that are in more than one directory
counter_global = Counter(all_ids)
dups_global = {k: v for k, v in counter_global.items() if v > 1}
print(f"IDs in multiple directories: {len(dups_global)}")
if dups_global:
    print(f"\nExamples:")
    for k, v in list(dups_global.items())[:5]:
        print(f"  ID {k}: appears {v} times")

# ID of the candidate you downloaded
print(f"\n=== Specific check: 00009883 ===")
for d in dirs:
    files = glob.glob(f'{d}/*00009883*')
    if files:
        print(f"  {d}:")
        for f in files:
            print(f"    {os.path.basename(f)}")
EOF