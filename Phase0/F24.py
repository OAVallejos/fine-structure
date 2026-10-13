python3 << 'EOF'
import os
import shutil
import pandas as pd
import glob

# Load candidates
df = pd.read_csv(os.path.expanduser('~/lrd_workspace/catalogs/lrd_candidates_z4_v2.csv'))
candidates = df[df['is_LRD_candidate'] == True]
candidate_ids = set(candidates['NIRSpec_ID'].astype(str).str.zfill(8))

print(f"LRD candidates: {len(candidate_ids)}")

# Source directories
source_dirs = [
    os.path.expanduser('~/todos_los_fits'),
    os.path.expanduser('~/lrd_workspace/spectra'),
    os.path.expanduser('~/lrd_workspace/spectra_lrd_candidates'),
]

target_dir = os.path.expanduser('~/lrd_workspace/spectra_lrd_candidates')
os.makedirs(target_dir, exist_ok=True)

# Move ALL candidate files to spectra_lrd_candidates/
copied = 0
for src_dir in source_dirs:
    if not os.path.exists(src_dir):
        continue
    for f in glob.glob(f'{src_dir}/*.fits'):
        basename = os.path.basename(f)
        # Extract ID
        import re
        m = re.search(r'-(\d{8})_', basename)
        if not m:
            continue
        id_ = m.group(1)
        if id_ in candidate_ids:
            dst = os.path.join(target_dir, basename)
            if not os.path.exists(dst):
                shutil.copy2(f, dst)
                copied += 1

print(f"Copied: {copied} new files")
print(f"In spectra_lrd_candidates: {len(glob.glob(f'{target_dir}/*.fits'))} files")
EOF