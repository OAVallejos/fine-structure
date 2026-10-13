python3 << 'EOF'
import pandas as pd
import os

BASE = os.path.expanduser('~/lrd_workspace')
csv_path = os.path.join(BASE, 'analysis', 'quality_metrics.csv')

df = pd.read_csv(csv_path)
print(f"Total spectra in CSV: {len(df)}")

# Assign tier
df['tier'] = 'C'
df.loc[df['snr_peak'] > 5, 'tier'] = 'B'
df.loc[df['snr_peak'] > 10, 'tier'] = 'A'

print("\n📊 Tier Distribution:")
print(df['tier'].value_counts().sort_index())
print()

# Save lists
for tier in ['A', 'B', 'C']:
    files = df[df['tier'] == tier]['file'].tolist()
    out = os.path.join(BASE, 'analysis', f'tier_{tier}_files.txt')
    with open(out, 'w') as f:
        for x in files:
            f.write(f"{x}\n")
    print(f"✅ Tier {tier}: {len(files)} files → {out}")
EOF