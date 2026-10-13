python3 << 'EOF'
import pandas as pd
import os

# Make sure this path is correct according to where you saved the CSV
csv_path = os.path.expanduser('~/lrd_workspace/analysis/lrd_classification.csv')

if not os.path.exists(csv_path):
    print("⚠️ Cannot find lrd_classification.csv. Did we generate it earlier?")
else:
    df = pd.read_csv(csv_path)
    print("Evaluating color (F444W - F277W) and magnitude thresholds:")
    for thresh in [0.5, 0.8, 1.0, 1.2, 1.5]:
        mask = (
            (df['color_F444W_F277W'] > thresh) &
            (df['mag_F444W'].between(20, 30)) &
            (df['mag_F444W'].notna())
        )
        print(f"  Color > {thresh}: {mask.sum()} sources")
EOF