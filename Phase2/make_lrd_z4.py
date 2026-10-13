"""
Emula el script perdido que generaba analysis/lrd_z4.csv
a partir de quality_metrics_lrd.csv (salida de P02b).

Columnas de salida: ID, file, z, snr_peak
"""
from pathlib import Path
import sys
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from config import ANALYSIS, ROOT

SPECTRA_LRD = ROOT / 'spectra_lrd_candidates'
SRC = ANALYSIS / 'quality_metrics_lrd.csv'      # salida de P02b
OUT = ANALYSIS / 'lrd_z4.csv'

# --- 1) Cargar la fuente ---
if not SRC.exists():
    raise SystemExit(f"No existe {SRC}. Ejecuta antes P02b_quality_lrd.py")

df = pd.read_csv(SRC)
print("Columnas fuente:", df.columns.tolist())

# --- 2) Normalizar columnas ---
df['ID'] = df['ID'].astype(str).str.zfill(8)
df['file'] = df['file'].astype(str)

# --- 3) Filtrar LRDs con S/N > 10 ---
df = df[df['snr_peak'] > 10].copy()
print(f"LRDs tras snr_peak > 10: {len(df)}")

# --- 4) Verificar que los FITS existen ---
missing = [f for f in df['file'] if not (SPECTRA_LRD / f).exists()]
if missing:
    print(f"AVISO: {len(missing)} FITS no encontrados, p.ej.:")
    for m in missing[:5]:
        print("   ", m)

# --- 5) Guardar con las 4 columnas que P19c espera ---
out = df[['ID', 'file', 'z', 'snr_peak']].copy()
out.to_csv(OUT, index=False)
print(f"Escrito {OUT} con {len(out)} filas")
print(out.head())
