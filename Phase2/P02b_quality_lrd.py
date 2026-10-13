"""
P02b: Quality control sobre LRDs.
Usa los catálogos de gratings para los z (no master_features).
"""
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table
from pathlib import Path
import re, sys, glob
sys.path.insert(0, str(Path(__file__).parent))
from config import *

SPECTRA_LRD = ROOT / 'spectra_lrd_candidates'


def build_z_map():
    """Construye z_map desde los catálogos de gratings."""
    z_map = {}
    catalogs = glob.glob(str(CATALOGS_JADES / '*gratings*.fits'))
    print(f"Catálogos: {len(catalogs)}")
    for c in catalogs:
        try:
            t = Table.read(c).to_pandas()
            print(f"  {Path(c).name}: {len(t)} filas")
            if 'NIRSpec_ID' in t.columns and 'z_Spec' in t.columns:
                t['ID'] = t['NIRSpec_ID'].astype(str).str.zfill(8)
                for _, row in t.iterrows():
                    z = row['z_Spec']
                    if pd.notna(z) and z > 0:
                        z_map[row['ID']] = z
        except Exception as e:
            print(f"  Error: {e}")
    return z_map


def peak_snr(fl, er):
    if len(fl) < 10: return 0.0
    snr = (fl - np.median(fl)) / er
    return float(np.max(snr))


def main():
    print("=" * 70)
    print("P02b: QUALITY CONTROL — LRDs")
    print("=" * 70)

    z_map = build_z_map()
    print(f"\nz_map: {len(z_map)} objetos")

    files = sorted(SPECTRA_LRD.glob('*.fits'))
    print(f"Archivos en {SPECTRA_LRD}: {len(files)}")

    results = []
    for f in files:
        m = re.search(r'(\d{8})', f.name)
        if not m: continue
        sid = m.group(1)
        z = z_map.get(sid)

        try:
            with fits.open(f, memmap=False) as hdul:
                d = hdul[1].data
                cn = hdul[1].columns.names
                wc = next((c for c in cn if 'wave' in c.lower()), cn[0])
                fc = next((c for c in cn if 'flux' in c.lower() and 'err' not in c.lower()), cn[1])
                ec = next((c for c in cn if 'err' in c.lower()), cn[2])
                w = d[wc].copy(); fl = d[fc].copy(); er = d[ec].copy()
            if np.nanmax(w) < 100: w = w * 1e4
            good = np.isfinite(fl) & np.isfinite(er) & (er > 0)
            if good.sum() < 50: continue
            fl = fl[good]; er = er[good]
        except Exception:
            continue

        snr = peak_snr(fl, er)
        results.append({
            'file': f.name, 'ID': sid, 'z': z, 'snr_peak': snr,
            'grating': 'G395H' if 'g395h' in f.name.lower() else 'G395M',
            'field': 'goods-s' if 'goods-s' in f.name else 'goods-n',
        })

    df = pd.DataFrame(results)
    print(f"\nTotal archivos: {len(df)}")
    print(f"IDs únicos: {df['ID'].nunique()}")
    print(f"Con z válido: {df['z'].notna().sum()}")
    if df['z'].notna().sum() > 0:
        print(f"z > 4: {(df['z'] > 4).sum()}")
        print(f"z > 5: {(df['z'] > 5).sum()}")
        print(f"z rango: {df['z'].min():.3f} – {df['z'].max():.3f}")
        print(f"z mediana: {df['z'].median():.3f}")
    print(f"Tier A (S/N > 10): {(df['snr_peak'] > 10).sum()}")
    print()
    print(f"Field: {df['field'].value_counts().to_dict()}")
    print(f"Grating: {df['grating'].value_counts().to_dict()}")

    df.to_csv(ANALYSIS / 'quality_metrics_lrd.csv', index=False)
    print(f"\nGuardado: quality_metrics_lrd.csv")


if __name__ == '__main__':
    main()
