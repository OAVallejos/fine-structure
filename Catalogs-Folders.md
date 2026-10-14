## 6 FITS files from JADES DR3:
 `hlsp_jades_jwst_nircam_goods-n_photometry_v1.0_catalog.fits` | NIRCam GOODS-N photometry (818 MB) |
| `hlsp_jades_jwst_nircam_goods-s-deep_photometry_v2.0_catalog.fits` | NIRCam GOODS-S deep photometry (674 MB) |
| `hlsp_jades_jwst_nirspec_goods-n_gratings-line-fluxes_v1.1_catalog.fits` | GOODS-N grating emission lines (1561 objects) |
| `hlsp_jades_jwst_nirspec_goods-n_prism-line-fluxes_v1.1_catalog.fits` | GOODS-N PRISM emission lines |
| `hlsp_jades_jwst_nirspec_goods-s_gratings-line-fluxes_v1.1_catalog.fits` | GOODS-S grating emission lines (2525 objects) |
| `hlsp_jades_jwst_nirspec_goods-s_prism-line-fluxes_v1.1_catalog.fits` | GOODS-S PRISM emission lines |



---

## Stack
- 212 galaxies at z ∈ [5, 8] from JADES DR3
- Grid: 4800–6700 Å rest-frame, 1900 points
- Columns: `WAVE_REST`, `FLUX`, `FLUX_ERR`
- Detected lines: Hβ, [OIII] 4959, [OIII] 5007, Hα, [NII] 6584
- z used: multi-line z (P13)


---

## LRDs

`lrd_workspace/spectra_lrd_candidates/` → 73 x1d files (62 unique IDs)

LRDs
1. Photometric
- Color: `F444W − F277W > 0.5` (red)
- Magnitude: `20 < mag_F444W < 29` (faint)
- Compact: `F444W_RHALF < 0.15"` (effective radius < 0.15 arcsec)

2. Spectroscopic:
- z > 4 (high-z candidates)
- Strong emission lines (Hβ, [OIII], Hα)
- V-shape at the Balmer break (3645 Å rest-frame) — defining feature (Setton+24)

3. Morphological:
- Compact sources
- Red in NIRCam
- No UV/optical counterpart

Origin: `catalogs/lrd_candidates_z4_v2.csv` (67 initial candidates, 62 with spectra).

Distribution:
- z range: 2.58 – 6.65
- z median: 4.41
- 57 with S/N > 10
- 9 with Hβ in range (z > 4.9)

---

## Visual summary

catalogs_jades/          → 6 JADES DR3 FITS catalogs
lrd_workspace/
├── analysis/ 
├── catalogs/            → lrd_candidates_z4_v2.csv (62 LRDs)
├── pipeline/v2/ 
├── spectra/             → 673 x1d (z > 5)
└── spectra_lrd_candidates/ → 73 x1d (LRDs, z ~ 4.4)

Main stack: (212 galaxies, z ∈ [5,8])

LRDs: `spectra_lrd_candidates/` (62 objects, red color + compact + z > 4 + V-shape)


The main stack uses `spectra/` (673 x1d from z > 5).

The LRDs are in `spectra_lrd_candidates/` (73 x1d). They are a separate folder.

Verification:

```bash
ls ~/lrd_workspace/spectra/ | grep -c "00000993\|00000997\|00001026"
# → 0
```

The LRDs are NOT in `spectra/`. They are disjoint samples.

Main stack: 212 galaxies from `spectra/` (z > 5, Hβ in range)
LRD stack: 57 LRDs from `spectra_lrd_candidates/` (z ~ 4.4, Hα in range)

They do not overlap.