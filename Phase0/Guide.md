### Phase 1 — `quality_metrics.py` (Fase1)
Purpose: Quality control of all x1d spectra.
- Scans `~/lrd_workspace/spectra/` for all `.fits` files.
- For each spectrum: cleans NaN/invalid pixels, trims edges, counts valid pixels, and computes the peak S/N.
- Outputs:
  - `quality_metrics.csv` — table with `file`, `n_valid_pixels`, `snr_peak`.
  - `data_quality_check.png` — histograms of valid pixels and peak S/N.
- Prints a summary (how many spectra pass S/N > 3, 5, 10) and the Top 20 by peak S/N.

---

### Phase 2 — `peak_inspection.py` (Fase2)
Purpose: Manual inspection of the 3 best spectra.
- Reads `quality_metrics.csv` and takes the Top 3 by `snr_peak`.
- For each one, finds the S/N peak pixel and estimates the line width (approximate FWHM in pixels and Å).
- Prints peak S/N, observed wavelength, peak flux, and peak error.
- Used to eyeball whether the strongest detections are real emission lines.

---

### Phase 3 — `line_identification.py` (Fase3)
Purpose: Identify which emission line produces the observed peak.
- Reads the NIRSpec GOODS-S grism catalog.
- For the 3 IDs from Phase 2 (`201249`, `202208`, `210003`) with their observed peak wavelengths, it looks up `z_Spec`.
- Tests candidate rest-frame lines (Hβ, [OIII] 4959/5007, Hα, [OII] 3727) and checks which one gives an implied redshift consistent with the catalog value (Δz < 0.05).
- Prints which line matches each source.

---

### Phase 4 — `stacking_restframe.py` (Fase4)
Purpose: Stack all Tier A spectra in the rest frame (the "Killer Figure").
- Loads `tier_A_files.txt` (highest S/N subsample).
- Builds a `{ID → z_Spec}` map from both GOODS-S and GOODS-N grism catalogs.
- For each Tier A spectrum: converts to rest-frame wavelengths, interpolates onto a common grid (4750–5150 Å), and weights by inverse variance.
- Produces:
  - `stacked_tierA_restframe.fits` — stacked flux, error, and N spectra per pixel.
  - `stacked_lrd_tierA.png` — final stacked spectrum with Hβ, [OIII] 4959, [OIII] 5007 marked.
- Prints how many spectra were stacked, skipped (no z), or failed.

---

### Phase 5 — `tier_assignment.py` (Fase5)
Purpose: Classify all spectra into quality tiers.
- Reads `quality_metrics.csv`.
- Assigns tiers based on peak S/N:
  - Tier A: S/N > 10
  - Tier B: S/N > 5
  - Tier C: the rest
- Writes three lists: `tier_A_files.txt`, `tier_B_files.txt`, `tier_C_files.txt`.
- These lists feed directly into Phase 4 (stacking).

---

### Pipeline flow (one line each):

| Phase | Role |
|-------|------|
| 1 | Measure quality → make `quality_metrics.csv` |
| 2 | Inspect Top 3 peaks in detail |
| 3 | Identify the emission line behind each peak |
| 4 | Stack Tier A spectra in rest frame → science figure |
| 5 | Split all spectra into Tier A/B/C lists |

In short: 1 → 5 (classify) → 4 (stack) is the main path; 2 and 3 are diagnostic steps to validate the detections before trusting the stack.

Here's a brief guide to what each script (Phase 6–10) does, in English:

---

### Phase 6 — `stacking_restframe.py` (Fase6)
Purpose: First clean version of the rest-frame stack of Tier A spectra.
- Loads `tier_A_files.txt` and builds a `{ID → z_Spec}` map from the GOODS-S and GOODS-N grism catalogs.
- Converts each spectrum to rest-frame wavelengths, interpolates onto a common grid (4750–5150 Å), and combines them with inverse-variance weighting.
- Outputs:
  - `stacked_tierA_restframe.fits` — stacked flux + error.
  - `stacked_lrd_tierA.png` — the "Killer Figure" with Hβ, [OIII] 4959, [OIII] 5007 marked.
- Limitation: Doesn't handle unit issues or NaN/`inf` edge cases carefully — this is why the next phases exist.

---

### Phase 7 — `stacking_aligned.py` (Fase7)
Purpose: Second stack attempt with per-spectrum peak alignment.
- Same as Phase 6, but adds a pre-alignment step: detects the [OIII] 5007 peak in each individual spectrum (via smoothed Gaussian detection), and shifts it to the expected rest-frame position before stacking.
- Only shifts if the detected peak is significant (>3σ) and the shift is small (<30 Å).
- Applies 4 fixes:
  1. Masks NaNs before interpolating.
  2. Uses safe `fill_value` (0 for flux, `inf` for error).
  3. Only aligns if the peak is real.
  4. Excludes `inf` values from the weight computation.
- Outputs: `stacked_tierA_aligned.fits` and `stacked_lrd_tierA_aligned.png`.

---

### Phase 8 — `stacking_final.py` (Fase8)
Purpose: Definitive stack with proper unit handling.
- Identical to Phase 7 but adds a robust unit conversion fix:
  - If `max(wave) < 1.0` → assume meters, convert to Å.
  - If `max(wave) < 100` → assume microns, convert to Å (JADES DR3 uses microns).
  - If `max(wave) > 100` → already in Å.
- Uses `np.interp` (faster than `interp1d`) with `left/right` fill values.
- Reports grid coverage (fraction of wavelength grid with valid data).
- Outputs: `stacked_tierA_final.fits` and `stacked_lrd_final.png`.
- This is the science-quality stack — the one that should show [OIII] peaks clearly.

---

### Phase 9 — `forecast_sensitivity.py` (Fase9)
Purpose: Forecast how sensitive the stack is to a varying fine-structure constant (Δα/α).
- Loads the final stack from Phase 8.
- Fits a Gaussian to the real [OIII] 5007 line to get baseline parameters (amplitude, center, σ, continuum).
- Injects a synthetic shift equivalent to Δα/α = 5×10⁻⁵ (via the relation Δλ/λ ≈ −2 Δα/α).
- Adds realistic Gaussian noise using the stack's own error array.
- Refits to recover the injected Δα/α and estimates its uncertainty from the covariance matrix.
- Outputs:
  - `forecast_sensitivity.png` — two-panel plot: top shows real vs mock vs fit; bottom shows residuals in σ.
- Goal: Demonstrate what level of Δα/α the JADES DR3 stack could detect.

---

### Phase 10 — `fit_results.txt` (Fase10)
Purpose: Not a script — a results summary from the multi-line Gaussian fit.
- Reports the fitted parameters for the three lines (Hβ, [OIII] 4959, [OIII] 5007):
  - Continuum (a0, a1), amplitude (A), center (μ), and width (σ) for each line, with uncertainties.
- Computes offsets vs laboratory wavelengths and converts them into Δα/α constraints:
  - Hβ: Δα/α = −1.17×10⁻⁵ ± 1.05×10⁻⁶
  - [OIII] 4959: Δα/α = −1.24×10⁻⁵ ± 7.03×10⁻⁷
  - [OIII] 5007: Δα/α = −6.69×10⁻⁶ ± 3.15×10⁻⁷
- Also reports line ratios ([OIII]5007/[OIII]4959 and [OIII]5007/Hβ) as a sanity check.

---

### Pipeline flow (one line each):

| Phase | Role |
|-------|------|
| 6 | First rest-frame stack (naive) |
| 7 | Stack with per-spectrum peak alignment + fixes |
| 8 | Final stack with robust unit handling → science product |
| 9 | Sensitivity forecast for Δα/α using the final stack |
| 10 | Fit results: line centers, offsets, and Δα/α constraints |

In short: 6 → 7 → 8 iteratively improve the stack (fixing units, alignment, NaN handling), 9 uses the final stack to forecast Δα/α sensitivity, and 10 is the actual measurement of Δα/α from the fitted line centers.

Science narrative: This is a fine-structure constant variation test using the [OIII] 4959/5007 and Hβ lines stacked from JADES DR3 LRD candidates.

---

### Phase 11 — `forecast_sensitivity_curve.py` (Fase11)
Purpose: First-pass sensitivity forecast for Δα/α.
- Loads the final stack from Phase 8.
- Fits the three lines (Hβ, [OIII] 4959, [OIII] 5007) with a Gaussian model to get baseline parameters.
- Injects a range of Δα/α values (1×10⁻⁶ → 1×10⁻³) into the mock spectrum using Δλ/λ = −2 Δα/α.
- Refits each mock and recovers Δα/α per line, then computes a weighted average.
- Writes `forecast_sensitivity_curve.csv` with injected, recovered, and error values.
- Limitation: Uses a perfect model (no systematics), so the errors are optimistic.

---

### Phase 12 — `forecast_conservative.py` (Fase12)
Purpose: Conservative forecast that includes realistic systematics.
- Same idea as Phase 11, but the injection model differs from the fit model:
  - Lines convolved with an extra Gaussian LSF (σ = 1 Å).
  - A z_Spec error (±0.001) that shifts the entire spectrum.
  - A continuum curvature term (a2) that the fit does not model.
- Runs 20 noise trials per injected Δα/α value and reports the mean, bias, and scatter.
- Writes `forecast_conservative.csv` and prints:
  - Statistical floor (~2.8×10⁻⁷).
  - Conservative error (~6.6×10⁻⁵).
- Key takeaway: The real sensitivity lies between these two, closer to the conservative value.

---

### Phase 13 — `offset_test.py` (Fase13)
Purpose: Diagnostic to decide whether the measured offsets are calibration or physics.
- Fits the stack and extracts Δλ for Hβ, [OIII] 4959, [OIII] 5007.
- Runs four tests:
  1. Consistency: Are all three offsets equal? (χ²_red test)
  2. Direction: Are all offsets negative (or positive)?
  3. Hβ vs [OIII]: Is Hβ shifted differently from [OIII]? (α-sensitivity test)
  4. Summary verdict: Calibration vs possible physical signal.
- Writes `offset_test.csv` with per-line Δλ, error, and significance.
- Key point: If all offsets are equal → NIRSpec calibration. If Hβ ≠ [OIII] → potential Δα/α signal.

---

### Phase 14 — `model_comparison.py` (Fase14)
Purpose: Compare four fit models with Δχ² to decide which restrictions are justified.
- Models tested:
  - A) Free (11 params): all lines independent.
  - B) Common widths (9 params): one σ for all lines.
  - C) Common shift (9 params): one shift applied to all three lines.
  - D) Fixed 5007/4959 ratio = 2.98 (10 params).
- Prints χ², χ²_red, and Δχ² vs the free model for each.
- Interpretation rules:
  - Δχ² > 25 (3σ) → restriction does not hold.
  - Δχ² < 6 (1σ) → restriction does hold.
  - χ²_red > 10 in all → the Gaussian model itself is inadequate.
- Key point: If the "common shift" model is as good as the free one, the offsets are consistent with a single systematic shift (calibration).

---

### Phase 15 — `generate_ml_dataset.py` (Fase15)
Purpose: Build a synthetic training set for a machine-learning redshift estimator.
- Generates 5000 mock spectra with:
  - Random z ∈ [5.0, 7.5].
  - Random Hβ, [OIII] 4959, [OIII] 5007 amplitudes (with realistic ratios).
  - Random line widths (σ = 2–4 Å) and continuum levels.
  - Gaussian noise (σ = 0.01) on a G395M-like wavelength grid (29000–53000 Å).
- Normalizes each spectrum (mean-subtracted, unit-variance).
- Saves:
  - `ml_data/spectra.npy` — shape (5000, 1500).
  - `ml_data/z_true.npy` — the true redshifts.
- Goal: Provide a labeled dataset to train a neural network that predicts z directly from the spectrum.

---

### Pipeline flow (one line each):

| Phase | Role |
|-------|------|
| 11 | Simple injection–recovery for Δα/α (idealized) |
| 12 | Conservative injection–recovery (with systematics) |
| 13 | Offset consistency test → calibration vs physics |
| 14 | Model comparison via Δχ² |
| 15 | Synthetic ML dataset for z estimation |

In short:
- 11 → 12 quantify the Δα/α sensitivity (optimistic → realistic).
- 13 → 14 test whether the measured offsets are calibration or a real signal.
- 15 pivots to machine learning: generate mock spectra to train a redshift predictor.

Here's a brief functional guide to Phase 16–20 in English:

---

### Phase 16 — `lrd_classification_v1.py` (Fase16)
Purpose: First attempt at classifying Tier A spectra as Little Red Dots (LRDs) using the GOODS-N photometric catalog.
- Loads the GOODS-N DR3 NIRCam photometric catalog (extensions 3 = SIZE/RHALF, 8 = KRON_CONV, 9 = PHOTOZ).
- Merges by `ID` to get RA, DEC, RHALF, F444W, F277W.
- Auto-detects whether the photometric values are fluxes (nJy) or magnitudes (by inspecting an example value): if > 100, converts fluxes → magnitudes using the AB zero point.
- Computes the color `F444W − F277W`.
- Cross-matches with the Tier A gratings list via `NIRCam_ID`.
- Applies a strict LRD cut: `color > 1.5` and `22 < mag_F444W < 29`.
- Outputs `lrd_classification.csv` with an `is_LRD` flag.
- Limitation: Only uses GOODS-N (no GOODS-S), so it misses half the sample.

---

### Phase 17 — `catalog_explorer.py` (Fase17)
Purpose: Debug tool — inspect the structure of both GOODS-N and GOODS-S photometric catalogs.
- Opens each FITS file and iterates through all extensions.
- For each table extension, prints:
  - Extension index and name.
  - Number of rows and columns.
  - Which columns contain `F444W`, `F277W`, and `RHALF`.
- Key point: Used to figure out the correct HDU indices (3 = SIZE, 8 = KRON_CONV) and exact column names before running the classification scripts.

---

### Phase 18 — `lrd_classification_v2.py` (Fase18)
Purpose: Corrected LRD classification using both GOODS-N and GOODS-S photometry.
- Loads photometry from both fields (extensions 3 and 8), concatenates them.
- Robust flux/magnitude detection using the median of `F444W_KRON` (instead of a single sample value):
  - If `0 < median < 1000` → fluxes (nJy), convert to magnitudes.
  - Otherwise → already magnitudes.
- Loads both gratings catalogs, tags them by field (`GOODS-S`, `GOODS-N`), and filters Tier A.
- Cross-matches by `NIRCam_ID` (converted to int).
- Applies the LRD cut: `color > 1.5`, `20 < mag < 30`.
- Prints:
  - Total Tier A, matched sources, red sources, valid magnitudes, confirmed LRDs.
  - Breakdown by field.
  - Top 10 LRDs by color.
- Outputs updated `lrd_classification.csv`.

---

### Phase 19 — `threshold_scan.py` (Fase19)
Purpose: Tune the color threshold for LRD selection.
- Reads `lrd_classification.csv`.
- Loops over a range of color thresholds: 0.5, 0.8, 1.0, 1.2, 1.5.
- For each threshold, counts how many sources pass `color > threshold` AND `20 < mag_F444W < 30`.
- Key point: Shows how sensitive the LRD count is to the color cut — used to decide whether 0.5 (relaxed) or 1.5 (strict) is more appropriate. This motivates the switch to a relaxed cut in Phase 20.

---

### Phase 20 — `lrd_candidates_z4.py` (Fase20)
Purpose: Final LRD candidate selection for all gratings spectra with z ≥ 4 (not just Tier A).
- Loads both gratings catalogs (GOODS-S + GOODS-N), tags by field, and filters `z_Spec ≥ 4`.
- Loads both photometric catalogs, converts fluxes → magnitudes, computes color.
- Cross-matches by `NIRCam_ID`.
- Applies a relaxed, multi-criteria LRD cut:
  - `color_F444W_F277W > 0.5` (relaxed vs 1.5 in Phase 16/18).
  - `20 < mag_F444W < 29`.
  - Compactness: `F444W_RHALF < 0.15"` (new — LRDs are compact).
- Prints:
  - Total z ≥ 4, matched, color-pass, mag-pass, compact-pass, final LRD candidates.
  - Breakdown by field.
  - z distribution of the candidates.
- Outputs the final catalog: `~/lrd_workspace/catalogs/lrd_candidates_z4.csv` — 67 candidates that feed the rest of the pipeline.

---

### Pipeline flow (one line each):

| Phase | Role |
|-------|------|
| 16 | First LRD cut (GOODS-N only, strict color > 1.5) |
| 17 | Debug: inspect FITS structure of both photometric catalogs |
| 18 | Corrected LRD cut using both fields + auto flux/mag detection |
| 19 | Threshold scan to tune the color cut |
| 20 | Final LRD catalog: all z ≥ 4 gratings + relaxed color + compactness → 67 candidates |

In short:
- 16 → 17 explore the photometric data and get the column structure right.
- 18 → 19 refine the LRD selection criteria (color + magnitude) on Tier A.
- 20 expands to all z ≥ 4 spectra and applies a relaxed + compactness-aware cut, producing the final `lrd_candidates_z4.csv` (67 sources) that drives the rest of the analysis.

Science narrative: This block transitions from spectral stacking (Phases 1–15) to photometric LRD identification. Phase 16–20 build the LRD sample by combining NIRSpec redshifts with NIRCam colors and sizes — the classic "little red dot" selection (red F444W−F277W color, compact morphology, faint magnitude).

Here's a brief functional guide to Phase 21–25 in English:

---

### Phase 21 — `lrd_candidates_v2.py` (Fase21)
Purpose: Corrected LRD candidate selection for all gratings spectra with z ≥ 4.
- Loads both gratings catalogs (GOODS-S + GOODS-N), tags them by field, and loads both photometric catalogs (GOODS-S deep v2.0 + GOODS-N v1.0).
- Key fix vs. Phase 20: merges each field separately by `NIRCam_ID` (avoiding ID collisions that occur when concatenating), then concatenates the two already-merged DataFrames.
- Converts fluxes → magnitudes, computes `color = F444W − F277W`.
- Applies the relaxed LRD cut:
  - `color > 0.5`
  - `20 < mag_F444W < 29`
  - `F444W_RHALF < 0.15"` (compactness)
- Prints total z ≥ 4, matched, color-pass, mag-pass, compact-pass, and final LRD candidates, with breakdown by field.
- Output: `~/lrd_workspace/catalogs/lrd_candidates_z4_v2.csv` — the definitive LRD candidate catalog used downstream.

---

### Phase 22 — `download_missing_spectra.py` (Fase22)
Purpose: Download any missing LRD spectra directly from MAST.
- Hardcodes a list of 11 missing NIRSpec IDs (e.g. `00005321`, `00006519`, `00009883`…).
- Queries MAST with `provenance_name='jades'` and `instrument_name='NIRSpec'`.
- Retrieves the full product list and filters for files where:
  - The filename contains any of the missing IDs.
  - The filename contains `x1d` (case-insensitive).
- Prints the matched products and downloads them into `~/lrd_workspace/spectra_lrd_candidates/`.
- Goal: Fill gaps in the LRD sample so stacking has maximum N.

---

### Phase 23 — `audit_spectra_inventory.py` (Fase23)
Purpose: Inventory and duplicate audit across all spectra directories.
- Scans three directories:
  - `~/todos_los_fits/` (full 3207-spectrum dataset).
  - `~/lrd_workspace/spectra/` (673 spectra at z > 5).
  - `~/lrd_workspace/spectra_lrd_candidates/` (LRD candidate spectra).
- For each directory:
  - Counts files and extracts 8-digit IDs from filenames.
  - Reports unique IDs and duplicate IDs (same ID in multiple files, e.g. G395M + G395H).
- Then does a global count:
  - Total files across all directories.
  - Global unique IDs.
  - IDs appearing in more than one directory (i.e. potential redundancy).
- Finally, checks specifically where ID `00009883` lives (as a sanity check on a candidate that was just downloaded).
- Key point: Diagnostic — confirms which candidates are present, duplicated, or entirely missing before consolidation.

---

### Phase 24 — `consolidate_lrd_spectra.py` (Fase24)
Purpose: Copy all LRD candidate spectra into a single folder.
- Reads `lrd_candidates_z4_v2.csv` and extracts the set of LRD candidate IDs (zfilled to 8 digits).
- Iterates over the three source directories:
  - `~/todos_los_fits/`
  - `~/lrd_workspace/spectra/`
  - `~/lrd_workspace/spectra_lrd_candidates/`
- For every `.fits` file, extracts its 8-digit ID via regex; if that ID is in the candidate set:
  - Copies the file to `~/lrd_workspace/spectra_lrd_candidates/` (using `shutil.copy2`, preserving metadata).
  - Skips files already present in the target (no overwrite).
- Prints how many new files were copied and the total count in `spectra_lrd_candidates/`.
- Goal: Produce a single consolidated folder with every LRD candidate spectrum, ready for stacking.

---

### Phase 25 — `check_lrd_completeness.py` (Fase25)
Purpose: Final completeness check of the LRD candidate sample.
- Lists all `.fits` files in `~/lrd_workspace/spectra_lrd_candidates/`, extracts 8-digit IDs, and reports:
  - Total files, unique IDs, IDs with multiple files (e.g. G395M + G395H).
- Loads the candidate list from `lrd_candidates_z4_v2.csv` and computes:
  - Candidates total.
  - Candidates present in the folder.
  - Candidates missing.
- Prints the missing IDs explicitly.
- Then searches all three spectra directories for each missing ID to see whether it exists somewhere else (and just wasn't copied), or is truly absent.
- Key point: Last QC step before stacking — tells you whether you have 100% coverage of the LRD sample, or whether Phase 22 (MAST download) needs to be run again for specific IDs.

---

### Pipeline flow (one line each):

| Phase | Role |
|-------|------|
| 21 | Final LRD candidate catalog (per-field merge, relaxed cuts) → `lrd_candidates_z4_v2.csv` |
| 22 | Download missing LRD spectra from MAST |
| 23 | Audit inventory + duplicates across spectra directories |
| 24 | Consolidate all LRD candidate spectra into `spectra_lrd_candidates/` |
| 25 | Verify completeness; flag any still-missing candidates |

In short:
- 21 → builds the definitive LRD catalog (67 candidates) with a correct per-field merge.
- 22 → 23 → download missing files and diagnose what exists where.
- 24 → 25 → consolidate all candidate spectra into one folder and verify coverage.

Science narrative: This block moves from catalog selection (Phase 21) to physical data assembly (Phase 22–25). The LRD sample has already been identified photometrically; now these scripts make sure every candidate has its NIRSpec x1d spectrum locally available — downloading missing ones from MAST, removing duplicates, consolidating into one folder, and confirming completeness — so that the final stacking and Δα/α analysis is done on a complete, clean, well-audited sample.
