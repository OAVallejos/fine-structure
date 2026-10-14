# Fine-Structure Constant at z > 5 from JWST/NIRSpec Stacked Spectra

Pipeline for measuring α at z > 5 using rest-frame stacked spectra of 212 JADES DR3 galaxies (JWST/NIRSpec G395M/G395H). We combine Hβ λ4862.68 Å and [OIII] λ5008.24 Å to constrain Δα/α, obtaining Δα/α = +3.02 ± 5.73 ppm (bootstrap, 100 iterations), consistent with no variation at 0.53σ, with a 2σ upper limit |Δα/α| < 1.15 × 10⁻⁵. The dominant systematic is per-object scatter (272 ppm) plus a z_Spec bias (−0.000626), corrected via multi-line redshift (Hβ, Hα, [OIII], [NII]). This is the first upper limit on Δα/α at z > 5 using emission lines. The repo contains the full pipeline (25+ modular Python scripts under `pipeline/v2/`, from quality cuts and stacking to bootstrap, subsamples, systematics and figures), the JADES DR3 and LRD candidate inputs, and the LaTeX preprint source. Built with numpy, scipy, astropy, matplotlib, pandas, scikit-learn, catboost and shap. To reproduce: prepare data, run the pipeline in order, then compile twice with pdflatex. Contact: arielvallejosok@gmail.com —


# Create folders 

mkdir lrd_workspace

cd lrd_workspace

mkdir analysis

mkdir pipeline/v2/ # Here are the *.py files 

mkdir spectra # data

mkdir spectra_lrd_candidates # data2


# Critical Phase 2:

F0.py

Phase1.py

extract_header_features.py


~/lrd_workspace/pipeline/v2

config.py

P01_load.py

P02_quality.py  
             
P02b_quality_lrd.py

P03_refine_z.py 
             
P03b_crosscorr.py

P04_stack.py

P13_multiline_z.py

P04b_stack_multiline.py  

VERIFY_P05.py 
   
P05_fit.py

P05c_fit_multiline.py   
     
P06_bootstrap.py

P06b_bootstrap_multiline.py

P07_subsamples.py

P08_systematics.py

P09_forecast.py

P10_figures.py

P14_validation.py

P15_per_object.py

P16_hb_profiles.py

P17_filter_absorption.py

P20_lrd_profiles.py

P21_lrd_corrected.py

P23_lrd_debug.py

P24_eda_lrd.py

P26_plot_lrd_examples.py

Vallejos, O. A. (2026). Constraint on the Fine Structure Constant at z > 5 from Stacked JWST/NIRSpec Spectra. Zenodo. https://doi.org/10.5281/zenodo.23373471

