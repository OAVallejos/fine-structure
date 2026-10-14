UNDER CONSTRUCTION 

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

P13_multiline_z.py

P14_validation.py

P15_per_object.py

P16_hb_profiles.py

P17_filter_absorption.py

P20_lrd_profiles.py

P21_lrd_corrected.py

P23_lrd_debug.py

P24_eda_lrd.py

P26_plot_lrd_examples.py


