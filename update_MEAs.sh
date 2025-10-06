#!/bin/bash

# Node start time, parent directory, and custom preferences for rereferencing,
# denoising, clustering

ROOT_DIR="/path/to/micros_pipeline/parent_directory"   # Parent directory containing 'MEA_pipeline' and 'MEA_database' folders
DO_COMBINATO=1                                         # (1/0) Whether to detect and cluster MEA spike data with Combinato
DO_COMBINATO_MAXWELL=1                                 # (1/0) Whether to cluster spikes detected by MEA with Combinato
DO_KILOSORT=1                                          # (1/0) Whether to detect and cluster MEA spike data with Kilosort
DO_PLOTS=0                                             # (1/0) Whether to generate plots of processing steps

cd $ROOT_DIR
cd MEA_pipeline
module load matlab

# Run MATLAB code without display
matlab -nodesktop -nodisplay -singleCompThread -r "update_MEA('$ROOT_DIR', $DO_COMBINATO, $DO_COMBINATO_MAXWELL, $DO_KILOSORT, $DO_PLOTS); exit;"