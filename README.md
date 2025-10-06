# MEA_pipeline
**This pipeline processes microelectrode recordings from MaxWell Biosystems to detect and cluster single unit neural activity.**

## Overview:
This project provides a Matlab-based pipeline that facilitates the processing and analysis of Maxwell Biosystems microelectrode array recordings using spike clustering methods [combinato](https://github.com/jniediek) and [Kilosort 3](https://github.com/MouseLand/Kilosort). It is tailored for cognitive electrophysiologists working on microelectrode data research. The data produced can readily be used for cognitive electrophysiology analyses using other Matlab-based tools developed by Texas Computational Memory Lab.


## Features:
- **Data Restructuring**: MaxWell Biosystems data is restructured so that it can be more easily manipulated and plotted using matrices.
- **Codebase Integration**: Pipeline expands analytic repertoire by making data compatible with Combinato and Kilosort. 
- **Spike Clustering Quality Metrics**: Spike clustering quality metrics in matrix format could aid in isolating distinct neurons within the microelectrode array. 
- **GPU Processing**: Faster processing of data from thousands of individual recording channels.


## Tech Stack:
- **Platform**: Linux.
- **Requirements**: Anaconda. GPU node.
- **Languages**: Matlab.
- **Key Packages**: Must install [combinato](https://github.com/jniediek). Must install [Kilosort](https://github.com/MouseLand/Kilosort).


## Installation:
1) Download repository into desired parent directory.
2) Install [combinato](https://github.com/jniediek) in MEA_pipeline/base_code folder. Anaconda environment for combinato is found in the same folder.
3) Install [Kilosort](https://github.com/MouseLand/Kilosort) in MEA_pipeline/base_code folder. Kilosort functions datashift and rezToPhy need to be modified to adapt to MEA data.
4) Edit paths for parent directory in all files.


## Usage:
1) Individual MaxWell Biosystems raw microelectrode array recordings (.raw.h5) should be stored in folders with the following format: /parent_directory/MEA_database/subject_code/yyyy-mm-dd_scan-name#/raw/MEA_recording.raw.h5

2) Execute MEA_pipeline/update_MEA.sh in terminal.


## Optional:
Install [Phy](https://github.com/cortex-lab/phy) to visualize and modify Kilosort clustering results.


## Future Improvements:
- Improve multirecording data processing and device methods for utilizing Combinato results for isolating distinct neurons in MEA. 
- Device methods for utilizing Kilosort results for extracting and refining spike time series data.
- Add more methods for referencing, denoising, and clustering.
