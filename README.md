# HRTF PCA

Performs a principal component analysis (PCA) on a set of head-related transfer functions (HRTFs). It reconstructs 128-pt HRTF's from the UC Davis CIPIC dataset from 1 (128x data reduction) base to 128 bases (no data reduction).

## Instructions
1. Download this repository.
2. Download the CIPIC dataset from the [UC Davis website](https://www.ece.ucdavis.edu/cipic/spatial-sound/hrtf-data/). Navigate to the "Full Database Download" and download the .zip file.
3. Unzip the file into the repository.
4. Make sure the HRTF directory paths point to the "standard_hrir_database/" folder in the CIPIC dataset. This needs to be checked in project_run.m (line 9, input.directory) and hrtf_load.m (line 12, input.dir).
5. Run run_project.m and resolve any runtime fixes that I may have forgotten.

**Note:** This project takes a while to run. It took about 10 minutes in total, of which 7 minutes were spent reconstructing the HRTFs and calculating loss functions, on my personal laptop (Intel i5-5300 4CPU @ 2.30GHz; 12GB RAM). It also spits out five plot figures.

**Note:** There is functionality to load the anthropomorphic data from CIPIC, but it is not used in the project.
