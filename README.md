# HRTF PCA

Performs a principal component analysis (PCA) on a set of head-related transfer functions (HRTFs). It reconstructs 128-pt HRTF's from the UC Davis CIPIC dataset from 1 (128x data reduction) base to 128 bases (no data reduction).

## Instructions
1. Download this repository.
2. Download the CIPIC dataset from the [UC Davis website](https://www.ece.ucdavis.edu/cipic/spatial-sound/hrtf-data/). Navigate to the "Full Database Download" section and download the .zip file (~170MB).
3. Unzip the file into the repository.
4. Make sure the HRTF directory paths point to the "standard_hrir_database/" folder in the CIPIC dataset. This needs to be checked in project_run.m (line 9, input.directory) and hrtf_load.m (line 12, input.dir).
5. Run run_project.m and resolve any runtime fixes that I may have forgotten.

**Note:** This project takes a while to run. It took about 10 minutes in total, of which 7 minutes were spent reconstructing the HRTFs and calculating loss functions, on my personal laptop (Intel i5-5300 4CPU @ 2.30GHz; 12GB RAM). It also spits out five plot figures.

**Note:** There is functionality to load the anthropomorphic data from CIPIC, but it is not used in the project.

## Pictures

My [write-up](https://github.com/alextongue/hrtf-pca/blob/master/writeup/Tung_HRTFPCA_report.pdf) could be found in this repo. It explains the photos below.

![Some ITD and IID slices](https://github.com/alextongue/hrtf-pca/blob/master/writeup/itd_iid.png?raw=true "Some ITD and IID slices")

![Preprocessing done to HRTF's (colloquially Directional Transfer Functions, or DTFs](https://github.com/alextongue/hrtf-pca/blob/master/writeup/hrtf_dtf.png?raw=true "Preprocessing done to HRTF's (colloquially Directional Transfer Functions, or DTFs")

![Covariance and correlation matrices of the mean-subtracted HRTFs](https://github.com/alextongue/hrtf-pca/blob/master/writeup/covar_corre.png?raw=true "Covariance and correlation matrices of the mean-subtracted HRTFs")

![The first five orthogonal bases that can be used to decompose HRTFs](https://github.com/alextongue/hrtf-pca/blob/master/writeup/components.png?raw=true "The first five orthogonal bases that can be used to decompose HRTFs")

![A reconstruction of an HRTF and loss functions](https://github.com/alextongue/hrtf-pca/blob/master/writeup/hrtf2_loss.png?raw=true "A reconstruction of an HRTF and loss functions")