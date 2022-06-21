## SECTIONS REQUIRING INPUT:
* Paths to input and output directories and models
* ROIs [earlier version of this code required all ROIs to be in the same directory; this does not, though it does require manual listing of ROI paths and names]
* Subjects
* Defining conditions in the model - - make sure that you include ALL task conditions / regressors in this list. Code is provided to read out regressor names from first level SPM.mat files.
* Edit contrasts in LOOP section (do not need to edit other settings in LOOP)


This notebook uses a Matlab kernel and is a wrapper of the gPPI toolbox in SPM12. 