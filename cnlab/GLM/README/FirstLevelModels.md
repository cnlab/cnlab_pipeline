# STEP 1 - Set up model specification json(s) 
## Model specification structure
* Three separate json files (analysis, task, job) contribute to model specifications. The analysis and task files are relatively stable; only the job json file should be created for each analysis instance. 
  * The **analysis template** containing components of model specification that have (lab-wide) defaults for a specific analysis type (e.g., regular GLM with 8mm smoothing, or single-trial beta with no smoothing for MVPA). This includes environment paths, model components like smoothing kernel, high pass filter cutoff and motion regressors, and design features like the HRF basis and modeling of serial correlations.
    * See data00/projects/megameta/scripts/jupyter_megameta/l1analysis/template-megameta.json
  * The **task template** should contain particular information for each task in a study. This includes information like the number of task runs, the TR, which subjects/runs to exclude from modeling, and the locations of the functional, event, and motion regressor files. 
    * See data00/projects/megameta/scripts/jupyter_megameta/l1analysis/darpa1/task-share.json
  * A **job specification** file. This should be edited for each analysis. This includes specification of mapping or melting components of the `events.tsv` into regressors, and specification of contrasts.
    * See data00/projects/megameta/scripts/jupyter_megameta/l1analysis/darpa1/task-share_model-beta.json
 * If a model differs from the **analysis** and **task** defaults (e.g., a different smoothing kernel is desired, or a different set of motion regressors), those features can be specified in the job json file and will override the defaults 
 * All specification options are described in data00/projects/megameta/scripts/jupyter_megameta/cnlab_pipeline/cnlab/GLM/l1analysis_description.json. Some of these are specific to our lab, but many refer back to nipype options
 * File organization (*under discussion*)
   * Analysis level json files live in data00/projects/megameta/scripts/jupyter_megameta/l1analysis
   * Task and job level json files live in a study specific sub-folder, e.g. data00/projects/megameta/scripts/jupyter_megameta/l1analysis/darpa1

## Overview of first-level analysis pipeline
The analysis pipeline is a lightweight wrapper of the `nipype` package. The analysis flow is (click on links to see details of the modules):
1. Smoothing (optional) - [IsotropicSmooth](https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.fsl.maths.html#isotropicsmooth)
2. Specify SPM model - [SpecifySPMModel](https://nipype.readthedocs.io/en/latest/api/generated/nipype.algorithms.modelgen.html#specifyspmmodel)
3. Generate design matrix - [Level1Design](https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.spm.model.html#level1design)
4. Estimate model - [EstimateModel](https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.spm.model.html#estimatemodel)
5. Estimate contrast - [EstimateContrast](https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.spm.model.html#estimateconstrast)

## Notes on specification options of particular interest
* Regressor specification
  * `map_event` groups multiple event categories into one regressor
  * `melt_event` generates individual events such as single-trial betas
  * `include_event` and `exclude_event` can be used when the user does not want to model all event types specified in the events file (for example, if fixation is specified as an event category, and the user wants to exclude this from explicit modeling)
* Contrast specification
  * The contrasts field expects as input: name of the contrast, statistic (T), the names of the regressors to be contrasted, and the contrast coding (click [here](https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.spm.model.html#estimatecontrast) for details).     

## Note on events.tsv files
* The BIDS default is for there to be a single set of events.tsv files for any task. Given the variety of adaptations on events.tsv we've seen in our lab, the current setup allows for flexibility in the events.tsv files and also stores a copy of the final events file for each model (*naming in discussion*). Events.tsv is specified in the model-task level json, but can be overridden by specifying a different file in the model-specific json. Further, a definitive events file is saved out within each model folder for later reference (e.g., /data00/projects/megameta/BA/models/task-walkstatement_model-message/events).
 
# STEP 2 - Create slurm jobs for running first level models
* Copy the codes in `Example script.ipynb` into a new notebook.
* This script takes in the model specification and template (analysis and task) jsons, processes the information from each, and creates and writes a slurm .job file for each participant. 
* The slurm job processes the first level model using /data00/projects/megameta/scripts/jupyter_megameta/cnlab_pipeline/cnlab/GLM/l1analysis_SPM.py 
* Options here are to create a short wrapper notebook around this first level model template notebook, or to make the template notebook an executable script (*under discussion*) 

# STEP 3 - Execute and monitor slurm jobs 
* To submit the jobs (and run the models), log on to the slurm cluster, and follow the notebook instructions - change directory to, e.g., /data00/projects/megameta/product/models/task-ARF_model-bin3/slurm. Each job can be pasted separately at the command line. 
* Logging on to slurm
  * Connect to the ASC VPN
  * Open a terminal window on your computer
  * Type: ssh <JANUS_UN>@asc.upenn.edu@cls000 - - replace your JANUS username, and enter your JANUS password when prompted
    * If you get a prompt “The authenticity of host 'cls000 (10.30.12.140)' can't be established. … Are you sure you want to continue connecting (yes/no)” → type yes
* Monitoring jobs in the queue
  * Type "squeue" to see what jobs are in the queue (currently being processed or waiting to be processed)
  * Errors will be written into files with a .err extension (e.g., sub-BA277.err) within an "out" directory nested in the relevant slurm directory (e.g., /data00/projects/megameta/BA/models/task-walkstatement_model-message/slurm/out)
  * View the .err files in any file viewer

# STEP 4 - Check model output 
[link to template notebook]









