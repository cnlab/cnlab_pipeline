# STEP 1 - Set up model specification json(s) 
* Three separate json files contribute to model specification; only one should be edited for most models. 
  * The base level is a lab-wide template (*the number of templates is under discussion*) containing components of model specification that have lab-wide defaults. This includes environment paths, model components like high pass filter cutoff and motion regressors, and design features like the HRF basis and modeling of serial correlations.
    * See data00/projects/megameta/scripts/jupyter_megameta/l1analysis/template-megameta.json
  * A study-task-level template builds on the base template. This should be created for each task in a study, at the start of processing for each new study. This includes information like the number of task runs, the TR, which subjects/runs to exclude from modeling, and the location of the functional, event.tsv, and motion regressor files. 
    * See data00/projects/megameta/scripts/jupyter_megameta/l1analysis/darpa1/task-share.json
  * A model-level specification file. This should be edited for each model. This includes specification of mapping or melting components of the events.tsv into regressors, and specification of contrasts.
    * See data00/projects/megameta/scripts/jupyter_megameta/l1analysis/darpa1/task-share_model-bin3.json
 * If a model differs from the study-task or base defaults (e.g., a different smoothing kernel is desired, or a different set of motion regressors), those features can be specified in the model-level json file and will override the defaults 
 * All specification options are described in data00/projects/megameta/scripts/jupyter_megameta/cnlab_pipeline/cnlab/GLM/description.json. Some of these are specific to our lab, but many refer back to nipype options
* File organization (*under discussion*)
  * The base level json files live in data00/projects/megameta/scripts/jupyter_megameta/l1analysis
  * Study-task and model level json files live in a study specific sub-folder, e.g. data00/projects/megameta/scripts/jupyter_megameta/l1analysis/darpa1
 
# STEP 2 - Create slurm jobs for running first level models
* First level model notebook: data00/projects/megameta/scripts/jupyter_megameta/cnlab_pipeline/cnlab/GLM/First Level Template - Convert and create slurm scripts.ipynb
* This script takes in the model specification json (model.json) as input, processes information from model.json and the template (base and study-task level) jsons, and writes a slurm .job file for each participant. 
* The slurm job processes the first level model using /data00/projects/megameta/scripts/jupyter_megameta/cnlab_pipeline/cnlab/GLM/first_level_HY.py 


# STEP 3 - Execute and monitor slurm jobs 
* To submit the jobs (and run the models), log on to the slurm cluster (INSTRUCTIONS), and follow the notebook instructions - change directory to, e.g., /data00/projects/megameta/product/models/task-ARF_model-bin3/slurm. Each job can be pasted separately at the command line. 
* INSTRUCTIONS for checking the queue and error output files 


# STEP 4 - Check model output 
[link to template notebook]









