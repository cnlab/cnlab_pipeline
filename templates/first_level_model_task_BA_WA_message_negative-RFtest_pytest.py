#!/usr/bin/env python
# coding: utf-8

# # Building first level models using _nipype_ and _SPM12_
# 
# ## model_task : WA,BA goal framing
# 
# -------
# #### History
# 
# * 8/31/2020 nc - edited .py to allow command line input of subject list 
# * 8/27/2020 nc - testing and export to py script
# * 8/21/2020 mbod - update master branch to do FAST and residuals
# * 8/21/2020 nc - update to new refactored from MURI including FAST 
# * 6/22/2020 nc - using refactored script on BA/WA models
# * 6/17/2020 mbod - testing refactored script to duplicate notebook setup 
# * 11/7/2019 jeesung - 1st level models without pmod
# * 3/4/19 cscholz - modify notebook for darpa 1 first-level model
# * 2/27/19 mbod  - modify example notebook to make template
# 
# -----
# 
# ### Description
# 
# * Set up a nipype workflow to use SPM12 to make first level models for _megameta_ task data (preprocessed using `batch8` SPM8 scripts) in BIDS derivative format   
# 

# -------------------
# ### Step 1: SET NEEDED PARAMETERS
# 
# * Checked out branch for `cnlab/GLM` is:
#     * /data00/tools/cnlab_pipeline/cnlab/GLM/
#     
#     
# 
# 
# * Import modules and the pipeline code
#  
# 
# 
# 

# In[2]:


import os
import sys
import argparse


sys.path.append('/data00/projects/megameta/scripts/jupyter_megameta/cnlab_pipeline/')
from cnlab.GLM import first_level


# ### Set params
# 
# * There are some parameters that need to be set
#     1. `MODEL_SPEC_FILE` is the name of the JSON model file
#     2. `MODEL_PATH` - at the moment the JSON files have been kept in a folder called `model_specifications` at the same level as the first_level notebook folder - but we could decide of another convention for where to place the JSON file and maybe it should live in the same folder as the notebook/script that runs the model pipeline.

# In[4]:


MODEL_SPEC_FILE = 'BA_negative_model_WALKSTATEMENT_message_RFfast.json' # replace with filename of JSON file

MODEL_PATH = os.path.abspath(
                    os.path.join('../model_specifications',
                                  MODEL_SPEC_FILE)
)



# -------------------------
# 
# ### Step 2: Run the `setup_pipeline` function
# 
# * This requires:
#     1. `MODEL_PATH` - the full/absolute path to the model JSON file (__REQUIRED__)
#     2. `include_subjects` - pass the list if you want to specify specific subjects to include (optional) - see below
#     3. `exclude_subjects` - pass the list if you want to specify specific subjects to exclude (optional) - see below
#     4. `DEBUG` - default is `False` whether to print out debugging info when setting up the pipeline (optional)

# ### Specify subjects to process
# 
# * By default the module will look for all the subjects that have BIDS data in the project folder.
# 
# 
# * But you can select specific subjects to:
#     1. be included with a list `include_subjects`
#     2. be excluded with a list `exclude_subjects`
# 
# * The two lists will be combined to give the resulting subject list passed to the pipeline

parser = argparse.ArgumentParser()
parser.add_argument('--include_subjects', default = "NA", nargs='*', help='optional: specify included subjects')
parser.add_argument('--exclude_subjects', default = "NA", nargs='*', help='optional: specify excluded subjects')
args = parser.parse_args()


## ncooper - you could refactor this to be simpler because the setup_pipeline function already
## does the logic of creating a union of include and exclude and if both or either is None
## it can cope with that
## so you could have them None by default and just one called to setup_pipeline

if args.include_subjects != "NA" and args.exclude_subjects != "NA":
    include_subjects = args.include_subjects
    exclude_subjects = args.exclude_subjects
    model_def=first_level.setup_pipeline(MODEL_PATH,
                          include_subjects=include_subjects,
                          exclude_subjects = exclude_subjects,
                          DEBUG=False)

elif args.include_subjects != "NA" and args.exclude_subjects == "NA":
    include_subjects = args.include_subjects
    model_def=first_level.setup_pipeline(MODEL_PATH,
                          include_subjects=include_subjects,
                          DEBUG=False)

elif args.exclude_subjects != "NA" and args.include_subjects == "NA":
    exclude_subjects = args.exclude_subjects
    model_def=first_level.setup_pipeline(MODEL_PATH,
                          exclude_subjects=exclude_subjects,
                          DEBUG=False)
    
else:
    model_def=first_level.setup_pipeline(MODEL_PATH,
                          DEBUG=False)

print("subject list:")
print(model_def['subject_list'])


# In[9]:  #### TEST IF THESE ARE NECESSARY


model_def['unzip_and_smooth']=False
model_def['resolutions'] = ['medium']
model_def['smoothing_list'] = [8]


# -------------------------
# 
# ### Step 3: CHECK WORKFLOW AND DIRECTORIES

# In[10]:


pipeline=first_level.build_pipeline(model_def)


# -------------------------
# 
# ### Step 4: RUN PIPELINE

# In[14]:


run_graph=pipeline.run(plugin='Linear', plugin_args={})

