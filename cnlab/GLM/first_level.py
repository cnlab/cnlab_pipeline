#!/usr/bin/env python
# coding: utf-8

# # Building first level models using _nipype_ and _SPM12_
# 
# ## Base functionality for _megameta_ project
# 
# -------
# #### History
# 
# * 3/28/19 mbod - update pipeline to include resampling to template & SPM path reference
# * 3/23/19 mbod - include contrast definition in the config JSON file
# * 3/9/19 mbod - updates from testing template with `darpa1`
# * 2/27/19 mbod  - modify example notebook to make base functionality notebook
# 
# -----
# 
# ### Description
# 
# * Set up a nipype workflow to use SPM12 to make first level models for _megameta_ task data (preprocessed using `batch8` SPM8 scripts) in BIDS derivative format   
# 

# -------------------
# 
# ### Template variables
# 
# * Specify the following values:
#     1. project name - should be name of folder under `/data00/project/megameta`, e.g. `project1`
#     2. filename for JSON model specification (should be inside `model_specification` folder), e.g. `p1_image_pmod_likeme.json`
#     3. TR value in seconds
#  
# 
# 

# -------------------
# 
# ### Setup
# 
# * import required modules and define parameters

# In[1]:


import os  # system functions

# NIYPE FUNCTIONS
import nipype.interfaces.io as nio           # Data i/o
import nipype.interfaces.spm as spm          # spm
import nipype.interfaces.matlab as mlab      # how to run matlab
import nipype.interfaces.utility as util     # utility
import nipype.pipeline.engine as pe          # pypeline engine
import nipype.algorithms.modelgen as model   # model specification
from nipype.interfaces.base import Bunch
from nipype.algorithms.misc import Gunzip

from itertools import combinations

from nilearn import plotting, image
from nistats import thresholding


from IPython.display import Image


import scipy.io as sio
import numpy as np
import json
import pandas as pd


# #### Matlab path
# 

# In[2]:


# Set the way matlab should be called
mlab.MatlabCommand.set_default_matlab_cmd("matlab -nodesktop -nosplash")
# If SPM is not in your MATLAB path you should add it here
mlab.MatlabCommand.set_default_paths(PATH_TO_SPM_FOLDER)


# ### Parameters
# 
# * These need to be reformatted to be consistent
# * as data is not smoothed commented out the `fwhm_size` param - but data probably has a value

# #### Load JSON model config

# In[5]:


JSON_MODEL_FILE = os.path.join('/data00/projects/megameta/scripts/jupyter_megameta/first_level_models',
                               PROJECT_NAME, 'model_specifications',
                               MODEL_SPEC_FILE)


# In[6]:


with open(JSON_MODEL_FILE) as fh:
    model_def = json.load(fh)


# In[7]:


TASK_NAME = model_def['TaskName']
RUNS = model_def['Runs']
MODEL_NAME = model_def['ModelName']
PROJECT_NAME = model_def['ProjectID']


# In[8]:


PROJECT_DIR = os.path.join('/data00/projects/megameta', PROJECT_NAME)
SUBJ_DIR = os.path.join(PROJECT_DIR, 'derivatives', 'batch8')


task_func_template = "{PID}_task-{TASK}_run-0{RUN}_space-MNI152-T1-1mm_desc-preproc_bold.nii.gz"


subject_list = [subj for subj in os.listdir(SUBJ_DIR) 
                   if os.path.exists(os.path.join(SUBJ_DIR,subj,'func',
                                        task_func_template.format(PID=subj, TASK=TASK_NAME, RUN=1)))]

output_dir = os.path.join(PROJECT_DIR,'derivatives', 'nipype','model_{}_{}'.format(TASK_NAME.upper(), MODEL_NAME))        # name of 1st-level output folder
working_dir = os.path.join(PROJECT_DIR, 'working', 
                           'nipype', 'workingdir_model_{}_{}'.format(TASK_NAME.upper(), MODEL_NAME))   # name of 1st-level working directory



# In[9]:


# check to see if output and work directories exist

if not os.path.exists(output_dir):
    os.makedirs(output_dir) 

if not os.path.exists(working_dir):
    os.makedirs(working_dir)


# In[4]:



try:
    subject_list = [ s for s in subject_list if s not in exclude_subjects ]
    print('\n\nApplied subject inclusion list:\n\t',' '.join(exclude_subjects))
except:
    print('\n\nNo subject exclusions applied')

try:
    subject_list = [ s for s in subject_list if s in include_subjects ]
    print('\n\nApplied subject inclusion list:\n\t',' '.join(include_subjects))
except:
    print('\n\nNo subject inclusions applied')

    
print('\n\nSUBJECT LIST IS:\n\t', ' '.join(subject_list))


# ### Utility functions for subject info and contrasts

# ### Setup design matrix data for subject
# 
# * need a function to set up the nipype `Bunch` format used
#     * https://nipype.readthedocs.io/en/latest/users/model_specification.html
# * read the onsets/dur/conditions from task logs and extract needed data

# In[9]:


def get_subject_info(subject_id, model_path, DEBUG=False):
    '''
    1. load model specification from JSON spec file
    2. get confound file for subject for task to add to design matrix     
    3. get task spec CSV for subject for task
    4. setup subject info structure
    '''
    
    import os
    import pandas as pd
    import json
    
    from nipype.interfaces.base import Bunch
    
    
    def make_pmod(df, conditions, pmods={}, normalize='mean'):
        
        pmod = []
        
        for cond in conditions:
        
            
            if not pmods.get(cond):
                pmod.append(None)
            else:
                df2 = df[df.trial_type==cond]
                
                pmod_name = pmods.get(cond)
                
                #pmod = [pmod] if not type(pmods) is list else pmod
                
                # MAKE SURE THERE IS VARIANCE IN PMOD VECTOR
                if df2[pmod_name].var()==0:
                    df2[pmod_name]+=0.001
                    
                # APPLY NORMALIZATION
                if normalize=='mean':
                    df2[pmod_name] = df2[pmod_name] - df2[pmod_name].mean()

                pmod.append(Bunch(name=[pmod_name],
                      param=[df2[pmod_name].values.tolist()
                            ],
                      poly=[1]
                     ))
        
        return pmod
    
    

    def map_spec_to_model(spec_df,model):
        """
        Maps spec trial names to model contrast trials.

        Args:
            spec: the events.tsv spec file
            model: the model.json file

        Returns:
            pandas dataframe object
        """

        spec=spec_df.copy()

        for con in model['Conditions']:
            if DEBUG:
                print('con is',con, model.get('ModelAsItem'))

            spec_trials = model['Conditions'][con]

            # --- mbod 9/18/19 added modelasitem 
            if model.get('ModelAsItem') and con in model['ModelAsItem']:
                pass
            else:
                spec.loc[spec.trial_type.isin(spec_trials),'trial_type'] = con
            spec.onset.sort_values()
            
        if DEBUG:
            print(spec) 
        
        return spec
    
    
    with open(model_path) as fh:
        model_def = json.load(fh)
    
    
    pmod = None if not model_def.get('Modulators') else []
    
    TASK_NAME = model_def['TaskName']
    TASK_RUNS = model_def['Runs']
    MODEL_NAME = model_def['ModelName']
    PROJECT_ID = model_def['ProjectID']
    
    condition_names = list(model_def['Conditions'].keys())
   
    
    PROJECT_DIR = os.path.join('/data00/projects/megameta', PROJECT_ID)
    SUBJ_DIR = os.path.join(PROJECT_DIR,'derivatives', 'batch8')
    
    realign_files = []
    subject_info = []
    

    
    # check to see which runs exist for subject
    # by looking for appropriate events.tsv files
    # this could (should?) also include looking for the nifti file?
    runs_for_subj = [run for run in TASK_RUNS
                    if 
                    os.path.exists(os.path.join(SUBJ_DIR, subject_id, 'func',
                                    '{}_task-{}_run-0{}_events.tsv'.format(subject_id, 
                                                                           TASK_NAME,
                                                                           run)))
                    ]
    
    if DEBUG:
        print("runs_for_subj", runs_for_subj)
        print("checked paths:")
        for run in TASK_RUNS:
            print('\t', os.path.join(SUBJ_DIR, subject_id, 'func',
                                    '{}_task-{}_run-0{}_events.tsv'.format(subject_id, 
                                                                           TASK_NAME,
                                                                           run)))
        print("TASK NAME", TASK_NAME)
        print("pmod", pmod)
        print("TASK_RUNS", TASK_RUNS)
        print("subject_id", subject_id)
    
    
    for run_num, _ in enumerate(runs_for_subj,1):
    
    
        events_df = pd.read_csv(os.path.join(SUBJ_DIR, subject_id, 'func',
                                             '{}_task-{}_run-0{}_events.tsv'.format(subject_id, 
                                                                                    TASK_NAME,
                                                                                    run_num)),
                               sep='\t')

        
        # mbod 9/18/2019 - adding ModelAsItem functionality
        # check to see if ModelAsItems field exists in the JSON and split conditions
        if model_def.get('ModelAsItem'):
            for cond_as_item in model_def['ModelAsItem']:
                if DEBUG:
                    print('***',cond_as_item)
                if condition_names.count(cond_as_item)>0:
                    condition_names.remove(cond_as_item)
            
                cond_trial_types = model_def['Conditions'][cond_as_item]
                condition_names.extend(events_df[events_df.trial_type.isin(cond_trial_types)]['trial_type'].unique())
                print(condition_names)


    
    
        onsets_df = map_spec_to_model(events_df, model_def)
        
                
        realign_file = os.path.join(PROJECT_DIR, 'working','nipype',
                                        'workingdir_model_{}_{}'.format(TASK_NAME.upper(),MODEL_NAME),
                                        '{}-run-0{}-realign.txt'.format(subject_id, run_num))

        confound_file=os.path.join(SUBJ_DIR, subject_id, 'func',
                                   '{}_task-{}_run-0{}_desc-confounds-regressors.tsv'.format(subject_id, 
                                                                                             TASK_NAME,
                                                                                             run_num)
                                   )

        confound_df = pd.read_csv(confound_file, sep='\t')

        cols_to_use = [ 'TransX','TransY', 'TransZ', 'RotX', 'RotY', 'RotZ']

        confound_df[cols_to_use].to_csv(realign_file, 
                                        header=False, 
                                        index=False,
                                        sep='\t')

        realign_files.append(realign_file)

        onsets = []
        dur = []

        for cond in condition_names:
            if DEBUG:
                print('val of cond',cond)
            onsets.append(onsets_df[onsets_df.trial_type==cond].onset.values)
            dur.append(onsets_df[onsets_df.trial_type==cond].duration.values)

                    
            
         
        #pmod = make_pmod(rdf, condition_names)   
        
        if model_def.get('Modulators'):
            pmod = make_pmod(onsets_df, condition_names, 
                         pmods=model_def['Modulators'])     

 
        
        
        subject_info.append(Bunch(conditions=condition_names,
                         onsets=onsets,
                         durations=dur,
                         amplitudes=None,
                         tmod=None,
                         pmod=pmod,
                         regressor_names=None,
                         regressors=None))
    
    
    DM_regressors = []
    for cond in condition_names:
        DM_regressors.append(cond)
        if pmod and model_def['Modulators'].get(cond):
            DM_regressors.append('{}x{}^1'.format(cond, model_def['Modulators'].get(cond)))
    
    
    return subject_info, realign_files, DM_regressors


# ### Set up contrasts
# 
# * This part of the template needs work to provide a cleaner way to specify contrasts
# * Could use the same vector contrasts approach as we have in batch8 and then have a function to convert this into the list of list data structure that nipype spm contrasts node looks for

# In[ ]:


def make_contrast_list(subject_id, condition_names, model_path, DEBUG=False):
    
    import json
    
    condition_names.append('constant')

    cont = []
    for idx, cname in enumerate(condition_names):
        ccode = [0 if pos!=idx else 1 for pos in range(len(condition_names))]
        cont.append([cname, 'T', condition_names, ccode])

    # add custom contrasts from the JSON model file
    with open(model_path) as fh:
        model_def = json.load(fh)
    
    contrasts = model_def.get('Contrasts')
    
    if not contrasts:
        return cont
    
    for contrast in contrasts:
            
        cname = contrast['name']
        
        
        
        pos_idx = [condition_names.index(p) for p in contrast['pos']]
        neg_idx = [condition_names.index(n) for n in contrast['neg']]
        
        pos_length = len(contrast['pos'])
        neg_length = len(contrast['neg'])
        
        ccode = []
        for idx, _ in enumerate(condition_names):
            if idx in pos_idx:
                ccode.append(1/pos_length)
            elif idx in neg_idx:
                ccode.append(-1/pos_length)
            else:
                ccode.append(0)

        cont.append([cname, 'T', condition_names, ccode])
        
                
        if DEBUG:
            print(contrast)
            print(ccode)
            
    return cont


# ## Set up processing nodes for modeling workflow

# #### Specify model node

# In[ ]:


# SpecifyModel - Generates SPM-specific Model
modelspec = pe.Node(model.SpecifySPMModel(concatenate_runs=False,
                                 input_units='secs',
                                 output_units='secs',
                                 time_repetition=TR,
                                 high_pass_filter_cutoff=128),
                                 output_units = 'scans',
                 name="modelspec")


# #### Level 1 Design node
# 
# ** TODO -- get the right matching template file for fmriprep **
# 
# * ??do we need a different mask than:
# 
#     `'/data00/tools/spm8/apriori/brainmask_th25.nii'`

# In[ ]:


# Level1Design - Generates an SPM design matrix
level1design = pe.Node(spm.Level1Design(bases={'hrf': {'derivs': [0, 0]}},
                                 timing_units='secs',
                                 interscan_interval=TR,
                                 model_serial_correlations='none', #'AR(1)',
                                 mask_image = '/data00/tools/spm8/apriori/brainmask_th25.nii',
                                 global_intensity_normalization='none'
                                       ),
                    name="level1design")


# #### Estimate Model node

# In[ ]:


# EstimateModel - estimate the parameters of the model
level1estimate = pe.Node(spm.EstimateModel(estimation_method={'Classical': 1}),
                      name="level1estimate")


# #### Estimate Contrasts node

# In[ ]:


# EstimateContrast - estimates contrasts
conestimate = pe.Node(spm.EstimateContrast(), name="conestimate")


# In[ ]:





# ## Setup pipeline workflow for level 1 model

# In[ ]:


# Initiation of the 1st-level analysis workflow
l1analysis = pe.Workflow(name='l1analysis')

# Connect up the 1st-level analysis components
l1analysis.connect([(modelspec, level1design, [('session_info',
                                                'session_info')]),
                    (level1design, level1estimate, [('spm_mat_file',
                                                     'spm_mat_file')]),
                    (level1estimate, conestimate, [('spm_mat_file',
                                                    'spm_mat_file'),
                                                   ('beta_images',
                                                    'beta_images'),
                                                   ('residual_image',
                                                    'residual_image')])
                   
                    ])


# In[ ]:





# ## Set up nodes for file handling and subject selection

# ### `getsubjectinfo` node 
# 
# * Use `get_subject_info()` function to generate spec data structure for first level model design matrix

# In[ ]:





# In[ ]:


# Get Subject Info - get subject specific condition information
getsubjectinfo = pe.Node(util.Function(input_names=['subject_id', 'model_path'],
                               output_names=['subject_info', 'realign_params', 'condition_names'],
                               function=get_subject_info),
                      name='getsubjectinfo')


# In[ ]:


makecontrasts = pe.Node(util.Function(input_names=['subject_id', 'condition_names', 'model_path'],
                                     output_names=['contrasts'],
                                      function=make_contrast_list),
                    name='makecontrasts')


# ### `infosource` node
# 
# * iterate over list of subject ids and generate subject ids and produce list of contrasts for subsequent nodes

# In[ ]:


# Infosource - a function free node to iterate over the list of subject names
infosource = pe.Node(util.IdentityInterface(fields=['subject_id', 'model_path']
                                   ),
                  name="infosource")

infosource.iterables = [('subject_id', subject_list),
                        ('model_path', [JSON_MODEL_FILE]*len(subject_list))
                        
                       ]



# ### `selectfiles` node
# 
# * match template to find source files (functional) for use in subsequent parts of pipeline

# In[ ]:


# SelectFiles - to grab the data (alternativ to DataGrabber)


## TODO: here need to figure out how to incorporate the run number and task name in call
templates = {'func': '{subject_id}/func/{subject_id}_task-'+TASK_NAME+'_run-0*_space-MNI152-T1-1mm_desc-preproc_bold.nii.gz'}          


selectfiles = pe.Node(nio.SelectFiles(templates,
                               base_directory='/data00/projects/megameta/{}/derivatives/batch8'.format(PROJECT_NAME)),
                      working_dir=working_dir,
                   name="selectfiles")


# ## Unzip and smoothing steps
# 
# * BIDS derivatives folders contain unsmoothed functional NIFTI files in zipped (.nii.gz) format
# * This subflow adds three nodes:
#     1. gunzip
#     2. resample
#     3. smooth

# #### Specify unzip node
# 
# * transform `.nii.gz` to `.nii`

# In[ ]:


gunzip = pe.MapNode(Gunzip(),name="gunzip", iterfield=['in_file'])


# #### Specify smoothing node

# In[ ]:


smooth = pe.Node(interface=spm.Smooth(), name="smooth")
#fwhmlist = [4,6,8]
fwhmlist = [8]
smooth.iterables = ('fwhm', fwhmlist)


# #### Specify resampling node

# In[ ]:


resample = pe.MapNode(interface=spm.utils.Reslice(), 
                      name='resample',
                     iterfield=['in_file'])


# In[1]:


resample.inputs.space_defining = '/data00/projects/megameta/templates/reference_medium_wad.nii'


# In[ ]:


unzip_resample_and_smooth = pe.Workflow(name='unzip_resample_and_smooth')

unzip_resample_and_smooth.base_dir = os.path.join(SUBJ_DIR, working_dir)

unzip_resample_and_smooth.connect(
    [
        (gunzip, resample, [('out_file', 'in_file')]),
        (resample, smooth, [('out_file', 'in_files')])
    ]
)


# ### Specify datasink node
# 
# * copy files to keep from various working folders to output folder for model for subject

# In[ ]:


# Datasink - creates output folder for important outputs
datasink = pe.Node(nio.DataSink(base_directory=SUBJ_DIR,
                         parameterization=True, 
                         #container=output_dir      
                               ),
                name="datasink")

datasink.inputs.base_directory = output_dir

# Use the following DataSink output substitutions
substitutions = []
subjFolders = [('_model_path.*subject_id_%s/_fwhm_%s' % (sub,f), 'fwhm_%s' % (f))
               for f in fwhmlist
               for sub in subject_list]
substitutions.extend(subjFolders)
datasink.inputs.regexp_substitutions = substitutions



# ---------

# ## Set up workflow for whole process

# In[ ]:


pipeline = pe.Workflow(name='first_level_model_{}_{}'.format(TASK_NAME.upper(),MODEL_NAME))
pipeline.base_dir = os.path.join(SUBJ_DIR, working_dir)


pipeline.connect([(infosource, selectfiles, [('subject_id', 'subject_id')]),
                  (infosource, getsubjectinfo, [('subject_id', 'subject_id'),
                                                ('model_path', 'model_path')
                                               ]),
                  (infosource, makecontrasts, [('subject_id', 'subject_id'),
                                               ('model_path', 'model_path')
                                              ]),
                  (getsubjectinfo, makecontrasts, [('condition_names', 'condition_names')]),
                  
                 (getsubjectinfo, l1analysis, [('subject_info',
                                                 'modelspec.subject_info'),
                                                ('realign_params',
                                                   'modelspec.realignment_parameters')]),
                  (makecontrasts, l1analysis, [('contrasts',
                                             'conestimate.contrasts')]),
                  
                  
                  (selectfiles, unzip_resample_and_smooth, [('func','gunzip.in_file')]),
                  
    
                  (unzip_resample_and_smooth, l1analysis, [('smooth.smoothed_files',
                                          'modelspec.functional_runs')]),
                  
                                    
                  (infosource, datasink, [('subject_id','container')]),
                  (l1analysis, datasink, [('conestimate.spm_mat_file','@spm'),
                                          ('level1estimate.beta_images','@betas'),
                                          ('level1estimate.mask_image','@mask'),
                                          ('conestimate.spmT_images','@spmT'),
                                          ('conestimate.con_images','@con'),
                                          ('conestimate.spmF_images','@spmF')
                                         ])
                 ]
)
                  


# In[ ]:




