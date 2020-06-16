#!/usr/bin/env python
# coding: utf-8
# %%

# # Building second level models using _nipype_ and _SPM12_
# 
# ## Base functionality for _megameta_ project
# 
# -------
# #### History
# * 5/4/19 cscholz - add datasink, incorporate mreg design, incorporate sampling of first-level contrast based on percentage of available first-level models per project
# * 4/15/19 mbod - incorporate function to read the 2nd level JSON model config
# * 4/9/19 mbod - modify template to work with fmriprep processed data
# * 3/20/19 mbod - initial setup for testing some simple one sample t-test models
# -----
# 
# ### Description
# 
# * Set up a nipype workflow to use SPM12 to make second level models for _megameta_ task data (preprocessed using `batch8` SPM8 scripts) in BIDS derivative format   
# 

# ### Setup

# %%


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

import scipy.io as sio
import numpy as np
import json
import pandas as pd

import random

from IPython.display import Image


from itertools import product


# #### Matlab path
# 

# %%


# Set the way matlab should be called
mlab.MatlabCommand.set_default_matlab_cmd("matlab -nodesktop -nosplash")
# If SPM is not in your MATLAB path you should add it here
mlab.MatlabCommand.set_default_paths(PATH_TO_SPM_FOLDER)


# %%


GROUP_DIR = '/data00/projects/megameta/group_models/'


# #### Load JSON model config

# %%


JSON_MODEL_FILE = os.path.join('/data00/projects/megameta/scripts/jupyter_megameta/second_level_models',
                               'model_specifications',
                               MODEL_SPEC_FILE)


# %%


with open(JSON_MODEL_FILE) as fh:
    model_def = json.load(fh)


# %%


MODEL_NAME = model_def['ModelName']

CONTRASTS = model_def['Contrasts']

ROOT_DIR = '/data00/projects/megameta'


# %%


l2_contrast_list = CONTRASTS # list of specific contrast files to use in 2nd level model (include .nii?)


output_dir = os.path.join(GROUP_DIR,'derivatives', 'nipype','model_2nd-level_{}'.format(MODEL_NAME))        
working_dir = os.path.join(GROUP_DIR, 'working', 
                           'nipype', 'workingdir_model_2nd-level_{}'.format(MODEL_NAME))   


# %%


if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    
if not os.path.exists(working_dir):
    os.makedirs(working_dir)


# ## Get list of contrast files

# %%


def process_project(project_name, model_def=model_def, exclude_subjects=exclude_subjects ,scan_all_subjs=False, DEBUG=False):
    
    project_spec = [pspec for pspec in model_def['Projects'] if pspec['Name']==project_name]
    
    if not project_spec:
        print('Cannot find specification for project: ', project_name)
        return None
    
    model_name = project_spec[0]['Model']
    cmap = project_spec[0]['ContrastMap']
    
    
    model_dir = os.path.join(ROOT_DIR, project_name, 
                             "derivatives", "nipype",
                             "model_{}".format(model_name)
                            )
    
    if not os.path.exists(model_dir):
        print('Cannot find first level model directory:', model_dir)
        return None
    
    subjs_with_models = [s for s in os.listdir(model_dir) if s.startswith('sub-')]
    #exclude_people
    subjs_with_models=[s for s in subjs_with_models if s not in exclude_subjects]
    
    if DEBUG:
        print("Found {} first level subject models\n".format(len(subjs_with_models)))
    
    
    contrast_lists = { cname: [] for cname in cmap}
    
    
    model_contrasts=None
    for sidx,subj in enumerate(subjs_with_models):
        
        if DEBUG:
            print('Processing',subj, '-',end='')
        
        first_level_dir = os.path.join(model_dir, subj, 'medium', 'fwhm_8')

        if scan_all_subjs or sidx==0:
            spm_mat_file = os.path.join(first_level_dir, 'SPM.mat')

            SPM = sio.loadmat(spm_mat_file, squeeze_me=True, struct_as_record=False)['SPM']

            model_contrasts = SPM.xCon

        if DEBUG:
            print(' found {} contrasts'.format(len(model_contrasts)))

        con_map = {con.name: 'con_{:0>4}.nii'.format(cidx) for cidx,con in enumerate(model_contrasts,1) }


        if DEBUG:
            print('\tContrasts are:', con_map)
        
        for model_con, proj_con in cmap.items():
            
            path_to_con = os.path.join(first_level_dir, con_map[proj_con])
            
            if os.path.exists(path_to_con):
                contrast_lists[model_con].append(path_to_con)
            
    return contrast_lists


# ## Define nodes

# %%


# Infosource - a function free node to iterate over the list of subject names
l2_infosource = pe.Node(util.IdentityInterface(fields=['contrast_id']),
                  name="infosource")

smoothing_kernels = [ 8 ]
resolutions = ['medium']

resolution_and_kernel_list = product(resolutions, smoothing_kernels)


l2_infosource.iterables = [('contrast_id', l2_contrast_list), 
                           ('resolution_and_smoothing', resolution_and_kernel_list)
                        ]


# %%


# SelectFiles - to grab the data (alternativ to DataGrabber)

subject_pattern='*'
OUTPUT_DIR = output_dir
l2_output_dir = output_dir

l2_templates = {'cons': os.path.join(output_dir, MODEL_NAME, subject_pattern, '{smoothing_ksize}',
                         '{contrast_id}.nii')}

l2_selectfiles = pe.Node(nio.SelectFiles(l2_templates,
                               base_directory=OUTPUT_DIR,
                               sort_filelist=True),
                   name="selectfiles")


# %%


def make_contrast_list(model_path, cname,exclude_subjects, sample_perc=80):
    #EDITED BY CHRISTIN to get randomly sample a given percentage of subjects for second-level model

    import json
    import random
    import os
    import scipy.io as sio
    
    ROOT_DIR = '/data00/projects/megameta'
    
    def process_project(project_name, model_def, scan_all_subjs=False, DEBUG=False):

        project_spec = [pspec for pspec in model_def['Projects'] if pspec['Name']==project_name]

        if not project_spec:
            print('Cannot find specification for project: ', project_name)
            return None

        model_name = project_spec[0]['Model']
        cmap = project_spec[0]['ContrastMap']


        model_dir = os.path.join(ROOT_DIR, project_name, 
                                 "derivatives", "nipype",
                                 "model_{}".format(model_name)
                                )

        if not os.path.exists(model_dir):
            print('Cannot find first level model directory:', model_dir)
            return None

        subjs_with_models = [s for s in os.listdir(model_dir) if s.startswith('sub-')]
        #Exclude people
        subjs_with_models=[s for s in subjs_with_models if s not in exclude_subjects]
        
        #Get a random sample of participants (based on a percentage)
        sample_size=(sample_perc/100)*len(subjs_with_models)
        subj_list=random.sample(subjs_with_models,int(sample_size))
        
        print('Project: {}, Sampling {} of {} participants with a model'.format(project_name, int(sample_size), len(subjs_with_models)))
        
        if DEBUG:
            print("Found {} first level subject models\n".format(len(subjs_with_models)))


        contrast_lists = { cname: [] for cname in cmap}


        model_contrasts=None
        for sidx,subj in enumerate(subj_list):

            if DEBUG:
                print('Processing',subj, '-',end='')

            first_level_dir = os.path.join(model_dir, subj, 'medium', 'fwhm_8')

            if scan_all_subjs or sidx==0:
                spm_mat_file = os.path.join(first_level_dir, 'SPM.mat')

                SPM = sio.loadmat(spm_mat_file, squeeze_me=True, struct_as_record=False)['SPM']

                model_contrasts = SPM.xCon

            if DEBUG:
                print(' found {} contrasts'.format(len(model_contrasts)))

            con_map = {con.name: 'con_{:0>4}.nii'.format(cidx) for cidx,con in enumerate(model_contrasts,1) }


            if DEBUG:
                print('\tContrasts are:', con_map)

            for model_con, proj_con in cmap.items():

                path_to_con = os.path.join(first_level_dir, con_map[proj_con])

                if os.path.exists(path_to_con):
                    contrast_lists[model_con].append(path_to_con)

        return contrast_lists

    with open(model_path) as fh:
        model_def = json.load(fh)
        
    conlist=[]
    for p in model_def['Projects']:
        print(p)
        conlist.extend(process_project(p['Name'], model_def)[cname])
        
    return conlist


# %%


l2_getcontrasts = pe.Node(util.Function(input_names=['model_path','cname','exclude_subjects'],
                                     output_names=['contrasts', 'covariates'],
                                      function=make_contrast_list),
                    name='makecontrasts')
MDIR = os.path.abspath('../model_specifications')
l2_getcontrasts.inputs.model_path=os.path.join(MDIR, MODEL_SPEC_FILE)
l2_getcontrasts.inputs.cname=CONTRAST_NAME
l2_getcontrasts.inputs.exclude_subjects=exclude_subjects


# %%


#EDITED BY CHRISTIN (ADDING DATASINK)
# Datasink - creates output folder for important outputs
datasink = pe.Node(nio.DataSink(base_directory=OUTPUT_DIR,
                         container=l2_output_dir),
                name="datasink")

# Use the following DataSink output substitutions
substitutions = [('_contrast_id_', '')]
datasink.inputs.substitutions = substitutions


# ## Model nodes

# %%


osttdesign = pe.Node(spm.model.OneSampleTTestDesign(),
                         name="osttdesign")

osttdesign.inputs.explicit_mask_file='/data00/tools/spm8/apriori/brainmask_th25.nii'
osttdesign.inputs.threshold_mask_none=True


# %%


#MODEL_SPEC_FILE = 'group_mreg_behav_nonavers.json'
#CONTRAST_NAME='puremessage'
#PATH_TO_SPM_FOLDER = '/data00/tools/spm12mega/'
#JSON_MODEL_FILE = os.path.join('/data00/projects/megameta/scripts/jupyter_megameta/second_level_models',
#                               'model_specifications',
#                               MODEL_SPEC_FILE)
#exclude_subjects=[]


# %%


#EDITED BY CHRISTIN TO IMPPLEMENT MREG

# Multiple Regression Design - creates mreg Design
mregdesign = pe.Node(spm.model.MultipleRegressionDesign(),
                         name="mregdesign")
# Add covariates
## Make a list of covariates based on the contrast list
covs=[]
contrast_list, subj_list=make_contrast_list(JSON_MODEL_FILE,CONTRAST_NAME,exclude_subjects)[0]
pjs=[c.split('/')[4] for c in contrast_list]
pjs=[s for s in set(pjs)]
print(pjs)
print(subj_list)

## Make dummy variables based on list of projects and add them to the covariate list of dictionaries
#for pj in set(pjs):
#    cur_cov_vector=[]
#    for idx, _ in enumerate(pjs):
##            if pjs[idx]==pj:
 #               cur_cov_vector.append(1)
 #           else:
 #               cur_cov_vector.append(0)
 #   #make dictionary for current covariate
 #   cur_dict={'name': pj, 'vector': cur_cov_vector}
 #   #add dictionary to list of covs
 #   covs.append(cur_dict)

##NOTE: THE CODE ABOVE CREATES ONE DUMMY PER PROJECT. NEED TO TAKE ONE OUT AND DECIDE WHICH PROJECT TO USE AS COMPARISON/REFERENCE. 
#BELOW ARE TWO VERSIONS OF DOING THAT. VERSIN 1 RANDOMLY CHOOSES (# OF PROJECTS)-1 COVARIATES TO INCLUDE - BUT WE PROBABLY WANT TO BE MORE STRATEGIC
#VERSION 1
#covs=random.sample(covs,(len(pjs)-1))
# VERSION 2 REMOVES DARPA1 TO MAKE IT THE REFERENCE PROJECT -- BUT I DON'T HAVE A CLEAR RATIONALE FOR WHY THAT OVER OTHERS RIGHT NOW...
#covs=[i for i in covs if i['name']!='darpa1']

# Intended covs format:
# covs = [
#    {'name':'alcohol','vector': []},
#    {'name':'darpa1','vector': []},
#    {'name':'darpa2','vector': []},
#    {'name':'cityyear','vector': []},
#    {'name':'project1','vector': []}
#]

# Add covariate of behaivor change and baseline
subj_list=[]
for pj in pjs:
    project_spec = [pspec for pspec in model_def['Projects'] if pspec['Name']==pj]

    

    model_name = project_spec[0]['Model']
    model_dir = os.path.join(ROOT_DIR, pj, 
                            "derivatives", "nipype",
                            "model_{}".format(model_name)
                        )
    subjs_with_models = [s for s in os.listdir(model_dir) if s.startswith('sub-')]
    #Exclude people
    subjs_with_models=[s for s in subjs_with_models if s not in exclude_subjects]
    subj_list=subj_list+subjs_with_models
    
subj_list=[s.replace('sub-','') for s in subj_list]

##make a new behavior vector for the people who are in subj_list
regressors=pd.read_csv('/data00/projects/megameta/scripts/jupyter_megameta/second_level_models/indbehav_data/behaviorchange_050919nc.csv')
behav_mreg=[]
for row_num, val in enumerate(regressors['change']):
    if regressors['pID'][row_num] in subj_list:
        behav_mreg.append(regressors['change'][row_num])

behav_mreg_dict={'name': 'behav_mreg', 'vector':behav_mreg}

behav_baseline=[]
for row_num, val in enumerate(regressors['baseline']):
    if regressors['pID'][row_num] in subj_list:
        behav_baseline.append(regressors['baseline'][row_num])

behav_baseline_dict={'name': 'behav_baseline', 'vector':behav_baseline}

covs=[behav_mreg_dict,behav_baseline_dict]

mregdesign.inputs.covariates=covs

mregdesign.inputs.explicit_mask_file='/data00/tools/spm8/apriori/brainmask_th25.nii'


# %%


# EstimateModel - estimate the parameters of the model
level2estimate = pe.Node(spm.model.EstimateModel(estimation_method={'Classical': 1}),
                      name="level2estimate")



# %%


# EstimateContrast - estimates simple group contrast
level2conestimate = pe.Node(spm.model.EstimateContrast(group_contrast=True),
                         name="level2conestimate")


# %%


'''
cont1 = ['QuitIntent', 'T', ['QuitIntent', 'FTND', 'mean_WC', 'mean'], [1, 0, 0, 0]]
cont2 = ['FTND', 'T', ['QuitIntent', 'FTND', 'mean_WC', 'mean'], [0, 1, 0, 0]]
cont3 = ['mean_WC', 'T', ['QuitIntent', 'FTND', 'mean_WC', 'mean'], [0, 0, 1, 0]]
cont4 = ['mean', 'T', ['QuitIntent', 'FTND', 'mean_WC', 'mean'], [0, 0, 0, 1]]
'''

cont = ['behav_mreg', 'T', ['behav_mreg','behav_baseline'], [1,0]]

level2conestimate.inputs.contrasts = [cont]


# ## Setup second level workflow

# %%


#l2_working_dir = os.path.join(PROJECT_DIR, 'nipype', 'workingdir_banner_2nd_level')
l2_working_dir = working_dir


# %%


# EDITED BY CHRISTIN (adding datasink to the workflow)
l2analysis = pe.Workflow(name='l2analysis')

l2analysis.base_dir = l2_working_dir

# Connect up the 2nd-level analysis components
l2analysis.connect(
                    [
                        
                    #(l2_infosource, l2_getcontrasts, [('contrast_id', 'contrast_id'),
                     #                                ('model_path')]),
                     
                     (l2_getcontrasts,  mregdesign, [('contrasts', 'in_files')]),
                     
                    (mregdesign, level2estimate, [('spm_mat_file',
                                                          'spm_mat_file')] ),
                    (level2estimate, level2conestimate, [('spm_mat_file',
                                                          'spm_mat_file'),
                                                         ('beta_images',
                                                          'beta_images'),
                                                         ('residual_image',
                                                          'residual_image')]),
                    (level2conestimate, datasink, [('spm_mat_file',
                                                    'contrasts.@spm_mat'),
                                                   ('spmT_images',
                                                    'contrasts.@T'),
                                                   ('con_images',
                                                    'contrasts.@con')])
                    ])

