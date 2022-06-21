
# coding: utf-8

# # Building first level models using `nipype` and `SPM` with parametric modulation
# 
# ## Parametric modulator setup from BIDS events tsv for _megameta_ tasks
# 
# ### Multiple run task setup (testing on P1 Image Task)
# 
# -------
# #### History
# 
# * 1/21/19 mbod  - modify notebook from P1 Banner task with pmod for megameta project
# * 1/22/19 mbod  - adjust functions for pmod and modelspec for two runs-single model for image_task
# * 2/2/19 mbod   - incorporate JSON model spec with `get_subject_info()` function
# * 2/6/19 mbod   - test with likeme rating pmod
# * 2/19/19 mbod  - extract code into a module for cleaner import
# -----
# 
# ### Description
# 
# * Set up a nipype workflow to use SPM12 to make first level models for _Project 1_ banner task data (preprocessed using `batch8` SPM8 scripts) in BIDS derivative format
# 
# * Workflow has steps (for each subject in subject list):
#     1. get subject NIFTI file
#     2. preparation steps
#         a. gunzip
#         b. resample to specified resolution (downsample or upsample)
#         c. smooth using specific kernel size(s)
#     3. call the `get_subject_info()` function to get spec data and the realignment confounds for run
#     4. create design matrix
#     5. estimate model (output betas)
#     6. estimate conditions (output conn files)
#     7. estimate constrasts (output spmT)
#     
#     
# 
# * **NOTES**
#  
#     * Mask specified in the model is `/data00/tools/spm8/apriori/brainmask_th25.nii` 
#     * Crashes of workflow produce a lot of logged data and usually a path to a `.pklz` file. In a terminal you can view this with:
#         * `nipypecli crash <PATH TO LOG>`
#         * e.g.
#         
#         ```
#         nipypecli crash /fmriNASTest/data00/projects/drisc/scripts/jupyter/GLM_models/crash-20180914-214850-mbod%40asc.upenn.edu-selectfiles.a0-f01ed26e-9f10-4520-8df5-94e2195b02a0.pklz        
#         
#         ```

# ### Setup
# 
# * import required modules and define parameters


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
import matplotlib.pyplot as plt

import scipy.io as sio
import numpy as np
import json
import pandas as pd


DEBUG=True

JSON_MODEL_PATH=None

# #### Matlab path
# 


# Set the way matlab should be called
mlab.MatlabCommand.set_default_matlab_cmd("matlab -nodesktop -nosplash")
# If SPM is not in your MATLAB path you should add it here
mlab.MatlabCommand.set_default_paths('/fmriNASTest/data00/tools/spm12')


# ### Parameters
# 
# * These need to be reformatted to be consistent
# * as data is not smoothed commented out the `fwhm_size` param - but data probably has a value

# #### Load JSON model config


def load_model_config(model_path):
    '''
    load the JSON model file and return dictionary
    '''
    JSON_MODEL_PATH=model_path
    
    with open(model_path) as fh:
        model_def = json.load(fh)
    return model_def



def setup_folders(output_dir, working_dir):
    '''
    check to see if output and work directories exist and create if needed
    '''

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if not os.path.exists(working_dir):
        os.makedirs(working_dir)


# ### Utility functions for subject info and contrasts

# ### Setup design matrix data for subject
# 
# * need a function to set up the nipype `Bunch` format used
#     * https://nipype.readthedocs.io/en/latest/users/model_specification.html
# * read the onsets/dur/conditions from task logs and extract needed data

def get_subject_info(subject_id, model_path):
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
    
    
    def make_pmod(df, conditions, pmods={}):
        
        pmod = []
        
        for cond in conditions:
        
            
            if not pmods.get(cond):
                pmod.append(None)
            else:
                df2 = df[df.trial_type==cond]
                
                pmod_name = pmods.get(cond)
                
                #pmod = [pmod] if not type(pmods) is list else pmod
                
                if df2[pmod_name].var()==0:
                    df2[pmod_name]+=0.001

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
            spec_trials = model['Conditions'][con]
            spec.loc[spec.trial_type.isin(spec_trials),'trial_type'] = con
            spec.onset.sort_values()

        cols_to_drop=[]
        for pmidx, (cond, colname) in enumerate(model['Modulators'].items(),1):
            #pmod="pmod{}_{}".format(pmidx,colname)
            pmod=colname

            spec.loc[spec.trial_type==cond,pmod] = spec[colname]
            cols_to_drop.append(colname)
        #spec.drop(columns=cols_to_drop,inplace=True)

        return spec
    
    
    with open(model_path) as fh:
        model_def = json.load(fh)
    
    
    PROJECT_DIR = '/data00/projects/megameta/project1'
    SUBJ_DIR = os.path.join(PROJECT_DIR,'derivatives', 'batch8')
    
    realign_files = []
    subject_info = []
    
    TASK_NAME = model_def['TaskName']
    TASK_RUNS = model_def['Runs']
    MODEL_NAME = model_def['ModelName']
    
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
    
    for run_num, _ in enumerate(runs_for_subj,1):
    
    
        events_df = pd.read_csv(os.path.join(SUBJ_DIR, subject_id, 'func',
                                             '{}_task-{}_run-0{}_events.tsv'.format(subject_id, 
                                                                                    TASK_NAME,
                                                                                    run_num)),
                               sep='\t')


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
        for cond in model_def['Conditions']:
            onsets.append(onsets_df[onsets_df.trial_type==cond].onset.values)
            dur.append(onsets_df[onsets_df.trial_type==cond].duration.values)

            
        #condition_names = list(onsets_df.trial_type.unique())
        condition_names = list(model_def['Conditions'].keys())
            
        #pmod = make_pmod(rdf, condition_names)   
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
        if model_def['Modulators'].get(cond):
            DM_regressors.append('{}x{}^1'.format(cond, model_def['Modulators'].get(cond)))
    
    
    return subject_info, realign_files, DM_regressors


# ### Set up contrasts
# 
# * This part of the template needs work to provide a cleaner way to specify contrasts
# * Could use the same vector contrasts approach as we have in batch8 and then have a function to convert this into the list of list data structure that nipype spm contrasts node looks for


def make_contrast_list(subject_id, condition_names):
    condition_names.append('constant')

    cont = []
    for idx, cname in enumerate(condition_names):
        ccode = [0 if pos!=idx else 1 for pos in range(len(condition_names))]
        cont.append([cname, 'T', condition_names, ccode])


    return cont



# ## Set up processing nodes for modeling workflow

def set_up_nodes(node_dict={}, TR=2.0):
    '''
    define needed nodes and return in a dictionary
    '''


    # #### Specify model node


    # SpecifyModel - Generates SPM-specific Model
    node_dict['modelspec'] = pe.Node(model.SpecifySPMModel(concatenate_runs=False,
                                     input_units='secs',
                                     output_units='secs',
                                     time_repetition=TR,
                                     high_pass_filter_cutoff=128),
                                     output_units = 'scans',
                                     name="modelspec")                                      

    if DEBUG:
        print('modelspec node defined\n\t', node_dict['modelspec'])

    # #### Level 1 Design node
    # 
    # ** TODO -- get the right matching template file for fmriprep **
    # 
    # * ??do we need a different mask than:
    # 
    #     `'/data00/tools/spm8/apriori/brainmask_th25.nii'`



    # Level1Design - Generates an SPM design matrix
    node_dict['level1design'] = pe.Node(spm.Level1Design(bases={'hrf': {'derivs': [0, 0]}},
                                     timing_units='secs',
                                     interscan_interval=TR,
                                     model_serial_correlations='none', #'AR(1)',
                                     mask_image = '/data00/tools/spm8/apriori/brainmask_th25.nii',
                                     global_intensity_normalization='none'
                                           ),
                        name="level1design")

    if DEBUG:
        print('level1design node defined\n\t', node_dict['level1design'])


    # EstimateModel - estimate the parameters of the model
    node_dict['level1estimate'] = pe.Node(spm.EstimateModel(estimation_method={'Classical': 1}),
                          name="level1estimate")

    if DEBUG:
        print('level1estimate node defined\n\t', node_dict['level1estimate'])


    # #### Estimate Contrasts node

    # EstimateContrast - estimates contrasts
    node_dict['conestimate'] = pe.Node(spm.EstimateContrast(), name="conestimate")

    if DEBUG:
        print('conestimate node defined\n\t', node_dict['conestimate'])

        
    # Initiation of the 1st-level analysis workflow
    l1analysis = pe.Workflow(name='l1analysis')

    # Connect up the 1st-level analysis components
    l1analysis.connect([(node_dict['modelspec'], node_dict['level1design'], [('session_info',
                                                    'session_info')]),
                        
                        (node_dict['level1design'], node_dict['level1estimate'], 
                                                 [('spm_mat_file', 'spm_mat_file')]),
                        
                        (node_dict['level1estimate'], node_dict['conestimate'], [('spm_mat_file',
                                                        'spm_mat_file'),
                                                       ('beta_images',
                                                        'beta_images'),
                                                       ('residual_image',
                                                        'residual_image')])

                        ])
    
    node_dict['l1analysis']=l1analysis
    
    if DEBUG:
        print('l1analysis workflow defined with:\n\tmodelspec\n\tlevel1design\n\tlevel1estimate\n\tconestimate\n\n\tnodes\n')

        
    # Get Subject Info - get subject specific condition information
    node_dict['getsubjectinfo'] = pe.Node(util.Function(input_names=['subject_id', 'model_path'],
                                   output_names=['subject_info', 'realign_params', 'condition_names'],
                                   function=get_subject_info),
                          name='getsubjectinfo')




    node_dict['makecontrasts'] = pe.Node(util.Function(input_names=['subject_id', 'condition_names'],
                                         output_names=['contrasts'],
                                          function=make_contrast_list),
                        name='makecontrasts')
        
        
    # Infosource - a function free node to iterate over the list of subject names
    node_dict['infosource'] = pe.Node(util.IdentityInterface(fields=['subject_id', 'model_path']
                                       ),
                      name="infosource")

    node_dict['infosource'].iterables = [('subject_id', subject_list),
                            ('model_path', [JSON_MODEL_FILE]*len(subject_list))

                           ]
        
        
        
    return node_dict



def OTHER_NODES():





    # ### `selectfiles` node
    # 
    # * match template to find source files (functional) for use in subsequent parts of pipeline

    # In[60]:

    # SelectFiles - to grab the data (alternativ to DataGrabber)


    ## TODO: here need to figure out how to incorporate the run number and task name in call
    templates = {'func': '{subject_id}/func/{subject_id}_task-image_run-0*_space-MNI152-T1-1mm_desc-preproc_bold.nii.gz'}          


    selectfiles = pe.Node(nio.SelectFiles(templates,
                                   base_directory='/data00/projects/megameta/project1/derivatives/batch8'),
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

    # In[61]:

    gunzip = pe.MapNode(Gunzip(),name="gunzip", iterfield=['in_file'])


    # #### Specify smoothing node

    # In[62]:

    smooth = pe.Node(interface=spm.Smooth(), name="smooth")
    fwhmlist = [4,6,8]
    smooth.iterables = ('fwhm', fwhmlist)


    # #### Specify resampling node

    # In[63]:

    resample = pe.Node(interface=spm.utils.ResliceToReference(), name='resample')


    # In[64]:

    resample.inputs.target = '/data00/projects/project1/data/subjs/P100/BOLD/argue_run01/wad_argue_run01.nii'


    # In[65]:

    unzip_resample_and_smooth = pe.Workflow(name='unzip_resample_and_smooth')

    unzip_resample_and_smooth.base_dir = os.path.join(SUBJ_DIR, working_dir)

    unzip_resample_and_smooth.connect(
        [
            #(gunzip, resample, [('out_file', 'in_files')]),
            #(resample, smooth, [('out_files', 'in_files')])
            (gunzip, smooth, [('out_file', 'in_files')])
        ]
    )


    # In[66]:

    Image(unzip_resample_and_smooth.write_graph(graph2use='flat', simple_form=False))


    # ### Specify datasink node
    # 
    # * copy files to keep from various working folders to output folder for model for subject



    # In[68]:

    # Datasink - creates output folder for important outputs
    datasink = pe.Node(nio.DataSink(base_directory=SUBJ_DIR,
                             parameterization=True, 
                             #container=output_dir      
                                   ),
                    name="datasink")

    datasink.inputs.base_directory = output_dir

    # Use the following DataSink output substitutions
    substitutions = []
    subjFolders = [('_subject_id_%s/_fwhm_%s' % (sub,f), 'fwhm_%s' % (f))
                   for f in fwhmlist
                   for sub in subject_list]
    substitutions.extend(subjFolders)
    datasink.inputs.substitutions = substitutions



def create_workflow(nd, workflow_name, working_dir):
    # ---------

    # ## Set up workflow for whole process


    workflow = pe.Workflow(name=workflow_name)
    workflow.base_dir = working_dir


    workflow.connect([(nd['infosource'], nd['selectfiles'], [('subject_id', 'subject_id')]),
                      (nd['infosource'], md['getsubjectinfo'], [('subject_id', 'subject_id'),
                                                    ('model_path', 'model_path')
                                                   ]),
                      (nd['infosource'], nd['makecontrasts'], [('subject_id', 'subject_id')]),
                      (nd['getsubjectinfo'], nd['makecontrasts'], [('condition_names', 'condition_names')]),

                     (nd['getsubjectinfo'], nd['l1analysis'], [('subject_info',
                                                     'modelspec.subject_info'),
                                                    ('realign_params',
                                                       'modelspec.realignment_parameters')]),
                      (nd['makecontrasts'], nd['l1analysis'], [('contrasts',
                                                 'conestimate.contrasts')]),


                      (nd['aselectfiles'], nd['gunzip'], [('func','in_file')]),

                      #(gunzip, resample, [('out_file', 'in_files')]),

                      #(resample, smooth, [('out_files', 'in_files')]),

                      (gunzip, smooth, [('out_file', 'in_files')]),


                      (smooth, l1analysis, [('smoothed_files',
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


    return workflow

    

def show_design_matrix(spm_mat_path):
    '''
    '''
    

    spm_mat=sio.loadmat(spm_mat_path,
               struct_as_record=False, squeeze_me=True)


    regressor_labels=spm_mat['SPM'].xX.name


    # use nistats.reporting plot function and a pandas data frame




