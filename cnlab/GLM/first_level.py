# coding: utf-8

# # Building first level models using _nipype_ and _SPM12_
# 
# ## Base functionality for _megameta_ project
# 
# -------
# #### History
# * 6/16/20 mbod - saved as a .py script for refactoring
# * 5/29/20 cscholz - include ModelAsItem functionality 
# * 5/18/20 hychan - include option for throwing dummy scans (`ExcludeDummyScans`) and excluding sessions (`ExcludeRuns`)
# * 9/18/19 hychan - include option for using custom event TSV
# * 4/9/19 cscholz - made small correction to make_contrast_list() (setting: -1/neg_length instead of -1/pos_length)
# * 4/2/19 mbod - split out processing pipeline for revised workflow
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

from nipype.interfaces.nipy.preprocess import Trim

from itertools import combinations

from nilearn import plotting, image
from nistats import thresholding



import scipy.io as sio
import numpy as np
import json
import pandas as pd
import glob


PATH_TO_SPM_FOLDER = '/data00/tools/spm12mega'
BRAIN_MASK_PATH = '/data00/tools/spm8/apriori/brainmask_th25.nii'



def setup_pipeline(MODEL_PATH, 
                   exclude_subjects=None,
                   include_subjects=None,
                   DEBUG=False):

    # Set the way matlab should be called
    mlab.MatlabCommand.set_default_matlab_cmd("matlab -nodesktop -nosplash")
    # If SPM is not in your MATLAB path you should add it here
    mlab.MatlabCommand.set_default_paths(PATH_TO_SPM_FOLDER)


    # ### Parameters
    # 
    # * These need to be reformatted to be consistent
    # * as data is not smoothed commented out the `fwhm_size` param - but data probably has a value

    # #### Load JSON model config


    with open(MODEL_PATH) as fh:
        model_def = json.load(fh)



    TASK_NAME = model_def['TaskName']
    RUNS = model_def['Runs']
    MODEL_NAME = model_def['ModelName']
    PROJECT_NAME = model_def['ProjectID']
    BASE_DIR = model_def['BaseDirectory']

    TR = model_def['TR']


    # TODO - add configurable project paths
    PROJECT_DIR = os.path.join(BASE_DIR, PROJECT_NAME)
    SUBJ_DIR = os.path.join(PROJECT_DIR, 'derivatives', 'nipype', 'resampled_and_smoothed')

    task_func_template = "sr{PID}_task-{TASK}_run-*_*preproc*.nii"

    subject_list = [subj for subj in os.listdir(SUBJ_DIR) 
                       if glob.glob(os.path.join(SUBJ_DIR,subj,'medium', 'fwhm_8',
                                            task_func_template.format(PID=subj, TASK=TASK_NAME)))
                    ]

    output_dir = os.path.join(PROJECT_DIR,'derivatives', 'nipype','model_{}_{}'.format(TASK_NAME.upper(), MODEL_NAME))        # name of 1st-level output folder
    working_dir = os.path.join(PROJECT_DIR, 'working', 
                               'nipype', 'workingdir_model_{}_{}'.format(TASK_NAME.upper(), MODEL_NAME))   # name of 1st-level working directory



    # check to see if output and work directories exist

    if not os.path.exists(output_dir):
        os.makedirs(output_dir) 

    if not os.path.exists(working_dir):
        os.makedirs(working_dir)



    try:
        subject_list = [ s for s in subject_list if s not in exclude_subjects ]
        
        if DEBUG:
            print('\n\nApplied subject exclusion list:\n\t',' '.join(exclude_subjects))
    except:
        if DEBUG:
            print('\n\nNo subject exclusions applied')

    try:
        subject_list = [ s for s in subject_list if s in include_subjects ]
        if DEBUG:
            print('\n\nApplied subject inclusion list:\n\t',' '.join(include_subjects))
    except:
        if DEBUG:
            print('\n\nNo subject inclusions applied')

    if DEBUG:
        print('\n\nSUBJECT LIST IS:\n\t', ' '.join(subject_list))


    # add directorys and subject list to dictionary
    # before return
    
    model_def['subject_list']=subject_list
    model_def['output_dir']=output_dir
    model_def['working_dir']=working_dir
    model_def['model_path']=MODEL_PATH
    model_def['SUBJ_DIR']=SUBJ_DIR
    model_def['PROJECT_DIR']=PROJECT_DIR
    
    return model_def

    
# ### Utility functions for subject info and contrasts    

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
    import glob
    
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
                    #df2[pmod_name]+=0.001
                    pmod.append(None)
                    continue
                    
                    
                    
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
    TASK_RUNS = model_def['Runs'].copy()
    MODEL_NAME = model_def['ModelName']
    PROJECT_ID = model_def['ProjectID']
    BASE_DIR = model_def['BaseDirectory']
    
    # Collect ModelAsItem conditions across all runs
    ALL_conditions=[]
    condition_names = list(model_def['Conditions'].keys())
   
    PROJECT_DIR = os.path.join(BASE_DIR, PROJECT_ID)
    FMRIPREP_SUBJ_DIR = os.path.join(PROJECT_DIR,'derivatives', 'fmriprep')
    BATCH8_SUBJ_DIR=os.path.join(PROJECT_DIR,'derivatives','batch8')
    
    # Determine which subject dir is the right one for this project
    if glob.glob(BATCH8_SUBJ_DIR):
        SUBJ_DIR=BATCH8_SUBJ_DIR
    else:
        SUBJ_DIR=FMRIPREP_SUBJ_DIR
    
    realign_files = []
    subject_info = []
    
    if model_def.get('CustomEventDir'):
        EVENT_DIR = model_def['CustomEventDir'].format(subject_id=subject_id)
    else:
        EVENT_DIR = os.path.join(SUBJ_DIR, subject_id, 'func')
        
        
    if model_def.get('CustomEventTemplate'):
        event_template=model_def['CustomEventTemplate']
    else:
        event_template='{SUBJ}_task-{TASK}_run-0{RUN}_events.tsv'
    
    # check to see which runs exist for subject
    # by looking for appropriate events.tsv files
    # this could (should?) also include looking for the nifti file?
    
    if model_def.get('ExcludeDummyScans'):
        ExcludeDummyScans = model_def['ExcludeDummyScans']
        TR = model_def['TR']

    else:
        ExcludeDummyScans = 0
        TR = 0

    if model_def.get('ExcludeRuns'):
        ExcludeRuns = model_def['ExcludeRuns']
        for er in ExcludeRuns:
            if er['subject'] == subject_id and er['run'] in TASK_RUNS:
                TASK_RUNS.remove(er['run'])
                    
    runs_for_subj = [run for run in TASK_RUNS
                    if os.path.exists(os.path.join(EVENT_DIR,
                                   event_template.format(SUBJ=subject_id,
                                                         TASK=TASK_NAME,
                                                         RUN=run)))
                    ]
    
    if DEBUG:
        print("runs_for_subj", runs_for_subj)
        print("checked paths:")
        for run in TASK_RUNS:
            print('\t', os.path.join(EVENT_DIR,
                                    event_template.format(SUBJ=subject_id,
                                                         TASK=TASK_NAME,
                                                         RUN=run)))
        print("TASK NAME", TASK_NAME)
        print("pmod", pmod)
        print("TASK_RUNS", TASK_RUNS)
        print("subject_id", subject_id)
    
    for _, run_num in enumerate(runs_for_subj,1):
        
        condition_names = list(model_def['Conditions'].keys())
    
        events_df = pd.read_csv(os.path.join(EVENT_DIR,
                                             event_template.format(SUBJ=subject_id,
                                                         TASK=TASK_NAME,
                                                         RUN=run_num)), sep='\t')

            
        events_df['onset'] = events_df['onset'] - (TR * ExcludeDummyScans)
        
        # cscholz 29/5/20 - adding ModelAsItem functionality written by mbod
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

        confound_file=glob.glob(os.path.join(SUBJ_DIR, subject_id, 'func',
                                   '{}_task-{}_run-0{}_*confounds*.tsv'.format(subject_id, 
                                                                               TASK_NAME,
                                                                               run_num)))[0]

        confound_df = pd.read_csv(confound_file, sep='\t')
        confound_df = confound_df[ExcludeDummyScans:]
        
        if 'TransX' in confound_df.columns:
            cols_to_use = [ 'TransX','TransY', 'TransZ', 'RotX', 'RotY', 'RotZ']
        else:
            cols_to_use = [ 'trans_x','trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z']

        confound_df[cols_to_use].to_csv(realign_file, 
                                        header=False, 
                                        index=False,
                                        sep='\t')

        realign_files.append(realign_file)

        onsets = []
        dur = []
        for cond in condition_names:
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
    
        # add run-specific conditions to ALL conditions
        ALL_conditions= ALL_conditions + condition_names
        
    
    ALL_conditions=set(ALL_conditions)
    
    DM_regressors = []
    for cond in ALL_conditions:
        DM_regressors.append(cond)
        if pmod and model_def['Modulators'].get(cond):
            DM_regressors.append('{}x{}^1'.format(cond, model_def['Modulators'].get(cond)))
    
    
    return subject_info, realign_files, DM_regressors





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
                ccode.append(-1/neg_length)
            else:
                ccode.append(0)

        cont.append([cname, 'T', condition_names, ccode])
        
                
        if DEBUG:
            print(contrast)
            print(ccode)
            
    return cont


# ## Set up processing nodes for modeling workflow


def build_pipeline(model_def):
    
    
    # create pointers to needed values from
    # the model dictionary
    # TODO - this could be refactored
    TR = model_def['TR']
    subject_list = model_def['subject_list']
    JSON_MODEL_FILE = model_def['model_path']

    working_dir = model_def['working_dir']
    output_dir = model_def['output_dir']
    
    SUBJ_DIR = model_def['SUBJ_DIR']
    PROJECT_DIR = model_def['PROJECT_DIR']
    TASK_NAME = model_def['TaskName']
    RUNS = model_def['Runs']
    MODEL_NAME = model_def['ModelName']
    PROJECT_NAME = model_def['ProjectID']
    BASE_DIR = model_def['BaseDirectory']
    
    SERIAL_CORRELATIONS = "AR(1)" if not model_def.get('SerialCorrelations') else model_def.get('SerialCorrelations')
    RESIDUALS = model_def.get('GenerateResiduals')

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




    # Level1Design - Generates an SPM design matrix
    level1design = pe.Node(spm.Level1Design(bases={'hrf': {'derivs': [0, 0]}},
                                     timing_units='secs',
                                     interscan_interval=TR,
                                     # model_serial_correlations='AR(1)', # [none|AR(1)|FAST]',
                                     # 8/21/20 mbod - allow for value to be set in JSON model spec
                                     model_serial_correlations=SERIAL_CORRELATIONS,   

                                     # TODO - allow for specified masks
                                     mask_image = BRAIN_MASK_PATH,
                                     global_intensity_normalization='none'
                                           ),
                        name="level1design")


    # #### Estimate Model node
    # EstimateModel - estimate the parameters of the model
    level1estimate = pe.Node(spm.EstimateModel(estimation_method={'Classical': 1}),
                          # 8/21/20 mbod - allow for value to be set in JSON model spec
                          write_residuals=RESIDUALS,
                             
                          name="level1estimate")


    # #### Estimate Contrasts node
    # EstimateContrast - estimates contrasts
    conestimate = pe.Node(spm.EstimateContrast(), name="conestimate")


    # ## Setup pipeline workflow for level 1 model
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


 
    # ## Set up nodes for file handling and subject selection
    # ### `getsubjectinfo` node 
    # 
    # * Use `get_subject_info()` function to generate spec data structure for first level model design matrix

 
    # Get Subject Info - get subject specific condition information
    getsubjectinfo = pe.Node(util.Function(input_names=['subject_id', 'model_path'],
                                   output_names=['subject_info', 'realign_params', 'condition_names'],
                                   function=get_subject_info),
                          name='getsubjectinfo')



    makecontrasts = pe.Node(util.Function(input_names=['subject_id', 'condition_names', 'model_path'],
                                         output_names=['contrasts'],
                                          function=make_contrast_list),
                        name='makecontrasts')




    if model_def.get('ExcludeDummyScans'):
        ExcludeDummyScans = model_def['ExcludeDummyScans']
    else:
        ExcludeDummyScans = 0

    print(f'Excluding {ExcludeDummyScans} dummy scans.')
    trimdummyscans = pe.MapNode(Trim(begin_index=ExcludeDummyScans),
                          name='trimdummyscans',
                          iterfield=['in_file'])


    # ### `infosource` node
    # 
    # * iterate over list of subject ids and generate subject ids and produce list of contrasts for subsequent nodes


    # Infosource - a function free node to iterate over the list of subject names
    infosource = pe.Node(util.IdentityInterface(fields=['subject_id', 'model_path', 'resolution', 'smoothing']
                                       ),
                      name="infosource")

    try:
        fwhm_list = model_def['smoothing_list']
    except:
        fwhm_list = [4,6,8]

    try:
        resolution_list = model_def['resolutions']
    except:
        resolution_list = ['low','medium','high']

    infosource.iterables = [('subject_id', subject_list),
                            ('model_path', [JSON_MODEL_FILE]*len(subject_list)),
                            ('resolution', resolution_list),
                            ('smoothing', ['fwhm_{}'.format(s) for s in fwhm_list])
                           ]




    # SelectFiles - to grab the data (alternativ to DataGrabber)

    ## TODO: here need to figure out how to incorporate the run number and task name in call
    templates = {'func': '{subject_id}/{resolution}/{smoothing}/sr{subject_id}_task-'+TASK_NAME+'_run-0*_*MNI*preproc*.nii'}          


    selectfiles = pe.Node(nio.SelectFiles(templates,
                                   base_directory='{}/{}/derivatives/nipype/resampled_and_smoothed'.format(
                                       BASE_DIR,
                                       PROJECT_NAME)),
                          working_dir=working_dir,
                       name="selectfiles")


    # ### Specify datasink node
    # 
    # * copy files to keep from various working folders to output folder for model for subject

    # Datasink - creates output folder for important outputs
    datasink = pe.Node(nio.DataSink(base_directory=SUBJ_DIR,
                             parameterization=True, 
                             #container=output_dir      
                                   ),
                    name="datasink")

    datasink.inputs.base_directory = output_dir


    # Use the following DataSink output substitutions
    substitutions = []
    subjFolders = [('_model_path.*resolution_(low|medium|high)_smoothing_(fwhm_\\d{1,2})_subject_id_sub-.*/(.*)$', '\\1/\\2/\\3')]
    substitutions.extend(subjFolders)
    datasink.inputs.regexp_substitutions = substitutions

    
     # datasink connections
    
    datasink_in_outs = [('conestimate.spm_mat_file','@spm'),
                        ('level1estimate.beta_images','@betas'),
                        ('level1estimate.mask_image','@mask'),
                        ('conestimate.spmT_images','@spmT'),
                        ('conestimate.con_images','@con'),
                        ('conestimate.spmF_images','@spmF')
                       ]
    
    if model_def.get('GenerateResiduals'):
        datasink_in_outs.append(('level1estimate.residual_images','@residuals'))   
    
    

    # ---------

    # ## Set up workflow for whole process




    pipeline = pe.Workflow(name='first_level_model_{}_{}'.format(TASK_NAME.upper(),MODEL_NAME))
    pipeline.base_dir = os.path.join(SUBJ_DIR, working_dir)


    pipeline.connect([(infosource, selectfiles, [('subject_id', 'subject_id'),
                                                 ('resolution', 'resolution'),
                                                 ('smoothing', 'smoothing')
                                                ]),
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


    #                  (selectfiles, l1analysis, [('func',
    #                                          'modelspec.functional_runs')]),
                      (selectfiles, trimdummyscans, [('func',
                                              'in_file')]),

                      (trimdummyscans, l1analysis, [('out_file',
                                              'modelspec.functional_runs')]),


                      
                     
                     
                     (infosource, datasink, [('subject_id','container')]),
                      
                     (l1analysis, datasink, datasink_in_outs)
                    ]
                     
    )

    
    return pipeline








