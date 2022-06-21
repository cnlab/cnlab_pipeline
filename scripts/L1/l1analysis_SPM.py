# coding: utf-8

import os, sys, json, copy, argparse  

import nipype.interfaces.io as nio           
import nipype.interfaces.spm as spm    
import nipype.interfaces.fsl as fsl
import nipype.interfaces.matlab as mlab     
import nipype.pipeline.engine as pe 
import nipype.algorithms.modelgen as model
from nipype.interfaces.base import Bunch
from nipype.algorithms.misc import Gunzip
from nipype.interfaces.utility import IdentityInterface

import pandas as pd
import numpy as np

import nibabel as nib

def ensure_list(obj):
    if type(obj) is not list:
        return [obj]
    else:
        return obj

def make_model(model_info, **env):
    
    modelspec_args = copy.deepcopy(model_info)
    
    subject_info = []
    
    for event_file, regressor_file, func_file in zip(modelspec_args['event_files'], modelspec_args['regressors'], modelspec_args['functional_runs']):
        
        func_img = nib.load(os.path.join(env['data_path'], func_file))
        tr_max = func_img.header.get_data_shape()[-1] * modelspec_args['time_repetition']
        
        event = (pd.read_csv(os.path.join(env['data_path'], event_file), sep='\t')
                 .query('onset >= 0')
                 .query(f'onset <= {tr_max}')
                 .sort_values(by='onset')
                 .reset_index(drop=True))

        bunch = { 'conditions': [],
                  'onsets': [],
                  'durations': [] }

        for trial_type, trial_df in event.groupby('trial_type'):
            bunch['conditions'].append(trial_type)
            bunch['onsets'].append(trial_df['onset'].tolist())
            bunch['durations'].append(trial_df['duration'].tolist())
            
        if regressor_file.endswith('txt'):
            regressors = np.loadtxt(os.path.join(env['data_path'], regressor_file)).T.tolist()
            regressor_names = [f"R{i}" for i in range(regressors.shape[1])]
        else:
            regressors_array = pd.read_csv(os.path.join(env['data_path'], regressor_file), sep='\t')
            regressor_names = modelspec_args.get('regressor_names', regressors_array.columns.tolist())
            
            regressors = regressors_array[regressor_names].fillna(0).values.T.tolist()
        
        bunch['regressors'] = regressors
        bunch['regressor_names'] = regressor_names
        
        if modelspec_args.get('pmod'):
            
            bunch['pmod'] = []
            
            for condition in bunch['conditions']:
                pmod_args = None
                
                if modelspec_args['pmod'].get(condition):
                    pmod_args = { 'name': [], 'param': [], 'poly': [] }
                    
                    pmod_variables = ensure_list(modelspec_args['pmod'][condition])

                    for pmod_variable in pmod_variables:
                        pmod_param = event.query(f'trial_type == "{condition}"')[pmod_variable]
                        
                        # Replace missing with mean and check variation
                        pmod_param = pmod_param.fillna(pmod_param.mean())
                        if pmod_param.var() == 0:
                            continue
                        
                        pmod_args['name'].append(pmod_variable)
                        pmod_args['param'].append(pmod_param.values.tolist())
                        pmod_args['poly'].append(1)
                        
                    pmod_args = Bunch(**pmod_args)
                        
                bunch['pmod'].append(pmod_args)       
        
        subject_info.append(Bunch(**bunch))
        
    modelspec_args['subject_info'] = subject_info
   
    outlier_files = []
    for outlier_file in modelspec_args.get('outlier_files',[]):
        outlier_files.append(os.path.join(env['data_path'], outlier_file))
    if len(outlier_files) > 0:
        modelspec_args['outlier_files'] = outlier_files
            
    for parameter in ['functional_runs', 'event_files', 'regressors', 'regressor_names', 'pmod']:
        if modelspec_args.get(parameter):
            del modelspec_args[parameter]        
    
    return modelspec_args

def build_pipeline(job, inSingularity=False):

    job_name = 'task-{}_model-{}'.format(job['Info']['task'], job['Info']['model'])
    
    env = job['Environment']
    
    if inSingularity:
        env['data_path'] = '/data'
        env['output_path'] = '/output'
        env['working_path'] = '/working'

    env['working_path'] = os.path.join(env['working_path'], job_name)
    env['output_path'] = os.path.join(env['output_path'], job_name)
    
    for parameter in ['functional_runs', 'event_files', 'regressors', 'outlier_files']:
        if job['SpecifySPMModel'].get(parameter):
            job['SpecifySPMModel'][parameter] = ensure_list(job['SpecifySPMModel'][parameter])
        
    # IsotropicSmooth - Smooth functional files
    if job.get("IsotropicSmooth",{}).get("fwhm",0) != 0:
        isotropicsmooth_args = job['IsotropicSmooth']
        preproc = pe.MapNode(interface=fsl.IsotropicSmooth(**isotropicsmooth_args), iterfield='in_file', name="smooth")
        preproc.inputs.in_file = [os.path.join(env['data_path'], func_file) for func_file in job['SpecifySPMModel']['functional_runs']]
        
    elif job['SpecifySPMModel']['functional_runs'][0].endswith('.gz'):
        preproc = pe.MapNode(interface=Gunzip(), iterfield='in_file', name="gunzip")
        preproc.inputs.in_file = [os.path.join(env['data_path'], func_file) for func_file in job['SpecifySPMModel']['functional_runs']]
        
    else:
        preproc = pe.MapNode(interface=IdentityInterface(fields=['out_file']), iterfield='out_file', name="identity")
        preproc.inputs.out_file = [os.path.join(env['data_path'], func_file) for func_file in job['SpecifySPMModel']['functional_runs']]

    # SpecifySPMModel - Generates SPM-specific model
    modelspec_args = make_model(job['SpecifySPMModel'], **env)
    modelspec = pe.Node(model.SpecifySPMModel(**modelspec_args), name="modelspec")

    # Level1Design - Generates an SPM design matrix
    level1design_args = job['Level1Design']
    if level1design_args.get('mask_image'):
        level1design_args['mask_image'] = os.path.join(env['data_path'], level1design_args['mask_image'])
    level1design = pe.Node(spm.Level1Design(**level1design_args), name="level1design")

    # EstimateModel - estimate the parameters of the model
    level1estimate_args = job['EstimateModel']
    level1estimate = pe.Node(spm.EstimateModel(**level1estimate_args), name="level1estimate")

    l1analysis_io = [(modelspec, level1design, [('session_info', 'session_info')]),
                     (level1design, level1estimate, [('spm_mat_file', 'spm_mat_file')])]

    datasink_io = [('level1estimate.beta_images','@betas'),
                   ('level1estimate.mask_image','@mask'),
                   ('level1estimate.residual_images','@residuals')]

    # EstimateContrast - estimates contrasts
    if len(job.get('EstimateContrast',{}).get('contrasts',[])) > 0:
        conestimate_args = job['EstimateContrast']
        conestimate_args['contrasts'] = [tuple(contrast) for contrast in conestimate_args['contrasts']]
        conestimate = pe.Node(spm.EstimateContrast(**conestimate_args), name="conestimate")
        
        l1analysis_io += [(level1estimate, conestimate, [('spm_mat_file', 'spm_mat_file'),
                                                         ('beta_images', 'beta_images'),
                                                         ('residual_image', 'residual_image')])]
                                                            
        datasink_io += [('conestimate.spm_mat_file','@spm'),
                        ('conestimate.con_images','@con'),
                        ('conestimate.spmT_images','@spmT'),
                        ('conestimate.spmF_images','@spmF')]
    else:
        datasink_io += [('level1estimate.spm_mat_file','@spm')]                                                 

    # Initiation of the 1st-level analysis workflow
    l1analysis = pe.Workflow(base_dir=env['working_path'], name='l1analysis')
    l1analysis.connect(l1analysis_io)

    # Datasink - creates output folder for important outputs
    datasink = pe.Node(nio.DataSink(base_directory=env['output_path'], container="sub-"+job['Info']["sub"]), name="datasink")
    
    # Putting everything together into a single pipeline
    pipeline = pe.Workflow(base_dir=env['working_path'], name="sub-"+job['Info']["sub"])
    
    pipeline.connect([(preproc, l1analysis, [('out_file', 'modelspec.functional_runs')]),
                      (l1analysis, datasink, datasink_io)])
    
    return pipeline

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="First level analysis")
    parser.add_argument("job", help="path to JSON file with job description")
    args = parser.parse_args()       

    if os.environ.get('SINGULARITY_CONTAINER'):
        inSingularity = True
    else:
        inSingularity = False

    with open(args.job, 'r') as f:
        job = json.load(f)

    # Set the way FSL and Matlab should be called
    mlab.MatlabCommand.set_default_matlab_cmd("matlab -nodesktop -nosplash")        
    if not inSingularity:
        
        spm_path = job['Environment']['spm_path']
        mlab.MatlabCommand.set_default_paths(spm_path)

        fsl_path = job['Environment']['fsl_path']
        if os.environ.get('FSLDIR') is None:            
            paths = os.environ.get('PATH',"").split(os.pathsep)
            paths.append(os.path.join(fsl_path, "bin"))
            os.environ['PATH'] = os.pathsep.join(paths)
            
    os.environ['FSLOUTPUTTYPE'] = 'NIFTI'        

    pipeline = build_pipeline(job, inSingularity)
    pipeline.run()