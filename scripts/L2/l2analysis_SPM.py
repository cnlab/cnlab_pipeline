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

def remove_key(keys, obj):
    for key in ensure_list(keys):
        if obj.get(key):
            del obj[key]

def build_pipeline(job, inSingularity=False):

    l2function = {
        "OneSampleTTestDesign": spm.OneSampleTTestDesign,
        "PairedTTestDesign": spm.PairedTTestDesign,
        "TwoSampleTTestDesign": spm.TwoSampleTTestDesign,
        "MultipleRegressionDesign": spm.MultipleRegressionDesign
    }

    centering_code = {
        'overall_mean': 1, 
        'factor1_mean': 2, 
        'factor2_mean': 3, 
        'factor3_mean': 4, 
        'no_centering': 5, 
        'user_specified': 6, 
        'as_implied_by_ancova': 7, 
        'gm': 8
    }

    job_name = 'task-{}_model-{}'.format(job['Info']['task'], job['Info']['model'])
    l2_name = 'l2-{}'.format(job['SecondLevel']['name'])
    l2_analysis = job['SecondLevel']['analysis']
    
    env = job['Environment']
    
    if inSingularity:
        env['data_path'] = '/data'
        env['output_path'] = '/output'
        env['working_path'] = '/working'

    env['working_path'] = os.path.join(env['working_path'], job_name)
    env['output_path'] = os.path.join(env['output_path'], job_name)

    l2_args = copy.deepcopy(job['SecondLevel'])

    if l2_args.get('covariate_file'):
        covariate_df = pd.read_csv(os.path.join(env['data_path'], l2_args['covariate_file']))

        if "sub" in covariate_df.columns:
            covariate_df = covariate_df.drop(columns='sub')

        covariates = []
        covariate_centering = l2_args.get('covariate_centering', {})

        for var in covariate_df.columns:
            centering = centering_code[covariate_centering.get(var, 'overall_mean')]
            covariates.append({'name': var, 'vector': covariate_df[var].to_list(), 'centering': centering})

        l2_args['covariates'] = covariates

    if l2_args.get('explicit_mask_file'):
        l2_args['explicit_mask_file'] = os.path.join(env['data_path'], l2_args['mask_image'])

    if l2_analysis in ['OneSampleTTestDesign', 'MultipleRegressionDesign']:
        l2_args['in_files'] = [os.path.join(env['data_path'], in_file) for in_file in l2_args['in_files']]

    elif l2_analysis == 'PairedTTestDesign':
        l2_args['paired_files'] = [[os.path.join(env['data_path'], c1), os.path.join(env['data_path'], c2)] for c1, c2 in l2_args['paired_files']]

    elif l2_analysis == 'TwoSampleTTestDesign':
        l2_args['group1_files'] = [os.path.join(env['data_path'], in_file) for in_file in l2_args['group1_files']]
        l2_args['group2_files'] = [os.path.join(env['data_path'], in_file) for in_file in l2_args['group2_files']]

    remove_key(['name', 'analysis', 'covariate_file', 'covariate_centering'], l2_args)

    # Level1Design - Generates an SPM design matrix
    level2design = pe.Node(l2function[l2_analysis](**l2_args), name="level2design")

    # EstimateModel - estimate the parameters of the model
    level2estimate = pe.Node(spm.EstimateModel(estimation_method={'Classical': 1}), name="level2estimate")

    # EstimateContrast - estimates contrasts
    conestimate_args = job['EstimateContrast']
    conestimate_args['group_contrast'] = True
    conestimate_args['contrasts'] = [tuple(contrast) for contrast in conestimate_args['contrasts']]
    conestimate = pe.Node(spm.EstimateContrast(**conestimate_args), name="conestimate")

    l2analysis_io = [(level2design, level2estimate, [('spm_mat_file', 'spm_mat_file')]),
                     (level2estimate, conestimate, [('spm_mat_file', 'spm_mat_file'),
                                                     ('beta_images', 'beta_images'),
                                                     ('residual_image', 'residual_image')])]
                                                        
    datasink_io = [('level2estimate.beta_images','@betas'),
                   ('level2estimate.mask_image','@mask'),
                   ('level2estimate.residual_images','@residuals'),
                   ('conestimate.spm_mat_file','@spm'),
                   ('conestimate.con_images','@con'),
                   ('conestimate.spmT_images','@spmT'),
                   ('conestimate.spmF_images','@spmF')]

    # Initiation of the 1st-level analysis workflow
    l2analysis = pe.Workflow(base_dir=env['working_path'], name='l2analysis')
    l2analysis.connect(l2analysis_io)

    # Datasink - creates output folder for important outputs
    datasink = pe.Node(nio.DataSink(base_directory=env['output_path'], container=l2_name), name="datasink")
    datasink.inputs.substitutions = [('_contrast_name_', '')]
    
    # Putting everything together into a single pipeline
    pipeline = pe.Workflow(base_dir=env['working_path'], name=l2_name)
    
    pipeline.connect([(l2analysis, datasink, datasink_io)])
    
    return pipeline

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Second level analysis")
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