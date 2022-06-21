
%========================================================================%
% GPPI_WRAPPER
%
%   This script allows you to run Donald McLaren's gPPI toolbox to build 
%   and estimate psychophysiological interaction (PPI) models for multiple
%   subjects and for multiple seed volumes-of-interest (VOI). Note that
%   before using this, you must first create image files corresponding to
%   your seed regions. 

%  SECTIONS REQUIRING INPUT:
%  * Paths to input and output directories and models
%  * ROIs [earlier version of this code required all ROIs to be in the same directory; this does not, though it does require manual listing of ROI paths and names]
%  * Subjects
%  * Defining conditions in the model - - make sure that you include ALL task conditions / regressors in this list. Code is provided to read out regressor names from first level SPM.mat files.
%  * Edit contrasts in LOOP section (do not need to edit other settings in LOOP)
%
%   History
%   Feb 8, 2020 - NC adapting for DARPA w/ RP and EB; using most recent version of gPPI toolbox (PPPIv13) and updated notes from McClaren and Mumford
%       Moving away from Bob Spunt's wrapper (as in data00/projects/socnets_darpa2/scripts/PPPIv13_test_nc.ipynb); now just using SPM toolbox
%=========================================================================%

clear all; home;

gPPIpath='/data00/tools/PPPIv13.1';
addpath(genpath(gPPIpath)); % path for gPPI toolbox
addpath(genpath('/data00/tools/spm12')) 

rmpath('/data00/tools/batch8');
rmpath('/data00/tools/spm8');


%---------PATHS TO INPUT AND OUTPUT DIRECTORIES AND MODELS------%
studypath = '/data00/projects/onr_fmri/derivatives/nipype'; % path for study directory
level1name = '/task-art_model-behav/'; % path of level1 analysis (relative to subject folder)
ppi_folder_affix = 'PPItest0315';


% -----------ROIs-----------
% Edit "regionfile" and "regions" lists

regionfile={'/data00/projects/socnets_darpa1/data/ROIs/PPI_ROIs/reward/reward_neurosynth.nii'};
regions={'reward_neurosynth'};

nmasks = length(regions);
domasks = 1:nmasks;


%-----------Subjects ------------------
% Edit subjectpattern and subTAG

subjectpattern = 'sub-OM12*'; % pattern for finding subject folders (use wildcards)
subTAG = 'all'; 

%%subTAG = 'all'; %  do which subjects? 'all' to do all, position indices to do subset, e.g., [1 3 7]
%subTAG = [ 1 2 3 4 5 6 8 9 10 11 12 13 14 15 16 17 18 20 21 22 23 24 26 27 28 29 30 32 33 34 37 38 39 40 ]; %exclude 7 & 19/ %25,31,35,36 do not have run 3s

cd([studypath, level1name]);
fprintf('\nSUBJECT LIST:\n');
d=dir(subjectpattern);
for i=1:length(d)
    subnam{i}=d(i).name;
    subdir{i} = [studypath, level1name,  subnam{i}];
    fprintf('\tAdding %s to subject list\n',subnam{i})
end
nsubs = length(subnam);
if strcmp(subTAG,'all')
    dosubs = 1:nsubs;
else
    dosubs = subTAG;
end

%-------- Defining conditions in the model -------%
% Visually confirm that the correct conditions are specified
% -- This code pulls the regressor list from the first run of the first subject
% -- Ensure that there are no conditions missing
% -- Use this list to create the conditions list in the next cell

load([[studypath, level1name, subnam{1}, '/','SPM.mat']]) 
SPM.Sess(1).U.name



% Set condition list here - use ALL regressors 
conditions = { 'final_stage' 'group_HIGHER_change' 'group_HIGHER_nochange' 'group_LOWER_change' 'group_LOWER_nochange' 'group_SAME_change' 'group_SAME_nochange' 'remind'}; % open SPM.mat file, Sess to see the regressor names. conditions to compute PPIs for (must match names used in level 1 analysis)


% LOOP OVER ALL SUBJECTS AND REGIONS TO RUN PPPI

for ii=1:numel(dosubs)
    for jj=1:numel(regions)
        cd(gPPIpath)
        clear P 
        load('/data00/tools/PPPIv13.1/PPPIv13/parameters.mat'); % load parameters template
%----------EDIT CONTRASTS------------%
        P.Contrasts(1).name = 'HIGHER_change_vs_nochange';
        P.Contrasts(1).left = {'group_HIGHER_change'};
        P.Contrasts(1).right = {'group_HIGHER_nochange'}; %right is the one that gets subtracted - if you want it to be rest put none

        %P.Contrasts(2) = P.Contrasts(1);
        %P.Contrasts(2).name = 'bCHANGEcLOWERvsbNOCHANGE';
        %P.Contrasts(2).left = {'bCHANGEcLOWER'};
        %P.Contrasts(2).right = {'bNOCHANGEgHIGHER' 'bNOCHANGEgLOWER'}; %right is the one that gets subtracted - if you want it to be rest put none


%------------------------------------%
    

        P.subject=subnam{ii};
        P.VOI.VOI = [regionfile{jj}];%filepath and name of the VOI
        P.VOI.masks = {[studypath, level1name, subnam{ii}, '/', 'mask.nii']}; %images to specify the subject specific ROI; % name of images to use to constrain definition of the seed region; default here is the mask.img file within the first level model
        P.VOI.thresh = 0.5; %threshold (statistic value) for the VOI
        P.VOI.min = 10; %minimum VOI size required
        P.VOImin = 10; %minimum VOI size required
        P.Region = [regions{jj} '_' ppi_folder_affix]; %string containing the basename of the output file(s)
        P.directory=[studypath level1name subnam{ii}]; %path to the first level SPM directory
        P.Estimate = 1; %estimate design (2 means already estimated, 0 means do not estimate)
        P.contrast = 1; %F contrast to adjust for (corresponds to number of the ess image in level 1 folder)
        P.extract='eig'; % extract timecourse as first eigenvariate ('eig') or mean  ('mean')
        P.Tasks = ['0' conditions]; %In the generalized context-dependent PPI, you need specify the tasks to include in the analyses, but put a ‘0’ or ‘1’ in front of them to specify if they must exist in all sessions
        P.Weights=[]; %only relevant for traditional PPI
        P.maskdir = []; %location to store VOI file if input VOI was a mat file
        P.equalroi = 0; %specifies that the ROIs must be the same size in all subjects (1 for true, 0 for false)
        P.FLmask = 1; %specifies that ROI should be restricted in each subject by the mask image (1 for true, 0 for false)
        P.analysis = 'psy'; %specifies psychophysiological interaction ('psy'), physiophysiological interaction ('phys'), or psychophysiophysiological interaction ('psyphy')
        P.method = 'cond'; %specifies traditional spm ppi ('trad') or generalized context-dependent ppi ('cond')
        P.CompContrasts = 1; %estimate contrasts
        P.Weighted = 0; %default is to not weight tasks by number of trials (0); to change this, specify which tasks should be weighted by trials
   
        spmfile=[studypath, level1name, subnam{ii}, '/','SPM.mat'];
        load(spmfile)
        
        PPPI(P);

    end
end
        

spmfile=[studypath, level1name, subnam{ii}, '/','SPM.mat'];
        load(spmfile)

for ii=1:numel(SPM.xY.VY)
    a{ii}=SPM.xY.VY(ii).fname;
end
a=unique(a);
filesgz={}; filesbz={};
for ii=1:numel(a)
    if ~exist(a{ii},'file')
       disp([a{ii} ' DOES NOT EXIST.'])
    else
       disp([a{ii} ' EXISTS']) 
    end
end


