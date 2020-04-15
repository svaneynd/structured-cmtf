% MAIN_PIPELINE runs the end-to-end pipeline for joint EEG-fMRI data
% analysis, which is described in 
%
% Van Eyndhoven, S., Dupont, P., Tousseyn, S., Vervliet, N., Van Paesschen, 
% W., Van Huffel, S. & Hunyadi, B., "Augmenting interictal mapping with 
% neurovascular coupling biomarkers by structured factorization of 
% epileptic EEG and fMRI data", submitted to NeuroImage
%
% Author: Simon Van Eyndhoven (simon.vaneyndhoven@gmail.com)
% Date: 15/4/2020


%% Initialization
clc
clearvars
close all;

% -- define the root/base directory
basedir = pwd; % base directory to which all subdirectories are referenced
    % (potentially: replace by your own directory after extracting all code)

% -- add some toolbox functionality to path
addpath(genpath( basedir ))

% -- make sure that SPM is on the path
assert(exist('spm.m')~=0,'You need have the SPM toolbox installed!')

%% Parameters 
ranks = 2 : 4 ; % select for which range of # components the sCMTF should be computed

%% Load data
% <!> Below, we generate (utterly meaningless!) EEG and fMRI data, to
% illustrate the format of the input data which is expected by the
% algorithm. With these dummy data, it is possible to run the algorithms
% below. Replace this block of code by a loading module for your own data.

patient = 1; % index of patient

nchans = 20; % number of EEG channels
fs = 250; % EEG sampling frequency (Hz)
TR = 2.5; % MR repetition time (s)
L = [500*fs,500*fs,800*fs]; % length of each acquisition run (in EEG samples)
nscans = L/(fs*TR); % length of each acquisition run (in MR scans)
eeg_sessiondata = cell(3,1); % cell array of EEG data (channels x time points) for each acquisition run ('session')
fmri_sessiondata = cell(3,1); % cell array of fMRI data (SPM volume info, not yet the full volumes) for each acquisition run ('session')
ied_sessiontimes = cell(3,1); % cell array of IED occurrences (time of IED peak) for each acquisition run ('session')
R = cell(3,1); % cell array of nuisance regressors for the fMRI data of each acquisition run ('session')

v = spm_vol(fullfile(basedir,'\utils\atlas\white_2mm.nii'));

% generate dummy 'data'  for a dataset of 3 acquisition runs
for i = 1 : 3
    eeg_sessiondata{i} = randn(nchans,L(i));
    fmri_sessiondata{i} = struct;
    fmri_sessiondata{i}.vol = repmat(v,nscans(i),1); % store the SPM volume info
    ied_sessiontimes{i} = randi(floor(L(i)/fs),20,1); % mark 20 random time points as IED peaks
    R{i} = randn(floor(nscans(i)),5); % generate 5 random time courses that serve as nuisance variables
end

%% Pipeline of computation

% -- signal enhancement of IEDs in the EEG data (2.3)
eeg_mwf_processing( basedir , patient , eeg_sessiondata , ied_sessiontimes , fs , TR , nscans );

% -- EEG data preprocessing (2.4.1)
eeg_TF_preprocessing( basedir , patient , fs , TR , nscans )

% -- BOLD data preprocessing (2.4.2)
fmri_preprocessing( basedir , patient , fmri_sessiondata , R , TR )

% -- EEG-only factorization for initialization (A.2)
init_eeg_TF_cpd( basedir , patient , ranks )

% -- sCMTF (2.5)
eegfmri_scmtf( basedir , patient , ranks );

% -- Stability analysis (B.2)
for r = ranks
    assess_stability_scmtf( basedir , patient , r , true );
end

fprintf('\n>>> Press any key to continue.\n')
pause;

% -- Statistical inference (2.6)
for r = ranks
    statistical_analysis( basedir , patient , r );
end

% After this step, directories for this patient's data have been created,
% for the sCMTF results of each rank r. These directories (under './inference/pXX')
% host the statistical (de)activation maps for all components, as well as
% the maps of the HRF variability.

