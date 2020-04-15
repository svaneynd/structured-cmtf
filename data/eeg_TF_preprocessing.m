function [ Z , params ] = eeg_TF_preprocessing( basedir , patient , fs , TR , nscans )
% EEG_TF_PREPROCESSING converts continuous (potentially filtered) EEG data
% into a multi-channel time-frequency spectrogram, using the Thomson
% multitaper method or the short-time Fourier transform.
%
% INPUTS
% - basedir = base directory of code execution
% - patient = index of the patient whose data need to be processed
% - ied_times = (#sessions x 1) cell array, containing the time points (in s)
% at which IEDs occurred in each session
% - fs = EEG sampling frequency in Hz
% - TR = repetition time of fMRI in s
% - nscans = number of scans (TRs) in each acquisition run ('session')
%
% OUTPUTS
% - Z = (#sessions x 1) cell array, containing a (time bins x frequencies x channels) 
% EEG spectrogram
%
% Author: Simon Van Eyndhoven (simon.vaneyndhoven@gmail.com)

%% Parameters
params = struct;

% -- files and directories
params.data.targetdir = basedir;

params.data.eegdir = 'data\eeg';
mkdir(params.data.targetdir,params.data.eegdir)

params.data.specfile = 'p%02d_data.mat'; % filename to save the spectrogram data
params.data.mwffile = 'p%02d_eegmwf.mat'; % file of filtered EEG data

% -- time-frequency transform
params.timefreq.foi = 0.5:0.5:40; % frequencies of interest
   

%% Load MWF filtered EEG data of the current patient
load( fullfile(params.data.targetdir,params.data.eegdir,...
        sprintf(params.data.mwffile,patient)) );  

% -- extract timing info
params.timefreq.TR = TR;
params.timefreq.fs = fs;
nsess = numel(eegmwf); % number of sessions
L = cellfun(@(x)size(x,2),eegmwf); % lenghts in samples of all EEG sessions

params.timefreq.L = []; % empty the timing info for the patient     

nchans = size(eegmwf{1},1);
    
%% Compute STFT for every of the sessions
Z = cell(nsess,1);
params.timefreq.L = cell(nsess,1);
for sess = 1:nsess   
    fprintf('\n/ session %d of %d/\n',sess,nsess)

    % -- determine the part of EEG that belongs to this session        
    ntbins = round( L(sess)/(params.timefreq.fs*params.timefreq.TR) ); % number of repetition times

    % -- check that the number of bins corresponds approximately to the
    % number of fMRI scans

    params.timefreq.L{sess}.fmri = nscans(sess);
    params.timefreq.L{sess}.eeg = ntbins;

    nt = round( params.timefreq.fs*params.timefreq.TR ); % number of EEG samples per scan

    % -- compute spectra 
    Z{sess}.pmtm = zeros( ntbins , length(params.timefreq.foi) , nchans );
    Z{sess}.stft = zeros( ntbins , length(params.timefreq.foi) , nchans );

    % 1. compute spectrum using Thomson multi-taper method
    fprintf('* computing spectrum using Thomson multi-taper method\n')
    for timebin = 1:ntbins
        t0 = (timebin-1)*params.timefreq.TR; % start second of time bin
        trange = find( 1/fs*(1:size(eegmwf{sess},2)) >= t0 , nt);
        Z{sess}.pmtm(timebin,:,:) = pmtm(eegmwf{sess}(:,trange)',[],params.timefreq.foi,params.timefreq.fs);
    end    

    % 2. compute spectrum using STFT
    fprintf('* computing spectrum using STFT\n')
    for chan = 1:nchans
        % compute the FFT using sliding Hamming windows on the
        % zero-padded signal
        [~,~,stfttimes,spectemp] = spectrogram( [ zeros(nt,1) ; eegmwf{sess}(chan,:)' ; zeros(nt,1) ] ,2*nt,nt,params.timefreq.foi,params.timefreq.fs);
        Z{sess}.stft(:,:,chan) = spectemp(:,1:ntbins)';
    end

end    

%% Save preprocessed data            
save( fullfile(params.data.targetdir,params.data.eegdir,...
        sprintf(params.data.specfile,patient) ) ,...
        'Z','params' )    

    
end