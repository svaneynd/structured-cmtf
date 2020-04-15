function [ eegmwf , eegmwfenv , params ] = eeg_mwf_processing( basedir , patient , eeg_sessiondata , ied_sessiontimes , fs , TR , nscans )
% EEG_MWF_PROCESSING filters continuous EEG data using a multi-channel
% Wiener Filter (MWF), to enhance the signal-to-noise ratio of interictal 
% epileptic discharges (IEDs) or spikes w.r.t. the background EEG waveforms.
%
% INPUTS 
% - basedir = base directory of code execution
% - eeg_sessiondata = (#sessions x 1) cell array, containing a (channels x
% time points) matrix of EEG data for each fMRI acquisition run
% - ied_times = (#sessions x 1) cell array, containing the time points (in s)
% at which IEDs occurred in each session
% - fs = EEG sampling frequency in Hz
% - TR = repetition time of fMRI in s
% - nscans = number of scans (TRs) in each acquisition run ('session')
%
% OUTPUTS
% - eegmwf = (#sessions x 1) cell array, containing a (channels x time points) 
% matrix of filtered EEG data for each fMRI acquisition run
% - eegmwfenv = (#sessions x 1) cell array, containing an average EEG envelope 
% of the filtered EEG data for each fMRI acquisition run, which can serve
% as a reference time course
%
% Author: Simon Van Eyndhoven (simon.vaneyndhoven@gmail.com)

%% Add the toolbox that implements multi-channel Wiener filtering
addpath(genpath(fullfile(basedir,'utils\MWF-artifact-removal')))

%% Parameters
params = struct;

% -- envelope computation
params.env.wsh = TR;
params.env.wlen = TR;
params.env.pow = 2;
params.env.prctile = 50;

% -- directories
params.data.targetdir = basedir;

params.data.eegdir = 'data\eeg';
mkdir(params.data.targetdir,params.data.eegdir)

params.data.eegmwffile = 'p%02d_eegmwf.mat';
params.data.eegmwfenvfile = 'p%02d_eegmwfenv.mat';

% -- start and end of the interval around the spike. this is best taken
% quite large, since the spikes are very sparse compared to the full data,
% so a bit of noise leaking into the signal covariance doesn't matter that
% much, since there are plenty of noise frames to estimate the noise
% covariance
params.interval.pre = 0.3; % in seconds
params.interval.post = 1.0; % in seconds

% -- high-pass filter
params.fs = fs;
params.filter.order = fs; % 2*fs
params.filter.fc = 0.5; % 0.75
hpf = fir1(params.filter.order,params.filter.fc/(fs/2),'high');
delay = grpdelay(hpf,1,fs);
    assert(mean(delay)/std(delay)>1e6) % group delay should be constant over frequencies for linear phase filter
delay = mean(delay);
params.filter.hpf = hpf;
params.filter.delay = delay;

% -- initialize the MWF parameters
params.mwf = mwf_params('rank',     'poseig', ...
                        'delay',    4       );
                  
%% Initialize the computation

% -- make sure that the sessions can be concatenated horizontally
if size(eeg_sessiondata,2) == 1
    eeg_sessiondata = eeg_sessiondata';
end

nsess = numel(eeg_sessiondata); % number of sessions
L = cellfun(@(x)size(x,2),eeg_sessiondata); % lenghts in samples of all sessions

%% Apply high-pass filtering
for sess = 1 : nsess
    % -- filter the data of the current session
    temp = fftfilt(params.filter.hpf,eeg_sessiondata{sess}');

    % -- apply a correction for end effects
    temp(1:params.filter.order+1,:) = 0; 

    % -- perform additional mean subtraction (not strictly necessary)
    temp = bsxfun(@minus,temp,mean(temp));

    % -- store the result back in the dataset
    eeg_sessiondata{sess} = temp';    

end
clear('temp')
    

%% Prepare training labels for the multi-channel Wiener filter
% -- initialize empty variables
ied_mask = cell(1,nsess); % binary mask that indicates which samples belong to IED periods (1) or to background EEG (0)
ied_epochs = []; % array in which the data segments around spikes will be stored

% -- find the times of the spikes (intra session)
ied_times = cellfun(@(x)round(x*fs),ied_sessiontimes,'UniformOutput',false);

% -- convert these times to 'absolute' time (i.e., w.r.t. to the
% beginning of the first session)
offset = 0;
for sess = 1:numel(ied_times)
    ied_times{sess} = ied_times{sess} + offset + params.filter.delay;
    offset = offset + size(eeg_sessiondata{sess},2);
end
ied_times = cell2mat(cellfun(@(x)x',ied_times,'UniformOutput',false));
params.data.ied_times = ied_times;

% -- segment the spikes and make spike mask
for sess = 1:numel(eeg_sessiondata)
    ied_mask{sess} = zeros(size(eeg_sessiondata{sess},2),1);
    for ied = 1:numel(ied_sessiontimes)
        % -- find the time of the current spike
        ied_time = ied_sessiontimes{sess}(ied);
        ied_time = params.filter.delay + round([ied_time-params.interval.pre,ied_time,ied_time+params.interval.post]*fs);
        % -- make sure the interictal interval lies within the recording
        ied_time(1) = max(ied_time(1),1);
        ied_time(3) = min(ied_time(3),size(eeg_sessiondata{sess},2));
        offset = round(params.interval.pre*fs) - (ied_time(2) - ied_time(1));
        if offset ~= 0, fprintf('offset = %d\n',offset), end
        % -- indicate the interval around the spike as interictal
        ied_mask{sess}(ied_time(1):ied_time(3)) = 1;
        % -- add the segmented spike
        ied_epoch = zeros(size(eeg_sessiondata{sess},1),round(params.interval.pre*fs)+round(params.interval.post*fs)+1);
        ied_epoch(:,offset + (1:ied_time(3)-ied_time(1)+1)) = eeg_sessiondata{sess}(:,ied_time(1):ied_time(3));
        ied_epochs = cat(3,ied_epochs,ied_epoch);   
    end
end
params.data.ied_mask = cat(1,ied_mask{:});
params.data.epochs = ied_epochs; % store also the epochs for potential later use

clear('offset','spikes','spike_epoch','ied','ied_time','sess')
    
%% Perform multi-channel Wiener filtering
% -- concatenate data from all sessions
eeg = cell2mat(eeg_sessiondata);

% -- normalize the variance of all channels
eeg = zscore(eeg,[],2); 

% -- compute the filter
W = mwf_compute(eeg,params.data.ied_mask,params.mwf);

% -- apply the filter to the data
[~,eegmwf] = mwf_apply(eeg,W);

% -- organize the EEG data back into the original cell structure
eegmwf = mat2cell(eegmwf,size(eegmwf,1),L);
    
%% Compute an average envelope of the filtered EEG as a reference IED time course
% -- initialize an empty variable
eegmwfenv = cell(nsess,1);

% -- calculate the envelope for every of the sessions
for sess = 1:nsess    

    % -- compute the envelope for each session
    [~,eegmwfenv{sess}] = multichan_envelope(eegmwf{sess},2,params.env.wsh,params.env.wlen,params.env.pow,fs);
    eegmwfenv{sess} = eegmwfenv{sess}(:);

    % -- correct for the noisefloor
    noisefloor = prctile(eegmwfenv{sess}(:,1),params.env.prctile);
    params.env.nf(sess) = noisefloor;
    eegmwfenv{sess}(:,1) = eegmwfenv{sess}(:,1) - params.env.nf(sess);

    % -- check the length
    if size(eegmwfenv{sess},1) > nscans(sess)
        result = size(eegmwfenv{sess},1) - nscans(sess);
        eegmwfenv{sess} = eegmwfenv{sess}(1:end-result,:);
    elseif size(eegmwfenv{sess},1) < nscans(sess)
        result = nscans(sess) - size(eegmwfenv{sess},1);
        eegmwfenv{sess} = [ eegmwfenv{sess} ; repmat(median(eegmwfenv{sess}(end-10+1:end,:),1),result,1) ]; % pad with the median of the last 10 samples
    end 
end

%% Save all data
% -- filtered EEG data 
mwffile = fullfile( params.data.targetdir , params.data.eegdir , ...
                    sprintf(params.data.eegmwffile,patient) );
save( mwffile , 'eegmwf' )

% -- envelope of the filterd EEG data
mwfenvfile = fullfile( params.data.targetdir , params.data.eegdir , ...
                    sprintf(params.data.eegmwfenvfile,patient) );
save( mwfenvfile , 'eegmwfenv' )
    
end