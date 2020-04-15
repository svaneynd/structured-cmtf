function [ sol , output , results , params  ] = init_eeg_TF_cpd( basedir , patient , ranks )
% INIT_EEG_TF_CPD fits a canonical polyadic decomposition (CPD),
% otherwise known as parallel factor analysis (PARAFAC) model, to a
% third-order EEG spectrogram tensor (time bins x frequencies x channels).
% Since the the algorithm to compute the CPD is non-convex, there's no
% guarantee that an optimal decomposition can be found. Therefore, the
% algorithm is run multiple times, to increase the odds of finding a
% 'good' solution. In each run, a different initialization is used for the
% CPD factors (either random or non-random).
% 
% INPUTS
% - basedir = base directory of code execution
% - patient = index of the patient whose data need to be processed
% - ranks = (range of) integer(s) indicating how many components the CPD
% model should contain
%
% OUTPUTS
% - sol = cell array, containing CPD factor sets ('solutions') after 
% repeated runs of of the algorithm for independent initialization
% - output = cell array, containing information on approximation error and
% convergence for each run of the algorithm
% - results = summary information on the rescaling that was performed on
% the EEG data tensor
%
% Author: Simon Van Eyndhoven (simon.vaneyndhoven@gmail.com)

%% Parameters
recompute = true; % if this decomposition has already been performed, you can choose
                    % to skip it
                    
params = struct;

% -- directories & files
params.basedir = basedir;

params.data.dir = 'data\eeg';
params.data.datafile = 'p%02d_data.mat';

params.results.dir = 'init';
    mkdir(params.basedir,params.results.dir)        
params.results.subdir = 'p%02d';        
params.results.solfile = 'normcpd_%s_R%d_Tnorm.mat';    
     
% -- spectrogram type
params.data.spectype = 'pmtm'; % 'pmtm' / 'stft'

% -- decomposition    
params.decomposition.Initialization = @cpd_rnd; % initialization strategy
params.decomposition.InitializationOptions.Slices = 'random';
params.decomposition.InitializationOptions.OptimalScaling = true;
params.decomposition.nrep = 50; % total number of repetitions of the (non-convex) optimization algorithm to fit a CPD
params.decomposition.crossval.K = 10; % number of repetitions for which a non-random initialization strategy is used
params.decomposition.TolFun = 1e-8; % convergence criterion for cost function
params.decomposition.TolX = 1e-10; % convergence criterion for step
params.decomposition.Display = 500; % print progress every ... iterations
params.decomposition.MaxIter = 2000; % 
    params.decomposition.CGMaxIter = 400;
    
%% Initialize parallel computation
if isempty(gcp('nocreate'))
    parpool
end        
    
%% Load the EEG data tensor            
% -- load data for the current patient 
eegfile = fullfile( params.basedir , params.data.dir , ...
                    sprintf(params.data.datafile,patient) );
eeg = load(eegfile);

% -- concatenate data across sessions
X = [];
for sess = 1:length(eeg.Z)
    sessdata = eeg.Z{sess}.(params.data.spectype);
    % adjust the length if needed: fMRI data dictates which length
    % should be taken
    excessL = eeg.params.timefreq.L{sess}.eeg - eeg.params.timefreq.L{sess}.fmri;
    if excessL > 0 % EEG data is longer than fMRI data : apply truncation
        sessdata = sessdata( 1:eeg.params.timefreq.L{sess}.eeg-excessL , : , : );
    elseif excessL < 0 % EEG data is shorter than fMRI data : apply zero-padding
        sessdata = cat( 1 , sessdata , zeros(-excessL,size(sessdata,2),size(sessdata,3)) );
    else % EEG and fMRI data are already of equal length for this session, which is good
    end
    X = cat( 1 , X , sessdata );
end    

% -- create a results directory
mkdir( fullfile(params.basedir,params.results.dir) , sprintf(params.results.subdir,patient) )

results = struct;

%% Preprocess the tensor data
% -- normalize the data over spectral and spatial mode
results.mu = mean(X,1);
X = bsxfun( @minus , X , results.mu );
[X,results.fnorms,results.scaling] = normalize_array( X , [2,3] );

% -- adjust for the norm of the dataset    
results.datanorm = frob(X);    
X = X / results.datanorm; % ensure that the final tensor has unit Frobenius norm

% -- create a subdivision in cross-validation folds for the non-random initialization
try
    cv_idx = crossvalind('Kfold',size(X,1),params.decomposition.crossval.K);    
catch
    cv_idx = randi(params.decomposition.crossval.K,size(X,1),1);
end

%% Compute the CPD for each desired number of components
for r = ranks
    % -- create the file where the results should be saved
    solfile = fullfile( params.basedir , params.results.dir , ...
                sprintf(params.results.subdir,patient) , ...
                sprintf(params.results.solfile,params.data.spectype,r) );

    sol_exists = exist(solfile,'file');

    if recompute | ~sol_exists

    % -- create empty results variables
    sol = cell(1,params.decomposition.nrep);
    output = cell(1,params.decomposition.nrep);

    % -- first part: non-random initialization using sliced tensor and GEVD
    parfor rep = 1:params.decomposition.crossval.K
        % -- initialize using GEVD for a part of the slices
        U0 = cpd_gevd(X(cv_idx==rep,:,:),r,params.decomposition.InitializationOptions);

        % -- compute the factor in the sliced mode by regression
        U0{1} = kr(U0{3},U0{2})\tens2mat(X,[2,3],1);
        U0{1} = U0{1}';

        % -- decompose the full tensor
        fprintf('\n%s\nRepetition nr. %d\n%s\n',repmat('-',1,30),rep,repmat('-',1,30))
        [sol{rep},output{rep}] = cpd(X,U0,params.decomposition);
    end

    % -- second part: using random initialization
    parfor rep = params.decomposition.crossval.K + 1 : params.decomposition.nrep

        % -- decompose the full tensor
        fprintf('\n%s\nRepetition nr. %d\n%s\n',repmat('-',1,30),rep,repmat('-',1,30))
        [sol{rep},output{rep}] = cpd(X,r,params.decomposition);

        assert(strcmp(output{rep}.Initialization.Name,'cpd_rnd')) % assert that random initialization was used
    end

    % -- save results
    save(   solfile , ...
            'sol','output','params','results') 
    end

    clc;

end
    
end