function [ sol , output ] = eegfmri_scmtf( basedir , patient , ranks )
% EEGfMRI_SCMTF computes a structured coupled matrix-tensor factorization
% (sCMTF) of simultaneously recorded EEG and fMRI data, which is
% represented in the form of a (time bins x frequencies x channels) tensor,
% and a (time bins x ROIs) matrix, respectively.
% The data of both modalities share the 'time' mode, which is the coupling
% mechanism in the factorization.
% Since the the algorithm to compute the sCTMF is non-convex, there's no
% guarantee that an optimal decomposition can be found. Therefore, the
% algorithm is run multiple times, to increase the odds of finding a
% 'good' solution. In each run, the algorithm starts from different
% initial factors, which result from the earlier computed factorization of
% the EEG tensor (from which all other initial factors can be derived).
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
%
% Author: Simon Van Eyndhoven (simon.vaneyndhoven@gmail.com)

%% Parameters
recompute = false; % if this decomposition has already been performed, you can choose
                    % to skip it

params = struct;

% -- directories & files
params.basedir = basedir;

params.data.dir = { 'data\eeg' , 'data\fmri' }; % directories where EEG and fMRI data are stored
params.data.datafile = { 'p%02d_data.mat' , 'p%02d_data.mat' }; % filenames for EEG and fMRI data

params.init.dir = 'init'; % directory where the EEG-only factorization is stored (used for initialization)    
params.init.subdir = 'p%02d';        
params.init.solfile = 'normcpd_%s_R%d_Tnorm.mat';

params.results.dir = 'cmtf';
    mkdir(params.basedir,params.results.dir)     
params.results.subdir = params.init.subdir;
params.results.solfile = 'normnorm3H_cmtf_l1%1.1e_%s_%s_R%d.mat';

% -- HRF initialization
params.data.Lnc = 4; % number of time samples that are 'non-causal' (negative lag)

hrfdisp = 1; % HRF dispersion ('width')
peak1 = [-2 4 10]; % time of first HRF peak (in s)
peak2 = peak1 + 10; % time of second (negative) HRF peak (in s)
relamp = 1/4; % relative amplitude of second peak

hrfidx = 2;
params.data.hrfinit_1 = [ peak1(hrfidx)*hrfdisp peak2(hrfidx)*hrfdisp hrfdisp hrfdisp relamp]; % mid-peaking HRF
hrfidx = 3;
params.data.hrfinit_2 = [ peak1(hrfidx)*hrfdisp peak2(hrfidx)*hrfdisp hrfdisp hrfdisp relamp]; % late-peaking HRF
hrfidx = 1;
params.data.hrfinit_3 = [ peak1(hrfidx)*hrfdisp peak2(hrfidx)*hrfdisp hrfdisp hrfdisp relamp]; % early-peaking HRF

params.data.Lhrf = 20; % number of time samples to evaluate the HRF on
params.data.snrhrf = 10; % signal-to-noise ratio of the random perturbations that are applied to the HRF parameters in each run of the algorithm

% -- compression
params.data.compress = false; % compress EEG tensor using MLSVD? (this can speed up computation, at the cost of accuracy)
params.data.varratio = 0.95; % percentage of variance to be retained in every unfolded mode

% -- data feature type
params.data.boldtype = 'avg'; % 'avg' / 'pca'
params.data.spectype = 'pmtm'; % 'stft' (short-time Fourier transform) / 'pmtm' (Thomson multitaper method)

% -- decomposition
L1pen = 1e-3; % L1 penalty on the components' amplitudes
add_fmriuniq = true; % add the 'free' components [[N,P]] that model residual structured variance in the fMRI?

params.decomposition.Q = 3; % number of rank-1 terms that solely model fMRI variance (per acquisition run)

params.decomposition.nrep = 50; % total number of repetitions of the (non-convex) optimization algorithm to fit the sCMTF
params.decomposition.options.Weights = [ 1 1 L1pen L1pen ]; % EEG - fMRI - component amplitudes EEG - component amplitudes fMRI 

do_sdf_nls = false; % NLS algorithm is very computationally demanding --> use a quasi-Newton algorithm (sdf_minf) instead
if do_sdf_nls % Gauss-Newton algorithm
    params.decomposition.options.TolFun = 1e-8;
    params.decomposition.options.TolX = 1e-12;
    params.decomposition.options.Display = 200;
    params.decomposition.options.MaxIter = 500;
    params.decomposition.options.CGMaxIter = 50;
    params.decomposition.algorithm = @(model,options)sdf_nls(model,options); 
else % quasi-Newton algorithm
    params.decomposition.options.TolFun = 1e-8;
    params.decomposition.options.TolX = 1e-12;
    params.decomposition.options.Display = 0;
    params.decomposition.options.MaxIter = 10;         
    params.decomposition.algorithm = @(model,options)sdf_minf(model,options);
end
    
%% Initialize parallel computation
if isempty(gcp('nocreate'))
    parpool
end 
fprintf('\n/// Processing patient %d ///\n\n',patient)

%% Load patient data
fprintf('* load spectrogram EEG data + ROI fMRI data\n')
eeg = load( fullfile( params.basedir , params.data.dir{1} , ...
                sprintf(params.data.datafile{1},patient) ) );
fmri = load( fullfile( params.basedir , params.data.dir{2} , ...
                sprintf(params.data.datafile{2},patient) ) );            

dataparams = struct;
dataparams.eeg = eeg.params;
dataparams.fmri = fmri.params;

% -- concatenate data across sessions    
    assert( length(eeg.Z) == length(fmri.Y) )  % ensure that equally many sessions are present in both datasets

X = [];
Y = [];

for sess = 1:length(eeg.Z)
    % 1. EEG session
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
    % append the data
    X = cat( 1 , X , sessdata );
    sessL(sess) = size(sessdata,1);
    % 2. fMRI session
    sessdata = fmri.Y{sess}.(params.data.boldtype);
    % append the data
    Y = cat( 1 , Y , sessdata );
    clearvars('sessdata')
end    
nsess = length(eeg.Z);

%% Preprocess the data
% -- create a results directory
mkdir( fullfile(params.basedir,params.results.dir) , sprintf(params.results.subdir,patient) )

results = struct;

% -- normalize the data if desired
results.mu_t = mean(X,1);
results.mu_m = mean(Y,1);
X = bsxfun( @minus , X , results.mu_t );
[X,results.fnorms_t,results.scaling_t] = normalize_array( X , [2,3] );
[Y,results.fnorms_m,results.scaling_m] = normalize_array( Y , 2 );


% -- adjust for the norm of the two datasets
results.datanorms = cell(1,2);

results.datanorms{1} = frob(X);
results.datanorms{2} = frob(Y);

X = X / results.datanorms{1};
Y = Y / results.datanorms{2};

% -- compression using MLSVD
if params.data.compress
    [Uml,Gml,sv] = mlsvd(X);
    for mode = 1:3
        cs = cumsum(sv{mode}.^2)./sum(sv{mode}.^2);
        cutoff{mode} = find(cs>params.data.varratio,1);
    end
    for mode = 1:3
        if isempty(cutoff{mode})
            cutoff{mode} = size(Uml{mode},2);
        end
        Uml{mode} = Uml{mode}(:,1:cutoff{mode});
    end
    Gml = Gml(1:cutoff{1},1:cutoff{2},1:cutoff{3});
    froblmlrares(X,Uml,Gml)
else
end

%% Create auxiliary optimization operators
% -- set the number of components that model BOLD variance which is not
% related to EEG
params.decomposition.Q = max(3,2*nsess); % allow at least as many additional rank-1 terms as sessions 

% -- auxiliary parameter to properly adjust the HRF parameters to model
% non-causal basis waveforms
tshift = params.data.Lnc * eeg.params.timefreq.TR

% -- define structure transformations to be used during optimization

    % HRF
    try
        timeax = (0 : params.data.Lhrf-1 ) * eeg.params.timefreq.TR ; % time axis
    catch
        timeax = (0 : size(Y,1)-1 ) * eeg.params.timefreq.TR ; % time axis
    end
    thrf = @(z,task)struct_hrf(z,task,timeax);

    % Khatri-Rao product
    tkr = @(z,task)struct_kr(z,task);

    % transpose
    ttransp = @(z,task)struct_transpose(z,task);

    % selection
    tsel14 = @(z,task)struct_selectN(z,task,[1 4]);
    tsel24 = @(z,task)struct_selectN(z,task,[2 4]);
    tsel34 = @(z,task)struct_selectN(z,task,[3 4]);
    tsel4 = @(z,task)struct_select(z,task,4);

    % non-negativity
    tnonneg = @(z,task)struct_abs(z,task,0.001); 

    % convolution
    toepsize = [ size(Y,1) , size(Y,1) ] ;
    toeppre = zeros( size(Y,1) - 1 - params.data.Lnc , 1 ) ;
    toeppost = zeros( size(Y,1) - length(timeax) + params.data.Lnc , 1 );
    ttoep = @(z,task)struct_toeplitz(z,task,toepsize,toeppre,toeppost);

    % constant
    tconst = @(z,task)struct_const(z,task);

    % diagonal
    tdiag = @(z,task)struct_diag(z,task);

    % L1-norm
    tsq = @(z,task)struct_nonneg(z,task); 
    tsum = @(z,task)struct_sum(z,task,1); 
    tsqrt = @(z,task)struct_sqrt(z,task);

    % normalization
    tnorm = @(z,task)struct_normalize(z,task);

%% Compute the sCMTF for each desired number of components
    
for r = ranks

    fprintf('\n%s\nComputing for rank %d and L1 penalty %3e\n%s\n',repmat('#',1,30),r,L1pen,repmat('#',1,30))

    solfile = fullfile( params.basedir , params.results.dir , ...
                sprintf(params.results.subdir,patient) , ...
                sprintf(params.results.solfile , ...
                    params.decomposition.options.Weights(3) , ...
                    params.data.spectype , ...
                    params.data.boldtype , ...
                    r) );

    sol_exists = exist(solfile,'file');

    if recompute | ~sol_exists

    % -- load the initialized factors
    initfile = fullfile( params.basedir , params.init.dir , ...
                sprintf(params.init.subdir,patient) , ...
                sprintf(params.init.solfile , ...
                    params.data.spectype , ...
                    r) );
    init = load(initfile);

    results.initfile = initfile;

    % -- create empty results variables
    sol = cell(1,params.decomposition.nrep);
    output = cell(1,params.decomposition.nrep);
    model = cell(1,params.decomposition.nrep);

    % -- loop over all pre-computed initializations
    parfor rep = 1:length(init.sol)
        model{rep} = struct;

        % -- initialize the EEG variables
        lambda_x = cell2mat( cellfun(@(x) sqrt(sum(x.^2)) , init.sol{rep}' , 'uniformoutput', false ) );
        model{rep} . variables . s = bsxfun( @rdivide , init.sol{rep}{1} , lambda_x(1,:) ); % variables for the temporal factor matrix
        model{rep} . variables . g = bsxfun( @rdivide , init.sol{rep}{2} , lambda_x(2,:) ); % variables for the spectral factor matrix
        model{rep} . variables . m = bsxfun( @rdivide , init.sol{rep}{3} , lambda_x(3,:) ); % variables for the spatial factor matrix

        model{rep} . variables . lambda_x = prod( lambda_x , 1 ); % amplitudes of the factors

        % -- initialize the HRF variables
        goodhrfinit = false;
        while ~goodhrfinit
            % randomly perturb the HRF parameters
            newpar = params.data.hrfinit_1 .* exp( 10^(-params.data.snrhrf/20) * randn(1,5) );            
            model{rep} . variables . hrfparams_1 = shiftlat(newpar,tshift);
            newpar = params.data.hrfinit_2 .* exp( 10^(-params.data.snrhrf/20) * randn(1,5) );            
            model{rep} . variables . hrfparams_2 = shiftlat(newpar,tshift);
            newpar = params.data.hrfinit_3 .* exp( 10^(-params.data.snrhrf/20) * randn(1,5) );            
            model{rep} . variables . hrfparams_3 = shiftlat(newpar,tshift);
            % verify whether no two basis functions are too correlated
            hrfsinit = [    thrf(model{rep}.variables.hrfparams_1,[]) , ...
                            thrf(model{rep}.variables.hrfparams_2,[]) , ...
                            thrf(model{rep}.variables.hrfparams_3,[]) ];
            hrfcorrs = triu(cossimil(hrfsinit),1);
            if max(abs(hrfcorrs(:)))<0.5
                goodhrfinit = true;
            else                
            end
        end

        % -- initialize the fMRI variables
        HSall = {   ttoep(tnorm(thrf(model{rep}.variables.hrfparams_1,[]),[]),[]) * model{rep}.variables.s , ...
                    ttoep(tnorm(thrf(model{rep}.variables.hrfparams_2,[]),[]),[]) * model{rep}.variables.s , ...
                    ttoep(tnorm(thrf(model{rep}.variables.hrfparams_3,[]),[]),[]) * model{rep}.variables.s };
        if false
            [ Bi , Vi ] = cmtf_multihrf_voxinit( Y , HSall );
        else
            [ Bi , Vi ] = cmtf_multihrf_voxinit_kpe( Y , HSall );
        end
        model{rep} . variables . spatial = { Bi(1,:) , Bi(2,:) , Bi(3,:) , Vi };      

        if add_fmriuniq
        % -- add an uncoupled (i.e., fMRI-specific) part based on the
        % residual
        Yres = Y - cell2mat(HSall)*kr(Bi,Vi);
            assert(frob(Yres)<1) % otherwise the initialization was performed very poorly
        [Umres,Smres,Vmres] = svds(Yres,params.decomposition.Q);

        model{rep} . variables . n = Umres * Smres; % principal component time courses of additional BOLD variance
        model{rep} . variables . p = Vmres; % principal component loadings of additional BOLD variance
        end

        % -- add a dummy variable that is needed for the LMLRA
        % specification
        model{rep} . variables . dummy = 1;

        % -- define the model factors
        model{rep} . factors . S = { 's' , tnorm }; % temporal factor matrix (shared) (Iv x R)
        model{rep} . factors . G = { 'g' , tnorm }; % spectral factor matrix (Ig x R)
        model{rep} . factors . M = { 'm' , tnorm }; % EEG spatial (leadfield) factor matrix (Im x R)
        model{rep} . factors . b1V = { 'spatial' , tsel14 , tkr , ttransp }; % fMRI spatial 'factor matrix' belonging to 1st HRF basis function (Iv x R)
        model{rep} . factors . b2V = { 'spatial' , tsel24 , tkr , ttransp }; % fMRI spatial 'factor matrix' belonging to 2nd HRF basis function (Iv x R)
        model{rep} . factors . b3V = { 'spatial' , tsel34 , tkr , ttransp }; % fMRI spatial 'factor matrix' belonging to 3rd HRF basis function (Iv x R)

        model{rep} . factors . HRFtoep_1 = { 'hrfparams_1' , tnonneg , thrf , tnorm , ttoep }; % Toeplitz matrix holding the 1st HRF basis function
        model{rep} . factors . HRFtoep_2 = { 'hrfparams_2' , tnonneg , thrf , tnorm , ttoep }; % Toeplitz matrix holding the 2nd HRF basis function
        model{rep} . factors . HRFtoep_3 = { 'hrfparams_3' , tnonneg , thrf , tnorm , ttoep }; % Toeplitz matrix holding the 3rd HRF basis function

        if add_fmriuniq
        model{rep} . factors . N = {'n'}; % principal component time courses of additional BOLD variance
        model{rep} . factors . eye_offset = sparse(diag(ones(size(Y,1),1)));
        model{rep} . factors . P = {'p'}; % principal component loadings of additional BOLD variance
        end

        model{rep} . factors . Lambda_x = { 'lambda_x' , tnonneg }; % amplitudes of the EEG components
            % Note: in the derivation in the paper, the amplitude is not
            % explicitly modeled as an additional factor, but both
            % formulations are equivalent.
        model{rep} . factors . Lambda_y = { {'spatial' , tsel14 , tkr , ttransp , tsq , tsum , tsqrt } ; ... % amplitudes of the fMRI components
                                        {'spatial' , tsel24 , tkr , ttransp , tsq , tsum , tsqrt } ; ...
                                        {'spatial' , tsel34 , tkr , ttransp , tsq , tsum , tsqrt } };

        model{rep} . factors . dummy = { 'dummy' , tconst }; % irrelevant mode, since actually a matrix is factorized (but necessary for correct implementation)

        % -- define factorizations
        if params.data.compress
            model{rep} . factorizations . eeg . data = { Uml , Gml };
        else
            model{rep} . factorizations . eeg . data = X;
        end
        model{rep} . factorizations . eeg . cpd = { 'S' , 'G' , 'M' , 'Lambda_x' };            
        model{rep} . factorizations . fmri . data = Y;
        model{rep} . factorizations . fmri . btd = { { 'HRFtoep_1' , 'b1V' , 'dummy' , 'S' } , ...
                                                    { 'HRFtoep_2' , 'b2V' , 'dummy' , 'S' } , ...
                                                    { 'HRFtoep_3' , 'b3V' , 'dummy' , 'S' } };
        if add_fmriuniq
            model{rep} . factorizations . fmri . btd = [ model{rep} . factorizations . fmri . btd , ...
                                                        {{ 'eye_offset' , 'P' , 'dummy' , 'N' }} ];
        end

        model{rep} . factorizations . groupsparseE . regL1 = 'Lambda_x' ;
        model{rep} . factorizations . groupsparseF . regL1 = 'Lambda_y' ;

        % -- solve the CMTF problem
        fprintf('\n%s\nRepetition nr. %d\n%s\n',repmat('-',1,30),rep,repmat('-',1,30))
        [sol{rep},output{rep}] = params.decomposition.algorithm(model{rep},params.decomposition.options);

        % -- compress the highly structured HRF factor and the
        % redundant block term factor
        sol{rep} . factors . HRFtoep_1 = sol{rep} . factors . HRFtoep_1(:,1+params.data.Lnc); % (1+Lnc)th column of the Toeplitz matrix contains the full waveform
        sol{rep} . factors . HRFtoep_2 = sol{rep} . factors . HRFtoep_2(:,1+params.data.Lnc); % (1+Lnc)th column of the Toeplitz matrix contains the full waveform
        sol{rep} . factors . HRFtoep_3 = sol{rep} . factors . HRFtoep_3(:,1+params.data.Lnc); % (1+Lnc)th column of the Toeplitz matrix contains the full waveform
        sol{rep} . factors . eye_offset = [];

        % -- remove redundant elements from the model
        model{rep} . factorizations . eeg . data = [];
        model{rep} . factorizations . fmri . data = [];

    end

    % -- save results
    Hi = struct('hrfparams_1',cell(params.decomposition.nrep,1),'hrfparams_2',cell(params.decomposition.nrep,1),'hrfparams_3',cell(params.decomposition.nrep,1));
    for rep = 1 : params.decomposition.nrep
        Hi(rep) . hrfparams_1 = model{rep} . variables . hrfparams_1;
        Hi(rep) . hrfparams_2 = model{rep} . variables . hrfparams_2;
        Hi(rep) . hrfparams_3 = model{rep} . variables . hrfparams_3;
    end
    results . Hi . hrfparams_1 = cell2mat({ Hi(:) . hrfparams_1 }');
    results . Hi . hrfparams_2 = cell2mat({ Hi(:) . hrfparams_2 }');
    results . Hi . hrfparams_3 = cell2mat({ Hi(:) . hrfparams_3 }');
    save(   solfile , ...
            'sol','output','params','results','dataparams') 
end

end
    
end

function [ phrf_new ] = shiftlat( phrf_orig , dt )
% SHIFTLAT transforms original HRF parameters into a new set of  parameters 
% that corresponds to a shift of the peak latencies with dt seconds.

phrf_new = phrf_orig;
phrf_new(1) = phrf_new(1) + dt*phrf_orig(3);
phrf_new(2) = phrf_new(2) + dt*phrf_orig(4);

    assert(all(phrf_new(1:4)>0))

end

function [ B , V , abserr , relerr , improve ] = cmtf_multihrf_voxinit_kpe( Y , HS )
% CMTF_MULTIHRF_VOXINIT_KPE computes initial factors in the fMRI's spatial mode
% of the CMTF with multiple HRF estimation. This requires a spatial loading
% that has a Khatri-Rao structure. A Kronecker-Product-Equation (KPE) is
% solved per spatial variable to iteratively find a Khatri-Rao structured
% solution.
%
% INPUTS
% - Y = fMRI data (time points x voxels/ROIs)
% - HS = (1 x K) cell, holding the temporal signatures convolved with every 
% of the K HRFs
% 
% OUTPUTS
% - B = coefficients of the HRF (basis functions x voxels/ROIs)
% - V = spatial signatures of the sources (number of sources x voxels/ROIs)
% - abserr , relerr = absolute , relative error of the approximation inevery voxel/ROI
% - improve = improvement factor of initial estimate to final KPE estimate

%% Check input arguments
    assert(iscell(HS))
    if size(HS,1)~=1 , HS = HS'; end
    
K = length(HS); % number of HRFs
[ Is , R ] = size(HS{1}); % number of time points / number of sources
Iv = size(Y,2); % number of voxels / ROIs

for k = 2 : K
    [ Is_h , R_h ] = size(HS{k});
    assert( Is_h == Is )
    assert( R_h == R )
end

%% Find an initial guess for the spatial factors using the pseudinverse
% (this does not obey the Khatri-Rao constraint)

% -- concatenate the convolved temporal signatures
HSn = cell2mat(HS);

% -- use regression to find the initial guess of the 'betas'
Beta_init = ( HSn \ Y ); % (nT x nB) x nV

%% Correct the initial guess by imposing Khatri-Rao structure
% every column of Si is theoretically the Khatri-Rao product of 2 vectors
% when unfolded into a matrix, this corresponds to rank-1 structure

% -- create empty factor matrices
B = zeros( K , Iv );
V = zeros( R , Iv );

abserr = zeros( 1 , Iv );
relerr = zeros( 1 , Iv );
improve = zeros( 1 , Iv );

% -- go over all columns of the initial guess to find a Khatri-Rao approximation
for iv = 1 : Iv
    % sample the v-th column of the initial guess
    x = Beta_init(:,iv);
    % unfold the column into a matrix of size nS x nB
    X = reshape( x , [ R , K ] );
    % find a rank-1 approximation of X using the SVD
    [ us , sig , ub ] = svds( X , 1 );
    % refine the estimates of the coefficients using KPE
    [ B( : , iv ) , V( : , iv ) , output ] = kpe_refine( Y(:,iv) , HSn , ub , us*sig );
    output.abserr = sqrt(2*output.fval);
    % calculate the errors
    abserr(iv) = output.abserr(end);
    relerr(iv) = abserr(iv) / frob(Y(:,iv));
    improve(iv) = output.abserr(1) / output.abserr(end);
end

end

function [ b , s , output ] = kpe_refine( y , X , bi , si )
% KPE_REFINE solves a Kronecker-Product-Equation (KPE), to obtain more 
% accurate estimates w of the parameters that approximately fit the linear 
% system X*w = y, where w is constrained to have a Kronecker structure

% -- optimization parameters
options = struct;
options.TolFun = 1e-8;
options.TolX = 1e-8;
options.CGMaxIter = length(bi) * length(si);
options.Display = 0;

% -- solve KPE to refine the bi and si variables
[ c , output ] = kpe_nls( X , y , {si,bi} , options );
b = c{2};
s = c{1};

% -- normalize to resolve scale ambiguities
s = s * frob(b);
b = b / frob(b);

end