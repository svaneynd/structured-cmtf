function [ statinf ] = statistical_analysis( basedir , patient , r )
% STATISTICAL_ANALYSIS performs a series of procedures to carry out model
% selection and statistical testing (for a specific sCMTF of rank r of the
% patient's data):
% - select the centroid component and centroid repetition among the 
% clustered solutions of the factorization
% - perform permutation testing
% - compute HRF variability metrics
% 
% INPUTS
% - basedir = base directory of code execution
% - patient = index of the patient whose data need to be processed
% - r = rank (i.e. number components) of the sCMTF
%
% OUTPUTS
% - statinf = structure containing a summary of the permutation-based
% inferece procedure
%
% Author: Simon Van Eyndhoven (simon.vaneyndhoven@gmail.com)

%% Parameters
params = struct;

% -- directories & files
params.basedir = basedir; 

params.data.dir = { 'data\eeg' , 'data\fmri' }; % to load the preprocessed data
params.data.datafile = { 'p%02d_data.mat' , 'p%02d_data.mat' };

params.results.dir = 'cmtf'; % to load the CMTF solution
params.stability.dir = 'stability';

params.results.subdir = 'p%02d';
params.results.version = '';

params.results.solfile = 'normnorm3H%s_cmtf_l1%1.1e_%s_%s_R%d.mat';

params.stat.cmtfdirmask = 'scmtf_l1%1.1e_R%d_%s';
params.stat.patientdir = params.results.subdir;
params.stat.resultdir = fullfile(params.basedir,'inference'); % to save CMTF image(s)
params.stat.basenames = {   {'snpmT' , '_thr', '.img' }, ... % spatial T signature
                            {'snpmF' , '_thr', '.img' }, ... % spatial F signature
                            'hrfB_homog.img', ... % HRF coefficients homogeneity (~distance to centroid)
                            'hrfTC_homog.img', ... % HRF time course homogeneity (~distance to centroid)
                            'hrfB_cluster.img', ...  % HRF coefficients clusters 
                            'hrfTC_cluster.img' };  % HRF time course clusters 

params.mwf.file = 'p%02d_eegmwfenv.mat'; % file containing the filtered EEG envelope

% -- data feature type
params.data.spectype = 'pmtm'; % 'stft' (short-time Fourier transform) / 'pmtm' (Thomson multitaper method)
params.data.boldtype = 'avg'; % 'avg' / 'pca'

% -- decomposition to investigate
params.decomposition.L1pen = 0.001; % component amplitudes (group sparsity)         
params.results.normp = 1; % exponent of the L(p)-norm to be used to normalize the HRF waveforms (e.g. p=2: Euclidean norm)

add_fmriuniq = true;

params.data.patient = patient;
params.decomposition.R = r;

params.stat.cmtfdir = sprintf(params.stat.cmtfdirmask,params.decomposition.L1pen,params.decomposition.R,params.results.version);

% -- results from stability analysis
params.clust.options = struct;
params.clust.solfile = [ 'stability_' , params.results.solfile ];

% -- inference
params.infer.alpha = 0.05; % significance level
params.infer.perm.compute = true;
params.infer.perm.n = 250; % number of surrogate data to generate to find null distribution
params.infer.solfile = [ 'summary_' , params.results.solfile ];

% -- atlas
params.atlas.dir = fullfile(basedir,'\utils\atlas');
params.atlas.file = 'brainnetome_2mm.nii'; 
    
%% Load the atlas
vatlas = spm_vol(fullfile(params.atlas.dir,params.atlas.file));
Vatlas = spm_read_vols(vatlas);
    
%% Initialize parallel computation
if isempty(gcp('nocreate'))
    parpool
end 

fprintf('\n\n--- Computing for patient %d, rank %d\n\n',params.data.patient,params.decomposition.R)
   
    
%% Load the solutions
% -- load the sCMTF factors
solfile = fullfile( params.basedir , params.results.dir , ...
                    sprintf(params.results.subdir,patient) , ...
                    sprintf(params.results.solfile , ...
                        params.results.version , ...
                        params.decomposition.L1pen , ...
                        params.data.spectype , ...
                        params.data.boldtype , ...
                        r) );
                    
eegfmrisol = load(solfile);   

% -- find out whether 'non-causal' HRFs were used
try
    lnc = eegfmrisol.params.data.Lnc; % lnc = number of negative lags at which the HRF differs from zero
catch
    lnc = 0;
end
                
%% Load the data and preprocess
% -- load the data
eeg = load( fullfile( params.basedir , params.data.dir{1} , ...
                    sprintf(params.data.datafile{1},params.data.patient) ) );
                
fmri = load( fullfile( params.basedir , params.data.dir{2} , ...
                    sprintf(params.data.datafile{2},params.data.patient) ) );
                
nsess = length(fmri.Y);
TR = eeg.params.timefreq.TR;
nscans = cellfun(@(x)size(x.avg,1),fmri.Y);
                
% -- preprocess data    
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

eeg = rmfield(eeg,'Z');
fmri = rmfield(fmri,'Y');
    
% -- normalize the data
results = struct;

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

%% Arrange the sCMTF solutions

sol = eegfmrisol.sol;
output = eegfmrisol.output;
nreps = length(sol); % total number of repetitions (runs) of the algorithm

% -- Extract spatial HRF information: basis functions + coefficients + full waveforms
fprintf('* Extracting spatial HRF information...\n')
Lhrf = ceil( 100 / eeg.params.timefreq.TR ); % maintain the samples of the HRF up till 50 seconds
hrfmaps1 = cell(1,nreps); % weighting coefficients B in all ROIs (K x Iv)
hrfmaps2 = cell(1,nreps); % full HRFs in all ROIs (Is x Iv)
hrfmaps3 = cell(1,nreps); % HRF basis functions (Is x K)
parfor solidx = 1 : nreps
    if mod(solidx-1,10)==0
        fprintf('--- solution %d of %d\n',solidx,nreps)
    end
    [ sol{solidx} , H , ~ , Hbasis ] = disambiguatefmri( sol{solidx} , [] , 'posfirst' , params.results.normp );
    hrfmaps1{1,solidx} = [ sol{solidx}.variables.spatial{1} ; sol{solidx}.variables.spatial{2} ; sol{solidx}.variables.spatial{3} ]; % basis coefficients
    hrfmaps2{1,solidx} = H(1:Lhrf,:); % full HRF
    hrfmaps3{1,solidx} = Hbasis(1:Lhrf,:); % full HRF basis function
end
hrfmaps = [hrfmaps1;hrfmaps2;hrfmaps3];
clearvars('hrfmaps1','hrfmaps2','hrfmaps3')
clearvars('sold','H','X')

% -- Compute spatial noise statistics for the fMRI data
fprintf('* Computing spatial noise statistics of fMRI data ...\n')
noiseest = cell(1,nreps);
parfor solidx = 1 : nreps
    if mod(solidx-1,10)==0
        fprintf('--- solution %d of %d\n',solidx,nreps)
    end
    % initialize the noise estimates
    noiseest{solidx}.recerr = zeros(1,size(Y,2));
    noiseest{solidx}.regvar = zeros(params.decomposition.R,size(Y,2));
    % compute the reconstruction errors' variance in every ROI
    % compute the component-wise variance in every ROI
    [recerr,regvar] = estnoise_fmri(Y,sol{solidx},output{solidx},lnc,Lhrf);
    noiseest{solidx}.recerr = recerr; % variance of the residual in each ROI
    noiseest{solidx}.regvar = regvar; % variance of each regressor in the 'local design matrix' Div in each ROI
end


% -- Perform sign correction of the factors     
parfor solidx = 1 : nreps
    % perform sign correction
    [ ~ , scaling , sinkmodes ] = cpdsigncorrect( ...
                                    {sol{solidx}.factors.S,sol{solidx}.factors.G,sol{solidx}.factors.M} , ...
                                    'energy' , 3 );  % use the EEG spatial mode as 'sink'

    sol{solidx}.factors.S = bsxfun( @times , sol{solidx}.factors.S , scaling{1} );
    sol{solidx}.factors.G = bsxfun( @times , sol{solidx}.factors.G , scaling{2} );
    sol{solidx}.factors.M = bsxfun( @times , sol{solidx}.factors.M , scaling{3} );
    
    % apply inverse scaling to the fMRI spatial factor 
        % (note: the norm of the S factor was constrained to be 1, so 'scaling' contains only +1 and -1. for this reason, multiplying is the same as dividing
    sol{solidx}.variables.spatial{4} = bsxfun( @times , sol{solidx}.variables.spatial{4} , scaling{1}' ); 
        
    % concatenate sign-corrected EEG factors with fMRI factors
    try % include the fMRI-only components if they exist
        fmriunc = { sol{solidx}.factors.N , sol{solidx}.factors.P };
    catch
        fmriunc = {};
    end
    sol{solidx} = { sol{solidx}.factors.S , ...
                    sol{solidx}.factors.G , ...
                    sol{solidx}.factors.M , ...
                    sol{solidx}.variables.spatial{4}' , ...
                    sol{solidx}.factors.Lambda_x , ...
                    mean(sol{solidx}.factors.Lambda_y,1) };
    sol{solidx} = [ sol{solidx} , fmriunc ];
end
clearvars('csol','scaling','sinkmodes')

%% Load the stability analysis results                     
clustfile = fullfile( params.basedir , params.stability.dir , ...
                    sprintf(params.results.subdir,params.data.patient) , ...
                    sprintf(params.clust.solfile , ...
                        params.results.version , ...
                        params.decomposition.L1pen , ...
                        params.data.spectype , ...
                        params.data.boldtype , ...
                        params.decomposition.R) );                    
                
if exist(clustfile)
    eegfmristability = load(clustfile);        
else
    warning('Clustering file does not exist yet! Will use a simplified heuristic')
end


%% Load the reference IED time course
% which is the envelope of EEG data that had been filtered using an MWF

% -- load the MWF envelope
mwfenv = importdata(fullfile(params.data.dir{1},sprintf(params.mwf.file,params.data.patient)));

    assert(sum(nscans)==sum(cellfun(@(x)length(x),mwfenv)))
 

%% Model selection: detect which component represents the IED the best

try   
ncomps = cellfun( @(x)length(x) , eegfmristability.similarity.cliquegroups );

regressors = zeros(sum(nscans),sum(ncomps));
regfeats = zeros(3,sum(ncomps));
runs = cell(1,sum(ncomps));

% -- loop over all 'cliques' and 'groups' to check which component has the
% greatest similarity
selmanual = false; % follow the automatic criterion, or manually select the best component to model the IED?

regidx = 0;
figh = figure('position',[ 200 0 1800 900 ]);
for clique = 1:numel(eegfmristability.similarity.cliquegroups)
    % -- find which groups belong to this clique
    groupidx = eegfmristability.similarity.cliquegroups{clique};
    
for g = 1:numel(groupidx)
    
    regidx = regidx + 1;

    % -- find out which in which runs this component appears
    group = groupidx(g);        
    m = eegfmristability.similarity.M(:,group);
    runs{regidx} = eegfmristability.similarity.cliqueruns{clique}{g};
    
    abscompidx = m(runs{regidx}) + (runs{regidx}-1)*params.decomposition.R; % absolute index of the involved components, in the similarity matrix
    
    % -- store the time course of this component
    TC = zeros(size(sol{1}{1},1),numel(runs{regidx})); % matrix holding all time courses    
    for ridx = 1 : numel(runs{regidx})
        TC(:,ridx) = sol{runs{regidx}(ridx)}{1}(:,m(runs{regidx}(ridx)));
    end
    
    % determine centroid component and centroid repetition
    Csim = eegfmristability.similarity.cont(abscompidx,abscompidx);
    [~,centroididx] = max(sum(Csim,2)-1); % minus 1: the diagonal term doesn't count (although subtraction doesn't change the ordering)
    regressors(:,regidx) = TC(:,centroididx);
    centroidrun(regidx) = runs{regidx}(centroididx);
    centroidcomp(regidx) = m(runs{regidx}(centroididx));

    
    % -- store several features of the candidate regressor    
    % 1. similarity with the MWF envelope
    simscore = zeros(1,nsess);
    for sess = 1 : nsess  
        range = sum(nscans(1:sess-1)) + 1 : sum(nscans(1:sess));
        simscore(sess) = corr( regressors(range,regidx) , mwfenv{sess} );
    end
    regfeats(1,regidx) = mean(simscore);
    
    % 2. number of runs in which the component appears
    regfeats(2,regidx) = length(runs{regidx});     
        
    % -- plot this candidate regressor
    figure(figh);
    subplot(ceil(sum(ncomps)/2),2,regidx)
    hold('on')
    if dohrfconvol
    regressor_temp = mat2cell( regressors(:,regidx) , cellfun(@(x)length(x),mwfenv) , 1 );
    regressor_temp = cellfun( @(x)fftfilt(hrf,x) , regressor_temp , 'uniformoutput', false);
    regressor_temp = cell2mat(regressor_temp);
    else
        regressor_temp = regressors(:,regidx);
    end
    plot(struct_normalize(regressor_temp,[]),'k','linewidth',1)
    plot(struct_normalize(cell2mat(mwfenv),[]),'color',[0.9, 0.7 , 0.7],'linewidth',1)
    titlestr = sprintf('comp %d : sim. = %2.2f, %d runs, L1 = %2.2f, %d ROIs',regidx,regfeats(1,regidx),regfeats(2,regidx),regfeats(3,regidx));
    title(titlestr,'fontsize',10,'interpreter','latex')
    axis('tight')
    
end
end


% -- select the component that best models the IEDs (based on highest correlation with EEG MWF envelope)    
allcomps = cell2mat(eegfmristability.similarity.cliquegroups);
minnruns = 10; % minimal number of runs for a component to be a candidate
% objf = regfeats(1,:) .* regfeats(2,:) .* (regfeats(2,:)>=minnruns) .* (regfeats(1,:)>max(regfeats(1,:))-0.10);
objf = regfeats(1,:) .* (regfeats(2,:)>=minnruns);
while ~any(objf>0)
    minnruns = minnruns-1;
    objf = regfeats(1,:) .* (regfeats(2,:)>=minnruns);
    % objf = regfeats(1,:) .* regfeats(2,:) .* (regfeats(2,:)>=minnruns) .* (regfeats(1,:)>max(regfeats(1,:))-0.10);
    if minnruns == 1 % condition to ensure that exiting from the loop is possible
        objf(objf==max(objf)) = 1;
    end
end
    [~,iedcomp] = max(objf);

    if selmanual
    pause;
    prompt = {sprintf('Enter a new index for the IED component\nEstimated IED component: %d',iedcomp)};
    prompttitle = 'Choose the IED component';
    promptdims = [1 70];
    definput = {sprintf('%d',iedcomp)};
    iedcomp = inputdlg(prompt,prompttitle,promptdims,definput);

    close;

iedcomp_relidx = str2num(iedcomp{1}); % index is relative within all evaluated components
else
    iedcomp_relidx = iedcomp;
end

catch % select the component and the run with the highest correlation to MWF envelope, 
        % regardless of clustering structure over runs of the algorithm
    
    centroidrun = 0;
    centroidcomp = 0;
    iedcomp_relidx = 1;
    
    all_timesignatures = zeros(size(sol{1}{1},1),nreps*r);
    
    for solidx = 1 : nreps
        all_timesignatures(:,(solidx-1)*r+1:solidx*r) = sol{solidx}{1};; % first factor is the temporal factor
    end
    
    [~,maxidx] = max(corr(all_timesignatures,cell2mat(mwfenv)));
    
    centroidrun(iedcomp_relidx) = ceil(maxidx/r);
    centroidcomp(iedcomp_relidx) = mod(maxidx-1,r)+1;
    
end

%% Perform permutation-based inference on the spatial activation

% -- conduct the inference on the coupled part of the decomposition only,
% i.e. without the fMRI-only factors [[N,P]]
if add_fmriuniq
    Yres = Y - (sol{ centroidrun(iedcomp_relidx) }{end-1} * sol{ centroidrun(iedcomp_relidx) }{end}');
    fprintf('\nsubtracted fMRI-only part\n')
else
    Yres = Y;
end

% -- create a structure holding all the necessary information in order to
% conduct permutation-based inference (using GLM-style notation: y ~ X*Beta)
regdata = struct;
regdata.X0 = sol{ centroidrun(iedcomp_relidx) }{1}; % original time signatures S
regdata.H0 = hrfmaps{2, centroidrun(iedcomp_relidx) }; % original HRFs
regdata.Beta0 = sol{ centroidrun(iedcomp_relidx) }{4}; % original spatial loadings
regdata.Beta0n = sqrt( noiseest{ centroidrun(iedcomp_relidx) }.recerr(:) ); % original spatial noise
regdata.output = output{ centroidrun(iedcomp_relidx) }; % original output summary
regdata.regvar = noiseest{ centroidrun(iedcomp_relidx) }.regvar;
regdata.recerr = noiseest{ centroidrun(iedcomp_relidx) }.recerr;
regdata.lnc = lnc; % number of non-causal HRF samples

% -- conduct permutation-based inference, to create a pseudo-T map and a
% pseudo-F map for each component in the centroid run
% (For a test on a single regressor: F = T.^2)
statinf = struct;
statinf.Tmap = struct;
statinf.Fmap = struct;

[ statinf.Tmap.res , statinf.Tmap.masks , statinf.Fmap.res , statinf.Fmap.masks ] = permutation_inference( Yres , regdata , params );

%% Quantify HRF variability using different metrics
statinf.hrf.res = struct;
[ entropies , enttypes ] = compute_hrfentropies( hrfmaps{1, centroidrun(iedcomp_relidx) } , hrfmaps{2, centroidrun(iedcomp_relidx) } );
for typeidx = 1 : numel(enttypes)
    statinf.hrf.res.(enttypes{typeidx}).e0 = entropies{typeidx}(:);
end

%% Create and save images of SnPMs, and HRF variability metrics
mkdir(fullfile(  params.stat.resultdir , ...
                        sprintf(params.stat.patientdir,params.data.patient) , ...
                        params.stat.cmtfdir ))
                    
% -- prepare a template volume
vtemplate = struct;
vtemplate.mat = vatlas.mat;
vtemplate.dim = vatlas.dim;
vtemplate.pinfo = [1;0;0];
vtemplate.dt = [16 0];
vtemplate.n = [1 1];
vtemplate.descrip = 'statistical image';

% -- construct a whole-brain image out of the pseudo-T maps
save_statmaps( statinf.Tmap , vtemplate , centroidcomp(iedcomp_relidx) , params , Vatlas , 1 , 'snpmT' ); % alternatively: clustsize = 350

% -- construct a whole-brain image out of the pseudo-F maps
try
save_statmaps( statinf.Fmap , vtemplate , centroidcomp(iedcomp_relidx) , params , Vatlas , 1 , 'snpmF' ); % alternatively: clustsize = 350
catch
end

% construct and save whole-brain images out of the HRF entropies
save_entropymaps( statinf.hrf , vtemplate , params , Vatlas , 1 ); % alternatively: clustsize = 350

%% Save the summary statistical metrics and other relevant data
% -- results of inference
save(fullfile(  params.stat.resultdir , ...
                        sprintf(params.stat.patientdir,params.data.patient) , ...
                        params.stat.cmtfdir , ...
                        'statresults.mat' ),...
        'statinf','params')
    
% -- other signatures and useful information    
summaryfile = fullfile( params.basedir , params.results.dir , ...
                    sprintf(params.results.subdir,params.data.patient) , ...
                    sprintf(params.infer.solfile , ...
                        params.results.version , ...
                        params.decomposition.L1pen , ...
                        params.data.spectype , ...
                        params.data.boldtype , ...
                        params.decomposition.R) );
summary = struct;

summary . hrf . coeff = hrfmaps(1,:);
summary . hrf . tc = hrfmaps(2,:);
summary . hrf . entropy = statinf.hrf.res;

summary . sol = sol;

summary . centroid . run = centroidrun;
summary . centroid . comp = centroidcomp;

summary . iedcomp = iedcomp_relidx;

try
summary . stability = eegfmristability;
catch
end

summary . inference = statinf;

summary . solparams = eegfmrisol.params;

save(summaryfile,'summary')

end