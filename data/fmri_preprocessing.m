function [ Y , params ] = fmri_preprocessing( basedir , patient , fmri_sessiondata , R , TR )
% fMRI_PREPROCESSING processes raw BOLD data by applying filtering,
% removing confounds, and averaging the time seriees within predefined 
% parcels of a brain atlas.
%
% INPUTS
% - basedir = base directory of code execution
% - patient = index of the patient whose data need to be processed
% - fmri_sessiondata = (#sessions x 1) cell array, containing the volume
% information of all MR scans per session / acquisition run
% - R = (#sessions x 1) cell array, containing a (time points x regressors)
% matrix containing time courses of nuisance regressors that should be
% removed from the BOLD data
% - TR = repetition time of fMRI in s
%
% OUTPUTS
% - Y = (#sessions x 1) cell array, containing a (time points x parcels)
% matrix of averaged BOLD signals in all of the atlas's parcels (ROIs)
%
% Dependency: SPM version 8 or higher needed
%
% Author: Simon Van Eyndhoven (simon.vaneyndhoven@gmail.com)

%% Parameters
params = struct;

% -- files & directories
params.basedir = basedir;
params.data.targetdir = basedir;

params.data.fmridir = 'data\fmri';
mkdir(params.data.targetdir,params.data.fmridir)

params.data.eegdir = 'data\eeg';
params.data.mwfenvfile = 'p%02d_eegmwfenv.mat';

params.data.filename = 'p%02d_data.mat'; % filename to save the parcellated BOLD data

% -- atlas
params.atlas.dir = fullfile(basedir,'utils\atlas');
params.atlas.file = 'brainnetome_2mm.nii'; % file of Brainnetome atlas for parcellation

% -- white matter (WM) and cerebrospinal fluid (CSF) masks
params.masks.dir = params.atlas.dir;
params.masks.wmfile = 'white_2mm.nii';
params.masks.csffile = 'csf_2mm.nii';
params.masks.thr = [ 0.99 , 0.750 ]; % threshold for tissue probability
params.masks.npc = 5; % number of principal components of the masks' time series to retain (per mask)
params.masks.include = [true true]; % wm / csf

% -- band-pass filter
params.filter.band = [1/128 0.2]; % high-pass cutoff of 1/128 Hz: default value in SPM

% -- confounds
params.confounds.inclderivsq = true; % include the derivative and squared signal of each confound and add these to the set of confounds
params.confounds.orthmwf = false; % if set to true: project the confounders' time courses on the null space of the MWF envelope, before regressing out of the BOLD data

% -- PCA compression
params.compress.dopca = false;
params.compress.varexp = 0.90; % percentage of the variance that needs to be explained by the retained singular vectors

%% Load the atlas and the masks
% -- atlas
params.atlas.vatlas = spm_vol(fullfile(params.atlas.dir,params.atlas.file));
[Vatlas,xyz.atlas] = spm_read_vols(params.atlas.vatlas);
params.atlas.rois = unique(Vatlas(:)); % find the indices of all regions of interest (ROIs), aka parcels
params.atlas.rois = params.atlas.rois(params.atlas.rois>0);

% -- masks
vcsf = spm_vol(fullfile(params.masks.dir,params.masks.csffile));
vwm = spm_vol(fullfile(params.masks.dir,params.masks.wmfile));
    assert(frob(vcsf.mat-params.atlas.vatlas.mat)/frob(vcsf.mat)<1e-6)
    assert(frob(vwm.mat-params.atlas.vatlas.mat)/frob(vwm.mat)<1e-6)
[Vcsf,xyz.csf] = spm_read_vols(vcsf);
[Vwm,xyz.wm] = spm_read_vols(vwm);

% -- apply a threshold to convert probabilities into binary masks
VMwm = zeros(size(Vwm));
VMwm(Vwm>params.masks.thr(1)) = 1;   
VMwm(Vatlas>0) = 0;
VMcsf = zeros(size(Vcsf));
VMcsf(Vcsf>params.masks.thr(2)) = 1;
VMcsf(Vatlas>0) = 0;

%% Load MWF envelope
fprintf('* load MWF-filtered EEG envelope\n')
load( fullfile(params.data.targetdir,params.data.eegdir,...
        sprintf(params.data.mwfenvfile,patient)) );     
  
%% Load and preprocess the BOLD data voxelwise
% Per session of the data, perform the following steps:
% - load session-specific BOLD data
% - regress out confounding variables
% - band-pass filtering
% - data compression: average BOLD signal / PCA compression in atlas ROIs

% -- extract timing-related info
nsess = numel(fmri_sessiondata); % number of sessions
L = cellfun(@(x)size(x(1).vol,1),fmri_sessiondata); % lenghts in samples of all sessions

% -- loop over all sessions to process the BOLD data
Y = cell(nsess,1);

for sess = 1:nsess
    fprintf('\n/ session %d /\n',sess)

    % -- load the data
    fprintf('* load data\n')
    Y{sess} = zeros( L(sess) , prod(fmri_sessiondata{sess}.vol(1).dim) );

    % load the data of every scan 
    for scan = 1:length(fmri_sessiondata{sess}.vol)
        yscan = spm_read_vols( fmri_sessiondata{sess}.vol(scan) );
        Y{sess}(scan,:) = tens2vec(yscan,1:3)'; % structure the 3D image as a long row vector
    end

    % -- compute confounding regressors that should be removed
    fprintf('* determine confounds\n')

        % extract the original regressors
        Rn = struct_normalize(R{sess},[]);
            assert(L(sess) == size(Y{sess},1))

        % compute WM and CSF regressors
        Rroi = [];
        if params.masks.npc > 0
            if params.confounds.orthmwf
                M = [ eegmwfenv{sess} , circshift(eegmwfenv{sess},1,1) , circshift(eegmwfenv{sess},-1,1) , ones(L(sess),1) ];
            else
                M = ones(L(sess),1);
            end
            [M,~] = qr(M,0);
            if params.masks.include(1) % include WM confounds
            Nr = Y{sess}(:,VMwm==1);
            Nr = Nr - M*M'*Nr;
            [Ur,sr,~] = svds(Nr,params.masks.npc);
            Rroi = [ Rroi , Ur ];
            end
            if params.masks.include(2) % include CSF confounds
            Nr = Y{sess}(:,VMcsf==1);
            Nr = Nr - M*M'*Nr;
            [Ur,sr,~] = svds(Nr,params.masks.npc);
            Rroi = [ Rroi , Ur ];
            end
        Rn = [ Rn , Rroi ];
        clearvars('Nr','Ur','sr')
        end

        % add derivative and squared signal of each confound and add these to the set of confounds
        if params.confounds.inclderivsq
            Rn = [ Rn , tderiv(Rn')' , Rn.^2 ];
        end        

        % add a constant confound to model baseline
        Rn = [ Rn , ones(L(sess),1) ];

        % store the extended regressor matrix
        Rn = struct_normalize(Rn,[]);
        Rext{sess} = Rn;


    % -- remove the confounding regressors
    fprintf('* removing confounds\n')
    [Rn,~] = qr(Rn,0);
    P_conf = eye(L(sess)) - Rn*Rn';
    Y{sess} = P_conf*Y{sess};

    % -- band-pass filtering      
    fprintf('* band-pass filtering\n')
    nparts = 100;
    L_part = ceil(size(Y{sess},2)/nparts);
    V_comp = cell(nparts,1);
    for part = 1:nparts % process the voxels in batches
        range = (part-1)*L_part+1 : part*L_part;
        range = range(range<=size(Y{sess},2));
        % <!> filtfilt produces large artifacts (especially at boundaries) --> don't use
        Y{sess}(:,range) = conn_filter(TR,params.filter.band,Y{sess}(:,range));
    end  

    % -- compression based on atlas
    fprintf('* atlas-based compression\n')
    Yc = struct;
    Yc.avg = zeros(L(sess),length(params.atlas.rois)); % average BOLD signal within ROI
    Yc.pca = zeros(L(sess),length(params.atlas.rois)); % 1st principal component of BOLD signal within ROI
    Yc.avghomog = zeros(length(params.atlas.rois),1);
    Yc.pcahomog = zeros(length(params.atlas.rois),1);
    for iv = 1:length(params.atlas.rois)
        % select the BOLD signals in the current ROI
        Yroi = Y{sess}(:,Vatlas==params.atlas.rois(iv));            
        params.atlas.nvox(iv) = size(Yroi,2);
        % average BOLD signal in the ROI
        Yc.avg(:,iv) = mean(Yroi,2);
        % compute the homogeneity of the average in the ROI
        Yc.avghomog(iv) = mean( abs(mean(Yroi,2))./std(Yroi,[],2) );
        % 1st principal component of the BOLD signals in the ROI
        [Uroi,Sigmaroi,Vroi] = svd(Yroi);
        Sigmaroi = Sigmaroi(1:min(size(Sigmaroi,1),size(Sigmaroi,2)),1:min(size(Sigmaroi,1),size(Sigmaroi,2)));
        Uroi = Uroi(:,1);
        Vroi = Vroi(:,1);
        % correct the sign of the 1st PC
        if median( cossimil( Yc.avg(:,iv) , Uroi*Sigmaroi(1,1) )) < 0 % / median(Vroi) < 0
            % scale the PC by -1 to correct the sign, and by
            % 1/sqrt(nvox) such that the scaling corresponds to the
            % scaling of the average BOLD signal
            Yc.pca(:,iv) = -Uroi*Sigmaroi(1,1)/sqrt(params.atlas.nvox(iv));
        else
            % scale the PC by 1/sqrt(nvox) such that the scaling 
            % corresponds to the scaling of the average BOLD signal
            Yc.pca(:,iv) = Uroi*Sigmaroi(1,1)/sqrt(params.atlas.nvox(iv));
        end
        % compute the homogeneity of principal component in the ROI
        Yc.pcahomog(iv) = Sigmaroi(1,1)^2 / sum( diag(Sigmaroi).^2 );
    end

    Y{sess} = Yc;     

end
    
%% (optional:) Data-driven compression of the (time points x ROIs) data using PCA
if params.compress.dopca
    % Note: this compression deals with the whole-brain multivariate BOLD data
    % and is distinct from the PCA-based extraction of a single time series
    % in a ROI, as done above.
    
fprintf('\n/ PCA compression (keeping %2.2f %% of variance) /\n',100*params.compress.varexp)

% -- compute the SVD
Z = [];
for sess = 1:nsess
    Z = [ Z ; Y{sess}.avg ];
    Y{sess} = rmfield(Y{sess},{'avg','pca'});
end
[Uy,Sy,Vy] = svd(Z,'econ'); 
clearvars('Z')
sy = diag(Sy);

% -- find an appropriate cutoff in the singular value spectrum 
% Note: due to the frequency-domain filtering, degrees of freedom are
% removed from the signal, and this can be seen (at least for
% single-session data) by the fact that the data are already lower
% rank)
idxc = find(cumsum(sy.^2)/sum(sy.^2)>=params.compress.varexp,1);
Uy = Uy(:,1:idxc);
Sy = Sy(1:idxc,1:idxc);
Vy = Vy(:,1:idxc);
fprintf('* compression: %2.2f %% of singular values kept\n', 100*idxc/length(sy))    

% -- rearrange the data into a cell array (for compatibility with other code)
Uy = mat2cell(Uy,L,idxc);
Y = cell(nsess,1);
for sess = 1:nsess
    Y{sess}.avg = Uy{sess};
    Y{sess}.Sig = sy(1:idxc)';
end

% -- store the necessary building blocks for decompression
Ysvd = struct;
Ysvd.V = Vy;
Ysvd.idxc = idxc;
Ysvd.sv = sy;

else
Ysvd = []; 

end
    
%% Save the processed BOLD data
save( fullfile(params.data.targetdir,params.data.fmridir,...
        sprintf(params.data.filename,patient) ) ,...
        'Y','R','Rext','eegmwfenv','params','Vatlas' )

end

function [ dq ] = tderiv( q , dt )
% TDERIV compute the temporal derivative (with time step sampling period
% dt) of the (variables x time points) data matrix q.

if nargin < 2 
    dt = 1;
end

dqbulk = (q(:,3:end)-q(:,1:end-2))/(2*dt); % bulk
dqstart = (q(:,2)-q(:,1))/dt; % begin
dqend = (q(:,end)-q(:,end-1))/dt; % end
dq = [dqstart,dqbulk,dqend]; 

end