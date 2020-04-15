function [ similarity ] = assess_stability_scmtf( basedir , patient , r , visualize )
% ASSESS_STABILITY_SCMTF performs a structured clustering procedure to
% analyze in which runs (repetitions) of the sCMTF optimization similar
% solutions for the factors are found.
% 
% INPUTS
% - basedir = base directory of code execution
% - patient = index of the patient whose data need to be processed
% - r = rank (i.e. number components) of the sCMTF
% - visualize = boolean, indicating whether or not to show all clustered
% components
%
% OUTPUTS
% - similarity = structure containing the results of the clustering
% procedure:
% - 
%
% Author: Simon Van Eyndhoven (simon.vaneyndhoven@gmail.com)

addpath('C:\Users\svaneynd\Documents\Studies\PhD\Research\Factorization clustering')

if nargin < 4
    visualize = false;
end

%% Parameters
params = struct;

% -- directories & files
params.basedir = basedir;

params.data.dir = { 'data\eeg' , 'data\fmri' }; % directories where EEG and fMRI data are stored
params.data.datafile = { 'p%02d_data.mat' , 'p%02d_data.mat' }; % filenames for EEG and fMRI data

params.stability.dir = 'stability';

params.scmtf.dir = 'cmtf';

mkdir( params.basedir , params.stability.dir )

params.results.subdir = 'p%02d';   
params.results.version = '';    
add_fmriuniq = true;

params.results.solfile = 'normnorm3H_cmtf_l1%1.1e_%s_%s_R%d.mat';%'3Hv4_cmtfBN_krgs%1.1e_%s_%s_randgevd_R%d.mat'; % CMTF filename

% -- data feature type
params.data.spectype = 'pmtm'; % 'stft' (short-time Fourier transform) / 'pmtm' (Thomson multitaper method)
params.data.boldtype = 'avg'; % 'avg' / 'pca'

% -- decomposition to investigate
params.decomposition.R = r;
params.decomposition.L1pen = 0.001; % L1 penalty on component amplitudes (group sparsity) 
params.decomposition.modes = 1 : 4; % factor modes to take into account for stability analysis (temporal/spectral/spatial/amplitudes)

% -- clustering of the solutions of the factorization
params.clust.options = struct;
params.clust.options.Xin = 'thresh';
params.clust.options.maxcost = 0.15; % similarity threshold = 1 - maxcost
params.clust.options.roundperm = 'avglink';
params.clust.solfile = [ 'stability_' , params.results.solfile ];
 

    
%% Load the solution
solfile = fullfile( params.basedir , params.scmtf.dir , ...
                    sprintf(params.results.subdir,patient) , ...
                    sprintf(params.results.solfile , ...
                        params.decomposition.L1pen , ...
                        params.data.spectype , ...
                        params.data.boldtype , ...
                        params.decomposition.R) );
                
eegfmrisol = load(solfile);    

% find out the degree of non-causality
try
    lnc = eegfmrisol.params.data.Lnc;
catch
    lnc = 0;
end
                
%% Load the data
eegdata = load( fullfile( params.basedir , params.data.dir{1} , ...
                    sprintf(params.data.datafile{1},patient) ) );

% -- concatenate data across sessions
X = [];
for sess = 1:length(eegdata.Z)
    sessdata = eegdata.Z{sess}.(params.data.spectype);
    % adjust the length if needed: fMRI data dictates which length
    % should be taken
    excessL = eegdata.params.timefreq.L{sess}.eeg - eegdata.params.timefreq.L{sess}.fmri;
    if excessL > 0 % EEG data is longer than fMRI data : apply truncation
        sessdata = sessdata( 1:eegdata.params.timefreq.L{sess}.eeg-excessL , : , : );
    elseif excessL < 0 % EEG data is shorter than fMRI data : apply zero-padding
        sessdata = cat( 1 , sessdata , zeros(-excessL,size(sessdata,2),size(sessdata,3)) );
    else % EEG and fMRI data are already of equal length for this session, which is good
    end
    X = cat( 1 , X , sessdata );
end    

% -- normalize the EEG data
results = struct;
results.mu = mean(X,1);
X = bsxfun( @minus , X , results.mu );
[X,results.fnorms,results.scaling] = normalize_array( X , [2,3] );

X = X / frob(X);

%% Arrange the sCMTF factors in the right format         

% -- arrange the solution
sol = eegfmrisol.sol;
output = eegfmrisol.output;
nsols = length(sol); % total number of runs

% -- Convert to cell array if needed (and perform dimension check), to use
% in stability analysis (keep the original solutions to have access to the
% HRFs)
[ ~, facnames ] = convert_cmtffac_kr( sol{1} , true );
sol = cellfun( @(x)convert_cmtffac_kr(x,true) , sol , 'uniformoutput',false );


%% Compute relative approximation error of the sCMTF in all runs
% -- compute summary metrics
approxerr = cell(length(sol),1);
ccd = zeros(length(sol),1);
algit = zeros(length(sol),1);
for s = 1:nsols % go over all runs of the algorithm
    approxerr{s} = output{s}.abserr; % approximation error
    ccd(s) = corcondia(X, { sol{s}{1}*diag(sol{s}{5}) , sol{s}{2} , sol{s}{3} } ); % core consistency diagnostic
    algit(s) = output{s}.iterations; % number of iterations until convergence of the algorithm
end
approxerr = cell2mat(approxerr);

% -- plot error
figure
subplot(2,1,1)
plot(approxerr(:,1:2)),ylim([0 1])
title('relative fitting error')
xlabel('runs of the algorithm')
legend('EEG','fMRI')
subplot(2,1,2)
plot(approxerr(:,3:4)),ylim([0 1.1*max(max(approxerr(:,3:4))+eps)])
title({'L1 penalty on';'component amplitudes'})
xlabel('runs of the algorithm')
legend('EEG','fMRI')

% -- plot CCD and number of iterations
figure
subplot(2,1,1)
plot(ccd)
title('Core Consistency Diagnostic')
xlabel('runs of the algorithm')
ylim([0 100])
subplot(2,1,2)
plot(algit)
title({'number of iterations';'until convergence'})
xlabel('runs of the algorithm')


fprintf('\n>>> Press any key to continue.\n')
pause;

%% Assess the stability of the decomposition
% -- compute a simple similarity matrix (for visualization only) using 'cpderr'
A = eye(length(sol)*params.decomposition.R);
for run1 = 1:length(sol)
    cols = (run1-1)*params.decomposition.R+1:run1*params.decomposition.R;
    for run2 = run1+1:length(sol)
        rows = (run2-1)*params.decomposition.R+1:run2*params.decomposition.R;
        [err,P,~] = cpderr(sol{run1}(params.decomposition.modes),sol{run2}(params.decomposition.modes));
            assert(size(P,1)==params.decomposition.R); assert(size(P,2)==params.decomposition.R)
        sc = 1 - prod(err)^(1/length(err));
        A(rows,cols) = sc*P;
        A(cols,rows) = A(rows,cols)';
    end
end
figure,imagesc(A),title({'Similarity matrix','(#runs x #components)x(#runs x #components)'})

fprintf('\n>>> Press any key to continue.\n')
pause;

close('all') 

% -- perform structured clustering
[ similarity , solinfo , csol ] = factorizationclustering( sol , params.decomposition.modes , params.clust.options );

%% Visualize the signatures of runs that were identified as belonging to the same clique (~local optimum)

if visualize

csol = cellfun( @(x) cpdsigncorrect(x(1:4),'energy',3) , csol , 'uniformoutput', false );

% -- disambiguate the spatial signatures of the fMRI factors
for solidx = 1 : length(eegfmrisol.sol)
    eegfmrisol.sol{solidx} = disambiguatefmri( eegfmrisol.sol{solidx} );
end

foi = eegdata.params.timefreq.foi; % frequencies of interest
toi = (0:size(csol{1}{1},1)-1)*eegfmrisol.dataparams.eeg.timefreq.TR;
Lhrf = 30; % round( 50 / eegfmrisol.dataparams.eeg.timefreq.RT ); % limit HRF to 30 seconds

    modenames = {   'spectral' , ...
                    'spatial (channels)' , ...
                    'time course' , ...
                    'spatial (ROIs)' , ...
                    'occ./ovl.' , ...
                    };
    relspace = [ 5 3 6 5 0.5 ];
    
% loop over supercliques = clusters of components that appear together in
% many runs
for supercliqueclique = 1:numel(similarity.cliquegroups)
    
    % -- find which groups belong to this clique
    groupidx = similarity.cliquegroups{supercliqueclique}; % group ~ clique
    
    % -- compute the minimal/median/maximal overlap between any two groups in the clique and report it
    A = similarity.cliqueruns{supercliqueclique};
    ov = cellfun(@(x,y)numel(intersect(x,y)),repmat(A',1,length(A)),repmat(A,length(A),1));
    
    % -- create a new figure to display this clique
    figparams = cmtf_figure_template( numel(groupidx) , modenames , relspace , ...
                    sprintf(['Superclique (cluster) %d: %d components (pairwise common runs = %d -- %d)'],supercliqueclique,numel(groupidx),min(ov(:)),max(ov(triu(ov,1)>0))) );
    
% -- collect the HRFs of all runs in the clique and plot
cidx = A{1}; 
for clique = 1:length(A)
    cidx = intersect(cidx,A{clique}); % i.e. find the runs in which _all_ components of this clique appear
end

subplot(cmtf_figure_gethandle( figparams , 'hrf' )), hold('on')
ylims = [];
hrfcolours = {'r','k','b'};
for hidx = 1 : 3 % plot all 3 HRF basis functions
    
    HRF = zeros(size(eegfmrisol.sol{1, 1}.factors.(sprintf('HRFtoep_%d',hidx)),1),length(cidx));
    for ridx = 1:length(cidx)
        HRF(:,ridx) = eegfmrisol.sol{cidx(ridx)}.factors.(sprintf('HRFtoep_%d',hidx));
        c = sign(cossimil(HRF(:,1),HRF));
        HRF = bsxfun(@times,HRF,c);
    end
    correctsign = sign(mean( sign(HRF) .* (HRF.^2) ,1)); % sign(mean(HRF))
    HRF = bsxfun(@times,HRF,correctsign);
    HRF = struct_normalize(HRF,[]);

    plot(toi,HRF,'color',hrfcolours{hidx}),axis('tight')    
    ylims(hidx,:) = ([min(HRF(:)) - 0.1*(max(HRF(:))-min(HRF(:))) , max(HRF(:)) + 0.1*(max(HRF(:))-min(HRF(:)))]);
    
end
xlim([0 , Lhrf])
ylim([min(ylims(:,1)) , max(ylims(:,2))])

% -- every clique = 1 template component in the superclique
for clique = 1:numel(groupidx)

    % -- find out which in which runs this component appears
    group = groupidx(clique);        
    m = similarity.M(:,group);
    runs = similarity.cliqueruns{supercliqueclique}{clique};

    % -- stack the signatures of all runs in this clique
    S = zeros(size(csol{1}{1},1),numel(runs));
    G = zeros(size(csol{1}{2},1),numel(runs));
    M = zeros(size(csol{1}{3},1),numel(runs));
    V = zeros(size(eegfmrisol.sol{1}.variables.spatial{4}',1),numel(runs)); %zeros(size(csol{1}{4},1),numel(runs));
    Lambda = zeros(2,numel(runs));
    
    for ridx = 1:numel(runs)
        S(:,ridx) = csol{runs(ridx)}{1}(:,m(runs(ridx)));
        G(:,ridx) = csol{runs(ridx)}{2}(:,m(runs(ridx)));
        M(:,ridx) = csol{runs(ridx)}{3}(:,m(runs(ridx)));
        V(:,ridx) = eegfmrisol.sol{runs(ridx)}.variables.spatial{4}(m(runs(ridx)),:)'; % csol{runs(ridx)}{4}(:,m(runs(ridx)));
        Lambda(1,ridx) = eegfmrisol.sol{runs(ridx)}.factors.Lambda_x(m(runs(ridx)));
        Lambda(2,ridx) = mean(eegfmrisol.sol{runs(ridx)}.factors.Lambda_y(m(runs(ridx))),1);
    end
    
    S = diag(1./diag(eegfmrisol.results.scaling_t{1}))*S;
    G = diag(1./diag(eegfmrisol.results.scaling_t{2}))*G;
    M = diag(1./diag(eegfmrisol.results.scaling_t{3}))*M;
    
    S = struct_normalize(S,[]);
    G = struct_normalize(G,[]);
    M = struct_normalize(M,[]);
    V = struct_normalize(V,[]);

    % -- plot the overlaid signatures 
    
    % 1. occurrence
    subplot(subplot(cmtf_figure_gethandle( figparams , 'occ' , clique )))

    occurrence = [  length(runs) ; min(ov(clique,:)) ];
    xlab = { '' , '' }; %{sprintf('%d runs',length(runs)),sprintf('min %d runs overlap',min(ov(g,:)))};
    bar(1:2,occurrence,0.7,'FaceColor',[1 0 0],'Linestyle','none')
    set(gca,'Xticklabel',xlab)
    ylim([0 sum(nsols)])

    % compute the coherence of this clique
    groupcomps = cumsum([0;similarity.ncomps(1:end-1)']) + m;
    groupcomps = groupcomps(runs);
    cluscoh = clustercoherence(similarity.cont,groupcomps); % ~ ICASSO terminology


    % add info about coherence and norms
    set( figparams.SourceTitlePane.handles{clique} , 'string' , ...
            sprintf(    'coh = %2.2f\nampE = %2.2f\n... +/-- %2.2f\nampF = %2.2f\n... +/-- %2.2f',cluscoh,median(Lambda(1,:)),std(Lambda(1,:)),median(Lambda(2,:)),std(Lambda(2,:))) , ...%{ sprintf('coh(%d) = %2.2f',g,cluscoh) , sprintf('med. $\\mid\\mid s_%d\\mid\\mid$ = %2.2f',g,prctile(compnorms_g,50)) } ,...
            'fontsize' , 10 )

    % 1. time
    subplot(cmtf_figure_gethandle( figparams , 'time' , clique ))
    plot(S),axis('tight')%set(gca,'Xtick',1:numel(EEG.params.data.cidx),'Xticklabel',EEG.params.data.chans(EEG.params.data.cidx),'Xticklabelrotation',45),xlabel('channels')

%     title({'norms = ';sprintf('%2.2f -- %2.2f -- %2.2f',prctile(compnorms_g,25),prctile(compnorms_g,50),prctile(compnorms_g,75))})

    % 2. frequency 
    subplot(cmtf_figure_gethandle( figparams , 'spec' , clique ))
    plot(foi,G),axis('tight'),%set(gca,'Xtick',EEG.params.freq.cfg.foi,'Xticklabel',EEG.params.freq.cfg.foi),xlabel('frequency (Hz)')

    % 3. EEG space 
    subplot(cmtf_figure_gethandle( figparams , 'spatial (channels)' , clique ))
    plot(M),axis('tight'),%set(gca,'Xtick',EEG.params.freq.cfg.foi,'Xticklabel',EEG.params.freq.cfg.foi),xlabel('frequency (Hz)')
   
    % 4. fMRI space 
    subplot(cmtf_figure_gethandle( figparams , 'spatial (ROIs)' , clique ))
    plot( V ),axis('tight'),%set(gca,'Xtick',EEG.params.freq.cfg.foi,'Xticklabel',EEG.params.freq.cfg.foi),xlabel('frequency (Hz)')
    
end

end    

end

%% Save the results of the stability analysis
mkdir(fullfile( params.basedir , params.stability.dir ),sprintf(params.results.subdir,patient))
clustfile = fullfile( params.basedir , params.stability.dir , ...
                    sprintf(params.results.subdir,patient) , ...
                    sprintf(params.clust.solfile , ...
                        params.decomposition.L1pen , ...
                        params.data.spectype , ...
                        params.data.boldtype , ...
                        params.decomposition.R) );
    
save(   clustfile , ...
        'similarity' , ...
        'solinfo' ,...
        'params' , ...
        'ccd' )
    
end
