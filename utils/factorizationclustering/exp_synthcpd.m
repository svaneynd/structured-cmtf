%--------------------------------------------------------------------------
%
% Experiment to check the stability of a factorization of synthetic data.
%
%--------------------------------------------------------------------------

%% Initialization
clc;
clearvars;
close('all');

%% Parameters
recompute = false; % if this decomposition has already been performed, you can choose
                    % to load the existing factorizations

params = struct;

    % -- synthetic data
    params.synth = struct;
    params.synth.s_t = [10 50 20]; % tensor dimensions, e.g. subjects x time course x frequency
    params.synth.R_t = 8; % true rank of the data
    params.synth.snr_t = -5; % signal-to-noise ratio of the data

    % -- decomposition
    params.decomposition.R = params.synth.R_t; % number of rank-1 terms to estimate
    params.decomposition.Initialization = @cpd_rnd;
    params.decomposition.Algorithm = @cpd_als;
    params.decomposition.nit = 50; % number of times to repeat the CPD computation
    params.decomposition.TolFun = 1e-6;
    params.decomposition.TolX = 1e-10;
    params.decomposition.Display = 0;
    params.decomposition.MaxIter = 200;
    
%     % -- initialization
%     params.init.angle = pi/6; % 
    
    % -- clustering
    params.clust = struct;
    params.clust.Xin = 'thresh';
    params.clust.maxcost = 0.10;
    params.clust.roundperm = 'avglink';

%% Generate the data

% function to generate exponentially distributed variables
exppdf = @(x)gamrnd(1,2,x);

Ut = cell(1,3); % allocate variable to hold the data factor matrices

% -- generate the factor matrices
    % mode 1: binary
    Ut{1} = double(rand(params.synth.s_t(1),params.synth.R_t)>0.5);
    
    % mode 2: exponential
    Ut{2} = exppdf([params.synth.s_t(2),params.synth.R_t]);
    
    % mode 3: Gaussian
    Ut{3} = randn(params.synth.s_t(3),params.synth.R_t);

% -- construct tensor
T = cpdgen(Ut);
        
% -- add noise
T = noisy(T,params.synth.snr_t,@randn);


%% Compute the CPD of the noisy data tensor   
% -- compute the CPD a number of times, from random initializations of the factors
for it = 1 : params.decomposition.nit
    fprintf('\n%s\nIteration nr. %d\n%s\n',repmat('-',1,20),it,repmat('-',1,20))
    [sol{it},output{it}] = cpd(T,params.decomposition.R,params.decomposition);
        assert(strcmp(output{it}.Initialization.Name,'cpd_rnd')) % verify whether random initialization was indeed used
end


%% Graph-clustering of the repeated CPD solutions
[ similarity , solinfo , csol ] = factorizationclustering( sol , 1:getorder(T) , params.clust );

% similarity.M = runs x cliques matrix, (i,j)-th entry is the index of the CPD
%   component of the i-th run that belongs to the j-th clique
% 
% similarity.cliquegroups = cell array, k-th cell holds the clique indices
%   of the k-th superclique (only sufficiently large cliques are kept)
%
% similarity.cliqueruns = cell array, k-th cell holds the runs in which
%   member components of each clique appear

supercliques = 1 : length(similarity.cliquegroups);

for k = supercliques
    % -- make one figure per superclique
    figure;
    sgtitle(sprintf('superclique %d',k))
    sc_size = numel(similarity.cliquegroups{k});
    for cliqueidx = 1 : sc_size
        % -- aggregate the CPD signatures of this clique (CPD component)
        clique = similarity.cliquegroups{k}(cliqueidx); % absolute index of this clique (CPD component)
        runs = similarity.cliqueruns{k}{cliqueidx}; % in which runs does this component appear?
        nruns = length(runs);
        A = zeros(size(T,1),nruns);
        B = zeros(size(T,2),nruns);
        C = zeros(size(T,3),nruns);
        for runidx = 1 : nruns
            % find a next run in which the current component appears
            cliquerun = similarity.cliqueruns{k}{cliqueidx}(runidx);
            % determine which index the CPD component has within this run
            compidx = similarity.M(cliquerun,clique);
            % store the signatures of this component (in all modes)
            A(:,runidx) = csol{cliquerun}{1}(:,compidx);
            B(:,runidx) = csol{cliquerun}{2}(:,compidx);
            C(:,runidx) = csol{cliquerun}{3}(:,compidx);
        end
        
        % -- plot the signatures in all modes
        subplot(getorder(T)+1,sc_size,cliqueidx) % mode-1 signatures
        hold('on');
        plot(A,'color',[0.2 0.2 0.2])
        plot(mean(A,2),'color',[1 0 0],'linewidth',2)
        title(sprintf('clique %d',clique))
        subplot(getorder(T)+1,sc_size,sc_size+cliqueidx) % mode-2 signatures
        hold('on');
        plot(B,'color',[0.2 0.2 0.2])
        plot(mean(B,2),'color',[1 0 0],'linewidth',2)
        subplot(getorder(T)+1,sc_size,2*sc_size+cliqueidx) % mode-3 signatures 
        hold('on');
        plot(C,'color',[0.2 0.2 0.2])
        plot(mean(C,2),'color',[1 0 0],'linewidth',2)     
        subplot(getorder(T)+1,sc_size,3*sc_size+cliqueidx) % visualization of the occurrence over runs
        occurrence = zeros(1,params.decomposition.nit);
        occurrence(runs) = 1;
%         imagesc(occurrence); colormap('gray')
        bar(1:params.decomposition.nit,occurrence,1.0,'k')
        set(gca,'xtick',0:10:params.decomposition.nit,'ylim',[0,1.5])
        xlabel('runs')
    end
end