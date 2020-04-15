function [ similarity , solinfo , csol ] = factorizationclustering( sol , modes , options )
% FACTORIZATION_CLUSTERING checks the stability of the estimated factors
% of a matrix-and/or-tensor factorization over runs of the optimization
% algorithm.
%
% INPUT
% - sol = cell array of solutions
% - modes = array of indices indicating which modes of the tensor should be
%           taken into account (one might choose to e.g. not consider LL modes)
% - options = structure with several hyperparameters (e.g. similarity criteria)
% 
% OUTPUT
% - similarity = structure holding the results of the clustering analysis
% - solinfo = meta-information on the solutions, i.e. info on the factor
% signatures, without the signatures themselves
% - csol = clustered solutions
%
% Author: Simon Van Eyndhoven (simon.vaneyndhoven@gmail.com)

%% 0. load data and initialize
    
    % define additional options
    if ~exist('options','var')
        options = struct;
        options.Xin = 'thresh';
        options.maxcost = 0.05;
        options.roundperm = 'avglink';
    end  
    
    % determine the tensor and factorization dimensionality    
    params = struct;
    params.decomposition.nfacs = numel(sol);
    params.decomposition.R = size(sol{1}{1},2);
    params.decomposition.nmodes = getorder(sol{1});
    
    % initialize
    params.thresh = zeros(1,params.decomposition.nmodes);
    params.TCthresh = -0.85;
    params.ll1thresh = 1; % similarity threshold above which 2 or more rank-1 terms are grouped into a block term (e.g. (L,L,1)-term)
    params.degthresh = 0.99; % similarity threshold above which 2 components are classified to be degenerate (i.e. nearly the same but opposite sign)
    params.normthresh = 0.90;
    
    if nargin < 2
        modes = 1:params.decomposition.nmodes;
    end

    if numel(params.thresh)==1
        params.thresh = repmat(params.thresh,1,params.decomposition.nmodes); % use the same threshold in all modes
    end

%% 1. sign correction
csol = sol;
for fac = 1:params.decomposition.nfacs
    csol{fac} = cpdsigncorrect(csol{fac},'energy');
end

%% 2. LL1 clustering: detect if some components form LL1 clusters. 
% If this is the case, no stability is expected for individual components.
% Rather, only stability of the whole 'LxL' matrix is expected.
% How? check if strong correlations in one mode exist (if it's in nmodes-1 modes:
% probably canceling/degenerate components)
% Also flag degenerate components.

solinfo = struct;
 
for fac = 1:params.decomposition.nfacs
    clear('*_temp')
    sol_temp = csol{fac};
    [cpd_terms_temp,ll1_terms_temp,deg_pairs_temp,norms_temp] = detect_ll1(sol_temp,params.ll1thresh,params.degthresh,params.normthresh);
    solinfo(fac).cpd_terms = cpd_terms_temp;
    solinfo(fac).ll1_terms = ll1_terms_temp;
    solinfo(fac).degenerateClus = deg_pairs_temp;
    solinfo(fac).norms = norms_temp;
end

for fac = 1:params.decomposition.nfacs
    keepl = true(1,size(solinfo(fac).ll1_terms,1)); % keep track of which LL1-terms can be taken off the list
    for l = 1:size(solinfo(fac).ll1_terms,1) % for all detected LL1-terms
        % a component is okay if it has strictly more than 2 'L-modes'
        % (cfr. CPD uniqueness)
        if params.decomposition.nmodes - numel(solinfo(fac).ll1_terms{l,2}) > 2 % L-modes form a CPD slice
            solinfo(fac).cpd_terms = [ solinfo(fac).cpd_terms , solinfo(fac).ll1_terms{l,1} ];
            solinfo(fac).cpd_terms = sort(solinfo(fac).cpd_terms,'ascend');
            keepl(l) = false;
        end
    end
    % only keep the real LL1-terms, i.e. those having a true ambiguity in
    % some modes
    solinfo(fac).ll1_terms = solinfo(fac).ll1_terms(keepl,:); 
end

%% 3. stability detection
% -- cluster the components over runs
similarity = detect_clusters(csol,solinfo,modes,options);

end