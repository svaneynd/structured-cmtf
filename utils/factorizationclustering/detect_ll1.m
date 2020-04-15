% DETECT_LL1 takes a CPD factorization of any higher order tensor as input 
% and categorizes each term according to the following possibilities:
% a) true rank-1 term
% b) belonging to a (L,L,...,L,1,1,...,1)-term
% c) rank-1 term that, together with another rank-1 term, cancels out or
%       forms a degenerate pair
%
% A greedy search is used to cluster the initial rank-1 terms in different
% (L,L,...,L,1,1,...,1)-terms (see below in the code). Clustering can be detected
% in more than 1 mode.
%
% Remark: In the remainder of the code, the notation 'll1' is used for
% simplicity, but without loss of generality.
%
% INPUT
% - factors = cell array of factor matrices, e.g. as output of the function cpd
% - ll1thresh = threshold above which the correlation between 2 signatures
%               from different terms is considered to be significant (indicating 
%               at least LL1 structure, and no longer simple rank-1 structure)
% - degthresh = threshold above which the correlation between 2 signatures
%               from different terms is considered to be highly significant 
%               (indicating a degenerate/canceling pair, if a significant correlation
%               is found in all modes, and if the signatures' signs differ in an 
%               odd number of modes)
%
% OUTPUT
% - cpd_terms = 1xP array of indices of true rank-1 terms
% - ll1_terms = Qx2 cell array, in which element (q,1) is an array of
%               indices of terms belonging to the q-th detected LL1-term, 
%               and element (q,2) is an array of indices of the mode(s) in
%               which the terms are 'clustered'
% - deg_pairs = Sx2 array of degenerate/canceling pairs, in which the s-th row
%               contains the indices of terms of the s-th canceling pair
%
% Remark 2: If R is the rank of the CPD solution, the following holds: 
% numel(cpd_terms) + sum(cellfun(@(x)numel(x),ll1terms(:,1))) + numel(deg_pairs) = R,
% which means that every rank-1 term is assigned to a single category
%
% Author: Simon Van Eyndhoven
% Date: June 2017

function [ cpd_terms , ll1_terms , deg_pairs , norms ] = detect_ll1( factors , ll1thresh , degthresh , normthresh )

%% Initialize
if nargin < 2
    ll1thresh = 1%0.95;
    degthresh = 0.97;
    normthresh = 0.90;
end
assert(ll1thresh>0)
assert(degthresh>0)
assert(iscell(factors))

nmodes = numel(factors);
R = size(factors{1},2);
for m = 2:nmodes
    assert(size(factors{m},2)==R)
end

%% compute the norm of every component
norms = ones(1,R);
for i = 1:R
    for m = 1:nmodes
        norms(i) = norms(i)*frob(factors{m}(:,i));
    end
end
[~,comps] = sort(norms,'descend'); % sort the component according to variance, because degenerate pairs should really be detected, and normally are outliers in variance

%% Compute similarity (correlation) of signatures in all modes
simil = cell(1,nmodes);
for m = 1:nmodes
    simil{m} = cossimil(factors{m},factors{m});
end

%% determine, for every component r = 1..R, whether it is a separate
% CPD-term, or if it belongs to an LL1-term or potentially a canceling pair
cpd_terms = [];
ll1_terms = {};
deg_pairs = [];
nll1 = 0;
while ~isempty(comps) % go over all components    
    r = comps(1);
    if numel(comps) == 1 % only 1 component left: automatic CPD assignment
        % assign r to the cpd list
        cpd_terms = [ cpd_terms , r ];
        % remove r from the component list
        comps(1) = [];
    else % more than 1 component left
        pairs_simil = zeros(nmodes,numel(comps)-1); % similarity between r and the remaining components, in all modes
        for m = 1:nmodes
            pairs_simil(m,:) = simil{m}(r,comps(2:end));
        end
        pairs_simil_tll1 = (abs(pairs_simil)>ll1thresh);
        pairs_simil_tdeg = (abs(pairs_simil)>degthresh);
        pairs_simil_sign = prod(sign(pairs_simil),1);
        % rationale: 
        % - CPD structure requires: no similarity in any mode
        % - LL1 structure requires: significant similarity in minimally 1,
        %   and maximally nmodes-1 modes. 
        % - degeneracy requires: significant similarity in all modes, and
        %   overall sign difference
        if isempty(find(pairs_simil_tll1)) % no similarity: CPD
            % assign r to the cpd list
            cpd_terms = [ cpd_terms , r ];
            % remove r from the component list
            comps(1) = [];
        else % there is significant similarity (i.e., at least LL1 similarity) with one or more components
            % check for degeneracy first?
            degidx = find(sum(pairs_simil_tdeg,1)==nmodes); % indices of potential degenerate partner terms
            degidx1 = degidx(find(pairs_simil_sign(degidx)==-1)); % find the candidate terms for which the sign is negative
            degnorms = norms(comps(degidx1+1));
            [degnorms,degorder] = sort(degnorms,'descend'); 
            degidx = degidx1;
            if ~isempty(degidx) & (min(degnorms(1),norms(r))>=normthresh*max(degnorms(1),norms(r)))
                % assign r and the found component to the list of degenerate pairs
                deg_pairs = [ deg_pairs ; r comps(degidx(degorder(1))+1) ];
                % remove r and the found component from the component list
                comps([1,degidx(1)+1]) = [];
            else
                nll1 = nll1 + 1; % additional LL1 term is present, and will now be characterized
                % check for LL1 structure alternatively?
                % greedy search: find the greatest possible subset of
                % components that can be clustered in an LL1-term (other
                % option: find the subset of components that correspond in 
                % as many modes as possible, but this might detect spurious
                % similarities in multiple modes, which might then neglect true
                % corresponding components that belong in LL1 term, but only correspond in fewer modes)
                simil_count = sum(pairs_simil_tll1,2);
                ll1modes = find(simil_count==max(simil_count)); % look for the modes in which the greatest possible number of components can be grouped in LL1 terms
                if numel(ll1modes) > 1
                    ll1comps = [];
                    for m = 1:numel(ll1modes)
                        ll1comps(m,:) = find(pairs_simil_tll1(ll1modes(m),:)); % list all the components that match in mode m
                    end
                    % correlation, like the old implementation below, is
                    % not an appropriate metric, since the correlation will
                    % always be one if the number of components is 2, and
                    % if there is only one component, it will be NaN
%                     ll1compcorr = (corr(ll1comps')>0.999); % check in which modes the same components match (not exactly 1 due to rounding errors)
                    ll1compdist = (squareform(pdist(ll1comps,'cosine'))<1e-4);
                    [~,ll1modesidx] = max(sum(ll1compdist,2)); % the mode which agrees with most other modes
                    ll1modes = ll1modes(find(ll1compdist(ll1modesidx,:))); % find all agreeing modes
                    compsll1idx = ll1comps(ll1modesidx,:)+1;                
                else
                    compsll1idx = find(pairs_simil_tll1(ll1modes,:))+1; % list all the components that match in mode m
                end
                % assign r and the found components to the list of LL1 terms
                ll1_terms{nll1,1} = sort(comps([1,compsll1idx]));
                ll1_terms{nll1,2} = ll1modes;
                % remove r and the found components from the component list
                comps([1,compsll1idx]) = [];
            end
        end
    end
end

% HACK : do not consider LL1 terms (makes everything more difficult...)
for ll = 1:size(ll1_terms,1)
    cpd_terms = [ cpd_terms , ll1_terms{ll,1} ];
end

ll1_terms = [];

cpd_terms = sort(cpd_terms);

end