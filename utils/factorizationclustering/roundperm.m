function [ M , X , compidx ] = roundperm( V , ncomps , mode )
%ROUNDPERM takes as input the leading eigenvectors of a similarity or 
% matching matrix over components from several runs, and uses this
% information to cluster the components into groups. Ideally, groups
% respect the cycle-consistency constraint, in which case they are cliques.
% 
% INPUTS
% - V = matrix with leading eigenvectors of adjacency/similarity matrix 
%           over components in all runs
% - ncomps = vector indicating the number of components in all runs
% - mode = string to indicate which rounding method should be used (see
%           notes below)
%
% OUTPUTS
% - M = cell matrix containing all groups of stable components
% - X = approximated adjacency matrix
% - compidx = indices of the cliqued components of all runs
%
% MODES
% 1. 'greedy' (default): Take the first non-matched row of V and find an
% orthogonal matrix that turns it (approximately) into a row of the
% identity matrix. Then, group this row with (at most) one row from every
% other run, if a similarity metric is sufficiently high. This might return
% more groups than the initial size of V (e.g. if the similarity threshold
% is high)
% 2. 'singlelink': apply single-linkage clustering on the rows of V
% directly, cfr. Bajaj, Chandrajit, et al. "SMAC: simultaneous mapping and 
% clustering using spectral decompositions." International Conference on 
% Machine Learning. 2018. This method returns exactly 
% 3. 'avglink': apply average-linkage clustering on the rows of V
% directly, cfr. ICASSO
% 4. 'hung': Sample rows from V such that all of the distinct components
% are observed, and use this to rotate and 'anchor' the matrix such that a
% linear assignment problem can be solved in every block, using the
% Hungarian algorithm.

if ~exist('mode','var')
    mode = 'greedy';
end

nruns = length(ncomps);
ngroups = size(V,2);
    
compidx = cell(2,nruns);
for run = 1:nruns
    compidx{1,run} = 1:ncomps(run); % relative index
    compidx{2,run} = compidx{1,run} + sum(ncomps(1:run-1)); % absolute index
end    

%% Group the components using one of several algorithms
switch mode
    case 'greedy'        
    %% Greedy clustering
    Tmatch = 0.75; % threshold for matching    
    
    compidxtemp = compidx; % save for later

    W = V;
    
    % Create empty output variable that holds the assignment groups
    % M is cell array that indicates which component of every run is placed
    % in which group
    M = nan(nruns,0);

    refset = 1; % at any time, refset is the index of the reference set, i.e. the set of rows that other rows will be compared to
    setisempty = repmat(false,nruns,1);
    while ~all(setisempty)

        % take the first of the remaining rows of the reference set
    %     v1 = W(compidx,:);
        v1 = W(compidx{2,refset}(1),:);
        m = nan(nruns,1);
        m(refset) = compidx{1,refset}(1);

        % find a unitary matrix P such that P*v1' = e1'
        % taking the outer product of both sides yields P*v1'*v1*P^T = e1'*e1 = diag([1 0 ... 0])
        % this means that the matrix P diagonalizes v1'*v1, and can hence be
        % computed as the set of (transposed) eigenvectors of that matrix    
        [P,D] = eig(v1'*v1);
        [d,order] = sort(diag(D),'descend');    
    %         assert(abs(log(d(1)))<log(1.5)) % the (only non-zero) eigenvalue
    %         should be reasonably close to one for well-behaved models 
    % THIS IS NOT NECESSARILY TRUE, THIS IS ONLY THE CASE IF THE INPUTS HAVE
    % BEEN SCALED BY THE DEGREE IN THE GRAPH!
        P = P(:,order)';

        % transform the matrix W such that its first row approximately equals e1'
        W = W*P';
        v1 = W(compidx{2,refset}(1),:);

        compidx{1,refset}(1) = []; % remove this row from further analysis
        compidx{2,refset}(1) = []; % remove this row from further analysis

        % cluster rows from other sets    
        for run = refset+1 : nruns
            if ~setisempty(run)%~isempty(compidx{2,run})
            % compute the cosine similarity with all the rows that are left in
            % this set        
            simil = cossimil(W(compidx{2,run},:)',v1'); % nrows x 1 vector
            [maxsimil,order] = sort(simil,'descend');

            % if a row meets the criterion: store it as a match
            if maxsimil(1) >= Tmatch
                m(run) = compidx{1,run}(order(1));
                compidx{1,run}(order(1)) = []; % remove this row from further analysis
                compidx{2,run}(order(1)) = []; % remove this row from further analysis
            end

            % check if the set is now empty
            if isempty(compidx{2,run})
                setisempty(run) = true;            
            end

            end

        end

        % check if the reference set is now empty
        if isempty(compidx{2,refset})
            setisempty(refset) = true;  
            refset = find(~setisempty,1);
        end

        % add the clustered rows
        M = [ M , m ];

        % if only one set remains, all those components end up in a separate
        % component
        if numel(find(setisempty==true)) == nruns - 1 % only the last set remains
            run = find(~setisempty);
            for comp = compidx{1,run}
                m = nan(nruns,1); % create an empty assignment column
                m(run) = comp;
                M = [ M , m ];
            end
            setisempty(run) = true;
        end

    end

    compidx = compidxtemp;
        
    case 'singlelink'
    %% Single-linkage clustering    
    
    % Create empty output variable that holds the assignment groups
    % M is cell array that indicates which component of every run is placed
    % in which group
    M = nan(nruns,ngroups); 
    
    % Perform the clustering
    T = cluster(linkage(V,'single','euclidean'),'maxclust',ngroups);
    
    % translate the clustering info into the grouping table
    for group = unique(T)'
        groupcomps = find(T==group);
        for c = 1:numel(groupcomps)
            crun = find(cellfun(@(x)~isempty(find(x==groupcomps(c))),compidx(2,:)));
            M(crun,group) = compidx{1,crun}(find(compidx{2,crun}==groupcomps(c)));
        end
    end
    
    case 'avglink'
    %% Single-linkage clustering    
    
    % Create empty output variable that holds the assignment groups
    % M is cell array that indicates which component of every run is placed
    % in which group
    M = nan(nruns,ngroups); 
    
    % Perform the clustering
    T = cluster(linkage(V,'average','euclidean'),'maxclust',ngroups);
    
    % translate the clustering info into the grouping table
    for group = unique(T)'
        groupcomps = find(T==group);
        for c = 1:numel(groupcomps)
            crun = find(cellfun(@(x)~isempty(find(x==groupcomps(c))),compidx(2,:)));
            M(crun,group) = compidx{1,crun}(find(compidx{2,crun}==groupcomps(c)));
        end
    end
    
    case 'hung'
    %% linear assignment via Hungarian algorithm
    
    % Create empty output variable that holds the assignment groups
    % M is cell array that indicates which component of every run is placed
    % in which group
    M = nan(nruns,ngroups);
    
    % Perform single linkage clustering to find initial, rough clusters
    T = cluster(linkage(V,'average','euclidean'),'maxclust',ngroups);
    groups = cell(ngroups,1);
    for group = 1:numel(groups)
        groups{group} = find(T==group);
    end
    
    % sample at random ngroups rows from V, taken into account that they 
    % should all belong to separate groups
    nit = min(1e3,round(0.5*size(V,1))); % number of iterations
    bestrows = cellfun(@(x)x(randi(numel(x))),groups);
    objf = Inf(nit,2); % objective function  % alternative: ratio of smallest/largest eigenvalue
    objf(1,1) = sum(sum(abs(eye(ngroups)-V(bestrows,:)*V(bestrows,:)')))/(ngroups^2);
    eigvals = sort(eig(V(bestrows,:)*V(bestrows,:)'),'descend');
    objf(1,2) = abs(eigvals(1)) - abs(eigvals(end));
    objfopt = objf(1,1);
    for it = 2:nit
        % sample 1 element from every group, at random
        randrows = cellfun(@(x)x(randi(numel(x))),groups);
        Vd = V(randrows,:);
        % compute the objective function
        objf(it,1) = sum(sum(abs(eye(ngroups)-Vd*Vd')))/(ngroups^2);
        eigvals = sort(eig(Vd*Vd'),'descend');
        objf(it,2) = abs(eigvals(1)) - abs(eigvals(end));
        if objf(it,1) < objfopt
            objfopt = objf(it);
            bestrows = randrows;
        end
    end
    
    % using the most representative rows, transform the eigenvectors into
    % cost matrices for every run
    D = V*V(bestrows,:)'; % similarity matrix
    D = max(D(:)) - D; % normally, the maximal value of the similarity should be close to 1
    D = max(D,0); % all cost values should be positive
    D = mat2cell(D,ncomps,ngroups); % convert to a cell per run
    
    % using the Hungarian algorithm / linear assignment, assign the
    % components to groups
    for run = 1:numel(ncomps)
        assignment = assignmentoptimal(D{run}); % for every row of this block: to which group does it belong?
        rungroups = find(assignment>0); % all the rows that are assigned in this run - this should be all the rows of the block
            assert(numel(rungroups)==numel(assignment))
        M(run,assignment(rungroups)) = rungroups';
    end
        
end

% sort the groups according to the components of the first run
run1 = find(M(1,:)); % indices of the groups in which components of run 1 appear
[~,run1] = sort(M(1,run1),'ascend');
M = M(:,[run1,setdiff(1:size(M,2),run1)]);

% construct an adjacency matrix from the grouping info in M
X = adjacency( M , ncomps );

end

function A = adjacency( M , ncomps )  
    [nruns,ncliques] = size(M);
        assert(nruns==length(ncomps))
    A = cell(nruns,nruns);
    % initialize the adjacency matrix
    for s1 = 1:nruns
        for s2 = 1:s1-1
            A{s1,s2} = zeros(ncomps(s1),ncomps(s2));
            A{s2,s1} = A{s1,s2}';
        end
        A{s1,s1} = eye(ncomps(s1));
    end
    % fill in the cliques in the adjacency matrix
    for cliqueidx = 1:ncliques
        clique = M(:,cliqueidx);
        cliqueruns = find(~isnan(clique))'; % (for loop needs a row vector)
        if length(cliqueruns) > 1
            for run1 = cliqueruns
                col = clique(run1);
                for run2 = setdiff(cliqueruns,run1)
                    row = clique(run2);
                    A{run2,run1}(row,col) = 1;
%                     A{run1,run2}(col,row) = 1; % adjacency is symmetric
                end
            end
        end
    end
    % convert into a matrix
    A = cell2mat(A);
end

