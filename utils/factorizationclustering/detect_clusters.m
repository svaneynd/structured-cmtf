% DETECT_CLUSTERS takes multiple CPD factorizations of the same higher 
% order tensor as input, and checks which terms are stable (i.e., are found
% in many of the factorizations).
%
% INPUT
% - sols = cell array of factorizations, each containing a cell array of 
%           factor matrices, e.g. as output of the function cpd
% - solsinfo = cell array of structures, which contain info on the nature
%               of the terms in the factorization (cpd-terms, ll1-terms, 
%               degenerate terms)
% - modes = indices of the modes in which similarity should be computed
%               not necessarily all modes should be considered, for example
%               if the different factorizations stem from 'split'
%               experiments along such mode
% - options = structure to specify additional options:
%               options.Xin = 'hung''/ 'hungN' / 'greedy' / 'round' : ways of binarizing the input matching matrix X
%               options.maxcost = value below which the match is validated  and retained (default: Inf)
%               options.roundperm = 'hung' / 'singlelink' : ways of rounding the eigenvectors into permutation matrices
%
% OUTPUT
% - similarity = structure holding the results of the clustering analysis

function [ similarity ] = detect_clusters( sols , solsinfo , modes , options )

%% Initialize
assert(iscell(sols))
assert(min(size(sols))==1)
assert(min(size(solsinfo))==1)
assert(isstruct(solsinfo));

% define additional options
if ~exist('options','var')
    options = struct;
    options.Xin = 'round';
    options.maxcost = Inf;
    options.roundperm = 'hung';
end

% determine the number of runs
nruns = numel(sols);
assert(nruns>1)
assert(numel(solsinfo)==nruns)

% determine the number of modes of the factorizations
nmodes = cellfun(@(x)numel(x),sols);
nmodes = unique(nmodes);
assert(numel(nmodes)==1); % all factorizations should have the same number of modes
assert(max(modes)<=nmodes)
modes = modes(:)';

% determine all possible pairs of runs
pairs = nchoosek(1:nruns,2);
npairs = size(pairs,1);

% create an array to store the similarity metric
simil = cell(npairs,nmodes);

% create an array to keep track of the number of components in every run
ncomps = cell(1,nruns);

%% Compare all pairs of runs
for p = 1:npairs
    
    fprintf('Processing pair of factorizations: %d out of %d\n',p,npairs)
    
    %% 0.1. pick two factorizations to compare
    solidx1 = pairs(p,1);
    solidx2 = pairs(p,2);
    sol1 = sols{solidx1};
    sol2 = sols{solidx2};
    
    %% 0.2. determine the 'true' number of terms in both factorizations:
    % discard degenerate/canceling terms and consider LL1-terms as one term only
    ncpd_1 = numel(solsinfo(solidx1).cpd_terms);
    nll1_1 = size(solsinfo(solidx1).ll1_terms,1);
    ncomps{solidx1} = [ncpd_1,nll1_1];
    ncpd_2 = numel(solsinfo(solidx2).cpd_terms);
    nll1_2 = size(solsinfo(solidx2).ll1_terms,1);
    ncomps{solidx2} = [ncpd_2,nll1_2];
    
    %% 0.3. initialize all similarity matrices
    for mode = modes
        if ~isempty(find(modes==mode)) % stability should also be checked in this mode
            % initialize a matrix of similarities between all components
            simil{p,mode} = NaN(ncpd_1+nll1_1,ncpd_2+nll1_2); 
        else % this mode is not of interest (e.g. if the different factorizations come from different splits in this mode)
            simil{p,mode} = [];
        end
    end
    
    %% 1.1. compare all cpd_1 terms to all cpd_2 terms
    if (ncpd_1>0) & (ncpd_2>0)
    for mode = 1:nmodes
        if ~isempty(find(modes==mode)) % stability should also be checked in this mode
            comps_1 = solsinfo(solidx1).cpd_terms;
            fac_1 = sol1{mode}(:,comps_1);
            comps_2 = solsinfo(solidx2).cpd_terms;        
            fac_2 = sol2{mode}(:,comps_2);
            simil{p,mode}(1:ncpd_1,1:ncpd_2) = cossimil(fac_1,fac_2);%corr(fac_1,fac_2);
        else % this mode is not of interest (e.g. if the different factorizations come from different splits in this mode)
        end
    end
    else % at least one of the factorizations has no true CPD terms
    end
    
    %% 1.2. compare all cpd_1 terms to all ll1_2 terms
    if ncpd_1 > 0
    for c2 = 1:nll1_2
        comps_1 = solsinfo(solidx1).cpd_terms;
        comps_2 = solsinfo(solidx2).ll1_terms{c2,1};
        cpdmodes = solsinfo(solidx2).ll1_terms{c2,2}; % modes in which terms could be clustered
        llmodes = 1:nmodes;
        llmodes(cpdmodes) = [];
        llmodes = intersect(llmodes,modes); % only modes of interest should be considered
        % treat all cpdmodes
        for mode = cpdmodes(:)'
            if ~isempty(find(modes==mode)) % stability should also be checked in this mode                
                fac_1 = sol1{mode}(:,comps_1);                        
                fac_2 = sol2{mode}(:,comps_2);
                simil_temp = abs(cossimil(fac_1,fac_2));%abs(corr(fac_1,fac_2));
                simil_temp = mean(simil_temp,2); % the vectors in fac_2 are clustered signatures in the LL1-terms            
                simil{p,mode}(1:ncpd_1,ncpd_2+c2) = simil_temp;
            else % this mode is not of interest (e.g. if the different factorizations come from different splits in this mode)
            end
        end
        % treat all llmodes
        for c1 = comps_1(:)'
            fac_1 = [];
            fac_2 = [];
            for m = 1:numel(llmodes)
                fac_1{m} = sol1{llmodes(m)}(:,c1);
                fac_2{m} = struct_normalize(sol2{llmodes(m)}(:,comps_2),[]);
            end
            fac_2{1} = fac_2{1}*diag(solsinfo(solidx2).norms(comps_2));
            if numel(llmodes) > 1 % normal case: at least 2 LL-modes
                cpd_slice = cpdgen(fac_1);
                ll_slice = cpdgen(fac_2);
                % compute the similarity between the CPD-like slice and the LL-slice
                % global scale difference possible! --> least-squares determination of scaling factor
                sc = dot(cpd_slice(:),ll_slice(:))/dot(cpd_slice(:),cpd_slice(:));
                cpd_slice = sc*cpd_slice;
                simil_temp = cossimil(cpd_slice(:),ll_slice(:));%corr(cpd_slice(:),ll_slice(:)); % other option: simil_temp = 1-frob(cpd_slice-ll_slice)/(max(frob(cpd_slice),frob(ll_slice)));
                for m = 1:numel(llmodes)
                    simil{p,llmodes(m)}(find(comps_1==c1),ncpd_2+c2) = simil_temp;
                end
            elseif numel(llmodes) == 1 % special case: only 1 LL-mode left
                [fac_2_est] = compare_subspaces(fac_2{1},fac_1{1});
                simil_temp = cossimil(fac_2{1}(:),fac_2_est(:));%corr(fac_2{1}(:),fac_2_est(:));
%                 simil_temp = abs(corr(fac_1{1},fac_2{1}));
%                 amp = sqrt(sum(fac_2{1}.^2));
%                 simil_temp = dot(amp,simil_temp)/sum(amp);
                simil{p,llmodes}(find(comps_1==c1),ncpd_2+c2) = simil_temp;
            else % all LL-modes are not of interest
            end
        end
    end
    else % first factorization has no true CPD terms
    end
    
    %% 2.1. compare all ll1_1 terms to all cpd_2 terms
    if ncpd_2 > 0
    for c1 = 1:nll1_1
        comps_1 = solsinfo(solidx1).ll1_terms{c1,1};
        comps_2 = solsinfo(solidx2).cpd_terms;
        cpdmodes = solsinfo(solidx1).ll1_terms{c1,2}; % modes in which terms could be clustered
        llmodes = 1:nmodes;
        llmodes(cpdmodes) = [];
        llmodes = intersect(llmodes,modes); % only modes of interest should be considered
        % treat all cpdmodes
        for mode = cpdmodes(:)'
            if ~isempty(find(modes==mode)) % stability should also be checked in this mode                
                fac_1 = sol1{mode}(:,comps_1);                        
                fac_2 = sol2{mode}(:,comps_2);
                simil_temp = abs(cossimil(fac_1,fac_2));%abs(corr(fac_1,fac_2));
                simil_temp = mean(simil_temp,1); % the vectors in fac_1 are clustered signatures in the LL1-terms            
                simil{p,mode}(ncpd_1+c1,1:ncpd_2) = simil_temp';
            else % this mode is not of interest (e.g. if the different factorizations come from different splits in this mode)
            end
        end
        % treat all llmodes
        for c2 = comps_2(:)'
            fac_1 = [];
            fac_2 = [];
            for m = 1:numel(llmodes)
                fac_1{m} = struct_normalize(sol1{llmodes(m)}(:,comps_1),[]);
                fac_2{m} = sol2{llmodes(m)}(:,c2);
            end
            fac_1{1} = fac_1{1}*diag(solsinfo(solidx1).norms(comps_1));
            if numel(llmodes) > 1 % normal case: at least 2 LL-modes
                cpd_slice = cpdgen(fac_2);
                ll_slice = cpdgen(fac_1);
                % compute the similarity between the CPD-like slice and the LL-slice
                % global scale difference possible! --> least-squares determination of scaling factor
                sc = dot(cpd_slice(:),ll_slice(:))/dot(cpd_slice(:),cpd_slice(:));
                cpd_slice = sc*cpd_slice;
                simil_temp = cossimil(cpd_slice(:),ll_slice(:));%corr(cpd_slice(:),ll_slice(:)); % other option: simil_temp = 1-frob(cpd_slice-ll_slice)/(max(frob(cpd_slice),frob(ll_slice)));
                for m = 1:numel(llmodes)
                    simil{p,llmodes(m)}(ncpd_1+c1,find(comps_2==c2)) = simil_temp;
                end
            elseif numel(llmodes) == 1 % special case: only 1 LL-mode left
                [fac_1_est] = compare_subspaces(fac_1{1},fac_2{1});
                simil_temp = cossimil(fac_1{1}(:),fac_1_est(:));%corr(fac_1{1}(:),fac_1_est(:));
%             simil_temp = abs(corr(fac_2{1},fac_1{1}));
%             amp = sqrt(sum(fac_1{1}.^2));
%             simil_temp = dot(amp,simil_temp)/sum(amp);
                simil{p,llmodes}(ncpd_1+c1,find(comps_2==c2)) = simil_temp;
            else % all LL-modes are not of interest
            end
        end
    end   
    else % second factorization has no true CPD terms
    end
    
    %% 2.2. compare all ll1_1 terms to all ll1_2 terms
    for c1 = 1:nll1_1
        for c2 = 1:nll1_2
            comps_1 = solsinfo(solidx1).ll1_terms{c1,1};
            comps_2 = solsinfo(solidx2).ll1_terms{c2,1};
            cpdmodes = intersect(solsinfo(solidx1).ll1_terms{c1,2},...
                        solsinfo(solidx2).ll1_terms{c2,2}); % modes in which terms could be clustered
            llmodes = 1:nmodes;
            llmodes(cpdmodes) = [];
            llmodes = intersect(llmodes,modes); % only modes of interest should be considered
            % treat all cpdmodes
            for mode = cpdmodes(:)'
                if ~isempty(find(modes==mode)) % stability should also be checked in this mode                
                    fac_1 = sol1{mode}(:,comps_1);                        
                    fac_2 = sol2{mode}(:,comps_2);
                    simil_temp = abs(cossimil(fac_1,fac_2));%abs(corr(fac_1,fac_2));
                    simil_temp = mean(simil_temp(:)); % the vectors in fac_1 and fac_2 are clustered signatures in the LL1-terms            
                    simil{p,mode}(ncpd_1+c1,ncpd_2+c2) = simil_temp;
                else % this mode is not of interest (e.g. if the different factorizations come from different splits in this mode)
                end
            end
            % treat all llmodes
            fac_1 = [];
            fac_2 = [];
            for m = 1:numel(llmodes)
                fac_1{m} = struct_normalize(sol1{llmodes(m)}(:,comps_1),[]);
                fac_2{m} = struct_normalize(sol2{llmodes(m)}(:,comps_2),[]);
            end
            fac_1{1} = fac_1{1}*diag(solsinfo(solidx1).norms(comps_1));
            fac_2{1} = fac_2{1}*diag(solsinfo(solidx2).norms(comps_2));
            if numel(llmodes) > 1 % normal case: at least 2 LL-modes
                ll_slice1 = cpdgen(fac_1);
                ll_slice2 = cpdgen(fac_2);
                % compute the similarity between the two LL-slices
                % global scale difference possible! --> least-squares determination of scaling factor
                sc = dot(ll_slice1(:),ll_slice2(:))/dot(ll_slice1(:),ll_slice1(:));
                ll_slice1 = sc*ll_slice1;
                simil_temp = cossimil(ll_slice1(:),ll_slice2(:)); %corr(ll_slice1(:),ll_slice2(:)); % other option: simil_temp = 1-frob(ll_slice1-ll_slice2)/(max(frob(ll_slice1),frob(ll_slice2)));
                for m = 1:numel(llmodes)
                    simil{p,llmodes(m)}(ncpd_1+c1,ncpd_2+c2) = simil_temp;
                end
            elseif numel(llmodes) == 1 % special case: only 1 LL-mode left
                if numel(comps_1) > numel(comps_2)
                    [fac_1_est] = compare_subspaces(fac_1{1},fac_2{1});
                    simil_temp = cossimil(fac_1{1}(:),fac_1_est(:));%corr(fac_1{1}(:),fac_1_est(:));
                else
                    [fac_2_est] = compare_subspaces(fac_2{1},fac_1{1});
                    simil_temp = cossimil(fac_2{1}(:),fac_2_est(:));%corr(fac_2{1}(:),fac_2_est(:));
                end
%                 simil_temp = abs(corr(fac_1{1},fac_2{1}));
%                 amp_1 = sqrt(sum(fac_1{1}.^2));
%                 amp_2 = sqrt(sum(fac_2{1}.^2));
%                 simil_temp = (amp_1'*amp_2).*simil_temp;
                simil{p,llmodes}(ncpd_1+c1,ncpd_2+c2) = simil_temp;
            else % all LL-modes are not of interest
            end
        end
    end
    
    %% 3. clear remaining temporary variables
    clear('*temp','*1','*2','*slice','fac*')
    
end

%% Determine the stability by looking at the reproducibility of components in different runs

    %% 0. assert that similarity was correctly computed
    % For modes that should not have been taken into account, the
    % similarity matrices should be zero. Note that this is different from
    % the initialized similarity values (NaN) in other modes.
    for m = 1:nmodes
        if find(modes==m)
            assert(all(cellfun(@(x)~isempty(x),simil(:,m))))
        else
            assert(all(cellfun(@(x)isempty(x),simil(:,m))))
        end
    end
    nnans = cell2mat(cellfun(@(x)numel(find(isnan(x))),simil,'UniformOutput',false));
    assert(sum(nnans(:))==0)
    
    %% 1. Construct a similarity graph of all components
    ncomps_tot = sum(cellfun(@(x)sum(x),ncomps));
    idx_start = downsample(cumsum(cell2mat(ncomps)),2,1);
    idx_start = [1,idx_start+1];
    Simil_tot = eye(ncomps_tot,ncomps_tot);
    Matched = eye(numel(ncomps),numel(ncomps));
    % calculate a large (n1+n2+...) x (n1+n2+...) matrix of the similarity
    % between all components over runs (averaged over modes)  
    Simil_totbin = Simil_tot; % binarized similarity graph
    method = options.Xin;
    for p = 1:npairs
        % calculate the similarity over modes for this pair of runs
        simil_p = 1;
        for m = modes
            simil_p = simil_p.*simil{p,m};%abs(simil{p,m});
        end
        simil_p = sign(simil_p).*abs(simil_p).^(1/numel(modes)); % geometric mean of the similarity over modes
        % store the similarity in the big matrix
        solidx1 = pairs(p,1);
        solidx2 = pairs(p,2);
        assert(solidx1<solidx2)
        Simil_tot(idx_start(solidx1):idx_start(solidx1+1)-1,idx_start(solidx2):idx_start(solidx2+1)-1) = simil_p;
        Simil_tot(idx_start(solidx2):idx_start(solidx2+1)-1,idx_start(solidx1):idx_start(solidx1+1)-1) = simil_p'; % symmetrical matrix
        
        switch method
            case 'hung'            
                % perform linear assignment to obtain a (pseudo-)permutation matrix in every block
                [M,cost] = assignmentoptimal(abs(1-abs(simil_p))); % matches
                midx = find(M);
                ass_p = zeros(size(simil_p));
                if cost < options.maxcost*length(midx)
                    for idx = midx'
                        ass_p(idx,M(idx)) = 1;
                    end
                    Matched(solidx1,solidx2) = 1;
                    Matched(solidx2,solidx1) = 1;
                else
                    Matched(solidx1,solidx2) = 0;
                    Matched(solidx2,solidx1) = 0;
                end
                
            case 'hungN' % perform Hungarian algorithm first in every mode, and aggregate then over modes (computationally expensive)
            case 'greedy' % use a similar greedy approach as in cpderr
                ass_p = zeros(size(simil_p));
                Ctemp = simil_p;
                if size(Ctemp,1) > size(Ctemp,2)
                    Ctemp = Ctemp';
                    ct = true;
                else
                    ct = false;
                end
                Ptemp = zeros(size(Ctemp));
                cost = 0;
                for r = 1:size(Ctemp,1)
                    [Cr,i] = max(Ctemp,[],1);
                    [maxCr,j] = max(Cr);
                    cost = cost + (1-maxCr);
                    Ptemp(i(j),j) = 1;
                    Ctemp(i(j),:) = 0;
                    Ctemp(:,j) = 0;
                end
                if cost < options.maxcost*size(Ctemp,1)
                    if ct
                        ass_p = Ptemp';
                    else
                        ass_p = Ptemp;
                    end
                    Matched(solidx1,solidx2) = 1;
                    Matched(solidx2,solidx1) = 1;
                else
                    Matched(solidx1,solidx2) = 0;
                    Matched(solidx2,solidx1) = 0;
                end
            case 'thresh'
                % use simple rounding (to 0 and 1) to binarize the matching matrix
%                 ass_p = round(max(simil_p,0));
                ass_p = double(simil_p > (1-options.maxcost));
                Matched(solidx1,solidx2) = 1;
                Matched(solidx2,solidx1) = 1;    
        end
        Simil_totbin(idx_start(solidx1):idx_start(solidx1+1)-1,idx_start(solidx2):idx_start(solidx2+1)-1) = ass_p;
        Simil_totbin(idx_start(solidx2):idx_start(solidx2+1)-1,idx_start(solidx1):idx_start(solidx1+1)-1) = ass_p'; % symmetrical matrix 
    end
    
    bl_cliqueruns = {}; % baseline groups (before map sync) of runs that belong together
    while ~all(Matched(:)==0)
        [~,row] = max(sum(Matched,2)); % index of run with the most matches
        m = find(Matched(row,:)); % find which other runs are matched with this run
        bl_cliqueruns = [ bl_cliqueruns , {m} ];
        for s = m
            Matched(s,:) = 0;
            Matched(:,s) = 0;
        end        
    end


    %% 3. Find clusters by low-rank approximation (using eigenvalue decomposition)
    [P,D] =  eig(Simil_totbin); % figure,imagesc(Simil_tot), figure,imagesc(Simil_totbin) , colormap('gray')

    [lambda,order] = sort(diag(D),'descend');
    
    % -- estimate the number of significant eigenvalues
    ncomps = cellfun(@(x)sum(x),ncomps); % count all the CPD and LL1-terms in every run
    Rmin = max(ncomps);
    Rest = universesize(lambda,'norm',Rmin);
    
    % -- ask the user to specify (adjust) the number of eigenvalues to
    % retain
    idx_tick = cumsum(ncomps);
    figure,plot(lambda,'k','linewidth',4) % linewidth 4 for exporting figures?
    xlabel('rank-1 terms of all runs','fontsize',28)
    title('eigenvalue spectrum of similarity matrix','fontsize',28)
    hold('on'),scatter(idx_tick,lambda(idx_tick),200,'r','filled')
    idx_tick = idx_tick([1 5:5:end]);
    set(gca,'Xtick',idx_tick,'fontsize',20)
    
    pause;
    
    prompt = {sprintf('Enter a new value for R (minimal value: %d)\nEstimated value: R = %d',Rmin,Rest)};
    prompttitle = 'Choose the number of significant eigenvalues';
    promptdims = [1 70];
    definput = {sprintf('%d',Rest)};
    newr = inputdlg(prompt,prompttitle,promptdims,definput);
    
    close(gcf)
    
    R = str2num(newr{1});    
    V = bsxfun(@times,P(:,order(1:R)),sqrt(max(lambda(1:R),0))');
    
    [ M , X , compidx ] = roundperm( V , ncomps , options.roundperm );
    [ M , cliquegroups , cliqueruns ] = sortcliquesauto(M,round(0.05*nruns)); % default: components that appear in less than 5% of runs are no longer considered
    
    %% 4. Find clusters as in ICASSO
    sR = struct;
    sR.cluster.similarity = max(0,Simil_tot);
    sR.cluster.strategy = 'AL';
%     [sR.cluster.partition,sR.cluster.dendrogram.Z,sR.cluster.dendrogram.order]=...
%                                     hcluster(abs( 1 - sR.cluster.similarity ),sR.cluster.strategy);
%                                 
%     % correct the representation of the clustering
%     maxR = max(ncomps);
%     sR.M = nan(length(ncomps),maxR);
%     
%         assert(max(sR.cluster.partition(maxR,:)==maxR))
%     
%     for run = 1:length(ncomps)
%         relidx = 1:ncomps(run);
%         absidx = sum(ncomps(1:run-1)) + relidx;
%         sR.M(run,sR.cluster.partition(maxR,absidx)) = relidx;
%     end    

    % -- create output variables
    similarity = struct;
    similarity.cont = Simil_tot;
    similarity.bin = Simil_totbin;
    similarity.optbin = X;
    similarity.ncomps = ncomps;
    similarity.cliquegroups = cliquegroups;
    similarity.cliqueruns = cliqueruns;
    similarity.bl_cliqueruns = bl_cliqueruns;
    similarity.icasso = sR;
    similarity.M = M;

end