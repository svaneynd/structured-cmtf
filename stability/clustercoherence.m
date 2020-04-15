function [ C ] = clustercoherence( S , comps )
% S = similarity matrix between components
% comps = indices of components that belong to one cluster

assert(max(comps)<=size(S,2))

nincl = length(comps);
nexcl = size(S,2) - length(comps);

partition = zeros(1,nincl+nexcl);
partition(comps) = 1;

withinsimil = S(partition==1,partition==1);
betweensimil = S(partition==1,partition==0);

    assert(numel(withinsimil)==nincl^2)
    assert(numel(betweensimil)==nincl*nexcl)

C = 1/(nincl^2)*sum(withinsimil(:)) - 1/(nincl*nexcl)*sum(betweensimil(:));

end