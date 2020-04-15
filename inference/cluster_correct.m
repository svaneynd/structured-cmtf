function [ Vo , frac ] = cluster_correct( Vi , Smin )
% INPUTS
% - Vi = binary image (mask) with (in)active voxels
% - Smin = threshold on the cluster extent
% 
% OUTPUTS
% - Vo = output binary image (mask) in which only voxels that appear in
%           clusters that meet the threshold are labeled 'active'
% - frac = percentage of active voxels of the original image that is kept

% -- Find the coordinates of all voxels
[i,j,k] = ind2sub(size(Vi),find(Vi)); % i,j,k = voxel coordinates in x,y,z directions

% -- Detect which voxels appear in which clusters
clustlabel = spm_clusters([i,j,k]');

% -- Count the number of voxels per cluster
Nv = histc(clustlabel,(0:max(clustlabel))+0.5);
Nv = Nv(1:end-1);

% -- Perform cluster extent correction
Vo = zeros(size(Vi));
validclust = find( Nv >= Smin ); % remove clusters that are too small

% loop over all clusters that meet the criterion
for cidx = validclust
    % determine which voxels are in this cluster
    vidx = find( clustlabel == cidx );
    % set those voxels to 'active' in the output image
    imidx = sub2ind(size(Vi),i(vidx),j(vidx),k(vidx));
    Vo( imidx ) = Vi( imidx );
end

% compute which fraction of original voxels is kept
frac = numel(find(Vo)) / numel(find(Vi));
assert( numel(find((Vo~=0)&(Vi==0))) == 0 ) % no new voxels should have been activated
assert( numel(find(Vo)) == sum(Nv(validclust)) )

end