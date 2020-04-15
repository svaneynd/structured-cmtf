function [ V ] = wholebrainfromrois( v , Vatlas ) 
% WHOLEBRAINFROMROIS takes as input a vector with every entry corresponding 
% to a BOLD intensity in one region of interest (ROI), and expands it again 
% into a whole brain volume with voxels at a finer resolution.
%
% INPUT
% - v = vector with activation for all ROIs 
% - Vatlas = 3D volume with an ROI label at every voxel
% 
% OUTPUT
% - V = 3D volume consisting of voxels at the resolution of the original
% atlas file
% 
% Author: Simon Van Eyndhoven (simon.vaneyndhoven@gmail.com)

%% Check input arguments
rois = unique(Vatlas(:));
rois = rois(rois>0);

%% Expand the ROI activation to full brain
% -- create an empty image
V = zeros(size(Vatlas));

% -- write an associated value of an ROI to all the voxels that constitute the ROI
for roiidx = 1:length(rois)
    voxidx = (Vatlas==rois(roiidx));
    V(voxidx) = v(roiidx);
%         assert(numel(find(voxidx))>1) % ROIs should have a size greater than 1 voxel
end

%% Check functionality:

% sz = [79 95 68];
% nslices = 15;
% nrows = 3;
% ncols = ceil(nslices/nrows);
% 
% zc = round(linspace(1,sz(3),nslices+2));
% zc = zc(2:end-1); 
% 
% idx = [ 1 : 30 ]; % ROIs to 'activate'
%
% v = zeros(120,1);
% v(idx) = 0.5 + 0.5*linspace(0,1,length(idx));
% 
% [ V ] = deROIfy( v , Vatlas );
% 
% figure
% for nf = 1 : nslices
%     subplot(nrows,ncols,nf)
%     image( 64 * V(:,:,zc(nf))); % scale factor 64 for standard color map
% end

end