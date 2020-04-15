function [ Xperm , permorder ] = waveletresample( X , dwtorder , adsync , permorder0 )
% WAVELETRESAMPLE resamples the data by shuffling the wavelet-domain 
% coefficients. First, the data X is transformed to the discrete wavelet 
% domain, where coefficients are normally quite decorrelated and potential 
% dependencies between samples are less severe than in the time domain. 
% Exchangeability of wavelet coefficients is then warranted, and coefficients 
% in every level can be resampled without replacement (i.e. permuted), 
% after which the data can be reconstructed in the time domain. 
% To preserve spatial correlations, the same permutation is performed 
% for all (spatial) variables in X.
%
% INPUTS 
% - X = data matrix (observations x variables)
% - dwtorder = maximal order of the wavelet decomposition (default: 5)
% - adsync = boolean to choose whether the approximation and the highest
% detail level coefficients should be reshuffled in the same way (default:
% false)
% - permorder0 = provided reshuffling order
% 
% OUTPUTS
% - Xperm = (wavelet-)resampled data matrix
% - permorder = reordering of the coefficients in every wavelet level
%
% Author: Simon Van Eyndhoven (simon.vaneyndhoven@gmail.com)

%% Check input arguments
if ~exist('dwtorder','var') | isempty(dwtorder)
    dwtorder = 5;
end
if ~exist('adsync','var') | isempty(adsync)
    adsync = false;
end

[nt,nv] = size(X);

Xperm = zeros(nt,nv);

wname = 'db4';

%% Compute a single wavelet transform and randomly pick a reshuffling order
% -- Pick the first voxel/ROI to compute a representative wavelet transform
% and obtain the sizes of all levels
v = 1;
[wcoeff,wl] = wavedec( X(:,v) , dwtorder , wname );

% -- Break the wavelet coefficients into a cell array
wcoeff = mat2cell(wcoeff,wl(1:end-1),1);

% -- Randomly determine a reshuffling for every level
if ~exist('permorder0','var') | isempty(permorder0)
permorder = cell(dwtorder+1,1);
for lev = 1 : dwtorder+1
    permorder{ lev } = randperm( wl(lev) );
end
else
    permorder = permorder0;
end

if adsync
    permorder{ 1 } = permorder{ 2 };
end

% -- Reshuffle the coefficients at each level
wcoeffperm = cellfun( @(x,y)x(y) , wcoeff , permorder , 'Uniformoutput',false );

% -- Reconstruct a time series from the reshuffled coefficients
Xperm(:,v) = waverec( cell2mat(wcoeffperm) , wl , wname );

%% Loop over spatial variables (voxels/ROIs)
for v = 2 : nv
    % -- compute the wavelet transform
    wcoeff = wavedec( X(:,v) , dwtorder , wname );
    wcoeff = mat2cell(wcoeff,wl(1:end-1),1);
    % -- reshuffle the coefficients using the previously chosen order
    wcoeffperm = cellfun( @(x,y)x(y) , wcoeff , permorder , 'Uniformoutput',false );
    % -- reconstruct the data in the time domain
    Xperm(:,v) = waverec( cell2mat(wcoeffperm) , wl , wname );
end

end