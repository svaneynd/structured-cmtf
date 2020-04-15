function [ sold , H , scaleB , Hbasis ] = disambiguatefmri( sol , scalenorm , signnorm , normp )
% DISAMBIGUATEFMRI takes as input a solution of the structured CMTF
% model with multiple HRF components and Khatri-Rao-structured spatial
% loading, and standardizes all (spatial) fMRI factors in order to remove 
% leftover sign and scale ambiguities.
%
% INPUTS
% - sol = SDF solution structure, with all variables and factors of the
%           structured CMTF model
% - scalenorm = scale normalization mode, i.e. based on what criterion to
%               adjust the scale of the HRFs (see below for options)
% - signnorm = sign normalization mode, i.e. based on what criterion to
%               adjust the sign of the HRFs (see below for options)
% - normp = p-parameter in p-norm computation (default: 2)
% 
% OUTPUTS
% - sold = disambiguated solution structure
% - H = matrix of time points x ROIs or voxels, holding all the adjusted HRFs
% - scaleB = the scaling factors that were applied to the spatial factors
% - Hbasis = the HRF basis functions that were used
%
% Author: Simon Van Eyndhoven (simon.vaneyndhoven@gmail.com)

%% Check input arguments
% -- ensure a valid scale normalization
normmodes = {   'maxampl' , ... % normalize based on maximal amplitude of the full HRF
                'energy' , ... % normalize based on the energy of the full HRF
                'stddev' }; % normalize based on the standard deviation of the full HRF
if nargin < 2 | isempty(scalenorm)
    scalenorm = 'energy';
else
    assert( ~isempty(find(strcmp(normmodes,lower(scalenorm)))) )
end

% -- ensure a valid sign normalization
normmodes = {   'posen' , ... % a valid HRF has most energy as overshoot (positive)
                'posfirst' }; % a valid HRF has overshoot before undershoot                
if nargin < 3 | isempty(signnorm)
    signnorm = 'posfirst';
else
    assert( ~isempty(find(strcmp(normmodes,lower(signnorm)))) )
end

if nargin < 4 | isempty(normp)
    normp = 2;
else
    assert( normp > 0 )
end

%% 0. Generate correction vectors to keep track of counterscaling factors
K = length(sol.variables.spatial) - 1; % number of HRF basis functions
[R,Iv] = size(sol.variables.spatial{K+1}); % number of sources / number of ROIs
    
scaleB = ones(1,Iv); % scaling factor for HRF weighting coefficients
scaleV = ones(1,Iv); % scaling factor for sources' spatial loadings

%% 1. Fix ambiguities for the HRF weighting coefficients
% The local HRF (i.e. in every ROI or voxel) should be unit norm (according
% to the chosen scale normalization mode!) and either have most of its energy as
% positive deviation (i.e. most of the HRF overshoots the horizontal axis)
% or have the HRF overshoot earlier than the undershoot (according to the
% chosen sign normalization).

sold = sol; % disambiguated solution
H = zeros(size(sol.factors.HRFtoep_1,1),Iv);

for iv = 1 : Iv % loop over all ROIs / voxels 
    
    % -- reconstruct the local HRF
    hv = zeros(size(sol.factors.HRFtoep_1,1),1);
    for hidx = 1 : K
        bf = sol.factors.(sprintf('HRFtoep_%d',hidx))(:,1); % basis function
        hv = hv + sol.variables.spatial{hidx}(iv)*bf ; % add weighted basis function to the local HRF
    end
    
    % -- appropriately correct the sign (if needed) of the local HRF
    posen = double(hv>0) .* (abs(hv).^normp); % energy of the HRF above the X-axis
    negen = double(hv<0) .* (abs(hv).^normp); % energy of the HRF above the X-axis
    switch lower(signnorm)
        case 'posen'
            if sum(negen) > sum(posen) % most energy lies below the X-axis
                hv = -hv;
                scaleB(iv) = -scaleB(iv);
            else % most energy lies above the X-axis, as desired
            end
        case 'posfirst'
            % compute the 'center of gravity' time of the positive/negative energy
            time = (1:length(hv))';
            if sum(time.*negen)/sum(negen) < sum(time.*posen)/sum(posen) % undershoot comes first
                hv = -hv;
                scaleB(iv) = -scaleB(iv);
            else % overshoot comes first, as desired
            end
        otherwise
            error('No valid sign normalization.')
    end    
    
    % -- appropriately correct the scale of the local HRF
    switch lower(scalenorm)
        case 'maxampl'
            sc = max(hv);
        case 'energy' % use variance as a surrogate for energy
            sc = (sum(abs(hv).^normp)).^(1/normp); % used to be: sqrt(sum(hv.^2)), but for hemodynamic signals the squaring is likely not meaningful?
        case 'stddev'
            sc = std(hv);
        otherwise
            error('No valid scale normalization.')
    end
    scaleB(iv) = (1/sc) * scaleB(iv);
    H(:,iv) = (1/sc) * hv;
    % -- correct all the involved HRF weighting coefficients
    for hidx = 1 : K
        sold.variables.spatial{hidx}(iv) = scaleB(iv) * sol.variables.spatial{hidx}(iv);
    end
end

%% 2. Store all aggregated corrections in the sources' spatial loadings
scaleV = 1./scaleB;
sold.variables.spatial{K+1} = bsxfun( @times , sol.variables.spatial{K+1} , scaleV );

for hidx = 1 : K
    Beta1 = ( kr( sold.variables.spatial{hidx} , sold.variables.spatial{K+1} ) )';
    Beta2 = sold.factors.(sprintf('b%dV',hidx));
        assert(frob(Beta1-Beta2)/frob(Beta1)<1e-9)  % since factors Bi and S from sol.variables.spatial are counterscaled,
                                            % the resulting products should not change!
    Hbasis(:,hidx) = sol.factors.(sprintf('HRFtoep_%d',hidx))(:,1);
end

end