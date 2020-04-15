function [ Xperm , deltaphi ] = fourierresample( X )
% FOURIERRESAMPLE resamples the data by shuffling the frequency-domain 
% coefficients. First, the data X is Fourier-transformed to the frequency
% domain, where coefficients are normally quite decorrelated and potential 
% dependencies between samples are less severe than in the time domain. 
% Exchangeability of Fourier coefficients is then warranted, and coefficients 
% in every level can be resampled without replacement (i.e. permuted), 
% after which the data can be reconstructed in the time domain. 
% To preserve spatial correlations, the same permutation is performed 
% for all (spatial) variables in X.
%
% INPUTS 
% - X = data matrix (observations x variables)
% 
% OUTPUTS
% - Xperm = (wavelet-)resampled data matrix
% - deltaphi = random deviations that were added to the phases of the
% original data
%
% Author: Simon Van Eyndhoven (simon.vaneyndhoven@gmail.com)

%% Check input arguments
[nt,nv] = size(X);

Xperm = zeros(nt,nv);

nfft = 2*ceil(nt/2);

%% Create random (but symmetric) phase differences
deltaphi = zeros( nfft , 1 );
deltaphi( 1 : nfft/2+1 ) = 2*pi * rand( nfft/2+1 , 1 );
deltaphi( nfft/2+2 : nfft ) = -flipud( deltaphi( 2 : nfft/2 ) );

%% Permute data by adding spatially synchronous phase noise to the Fourier coefficients
% -- Compute the Fourier transform of the data
Xf = fft( X , nfft );

% -- Distill the phase info
phiX = angle(Xf);

% -- Add random deviations to the phase
phiX = bsxfun( @plus , phiX , deltaphi );
Xf = abs(Xf) .* exp(i*phiX) ;

% -- Reconstruct the data
Xperm = ifft( Xf );
Xperm = real(Xperm(1:nt,:));

end