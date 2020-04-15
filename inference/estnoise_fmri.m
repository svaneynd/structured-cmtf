function [ recerr , regvar ] = estnoise_fmri( Y , sol , output , lnc , Lhrf )
% ESTNOISE_fMRI computes the residuals, and some second order statistics
% about them, for given sCMTF factors that fit a given fMRI data matrix Y.
%
% INPUTS
% - Y = (time bins x ROIs) matrix holding BOLD time series in different
% parcels of an atlas
% - sol = Tensorlab solution structure holding the sCMTF factors to approximate Y
% - output = Tensorlab output structure holding elementary information
% about the algorithm's convergence 
% - lnc = number of 'non-causal' samples in the HRFs
% - Lhrf = number of HRF samples to retain
% 
% OUTPUTS
% - recerr = vector holding the squared reconstruction/approximation error in each ROI
% - regvar = matrix holding, for each ROI, the variances of all
% 'regressors' (=temporal signatures in S)
%
% Author: Simon Van Eyndhoven (simon.vaneyndhoven@gmail.com)

if nargin < 4 | isempty(lnc)
    lnc = 0;
end
if nargin < 5 | isempty(Lhrf)
    Lhrf = 20;
end

[Is,Iv] = size(Y);
R = size(sol.factors.S,2);

h1 = sol . factors . HRFtoep_1(1:Lhrf);
h2 = sol . factors . HRFtoep_2(1:Lhrf);
h3 = sol . factors . HRFtoep_3(1:Lhrf);

truncateT = @(x)x(1+lnc:end,:);
convolveT = @(h,x)fftfilt(h,cat(1,x,zeros(lnc,size(x,2))));
ncconv = @(h,x)truncateT(convolveT(h,x));

% reconstruct the full approximation of the BOLD data
Yhat =  ncconv(h1,sol.factors.S) * (kr(sol.variables.spatial{1},sol.variables.spatial{4})) + ...
        ncconv(h2,sol.factors.S) * (kr(sol.variables.spatial{2},sol.variables.spatial{4})) + ...
        ncconv(h3,sol.factors.S) * (kr(sol.variables.spatial{3},sol.variables.spatial{4})) ;      
    
try % if the '_fmriuniq' model was used, with additional rank-1 terms for fMRI
    Yhat = Yhat + (sol.factors.N * sol.factors.P');
catch
end
 
Yres = Y - Yhat;        


recerr = zeros( 1 , Iv );
regvar = zeros( R , Iv );

% -- process all ROIs
for v = 1 : Iv
    
    % compute the (unnormalized) variance of the residual
    recerr(v) = sum(Yres(:,v).^2);
    
    hv = 0;
    for hidx = 1 : 3
        bf = sol.factors.(sprintf('HRFtoep_%d',hidx))(:,1); % basis function
        hv = hv + sol.variables.spatial{hidx}(v) * bf; % add weighted basis function to the local HRF
    end
    
    % construct a local 'design matrix'
    if lnc == 0
        Div = fftfilt(hv,sol.factors.S);
    else
        Div = ncconv(hv,sol.factors.S);
    end
    
    % compute the variance of the 'regressors' using standard result from
    % linear regression 
    % Note: to obtain the final variance, scaling with sqrt(recerr) is
    % still needed. This is performed in the later inference function.
    varbetav = ( 1 / (Is - R)) * inv( Div'*Div );    
    regvar( : , v ) = diag(varbetav);
    
end

tolerance = abs(sqrt(sum(recerr)) - output.relerr(2));
assert( tolerance < 1e-4 , sprintf('Correspondence within tolerance %2.2e',tolerance) )

end