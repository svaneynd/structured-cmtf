function [ pval , Hval , kdepar , pdf ] = kde_entropy( X , ng , verbose , bw )
% KDE_ENTROPY_FIXCENTERS computes a multivariate PDF for the data based on
% kernel density estimation (KDE). Afterwards, this function 
% computes the probability (p-value) for each point (which leads to its 
% entropy and p-value).
%
% INPUT
% - X = input data (observations x variables)
% - ng = number of grid points in each dimension
%
% OUTPUT
% - pval = probability for each input data point, computed in leave-one-out fashion
% - Hval = entropy for each input data point, computed in leave-one-out fashion
%
% Reference:
%  Kernel density estimation via diffusion
%  Z. I. Botev, J. F. Grotowski, and D. P. Kroese (2010)
%  Annals of Statistics, Volume 38, Number 5, pages 2916-2957.
%  https://nl.mathworks.com/matlabcentral/fileexchange/58312-kernel-density-estimator-for-high-dimensions?s_tid=prof_contriblnk

%% Initialize
try
addpath('C:\Matlab toolboxes\akde')
catch
end

[n,d] = size(X);
    
nk = n; % number of kernels = number of observations

if nargin < 2
    ng = 100;
end

if nargin < 3 | isempty(verbose)
    verbose = false;
end

if nargin < 4 | isempty(bw)
    silverman = true;
else
    silverman = false;
end

if d==3
%% Create an evalution grid
% -- find the bounding box dimensions
MAX = max(X,[],1); 
MIN = min(X,[],1); 
scaling = MAX-MIN;

% -- add a bit more space to the bounding box, in all dimensions
MAX = MAX + scaling/10;
MIN = MIN - scaling/10;
scaling = MAX-MIN;

% -- construct the grid                    
[X1,X2,X3] = meshgrid(  linspace(MIN(1),MAX(1),ng) ,...
                        linspace(MIN(2),MAX(2),ng) ,...
                        linspace(MIN(3),MAX(3),ng) );

grid = reshape([X1(:),X2(:),X3(:)],ng^d,d);
dV = prod( scaling / (ng-1) ); % elementary cubic volume of each grid point
else
    grid = [];
end

%% Compute PDF by kernel density estimation
% -- define a bandwidth (variance) parameter for all dimensions
if silverman % based on Silverman's heuristic rule
    % cfr. https://en.wikipedia.org/wiki/Multivariate_kernel_density_estimation
    ss = 1 * (4/(d+2))^(1/(d+4)) * n^(-1/(d+4)) * std(X,[],1);
else
    ss = bw * std(X,[],1); 
end         

if ~isempty(grid)
[kdepar,pdf] = kde(X,ss,grid);  
pdf = reshape(pdf,size(X1)); % reshape the PDF into a 3D volume  
kdepar.grid.X1 = X1;
kdepar.grid.X2 = X2;
kdepar.grid.X3 = X3;
kdepar.grid.dV = dV;

if verbose
figure, hold('on'), colormap('jet')
sz = 10;
handkc = scatter3(  X(:,1) , ...
                    X(:,2) , ...
                    X(:,3) , ...
                    sz , ...
                    'k' , ...
                    'filled' );
isoval = max(pdf(:)) * [ 0.5 0.25 0.10 0.05 ];
for iso=[isoval] % isosurfaces with pdf = 0.005,0.01,0.015
    isosurface(X1,X2,X3,pdf,iso),view(3),alpha(.3)
end 
else
end

else
    kdepar = kde(X,ss); 
    pdf = [];
end

%% Evaluate leave-one-out probabilities of each point
pval = zeros(n,1);
Hval = zeros(n,1);
pvalnorm = zeros(n,1);
Hvalnorm = zeros(n,1);

kernelcenters = kdepar . mu;

% -- compute probability for each point in a leave-one-out fashion
% since all kernels lie at the points themselves, LOO can be emulated post
% hoc, by subtracting a (always the same) incremental bias probability that
% results from the kernel at the current data point
pval = probfun( X , ...
                kdepar.w/sum(kdepar.w) , ...
                kdepar.mu , ...
                kdepar.Sig );    
            
pbias = probfun(    zeros(1,d) , ...
                    1/n , ...
                    zeros(1,d) , ...
                    kdepar.Sig(:,:,1) );
 
% -- compute P values
pval = bsxfun(@minus,pval,pbias)*n/(n-1);   

% -- compute the entropy values
Hval = -log2(pval);

end