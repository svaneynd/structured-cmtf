function [ kdepar , pdf ] = kde( X , sigma , grid )

% Reference:
%  Kernel density estimation via diffusion
%  Z. I. Botev, J. F. Grotowski, and D. P. Kroese (2010)
%  Annals of Statistics, Volume 38, Number 5, pages 2916-2957.
%  https://nl.mathworks.com/matlabcentral/fileexchange/58312-kernel-density-estimator-for-high-dimensions?s_tid=prof_contriblnk


[n,d] = size(X);

if numel(sigma) == 1
    sigma = repmat(sigma,n,1);
else
    assert(numel(sigma) == d)
end

% KDE parameters
kdepar = struct;
kdepar . w = ones(1,n)/n;
kdepar . mu = X;
kdepar . Sig = repmat( diag(sigma.^2) , 1 , 1 , n );

% evaluate on grid
if (nargin > 2) & ~isempty(grid)
    pdf = probfun( grid , kdepar.w , kdepar.mu , kdepar.Sig );
else
    pdf = [];
end


end