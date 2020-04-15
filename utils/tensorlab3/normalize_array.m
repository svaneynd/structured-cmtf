% NORMALIZE_ARRAY normalizes a multiway array in one or more modes. If more
% than 1 mode should be normalized, the procedure is carried out
% iteratively. If no modes are specified, the norms in all modes are
% computed.
%
% INPUT
% - X = multiway array of arbitrary order M (assumed to be in 'full' format)
% - modes = the modes over which normalization should be carried out. if a 
%           mode m is normalized: all array slabs of order M-1 are scaled to have
%           the same frobenius norm, for varying index k = 1:Jm, in which Jm is the
%           dimension of mode m
% - nit = maximum number of iterations, if more than one mode is normalized
%
% OUTPUT
% - Y = normalized array
% - f = 1xM cell array containing the frobenius norms of Y in all modes 
% - S = cell array of scaling matrices in all modes
%
% Author: Simon Van Eyndhoven

function [ Y , f , S ] = normalize_array( X , modes , nit )

ordX = getorder(X);

%% Check input arguments
if nargin < 3
    nit = 20;
elseif nargin < 2
    modes = repmat(false,1,ordX);
end

if numel(modes) == 1
    if modes == 0
        modes = repmat(false,1,ordX);
    else
        assert((modes>0)&(modes<=ordX),'Choose a valid one mode to normalize.')
        temp = false(1,ordX);
        temp(modes) = true;
        modes = temp;
    end
elseif numel(modes) == ordX
    modes = boolean(modes);
else
    modes = unique(modes);
    assert(all(modes<=ordX))
    assert(all(modes>0))
    temp = false(1,ordX);
    temp(modes) = true;
    modes = temp;
end
clear('temp')

%% Normalize the array
% -- Initialize all scaling matrices
S = cell(1,ordX);
for m = 1:ordX
    S{m} = sparse(diag(ones(size(X,m),1)));
end

% -- If only one mode needs to be normalized: no iterative procedure needed
if numel(find(modes)) == 1
    m = find(modes);
    % compute the frobenius norm of all array slabs for varying index in mode m
    N = tens2mat(X,[],m);
    N = sqrt(sum(N.^2,1));
    % find the scaling matrix by inverting the norms
    S{m} = sparse(diag(1./N));
    Y = tmprod(X,S{m},m);
else
% -- If multiple modes need to be normalized: iterate
Xtemp = X;
for it = 1:nit
    % loop over all modes
    for m = 1:ordX
        % check whether this mode needs to be normalized
        if modes(m)
            % compute the frobenius norm of all array slabs for varying index in mode m
            N = tens2mat(Xtemp,[],m);
            N = sqrt(sum(N.^2,1));
            N = N / mean(N); % normalization: same slab norm for all indices (unit norm is not attainable for more than 1 mode)       
            % find the scaling matrix by inverting the norms
            S{m} = S{m}*sparse(diag(1./N));
            Xtemp = tmprod(Xtemp,sparse(diag(1./N)),m);
        else
        end        
    end
end
Y = Xtemp;
clear('Xtemp')
end

% -- compute the norms in all modes
f = cell(1,ordX);
for m = 1:ordX
    % compute the frobenius norm of all array slabs for varying index in mode m
    N = tens2mat(Y,[],m);
    N = sqrt(sum(N.^2,1));
    f{m} = N(:);        
end

end

