function [ r , f0 ] = universesize( lambda , mode , lims )
%UNIVERSESIZE estimates the number of significant eigenvalues

if ~exist('mode','var') | isempty(mode)
    mode = 'norm';
else
    if ~strcmp(mode,'norm')
        mode = 'abs';
    end
end

% -- check if limits for the universe size have been provided
if ~exist('lims','var')
    lims = repmat(length(lambda),1,2);
elseif numel(lims) == 1;
    lims(2) = length(lambda);
else
end

% -- apply a correction: the eigenvalues come in principle from a positive
% semi-definite matrix and should hence be positive
if ~isempty(find(lambda<0))
    lims(2) = min(find(lambda<0)) - 1;
end

[lambda,order] = sort(lambda,'descend');

f = abs(lambda(1:end-1)-lambda(2:end));

if strcmp(mode,'norm')
    f = f./(abs(lambda(1:end-1))+abs(lambda(2:end))+(1e-9)*lambda(1)); % regularization by a fraction of the largest eigenvalue to prevent division by (nearly) zero
elseif strcmp(mode,'abs')
else
    error('Mode not valid')
end

f0 = f;
f(setdiff(1:length(f),lims(1):lims(2))) = 0;

[~,r] = max(f);

end

