function [ b , LL , pen ] = bic( T , model , makefig )
% BIC computes the Bayesian Information Criterion for tensor-valued data T
% and (potentially more than 1) model in CPD format
% 
% INPUTS
% - T = Dth order tensor
% - model = 1 or more CPD approxiations of T, in a cell array
% - makefig = create a figure
%
% OUTPUTS
% - b = value(s) for the BIC
% - LL = value(s) for the log-likelihood, multiplied with -2
% - pen = value(s) for the penalty term that relates to the number of
% parameters
% 
% Reference:
% M. Morup and L. K. Hansen, "Automatic relevance determination for multi-way models"

%% Initialize
if nargin < 3
    makefig = false;
end

nmodels = length(model);
S = numel(T); % number of observations

K = zeros(nmodels,1);
for m = 1:nmodels
    R = size(model{1}{1},2); % CP rank
    K(m) = sum(cellfun(@(x)numel(x),model{m})) - R*(getorder(T)-1);% number of parameters
end

%% Compute BIC
% Compute -2*LL
LL = zeros(nmodels,1);
% for m = 1:nmodels
%     sse = sum(tens2vec(cpdres(T,model{m})).^2);    
%     LL(m) = S*log(sse/S);
% end


prob = zeros(nmodels,1);
for m = 1:nmodels
    Res = cpdres(T,model{m});
    sse = sum(Res(:).^2); % variance = sse/S
%     probi = 1/sqrt(2*pi*sse/S)*exp(-1/(2*sse/S)*(Res.^2));
%     assert(all(probi(:)>0)),assert(all(probi(:)<1))
%     prob(m) = prod(probi(:));
% remark: this is a probability density function; as such the individual
% values may be higher than 1 if sigma^2 is very small. if the data are not
% truly normally distributed, the product of the individual values may also
% be larger than one 
%     sig2 = var(Res(:));
%     LL(m) = S*log(2*pi*sse/S) + S;
%     LL(m) = S*log(2*pi*sig2) + sse/sig2;

    LL(m) = S*log(sse/S) + S*log(2*pi) + S;
end

% Compute pen
pen = K.*log(S);

% Compute BIC
b = LL + pen;

%% Make a figure
if makefig
    figure,hold('on')
    plot(LL,'--r','linewidth',2)
    plot(pen,'--b','linewidth',2)
    plot(b,'-k','linewidth',4)
    legend('log-likelihood','number of parameters','BIC')
    xlabel('models')
end


end