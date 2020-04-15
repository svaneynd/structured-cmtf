function [ csol , scaling , sinkmodes , compnorms , C ] = cpdsigncorrect( sol , criterion , sinkmode )
% CPDSIGNCORRECT takes as input a CPD solution, and applies sign correction to
% the factors, in order to obtain signatures that are as non-negative as
% possible (according to the chosen criterion).
% 
% INPUTS
% - sol = CPD solution as a cell array
% - criterion = the criterion which is used to correct the sign of signatures (see below)
%               'energy': 'most energy should be due to positive samples'
%               'mean': 'the mean of the signature should be positive'
%               'maxampl': 'the maximal amplitude should be attained by a positive sample'
% - sinkmode = to maintain the same solution, one mode (the 'sink mode') has to compensate
%               for the accumulative sign corrections in all other modes (hence, acts as
%               a sink for the aggregate sign). When sinkmode is specified as zero, the
%               sink mode is automatically determined (for every component!) to be the 
%               mode that has the most ambiguous score for the chosen criterion. 
%               (This mode also receives the full amplitude of the
%               component)
%
% OUTPUTS
% - csol = CPD solution to which the necessary corrections have been applied
% - scaling = the applied scaling factors for all modes and signatures
% - sinkmodes = the selected sink mode for every component
% - compnorms = norms of all signatures
% - C = scores of all signatures on the criterion
%
% Author: Simon Van Eyndhoven (simon.vaneyndhoven@gmail.com)

criteria = {    'energy' , ... % 'most energy should be due to positive samples'
                'mean' , ... % 'the mean of the signature should be positive'
                'maxampl', ... % 'the maximal amplitude should be attained by a positive sample'
                };

%% Check input arguments
ordT = length(sol); % tensor order
R = size(sol{1},2); % CPD rank

if ~exist('criterion','var') | isempty(criterion)
    criterion = 'energy';
else
    criterion = lower(criterion);
end
if ~exist('sinkmode','var') | isempty(sinkmode)
    sinkmode = 0;
end

%% Define criteria 
% (all should be independent of the signatures' length, so that comparisons
% across modes are possible)
switch criterion
    case 'energy'
        critfun = @(x) sum(double(x>0) .* (x.^2)) - sum(double(x<0) .* (x.^2)); % difference of 'positive' energy to 'negative' energy
    case 'mean'
        critfun = @(x) mean(x); % mean of all samples
    case 'maxampl'
        critfun = @(x) max( x(x>0) ) - max( abs(x(x<0)) ); % difference of maximal amplitudes of positive samples and negative samples
    otherwise
        error('Choose a valid correction criterion')
end

%% Normalize each signature to remove scaling ambiguity, and compute the criterion score
compnorms = nan(ordT,R); % Frobenius norms of signatures
C = nan(ordT,R); % criterion scores for every mode and every component
scaling = repmat({ones(1,R)},ordT,1); % applied scaling factors to every mode and component

csol = sol;

for r = 1 : R
    for mode = 1 : ordT
        % -- compute norm of the signature and normalize
        compnorms(mode,r) = sqrt(sum(sol{mode}(:,r).^2));
        csol{mode}(:,r) = sol{mode}(:,r) / compnorms(mode,r);
        scaling{mode}(r) = 1/compnorms(mode,r);
        % -- compute the criterion score
        C(mode,r) = critfun( csol{mode}(:,r) );
    end
end

%% Determine, for each component, how the sign should be corrected
if sinkmode == 0
    sinkmodes = zeros(1,R);
    for r = 1 : R
        % -- determine the sink mode
        [~,sinkmodes(r)] = min(abs(C(:,r)));
        % -- apply corrections for the non-sink modes
        freemodes_r = setdiff(1:ordT,sinkmodes(r));
        for mode = freemodes_r
            csol{mode}(:,r) = csol{mode}(:,r) * sign(C(mode,r)) ; % sign
            csol{mode}(:,r) = csol{mode}(:,r) * compnorms(mode,r) ; % amplitude
            scaling{mode}(r) = scaling{mode}(r) * sign(C(mode,r)) * compnorms(mode,r);
        end
        % -- apply corrections for the sink mode
        aggsign = prod( sign(C(freemodes_r,r)) ,1);
        csol{sinkmodes(r)}(:,r) = csol{sinkmodes(r)}(:,r) * aggsign ; % sign
        csol{sinkmodes(r)}(:,r) = csol{sinkmodes(r)}(:,r) * compnorms(sinkmodes(r),r) ; % amplitude
        scaling{sinkmodes(r)}(r) = scaling{sinkmodes(r)}(r) * aggsign * compnorms(sinkmodes(r),r);
    end
elseif sinkmode <= ordT
    sinkmodes = repmat(sinkmode,1,R);
    freemodes = setdiff(1:ordT,sinkmode);
    % -- apply corrections for the non-sink modes
    for mode = freemodes
        csol{mode} = bsxfun( @times , csol{mode} , sign(C(mode,:)) ); % sign
        scaling{mode} = bsxfun( @times , scaling{mode} , sign(C(mode,:)) );
    end
    % -- apply corrections for the sink mode
    aggsign = prod( sign(C(freemodes,:)) ,1);
    csol{sinkmode} = bsxfun( @times , csol{sinkmode} , aggsign ); % sign
    aggnorm = prod( compnorms ,1);
    csol{sinkmode} = bsxfun( @times , csol{sinkmode} , aggnorm ); % amplitude
    scaling{sinkmode} = bsxfun( @times , scaling{sinkmode} , aggsign .* aggnorm );
else
    error('Choose a valid mode as sink.')
end

%% Perform a sanity check to validate the applied changes
[err,P,D] = cpderr( sol , csol );
Dtot = D{1};
for mode = 2 : ordT
    Dtot = Dtot * D{mode};
end    
    
    assert( max(err) < 1e-12 )
    assert( frob(P-eye(R))/sqrt(R) < 1e-12 )
    assert( frob(Dtot-eye(R))/sqrt(R) < 1e-12 )

end