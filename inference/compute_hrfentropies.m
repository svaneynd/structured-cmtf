function [ entropies , enttypes ] = compute_hrfentropies( hB , htc )
% COMPUTE_HRFENTROPIES computes different metrics that quantify how unusual
% each ROI's HRF is (~outlier HRF), given the HRFs of all other ROIs, based 
% on the weighting coefficients B, or the full HRF time courses. (The term
% 'entropies' is used as an umbrella term for all these kinds of metrics.)
% 
% <!> We have empirically found that the full HRF waveform in a ROI is
% reproducible, but that the basis functions and associated weighting
% coefficients B are not. Hence, there may be an unresolved ambiguity
% between these quantities, making them unreliable to compute metrics.
% 
% INPUTS
% - hB = (K x Iv) matrix holding the HRF basis function weight coefficients
% - htc = (time points x Iv) matrix holding the full HRF in each ROI
%
% OUTPUTS
% - entropies = cell array containing different metrics that quantify which
% HRFs may be outliers (each metric as a vector of length Iv)
% - enttypes = names of all metrics
% 
% Author: Simon Van Eyndhoven (simon.vaneyndhoven@gmail.com)

enttypes = {    'Bmahal' ,          false ; ... % Mahalanobis distance for weight coefficients (sign-sensitive)
                'Bmahalmin' ,       false ; ... % Mahalanobis distance for weight coefficients (sign-insensitive)
                'Bentropy' ,        false ; ... % entropy for weight coefficients (sign-sensitive)
                'Bentropymin' , 	false ; ... % entropy for weight coefficients (sign-insensitive)
                'TCextremitysigned' , true ; ... % extremity for HRF time course (sign-sensitive)
                'TCextremityunsign' , true ; ... % extremity for HRF time course (sign-insensitive)
                'TCentropysigned' , true ; ... % entropy for HRF time course (sign-sensitive)
                'TCentropyunsign' , true , ... % entropy for HRF time course (sign-insensitive)
            };
ntypes = size(enttypes,1);
entropies = cell(1,ntypes);

[K,Iv] = size(hB);
nt = size(htc,1);

% truncate time course if needed (samples at which there is no variation
% over the HRFs, need to be removed for numerical reasons)
tincl = true(nt,1);
for tidx = 1 : nt
    if std(htc(tidx,:)) < 1e-6
        tincl(tidx) = false;
    else
    end
end
htc = htc(tincl,:);

if enttypes{1,2}
    Hval = sqrt(mahal(hB',hB'));
    entropies{1} =  Hval;
end
if enttypes{2,2}
    Hval = sqrt(mahal(hB',[hB';-hB']));
    entropies{2} =  Hval;
end
if enttypes{3,2}
    [ ~ , Hval , ~ , ~ ] = kde_entropy_fixcenters( hB' , 2 );  
    entropies{3} = Hval - min(Hval);
end
if enttypes{4,2}
    [ ~ , Hval , ~ , ~ ] = kde_entropy_fixcenters( [hB';-hB'] , 2 );
    Hval = Hval(1:Iv);  
    entropies{4} = Hval - min(Hval);
end
if enttypes{5,2}
    hrfdist = sqrt(max(2-2*cossimil(htc),0))/2;
    entropies{5} = 10 * mean(hrfdist,2);
end
if enttypes{6,2}
    hrfdist = (1-abs(cossimil(htc)));
    entropies{6} = 10 * mean(hrfdist,2);
end
if enttypes{7,2}
    [~,Hval, ~ , ~ ] = kde_entropy_fixcenters( htc(2:end,:)' , 2 );
    entropies{7} = Hval - min(Hval);
end
if enttypes{8,2}
    [~,Hval, ~ , ~ ] = kde_entropy_fixcenters( [htc(2:end,:)';-htc(2:end,:)'] , 2 );
    Hval = Hval(1:Iv);
    entropies{8} = Hval - min(Hval);
end

include = cell2mat(enttypes(:,2));
enttypes = enttypes(include,1);
entropies = entropies(find(include));

end