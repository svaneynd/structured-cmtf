function [ csol , facnames ] = convert_cmtffac( sol , withampl )
% CONVERTFAC takes as input an instance/solution of an sCMTF model with 
% multiple structured factors, and converts it to "CPD format", i.e. arranging all 
% factors in a cell array, with the same number of columns per factor.
%
% In the converted solution, the EEG CPD-structured factors are included, 
% as well as the fMRI's spatial loadings of the components (it is assumed 
% that the local HRFs, which are then the only variables left out of the 
% analysis, are uniquely determined from the other factors, as the data are fixed).
%
% INPUTS
% - sol = cell array containing the sCMTF solutions in different
% repetitions of the optimization algorithm
% - withampl = boolean indicating whether the components' amplitudes are
% relevant and should be included
%
% OUTPUTS
% - csol = cell array containing the sCMTF solutions, converted to CPD
% format
% - facnames = cell array containing the names of all factors
%
% Author: Simon Van Eyndhoven (simon.vaneyndhoven@gmail.com)

if nargin < 2 | isempty(withampl)
    withampl = true;
end

if nargin < 3 | isempty(lnc)
    lnc = 0;
else
    assert( lnc>=0 )
end

if iscell( sol ) 
    % -- The factors are already in CPD format: check whether dimensions match
    nmodes = length(sol);
    
    Rm = cellfun( @(x)size(x,2) , sol );
    
    uRm = unique(Rm);
    if length(uRm) == 2
        error('Convert solution from (L,L,1)-format to CDP format.')
    elseif length(uRm) > 2
        error('The provided factors should have the same column dimension.')
    else
        csol = sol;
    end
    
elseif isstruct( sol )
    % -- The factors are already in a structure format: convert to cell array
    facnames = fieldnames( sol.factors );
    
    % find the EEG factors
    eegnames = { 'S' , 'G', 'M' };
    eegfacs = cell(1,length(eegnames));
    for f = 1:length(eegnames)
        eegfacs{f} = sol.factors.(eegnames{f});
    end
    
    % find the HRF factors
    hrfname = 'HRFtoep';
    hrffacidx = find( cellfun( @(x)~isempty(x) , strfind( facnames , hrfname )));
    hrfnames = facnames( hrffacidx ); 
    
    % find (and compose) the fMRI factors
    fmrinames = {'V'};
    
    sol = disambiguatefmri(sol);
    
    fmrifacs = sol.variables.spatial(end);
    fmrifacs{1} = fmrifacs{1}';
    
    % find amplitudes
    if withampl
        amplnames = { 'Lambda_x' , 'Lambda_y' };
        amplfacs = cell(1,length(amplnames));
        for f = 1:length(amplnames)
            amplfacs{f} = sol.factors.( amplnames{f} );
        end
    else
        amplnames = {};
        amplfacs = {};
    end
    
    % concatenate all factors
    csol = [ eegfacs , fmrifacs , amplfacs ];
    facnames = [ eegnames , fmrinames , amplnames ];
    
    % do dimension check
    csol = convert_cmtffac_kr( csol);
    
else
    error('The provided solution should be a cell array or structure containing all factors.')
end


end