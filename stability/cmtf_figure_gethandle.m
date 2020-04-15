function [ plothandle ] = cmtf_figure_gethandle( figparams , mode , r )
% CMTF_FIGURE_GETHANDLE adds a factor or set of signatures of a desired
% source and modes to an existing figure.

%% Check input arguments
% -- find the correct source index
if ~exist('r','var') | isempty(r) | (r<1), r = 0; end % no proper source index given
    assert(r<=figparams.sources.N)

%% Get the handle of the correct subplot
% -- determine which mode the signatures belong to
modeidx = find( cellfun( @(x)~isempty(strfind(x,mode)) , figparams.modes.names ));

if ~isempty(modeidx)
    plothandle = figparams.SourcePane.axhandles{ r , modeidx };
elseif ~isempty(strfind(lower(mode),'hrf'))
    plothandle = figparams.HrfPane.axhandles{1+r}; % time courses / nhrf spatial modulations
else
    warning('No plot handle for this factor was found')
end

end