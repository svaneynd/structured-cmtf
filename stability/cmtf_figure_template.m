function [ figparams ] = cmtf_figure_template( R , modes , relspace , figtitle )
% CMTF_FIGURE_TEMPLATE creates a blank frame wherein factor signatures of
% different modes can be plotted later. A frame for the HRF is always
% created (this is not a real 'mode').
% 
% INPUTS

figparams = struct;

%% Check input arguments
if ~exist('modes','var'), modes = []; end % modes not specified
if ~exist('relspace','var'), relspace = []; end % relspace not specified
if ~exist('figtitle','var'), figparams.figtitle = 'CMTF signatures'; else figparams.figtitle = figtitle; end % title not specified

% -- number of sources
figparams.sources.N = R;

% -- modes and attributes
figparams.modes = struct;
if iscell(modes)
    figparams.modes.names = modes;
    figparams.modes.N = length(modes); % number of modes
    if (length(relspace)~=length(modes)) | isempty(relspace)
        figparams.dims.relspace = ones(1,figparams.modes.N);
    else
        figparams.dims.relspace = relspace;
    end
elseif isempty(modes)
    figparams.modes.names = {   'spectral' , ...
                                'spatial (channels)' , ...
                                'time course' , ...
                                'spatial (ROIs)' , ...
                                'occ./ovl.' , ...
                                };
    figparams.modes.N = 5;
    figparams.dims.relspace = [ 5 3 6 5 0.5 ];
elseif isstr(modes)
    figparams.modes.names = modes;
    figparams.modes.N = 1;
    figparams.dims.relspace = 1;
else
    error('Use a valid format for the modes')
end

%% Define default figure parameters
% -- size of the figure
figparams.dims.figsize = [ 200 0 1800 900 ];

% -- portion of the figure devoted to whitespace
figparams.dims.hwhite = 0.025; % horizontal / columnwise: between modes
figparams.dims.vwhite = 0.03; % vertical / rowwise: between sources

% -- figure title pane
figparams.FigTitlePane = struct;
figparams.FigTitlePane.height = 0.05;

% -- mode title pane
figparams.ModeTitlePane = struct;
figparams.ModeTitlePane.height = 0.025;
% figparams.ModeTitlePane.pos = [ 0 , 1-figparams.FigTitlePane.height-figparams.ModeTitlePane.height , 1 , figparams.ModeTitlePane.height ];

% -- source title pane
figparams.SourceTitlePane = struct;
figparams.SourceTitlePane.width = 0.08;
% figparams.SourceTitlePane.pos = [ 0 , 1-figparams.FigTitlePane.height-figparams.ModeTitlePane.height , 1 , figparams.ModeTitlePane.height ];

% -- source pane
figparams.SourcePane = struct;
figparams.SourcePane.height = 0.65;
figparams.SourcePane.width = 0.92;

% -- hrf pane
figparams.HrfPane = struct;
% figparams.HrfPane.height = 0.2;
figparams.HrfPane.width = 0.7;

%% Generate the figure with panels
figparams.fighandle = figure('name',figparams.figtitle,'position',figparams.dims.figsize);

% -- add the figure title
figparams.FigTitlePane.handle = ...
    annotation(	'textbox',[ 0 , 1-figparams.FigTitlePane.height , 1 , figparams.FigTitlePane.height ], ...
                'string',figparams.figtitle , ...
                'fontsize',18 , ...
                'horizontalalignment','center' , ...
                'verticalalignment','top' , ...
                'interpreter','latex', ...
                'linestyle','none');
            
% --  generate panes for sources and modes
figparams.modes.space = (figparams.dims.relspace ./ sum(figparams.dims.relspace)) * ...
                            (figparams.SourcePane.width - figparams.dims.hwhite*(figparams.modes.N+1)) ;
                        
figparams.sources.space = (1 / figparams.sources.N) * ...
                            (figparams.SourcePane.height - figparams.dims.vwhite*(figparams.sources.N+1)) ;
                        
figparams.SourcePane.upperleft = [  figparams.SourceTitlePane.width , ...
                                        1 - figparams.FigTitlePane.height - figparams.ModeTitlePane.height  ];

figparams.SourcePane.axhandles = cell(figparams.sources.N,figparams.modes.N);   

% make subplots for the standard modes                                 
modeidx = 1 : figparams.modes.N; %find( cellfun( @(x)isempty(strfind(x,'occur')) , figparams.modes.names ));
for m = modeidx
   for s = 1 : figparams.sources.N
       % determine where the subplot should be placed
       pos = zeros(1,4);
       pos(1) = figparams.SourcePane.upperleft(1) + m*figparams.dims.hwhite + sum(figparams.modes.space(1:m-1));
       pos(2) = figparams.SourcePane.upperleft(2) - s * (figparams.dims.vwhite + figparams.sources.space);
       pos(3) = figparams.modes.space(m);
       pos(4) = figparams.sources.space;
       % store the handle 
       figparams.SourcePane.axhandles{ s , m } = subplot('position',pos);
   end
end

% % make subplots for the occurrence modes  
% m = find( cellfun( @(x)~isempty(strfind(x,'occur')) , figparams.modes.names ));
%     for s = 1 : figparams.sources.N
%        % determine where the subplot should be placed
%        pos = zeros(1,4);
%        pos(1) = figparams.SourcePane.upperleft(1) + m*figparams.dims.hwhite + sum(figparams.modes.space(1:m-1));
%        pos(2) = figparams.SourcePane.upperleft(2) - figparams.sources.N * (figparams.sources.space + figparams.dims.vwhite);
%        pos(3) = figparams.modes.space(m);
%        pos(4) = figparams.sources.N * (figparams.sources.space + figparams.dims.vwhite) - figparams.dims.vwhite;
%        % store the handle 
%        figparams.SourcePane.axhandles{ s , m } = subplot('position',pos);
%    end

% add titles in the ModeTitlePane
modetitletemplate = '\\textbf{%s}';
for m = 1 : figparams.modes.N
    % determine where the subplot should be placed
    pos = zeros(1,4);
    pos(1) = figparams.SourcePane.upperleft(1) + m*figparams.dims.hwhite + sum(figparams.modes.space(1:m-1));
    pos(2) = figparams.SourcePane.upperleft(2) - 0.5*figparams.dims.vwhite;
    pos(3) = figparams.modes.space(m);
    pos(4) = figparams.ModeTitlePane.height;  
    % add title
    figparams.ModeTitlePane.handles{m} = ...
    annotation(	'textbox',pos, ...
                'string',sprintf(modetitletemplate,figparams.modes.names{m}) , ...
                'fontsize',14 , ...
                'horizontalalignment','center' , ...
                'verticalalignment','bottom' , ...
                'interpreter','latex', ...
                'linestyle','none');
end

% add titles in the SourceTitlePane
sourcetitletemplate = '\\textbf{s%d}';
for s = 1 : figparams.sources.N
    % determine where the subplot should be placed
    pos = zeros(1,4);
    pos(1) = 0;
    pos(2) = figparams.SourcePane.upperleft(2) - s * (figparams.dims.vwhite + figparams.sources.space);
    pos(3) = figparams.SourceTitlePane.width;
    pos(4) = figparams.sources.space;  
    % add title
    figparams.SourceTitlePane.handles{s} = ...
        annotation(	'textbox',pos, ...
                    'string',sprintf(sourcetitletemplate,s) , ...
                    'fontsize',14 , ...
                    'fontweight','normal' , ...
                    'horizontalalignment','center' , ...
                    'verticalalignment','middle' , ...
                    'interpreter','latex', ...
                    'linestyle','none');
end

% --  generate hrf pane
figparams.HrfPane.axhandles = cell(1,1);
figparams.HrfPane.titlehandle = cell(1,1);

% determine where the subplot should be placed
pos = zeros(1,4);
pos(1) = figparams.SourcePane.upperleft(1) + figparams.dims.hwhite;
pos(2) = figparams.dims.vwhite;
pos(3) = figparams.HrfPane.width;  
pos(4) = figparams.SourcePane.upperleft(2) - figparams.sources.N * (figparams.dims.vwhite + figparams.sources.space)- 2*figparams.dims.vwhite; % count from the bottom of the SourcePane until the bottom of the figure

% store the handle 
figparams.HrfPane.axhandles{1} = subplot('position',pos);

% determine where the title should be placed
pos = pos + [ pos(3) + figparams.dims.hwhite , 0 , figparams.SourceTitlePane.width-pos(3) , 0 ]; 

% add the title
hrftitle = '\textbf{HRF basis function(s)}';
figparams.HrfPane.titlehandle{1} = ...
        annotation(	'textbox',pos , ...
                    'string',hrftitle , ...
                    'fontsize',14 , ...
                    'fontweight','normal' , ...
                    'horizontalalignment','center' , ...
                    'verticalalignment','middle' , ...
                    'interpreter','latex', ...
                    'linestyle','none');

end