function [ Xenv , Xenvtot , t_ax , pow ] = envelope( X , tmode , wsh , wlen , pow , fs , im )
% ENVELOPE computes a p-norm envelope for every variable of a higher-order 
% dataset (potentially followed by sliding-window averaging).
%
% Author: Simon Van Eyndhoven (simon.vaneyndhoven@gmail.com)

%% Initialization
OrdX = getorder(X);
assert(tmode>0&tmode<=OrdX)
sX = getsize(X);

X = tens2mat(X,[],tmode);

if isempty(wsh)
    wsh = 1;
else
    wsh = round(wsh*fs);
end
if isempty(wlen)
    wlen = wsh;
else
    wlen = round(wlen*fs);
end
assert(wlen>=wsh,'window length should be greater than the step size')
if isempty(pow)
    pow = 2;
end

if nargin < 7
    im = false;
end

[nvars,ntimes] = size(X);

nw = floor((ntimes-wlen)/wsh)+1;

%% Compute the p-norm envelope
Xenv = zeros(nvars,nw);
for w = 1:nw
    start = (w-1)*wsh+1;
    stop = start + wlen - 1;
    Xenv(:,w) = squeeze(mean(abs(X(:,start:stop,:)).^pow,2));    
end

% -- average the envelopes of all variables in the data to obtain a global envelope
Xenvtot = mean(Xenv,1);

% -- reconstruct the envelope to tensor format
sXenv = sX;
sXenv(tmode) = nw;
modes = 1:OrdX;
modes(tmode) = [];
Xenv = mat2tens(Xenv,sXenv,modes,2);

% -- construct a corresponding time axis (in seconds)
if nargin > 2
    fsenv = fs/wsh;
    assert(logical(exist('fs','var')))
    t_ax = wlen/(fs*2) + (0:nw-1)/fsenv;
end

%% Plot the resulting envelopes
if im
    % -- determine sizes for plotting
    height = 1/(nvars + 1 + 3);
    width = 0.8;
    
    figure
    h1 = subplot('position',[(1-width)/2,3*height,width,nvars*height]);
    imagesc(Xenv)
    h2 = subplot('position',[(1-width)/2,1*height,width,1*height]);
    
    % -- channelwise envelopes
    subplot(h1)
    imagesc(Xenv)
    xt = get(gca,'Xtick');
    set(gca,'Xticklabel',t_ax(xt))

    
    % -- global envelope
    subplot(h2)
    imagesc(Xenvtot)
    set(gca,'Xtick',xt,'Xticklabel',t_ax(xt))
    
end

end