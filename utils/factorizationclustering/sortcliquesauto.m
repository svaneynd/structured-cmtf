function [ M , cliquegroups , cliqueruns ] = sortcliquesauto( M , minN , thresh )
% SORTCLIQUESAUTO identifies runs which share 1 or more common components. This
% is done by taking as input a set of 'groups' (corresponding to one underlying
% template components), which are lists that indicate in which runs the
% component appears. These groups can then be categorized in a smaller 
% set of 'cliques': one clique hence consists of a number of template
% components (groups) that occur frequently together in several runs.
% Creating sets of cliques is done by checking which groups overlap in many
% runs, and applying average linkage clustering.
%
% Note 1: Groups can hence be considered as component-level/point-level cliques, as the
% runs in which a component appears are seen as fully mutually
% interconnected.
% Note 2: Cliques may be seen as corresponding to local solutions (one
% solution is a set of components) (object-level). These cliques may not be
% perfect (pseudo-cliques), if the participation of runs differs over
% components in the clique.

if ~exist('minN','var')
    minN = 1;
end
if ~exist('thresh','var')    
    thresh = [0.5:0.1:0.9];
end

%% Arrange the input matrix M (= a set of groups)
% every column/group has a length equal to the number of runs
% there are 'ngroups' groups
[nruns,ngroups] = size(M);
M = mat2cell(M,nruns,ones(1,ngroups));
groupruns = cellfun(@(x)find(~isnan(x)),M,'UniformOutput',false);
ngroupruns = cellfun(@(x)numel(x),groupruns);

% find which column/group has the most assigned runs ans sort the
% input accordingly
[ngroupruns,order] = sort(ngroupruns,'descend');
M = M(:,order);
groupruns = groupruns(order);

%% Compute the similarity between different groups
% limit to the number of groups that have more than one associated run
% (more meaningful clustering)
Mred = M(ngroupruns>=minN);
groupruns = groupruns(ngroupruns>=minN);

% compute the Dice similarity
S = cellfun(@(x,y)2*numel( intersect(x,y) )/(numel(x)+numel(y)),...
                    repmat(groupruns',1,length(Mred)),repmat(groupruns,length(Mred),1));

% compute dissimilarity
assert(max(S(:))<=1)

D = sqrt( 1 - S ).^2;    

M = cell2mat(M);    
    
%% Use average linkage clustering to obtain partitions     
strategy = 'AL';
[partition,dendrogram.Z,dendrogram.order]=...
    hcluster(D,strategy);

%% Find the 'best' number of cliques
Rmin = 1;
L = size(D,1);  
% compute the rank of the binarized distance matrix
hf = figure; 
subplot(2,3,4:6), hold('on')
cm = colormap(jet(length(thresh)));
idxt = 1:length(thresh);% round(linspace(1,size(cm,1)/2,length(thresh)));
for t = 1:length(thresh)
    legendlab{t} = sprintf('threshold %2.2f',thresh(t));
    Sbin = double(S>=thresh(t));
    lambda = eig(Sbin);
    lambda = sort(lambda,'descend');
    plot(lambda,'color',cm(idxt(t),:),'linewidth',2);
end
title('Eigenvalues of thresholded similarity matrix')
legend(legendlab)

Rest = universesize(lambda,'norm',Rmin);

pause;
    
prompt = {sprintf('Enter a new value for R (minimal value: %d)\nEstimated value: R = %d at threshold %2.2f',Rmin,Rest,thresh(end))};
prompttitle = 'Choose the number of significant eigenvalues';
promptdims = [1 70];
definput = {sprintf('%d',Rest)};
newr = inputdlg(prompt,prompttitle,promptdims,definput);

ncliquesopt = str2num(newr{1});

figure(hf);  colormap('gray');
subplot(2,3,1),imagesc(~isnan(cell2mat(Mred)')),title('distribution of runs over cliques')
% subplot(2,3,2),barh(1:L,partition(ncliques(1),:),0.7,'FaceColor',[1 0 0],'Linestyle','none'),title('suggested partition')
subplot(2,3,3),imagesc(D),title('dissimilarity between cliques'),colorbar;

ready = false;
while ~ready
    subplot(2,3,2)
    P = zeros(L,L);
    for idx = 1:L
        P(idx,1:partition(ncliquesopt,idx)) = 1;
    end
    imagesc(P),title('suggested partition')
%     barh(1:L,partition(ncliquesopt,:),0.7,'FaceColor',[1 0 0],'Linestyle','none'),title('suggested partition')
    pause;
    prompt = {sprintf('Choose the number of cliques, or enter 0 to stop\nEstimated value: n = %d',ncliquesopt)};
    prompttitle = 'Choose the number of cliques';
    promptdims = [1 70];
    definput = {sprintf('%d',ncliquesopt)};
    newn = inputdlg(prompt,prompttitle,promptdims,definput);
    newn = str2num(newn{1});
    if newn == 0
        ready = true;
    else
        ncliquesopt = newn;
    end
end
close(hf);

%% Arrange the output
cliquegroups = cell(1,max(partition(ncliquesopt,:)));
cliqueruns = cell(1,max(partition(ncliquesopt,:)));

nruns = zeros(ncliquesopt,1);
for c = 1:ncliquesopt
    cliquegroups{c} = find(partition(ncliquesopt,:)==c);
    cliqueruns{c} = groupruns(find(partition(ncliquesopt,:)==c));
    nruns(c) = mean(cellfun(@(x)length(x),groupruns(cliquegroups{c})));
end

[nruns,order] = sort(nruns,'descend');

cliquegroups = cliquegroups(order);
cliqueruns = cliqueruns(order);

end

function [ ri , stats ] =rindex(dist,partition)
%function ri=rindex(dist,partition)

N=size(partition,1);

% Number of clusters in each partition
Ncluster=max(partition');

disp('Computing R-index...');

for k=1:N,
    % compute cluster statistics
    s=clusterstat(dist,partition(k,:),1);
    % Compute R-index: (set diagonal to Inf to exclude self-distance!
    s.between.avg(eye(size(s.between.avg))==1)=0;  
    ri(k,3)=mean(s.internal.avg'./max(s.between.avg));    
    ri(k,9)=median(s.internal.avg'./max(s.between.avg)); 
    
    ri(k,2)=mean(s.internal.avg'./mean(s.between.avg)); 
    ri(k,8)=median(s.internal.avg'./mean(s.between.avg)); 
    
    ri(k,10)=mean(s.internal.avg'./median(s.between.avg)); 
    ri(k,11)=median(s.internal.avg'./median(s.between.avg)); 
    
    % a good cluster is one in which the internal distance is low and
    % the external distance is high

    s.between.avg(eye(size(s.between.avg))==1)=Inf;
    ri(k,1)=mean(s.internal.avg'./min(s.between.avg));     
    ri(k,7)=median(s.internal.avg'./min(s.between.avg));     
     
    
    ri(k,4) = min(s.internal.avg'./s.external.avg'); 
    ri(k,5) = mean(s.internal.avg'./s.external.avg'); 
    ri(k,6) = max(s.internal.avg'./s.external.avg'); 
    
    ri(k,12) = mean(s.internal.avg') - mean(min(s.between.avg));
    ri(k,13) = median(s.internal.avg') - median(min(s.between.avg));
    
    stats{k,1} = s.internal.avg;
    stats{k,2} = s.external.avg;
    stats{k,3} = s.between.avg;
end
end

function Stat=clusterstat(S,partition,between)
%function Stat=clusterstat(S,partition,[between])

if nargin<2,
   error('You must give at least two input arguments.');
end

if nargin<3|isempty(between)
  between=0;
end

% Number of clusters
Ncluster=max(partition);

%Initialize the struct
Stat.internal.sum(1:Ncluster,1)=NaN;
Stat.internal.min(1:Ncluster,1)=NaN;
Stat.internal.avg(1:Ncluster,1)=NaN;
Stat.internal.max(1:Ncluster,1)=NaN;
Stat.external.sum(1:Ncluster,1)=NaN;
Stat.external.min(1:Ncluster,1)=NaN;
Stat.external.avg(1:Ncluster,1)=NaN;
Stat.external.max(1:Ncluster,1)=NaN;

for cluster=1:Ncluster,
  thisPartition=(partition==cluster);
  S_=S(thisPartition,thisPartition);
  Stat.N(cluster)=size(S_,1);
  S_(eye(size(S_))==1)=[];
  if isempty(S_),
          S_ = 0.0; % self-distance
          includeclust(cluster) = true;
  else
      includeclust(cluster) = true;
  end
    Stat.internal.sum(cluster)=sum(S_);
    Stat.internal.min(cluster)=min(S_);
    Stat.internal.avg(cluster)=mean(S_);
    Stat.internal.max(cluster)=max(S_);
  if Ncluster>1,
    S_=S(thisPartition,~thisPartition);
    Stat.external.sum(cluster)=sum(S_(:));
    Stat.external.min(cluster)=min(S_(:));
    Stat.external.avg(cluster)=mean(S_(:));
    Stat.external.max(cluster)=max(S_(:));
  end
  
end

Stat.internal.sum = Stat.internal.sum(includeclust);
  Stat.internal.min = Stat.internal.min(includeclust);
  Stat.internal.avg = Stat.internal.avg(includeclust);
  Stat.internal.max = Stat.internal.max(includeclust);
  
  Stat.external.sum = Stat.external.sum(includeclust);
  Stat.external.min = Stat.external.min(includeclust);
  Stat.external.avg = Stat.external.avg(includeclust);
  Stat.external.max = Stat.external.max(includeclust);

if between,
  Stat.between.min=zeros(Ncluster,Ncluster);
  Stat.between.max=zeros(Ncluster,Ncluster);
  Stat.between.avg=zeros(Ncluster,Ncluster);

  % perform computation for upper triangular part
  for i=1:Ncluster,
    Pi=find(i==partition);
    for j=i+1:Ncluster,
      Pj=find(j==partition);
      d_=S(Pi,Pj); 
      Stat.between.min(i,j)=min(d_(:));
      Stat.between.avg(i,j)=mean(d_(:));		
      Stat.between.max(i,j)=max(d_(:));				
    end
  end  
  
  Stat.between.min = Stat.between.min(includeclust,includeclust);
  Stat.between.avg = Stat.between.avg(includeclust,includeclust);
  Stat.between.max = Stat.between.max(includeclust,includeclust);
  
  % make the matrix symmetric
  Stat.between.min=Stat.between.min+ ...
    Stat.between.min';
  Stat.between.max=Stat.between.max+ ...
      Stat.between.max';
  Stat.between.avg=Stat.between.avg+ ...
      Stat.between.avg';
end
end