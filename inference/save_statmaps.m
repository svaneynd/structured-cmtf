function [ ] = save_statmaps( statinfo , vtemplate , iedcomp , params , Vatlas , clustsize , statname )
% SAVE_STATMAPS takes as input the results of a statistical testing
% procedure on the fMRI spatial signatures, and creates whole-brain images
% of the statistical (de)activation maps for all components in the
% factorization.
% 
% INPUTS
% - statinfo = structure containing the output of a statistical testing
% procedure
% - vtemplate = template volume information to store images ('spm_vol') 
% - iedcomp = index of the component within the factorization that represents the IEDs
% - params = structure with several parameters regarding directories etc.
% - Vatlas = 3D volume of the brain atlas which is used to transform the
% spatial signatures (over atlas ROIs) back to voxel space
% - clustsize = cluster extent threshold (default: 1 voxel)
% - statname = prefix which is added to the filenames of the images
%
% Author: Simon Van Eyndhoven (simon.vaneyndhoven@gmail.com)

%% 0 . Initialize
if isempty(clustsize) | (clustsize<0)
    clustsize = 1;
end


basedir = fullfile( params.stat.resultdir , ...
                    sprintf(params.stat.patientdir,params.data.patient) , ...
                    params.stat.cmtfdir );

% -- define suffixes to the file names to indicate modifications
suffix_ied = '_ied';
suffix_thr = '_thr';
suffix_neg = '_deact';
suffix_im = '.img';

% -- define template name for directory of non-IED components
dirnoied = 'comp%d_nonied';

% -- find the indices of the non-IED components
nr = params.decomposition.R;
    assert(size(statinfo.res.T0,2)==nr)
idxnoied = setdiff(1:nr,iedcomp);

%% 1 . Store map + thresholded map for _activations_ of IED component

% construct a whole-brain image out of the statistical signature
Vied = wholebrainfromrois( statinfo.res.T0(:,iedcomp) , Vatlas );
vtemplate.fname = fullfile(  basedir , ...
                        [ statname , suffix_ied , suffix_im ] );
spm_write_vol( vtemplate , Vied );

    % ... and also store the thresholded T-map
    % iedmask = statinf.Tglm.masks.gmm.pos(:,centroidcomp(iedcomp_relidx)) & statinf.Tglm.masks.perm.pos(:,centroidcomp(iedcomp_relidx)) ;
iedmask = statinfo.masks.perm.pos(:,iedcomp) ;
[ Vied_thr , frac ] = cluster_correct( wholebrainfromrois( statinfo.res.T0(:,iedcomp) .* iedmask , Vatlas ) , clustsize );
% Vied_thr = wholebrainfromrois( statinf.Tglm.res.T0(:,centroidcomp(iedcomp_relidx)) .* iedmask , Vatlas );
vtemplate.fname = fullfile(  basedir , ...
                        [ statname , suffix_thr , suffix_ied , suffix_im ] );
spm_write_vol( vtemplate , Vied_thr );

%% 2 . Store map + thresholded map for _deactivations_ of IED component

% construct a whole-brain image out of the statistical signature
Vied = wholebrainfromrois( -1 * statinfo.res.T0(:,iedcomp) , Vatlas );
vtemplate.fname = fullfile(  basedir , ...
                        [ statname , suffix_neg , suffix_im ] );
spm_write_vol( vtemplate , Vied );

    % ... and also store the thresholded T-map
    % iedmask = statinf.Tglm.masks.gmm.pos(:,centroidcomp(iedcomp_relidx)) & statinf.Tglm.masks.perm.pos(:,centroidcomp(iedcomp_relidx)) ;
iedmask = statinfo.masks.perm.neg(:,iedcomp) ;
[ Vied_thr , frac ] = cluster_correct( -1 * wholebrainfromrois( statinfo.res.T0(:,iedcomp) .* iedmask , Vatlas ) , clustsize );
% Vied_thr = wholebrainfromrois( statinf.Tglm.res.T0(:,centroidcomp(iedcomp_relidx)) .* iedmask , Vatlas );
vtemplate.fname = fullfile(  basedir , ...
                        [ statname , suffix_neg, suffix_thr , suffix_im ] );
spm_write_vol( vtemplate , Vied_thr );

%% 3 . Store maps + thresholded maps for non-IED components
for k = 1 : length(idxnoied)
    noiedcomp = idxnoied(k);
    subdir = sprintf(dirnoied,k+1);
    mkdir(basedir,subdir)
    
    % -- 3.1. activation
    Vnoied = wholebrainfromrois( statinfo.res.T0(:,noiedcomp) , Vatlas );
    vtemplate.fname = fullfile(  basedir , subdir , ...
                            [ statname , suffix_im ] );
    spm_write_vol( vtemplate , Vnoied );    
    
    % -- 3.2. thresholded activation
    noiedmask = statinfo.masks.perm.pos(:,noiedcomp) ;
    [ Vnoied_thr , frac ] = cluster_correct( wholebrainfromrois( statinfo.res.T0(:,noiedcomp) .* noiedmask , Vatlas ) , clustsize );
    vtemplate.fname = fullfile(  basedir , subdir , ...
                            [ statname , suffix_thr , suffix_im ] );
    spm_write_vol( vtemplate , Vnoied_thr );
    
    % -- 3.3. deactivation
    Vnoied = wholebrainfromrois( -1 * statinfo.res.T0(:,noiedcomp) , Vatlas );
    vtemplate.fname = fullfile(  basedir , subdir , ...
                            [ statname , suffix_neg , suffix_im ] );
    spm_write_vol( vtemplate , Vnoied );    
    
    % -- 3.4. thresholded deactivation
    noiedmask = statinfo.masks.perm.neg(:,noiedcomp) ;
    [ Vnoied_thr , frac ] = cluster_correct( wholebrainfromrois( -1 * statinfo.res.T0(:,noiedcomp) .* noiedmask , Vatlas ) , clustsize );
    vtemplate.fname = fullfile(  basedir , subdir , ...
                            [ statname , suffix_neg , suffix_thr , suffix_im ] );
    spm_write_vol( vtemplate , Vnoied_thr );
    
end

end
