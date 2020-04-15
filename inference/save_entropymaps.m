function [ ] = save_entropymaps( hrfinfo , vtemplate , params , Vatlas , clustsize )
% SAVE_ENTROPYMAPS takes as input the parcelwise HRF variability, computed
% using differnet metrics, and converts these to whole-brain images.
% 
% INPUTS
% - statinfo = structure containing the output of a statistical testing
% procedure
% - vtemplate = template volume information to store images ('spm_vol') 
% - params = structure with several parameters regarding directories etc.
% - Vatlas = 3D volume of the brain atlas which is used to transform the
% spatial signatures (over atlas ROIs) back to voxel space
% - clustsize = cluster extent threshold (default: 1 voxel)
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
suffix_thr = '_thr';
suffix_im = '.img';

%% Store maps for all computed HRF variability metrics
try
    hrfinfo.res = rmfield(hrfinfo.res,'perm');
catch
end

% -- extract the names of all HRF variability types
enttypes = fieldnames(hrfinfo.res);

% -- create and save an image for each metric
for typeidx = 1 : numel(enttypes)
    e0 = hrfinfo.res.(enttypes{typeidx}).e0;
    
    % construct a whole-brain image out of the statistical signature
    Vied = wholebrainfromrois( e0 , Vatlas );
    vtemplate.fname = fullfile(  basedir , ...
                            [ enttypes{typeidx} , suffix_im ] );
    spm_write_vol( vtemplate , Vied );

    % ... and also store the thresholded maps if available
    try
    iedmask = hrfinfo.masks.(enttypes{typeidx}).pos;
    [ Vied_thr , frac ] = cluster_correct( wholebrainfromrois( e0 .* iedmask , Vatlas ) , clustsize );
    vtemplate.fname = fullfile(  basedir , ...
                            [ enttypes{typeidx} , suffix_thr , suffix_im ] );
    spm_write_vol( vtemplate , Vied_thr );
    catch % no thresholded map present
    end
end

end
