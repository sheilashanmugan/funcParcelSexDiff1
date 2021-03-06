
clear
addpath(genpath('/cbica/projects/pncSingleFuncParcel/Replication/Toolbox/Code_mvNMF_l21_ard_v3_release/'));
load('/cbica/projects/funcParcelSexDiff/results/genetics/parc_mat/w_Brain_Sex_Abs_sum_schaefer400_lh.mat');

% Create file for workbench visualization
SubjectsFolder = '/cbica/software/external/freesurfer/centos7/5.3.0/subjects/fsaverage5';
% for surface data
surfML = '/cbica/projects/pncSingleFuncParcel/Replication/data/SNR_Mask/fsaverage5/lh.Mask_SNR.label';
mwIndVec_l = read_medial_wall_label(surfML);
Index_l = setdiff([1:10242], mwIndVec_l);
surfMR = '/cbica/projects/pncSingleFuncParcel/Replication/data/SNR_Mask/fsaverage5/rh.Mask_SNR.label';
mwIndVec_r = read_medial_wall_label(surfMR);
Index_r = setdiff([1:10242], mwIndVec_r);

V_lh = gifti;
V_lh.cdata = parc_vert;
outdir = '/cbica/projects/funcParcelSexDiff/results/genetics/parc_mat/'
V_lh_File = [outdir 'w_Brain_Sex_Abs_sum_schaefer400_lh.func.gii'];
save(V_lh, V_lh_File);