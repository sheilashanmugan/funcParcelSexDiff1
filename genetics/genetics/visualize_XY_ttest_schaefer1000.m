
clear
addpath(genpath('/cbica/projects/pncSingleFuncParcel/Replication/Toolbox/Code_mvNMF_l21_ard_v3_release/'));
load('/cbica/projects/funcParcelSexDiff/results/genetics/schaefer1000/ttest/XY_ttest_schaefer1000_net7.mat');

% Create file for workbench visualization
SubjectsFolder = '/cbica/projects/funcParcelSexDiff/programs/freesurfer/5.3.0/subjects/fsaverage5';
% for surface data
surfML = '/cbica/projects/pncSingleFuncParcel/Replication/data/SNR_Mask/fsaverage5/lh.Mask_SNR.label';
mwIndVec_l = read_medial_wall_label(surfML);
Index_l = setdiff([1:10242], mwIndVec_l);
surfMR = '/cbica/projects/pncSingleFuncParcel/Replication/data/SNR_Mask/fsaverage5/rh.Mask_SNR.label';
mwIndVec_r = read_medial_wall_label(surfMR);
Index_r = setdiff([1:10242], mwIndVec_r);

V_lh = gifti;
V_lh.cdata = XY_ttest_vert.XY_ttest_schaefer1000;
outdir = '/cbica/projects/funcParcelSexDiff/results/genetics/schaefer1000/ttest/'
V_lh_File = [outdir 'XY_ttest_net7_lh.func.gii'];
save(V_lh, V_lh_File);


%%Xgenes
load('/cbica/projects/funcParcelSexDiff/results/genetics/schaefer1000/ttest/Xgenes_mean_schaefer1000_net7.mat');
V_lh = gifti;
V_lh.cdata = all_X_genes_mean_vert.mean_X_schaefer1000;
outdir = '/cbica/projects/funcParcelSexDiff/results/genetics/schaefer1000/ttest/'
V_lh_File = [outdir 'mean_X_net7_lh.func.gii'];
save(V_lh, V_lh_File);

%%Ygenes
load('/cbica/projects/funcParcelSexDiff/results/genetics/schaefer1000/ttest/Ygenes_mean_schaefer1000_net7.mat');
V_lh = gifti;
V_lh.cdata = all_Y_genes_mean_vert.mean_Y_schaefer1000;
outdir = '/cbica/projects/funcParcelSexDiff/results/genetics/schaefer1000/ttest/'
V_lh_File = [outdir 'mean_Y_net7_lh.func.gii'];
save(V_lh, V_lh_File);