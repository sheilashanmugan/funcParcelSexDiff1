clear
addpath(genpath('/cbica/projects/pncSingleFuncParcel/Replication/Toolbox/Code_mvNMF_l21_ard_v3_release/'))

ResultFolder= '/cbica/projects/funcParcelSexDiff/results';
ResultantFolder = [ResultFolder '/GamAnalysis/AtlasLoading'];

% Create file for workbench visualization
SubjectsFolder = '/cbica/software/external/freesurfer/centos7/5.3.0/subjects/fsaverage5';
% for surface data
surfML = '/cbica/projects/pncSingleFuncParcel/Replication/data/SNR_Mask/fsaverage5/lh.Mask_SNR.label';
mwIndVec_l = read_medial_wall_label(surfML);
Index_l = setdiff([1:10242], mwIndVec_l);
surfMR = '/cbica/projects/pncSingleFuncParcel/Replication/data/SNR_Mask/fsaverage5/rh.Mask_SNR.label';
mwIndVec_r = read_medial_wall_label(surfMR);
Index_r = setdiff([1:10242], mwIndVec_r);


%% Sex
ResultantFolder_Sex = [ResultantFolder '/SexEffects'];
mkdir(UnthreshAbsSumFolder);
for i = 1:17
    SexEffects = load([ResultantFolder_Sex '/SexEffect_AtlasLoading_17_Network_' num2str(i) '.mat']);
    SexEffects_Matrix(i, :) = SexEffects.Gam_Z_FDR_Sig_Vector_All;
end

 save([ResultantFolder_Sex '/SexEffects_Matrix_Gam_17NetAll_FDR_Sig.mat'], 'SexEffects_Matrix');
