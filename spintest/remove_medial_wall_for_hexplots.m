clear
addpath(genpath('/cbica/projects/pncSingleFuncParcel/Replication/Toolbox/Code_mvNMF_l21_ard_v3_release/'))


% Create file for workbench visualization
SubjectsFolder = '/cbica/projects/funcParcelSexDiff/programs/freesurfer/5.3.0/subjects/fsaverage5';
% for surface data
surfML = '/cbica/projects/pncSingleFuncParcel/Replication/data/SNR_Mask/fsaverage5/lh.Mask_SNR.label';
mwIndVec_l = read_medial_wall_label(surfML);
Index_l = setdiff([1:10242], mwIndVec_l);
surfMR = '/cbica/projects/pncSingleFuncParcel/Replication/data/SNR_Mask/fsaverage5/rh.Mask_SNR.label';
mwIndVec_r = read_medial_wall_label(surfMR);
Index_r = setdiff([1:10242], mwIndVec_r);



load('/cbica/projects/funcParcelSexDiff/results/PredictionAnalysis/SVM/2fold_CSelect_Cov_SubIndex/Sex_CovAgeMotion/Permutation/res_MultiTimes/AtlasLoading/WeightVisualize_Sex_SVM_2fold_CSelect_Cov_MultiTimes/w_Brain_Sex_Abs_sum.mat'); % NonZeroIndex was here
load('/cbica/projects/funcParcelSexDiff/results/GamAnalysis/AtlasLoading/SexEffects/UnthreshAbsSum/Gam_Sex_Abs_sum.mat'); 


w_Brain_Sex_Abs_sum_rh_NoMedialWall = w_Brain_Sex_Abs_sum_rh(Index_r);
w_Brain_Sex_Abs_sum_lh_NoMedialWall = w_Brain_Sex_Abs_sum_lh(Index_l);
w_Brain_Sex_Abs_sum_All_NoMedialWall = [w_Brain_Sex_Abs_sum_lh(:, Index_l) w_Brain_Sex_Abs_sum_rh(:, Index_r)];


Gam_Sex_Abs_sum_rh_NoMedialWall = Gam_Sex_Abs_sum_rh(Index_r);
Gam_Sex_Abs_sum_lh_NoMedialWall = Gam_Sex_Abs_sum_lh(Index_l);
Gam_Sex_Abs_sum_All_NoMedialWall = [Gam_Sex_Abs_sum_lh(:, Index_l) Gam_Sex_Abs_sum_rh(:, Index_r)];



 save(['/cbica/projects/funcParcelSexDiff/results/PredictionAnalysis/SVM/2fold_CSelect_Cov_SubIndex/Sex_CovAgeMotion/Permutation/res_MultiTimes/AtlasLoading/WeightVisualize_Sex_SVM_2fold_CSelect_Cov_MultiTimes/w_Brain_Sex_Abs_sum_NoMedialWall.mat'], 'w_Brain_Sex_Abs_sum_All_NoMedialWall', 'w_Brain_Sex_Abs_sum_rh_NoMedialWall', 'w_Brain_Sex_Abs_sum_lh_NoMedialWall');

 save(['/cbica/projects/funcParcelSexDiff/results/GamAnalysis/AtlasLoading/SexEffects/UnthreshAbsSum/Gam_Sex_Abs_sum_NoMedialWall.mat'], 'Gam_Sex_Abs_sum_All_NoMedialWall', 'Gam_Sex_Abs_sum_rh_NoMedialWall', 'Gam_Sex_Abs_sum_lh_NoMedialWall');
