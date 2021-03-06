addpath(genpath('/cbica/projects/pncSingleFuncParcel/Replication/Toolbox/Code_mvNMF_l21_ard_v3_release/'))
Weight_lh = gifti('/cbica/projects/funcParcelSexDiff/results/PredictionAnalysis/SVM/2fold_CSelect_Cov_SubIndex/Sex_CovAgeMotion/Permutation/res_MultiTimes/AtlasLoading/WeightVisualize_Sex_SVM_2fold_CSelect_Cov_MultiTimes/w_Brain_Sex_Abs_lh.func.gii');
Weight_lh_cdata = Weight_lh.cdata;

Weight_rh = gifti('/cbica/projects/funcParcelSexDiff/results/PredictionAnalysis/SVM/2fold_CSelect_Cov_SubIndex/Sex_CovAgeMotion/Permutation/res_MultiTimes/AtlasLoading/WeightVisualize_Sex_SVM_2fold_CSelect_Cov_MultiTimes/w_Brain_Sex_Abs_sum_rh.func.gii');
Weight_rh_cdata = Weight_rh.cdata;
