clear
addpath(genpath('/cbica/projects/pncSingleFuncParcel/Replication/Toolbox/Code_mvNMF_l21_ard_v3_release/'))
addpath(genpath('/cbica/projects/funcParcelSexDiff/scripts/SVM_scripts/roc/'))

mat_yval = load(['/cbica/projects/funcParcelSexDiff/results/PredictionAnalysis/SVM/2fold_CSelect_Cov_SubIndex/Sex_CovAgeMotion/Permutation/res_MultiTimes/ROC/mat_yval.mat']);
mat_sex = load(['/cbica/projects/funcParcelSexDiff/results/PredictionAnalysis/SVM/2fold_CSelect_Cov_SubIndex/Sex_CovAgeMotion/Permutation/res_MultiTimes/ROC/mat_sex.mat']);
%vec_sex=mat_sex.mat_sex;
%vec_yval=mat_yval.mat_yval;

Label=mat_sex.mat_sex;
DecisionValues=mat_yval.mat_yval;
ROC_Draw_Flag=1

AUC_Calculate_ROC_Draw2(DecisionValues, Label, ROC_Draw_Flag)

%AUC_Calculate_ROC_Draw2(vec_yval, vec_sex, 1);