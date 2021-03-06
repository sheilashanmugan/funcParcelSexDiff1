
clear
addpath(genpath('/cbica/projects/funcParcelSexDiff/scripts/SVM_scripts/scripts_from_zaixu/'))
ProjectFolder = '/cbica/projects/pncSingleFuncParcel/Replication/Revision';
DataFolder = '/cbica/projects/pncSingleFuncParcel/Replication/data';
ResultantFolder = '/cbica/projects/funcParcelSexDiff/results/PredictionAnalysis/SVM/2fold_CSelect_Cov_SubIndex/Sex_CovAgeMotion/Permutation';
AtlasLoading_Folder = '/cbica/projects/funcParcelSexDiff/data/Revision/SingleParcellation/SingleAtlas_Analysis/FinalAtlasLoading';
PredictionFolder = '/cbica/projects/pncSingleFuncParcel/Replication/Revision/PredictionAnalysis';
Behavior_Mat = load([PredictionFolder '/Behavior_693.mat']);
Behavior_Mat.sex_new= zeros(693,1);
Covariate_mat=zeros(693,2);
Covariate_mat(:,1)=Behavior_Mat.AgeYears;
Covariate_mat(:,2)=Behavior_Mat.Motion;

BBLID = Behavior_Mat.BBLID;
load('/cbica/projects/funcParcelSexDiff/results/AtlasLoading/AtlasLoading_All_RemoveZero.mat'); % NonZeroIndex is here

for i = 1:length(BBLID)
    if Behavior_Mat.Sex(i)>1.5
       Behavior_Mat.sex_new(i) = 1;
    else
       Behavior_Mat.sex_new(i) = -1;
    end
end

Subjects_Data = AtlasLoading_All_RemoveZero;
Subjects_Label = Behavior_Mat.sex_new';
Fold_Quantity= 2;
Pre_Method = 'Scale';
CVRepeatTimes = 10;
C_Range = 2.^(-5:10);
%C_Range = 2.^(-2:2);
Covariates= Covariate_mat;

% 100 Repeat
Subjects_Data_Path = '/cbica/projects/pncSingleFuncParcel/Chead_Backup/Sheila/runSVM/Subjects_Data.mat';
save(Subjects_Data_Path, 'Subjects_Data');
CVRepeatTimes_Range = [1:100];
ResultantFolder = '/cbica/projects/funcParcelSexDiff/results/PredictionAnalysis/SVM/2fold_CSelect_Cov_SubIndex/Sex_CovAgeMotion/Permutation/res_MultiTimes';
SVM_2group_NFolds_CSelect_Cov_MultiTimes(Subjects_Data_Path, Subjects_Label, Fold_Quantity, C_Range, Covariates, CVRepeatTimes_Range, 0, ResultantFolder)

% Permutation, 1000 times
Subjects_Data_Path = '/cbica/projects/pncSingleFuncParcel/Chead_Backup/Sheila/runSVM/Subjects_Data.mat';
save(Subjects_Data_Path, 'Subjects_Data');
CVRepeatTimes_Range = [1:1000];
ResultantFolder = '/cbica/projects/funcParcelSexDiff/results/PredictionAnalysis/SVM/2fold_CSelect_Cov_SubIndex/Sex_CovAgeMotion/Permutation/res_Permutation';
SVM_2group_NFolds_CSelect_Cov_MultiTimes(Subjects_Data_Path, Subjects_Label, Fold_Quantity, C_Range, Covariates, CVRepeatTimes_Range, 1, ResultantFolder)

