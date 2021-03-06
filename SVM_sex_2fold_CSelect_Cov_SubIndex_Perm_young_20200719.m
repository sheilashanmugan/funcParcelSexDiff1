
clear
addpath(genpath('/cbica/projects/funcParcelSexDiff/scripts/SVM_scripts/scripts_from_zaixu/'))
ProjectFolder = '/cbica/projects/pncSingleFuncParcel/Replication/Revision';
DataFolder = '/cbica/projects/pncSingleFuncParcel/Replication/data';
AtlasLoading_Folder = '/cbica/projects/funcParcelSexDiff/data/Revision/SingleParcellation/SingleAtlas_Analysis/FinalAtlasLoading';
PredictionFolder = '/cbica/projects/pncSingleFuncParcel/Replication/Revision/PredictionAnalysis';
Behavior_Mat = load('/cbica/projects/funcParcelSexDiff/inputData/behavior/Behavior_data_young.mat');
Behavior_Mat.Behavior_data_young.sex_new= zeros(231,1);
Covariate_mat=zeros(231,2);
Covariate_mat(:,1)=Behavior_Mat.Behavior_data_young.AgeYears;
Covariate_mat(:,2)=Behavior_Mat.Behavior_data_young.Motion;

BBLID = Behavior_Mat.Behavior_data_young.BBLID;
load('/cbica/projects/funcParcelSexDiff/results/AtlasLoading/tertiles/young/AtlasLoading/AtlasLoading_All_RemoveZero.mat'); % NonZeroIndex is here

for i = 1:length(BBLID)
    if Behavior_Mat.Behavior_data_young.Sex(i)>1.5
       Behavior_Mat.Behavior_data_young.sex_new(i) = 1;
    else
       Behavior_Mat.Behavior_data_young.sex_new(i) = -1;
    end
end

Subjects_Data = AtlasLoading_All_RemoveZero;
Subjects_Label = Behavior_Mat.Behavior_data_young.sex_new';
Fold_Quantity= 2;
Pre_Method = 'Scale';
CVRepeatTimes = 10;
C_Range = 2.^(-5:10);
%C_Range = 2.^(-2:2);
Covariates= Covariate_mat;

% 100 Repeat
Subjects_Data_Path = '/cbica/projects/funcParcelSexDiff/results/AtlasLoading/tertiles/young/AtlasLoading/Subjects_Data.mat'; % is this the same thing?'/cbica/projects/pncSingleFuncParcel/Chead_Backup/Sheila/runSVM/Subjects_Data.mat';
save(Subjects_Data_Path, 'Subjects_Data');
CVRepeatTimes_Range = [1:100];
ResultantFolder = '/cbica/projects/funcParcelSexDiff/results/PredictionAnalysis/SVM/2fold_CSelect_Cov_SubIndex/Sex_CovAgeMotion/tertiles/young/Permutation/res_MultiTimes';
SVM_2group_NFolds_CSelect_Cov_MultiTimes(Subjects_Data_Path, Subjects_Label, Fold_Quantity, C_Range, Covariates, CVRepeatTimes_Range, 0, ResultantFolder)

% Permutation, 1000 times
Subjects_Data_Path = '/cbica/projects/funcParcelSexDiff/results/AtlasLoading/tertiles/young/AtlasLoading/Subjects_Data.mat'; % is this the same thing?'/cbica/projects/pncSingleFuncParcel/Chead_Backup/Sheila/runSVM/Subjects_Data.mat';
save(Subjects_Data_Path, 'Subjects_Data');
CVRepeatTimes_Range = [1:1000];
ResultantFolder = '/cbica/projects/funcParcelSexDiff/results/PredictionAnalysis/SVM/2fold_CSelect_Cov_SubIndex/Sex_CovAgeMotion/tertiles/young/Permutation/res_Permutation';
SVM_2group_NFolds_CSelect_Cov_MultiTimes(Subjects_Data_Path, Subjects_Label, Fold_Quantity, C_Range, Covariates, CVRepeatTimes_Range, 1, ResultantFolder)

