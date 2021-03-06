clear
PredictionFolder = '/cbica/projects/funcParcelSexDiff/results/PredictionAnalysis/SVM/2fold_CSelect_Cov_SubIndex/Sex_CovAgeMotion/Permutation/res_MultiTimes/';
for i = 1:100
tmp = load([PredictionFolder '/Time_' num2str(i) '/Accuracy.mat']);
   acc_AllModels1(i, :) = tmp;
tmp = load([PredictionFolder '/Time_' num2str(i) '/Sensitivity.mat']);
   sensitivity_AllModels1(i, :) = tmp;
tmp = load([PredictionFolder '/Time_' num2str(i) '/Specificity.mat']);
   specificity_AllModels1(i, :) = tmp;
end
acc_avg = mean(cell2mat(struct2cell(acc_AllModels1)))
sens_avg = mean(cell2mat(struct2cell(sensitivity_AllModels1)))
spec_avg = mean(cell2mat(struct2cell(specificity_AllModels1)))


PredictionFolder_perm = '/cbica/projects/funcParcelSexDiff/results/PredictionAnalysis/SVM/2fold_CSelect_Cov_SubIndex/Sex_CovAgeMotion/Permutation/res_Permutation/';
for i = 1:1000
tmp = load([PredictionFolder_perm '/Time_' num2str(i) '/Accuracy.mat']);
   acc_AllModels1_perm(i, :) = tmp;
end


P = length(find((cell2mat(struct2cell(acc_AllModels1_perm))) > acc_avg)) / 1000
