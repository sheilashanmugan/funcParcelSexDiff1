clear
PredictionFolder_young = '/cbica/projects/funcParcelSexDiff/results/PredictionAnalysis/SVM/2fold_CSelect_Cov_SubIndex/Sex_CovAgeMotion/tertiles/young/Permutation/res_MultiTimes/';
for i = 1:100
tmp = load([PredictionFolder_young '/Time_' num2str(i) '/Accuracy.mat']);
   acc_AllModels1_young(i, :) = tmp;
tmp = load([PredictionFolder_young '/Time_' num2str(i) '/Sensitivity.mat']);
   sensitivity_AllModels1_young(i, :) = tmp;
tmp = load([PredictionFolder_young '/Time_' num2str(i) '/Specificity.mat']);
   specificity_AllModels1_young(i, :) = tmp;
end
acc_avg_young = mean(cell2mat(struct2cell(acc_AllModels1_young)))
sens_avg_young = mean(cell2mat(struct2cell(sensitivity_AllModels1_young)))
spec_avg_young = mean(cell2mat(struct2cell(specificity_AllModels1_young)))



clear
PredictionFolder_middle = '/cbica/projects/funcParcelSexDiff/results/PredictionAnalysis/SVM/2fold_CSelect_Cov_SubIndex/Sex_CovAgeMotion/tertiles/middle/Permutation/res_MultiTimes/';
for i = 1:100
tmp = load([PredictionFolder_middle '/Time_' num2str(i) '/Accuracy.mat']);
   acc_AllModels1_middle(i, :) = tmp;
tmp = load([PredictionFolder_middle '/Time_' num2str(i) '/Sensitivity.mat']);
   sensitivity_AllModels1_middle(i, :) = tmp;
tmp = load([PredictionFolder_middle '/Time_' num2str(i) '/Specificity.mat']);
   specificity_AllModels1_middle(i, :) = tmp;
end
acc_avg_middle = mean(cell2mat(struct2cell(acc_AllModels1_middle)))
sens_avg_middle = mean(cell2mat(struct2cell(sensitivity_AllModels1_middle)))
spec_avg_middle = mean(cell2mat(struct2cell(specificity_AllModels1_middle)))




clear
PredictionFolder_old = '/cbica/projects/funcParcelSexDiff/results/PredictionAnalysis/SVM/2fold_CSelect_Cov_SubIndex/Sex_CovAgeMotion/tertiles/old/Permutation/res_MultiTimes/';
for i = 1:100
tmp = load([PredictionFolder_old '/Time_' num2str(i) '/Accuracy.mat']);
   acc_AllModels1_old(i, :) = tmp;
tmp = load([PredictionFolder_old '/Time_' num2str(i) '/Sensitivity.mat']);
   sensitivity_AllModels1_old(i, :) = tmp;
tmp = load([PredictionFolder_old '/Time_' num2str(i) '/Specificity.mat']);
   specificity_AllModels1_old(i, :) = tmp;
end
acc_avg_old = mean(cell2mat(struct2cell(acc_AllModels1_old)))
sens_avg_old = mean(cell2mat(struct2cell(sensitivity_AllModels1_old)))
spec_avg_old = mean(cell2mat(struct2cell(specificity_AllModels1_old)))

