function [Accuracy, Sensitivity, Specificity] = SVM_2group_NFolds_CSelect_Cov(Subjects_Data, Subjects_Label, Fold_Quantity, C_Range, Covariates, Permutation_Flag, ResultantFolder)
%
% Subject_Data:
%           m*n matrix
%           m is the number of subjects
%           n is the number of features
%
% Subject_Label:
%           array of -1 or 1
%
% Fold_Quantity:
%           number of folds
%
% CVRepeatTimes:
%           Number of repetation of cross-validations
%           Because splitting the subjects into N groups is random, we generally repeat several times, i.e., 100
%
% Alpha_Range:
%           Range of alpha, the regularization parameter balancing the training error and regularization penalty
%           Our previous paper used (2^(-5), 2^(-4), ..., 2^9, 2^10), see Cui and Gong (2018), NeuroImage
%
% Pre_Method:
%           'Normalize' or 'Scale'
%
% Covariates:
%           m * n matrix
%           m is the number of subjects
%           n is the number of covariates
%           if no covariates, please set it []
%
% PermutationFlag:
%           1: doing permutation testing
%           Set 0 if not doing permutation testing
%
% ResultantFolder:
%           the path of folder storing resultant files
%

if ~exist(ResultantFolder, 'dir')
    mkdir(ResultantFolder);
end

[Subjects_Quantity, Feature_Quantity] = size(Subjects_Data);

[Splited_Data, Splited_Data_Label, Origin_ID_Cell] = Split_NFolds(Subjects_Data, Subjects_Label, Fold_Quantity);

predicted_labels = [];
decision_values = [];
for i = 1:Fold_Quantity
    
  disp(['The ' num2str(i) ' iteration!']);
    
  % Select training data and testing data
  test_label = Splited_Data_Label{i};
  test_data = Splited_Data{i};
  if ~isempty(Covariates)
      Covariates_test = Covariates(Origin_ID_Cell{i}, :);
  end
    
  Training_all_data = [];
  Label = [];
  Covariates_training = [];
  for j = 1:Fold_Quantity
      if j == i
          continue;
      end
      Training_all_data = [Training_all_data; Splited_Data{j}];
      Label = [Label; Splited_Data_Label{j}];
      if ~isempty(Covariates)
          Covariates_training = [Covariates_training; Covariates(Origin_ID_Cell{j}, :)];
      end
  end

  if Permutation_Flag
      Training_Quantity = length(Label);
      rng('shuffle')
      RandIndex{i} = randperm(Training_Quantity);
      Label = Label(RandIndex{i});
      save([ResultantFolder filesep 'RandIndex.mat'], 'RandIndex');
  end

  if ~isempty(Covariates)
      [Training_quantity, Covariates_quantity] = size(Covariates_training);
      M = 1;
      for j = 1:Covariates_quantity
          M = M + term(Covariates_training(:, j));
      end
      slm = SurfStatLinMod(Training_all_data, M);

      Training_data_residual = Training_all_data - repmat(slm.coef(1, :), Training_quantity, 1);
      for j = 1:Covariates_quantity
          Training_data_residual = Training_data_residual - ...
              repmat(Covariates_training(:, j), 1, Feature_Quantity) .* repmat(slm.coef(j + 1, :), Training_quantity, 1);
      end
      Training_all_data = Training_data_residual;
  end

  % Select optimal C
  for m = 1:length(C_Range)
    for k = 1:10
      Accuracy_tmp(k) = SVM_2group_NFolds(Training_all_data, Label', Fold_Quantity, C_Range(m));
    end
    Accuracy_Inner(m) = mean(Accuracy_tmp);
  end
  [~, Max_Index] = max(Accuracy_Inner);
  C_Optimal(i) = C_Range(Max_Index);
  disp(Accuracy_Inner);
  disp(C_Optimal(i));
  disp('*****************************************');
  % Scaling to [0 1]
  MinValue = min(Training_all_data);
  MaxValue = max(Training_all_data);
  [rows, columns_quantity] = size(Training_all_data);
  for j = 1:columns_quantity
      Training_all_data(:, j) = (Training_all_data(:, j) - MinValue(j)) / (MaxValue(j) - MinValue(j));
  end

  % classification
  Label = reshape(Label, length(Label), 1);
  Training_all_data = double(Training_all_data);
  model(i) = svmtrain(Label, Training_all_data, ['-t 0 -c ' num2str(C_Optimal(i))]);
  if ~isempty(Covariates)
      [test_quantity, ~] = size(Covariates_test);
      test_data_residual = test_data - repmat(slm.coef(1, :), test_quantity, 1);
      for j = 1:Covariates_quantity
          test_data_residual = test_data_residual - ...
              repmat(Covariates_test(:, j), 1, Feature_Quantity) .* repmat(slm.coef(j + 1, :), test_quantity, 1);
      end
      test_data = test_data_residual;
  end
  % Scale
  MaxValue_New = repmat(MaxValue, length(test_label), 1);
  MinValue_New = repmat(MinValue, length(test_label), 1);
  test_data = (test_data - MinValue_New) ./ (MaxValue_New - MinValue_New);

  % predicts
  test_data = double(test_data);
  [predicted_labels_tmp, ~, ~] = svmpredict(test_label, test_data, model(i));
  predicted_labels = [predicted_labels predicted_labels_tmp'];
    
  w{i} = zeros(1, Feature_Quantity);
  for j = 1 : model(i).totalSV
      w{i} = w{i} + model(i).sv_coef(j) * model(i).SVs(j, :);
  end
  decision_values_tmp = w{i} * test_data' - model(i).rho;
  decision_values = [decision_values decision_values_tmp];
    
end

Origin_ID = [];
for i = 1:length(Origin_ID_Cell)
  Origin_ID = [Origin_ID; Origin_ID_Cell{i}];
end

Group1_Index = find(Subjects_Label(Origin_ID) == 1);
Group0_Index = find(Subjects_Label(Origin_ID) == -1);
Category_group1 = predicted_labels(Group1_Index);
Y_group1 = decision_values(Group1_Index);
Y_group1_subj_index = Origin_ID(Group1_Index);
Category_group0 = predicted_labels(Group0_Index);
Y_group0 = decision_values(Group0_Index);
Y_group0_subj_index = Origin_ID(Group0_Index);

mkdir(ResultantFolder);
save([ResultantFolder filesep 'Y.mat'], 'Y_group1', 'Y_group0', 'C_Optimal', 'Y_group1_subj_index', 'Y_group0_subj_index');
save([ResultantFolder filesep 'Category.mat'], 'Category_group1', 'Category_group0');

group0_Wrong_ID = find(Category_group0 == 1);
group0_Wrong_Quantity = length(group0_Wrong_ID);
group1_Wrong_ID = find(Category_group1 == -1);
group1_Wrong_Quantity = length(group1_Wrong_ID);
disp(['group0: ' num2str(group0_Wrong_Quantity) ' subjects are wrong ' mat2str(group0_Wrong_ID) ]);
disp(['group1: ' num2str(group1_Wrong_Quantity) ' subjects are wrong ' mat2str(group1_Wrong_ID) ]);
save([ResultantFolder filesep 'WrongInfo.mat'], 'group0_Wrong_Quantity', 'group0_Wrong_ID', 'group1_Wrong_Quantity', 'group1_Wrong_ID');
Accuracy = (Subjects_Quantity - group0_Wrong_Quantity - group1_Wrong_Quantity) / Subjects_Quantity;
disp(['Accuracy is ' num2str(Accuracy) ' !']);
save([ResultantFolder filesep 'Accuracy.mat'], 'Accuracy');
Sensitivity = (length(Group0_Index) - group0_Wrong_Quantity) / length(Group0_Index);
Specificity = (length(Group1_Index) - group1_Wrong_Quantity) / length(Group1_Index);
disp(['Sensitivity is ' num2str(Sensitivity) ' !']);
save([ResultantFolder filesep 'Sensitivity.mat'], 'Sensitivity');
disp(['Specificity is ' num2str(Specificity) ' !']);
save([ResultantFolder filesep 'Specificity.mat'], 'Specificity');
% save contribution map
w{1} = w{1} / norm(w{1});
w{2} = w{2} / norm(w{2});
save([ResultantFolder filesep 'w_Brain.mat'], 'w');
save([ResultantFolder filesep 'Origin_ID_Cell.mat'], 'Origin_ID_Cell');

