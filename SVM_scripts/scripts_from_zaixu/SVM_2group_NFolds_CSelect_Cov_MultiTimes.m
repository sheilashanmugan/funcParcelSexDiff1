function SVM_2group_NFolds_CSelect_Cov_MultiTimes(Subjects_Data_Path, Subjects_Label, Fold_Quantity, C_Range, Covariates, CVRepeatTimes_Range, Permutation_Flag, ResultantFolder) 
%
% Subject_Data_Path:
%           Path of a .mat file that stores a m*n matrix
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
% ResultantFolder:
%           the path of folder storing resultant files
%

if ~exist(ResultantFolder, 'dir')
    mkdir(ResultantFolder);
end

LogFolder = [ResultantFolder '/logs'];
mkdir(LogFolder);
for i = CVRepeatTimes_Range
    ResultantFolder_I = [ResultantFolder '/Time_' num2str(CVRepeatTimes_Range(i))];
    ResultantFile = [ResultantFolder_I '/Accuracy.mat'];
    if ~exist(ResultantFile, 'file')
        mkdir(ResultantFolder_I);
        ScriptPath = [ResultantFolder_I '/script.sh'];
        ParametersPath = [ResultantFolder_I '/Parameters.mat'];
        save(ParametersPath, 'Subjects_Label', 'Fold_Quantity', 'C_Range', 'Covariates', 'Permutation_Flag', 'ResultantFolder_I'); 
        Cmd = ['/cbica/software/external/matlab/R2018A/bin/matlab -nosplash -nodesktop -r ' ...
           '"addpath(genpath(''/cbica/projects/pncSingleFuncParcel/Chead_Backup/Sheila''));' ...
           'load(''' Subjects_Data_Path ''');' ...
           'Parameters = load(''' ParametersPath ''');' ...
           'SVM_2group_NFolds_CSelect_Cov(Subjects_Data, Parameters.Subjects_Label, Parameters.Fold_Quantity, Parameters.C_Range, Parameters.Covariates, Parameters.Permutation_Flag, Parameters.ResultantFolder_I);' ...
           'exit(1)">"' ...
           LogFolder '/log_' num2str(CVRepeatTimes_Range(i)) '" 2>&1'];
        fid = fopen(ScriptPath, 'w');
        fprintf(fid, Cmd);
        system(['qsub -l h_vmem=30G ' ScriptPath]);
    end
end


