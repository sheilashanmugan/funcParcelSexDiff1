clear
addpath(genpath('/cbica/projects/pncSingleFuncParcel/Replication/Toolbox/Code_mvNMF_l21_ard_v3_release/'))
ReplicationFolder = '/cbica/projects/pncSingleFuncParcel/sheila/pncSingleFuncParcel/Replication';
PredictionFolder = '/cbica/projects/funcParcelSexDiff/results/PredictionAnalysis/SVM/2fold_CSelect_Cov_SubIndex/Sex_CovAgeMotion/Permutation/res_MultiTimes';
%Results_Cell = g_ls([PredictionFolder '/Time_*/w_Brain.mat']);
Results_Cell_a = ls([PredictionFolder '/Time_*/w_Brain.mat']);
Results_Cell= textscan( Results_Cell_a, '%s', 'delimiter', '\n' );
%for i = 1:length(Results_Cell)
for i = 1:100
tmp = load([PredictionFolder '/Time_' num2str(i) '/w_Brain.mat']);
 % tmp = load(Results_Cell{i});
  w_Brain_AllModels1(i, :) = tmp.w{1,1};
  w_Brain_AllModels2(i, :) = tmp.w{1,2};
end
w_Brain_AllModels=vertcat(w_Brain_AllModels1,w_Brain_AllModels2);
w_Brain_Sex = mean(w_Brain_AllModels);
VisualizeFolder = [PredictionFolder '/AtlasLoading/WeightVisualize_Sex_SVM_2fold_CSelect_Cov_MultiTimes'];
mkdir(VisualizeFolder);
save([VisualizeFolder '/w_Brain_Sex.mat'], 'w_Brain_Sex');

load('/cbica/projects/funcParcelSexDiff/results/AtlasLoading/AtlasLoading_All_RemoveZero.mat'); % NonZeroIndex was here


% Create file for workbench visualization
SubjectsFolder = '/cbica/software/external/freesurfer/centos7/5.3.0/subjects/fsaverage5';
% for surface data
surfML = '/cbica/projects/pncSingleFuncParcel/Replication/data/SNR_Mask/fsaverage5/lh.Mask_SNR.label';
mwIndVec_l = read_medial_wall_label(surfML);
Index_l = setdiff([1:10242], mwIndVec_l);
surfMR = '/cbica/projects/pncSingleFuncParcel/Replication/data/SNR_Mask/fsaverage5/rh.Mask_SNR.label';
mwIndVec_r = read_medial_wall_label(surfMR);
Index_r = setdiff([1:10242], mwIndVec_r);

%%%%%%%%%%%%%%%%%%
% Sex Prediction %
%%%%%%%%%%%%%%%%%%
VertexQuantity = 17734;
%% Display the weight of the first 25% regions with the highest absolute weight
mkdir([VisualizeFolder '/First25Percent']);
w_Brain_Sex_Abs = abs(w_Brain_Sex);
[~, Sorted_IDs] = sort(w_Brain_Sex_Abs);
w_Brain_Sex_FirstPercent = w_Brain_Sex;
w_Brain_Sex_FirstPercent(Sorted_IDs(1:round(length(Sorted_IDs) * 0.75))) = 0;
w_Brain_Sex_FirstPercent_All = zeros(1, 17734*17);
w_Brain_Sex_FirstPercent_All(NonZeroIndex) = w_Brain_Sex_FirstPercent;
%SexEffect_Folder = [ReplicationFolder '/Revision/GamAnalysis/AtlasLoading/SexEffects'];
for i = 1:17
    i
    w_Brain_Sex_FirstPercent_I = w_Brain_Sex_FirstPercent_All([(i - 1) * VertexQuantity + 1 : i * VertexQuantity]);
    % left hemi
    w_Brain_Sex_FirstPercent_lh = zeros(1, 10242);
    w_Brain_Sex_FirstPercent_lh(Index_l) = w_Brain_Sex_FirstPercent_I(1:length(Index_l));
    V_lh = gifti;
    V_lh.cdata = w_Brain_Sex_FirstPercent_lh';
    V_lh_File = [VisualizeFolder '/First25Percent/w_Brain_Sex_First25Percent_lh_Network_' num2str(i) '.func.gii'];
    save(V_lh, V_lh_File);

    % right hemi
    w_Brain_Sex_FirstPercent_rh = zeros(1, 10242);
    w_Brain_Sex_FirstPercent_rh(Index_r) = w_Brain_Sex_FirstPercent_I(length(Index_l) + 1:end);
    V_rh = gifti;
    V_rh.cdata = w_Brain_Sex_FirstPercent_rh';
    V_rh_File = [VisualizeFolder '/First25Percent/w_Brain_Sex_First25Percent_rh_Network_' num2str(i) '.func.gii'];
    save(V_rh, V_rh_File);
    % convert into cifti file
    cmd = ['wb_command -cifti-create-dense-scalar ' VisualizeFolder '/First25Percent/w_Brain_Sex_First25Percent_Network_' ...
           num2str(i) '.dscalar.nii -left-metric ' V_lh_File ' -right-metric ' V_rh_File];
    system(cmd);
end

%% Display sum absolute weight of the 17 maps
w_Brain_Sex_All = zeros(1, 17734*17);
w_Brain_Sex_All(NonZeroIndex) = w_Brain_Sex;
%% Display weight of all regions
for i = 1:17
    w_Brain_Sex_Matrix(i, :) = w_Brain_Sex_All([(i - 1) * VertexQuantity + 1 : i * VertexQuantity]);
end
save([VisualizeFolder '/w_Brain_Sex_Matrix.mat'], 'w_Brain_Sex_Matrix');

w_Brain_Sex_Abs_sum = sum(abs(w_Brain_Sex_Matrix));
w_Brain_Sex_Abs_sum_lh = zeros(1, 10242);
w_Brain_Sex_Abs_sum_lh(Index_l) = w_Brain_Sex_Abs_sum(1:length(Index_l));
w_Brain_Sex_Abs_sum_rh = zeros(1, 10242);
w_Brain_Sex_Abs_sum_rh(Index_r) = w_Brain_Sex_Abs_sum(length(Index_l) + 1:end);
save([VisualizeFolder '/w_Brain_Sex_Abs_sum.mat'], 'w_Brain_Sex_Abs_sum', ...
                         'w_Brain_Sex_Abs_sum_lh', 'w_Brain_Sex_Abs_sum_rh');

V_lh = gifti;
V_lh.cdata = w_Brain_Sex_Abs_sum_lh';
V_lh_File = [VisualizeFolder '/w_Brain_Sex_Abs_lh.func.gii'];
save(V_lh, V_lh_File);
pause(1);
V_rh = gifti;
V_rh.cdata = w_Brain_Sex_Abs_sum_rh';
V_rh_File = [VisualizeFolder '/w_Brain_Sex_Abs_sum_rh.func.gii'];
save(V_rh, V_rh_File);
% combine 
cmd = ['wb_command -cifti-create-dense-scalar ' VisualizeFolder '/w_Brain_Sex_Abs_sum' ...
         '.dscalar.nii -left-metric ' V_lh_File ' -right-metric ' V_rh_File];
system(cmd);
pause(1);
system(['rm -rf ' V_lh_File ' ' V_rh_File]);