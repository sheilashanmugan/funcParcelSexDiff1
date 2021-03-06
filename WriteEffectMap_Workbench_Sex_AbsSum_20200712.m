clear
addpath(genpath('/cbica/projects/pncSingleFuncParcel/Replication/Toolbox/Code_mvNMF_l21_ard_v3_release/'))

ResultFolder= '/cbica/projects/funcParcelSexDiff/results';
ResultantFolder = [ResultFolder '/GamAnalysis/AtlasLoading'];

% Create file for workbench visualization
SubjectsFolder = '/cbica/software/external/freesurfer/centos7/5.3.0/subjects/fsaverage5';
% for surface data
surfML = '/cbica/projects/pncSingleFuncParcel/Replication/data/SNR_Mask/fsaverage5/lh.Mask_SNR.label';
mwIndVec_l = read_medial_wall_label(surfML);
Index_l = setdiff([1:10242], mwIndVec_l);
surfMR = '/cbica/projects/pncSingleFuncParcel/Replication/data/SNR_Mask/fsaverage5/rh.Mask_SNR.label';
mwIndVec_r = read_medial_wall_label(surfMR);
Index_r = setdiff([1:10242], mwIndVec_r);


%% Sex
ResultantFolder_Sex = [ResultantFolder '/SexEffects'];
UnthreshAbsSumFolder = [ResultantFolder_Sex '/UnthreshAbsSum'];
mkdir(UnthreshAbsSumFolder);
for i = 1:17
    SexEffects = load([ResultantFolder_Sex '/SexEffect_AtlasLoading_17_Network_' num2str(i) '.mat']);
    SexEffects_Matrix(i, :) = SexEffects.Gam_Z_Vector_All;
end

Gam_Sex_Abs_sum = sum(abs(SexEffects_Matrix));
Gam_Sex_Abs_sum_lh = zeros(1, 10242);
Gam_Sex_Abs_sum_lh(Index_l) = Gam_Sex_Abs_sum(1:length(Index_l));
Gam_Sex_Abs_sum_rh = zeros(1, 10242);
Gam_Sex_Abs_sum_rh(Index_r) = Gam_Sex_Abs_sum(length(Index_l) + 1:end);
save([UnthreshAbsSumFolder '/Gam_Sex_Abs_sum.mat'], 'Gam_Sex_Abs_sum', ...
                         'Gam_Sex_Abs_sum_lh', 'Gam_Sex_Abs_sum_rh');

V_lh = gifti;
V_lh.cdata = Gam_Sex_Abs_sum_lh';
V_lh_File = [UnthreshAbsSumFolder '/Gam_Sex_Abs_sum_lh.func.gii'];
save(V_lh, V_lh_File);
pause(1);
V_rh = gifti;
V_rh.cdata = Gam_Sex_Abs_sum_rh';
V_rh_File = [UnthreshAbsSumFolder '/Gam_Sex_Abs_sum_rh.func.gii'];
save(V_rh, V_rh_File);
% combine 
cmd = ['wb_command -cifti-create-dense-scalar ' UnthreshAbsSumFolder '/Gam_Sex_Abs_sum' ...
         '.dscalar.nii -left-metric ' V_lh_File ' -right-metric ' V_rh_File];
system(cmd);
pause(1);
system(['rm -rf ' V_lh_File ' ' V_rh_File]);


