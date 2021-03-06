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
mkdir(ResultantFolder_Sex);
for i = 1:17

  i

  SexEffects = load([ResultantFolder_Sex '/SexEffect_AtlasLoading_17_Network_' num2str(i) '.mat']);

  %% FDR
  %% display Z value (FDR corrected q < 0.05)
  Gam_Z_FDR_Sig_Vector_All_lh = zeros(10242, 1);
  Gam_Z_FDR_Sig_Vector_All_lh(Index_l) = SexEffects.Gam_Z_FDR_Sig_Vector_All(1:length(Index_l));
  Gam_Z_FDR_Sig_Vector_All_rh = zeros(10242, 1);
  Gam_Z_FDR_Sig_Vector_All_rh(Index_r) = SexEffects.Gam_Z_FDR_Sig_Vector_All(length(Index_l) + 1:end);
  % left hemi
  V_lh = gifti;
  V_lh.cdata = Gam_Z_FDR_Sig_Vector_All_lh;
  V_lh_File = [ResultantFolder_Sex '/Gam_Z_FDR_Sig_Vector_All_Sex_lh_Network_' num2str(i) '.func.gii'];
  save(V_lh, V_lh_File);
  % right hemi
  V_rh = gifti;
  V_rh.cdata = Gam_Z_FDR_Sig_Vector_All_rh;
  V_rh_File = [ResultantFolder_Sex '/Gam_Z_FDR_Sig_Vector_All_Sex_rh_Network_' num2str(i) '.func.gii'];
  save(V_rh, V_rh_File);
  % convert into cifti file
  cmd = ['wb_command -cifti-create-dense-scalar ' ResultantFolder_Sex '/Gam_Z_FDR_Sig_Vector_All_Sex_Network_' ...
         num2str(i) '.dscalar.nii -left-metric ' V_lh_File ' -right-metric ' V_rh_File];
  system(cmd);

  Gam_Z_Vector_All_lh = zeros(10242, 1);
  Gam_Z_Vector_All_lh(Index_l) = SexEffects.Gam_Z_Vector_All(1:length(Index_l));
  Gam_Z_Vector_All_rh = zeros(10242, 1);
  Gam_Z_Vector_All_rh(Index_r) = SexEffects.Gam_Z_Vector_All(length(Index_l) + 1:end);
  save([ResultantFolder_Sex '/SexEffect_AtlasLoading_17_Network_' num2str(i) '.mat'], 'Gam_Z_FDR_Sig_Vector_All_lh', 'Gam_Z_FDR_Sig_Vector_All_rh', 'Gam_Z_Vector_All_lh', 'Gam_Z_Vector_All_rh', '-append');

end