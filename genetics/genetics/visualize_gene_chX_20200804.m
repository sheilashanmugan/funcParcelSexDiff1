
clear
addpath(genpath('/cbica/projects/pncSingleFuncParcel/Replication/Toolbox/Code_mvNMF_l21_ard_v3_release/'));
load('/cbica/projects/funcParcelSexDiff/results/genetics/gene_mat/all_genes_mean.mat');
gene_all_vert = all_genes_mean.gene_avg;

% Create file for workbench visualization
SubjectsFolder = '/cbica/software/external/freesurfer/centos7/5.3.0/subjects/fsaverage5';
% for surface data
surfML = '/cbica/projects/pncSingleFuncParcel/Replication/data/SNR_Mask/fsaverage5/lh.Mask_SNR.label';
mwIndVec_l = read_medial_wall_label(surfML);
Index_l = setdiff([1:10242], mwIndVec_l);
surfMR = '/cbica/projects/pncSingleFuncParcel/Replication/data/SNR_Mask/fsaverage5/rh.Mask_SNR.label';
mwIndVec_r = read_medial_wall_label(surfMR);
Index_r = setdiff([1:10242], mwIndVec_r);

V_lh = gifti;
V_lh.cdata = gene_all_vert;
outdir = '/cbica/projects/funcParcelSexDiff/results/genetics/gene_mat/'
V_lh_File = [outdir 'gene_all_lh.func.gii'];
save(V_lh, V_lh_File);










% 
% 
% sigX = load('/cbica/projects/funcParcelSexDiff/results/genetics/sigX.mat')
% genelist = sigX.sigX.gene
% 
% 
% for i = 1:length(genelist)
% matPath = '/cbica/projects/funcParcelSexDiff/results/genetics/gene_mat/';
% gene = genelist{i}
% gene_vert = load([matPath gene '_vert.mat']);
% all_genes(i, :) = gene_vert.gene_vert;
% end
% 
% avg_gene = mean(cell2mat(struct2cell(all_genes)))
% 
% % Create file for workbench visualization
% SubjectsFolder = '/cbica/software/external/freesurfer/centos7/5.3.0/subjects/fsaverage5';
% % for surface data
% surfML = '/cbica/projects/pncSingleFuncParcel/Replication/data/SNR_Mask/fsaverage5/lh.Mask_SNR.label';
% mwIndVec_l = read_medial_wall_label(surfML);
% Index_l = setdiff([1:10242], mwIndVec_l);
% surfMR = '/cbica/projects/pncSingleFuncParcel/Replication/data/SNR_Mask/fsaverage5/rh.Mask_SNR.label';
% mwIndVec_r = read_medial_wall_label(surfMR);
% Index_r = setdiff([1:10242], mwIndVec_r);
% 
% V_lh = gifti;
% V_lh.cdata = gene_vert;
% outdir = '/cbica/projects/pncSingleFuncParcel/sheila/pfn_sex_diff/results/Allen/VisualizeMap/'
% V_lh_File = [outdir gene '_lh.func.gii'];
% save(V_lh, V_lh_File);
% 
% 
% %to read in label files
% [v, L, ct] = read_annotation('/cbica/projects/pncSingleFuncParcel/sheila/pfn_sex_diff/inputData/genetics/lh.V2.label');
% 
% 
% %to read in label files
% [v, L, ct] = read_annotation('/cbica/projects/pncSingleFuncParcel/sheila/pfn_sex_diff/inputData/genetics/lh.500.aparc.annot');