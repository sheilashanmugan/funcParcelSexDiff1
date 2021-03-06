addpath(genpath('/cbica/projects/pncSingleFuncParcel/Replication/Toolbox/Code_mvNMF_l21_ard_v3_release/'))
addpath(genpath('/cbica/projects/funcParcelSexDiff/scripts/spintest/scripts/'))
fshome = '/cbica/projects/funcParcelSexDiff/programs/freesurfer/5.3.0';

% Step 1: SpinPermuFs.m to obtain 'spins' of the data 
% (or use SpinPermuCIVET.m):

% left and right surfaces (group-averaged at every vertex):
readleft = '/cbica/projects/funcParcelSexDiff/inputData/spintest/wBrainSexAbssum_lh.csv';
readright = '/cbica/projects/funcParcelSexDiff/inputData/spintest/wBrainSexAbssum_rh.csv';
permno = 1000; % how many spins
wsname = sprintf('/cbica/projects/funcParcelSexDiff/results/spintest/SVMwBrainSexAbsSum_perm1000.mat');

SpinPermuFS(readleft,readright,permno,wsname)

% Step 2: pvalvsNull.m to run the spin test
% left and right hemispheres for the second modality:
readleft1 = readleft;
readright1 = readright;
readleft2 = '/cbica/projects/funcParcelSexDiff/inputData/spintest/GamSexAbssum_lh.csv'; 
readright2 = '/cbica/projects/funcParcelSexDiff/inputData/spintest/GamSexAbssum_rh.csv';

% indicate (with 0's and 1's) which vertices in the left and right
% hemispheres are part of the medial wall
[vl, left_labels, ctl] = read_annotation(fullfile(fshome,'/subjects/fsaverage5/label/lh.aparc.a2009s.annot'));
v_exclude_left = left_labels==1644825; % label of vertices in the medial wall is 1644825
[vr,right_labels,ctr] = read_annotation(fullfile(fshome,'/subjects/fsaverage5/label/rh.aparc.a2009s.annot'));
v_exclude_right = right_labels==1644825;

pval=pvalvsNull(readleft1,readright1,readleft2,readright2,permno,wsname, v_exclude_left, v_exclude_right);
pval
csvwrite('/cbica/projects/funcParcelSexDiff/results/spintest/pval_1000.csv',pval);