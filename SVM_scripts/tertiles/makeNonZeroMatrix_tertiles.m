%young
clear
ResultantFolder = '/cbica/projects/funcParcelSexDiff/results/AtlasLoading/tertiles/young';
AtlasLoading_Folder = '/cbica/projects/funcParcelSexDiff/data/Revision/SingleParcellation/SingleAtlas_Analysis/FinalAtlasLoading';
PredictionFolder = '/cbica/projects/pncSingleFuncParcel/Replication/Revision/PredictionAnalysis';
Behavior_Mat = load('/cbica/projects/funcParcelSexDiff/inputData/behavior/Behavior_data_young.mat');
BBLID = Behavior_Mat.Behavior_data_young.BBLID;

for i = 1:length(BBLID)
    i
    tmp = load([AtlasLoading_Folder '/' num2str(BBLID(i)) '.mat']); 
    sbj_AtlasLoading_NoMedialWall_Tmp = tmp.sbj_AtlasLoading_NoMedialWall;
    [rowQuantity, colQuantity] = size(sbj_AtlasLoading_NoMedialWall_Tmp);
    AtlasLoading_All(i, :) = reshape(sbj_AtlasLoading_NoMedialWall_Tmp, 1, rowQuantity * colQuantity);
end

AtlasLoading_Sum = sum(AtlasLoading_All);
NonZeroIndex = find(AtlasLoading_Sum);
AtlasLoading_All_RemoveZero = AtlasLoading_All(:, NonZeroIndex);
mkdir([ResultantFolder '/AtlasLoading']);
save([ResultantFolder '/AtlasLoading/AtlasLoading_All.mat'], 'AtlasLoading_All');
save([ResultantFolder '/AtlasLoading/AtlasLoading_All_RemoveZero.mat'], 'AtlasLoading_All_RemoveZero', 'NonZeroIndex');


%middle
clear
ResultantFolder = '/cbica/projects/funcParcelSexDiff/results/AtlasLoading/tertiles/middle';
AtlasLoading_Folder = '/cbica/projects/funcParcelSexDiff/data/Revision/SingleParcellation/SingleAtlas_Analysis/FinalAtlasLoading';
PredictionFolder = '/cbica/projects/pncSingleFuncParcel/Replication/Revision/PredictionAnalysis';
Behavior_Mat = load('/cbica/projects/funcParcelSexDiff/inputData/behavior/Behavior_data_middle.mat');
BBLID = Behavior_Mat.Behavior_data_middle.BBLID;

for i = 1:length(BBLID)
    i
    tmp = load([AtlasLoading_Folder '/' num2str(BBLID(i)) '.mat']); 
    sbj_AtlasLoading_NoMedialWall_Tmp = tmp.sbj_AtlasLoading_NoMedialWall;
    [rowQuantity, colQuantity] = size(sbj_AtlasLoading_NoMedialWall_Tmp);
    AtlasLoading_All(i, :) = reshape(sbj_AtlasLoading_NoMedialWall_Tmp, 1, rowQuantity * colQuantity);
end

AtlasLoading_Sum = sum(AtlasLoading_All);
NonZeroIndex = find(AtlasLoading_Sum);
AtlasLoading_All_RemoveZero = AtlasLoading_All(:, NonZeroIndex);
mkdir([ResultantFolder '/AtlasLoading']);
save([ResultantFolder '/AtlasLoading/AtlasLoading_All.mat'], 'AtlasLoading_All');
save([ResultantFolder '/AtlasLoading/AtlasLoading_All_RemoveZero.mat'], 'AtlasLoading_All_RemoveZero', 'NonZeroIndex');


%old
clear
ResultantFolder = '/cbica/projects/funcParcelSexDiff/results/AtlasLoading/tertiles/old';
AtlasLoading_Folder = '/cbica/projects/funcParcelSexDiff/data/Revision/SingleParcellation/SingleAtlas_Analysis/FinalAtlasLoading';
PredictionFolder = '/cbica/projects/pncSingleFuncParcel/Replication/Revision/PredictionAnalysis';
Behavior_Mat = load('/cbica/projects/funcParcelSexDiff/inputData/behavior/Behavior_data_old.mat');
BBLID = Behavior_Mat.Behavior_data_old.BBLID;

for i = 1:length(BBLID)
    i
    tmp = load([AtlasLoading_Folder '/' num2str(BBLID(i)) '.mat']); 
    sbj_AtlasLoading_NoMedialWall_Tmp = tmp.sbj_AtlasLoading_NoMedialWall;
    [rowQuantity, colQuantity] = size(sbj_AtlasLoading_NoMedialWall_Tmp);
    AtlasLoading_All(i, :) = reshape(sbj_AtlasLoading_NoMedialWall_Tmp, 1, rowQuantity * colQuantity);
end

AtlasLoading_Sum = sum(AtlasLoading_All);
NonZeroIndex = find(AtlasLoading_Sum);
AtlasLoading_All_RemoveZero = AtlasLoading_All(:, NonZeroIndex);
mkdir([ResultantFolder '/AtlasLoading']);
save([ResultantFolder '/AtlasLoading/AtlasLoading_All.mat'], 'AtlasLoading_All');
save([ResultantFolder '/AtlasLoading/AtlasLoading_All_RemoveZero.mat'], 'AtlasLoading_All_RemoveZero', 'NonZeroIndex');

