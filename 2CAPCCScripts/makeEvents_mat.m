%% Pull Events Into MATLAB, Save File
addpath(genpath('H:\dissDat\analysis-tools-master'))
path = 'H:\dissDat\Data\';
load(['H:\dissDat\masterTable.mat'])
dataSetIDs = masterTbl.DataSetName;

for i = 1:length(dataSetIDs)
    helper = dir([path dataSetIDs{i} '/openEphysData/*.events']);
    eventsName = helper(1).name;
    [maEvents,maTimestamps] = load_open_ephys_data_faster([path dataSetIDs{i} '/openEphysData/' eventsName]);
    save([path dataSetIDs{i} '/maEvents.mat'],'maEvents','maTimestamps')
end