function [maEvents,maTimestamps] = createEvents(path,dataSetIDs)
    addpath(genpath('H:\dissDat\analysis-tools-master'))
    helper = dir([path dataSetIDs '/openEphysData/*.events']);
    eventsName = helper(1).name;
    [maEvents,maTimestamps] = load_open_ephys_data_faster([path dataSetIDs '/openEphysData/' eventsName]);
    return