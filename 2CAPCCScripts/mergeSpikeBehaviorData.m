%% Merge spikes and other data into one loadable variable
% Paths for existing data
addpath(genpath('F:\dissDat\analysis-tools-master'))
path = 'F:\dissDat';
load(['F:\dissDat\masterTable.mat'])
% Paths for spike data

ksresdirpath = 'F:\dissDat\carDats';
ksdir = dir(ksresdirpath);
ksfolders = [ksdir.isdir];
allKsFolders = ksdir(ksfolders);
allKsFoldersNames = {allKsFolders(3:end).name};

addpath(genpath('F:\dissDat\analysis-tools-master'))
path = 'F:\dissDat';
load(['F:\dissDat\restoredScripts\masterTable.mat'])
addpath(genpath('F:\dissDat\rawData'))
fpath = 'F:/dissDat'
NewDataDir = [fpath,filesep,'res'];
DLCDir = [fpath,filesep,'dlcRes'];
DLCsuf = '_videoDLC_resnet50_2CAPFullRun1Mar28shuffle1_1030000.csv';
dirList = dir([fpath,filesep,'rawData']);
folderNames = {dirList([dirList.isdir]).name};
folderNames = folderNames(~ismember(folderNames ,{'.','..'}));
dataSetNames = folderNames;
dataDirs = dataSetNames;
for i = 1:length(dataSetNames)
    splitName = strsplit(dataSetNames{i},'-');
    DLCNames{i} = [splitName{1} '-' splitName{2} '-' splitName{3} splitName{4}];
end
dataSetIDs = DLCNames;
for i = 1:length(dataSetIDs)
    load([path,filesep,'res',filesep,dataSetIDs{i} '_AllData.mat'])
    if exist([ksresdirpath filesep allKsFoldersNames{i} filesep 'spikes.mat'],'file')
        load([ksresdirpath filesep allKsFoldersNames{i} filesep 'spikes.mat'])
    else
        stMtx = []; fr = [];
    end
    save([path,filesep,'res',filesep,dataSetIDs{i},'_AllData.mat'],'sessionInfo','trialTimes','trialTypes','approach','approachParams','boxCoords','bodyCoords','stMtx','fr')

end