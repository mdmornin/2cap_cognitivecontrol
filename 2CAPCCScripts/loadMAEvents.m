fpath = 'F:/dissDat/rawData';
addpath(genpath([fpath,filesep,'analysis-tools-master']))
dirList = dir([fpath]);
folderNames = {dirList([dirList.isdir]).name};
folderNames = folderNames(~ismember(folderNames ,{'.','..'}));
dataDirs = folderNames;


for i = 1:length(dataDirs)
    eventsPath = [fpath,filesep,dataDirs{i},filesep,'openEphysData/all_channels.events'];
    [maEvents,maTimestamps] = load_open_ephys_data_faster(eventsPath);
    savePath = [fpath,filesep,dataDirs{i},filesep,'maEvents.mat'];
    save(savePath,'maEvents','maTimestamps');
end

