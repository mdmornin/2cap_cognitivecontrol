fpath = 'F:\dissDat';
dataPath = 'F:/dissDat/dats';    dirList = dir([dataPath]);  folderNames = {dirList.name}; folderNames = folderNames(~ismember(folderNames ,{'.','..'}));
dataDirs = folderNames;
saveDir = 'F:\dissDat\dats';
addpath(genpath('F:\dissDat'))


for i = 1:length(dataDirs)
    applyCARtoDat([saveDir filesep dataDirs{i}],64);
end

