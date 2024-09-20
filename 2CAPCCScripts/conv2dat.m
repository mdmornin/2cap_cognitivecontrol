fpath = 'F:\dissDat';
dataPath = 'F:/dissDat/rawData';    dirList = dir([dataPath]);  folderNames = {dirList([dirList.isdir]).name}; folderNames = folderNames(~ismember(folderNames ,{'.','..'}));
dataDirs = folderNames;
saveDir = 'F:\dissDat\dats';
addpath(genpath('F:\dissDat'))
for i = 2:length(dataDirs)

    oe2dat([dataPath filesep dataDirs{i} filesep 'openEphysData'],saveDir);

end