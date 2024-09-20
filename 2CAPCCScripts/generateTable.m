% Loaded DataLog3.xlsx and manually removed columns and datasets not needed
% Saved this manually cleaned data as masterTable.mat

load('F:\dissDat\masterTable.mat')

% Add columns to masterTable related to drinking intake and other useful
% information from each day's session

fpath = 'F:/dissDat/rawData';
dirList = dir([fpath]);
folderNames = {dirList([dirList.isdir]).name};
folderNames = folderNames(~ismember(folderNames ,{'.','..'}));
dataDirs = folderNames;

for i = 1:length(dataDirs)
    loadPath = [fpath filesep dataDirs{i} filesep 'sessionInfo.xlsx'];
    sessionTable = readtable(loadPath);
    bottleLabels = {'Initial Left Bottle Weight (g)','Initial Right Bottle Weight (g)','Final Left Bottle Weight (g)','Final Right Bottle Weight (g)'};
    bottleIndex = find(contains(sessionTable.Var1,bottleLabels));
    bottleValues = str2double(sessionTable.Var2(bottleIndex))';
    weightIndex = find(contains(sessionTable.Var1,'Animal Weight (g)'));
    weightValue = str2double(sessionTable.Var2(weightIndex));
    % Extraneous code for inputting values into table - could be improved.
    masterTable.animalWeight(i) = weightValue;
    masterTable.initLWeight(i) = bottleValues(1);
    masterTable.initRWeight(i) = bottleValues(2);
    masterTable.finlLWeight(i) = bottleValues(3);
    masterTable.finlRweight(i) = bottleValues(4);
    masterTable.diffL(i) = bottleValues(1) - bottleValues(3);
    masterTable.diffR(i) = bottleValues(2) - bottleValues(4);
    % Dose (g/kg) of EtOH calcluated by summing the two difference values
    % (g), multiplying by EtOH concentration (10%, 0.1), and then dividing
    % by the animals weight in kg (g/1000)
    masterTable.dose(i) = (sum([masterTable.diffL(i),masterTable.diffR(i)])*0.1) / (masterTable.animalWeight(i)/1000);
end

% All relevant information from Excel sheets has been gained, save
% variable.

save('masterTable.mat','masterTable');

