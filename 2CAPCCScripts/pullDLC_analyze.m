clc
%% 
leftLight = 0;
leftSipper = 2;
rightLight = 1;
rightSipper = 3;
addpath(genpath('H:\dissDat\analysis-tools-master'))
path = 'H:\dissDat\Data\';
load(['H:\dissDat\masterTable.mat'])
dataSetIDs = masterTbl.DataSetName;
affix = 'DLC_resnet50_2CAPFullRun1Mar28shuffle1_1030000.csv';
for j = 1:length(dataSetIDs)
    strSplit = split(dataSetIDs{j});
    prefixDLC = [strSplit{1} strSplit{3} '_video'];
    dlcFileName = [path dataSetIDs{j} '\' prefixDLC affix];
    opts = detectImportOptions(dlcFileName,'NumHeaderLines',1);
    opts.VariableNamesLine = 2;
    opts.VariableUnitsLine = 3;
    opts.DataLine = 4;
    dlcTable = readtable(dlcFileName,opts);
    dlcTable.bodyparts = round(((dlcTable.bodyparts+1)/30),3);
    dlcTimeStamps = dlcTable.bodyparts;
    firstLight = dlcTable.LeftLight_2 > 0.5 | dlcTable.RightLight_2 > 0.5;
    firstLightTime = dlcTimeStamps(find(firstLight,1));
    % Load OE Data & Create Offset Data to Match DLC
    helper = dir([path dataSetIDs{j} '/openEphysData/*.events']);
    eventsName = helper(1).name;
    [events,OEtimestamps] = load_open_ephys_data_faster([path dataSetIDs{j} '/openEphysData/' eventsName]);
    offset = diff([OEtimestamps(1),firstLightTime]);
    OEtimestamps = OEtimestamps + offset;
    OEevents =  [events,OEtimestamps];
    % Create PSTH of DLC Data
    rightSipperTimes = OEevents(OEevents(:,1) == rightSipper,2);
    leftSipperTimes = OEevents(OEevents(:,1) == leftSipper,2);
    rightLightTimes = OEevents(OEevents(:,1) == rightLight,2);
    leftLightTimes = OEevents(OEevents(:,1) == leftLight,2);
    % Reshape Data
    numTrials = length(rightSipperTimes)/2;
    rightSipperTimes = reshape(rightSipperTimes,2,numTrials);
    leftSipperTimes = reshape(leftSipperTimes,2,numTrials);

    if strcmp(masterTbl.CSType{j},'B')
        divLight = 8;
    else
        divLight = 2;
    end
    rightLightTimes = reshape(rightLightTimes,2,length(rightLightTimes)/2);
    leftLightTimes = reshape(leftLightTimes,2,length(leftLightTimes)/2);

    diffRightSipper = diff(rightSipperTimes);
    diffLeftSipper = diff(leftSipperTimes);

    diffRightLight = diff(rightLightTimes);
    diffLeftLight = diff(leftLightTimes);

    offsetSipperLight = 5;

    idxRightSipperTrials = rightSipperTimes(1,round(diffRightSipper) == 8);
    idxLeftSipperTrials = leftSipperTimes(1,round(diffLeftSipper) == 8);
    
    idxRightSipperTrialsCSMIN = rightSipperTimes(1,round(diffRightSipper) ~= 8);
    idxLeftSipperTrialsCSMIN = leftSipperTimes(1,round(diffLeftSipper) ~= 8);
  
    CSPlusON = [idxRightSipperTrials - 5 idxLeftSipperTrials - 5];
    CSMinON = [idxRightSipperTrialsCSMIN - 5 idxLeftSipperTrialsCSMIN - 5];

    CSPlusON = round(CSPlusON,3);
    CSMinON = round(CSMinON,3);

    % Find Positions, Orient to Trial Starts

    dlcleftSipperLocation = mean([dlcTable.LeftSipper dlcTable.LeftSipper_1]);
    dlcrightSipperLocation = mean([dlcTable.RightSipper dlcTable.RightSipper_1]);

    animalSnoutXY = [dlcTable.SnoutTip dlcTable.SnoutTip_1];
    animalHeadXY = [dlcTable.HeadCap dlcTable.HeadCap_1];

    snoutCSPlusPSTH = []; headsCSPlusPSTH = [];
    snoutCSMinPSTH = []; headsCSMinPSTH = [];
    CorrectSipper = [];

for k = 1:length(CSPlusON)
    [~,csPlusIdxMin] = min(abs(CSPlusON(k)-5 - dlcTimeStamps));
    [~,csPlusIdxMax] = min(abs(CSPlusON(k)+20 - dlcTimeStamps));
    snoutCSPlusPSTH(:,:,k) = animalSnoutXY(csPlusIdxMin:csPlusIdxMax,:);
    headsCSPlusPSTH(:,:,k) = animalHeadXY(csPlusIdxMin:csPlusIdxMax,:);
    [~,csMinIdxMin] = min(abs(CSMinON(k)-5 - dlcTimeStamps));
    [~,csMinIdxMax] = min(abs(CSMinON(k)+20 - dlcTimeStamps));
    snoutCSMinPSTH(:,:,k) = animalSnoutXY(csMinIdxMin:csMinIdxMax,:);
    headsCSMinPSTH(:,:,k) = animalSnoutXY(csMinIdxMin:csMinIdxMax,:);
end

if strcmp(masterTbl.SessionType{j},'Regular')
    CorrectSipper = [repmat(2,length(idxRightSipperTrials),1); repmat(3,length(idxLeftSipperTrials),1)]';
else
    CorrectSipper = [repmat(3,length(idxRightSipperTrials),1); repmat(2,length(idxLeftSipperTrials),1)]';
end

trialStructure(j).snoutCSPlusPSTH = snoutCSPlusPSTH;
trialStructure(j).headCSPlusPSTH = headsCSPlusPSTH;
trialStructure(j).snoutCSMinPSTH = snoutCSMinPSTH;
trialStructure(j).headsCSMinPSTH = headsCSMinPSTH;
trialStructure(j).leftSipperLocation = dlcleftSipperLocation;
trialStructure(j).rightSipperLocation = dlcrightSipperLocation;
trialStructure(j).CorrectSipper = CorrectSipper;
trialStructure(j).events = OEevents;



end