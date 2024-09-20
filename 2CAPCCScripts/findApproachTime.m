%% Loads & Paths
clear all; close all
if ispc 
    addpath(genpath('F:\dissDat\analysis-tools-master'))
    path = 'F:\dissDat';
    load(['F:\dissDat\restoredScripts\masterTable.mat'])
    addpath(genpath('F:\dissDat\res'))
    NewDataDir = [path,filesep,'res'];
    DLCDir = [path,filesep,'dlcRes'];
    DLCsuf = '_videoDLC_resnet50_2CAPFullRun1Mar28shuffle1_1030000.csv';
    dirList = dir([path,filesep,'res']);
    folderNames = {dirList.name};
    folderNames = folderNames(~ismember(folderNames ,{'.','..'}));
    dataSetIDs = folderNames;
    structurePath = [path filesep 'trlStruct.mat'];
elseif isunix
    figuredir = '/research3/figures';
    
    addpath(genpath('/research3/analysis-tools-master'))
    path = '/research3/';
    load(['/research3/restoredScripts/masterTable.mat'])
    addpath(genpath('/research3/res'))
    fpath = '/research3/';
    NewDataDir = [fpath,filesep,'res'];
    dirList = dir([fpath,filesep,'res']);
    folderNames = {dirList.name};
    folderNames = folderNames(~ismember(folderNames ,{'.','..'}));
    dataSetIDs = folderNames;
    stMtx = []; % For now
    structurePath = '/research3/trlStruct.mat';
end
    
load(structurePath)
    
%%
% Box Coordinate Key
% left sipper, right sipper, upper left corner, lower left corner, upper right corner, lower right corner

posRegCSPos = []; posRegCSMin = [];
posRevCSPos = []; posRevCSMin = [];
for i = 1:length(trlStruct)

    [nmBoxCoords,nmCenter,nmScale] = normalize(trlStruct(i).boxCoords);
    for k = 1:length(trlStruct(i).trialTimes)
        nmNosePos(k,:,:) = normalize(squeeze(trlStruct(i).trlnosePos(k,:,:)),'center',nmCenter,'scale',nmScale);
    end

    trlStruct(i).nmNosePos = nmNosePos;
    trlStruct(i).nmBoxCoords = nmBoxCoords;
end

for i = 1:length(dataSetIDs)
    if regexp(masterTbl.SessionType{i},regexptranslate('wildcard','Regular*'))
        posRegCSPos = [posRegCSPos; trlStruct(i).nmNosePos(1:48,:,:)];
        posRegCSMin = [posRegCSMin; trlStruct(i).nmNosePos(49:end,:,:)];
    elseif regexp(masterTbl.SessionType{i},regexptranslate('wildcard','Reversal*'))
        posRevCSPos = [posRevCSPos; trlStruct(i).nmNosePos(1:48,:,:)];
        posRevCSMin = [posRevCSMin; trlStruct(i).nmNosePos(49:end,:,:)];
    end
end

%% 
for i = 1:length(trlStruct)
    tMin = 150;
    tMax = 540;
    idx = tMax - tMin;
    
    testData = [trlStruct(i).nmNosePos(:,tMin:tMax-1,1), trlStruct(i).nmNosePos(:,tMin:tMax-1,2), diff(trlStruct(i).nmNosePos(:,tMin:tMax,1),1,2), diff(trlStruct(i).nmNosePos(:,tMin:tMax,2),1,2)];
    testData = reshape(testData,96,idx,4);
    
    data = testData(1:48,:,:);
    

    lType = trlStruct(1).approach(1:48,1,1) == 1 | trlStruct(1).approach(1:48,2,1) == 1;
    lickLikelihood = [];
    for iTrial = 1:48
        % Calculate the distance in 4-D space between this trial and all other
        % trials for each time step
        dis = sqrt(sum((data - repmat(data(iTrial,:,:),[size(data,1),1,1])).^2,3));
%         dis = diag(pdist2(squeeze(data(iTrial,:,:)),squeeze(data(iTrial,:,:)),'euclidean'),1);

        % The weights are the inverse distance
        weights = 1./dis;
        
        % Normalize the weights. Note, remove this trial because its weight is
        % infinity
        weights = weights./repmat(sum(weights(setdiff(1:size(data,1),iTrial),:),1),[size(data,1),1]);
        
        % Calculate the likelihood using the weighted average
        lickLikelihood(iTrial,:) = sum(weights(setdiff(1:size(data,1),iTrial),:).*repmat(lType(setdiff(1:size(data,1),iTrial)),[1,size(data,2)]),1);
        
    end

    for iTrial = 1:size(lickLikelihood,1)
        changePts(iTrial) = findchangepts(lickLikelihood(iTrial,:));
    end
    trlStruct(i).approachChangePoints = changePts./30;

end
%% Plot the lick likelihoods
figure
hold on

for iTrial = 1:48
    if lType(iTrial) == 1
        plot(1:size(data,2),smooth(lickLikelihood(iTrial,:),3),'r');
    else
        plot(1:size(data,2),smooth(lickLikelihood(iTrial,:),3),'b');
    end
end

xlabel('Time Steps')
ylabel('Lick Likelihood')
legend('Lick Trials','No Lick Trials')
title('Lick Likelihoods')

%% Change Points

for iTrial = 1:size(lickLikelihood,1)
    changePts(iTrial) = findchangepts(lickLikelihood(iTrial,:));
end


figure
hold on
scatter(1:48,changePts)
scatter(1:48,lType*100)
yline(4*30);
yline(9*30);
yline(16*30)

%% Plot a mean


