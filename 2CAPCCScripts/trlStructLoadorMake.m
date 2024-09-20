function trlStruct = trlStructLoadorMake(dataPath,checkPath,savePath,tablePath)
    
% First determine if trlStruct already exists at checkPath
% If trlStruct already exists, load it and end function. If not, continue.
if exist(checkPath,'file')
    tic;
    load(checkPath)
    toc;
    disp(['File exists and has been loaded into the workspace after ' num2str(toc) 'seconds'])
    trlStruct = trlStruct;
    return 
end

% 
if exist(tablePath,'file')
    load(tablePath);
else
    warning('No data table found at passed location.')
end


% Generate dataset IDs
if exist(dataPath,'dir')
    dirList = dir(dataPath);
    fileNames = {dirList.name};
    fileNames = fileNames(~ismember(fileNames, {'.','..'}));
    dataSetIDs = fileNames;
else
    warning('No data found')
    return
end

trlStruct = struct();

for i = 1:length(dataSetIDs)
    load([dataPath,filesep,dataSetIDs{i}])
    dlcTimestamps = bodyCoords{1};
    nosePos = bodyCoords{2};
    headPos = bodyCoords{3};
    headPos = headPos(:,1:2);
    nosePos = nosePos(:,1:2);
    LeftSipperPos = boxCoords(1,:);
    RightSipperPos = boxCoords(2,:);
    LSipDist = zeros(1,length(dlcTimestamps));
    RSipDist = zeros(1,length(dlcTimestamps));
    LeftTrials = 1:24;
    RightTrials = 25:48;

    for j = 1:length(dlcTimestamps)
        LSipDist(j) = pdist([nosePos(j,:); RightSipperPos]);    % Invert
        RSipDist(j) = pdist([nosePos(j,:); LeftSipperPos]);     % Invert
    end
    
    tRange = 25; % In Seconds
    tRangeAdj = tRange*30; % Value converted to 30 fps

    trlLSipDist = zeros(length(trialTimes),tRangeAdj+1);
    trlRSipDist = zeros(length(trialTimes),tRangeAdj+1);
    
    

    for k = 1:length(trialTimes)
        [~,tMin] = min(abs(trialTimes(k) - 5 - dlcTimestamps));
        [~,tMax] = min(abs(trialTimes(k) + 20 - dlcTimestamps)); 
        if diff([tMin,tMax]) == 749
            tMax = tMax+1;
        end
        if diff([tMin,tMax]) == 750
            trlLSipDist(k,:) = LSipDist(tMin:tMax);
            trlRSipDist(k,:) = RSipDist(tMin:tMax);
            trlnosePos(k,:,:) = nosePos(tMin:tMax,:);
            trlheadPos(k,:,:) = headPos(tMin:tMax,:);
        end
    end
    
    % Update structure

    trlStruct(i).Session = dataSetIDs{i};
    trlStruct(i).trlLSipDist = trlLSipDist;
    trlStruct(i).trlRSipDist = trlRSipDist;
    trlStruct(i).trlnosePos = trlnosePos;
    trlStruct(i).trlheadPos = trlheadPos;
    % Add Session Type Conditional

    if startsWith(masterTbl.SessionType{i},'Regular')
        LcorrectIdx = (approach(LeftTrials,1,1)) == 1 & approach(LeftTrials,2,1) == 0;
        LincorrectIdx = (approach(LeftTrials,2,1)) == 1;
        RcorrectIdx = (approach(RightTrials,1,1)) == 1 & approach(RightTrials,2,1) == 0;
        RincorrectIdx = (approach(RightTrials,2,1)) == 1;
        LcorrectApproachTrls = trlLSipDist(LeftTrials,:);
        RcorrectApproachTrls = trlRSipDist(RightTrials,:);
    elseif startsWith(masterTbl.SessionType{i},'Reversal')
        LcorrectIdx = (approach(LeftTrials,2,1)) == 1 & approach(LeftTrials,1,1) == 0;
        LincorrectIdx = (approach(LeftTrials,1,1)) == 1;

        RcorrectIdx = (approach(RightTrials,2,1)) == 1 & approach(RightTrials,1,1) == 0;
        RincorrectIdx = (approach(RightTrials,1,1)) == 1;

        LcorrectApproachTrls = trlRSipDist(LeftTrials,:);
        RcorrectApproachTrls = trlLSipDist(RightTrials,:);
    else
        LcorrectIdx = [];
        LincorrectIdx = [];

        RcorrectIdx = [];
        RincorrectIdx = [];

        LcorrectApproachTrls = [];
        RcorrectApproachTrls = [];
    end

    % Clean vars
    LcorrectIdx(isnan(LcorrectIdx))=0; LincorrectIdx(isnan(LincorrectIdx))=0;
    RcorrectIdx(isnan(RcorrectIdx))=0; RincorrectIdx(isnan(RincorrectIdx))=0;


    LcorrectIdx = logical(LcorrectIdx); LincorrectIdx = logical(LincorrectIdx);
    RcorrectIdx = logical(RcorrectIdx); RincorrectIdx = logical(RincorrectIdx);

    % Index Vars
    % Correct
    LDistCorrectAppr = LcorrectApproachTrls(LcorrectIdx,:);
    RDistCorrectAppr = RcorrectApproachTrls(RcorrectIdx,:);
    % Incorrect
    LDistInCorrectAppr = LcorrectApproachTrls(LincorrectIdx,:);
    RDistInCorrectAppr = RcorrectApproachTrls(RincorrectIdx,:);

    % Pull into structure
    trlStruct(i).LcorrectApproachTrls = LcorrectApproachTrls;
    trlStruct(i).RcorrectApproachTrls = RcorrectApproachTrls;
    %
    trlStruct(i).LDistCorrectAppr = LDistCorrectAppr;
    trlStruct(i).RDistCorrectAppr = RDistCorrectAppr;
    %
    trlStruct(i).LDistInCorrectAppr = LDistInCorrectAppr;
    trlStruct(i).RDistInCorrectAppr = RDistInCorrectAppr;
    %
    trlStruct(i).LcorrectIdx = LcorrectIdx;
    trlStruct(i).RcorrectIdx = RcorrectIdx;
    trlStruct(i).LincorrectIdx = LincorrectIdx;
    trlStruct(i).RincorrectIdx = RincorrectIdx;
    %
    trlStruct(i).approach = approach;
    trlStruct(i).bodyCoords = bodyCoords;
    trlStruct(i).boxCoords = boxCoords;
    trlStruct(i).trialTimes = trialTimes;
    trlStruct(i).trialTypes = trialTypes;

    trlStruct(i).NeuronYield = min(size(stMtx));
    
    progress = i/length(dataSetIDs) * 100;

    disp(['Processed ' num2str(progress) '% of datasets'])
end

    if exist(savePath,'dir')
        save([savePath filesep 'trlStruct.mat'],'trlStruct')
    else
        mkdir(savePath)
        save([savePath filesep 'trlStruct.mat'],'trlStruct')
    end

    return

end
