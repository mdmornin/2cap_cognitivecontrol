% pWistDataOrganizerVer2
%
%   This script organizes data from the P rat and Wistar quinine ephys
%   experiment that included recordings in mPFC. Version 2 operated on
%   version 3 of the data which had been sorted with kilosort1 and reviewed
%   using Nicks GUI. It takes those data and replaces the original spiking
%   data with new spiking data from batch 7 (kilosort 2, reviewed in phy).
%   Version 3 combines batches 7 and 8 into one data set and incorporates
%   results from DeepLabCut.

%% Loads
addpath(genpath('H:\dissDat\analysis-tools-master'))
path = 'H:\dissDat\Data\';
load(['H:\dissDat\masterTable.mat'])
dataSetIDs = masterTbl.DataSetName;
%% Settings
affix = 'DLC_resnet50_2CAPFullRun1Mar28shuffle1_1030000.csv';
% Set the location of the new data sets
NewDataDir = 'H:\dissDat\curated';

% Set the location of the DeepLabCut results
DLCDir = path;

% Set the file name suffix for the DLC results
DLCsuf = 'videoDLC_resnet50_2CAPFullRun1Mar28shuffle1_1030000.csv';

% Set the locations of the data sets on the data capacitor

% Set the names of the data sets. 
dataSetNames = dataSetIDs;
dataDirs = dataSetNames;
%% Organize the data

% If necessary, create the directory to hold the organized data
if exist(NewDataDir,'dir') ~= 7
    mkdir(NewDataDir)
end

for iFile = 1:length(dataSetNames)
    
    % Load the results from DeepLabCut
    % !!!! Recall, the RGB videos (used in DLC) are upside down and flipped
    % from labels in the rest of the recording system !!!!
    strSplit = split(dataSetIDs{iFile});
    prefixDLC = [strSplit{1} strSplit{3} '_video'];
    dlcFileName = [path dataSetIDs{iFile} '\' prefixDLC affix];

    M = readmatrix(dlcFileName);
    C = readcell(dlcFileName,'Range',[1,1,3,size(M,2)]);
    frame = M(:,1) + 1;
    SnoutTip = M(:,2:4);
    HeadCap = M(:,5:7);
    MidBack = M(:,8:10);
    BackTail = M(:,11:13);
    lightLDLC = M(:,16); % We only care about the probabilities
    lightRDLC = M(:,19); % We only care about the probabilities
    sipLDLC = M(:,20:22);
    sipRDLC = M(:,23:25);
    ULC = M(:,26:28);
    LLC = M(:,29:31);
    URC = M(:,32:34);
    LRC = M(:,35:37);
    
    % Replace outliers in SnoutTip, HeadCap, MidBack, and BackTail data
    SnoutTip(:,1) = filloutliers(SnoutTip(:,1),'linear','movmedian',10);
    SnoutTip(:,2) = filloutliers(SnoutTip(:,2),'linear','movmedian',10);
    HeadCap(:,1) = filloutliers(HeadCap(:,1),'linear','movmedian',10);
    HeadCap(:,2) = filloutliers(HeadCap(:,2),'linear','movmedian',10);
    MidBack(:,1) = filloutliers(MidBack(:,1),'linear','movmedian',10);
    MidBack(:,2) = filloutliers(MidBack(:,2),'linear','movmedian',10);
    BackTail(:,1) = filloutliers(BackTail(:,1),'linear','movmedian',10);
    BackTail(:,2) = filloutliers(BackTail(:,2),'linear','movmedian',10);
    
    % Organize body coordinates (column 1: x coord, column 2: y coord,
    % column 3: DLC probability)
    bodyCoords = {SnoutTip,HeadCap,MidBack,BackTail};
    
    % Find the coordinates of the sippers and corners
    sipLDLC(sipLDLC(:,3) < 0.9,:) = [];
    sipLDLC = mean(sipLDLC(:,1:2));
    sipRDLC(sipRDLC(:,3) < 0.9,:) = [];
    sipRDLC = mean(sipRDLC(:,1:2));
    ULC(ULC(:,3) < 0.9,:) = [];
    ULC = mean(ULC(:,1:2));
    LLC(LLC(:,3) < 0.9,:) = [];
    LLC = mean(LLC(:,1:2));
    URC(URC(:,3) < 0.9,:) = [];
    URC = mean(URC(:,1:2));
    LRC(LRC(:,3) < 0.9,:) = [];
    LRC = mean(LRC(:,1:2));
    
    % Organize box landmark coordinates
    % Note, these are from the RGB video which is flipped from the standard
    % labeling 
    % left sipper, right sipper, upper left corner, lower left corner, upper right corner, lower right corner
    boxCoords = [sipLDLC;sipRDLC;ULC;LLC;URC;LRC]; 
    
    % Load the med associates events
    load([dataSetNames{iFile},filesep,'maEvents.mat'],'maEvents','maTimestamps')
    
    % Reorganize the med associates events !!!! Double Check Assignments !!!!
    lightL = maTimestamps(maEvents == 0);
    lightL = reshape(lightL,[2,length(lightL)/2]); % Top row is on, bottom row is off
    lightR = maTimestamps(maEvents == 1);
    lightR = reshape(lightR,[2,length(lightR)/2]);
    sipL = maTimestamps(maEvents == 2);
    sipL = reshape(sipL,[2,length(sipL)/2]); % Top row is insertion, bottom row is retraction
    sipR = maTimestamps(maEvents == 3);
    sipR = reshape(sipR,[2,length(sipR)/2]);
    
    % Find closest match of med associates light signals to DLC light
    % signals to get a first pass linear transformation from frames to
    % openEphys time
    temp = lightL - lightL(1);
    temp = round((temp*30) + 1);
    maOn = zeros([1,temp(end)]);
    for iLight = 1:size(temp,2)
        maOn(temp(1,iLight):temp(2,iLight)) = 1;
    end
    simScore = zeros([1,length(lightRDLC) - length(maOn)]);
    for i = 1:(length(lightRDLC) - length(maOn))
        simScore(i) = sum(lightRDLC(i:(i + length(maOn) - 1))'.*maOn);
    end
    frameIR = find(simScore == max(simScore));
    P0 = [1/30,lightL(1,1) - (frameIR*(1/30))]; % Coefficients for first pass linear fit
    
    % Find light on frames closest to expected light on frames based on
    % first pass linear fit and med associates recorded times
    lightLDLC(lightLDLC > 0.95) = 1;
    lightLDLC(lightLDLC < 1) = 0;
    temp = find([0,diff(lightLDLC)'] == 1);
    lightLDLCon = NaN([1,size(lightR,2)]);
    for iOn = 1:size(lightR,2)
        predFrame = (lightR(1,iOn) - P0(2))/P0(1);
        lightLDLCon(iOn) = temp(find(abs(temp - predFrame) == min(abs(temp - predFrame)),1,'first'));
    end
    lightRDLC(lightRDLC > 0.95) = 1;
    lightRDLC(lightRDLC < 1) = 0;
    temp = find([0,diff(lightRDLC)'] == 1);
    lightRDLCon = NaN([1,size(lightL,2)]);
    for iOn = 1:size(lightL,2)
        predFrame = (lightL(1,iOn) - P0(2))/P0(1);
        lightRDLCon(iOn) = temp(find(abs(temp - predFrame) == min(abs(temp - predFrame)),1,'first'));
    end
    
    % Refine the linear fit with the updated light on frames
    if ~isequal(size(lightRDLCon),size(lightL(1,:))) || ~isequal(size(lightLDLCon),size(lightR(1,:)))
        error('Different number of light on/off times')
    end
    x = [lightRDLCon';lightLDLCon']; 
    y = [lightL(1,:)';lightR(1,:)'];
    P = polyfit(x,y,1);
    frameT = P(1)*frame + P(2);
    
    % Get the med associates task type from the sessionInfo spreadsheet
    [~,~,sessionInfo] = xlsread([dataDirs{iFile},filesep,'sessionInfo.xlsx']);
    if strcmp(sessionInfo{8,2},'TimmeVisual2CapVer3A')
        % Solid CS+, Blinking CS-
        CSpL = lightL(1,(lightL(2,:) - lightL(1,:)) > 3);
        CSpR = lightR(1,(lightR(2,:) - lightR(1,:)) > 3);
        CSn = lightL(1,(lightL(2,:) - lightL(1,:)) < 3);
        CSn = CSn(1:4:end);
        % Error Check
        if nnz(ismember(CSpL,CSpR)) > 0
            error('Task name is incorrect in sessionInfo')
        end
        if length(CSpR) ~= length(CSpL)
            error('Unequal numbers of CS positive trials.')
        end
        if length(CSn) ~= (length(CSpL) + length(CSpR))
            error('CS negative trial count error.')
        end
    elseif strcmp(sessionInfo{8,2},'TimmeVisual2CapVer3B')
        % Blinking CS+, Solid CS-
        CSpL = lightL(1,(lightL(2,:) - lightL(1,:)) < 3);
        CSpL = CSpL(1:4:end);
        CSpR = lightR(1,(lightR(2,:) - lightR(1,:)) < 3);
        CSpR = CSpR(1:4:end);
        CSn = [];
        iLight = 1;
        while iLight <= size(lightL,2)
            if (lightL(2,iLight) - lightL(1,iLight)) < 3
                iLight = iLight + 1;
            else
                nCSs = round((lightL(2,iLight) - lightL(1,iLight))/4);
                CSn = [CSn,lightL(1,iLight) + 4*(0:(nCSs - 1))];
                iLight = iLight + 1;
            end
        end
        % Error Check
        if nnz(ismember(CSpL,CSpR)) > 0
            error('Task name is incorrect in sessionInfo')
        end
        if length(CSpR) ~= length(CSpL)
            error('Unequal numbers of CS positive trials.')
        end
        if length(CSn) ~= (length(CSpL) + length(CSpR))
            error('CS negative trial count error.')
        end
    end
    
    % Make trial time master list (in open ephys time)
    trialTimes = [CSpL';CSpR';CSn'];
    trialTypes = [1*ones(size(CSpL'));2*ones(size(CSpR'));3*ones(size(CSn'))];
    nTrials = length(trialTimes);
    
    % Calculate the minimum distance between the sippers and the snout tip at
    % several key points (access on CS+, stimulus on CS+, stimulus on CS-)
    approach = NaN([nTrials,4,2]);
    disThresh = 9; % Pixel threshold for approaches based on comparison to manually coded data sets
    disRun = 3; % Required number of frames under the pixel threshold to address momentary jumps
    for iTrial = 1:nTrials
        if trialTypes(iTrial) == 1
            accessTrack = SnoutTip((frameT > (trialTimes(iTrial) + 5.5)) & (frameT < (trialTimes(iTrial) + 13.5)),1:2);
            accessTime = frameT((frameT > (trialTimes(iTrial) + 5.5)) & (frameT < (trialTimes(iTrial) + 13.5)));
            stimTrack = SnoutTip((frameT > (trialTimes(iTrial) + 0)) & (frameT < (trialTimes(iTrial) + 4)),1:2);
            stimTime = frameT((frameT > (trialTimes(iTrial) + 0)) & (frameT < (trialTimes(iTrial) + 4)));
            correctSip = sipRDLC; % Recall flip from RGB in DLC
            incorrectSip1 = sipLDLC;
            incorrectSip2 = [];
        elseif trialTypes(iTrial) == 2
            accessTrack = SnoutTip((frameT > (trialTimes(iTrial) + 5.5)) & (frameT < (trialTimes(iTrial) + 13.5)),1:2);
            accessTime = frameT((frameT > (trialTimes(iTrial) + 5.5)) & (frameT < (trialTimes(iTrial) + 13.5)));
            stimTrack = SnoutTip((frameT > (trialTimes(iTrial) + 0)) & (frameT < (trialTimes(iTrial) + 4)),1:2);
            stimTime = frameT((frameT > (trialTimes(iTrial) + 0)) & (frameT < (trialTimes(iTrial) + 4)));
            correctSip = sipLDLC; % Recall flip from RGB in DLC
            incorrectSip1 = sipRDLC;
            incorrectSip2 = [];
        elseif trialTypes(iTrial) == 3
            accessTrack = [];
            accessTime = [];
            stimTrack = SnoutTip((frameT > (trialTimes(iTrial) + 0)) & (frameT < (trialTimes(iTrial) + 4)),1:2);
            stimTime = frameT((frameT > (trialTimes(iTrial) + 0)) & (frameT < (trialTimes(iTrial) + 4)));
            correctSip = []; % Recall flip from RGB in DLC
            incorrectSip1 = sipLDLC;
            incorrectSip2 = sipRDLC;
        end
        
        % Access Time
        if ~isempty(accessTrack)
            dis = sqrt(sum((accessTrack - repmat(correctSip,[size(accessTrack,1),1])).^2,2))';
            close = (dis < disThresh);
            hits = strfind(close,ones([1,disRun]));
            if ~isempty(hits)
                approach(iTrial,1,1) = 1; % Correct approach during access
                approach(iTrial,1,2) = accessTime(hits(1)) - trialTimes(iTrial); % Time for correct approach during access
            else
                approach(iTrial,1,1) = 0;
            end
%             dis(iTrial,1) = min(sqrt(sum((accessTrack - repmat(correctSip,[size(accessTrack,1),1])).^2,2)));
            dis = sqrt(sum((accessTrack - repmat(incorrectSip1,[size(accessTrack,1),1])).^2,2))';
            close = (dis < disThresh);
            hits = strfind(close,ones([1,disRun]));
            if ~isempty(hits)
                approach(iTrial,2,1) = 1; % Incorrect approach during access
                approach(iTrial,2,2) = accessTime(hits(1)) - trialTimes(iTrial); % Time for incorrect approach during access
            else
                approach(iTrial,2,1) = 0;
            end
%             dis(iTrial,2) = min(sqrt(sum((accessTrack - repmat(incorrectSip1,[size(accessTrack,1),1])).^2,2)));
        end
        
        % Stimulus Time
        if ~isempty(correctSip)
            dis = sqrt(sum((stimTrack - repmat(correctSip,[size(stimTrack,1),1])).^2,2))';
            close = (dis < disThresh);
            hits = strfind(close,ones([1,disRun]));
            if ~isempty(hits)
                approach(iTrial,3,1) = 1; % Correct approach during stimulus
                approach(iTrial,3,2) = stimTime(hits(1)) - trialTimes(iTrial); % Time for correct approach during stimulus
            else
                approach(iTrial,3,1) = 0;
            end
%             dis(iTrial,3) = min(sqrt(sum((stimTrack - repmat(correctSip,[size(stimTrack,1),1])).^2,2)));
        end
        dis = sqrt(sum((stimTrack - repmat(incorrectSip1,[size(stimTrack,1),1])).^2,2))';
        close = (dis < disThresh);
        hits = strfind(close,ones([1,disRun]));
        if ~isempty(hits)
            approach(iTrial,4,1) = 1; % Incorrect approach during stimulus
            approach(iTrial,4,2) = stimTime(hits(1)) - trialTimes(iTrial); % Time for incorrect approach during stimulus
        else
            approach(iTrial,4,1) = 0;
        end
%         dis(iTrial,4) = min(sqrt(sum((stimTrack - repmat(incorrectSip1,[size(stimTrack,1),1])).^2,2)));
        if ~isempty(incorrectSip2)
            dis = sqrt(sum((stimTrack - repmat(incorrectSip2,[size(stimTrack,1),1])).^2,2))';
            close = (dis < disThresh);
            hits = strfind(close,ones([1,disRun]));
            if ~isempty(hits)
                approach(iTrial,4,1) = 1; % Incorrect approach during stimulus (other side for CS-)
                approach(iTrial,4,2) = min([approach(iTrial,4,2),stimTime(hits(1)) - trialTimes(iTrial)]); % Time for incorrect approach during stimulus (other side for CS-, first approach)
            end
%             dis(iTrial,4) = min([dis(iTrial,4),min(sqrt(sum((stimTrack - repmat(incorrectSip2,[size(stimTrack,1),1])).^2,2)))]);
        end
    end
    
    % Record the parameters used for determining approaches
    approachParams = [disThresh,disRun];
    
    % Organize body coordinates (column 1: x coord, column 2: y coord,
    % column 3: DLC probability)
    bodyCoords = {frameT,SnoutTip,HeadCap,MidBack,BackTail};
    
%     % Load other necessary data
%     load([dataDirs{iFile},filesep,'spkData.mat'])
%     load([dataDirs{iFile},filesep,'clusType.mat'])
%     load([dataDirs{iFile},filesep,'shankNum.mat'])
%     
%     % Fix time offset in spike times (open ephys voltage recording offset)
%     [~,timestamps,~] = load_open_ephys_data_faster([dataDirs{iFile},filesep,'100_CH1.continuous']);
%     if max(abs(diff(timestamps) - diff(timestamps(1:2)))) > 10^(-8)
%         disp('Unstable time stamps.')
%     end
%     for iNeuron = 1:length(spkData)
%         spkData{iNeuron} = timestamps(round(spkData{iNeuron}*30000));
%     end
    
    % Save the data set
    save([NewDataDir,filesep,dataSetNames{iFile},'_AllData.mat'],'sessionInfo','trialTimes','trialTypes','approach','approachParams','boxCoords','bodyCoords')
    
    % Display progress
    disp([num2str(100*iFile/length(dataDirs),3),'% Done'])
end
