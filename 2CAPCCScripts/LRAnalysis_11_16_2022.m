%
%% Load and check data
clear all; close all % Refresh workspace
rng(1)
% PATHS
% For brain3 Linux Machine:
if ispc
    parentPath = 'F:/dissDat';
    figPath = 'F:/dissDat/lrfigs';
else
    parentPath = '/research/lapishla/dissDat';
    figPath = '/research/lapishla/dissDat/figs';
end
% LOAD
load([parentPath filesep 'ephysStruct.mat']);
load([parentPath filesep 'trlStruct.mat']);
load([parentPath filesep 'masterTable.mat']);

addpath(genpath(parentPath))

sessionType = 'Incongruent';

%% Set indexing based on Congruent versus Congruent Sessions
regPIdx = startsWith(masterTbl.SessionType,'Regular') & startsWith(masterTbl.Strain,'P');
revPIdx = startsWith(masterTbl.SessionType,'Reversal') & startsWith(masterTbl.Strain,'P');
    
regWIdx = startsWith(masterTbl.SessionType,'Regular') & startsWith(masterTbl.Strain,'W');
revWIdx = startsWith(masterTbl.SessionType,'Reversal') & startsWith(masterTbl.Strain,'W');

%% Timepoints and other parameters
numTrials = 48;


sipDescent = 6;
sipAscent = 14;
cueOn = 0;

LeftTrials = 1:24;
RightTrials = 25:48;
colors = [{'#EDB120'}, {'#7E2F8E'}, {'#0072BD'},  {'#D95319'} ];

%%
% In order to create a PSTH with approach, it will be necessary to create a
% new variable akin to 'trialTimes'.
% This new variable will be trialTimes + time of approach.
% Time of approach is in the 'approach' variable. 
% First, there needs to be a determination of correct/incorrect choices,
% and then based on that determination add the time each correct/incorrect
% choice occurs to the time of 'trialTimes'

%% ApproachTime PSTH creation
% 
rng(1)
numTrials = 13;
PSTH = [];
% Pull ephys and trial data, create PSTH
    pBin = 1/10; % 33 ms bins for fr variable, match video framerate
    PSTH_Approach = [];
    if ~isfield(ephysStruct,'PSTH_Approach')
        
        for i = 1:length(trlStruct)
    
            % Load trial times from trlStruct
            trialTimes = trlStruct(i).trialTimes(1:48);
            [trialTimes, trialTimesIdx] = sort(trialTimes);


%             approachVar = trlStruct(i).approachChangePoints;           
%             sortedApproachVar = approachVar(trialTimesIdx);

            approachVar = trlStruct(i).approach(1:48,1,2);
            approachVar(isnan(approachVar)) = datasample(approachVar(~isnan(approachVar)),1);
            sortedApproachVar = approachVar(trialTimesIdx);

            approachIdx = trlStruct(i).approach(1:48,1,1) == 1 | trlStruct(i).approach(1:48,2,1) == 1;
            approachIdx = approachIdx(trialTimesIdx);
%             sortedApproachVar(approachIdx == 0) = NaN;
%             sortedApproachVar(approachIdx == 0) = randsample(sortedApproachVar(~isnan(sortedApproachVar)),1);
            
%             correctionMidPoints = (( trlStruct(i).approach(1:48,1,2) + trlStruct(i).approach(1:48,2,2) ) / 2) - 5;
%             correctionMidPoints = correctionMidPoints(trialTimesIdx);
% % 
%             maxVar = max(trlStruct(i).approach(1:48,1:2,2),[],2);
%             maxVar = maxVar(trialTimesIdx);
%             maxVar(isnan(maxVar)) = datasample(maxVar(~isnan(maxVar)),1);
%             
%             sortedApproachVar = maxVar;

            if startsWith(masterTbl.SessionType(i),'Regular')            
                correctIndex = trlStruct(i).approach(1:48,1,1) == 1 & trlStruct(i).approach(1:48,2,1) == 0;
                correctionIndex = trlStruct(i).approach(1:48,1,1) == 1 & trlStruct(i).approach(1:48,2,1) == 1;
                incorrectIndex = trlStruct(i).approach(1:48,2,1) == 1;
                
                correctIndexL = trlStruct(i).approach(1:24,1,1) == 1 & trlStruct(i).approach(1:24,2,1) == 0;
                correctIndexR = trlStruct(i).approach(25:48,1,1) == 1 & trlStruct(i).approach(25:48,2,1) == 0;
                correctionL = trlStruct(i).approach(1:24,1,2) > trlStruct(i).approach(1:24,2,2); 
                correctionR = trlStruct(i).approach(25:48,1,2) > trlStruct(i).approach(25:48,2,2);
             
            elseif startsWith(masterTbl.SessionType(i),'Reversal')
                correctIndex = trlStruct(i).approach(1:48,1,1) == 0 & trlStruct(i).approach(1:48,2,1) == 1;
                correctionIndex = trlStruct(i).approach(1:48,1,1) == 1 & trlStruct(i).approach(1:48,2,1) == 1;

                incorrectIndex = trlStruct(i).approach(1:48,1,1) == 1;
                
                correctIndexL = trlStruct(i).approach(25:48,1,1) == 0 & trlStruct(i).approach(25:48,2,1) == 1;
                correctIndexR = trlStruct(i).approach(1:24,1,1) == 0 & trlStruct(i).approach(1:24,2,1) == 1;

            end

            approachIndex = trlStruct(i).approach(1:48,1,1) == 1 | trlStruct(i).approach(1:48,2,1) == 1;
            noApproachIndex = approachIndex(approachIndex == 0);

            correctIndex = correctIndex(trialTimesIdx);
            incorrectIndex = incorrectIndex(trialTimesIdx);
            correctionIndex = correctionIndex(trialTimesIdx);

            LRIdx = double([correctIndexL; correctIndexR]);
            LRIdx(1:24,2) = 1;
            LRIdx = LRIdx(trialTimesIdx,:);
            correctIndexL = (LRIdx(:,1) == 1 & LRIdx(:,2) == 1);
            correctIndexR = (LRIdx(:,1) == 1 & LRIdx(:,2) == 0);
            
            correctionLRIdx = double([correctionL; correctionR]);
            correctionLRIdx(1:24,2) = 1;
            correctionLRIdx = correctionLRIdx(trialTimesIdx,:);
            correctionIndexL = (correctionLRIdx(:,1) == 1 & correctionLRIdx(:,2) == 1);
            correctionIndexR = (correctionLRIdx(:,1) == 1 & correctionLRIdx(:,2) == 0);
            
%             sortedApproachVar(correctionIndex) = correctionMidPoints(correctionIndex);
            % Load ephys related data
            stMtx = ephysStruct(i).stmtx;
            if ~isempty(stMtx)
                stMtx = bootISI(stMtx,0.05,0.003);  % Custom function to boot neurons where ISI Violations are a greater than 5% occurence 
            end
            % 
            % Create fr based on pBin above
%             h=histc(stMtx,[min(stMtx(:)):pBin:max(stMtx(:))]);
            fr = gaussconv(stMtx,pBin);  % Gaussconv function is preferred
    
            % Find time, in seconds, that each fr index corresponds to.
            frTime = (1:length(fr))*pBin;   % Denominator corresponds with pBin
            if ~isempty(frTime)
                if frTime(end) < trialTimes(end)
                    zeroPad = frTime(end):pBin:trialTimes(end)+20;
                    fr = [fr; zeros(length(zeroPad),size(fr,2))];
                    frTime = [frTime zeroPad];
                end
            end
            dlcTimestamps = trlStruct(i).bodyCoords{1};
            dlcTimestamps = dlcTimestamps(dlcTimestamps >= 0);
            dlcTimestamps = (dlcTimestamps);
            trialTimes = trialTimes + sortedApproachVar';
        
        for k = 1:length(trialTimes(1:48))
            [~,tMin] = min(abs((trialTimes(k) - 2) - dlcTimestamps));
            [~,tMax] = min(abs((trialTimes(k) + 2) - dlcTimestamps)); 

            [~,tMin_fr] = min(abs((dlcTimestamps(tMin)) - frTime));
            [~,tMax_fr] = min(abs((dlcTimestamps(tMax)) - frTime));
            
            trlSize(k) = size(fr(tMin_fr:tMax_fr,:),1);

            if size(fr(tMin_fr:tMax_fr,:),1) == 41
                PSTH(k,:,:) = fr(tMin_fr:tMax_fr,:);
            elseif size(fr(tMin_fr:tMax_fr,:),1) == 42
                PSTH(k,:,:) = fr(tMin_fr:tMax_fr-1,:);
            elseif size(fr(tMin_fr:tMax_fr,:),1) == 40
                PSTH(k,:,:) = fr(tMin_fr:tMax_fr+1,:);          
            end
            
                
               
        end

        

        % Remove low firing rate neurons. Not sure what the best value is
        % here.
        minFR = 0.01;
        if ~isempty(PSTH)            
            meanFR = squeeze(mean(mean(PSTH,1),2)); 
            lowFR = meanFR > minFR; % Find values greater than minFR
            PSTH = PSTH(:,:,lowFR);
            if length(size(PSTH)) == 2
                S = size(PSTH);

                PSTH = reshape(PSTH,S(1),1,S(2));

            end
        end
        if ~isempty(PSTH)
            ephysStruct(i).PSTH_Approach = PSTH(approachIdx,:,:);
            ephysStruct(i).PSTH_NoApproach = PSTH(approachIdx == 0,:,:);
            ephysStruct(i).PSTH = PSTH;
        else 
            ephysStruct(i).PSTH_Approach = [];
        end
        if ~isempty(PSTH)
            ephysStruct(i).neuronNumApproach = size(PSTH,3);
        else
            ephysStruct(i).neuronNumApproach = 0;
        end
        ephysStruct(i).maxFR = max(fr);
        
        if ~isempty(PSTH)
            
            ephysStruct(i).PSTH_Correct = PSTH(correctIndex(1:15),:,:);
            ephysStruct(i).PSTH_Incorrect = PSTH(incorrectIndex(1:15),:,:);

            ephysStruct(i).PSTH_Left = PSTH(correctIndexL,:,:);
            ephysStruct(i).PSTH_Right = PSTH(correctIndexR,:,:);
        
%             ephysStruct(i).PSTH_Approach(:,:,:);
%             ephysStruct(i).PSTH_NoApproach = [];
   
        else

            ephysStruct(i).PSTH_Correct = [];
            ephysStruct(i).PSTH_Incorrect = [];

        end

        ephysStruct(i).PSTH_ApproachMean = squeeze(mean(ephysStruct(i).PSTH_Approach,1));
        ephysStruct(i).PSTH_NoApproachMean = squeeze(mean(ephysStruct(i).PSTH_NoApproach,1));

        ephysStruct(i).PSTH_ApproachMean_Trial = squeeze(mean(ephysStruct(i).PSTH_Approach,3));
        ephysStruct(i).PSTH_NoApproachMean_Trial = squeeze(mean(ephysStruct(i).PSTH_NoApproach,3));

        ephysStruct(i).NumApproaches = size(ephysStruct(i).PSTH_ApproachMean_Trial,1);
        ephysStruct(i).NumNoApproaches = size(ephysStruct(i).PSTH_NoApproachMean_Trial,1);

        ephysStruct(i).correctIndex = correctIndex;
        ephysStruct(i).incorrectIndex = incorrectIndex;

        ephysStruct(i).PSTH_CorrectMean = squeeze(mean(ephysStruct(i).PSTH_Correct,1));
        ephysStruct(i).PSTH_IncorrectMean = squeeze(mean(ephysStruct(i).PSTH_Incorrect,1));

        ephysStruct(i).LeftIndex = correctIndexL;
        ephysStruct(i).RightIndex = correctIndexR;
        
        ephysStruct(i).correctionLeftIndex = correctionIndexL;
        ephysStruct(i).correctionRightIndex = correctionIndexR;

        ephysStruct(i).PSTH_LeftMean= squeeze(mean(ephysStruct(i).PSTH_Left,1));
        ephysStruct(i).PSTH_RightMean = squeeze(mean(ephysStruct(i).PSTH_Right,1));

        ephysStruct(i).correctionsIndex = correctionIndex;


        PSTH_Approach = [];
        PSTH_NoApproach = [];
        PSTH = [];

        

    end

end

%%
for i = 1:length(ephysStruct)
    if size(ephysStruct(i).PSTH_ApproachMean,1) == 1
       ephysStruct(i).PSTH_ApproachMean = ephysStruct(i).PSTH_ApproachMean';
       ephysStruct(i).PSTH_NoApproachMean = ephysStruct(i).PSTH_NoApproachMean';
       ephysStruct(i).PSTH_CorrectMean = ephysStruct(i).PSTH_CorrectMean';
       ephysStruct(i).PSTH_IncorrectMean = ephysStruct(i).PSTH_IncorrectMean';
       ephysStruct(i).PSTH_LeftMean= ephysStruct(i).PSTH_LeftMean';
       ephysStruct(i).PSTH_RightMean = ephysStruct(i).PSTH_RightMean';
    end
end
ephysStruct(19).PSTH_LeftMean = [];
ephysStruct(19).PSTH_RightMean = [];
ephysStruct(19).neuronNumApproach = []; 
%% Based on data above, there are no real differences within Genotype (WCon,WInc OR PCon,PInc) coefficient loadings. However, there are consistent differences between genotype loadings
% Therefore, let us analyze P and Wistars separate from one another. 
% This will allow within-session type comparisons but not necessarily allow
% between genotype comparisons except in certain circumstances such as
% comparing the distance metrics between two PC subspaces. 

% Depending on the session types we would like to analyze, we can either
% take the index of all P and W rats or we can take the index of ONLY P / W
% rats in Congruent Sessions
if strcmp(sessionType,'Congruent')
    PIdx = startsWith(masterTbl.Strain,'P') & startsWith(masterTbl.SessionType,'Regular');
    WIdx = startsWith(masterTbl.Strain,'W') & startsWith(masterTbl.SessionType,'Regular');
elseif strcmp(sessionType,'Incongruent')
    PIdx = startsWith(masterTbl.Strain,'P') & startsWith(masterTbl.SessionType,'Reversal');
    WIdx = startsWith(masterTbl.Strain,'W') & startsWith(masterTbl.SessionType,'Reversal'); 
    PIdx(19) = 0; WIdx(19) = 0;
end
% PIdx = startsWith(masterTbl.Strain,'P');
% WIdx = startsWith(masterTbl.Strain,'W');

% A lot of this will repeat code from above. 
% Set indexing based on Genotype
pMat = [];
wMat = [];
%


% %% Index Wistar, P
% % Congruent
pMatLC = [ephysStruct(regPIdx).PSTH_LeftMean];
wMatLC = [ephysStruct(regWIdx).PSTH_LeftMean];

pMatRC = [ephysStruct(regPIdx).PSTH_RightMean];
wMatRC = [ephysStruct(regWIdx).PSTH_RightMean];
% Congruent
pMatLIC = [ephysStruct(revPIdx).PSTH_LeftMean];
wMatLIC = [ephysStruct(revWIdx).PSTH_LeftMean];

pMatRIC = [ephysStruct(revPIdx).PSTH_RightMean];
wMatRIC = [ephysStruct(revWIdx).PSTH_RightMean];

conditionIndexes = [length(wMatLC) length(wMatLIC) length(pMatLC) length(pMatLIC)];
conditionLabels = [{'Congruent W'}, {'Congruent W'}, {'Congruent P'}, {'Congruent P'}];


%% Reshape Matrices 
% For the purposes of the coefficient matrix:
% allMat contains all Left choices followed in time by a matrix of all of
% the Right choices
% Each individual matrix itself repeats the following pattern:
% Wistar Congruent, Wistar Congruent, P Congruent, P Congruent
% This will be important to pull out the coefficients pertaining to each
% trial type 
% From there, we can take the coefficients of known datasets and determine
% if there are differences between the L/R
% 
% 08/22/2022 - - First, determine if there are baseline differences between
% L and R choices in the Congruent sessions prior to the addition of
% Congruent choices.
% In order to do this, we just take the following modification of the above
% pattern:
% Congruent Wistar Left Choices, Congruent P Rat Left Choices, 

% Wistar Right Choices, Congruent P Rats Right Choices.

allMat = [];

% Below contain ALL session types
% allMatL = [wMatLC wMatLIC pMatLC pMatLIC];
% allMatR = [wMatRC wMatRIC pMatRC pMatRIC];

% Below contain Congruent only session types
if strcmp(sessionType,'Congruent')

    allMatL = [wMatLC pMatLC];
    allMatR = [wMatRC pMatRC];
    conditionIndexes = [length(wMatLC) length(pMatLC)];
    conditionLabels = [{'Congruent W'}, {'Congruent P'}];
    % Format matrix such that it expands 'in time / observations'
    allMat = [allMatL; allMatR];

elseif strcmp(sessionType,'Incongruent')

    allMatL = [wMatLIC pMatLIC];
    allMatR = [wMatRIC pMatRIC];
    conditionIndexes = [length(wMatLIC) length(pMatLIC)];
    conditionLabels = [{'Incongruent W'}, {'Incongruent P'}];
    % Format matrix such that it expands 'in time / observations'
    allMat = [allMatL; allMatR];

end


%% Regress allMat

%% PCA Preprocessing | Separation between P, W groups | LEFT TRIALS



[ConditionsCoeff,ConditionsScore,ConditionsLatent,~,ConditionsExplained] = pca(allMat);

% Need to find an index of neurons per dataset, be able to sort that index
% into our conditions (W/P; Con/Inc). Should ideally be a 1x54 matrix
% containing N # of neurons.

PneuronIndex = [ephysStruct(PIdx).neuronNumApproach];
WneuronIndex = [ephysStruct(WIdx).neuronNumApproach];

%% Basic Figures First
% Plot mean of allMat
% figure('Units','normalized','Position',[0 0 1 1])
% plot(mean(allMat(1:71,:),2))
% hold on
% plot(mean(allMat(72:end,:),2))

%% Plot PCs over Time without Separations
timeVec = [-2:0.1:2];


figure('Units','normalized','Position',[0 0 1 1])
% Left
plot(timeVec,ConditionsScore(1:41,1),'b-','LineWidth',3)
hold on
plot(timeVec,ConditionsScore(1:41,2),'r-','LineWidth',3)
plot(timeVec,ConditionsScore(1:41,3),'k-','LineWidth',3)
% Right
plot(timeVec,ConditionsScore(42:end,1),'b--','LineWidth',3)
plot(timeVec,ConditionsScore(42:end,2),'r--','LineWidth',3)
plot(timeVec,ConditionsScore(42:end,3),'k--','LineWidth',3)

legend({'PC1 Left','PC2 Left','PC3 Left','PC1 Right','PC2 Right','PC3 Right'},'AutoUpdate','off')
xlabel('Time Bins')
ylabel('PC Score')
xline(0,'--','Approach','LineWidth',5)
title('Top 3 PCs')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

saveas(gca,[figPath filesep 'top3PCs_Time'     '_100ms_' sessionType ],'svg');     saveas(gca,[figPath filesep 'top3PCs_Time'    '_100ms_' sessionType ],'png')

%% Variance explained
p = length(ConditionsExplained);
pVec = 1:p;
pExpected = zeros(length(pVec),1);

for i = pVec
    pExpected(i) = sum( 1 ./ pVec(i:end)) / p;
end
pExpected = pExpected .* 100; % Convert to percentage
stickPoint = find(pExpected >= ConditionsExplained,1);
figure('Units','normalized','Position',[0 0 1 1])
plot((ConditionsExplained),'ko','MarkerSize',5);
hold on
plot(pExpected,'r--','LineWidth',2)
xline(stickPoint,'k--','Broken Stick Point','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
xlabel('PC Number')
ylabel('Cumulative Explained Variance')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,[figPath filesep 'varianceExplained'     '_100ms_' sessionType],'svg');     saveas(gca,[figPath filesep 'varianceExplained'    '_100ms_' sessionType],'png')

%% Pull Individual Scores
% Calculate scores for reach 
allCond = 0;
% First determine boundaries for each index 
if allCond == 1
    wConIdx = 1:conditionIndexes(1);
    wIncIdx = conditionIndexes(1)+1:conditionIndexes(1)+conditionIndexes(2);
    pConIdx = conditionIndexes(1)+conditionIndexes(2)+1:conditionIndexes(1)+conditionIndexes(2)+conditionIndexes(3);
    pIncIdx = conditionIndexes(1)+conditionIndexes(2)+conditionIndexes(3)+1:conditionIndexes(1)+conditionIndexes(2)+conditionIndexes(3)+conditionIndexes(4);
else
    wConIdx = 1:conditionIndexes(1);
    pConIdx = conditionIndexes(1)+1:conditionIndexes(1)+conditionIndexes(2);
end

%% Generate scores using indices 
% Will require the mean centered data 
rawData = ConditionsCoeff*ConditionsScore';
rawData = allMat';

wConScore = rawData(wConIdx,:)'*ConditionsCoeff(wConIdx,:);
% wIncScore = rawData(wIncIdx,:)'*ConditionsCoeff(wIncIdx,:);

pConScore = rawData(pConIdx,:)'*ConditionsCoeff(pConIdx,:);
% pIncScore = rawData(pIncIdx,:)'*ConditionsCoeff(pIncIdx,:);

%%
LeftTime = 1:41;
RightTime = 42:82;

wConL = wConScore(LeftTime,1:stickPoint);
wConR = wConScore(RightTime,1:stickPoint);

pConL = pConScore(LeftTime,1:stickPoint);
pConR = pConScore(RightTime,1:stickPoint);

[test,testvol] = convhulln(wConR);

%% Project PCs onto Trials 
% One way of doing this is to take the coefficients from each session and
% project those onto the full PSTH for that session. From there, we can
% take the average of all trial types for that sessions and concatenate
% those into one variable per genotype / session type
%
%
%  First statement is to make sure our indices align!
if sum(PneuronIndex) + sum(WneuronIndex) == size(ConditionsCoeff,1)
    disp('Sizes of individual indices match total coeff')
else
    error('Indices mismatch')
end

indexedCoeff = [];
% Make sure that W comes before P since this is how we have it set up above
bothIndices = [WneuronIndex PneuronIndex];
indexHelper = zeros(1,length(bothIndices)+1);
indexHelper(2:end) = cumsum(bothIndices);

for i = 1:length(bothIndices)
    indexedCoeff{i} = ConditionsCoeff(indexHelper(i)+1:indexHelper(i+1),:);    
end

% Create a 4D PSTH Housing All Relevant Data
% Find position of each dataset in the order presented
dataPos = [find(WIdx == 1); find(PIdx == 1)];
% Monstrous For Loop Structure to Find Trial Projections
% i = dataset, k = trial, j = PC#
for i = 1:length(dataPos)
    numPC = 1:length(ConditionsExplained);
    PSTH = ephysStruct(dataPos(i)).PSTH;
    % Mean Center PSTH
    PSTH = (mean(PSTH,3) - PSTH);
    if ~isempty(PSTH)
        for k = 1:size(PSTH,1)
            for j = 1:length(numPC)
                trlProj{i,j}(k,:,:) = squeeze(PSTH(k,:,:)) * indexedCoeff{i}(:,j);
            end
        end    
    end
end
% This leaves us with a Number of Session x Number of PC Cell Array
% Cell Aray has 48 Trials x 71 Observations
% Now, we can index on the cell aray to obtain L/R Trial Types

for i = 1:length(dataPos)
    correctL = ephysStruct(dataPos(i)).LeftIndex;
    correctR = ephysStruct(dataPos(i)).RightIndex;
    correctionsL = ephysStruct(dataPos(i)).correctionLeftIndex;
    correctionsR = ephysStruct(dataPos(i)).correctionRightIndex;
    if ~isempty(trlProj{i})
        for j = 1:length(numPC)
            leftProj{i,j} = trlProj{i,j}(correctL,:,:);
            rightProj{i,j} = trlProj{i,j}(correctR,:,:);
            correctionsLProj{i,j} = trlProj{i,j}(correctionsL,:,:);
            correctionsRProj{i,j} = trlProj{i,j}(correctionsR,:,:);
        end
    end
end

% 2 Cell Arrays now exist with Left and Right Projections
% Can now take the mean for each in each genotype and each PC#

for j = 1:length(numPC)
    wLProj{j} = vertcat(leftProj{1:13,j}); 
    wRProj{j} = vertcat(rightProj{1:13,j}); 
    wLCrct{j} = vertcat(correctionsLProj{1:13,j});
    wRCrct{j} = vertcat(correctionsRProj{1:13,j});
    

    pLProj{j} = vertcat(leftProj{14:end,j});
    pRProj{j} = vertcat(rightProj{14:end,j});
    pLCrct{j} = vertcat(correctionsLProj{14:end,j});
    pRCrct{j} = vertcat(correctionsRProj{14:end,j});

end
% %% Clear some big variables
% ephysStruct = [];
% trlStruct = [];
% masterTbl = [];

%% Calculate volume, distance per trial in an ascending manner.
% I.e. 2PCs, 3PCs, 4PCs, usw
% Calculate volume, distance per trial and then average these 
% Variables of interest are wLProj, wRProj, pLProj pRProj
% Distances can be between Left and Right
% Volumes will be compared between Left and Right later
% Therefore volumes will be 2x2 (LR PW)
% Distances 1x2 (PW)

numPC = stickPoint;

% Init matrices so they do not stack when rerunning 
wLData = []; wRData = []; pLData = []; pRData = [];
wLCData = []; wRCData = []; pLCData = []; pRCData = [];


% Convert all cell variables into 3D Matrices
wLData = cat(3,wLProj{:});
wRData = cat(3,wRProj{:});
wLCData = cat(3,wLCrct{:});
wRCData = cat(3,wRCrct{:});

pLData = cat(3,pLProj{:});
pRData = cat(3,pRProj{:});
pLCData = cat(3,pLCrct{:});
pRCData = cat(3,pRCrct{:});

wLVol = [];     pLVol = [];
wRVol = [];     pRVol = [];


for k = 1:size(wLData,1)    % Trials
    wLData_helper = [];
    for i = 1:numPC         % Must start at 2 
        % Generate trial-by-trial data structure containing N PCs
        wLData_helper = [wLData_helper; squeeze(wLData(k,:,i))];
        if i ~= 1
            try
                [~,wLVol(i,k)] = convhulln(zscore(wLData_helper'),{'Qt'});
            catch
                warning('Not enough dimensionality')
                wLVol(i,k) = NaN;
            end
        end
    end
end

for k = 1:size(wRData,1)    % Trials
    wRData_helper = [];
    for i = 1:numPC         % Must start at 2 
        % Generate trial-by-trial data structure containing N PCs
        wRData_helper = [wRData_helper; squeeze(wRData(k,:,i))];
        if i ~= 1
            try
                [~,wRVol(i,k)] = convhulln(zscore(wRData_helper'),{'Qt'});
            catch
                warning('Not enough dimensionality')
                wRVol(i,k) = NaN;
            end
        end
    end
end

for k = 1:size(pLData,1)    % Trials
    pLData_helper = [];
    for i = 1:numPC         % Must start at 2 
        % Generate trial-by-trial data structure containing N PCs
        pLData_helper = [pLData_helper; squeeze(pLData(k,:,i))];
        if i ~= 1
            try
                [~,pLVol(i,k)] = convhulln(zscore(pLData_helper'),{'Qt'});
            catch
                warning('Not enough dimensionality')
                pLVol(i,k) = NaN;
            end
        end
    end
end

for k = 1:size(pRData,1)    % Trials
    pRData_helper = [];
    for i = 1:numPC         % Must start at 2 
        % Generate trial-by-trial data structure containing N PCs
        pRData_helper = [pRData_helper; squeeze(pRData(k,:,i))];
        if i ~= 1
            try
                [~,pRVol(i,k)] = convhulln(zscore(pRData_helper'),{'Qt'});
            catch
                warning('Not enough dimensionality')
                pRVol(i,k) = NaN;
            end
        end
    end
end

% Plotting of volumes
% Obtain mean, SEM
wLVolM = nanmean(wLVol,2);     wLVolSEM = nanstd(wLVol,[],2)/sqrt(length(wLVol));
wRVolM = nanmean(wRVol,2);     wRVolSEM = nanstd(wRVol,[],2)/sqrt(length(wRVol));

pLVolM = nanmean(pLVol,2);     pLVolSEM = nanstd(pLVol,[],2)/sqrt(length(pLVol));
pRVolM = nanmean(pRVol,2);     pRVolSEM = nanstd(pRVol,[],2)/sqrt(length(pRVol));

% Plotting
figure('Units','normalized','Position',[0 0 1 1])
plot(wLVolM,'b-o','LineWidth',3); hold on; plot(wRVolM,'b--o','LineWidth',3)
plot(pLVolM,'r-o','LineWidth',3); plot(pRVolM,'r--o','LineWidth',3)
legend([{'WL','WR','PL','PR'}],'AutoUpdate','off')

er = errorbar(wLVolM, wLVolSEM);
er.Color = [0 0 0];
er.LineStyle = 'none';
er.LineWidth = 2;

er = errorbar(wRVolM, wRVolSEM);
er.Color = [0 0 0];
er.LineStyle = 'none';
er.LineWidth = 2;

er = errorbar(pLVolM, pLVolSEM);
er.Color = [0 0 0];
er.LineStyle = 'none';
er.LineWidth = 2;

er = errorbar(pRVolM, pRVolSEM);
er.Color = [0 0 0];
er.LineStyle = 'none';
er.LineWidth = 2;

xlabel('PC Number')
ylabel('Volume of State Space')
title('Relationship between Volume and PC Space')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,[figPath filesep 'volume_v_pcspaces'     '_100ms_' sessionType],'svg');     saveas(gca,[figPath filesep 'volume_v_pcspaces'    '_100ms_' sessionType],'png')

%% Calculate volume as a fcn of permuted 3D spaces
% For whatever reason, adding dimensionality reduced the total volume found
% From now, we can shuffle through some amount of dimensions and take the
% mean of those to see if there are meaningful differences between LR PW
numDim = stickPoint;
allPerms = nchoosek([1:numDim],3);

% Init variables
wLVol3D = []; wRVol3D = [];     wLCVol3D = []; wRCVol3D = [];
pLVol3D = []; pRVol3D = [];     pLCVol3D = []; pRCVol3D = [];

for i = 1:length(allPerms)

    wLData_helper = [squeeze(mean(wLData(:,:,allPerms(i,1)))); squeeze(mean(wLData(:,:,allPerms(i,2)))); squeeze(mean(wLData(:,:,allPerms(i,3))))];  
    wRData_helper = [squeeze(mean(wRData(:,:,allPerms(i,1)))); squeeze(mean(wRData(:,:,allPerms(i,2)))); squeeze(mean(wRData(:,:,allPerms(i,3))))];  
    pLData_helper = [squeeze(mean(pLData(:,:,allPerms(i,1)))); squeeze(mean(pLData(:,:,allPerms(i,2)))); squeeze(mean(pLData(:,:,allPerms(i,3))))];  
    pRData_helper = [squeeze(mean(pRData(:,:,allPerms(i,1)))); squeeze(mean(pRData(:,:,allPerms(i,2)))); squeeze(mean(pRData(:,:,allPerms(i,3))))];

    wLCData_helper = [squeeze(mean(wLCData(:,:,allPerms(i,1)))); squeeze(mean(wLCData(:,:,allPerms(i,2)))); squeeze(mean(wLCData(:,:,allPerms(i,3))))]; 
    wRCData_helper = [squeeze(mean(wRCData(:,:,allPerms(i,1)))); squeeze(mean(wRCData(:,:,allPerms(i,2)))); squeeze(mean(wRCData(:,:,allPerms(i,3))))]; 
    pLCData_helper = [squeeze(mean(pLCData(:,:,allPerms(i,1)))); squeeze(mean(pLCData(:,:,allPerms(i,2)))); squeeze(mean(pLCData(:,:,allPerms(i,3))))]; 
    pRCData_helper = [squeeze(mean(pRCData(:,:,allPerms(i,1)))); squeeze(mean(pRCData(:,:,allPerms(i,2)))); squeeze(mean(pRCData(:,:,allPerms(i,3))))]; 

    [points,wLVol3D(i)] = convhulln(zscore(wLData_helper'),{'Qt'});
    [points2,wRVol3D(i)] = convhulln(zscore(wRData_helper'),{'Qt'});
    [~,pLVol3D(i)] = convhulln(zscore(pLData_helper'),{'Qt'});
    [~,pRVol3D(i)] = convhulln(zscore(pRData_helper'),{'Qt'});

    [~,wLCVol3D(i)] = convhulln(zscore(wLCData_helper'),{'Qt'});
    [~,wRCVol3D(i)] = convhulln(zscore(wRCData_helper'),{'Qt'});
    [~,pLCVol3D(i)] = convhulln(zscore(pLCData_helper'),{'Qt'});
    [~,pRCVol3D(i)] = convhulln(zscore(pRCData_helper'),{'Qt'});
end

% Quick boxplot between the four conditions
data = [wLVol3D; wRVol3D; pLVol3D; pRVol3D];
figure
boxplot(data','notch','on')
ylabel('State Space Volume')
xlabel('Group')
xticklabels({'WL','WR','PL','PR'})
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

% Corrections
data2 = [wLCVol3D; wRCVol3D; pLCVol3D; pRCVol3D];
figure
boxplot(data2','notch','on')
ylabel('State Space Volume')
xlabel('Group')
xticklabels({'WLC','WRC','PLC','PRC'})
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

% Try to plot side-by-side
figure('Units','normalized','Position',[0 0 1 1])
tiledlayout(1,2)
ax1 = nexttile;
boxchart(ax1,data','notch','on')
ylabel(ax1,'State Space Volume')
xticklabels(ax1,{'WL','WR','PL','PR'})
ylim(ax1,[0 30])
set(ax1,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

ax2 = nexttile;
boxchart(ax2,data2','notch','on')
ylabel(ax2,'State Space Volume')
xticklabels(ax2,{'WLC','WRC','PLC','PRC'})
ylim(ax2,[0 30])
set(ax2,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,[figPath filesep '3DVolumes'     '_100ms_' sessionType],'svg');     saveas(gca,[figPath filesep '3DVolumes'    '_100ms_' sessionType],'png')
%%  Plot example of volume space
figure('Units','normalized','Position',[0 0 1 1])
plot3(wLData_helper(1,:),wLData_helper(2,:),wLData_helper(3,:),'k-o','MarkerSize',10,'LineWidth',3);
hold on
plot3(wLData_helper(1,1),wLData_helper(2,1),wLData_helper(3,1),'k^','MarkerSize',20,'LineWidth',3);
plot3(wLData_helper(1,end),wLData_helper(2,end),wLData_helper(3,end),'kv','MarkerSize',20,'LineWidth',3);
trisurf(points,wLData_helper(1,:),wLData_helper(2,:),wLData_helper(3,:),'FaceColor','cyan');
xlabel('PC 2')
ylabel('PC 3')
zlabel('PC 5')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,[figPath filesep 'volumeExample'     '_100ms_' sessionType],'svg');     saveas(gca,[figPath filesep 'volumeExample'    '_100ms_' sessionType],'png')

figure('Units','normalized','Position',[0 0 1 1])
plot3(wRData_helper(1,:),wRData_helper(2,:),wRData_helper(3,:),'k-o','MarkerSize',10,'LineWidth',3);
hold on
plot3(wRData_helper(1,1),wRData_helper(2,1),wRData_helper(3,1),'k^','MarkerSize',20,'LineWidth',3);
plot3(wRData_helper(1,end),wRData_helper(2,end),wRData_helper(3,end),'kv','MarkerSize',20,'LineWidth',3);
trisurf(points2,wRData_helper(1,:),wRData_helper(2,:),wRData_helper(3,:),'FaceColor','cyan');
xlabel('PC 2')
ylabel('PC 3')
zlabel('PC 5')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,[figPath filesep 'volumeExample2'     '_100ms_' sessionType],'svg');     saveas(gca,[figPath filesep 'volumeExample2'    '_100ms_' sessionType],'png')
%% Iteratively calculate distance between L and R as a fcn of PCs
wLData_helper = [];
wRData_helper = [];
wDist = [];
% Bootstrap 100 random trials for each condition due to unequal trials
numBootStraps = 500;
rng(1)
numPC = length(ConditionsExplained);
method = 'euclidean';

wSamp = [randi(size(wLProj{1},1), 1, numBootStraps); randi(size(wRProj{1},1), 1, numBootStraps)];
pSamp = [randi(size(pLProj{1},1), 1, numBootStraps); randi(size(pRProj{1},1), 1, numBootStraps)];

for k = 1:numBootStraps
    wLData_helper = [];
    wRData_helper = [];
    for i = 1:numPC         % Must start at 2 
        % Generate trial-by-trial data structure containing N PCs
        wLData_helper = [wLData_helper; squeeze(wLData(wSamp(1,k),:,i))];
        wRData_helper = [wRData_helper; squeeze(wRData(wSamp(2,k),:,i))];
        
        wDist(i,k,:) = diag(pdist2(wLData_helper', wRData_helper',method));
    end
end

wDistM = squeeze(mean(wDist,2));
wDistM_PC = [0; mean(wDistM,2)];
wDistM_PCSEM = [0; std(wDistM,[],2)/sqrt(numBootStraps)];

figure('Units','normalized','Position',[0 0 1 1])

shadedErrorBar(0:numPC,wDistM_PC,wDistM_PCSEM,'lineprops',{'b-','LineWidth',3})
hold on
%
pLData_helper = [];
pRData_helper = [];
pDist = [];

% For K random pairs of trials, calculate distance across time for every
% Ith principal componenet 
for k = 1:numBootStraps
    pLData_helper = [];
    pRData_helper = [];
    for i = 1:numPC         % Must start at 2 
        % Generate trial-by-trial data structure containing N PCs
        pLData_helper = [pLData_helper; squeeze(pLData(pSamp(1,k),:,i))];
        pRData_helper = [pRData_helper; squeeze(pRData(pSamp(2,k),:,i))];
    
        pDist(i,k,:) = diag(pdist2(pLData_helper', pRData_helper',method));
    end
end

pDistM = squeeze(mean(pDist,2));
pDistM_PC = [0; mean(pDistM,2)];
pDistM_PCSEM = [0; std(pDistM,[],2)/sqrt(numBootStraps)];

shadedErrorBar(0:numPC,pDistM_PC,pDistM_PCSEM,'lineprops',{'r-','LineWidth',3})
xlabel('PC Number')
ylabel('Distance')
title('Left vs. Right')
legend('Wistars','P rats')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,[figPath filesep 'LRDist_v_PC'     '_100ms_' sessionType],'svg');     saveas(gca,[figPath filesep 'LRDist_v_PC'    '_100ms_' sessionType],'png')

% Repeat analysis with different combinations of Corrections etc.
wLData_helper = [];
wLCData_helper = [];
wDist = [];
% Bootstrap 100 random trials for each condition due to unequal trials
numBootStraps = 500;
rng(1)
numPC = length(ConditionsExplained);
method = 'euclidean';

wSamp = [randi(size(wLProj{1},1), 1, numBootStraps); randi(size(wLCrct{1},1), 1, numBootStraps)];
pSamp = [randi(size(pLProj{1},1), 1, numBootStraps); randi(size(pLCrct{1},1), 1, numBootStraps)];

for k = 1:numBootStraps
    wLData_helper = [];
    wLCData_helper = [];
    for i = 1:numPC         % Must start at 2 
        % Generate trial-by-trial data structure containing N PCs
        wLData_helper = [wLData_helper; squeeze(wLData(wSamp(1,k),:,i))];
        wLCData_helper = [wLCData_helper; squeeze(wLCData(wSamp(2,k),:,i))];
        
        wDist(i,k,:) = diag(pdist2(wLData_helper', wLCData_helper',method));
    end
end

wDistM = squeeze(mean(wDist,2));
wDistM_PC = [0; mean(wDistM,2)];
wDistM_PCSEM = [0; std(wDistM,[],2)/sqrt(numBootStraps)];

figure('Units','normalized','Position',[0 0 1 1])

shadedErrorBar(0:numPC,wDistM_PC,wDistM_PCSEM,'lineprops',{'b-','LineWidth',3})
hold on
%
pLData_helper = [];
pLCData_helper = [];
pDist = [];

% For K random pairs of trials, calculate distance across time for every
% Ith principal componenet 
for k = 1:numBootStraps
    pLData_helper = [];
    pLCData_helper = [];
    for i = 1:numPC         % Must start at 2 
        % Generate trial-by-trial data structure containing N PCs
        pLData_helper = [pLData_helper; squeeze(pLData(pSamp(1,k),:,i))];
        pLCData_helper = [pLCData_helper; squeeze(pLCData(pSamp(2,k),:,i))];
    
        pDist(i,k,:) = diag(pdist2(pLData_helper', pLCData_helper',method));
    end
end

pDistM = squeeze(mean(pDist,2));
pDistM_PC = [0; mean(pDistM,2)];
pDistM_PCSEM = [0; std(pDistM,[],2)/sqrt(numBootStraps)];

shadedErrorBar(0:numPC,pDistM_PC,pDistM_PCSEM,'lineprops',{'r-','LineWidth',3})
xlabel('PC Number')
ylabel('Distance')
title('Left vs. Left Corrections')
legend('Wistars','P rats')

set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,[figPath filesep 'LLCDist_v_PC'     '_100ms_' sessionType],'svg');     saveas(gca,[figPath filesep 'LLCDist_v_PC'    '_100ms_' sessionType],'png')

% Right vs Right Corrections

wRData_helper = [];
wRCData_helper = [];
wDist = [];
% Bootstrap 100 random trials for each condition due to unequal trials
numBootStraps = 500;
rng(1)
numPC = length(ConditionsExplained);
method = 'euclidean';

wSamp = [randi(size(wRProj{1},1), 1, numBootStraps); randi(size(wRCrct{1},1), 1, numBootStraps)];
pSamp = [randi(size(pRProj{1},1), 1, numBootStraps); randi(size(pRCrct{1},1), 1, numBootStraps)];

for k = 1:numBootStraps
    wRData_helper = [];
    wRCData_helper = [];
    for i = 1:numPC         % Must start at 2 
        % Generate trial-by-trial data structure containing N PCs
        wRData_helper = [wRData_helper; squeeze(wRData(wSamp(1,k),:,i))];
        wRCData_helper = [wRCData_helper; squeeze(wRCData(wSamp(2,k),:,i))];
        
        wDist(i,k,:) = diag(pdist2(wRData_helper', wRCData_helper',method));
    end
end

wDistM = squeeze(mean(wDist,2));
wDistM_PC = [0; mean(wDistM,2)];
wDistM_PCSEM = [0; std(wDistM,[],2)/sqrt(numBootStraps)];

figure('Units','normalized','Position',[0 0 1 1])

shadedErrorBar(0:numPC,wDistM_PC,wDistM_PCSEM,'lineprops',{'b-','LineWidth',3})
hold on
%
pRData_helper = [];
pRCData_helper = [];
pDist = [];

% For K random pairs of trials, calculate distance across time for every
% Ith principal componenet 
for k = 1:numBootStraps
    pRData_helper = [];
    pRCData_helper = [];
    for i = 1:numPC         % Must start at 2 
        % Generate trial-by-trial data structure containing N PCs
        pRData_helper = [pRData_helper; squeeze(pRData(pSamp(1,k),:,i))];
        pRCData_helper = [pRCData_helper; squeeze(pRCData(pSamp(2,k),:,i))];
    
        pDist(i,k,:) = diag(pdist2(pRData_helper', pRCData_helper',method));
    end
end

pDistM = squeeze(mean(pDist,2));
pDistM_PC = [0; mean(pDistM,2)];
pDistM_PCSEM = [0; std(pDistM,[],2)/sqrt(numBootStraps)];

shadedErrorBar(0:numPC,pDistM_PC,pDistM_PCSEM,'lineprops',{'r-','LineWidth',3})
xlabel('PC Number')
ylabel('Distance')
title('Right vs. Right Corrections')
legend('Wistars','P rats')

set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,[figPath filesep 'RRCDist_v_PC'     '_100ms_' sessionType],'svg');     saveas(gca,[figPath filesep 'RRCDist_v_PC'    '_100ms_' sessionType],'png')

% Left Corrections vs Right Corrections

wLCData_helper = [];
wRCData_helper = [];
wDist = [];
% Bootstrap 100 random trials for each condition due to unequal trials
numBootStraps = 500;
rng(1)
numPC = length(ConditionsExplained);
method = 'euclidean';

wSamp = [randi(size(wLCrct{1},1), 1, numBootStraps); randi(size(wRCrct{1},1), 1, numBootStraps)];
pSamp = [randi(size(pLCrct{1},1), 1, numBootStraps); randi(size(pRCrct{1},1), 1, numBootStraps)];

for k = 1:numBootStraps
    wLCData_helper = [];
    wRCData_helper = [];
    for i = 1:numPC         % Must start at 2 
        % Generate trial-by-trial data structure containing N PCs
        wLCData_helper = [wLCData_helper; squeeze(wLCData(wSamp(1,k),:,i))];
        wRCData_helper = [wRCData_helper; squeeze(wRCData(wSamp(2,k),:,i))];
        
        wDist(i,k,:) = diag(pdist2(wLCData_helper', wRCData_helper',method));
    end
end

wDistM = squeeze(mean(wDist,2));
wDistM_PC = [0; mean(wDistM,2)];
wDistM_PCSEM = [0; std(wDistM,[],2)/sqrt(numBootStraps)];

figure('Units','normalized','Position',[0 0 1 1])

shadedErrorBar(0:numPC,wDistM_PC,wDistM_PCSEM,'lineprops',{'b-','LineWidth',3})
hold on
%
pLCData_helper = [];
pRCData_helper = [];
pDist = [];

% For K random pairs of trials, calculate distance across time for every
% Ith principal componenet 
for k = 1:numBootStraps
    pLCData_helper = [];
    pRCData_helper = [];
    for i = 1:numPC         % Must start at 2 
        % Generate trial-by-trial data structure containing N PCs
        pLCData_helper = [pLCData_helper; squeeze(pLCData(pSamp(1,k),:,i))];
        pRCData_helper = [pRCData_helper; squeeze(pRCData(pSamp(2,k),:,i))];
    
        pDist(i,k,:) = diag(pdist2(pLCData_helper', pRCData_helper',method));
    end
end

pDistM = squeeze(mean(pDist,2));
pDistM_PC = [0; mean(pDistM,2)];
pDistM_PCSEM = [0; std(pDistM,[],2)/sqrt(numBootStraps)];

shadedErrorBar(0:numPC,pDistM_PC,pDistM_PCSEM,'lineprops',{'r-','LineWidth',3})
xlabel('PC Number')
ylabel('Distance')
title('Left Corrections vs. Right Corrections')
legend('Wistars','P rats')

set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,[figPath filesep 'LCRCDist_v_PC'     '_100ms_' sessionType],'svg');     saveas(gca,[figPath filesep 'LCRCDist_v_PC'    '_100ms_' sessionType],'png')
