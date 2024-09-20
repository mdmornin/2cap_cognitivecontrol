%% Load and check data
clear all; close all % Refresh workspace
rng(1)
% PATHS
% For brain3 Linux Machine:
if ispc
    parentPath = 'F:/dissDat';
    figPath = 'F:/dissDat/figs/redirection';
else
    parentPath = '/research/dissDat';
    figPath = '/research/dissDat/figs';
end
% LOAD
if ~exist(figPath)
    mkdir(figPath)
end
load([parentPath filesep 'ephysStruct.mat']);
load([parentPath filesep 'trlStruct.mat']);
load([parentPath filesep 'masterTable.mat']);

addpath(genpath([parentPath filesep 'analysisScripts']))
%% Set indexing based on Congruent versus Incongruent Sessions
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
    pBin = 1/20; % 33 ms bins for fr variable, match video framerate
    PSTH_Approach = [];
    if ~isfield(ephysStruct,'PSTH_Approach')
        
        for i = 1:length(trlStruct)
    
            % Load trial times from trlStruct
            trialTimes = trlStruct(i).trialTimes(1:48);

            approachTimes = trlStruct(i).approach(1:48,1:2,2);
            approachTimes = approachTimes + trialTimes;

            [trialTimes, trialTimesIdx] = sort(trialTimes);
            approachTimesSorted = approachTimes(trialTimesIdx,:);
            approachIdx = trlStruct(i).approach(1:48,1,1) == 1 | trlStruct(i).approach(1:48,2,1) == 1;
            approachIdx = approachIdx(trialTimesIdx);
            Events(1,:) = trialTimes - 5;
            Events(2,:) = trialTimes;
            Events(3,:) = min(approachTimesSorted,[],2);
            Events(4,:) = max(approachTimesSorted,[],2);

            diffEvents = diff([Events(3,:); Events(4,:)]);
            idxTrl(diffEvents == 0) = 1;    % Correct Approaches
            idxTrl(diffEvents > 0) = 2;     % Redirections
            idxTrl(isnan(diffEvents)) = 3;  % No Approaches

             if startsWith(masterTbl.SessionType(i),'Regular')            
                correctIndex = trlStruct(i).approach(1:48,1,1) == 1 & trlStruct(i).approach(1:48,2,1) == 0;
                correctionIndex = trlStruct(i).approach(1:48,1,1) == 1 & trlStruct(i).approach(1:48,2,1) == 1;
                incorrectIndex = trlStruct(i).approach(1:48,2,1) == 1;
                
                correctIndexL = trlStruct(i).approach(1:24,1,1) == 1 & trlStruct(i).approach(1:24,2,1) == 0;
                correctIndexR = trlStruct(i).approach(25:48,1,1) == 1 & trlStruct(i).approach(25:48,2,1) == 0;
                correctionL = trlStruct(i).approach(:,1,2) > trlStruct(i).approach(:,2,2); 
                correctionR = trlStruct(i).approach(:,2,2) > trlStruct(i).approach(:,1,2);
             
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
            correctIndexL = (LRIdx(:,2) == 1);
            correctIndexR = (LRIdx(:,2) == 0);
            
            correctionLRIdx = double([correctionL; correctionR]);
            correctionLRIdx(1:24,2) = 1;
            correctionLRIdx = correctionLRIdx(trialTimesIdx,:);
            correctionIndexL = (correctionLRIdx(:,1) == 1 & correctionLRIdx(:,2) == 1);
            correctionIndexR = (correctionLRIdx(:,1) == 1 & correctionLRIdx(:,2) == 0);
            
            %
            epochBins = [11 41 41 11];
            
            for k = 1:length(idxTrl)
                if idxTrl(k) == 3
                    Events(3,k) = (rand(1,1) + randi([5,7],1)) + trialTimes(k);
                    Events(4,k) = Events(3,k);                   
                    rateVector(k,:,:) = makeTimeNormPSTH(ephysStruct(i).stmtx,Events(:,k),epochBins);
                else
                    rateVector(k,:,:) = makeTimeNormPSTH(ephysStruct(i).stmtx,Events(:,k),epochBins);
                end
            end
            PSTH = rateVector;
            rateVector = [];
            
            neurNum = size(PSTH,3);

            for k = 1:neurNum
                isimean = 1./(nanmean(PSTH(:,:,k)));
                isimean(isinf(isimean)) = NaN;
                isistd = nanstd(isimean);
                if isistd == 0
                    isistd = 0.000001;
                end
%                 sigma = nanmean(isimean)^(1/2) *  (1 / (isistd / nanmean(isimean)));
                sigma = 4;
                PSTH(:,:,k) = smoothdata(PSTH(:,:,k),2,'gaussian',sigma);
                
            end

%             sigma = 5;
%             PSTH = smoothdata(PSTH,2,'gaussian',sigma);

            PSTH = zscore(PSTH,0,[1 3]);


            ephysStruct(i).PSTH = PSTH;
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

        ephysStruct(i).PSTH_LeftMean = squeeze(mean(ephysStruct(i).PSTH_Left,1));
        ephysStruct(i).PSTH_RightMean = squeeze(mean(ephysStruct(i).PSTH_Right,1));

        ephysStruct(i).correctionsIndex = correctionIndex;


        PSTH_Approach = [];
        PSTH_NoApproach = [];
        PSTH = [];
    end

    end
   

    %% Based on data above, there are no real differences within Genotype (WCon,WInc OR PCon,PInc) coefficient loadings. However, there are consistent differences between genotype loadings
% Therefore, let us analyze P and Wistars separate from one another. 
% This will allow within-session type comparisons but not necessarily allow
% between genotype comparisons except in certain circumstances such as
% comparing the distance metrics between two PC subspaces. 

% Depending on the session types we would like to analyze, we can either
% take the index of all P and W rats or we can take the index of ONLY P / W
% rats in Congruent Sessions

PIdx = startsWith(masterTbl.Strain,'P') & startsWith(masterTbl.SessionType,'Regular');
WIdx = startsWith(masterTbl.Strain,'W') & startsWith(masterTbl.SessionType,'Regular');

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
% Incongruent
pMatLIC = [ephysStruct(revPIdx).PSTH_LeftMean];
wMatLIC = [ephysStruct(revWIdx).PSTH_LeftMean];

pMatRIC = [ephysStruct(revPIdx).PSTH_RightMean];
wMatRIC = [ephysStruct(revWIdx).PSTH_RightMean];

conditionIndexes = [length(wMatLC) length(wMatLIC) length(pMatLC) length(pMatLIC)];
conditionLabels = [{'Congruent W'}, {'Incongruent W'}, {'Congruent P'}, {'Incongruent P'}];


%% Reshape Matrices 
% For the purposes of the coefficient matrix:
% allMat contains all Left choices followed in time by a matrix of all of
% the Right choices
% Each individual matrix itself repeats the following pattern:
% Wistar Congruent, Wistar Incongruent, P Congruent, P Incongruent
% This will be important to pull out the coefficients pertaining to each
% trial type 
% From there, we can take the coefficients of known datasets and determine
% if there are differences between the L/R
% 
% 08/22/2022 - - First, determine if there are baseline differences between
% L and R choices in the Congruent sessions prior to the addition of
% incongruent choices.
% In order to do this, we just take the following modification of the above
% pattern:
% Congruent Wistar Left Choices, Congruent P Rat Left Choices, Congruent
% Wistar Right Choices, Congruent P Rats Right Choices.

allMat = [];

% Below contain ALL session types
% allMatL = [wMatLC wMatLIC pMatLC pMatLIC];
% allMatR = [wMatRC wMatRIC pMatRC pMatRIC];

% Below contain CONGRUENT only sessino types
allMatL = [wMatLC pMatLC];
allMatR = [wMatRC pMatRC];
conditionIndexes = [length(wMatLC) length(pMatLC)];
conditionLabels = [{'Congruent W'}, {'Congruent P'}];
% Format matrix such that it expands 'in time / observations'
allMat = [allMatL; allMatR];

%% PCA Preprocessing | Separation between P, W groups | LEFT TRIALS



[ConditionsCoeff,ConditionsScore,ConditionsLatent,~,ConditionsExplained] = pca(allMat);

% Need to find an index of neurons per dataset, be able to sort that index
% into our conditions (W/P; Con/Inc). Should ideally be a 1x54 matrix
% containing N # of neurons.

PneuronIndex = [ephysStruct(PIdx).neuronNumApproach];
WneuronIndex = [ephysStruct(WIdx).neuronNumApproach];
%% Plot PCs over Time without Separations
timeVec = [-5:0.1:2.1];


figure('Units','normalized','Position',[0 0 1 1])
plot(ConditionsScore(:,1),'b-','LineWidth',3)
hold on
plot(ConditionsScore(:,2),'r-','LineWidth',3)
plot(ConditionsScore(:,3),'k-','LineWidth',3)

xlabel('Time Bins')
ylabel('PC Score')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
% xline(11,'--','LApproach')
% xline(51,'--','RApproach')
saveas(gca,[figPath filesep 'PCs_Time'     '_100ms'],'svg');     saveas(gca,[figPath filesep 'PCs_Time'    '_100ms'],'png')

%% Variance explained
figure('Units','normalized','Position',[0 0 1 1])
plot(cumsum(ConditionsExplained),'ko');
xlabel('PC Number')
ylabel('Cumulative Explained Variance')
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

wConScore = rawData(wConIdx,:)'*ConditionsCoeff(wConIdx,:);
% wIncScore = rawData(wIncIdx,:)'*ConditionsCoeff(wIncIdx,:);

pConScore = rawData(pConIdx,:)'*ConditionsCoeff(pConIdx,:);
% pIncScore = rawData(pIncIdx,:)'*ConditionsCoeff(pIncIdx,:);

%%
if allCond == 1
    allPCProjections = [{wConScore} {wIncScore} {pConScore} {pIncScore}];
    allPCCoeff = [{ConditionsCoeff(wConIdx,:)} {ConditionsCoeff(wIncIdx,:)} {ConditionsCoeff(pConIdx,:)} {ConditionsCoeff(pIncIdx,:)}];
    titles = conditionLabels;
else
    allPCProjections = [{wConScore} {pConScore} ];
    allPCCoeff = [{ConditionsCoeff(wConIdx,:)} {ConditionsCoeff(pConIdx,:)}];
    titles = conditionLabels;
end
% %% Plot Separated PCs Over Time
pcSort = 1:5;
colors = [{'#0072BD'} {'#4DBEEE'} {'#77AC30'} {'#D95319'} {'#A2142F'} {'#0072BD'} {'#4DBEEE'} {'#77AC30'} {'#D95319'} {'#A2142F'} {'#0072BD'} {'#4DBEEE'} {'#77AC30'} {'#D95319'} {'#A2142F'} {'#0072BD'} {'#4DBEEE'} {'#77AC30'} {'#D95319'} {'#A2142F'} {'#0072BD'} {'#4DBEEE'} {'#77AC30'} {'#D95319'} {'#A2142F'} {'#0072BD'} {'#4DBEEE'} {'#77AC30'} {'#D95319'} {'#A2142F'} {'#0072BD'} {'#4DBEEE'} {'#77AC30'} {'#D95319'} {'#A2142F'} {'#0072BD'} {'#4DBEEE'} {'#77AC30'} {'#D95319'} {'#A2142F'} {'#0072BD'} {'#4DBEEE'} {'#77AC30'} {'#D95319'} {'#A2142F'} {'#0072BD'} {'#4DBEEE'} {'#77AC30'} {'#D95319'} {'#A2142F'}];
time = 1:length(wConScore)/2;
timeIdx = length(wConScore)/2;
cueOn = 0;

figure('Units','normalized','Position',[0 0 1 1])

for j = 1:length(allPCProjections)
    for i=1:length(pcSort)
        
        subplot(1,2,j);
        plot(time,allPCProjections{j}(1:timeIdx,pcSort(i)),'LineWidth',3,'Color',colors{i},'LineStyle','-')

        hold on
        plot(time,allPCProjections{j}(timeIdx+1:end,pcSort(i)),'LineWidth',3,'Color',colors{i},'LineStyle','--')
        xlabel('Time (s)');
        title([titles{j}]);
        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',10,'FontWeight','bold','LineWidth',4)
%         xline(cueOn,'k--','Approach','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
    end
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    if j == 1; legend([{'PC1'},{'PC2'},{'PC3'},{'PC4'},{'PC5'}],'AutoUpdate','off','location','best'); end

end

saveas(gca,[figPath filesep 'allNeurons_split_top5pcs_CongruentOnly'  num2str(pcSort)  '_100ms'],'svg');     saveas(gca,[figPath filesep 'allNeurons_split_top5pcs_CongruentOnly'  num2str(pcSort)  '_100ms'],'png')
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
    numPC = 1:10;
    PSTH = ephysStruct(dataPos(i)).PSTH;
    % Mean Center PSTH
    PSTH = (mean(PSTH,3) - PSTH);
    for k = 1:size(PSTH,1)
        for j = 1:length(numPC)
            trlProj{i,j}(k,:,:) = squeeze(PSTH(k,:,:)) * indexedCoeff{i}(:,j);
        end
    end    
end
% This leaves us with a Number of Session x Number of PC Cell Array
% Cell Aray has 48 Trials x 71 Observations
% Now, we can index on the cell aray to obtain L/R Trial Types

for i = 1:length(dataPos)
    correctL = ephysStruct(dataPos(i)).LeftIndex;
    correctR = ephysStruct(dataPos(i)).RightIndex;
    correctionsIndex = ephysStruct(dataPos(i)).correctionsIndex;
    correctionsL = ephysStruct(dataPos(i)).correctionLeftIndex;
    correctionsR = ephysStruct(dataPos(i)).correctionRightIndex;

    for j = 1:length(numPC)
        leftProj{i,j} = trlProj{i,j}(correctL,:,:);
        rightProj{i,j} = trlProj{i,j}(correctR,:,:);
        corrections{i,j} = trlProj{i,j}(correctionsIndex,:,:);
        correctionsLProj{i,j} = trlProj{i,j}(correctionsL,:,:);
        correctionsRProj{i,j} = trlProj{i,j}(correctionsR,:,:);

    end
end

% 2 Cell Arrays now exist with Left and Right Projections
% Can now take the mean for each in each genotype and each PC#

for j = 1:length(numPC)
    wLProj{j} = vertcat(leftProj{1:13,j}); 
    wRProj{j} = vertcat(rightProj{1:13,j}); 
    wCrctn{j} = vertcat(corrections{1:13,j});
    wLCrct{j} = vertcat(correctionsLProj{1:13,j});
    wRCrct{j} = vertcat(correctionsRProj{1:13,j});
    

    pLProj{j} = vertcat(leftProj{14:end,j});
    pRProj{j} = vertcat(rightProj{14:end,j});
    pCrctn{j} = vertcat(corrections{14:end,j});
    pLCrct{j} = vertcat(correctionsLProj{14:end,j});
    pRCrct{j} = vertcat(correctionsRProj{14:end,j});

end

%% Quick plots for simple data vis
% Wistars
allPerms = nchoosek([1:5],2);

for i = 1:length(allPerms)
    figure('Units','normalized','Position',[0 0 1 1])

    plot(mean(wLProj{i}),'LineWidth',3); hold on; plot(mean(wRProj{i}),'LineWidth',3,'LineStyle','--'); plot(mean(wLCrct{i}),'LineWidth',3); plot(mean(wRCrct{i}),'LineWidth',3,'LineStyle','--');
    title(['Ws | PC Number ' num2str(i)])
    legend([{'Left'},{'Right'},{'Correction to Left'},{'Correction to Right'}],'AutoUpdate','off','location','best');
    xlabel('Time (Bins)')
    ylabel('PC Score')
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figPath filesep  'WPCTime_CongruentOnly_PC_' num2str(i)  '_100ms'],'svg');     saveas(gca,[figPath filesep 'WPCTime_CongruentOnly_PC_' num2str(i) '_100ms'],'png')

    % P Rats
    figure('Units','normalized','Position',[0 0 1 1])

    plot(mean(pLProj{i}),'LineWidth',3); hold on; plot(mean(pRProj{i}),'LineWidth',3,'LineStyle','--'); plot(mean(pLCrct{i}),'LineWidth',3); plot(mean(pRCrct{i}),'LineWidth',3,'LineStyle','--');
    legend([{'Left'},{'Right'},{'Correction to Left'},{'Correction to Right'}],'AutoUpdate','off','location','best');
    title(['Ps | PC Number ' num2str(i)])
    xlabel('Time (Bins)')
    ylabel('PC Score')
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
%     saveas(gca,[figPath filesep  'PPCTime_CongruentOnly_PC_' num2str(i)  '_100ms'],'svg');     saveas(gca,[figPath filesep  'PPCTime_CongruentOnly_PC_' num2str(i)  '_100ms'],'png')
end
%% 3D Plots
close all
allPerms = nchoosek([1:5],3);
% Wistars
for i = 1:length(allPerms)
    figure('Units','normalized','Position',[0 0 1 1])
    plot3(mean(wLProj{allPerms(i,1)}),mean(wLProj{allPerms(i,2)}),mean(wLProj{allPerms(i,3)}),'LineWidth',3,'Color',colors{i},'LineStyle','-')
    hold on
    plot3(mean(wLProj{allPerms(i,1)}(:,1)),mean(wLProj{allPerms(i,2)}(:,1)),mean(wLProj{allPerms(i,3)}(:,1)),'ko','MarkerSize',15,'MarkerFaceColor','k')

    plot3(mean(wRProj{allPerms(i,1)}),mean(wRProj{allPerms(i,2)}),mean(wRProj{allPerms(i,3)}),'LineWidth',3,'Color',colors{i},'LineStyle','--')
    plot3(mean(wRProj{allPerms(i,1)}(:,1)),mean(wRProj{allPerms(i,2)}(:,1)),mean(wRProj{allPerms(i,3)}(:,1)),'ko','MarkerSize',15,'MarkerFaceColor','k')

    xlabel(['PC' num2str(allPerms(i,1))]);
    ylabel(['PC' num2str(allPerms(i,2))]);
    zlabel(['PC' num2str(allPerms(i,3))]);
    title('Wistar Congruent')

    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figPath filesep '3DPROJ_TrialMeans_Wistars'  num2str(allPerms(i,:))  '_100ms'],'svg');     saveas(gca,[figPath filesep '3DPROJ_TrialMeans_Wistars' num2str(allPerms(i,:))  '_100ms'],'fig')
end
% Left and Redirect Left
for i = 1:length(allPerms)
    figure('Units','normalized','Position',[0 0 1 1])
    plot3(mean(wLProj{allPerms(i,1)}),mean(wLProj{allPerms(i,2)}),mean(wLProj{allPerms(i,3)}),'LineWidth',3,'Color',colors{i},'LineStyle','-')
    hold on
    plot3(mean(wLProj{allPerms(i,1)}(:,1)),mean(wLProj{allPerms(i,2)}(:,1)),mean(wLProj{allPerms(i,3)}(:,1)),'ko','MarkerSize',15,'MarkerFaceColor','k')

    plot3(mean(wLCrct{allPerms(i,1)}),mean(wLCrct{allPerms(i,2)}),mean(wLCrct{allPerms(i,3)}),'LineWidth',3,'Color',colors{i},'LineStyle','--')
    plot3(mean(wLCrct{allPerms(i,1)}(:,1)),mean(wLCrct{allPerms(i,2)}(:,1)),mean(wLCrct{allPerms(i,3)}(:,1)),'ko','MarkerSize',15,'MarkerFaceColor','k')

    xlabel(['PC' num2str(allPerms(i,1))]);
    ylabel(['PC' num2str(allPerms(i,2))]);
    zlabel(['PC' num2str(allPerms(i,3))]);
    title('Wistar Congruent | Left + Redirect Left')

    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figPath filesep '3DPROJ_TrialMeans_Wistars_Left_RedirectLeft'  num2str(allPerms(i,:))  '_100ms'],'svg');     saveas(gca,[figPath filesep '3DPROJ_TrialMeans_Wistars_Left_RedirectLeft' num2str(allPerms(i,:))  '_100ms'],'fig')
end

for i = 1:length(allPerms)
    figure('Units','normalized','Position',[0 0 1 1])
    plot3(mean(wRProj{allPerms(i,1)}),mean(wRProj{allPerms(i,2)}),mean(wRProj{allPerms(i,3)}),'LineWidth',3,'Color',colors{i},'LineStyle','-')
    hold on
    plot3(mean(wRProj{allPerms(i,1)}(:,1)),mean(wRProj{allPerms(i,2)}(:,1)),mean(wRProj{allPerms(i,3)}(:,1)),'ko','MarkerSize',15,'MarkerFaceColor','k')

    plot3(mean(wRCrct{allPerms(i,1)}),mean(wRCrct{allPerms(i,2)}),mean(wRCrct{allPerms(i,3)}),'LineWidth',3,'Color',colors{i},'LineStyle','--')
    plot3(mean(wRCrct{allPerms(i,1)}(:,1)),mean(wRCrct{allPerms(i,2)}(:,1)),mean(wRCrct{allPerms(i,3)}(:,1)),'ko','MarkerSize',15,'MarkerFaceColor','k')

    xlabel(['PC' num2str(allPerms(i,1))]);
    ylabel(['PC' num2str(allPerms(i,2))]);
    zlabel(['PC' num2str(allPerms(i,3))]);
    title('Wistar Congruent | Right + Redirect Right')

    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figPath filesep '3DPROJ_TrialMeans_Wistars_Right_RedirectRight'  num2str(allPerms(i,:))  '_100ms'],'svg');     saveas(gca,[figPath filesep '3DPROJ_TrialMeans_Wistars_Right_RedirectRight' num2str(allPerms(i,:))  '_100ms'],'fig')
end

% P rats
for i = 1:length(allPerms)
    figure('Units','normalized','Position',[0 0 1 1])
    plot3(mean(pLProj{allPerms(i,1)}),mean(pLProj{allPerms(i,2)}),mean(pLProj{allPerms(i,3)}),'LineWidth',3,'Color',colors{i},'LineStyle','-')
    hold on
    plot3(mean(pLProj{allPerms(i,1)}(:,1)),mean(pLProj{allPerms(i,2)}(:,1)),mean(pLProj{allPerms(i,3)}(:,1)),'ko','MarkerSize',15,'MarkerFaceColor','k')

    plot3(mean(pRProj{allPerms(i,1)}),mean(pRProj{allPerms(i,2)}),mean(pRProj{allPerms(i,3)}),'LineWidth',3,'Color',colors{i},'LineStyle','--')
    plot3(mean(pRProj{allPerms(i,1)}(:,1)),mean(pRProj{allPerms(i,2)}(:,1)),mean(pRProj{allPerms(i,3)}(:,1)),'ko','MarkerSize',15,'MarkerFaceColor','k')

    xlabel(['PC' num2str(allPerms(i,1))]);
    ylabel(['PC' num2str(allPerms(i,2))]);
    zlabel(['PC' num2str(allPerms(i,3))]);
    title('P Congruent')

    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figPath filesep '3DPROJ_TrialMeans_PRats'  num2str(allPerms(i,:))  '_100ms'],'svg');     saveas(gca,[figPath filesep '3DPROJ_TrialMeans_PRats' num2str(allPerms(i,:))  '_100ms'],'fig')
end

for i = 1:length(allPerms)
    figure('Units','normalized','Position',[0 0 1 1])
    plot3(mean(pLProj{allPerms(i,1)}),mean(pLProj{allPerms(i,2)}),mean(pLProj{allPerms(i,3)}),'LineWidth',3,'Color',colors{i},'LineStyle','-')
    hold on
    plot3(mean(pLProj{allPerms(i,1)}(:,1)),mean(pLProj{allPerms(i,2)}(:,1)),mean(pLProj{allPerms(i,3)}(:,1)),'ko','MarkerSize',15,'MarkerFaceColor','k')

    plot3(mean(pLCrct{allPerms(i,1)}),mean(pLCrct{allPerms(i,2)}),mean(pLCrct{allPerms(i,3)}),'LineWidth',3,'Color',colors{i},'LineStyle','--')
    plot3(mean(pLCrct{allPerms(i,1)}(:,1)),mean(pLCrct{allPerms(i,2)}(:,1)),mean(pLCrct{allPerms(i,3)}(:,1)),'ko','MarkerSize',15,'MarkerFaceColor','k')

    xlabel(['PC' num2str(allPerms(i,1))]);
    ylabel(['PC' num2str(allPerms(i,2))]);
    zlabel(['PC' num2str(allPerms(i,3))]);
    title('P Rats Congruent | Left + Redirect Left')

    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figPath filesep '3DPROJ_TrialMeans_P_Left_RedirectLeft'  num2str(allPerms(i,:))  '_100ms'],'svg');     saveas(gca,[figPath filesep '3DPROJ_TrialMeans_P_Left_RedirectLeft' num2str(allPerms(i,:))  '_100ms'],'fig')
end

for i = 1:length(allPerms)
    figure('Units','normalized','Position',[0 0 1 1])
    plot3(mean(pRProj{allPerms(i,1)}),mean(pRProj{allPerms(i,2)}),mean(pRProj{allPerms(i,3)}),'LineWidth',3,'Color',colors{i},'LineStyle','-')
    hold on
    plot3(mean(pRProj{allPerms(i,1)}(:,1)),mean(pRProj{allPerms(i,2)}(:,1)),mean(pRProj{allPerms(i,3)}(:,1)),'ko','MarkerSize',15,'MarkerFaceColor','k')

    plot3(mean(pRCrct{allPerms(i,1)}),mean(pRCrct{allPerms(i,2)}),mean(pRCrct{allPerms(i,3)}),'LineWidth',3,'Color',colors{i},'LineStyle','--')
    plot3(mean(pRCrct{allPerms(i,1)}(:,1)),mean(pRCrct{allPerms(i,2)}(:,1)),mean(pRCrct{allPerms(i,3)}(:,1)),'ko','MarkerSize',15,'MarkerFaceColor','k')

    xlabel(['PC' num2str(allPerms(i,1))]);
    ylabel(['PC' num2str(allPerms(i,2))]);
    zlabel(['PC' num2str(allPerms(i,3))]);
    title('P Rats Congruent | Right + Redirect Right')

    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figPath filesep '3DPROJ_TrialMeans_P_Right_RedirectRight'  num2str(allPerms(i,:))  '_100ms'],'svg');     saveas(gca,[figPath filesep '3DPROJ_TrialMeans_P_Right_RedirectRight' num2str(allPerms(i,:))  '_100ms'],'fig')
end


%% Calculate Pairwise Distances Between Trial Types 
% Find the distance between L and R trial types for a given PC subspace.
% Determine if this distance is greater or lesser in P rats versus Wistars
% in Congruent sessions. 

% Therea are inequal trial amounts for L and R. One remedy to this is to
% subsample across the population of trials. The minimum amount of trials
% for a given type is 15 within the Wistar Right Corrections. This may be
% used for comparing all 4 trial types (Left, Right, Correct Left, Correct
% Right). 

% For just the Left and Right, we can use 50 trials per bootstrapped sample
% This gives a decent amount of sample size without sacrificing data.
close all
numBootStraps = 1000;
distanceMeasurement = 'euclidean';
wDistance = []; pDistance = [];
for i = 1:numBootStraps

    wSamp(i) = datasample(1:min(size(wLProj{1},1), size(wRProj{1},1)),1);
    pSamp(i) = datasample(1:min(size(pLProj{1},1), size(pRProj{1},1)),1);
     
    wDistance(i,:) = diag(pdist2([wLProj{1}(wSamp(i),:); wLProj{2}(wSamp(i),:); wLProj{3}(wSamp(i),:); wLProj{4}(wSamp(i),:); wLProj{5}(wSamp(i),:)]',[wRProj{1}(wSamp(i),:); wRProj{2}(wSamp(i),:); wRProj{3}(wSamp(i),:); wRProj{4}(wSamp(i),:); wRProj{5}(wSamp(i),:)]',distanceMeasurement));
    pDistance(i,:) = diag(pdist2([pLProj{1}(pSamp(i),:); pLProj{2}(pSamp(i),:); pLProj{3}(pSamp(i),:); pLProj{4}(pSamp(i),:); pLProj{5}(pSamp(i),:)]',[pRProj{1}(pSamp(i),:); pRProj{2}(pSamp(i),:); pRProj{3}(pSamp(i),:); pRProj{4}(pSamp(i),:); pRProj{5}(pSamp(i),:)]',distanceMeasurement));

end
figure('Units','normalized','Position',[0 0 1 1])
shadedErrorBar(time,mean(wDistance),std(wDistance)/sqrt(numBootStraps),'lineprops',{'b-','LineWidth',3});
hold on
shadedErrorBar(time,mean(pDistance),std(pDistance)/sqrt(numBootStraps),'lineprops',{'r-','LineWidth',3});
xlabel('Time (s)')
ylabel(['Distance ' distanceMeasurement])

set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,[figPath filesep 'LR_Distance' '_100ms'],'svg');     saveas(gca,[figPath filesep 'LR_Distance'  '_100ms'],'fig')

figure
boxplot([mean(wDistance,2) mean(pDistance,2)],'notch','on')

%
wDistance = []; pDistance = [];
for i = 1:numBootStraps

    wSamp(i) = datasample(1:min(size(wLProj{1},1), size(wLCrct{1},1)),1);
    pSamp(i) = datasample(1:min(size(pLProj{1},1), size(pLCrct{1},1)),1);
     
    wDistance(i,:) = diag(pdist2([wLProj{1}(wSamp(i),:); wLProj{2}(wSamp(i),:); wLProj{3}(wSamp(i),:); wLProj{4}(wSamp(i),:); wLProj{5}(wSamp(i),:)]',[wLCrct{1}(wSamp(i),:); wLCrct{2}(wSamp(i),:); wLCrct{3}(wSamp(i),:); wLCrct{4}(wSamp(i),:); wLCrct{5}(wSamp(i),:)]',distanceMeasurement));
    pDistance(i,:) = diag(pdist2([pLProj{1}(pSamp(i),:); pLProj{2}(pSamp(i),:); pLProj{3}(pSamp(i),:); pLProj{4}(pSamp(i),:); pLProj{5}(pSamp(i),:)]',[pLCrct{1}(pSamp(i),:); pLCrct{2}(pSamp(i),:); pLCrct{3}(pSamp(i),:); pLCrct{4}(pSamp(i),:); pLCrct{5}(pSamp(i),:)]',distanceMeasurement));

end
figure('Units','normalized','Position',[0 0 1 1])
shadedErrorBar(time,mean(wDistance),std(wDistance)/sqrt(numBootStraps),'lineprops',{'b-','LineWidth',3});
hold on
shadedErrorBar(time,mean(pDistance),std(pDistance)/sqrt(numBootStraps),'lineprops',{'r-','LineWidth',3});
xlabel('Time (s)')
ylabel(['Distance ' distanceMeasurement])

set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,[figPath filesep 'L_LCorrection_Distance' '_100ms'],'svg');     saveas(gca,[figPath filesep 'L_LCorrection_Distance'  '_100ms'],'fig')

figure
boxplot([mean(wDistance,2) mean(pDistance,2)],'notch','on')

%
wDistance = []; pDistance = [];
for i = 1:numBootStraps

    wSamp(i) = datasample(1:min(size(wRProj{1},1), size(wRCrct{1},1)),1);
    pSamp(i) = datasample(1:min(size(pRProj{1},1), size(pRCrct{1},1)),1);
     
    wDistance(i,:) = diag(pdist2([wRCrct{1}(wSamp(i),:); wRCrct{2}(wSamp(i),:); wRCrct{3}(wSamp(i),:); wRCrct{4}(wSamp(i),:); wRCrct{5}(wSamp(i),:)]',[wRProj{1}(wSamp(i),:); wRProj{2}(wSamp(i),:); wRProj{3}(wSamp(i),:); wRProj{4}(wSamp(i),:); wRProj{5}(wSamp(i),:)]',distanceMeasurement));
    pDistance(i,:) = diag(pdist2([pRCrct{1}(pSamp(i),:); pRCrct{2}(pSamp(i),:); pRCrct{3}(pSamp(i),:); pRCrct{4}(pSamp(i),:); pRCrct{5}(pSamp(i),:)]',[pRProj{1}(pSamp(i),:); pRProj{2}(pSamp(i),:); pRProj{3}(pSamp(i),:); pRProj{4}(pSamp(i),:); pRProj{5}(pSamp(i),:)]',distanceMeasurement));


end
figure('Units','normalized','Position',[0 0 1 1])
shadedErrorBar(time,mean(wDistance),std(wDistance)/sqrt(numBootStraps),'lineprops',{'b-','LineWidth',3});
hold on
shadedErrorBar(time,mean(pDistance),std(pDistance)/sqrt(numBootStraps),'lineprops',{'r-','LineWidth',3});
xlabel('Time (s)')
ylabel(['Distance ' distanceMeasurement])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,[figPath filesep 'R_RCorrection_Distance' '_100ms'],'svg');     saveas(gca,[figPath filesep 'R_RCorrection_Distance'  '_100ms'],'fig')

figure
boxplot([mean(wDistance,2) mean(pDistance,2)],'notch','on')

%% Create digestible graphs that tell the distance between each epoch of interest
% Each epoch is 25 bins long
% Create pariwise distance graphs for each epoch 
% 
close all
numBootStraps = 500;
distanceMeasurement = 'euclidean';
wDistanceLR = []; pDistanceLR = [];
wVelLR = []; pVelLR = [];
wDistanceLLCR = []; pDistanceLLCR = [];
wVelLLCR = []; pVelLLCR = [];
wDistanceRRCR = []; pDistanceRRCR = [];
wVelRRCR = []; pVelRRCR = [];

for i = 1:numBootStraps

    wSamp(i) = datasample(1:min(size(wLProj{1},1), size(wRProj{1},1)),1);
    pSamp(i) = datasample(1:min(size(pLProj{1},1), size(pRProj{1},1)),1);
     
    wDistanceLR(i,:) = diag(pdist2([wLProj{1}(wSamp(i),:); wLProj{2}(wSamp(i),:); wLProj{3}(wSamp(i),:); wLProj{4}(wSamp(i),:); wLProj{5}(wSamp(i),:)]',[wRProj{1}(wSamp(i),:); wRProj{2}(wSamp(i),:); wRProj{3}(wSamp(i),:); wRProj{4}(wSamp(i),:); wRProj{5}(wSamp(i),:)]',distanceMeasurement));
    pDistanceLR(i,:) = diag(pdist2([pLProj{1}(pSamp(i),:); pLProj{2}(pSamp(i),:); pLProj{3}(pSamp(i),:); pLProj{4}(pSamp(i),:); pLProj{5}(pSamp(i),:)]',[pRProj{1}(pSamp(i),:); pRProj{2}(pSamp(i),:); pRProj{3}(pSamp(i),:); pRProj{4}(pSamp(i),:); pRProj{5}(pSamp(i),:)]',distanceMeasurement));
    
    wVelLR(i,:) = ([0 diff(wDistanceLR(i,:),[],2)]);
    pVelLR(i,:) = ([0 diff(pDistanceLR(i,:),[],2)]);

    wSamp(i) = datasample(1:min(size(wLProj{1},1), size(wLCrct{1},1)),1);
    pSamp(i) = datasample(1:min(size(pLProj{1},1), size(pLCrct{1},1)),1);
     
    wDistanceLLCR(i,:) = diag(pdist2([wLProj{1}(wSamp(i),:); wLProj{2}(wSamp(i),:); wLProj{3}(wSamp(i),:); wLProj{4}(wSamp(i),:); wLProj{5}(wSamp(i),:)]',[wLCrct{1}(wSamp(i),:); wLCrct{2}(wSamp(i),:); wLCrct{3}(wSamp(i),:); wLCrct{4}(wSamp(i),:); wLCrct{5}(wSamp(i),:)]',distanceMeasurement));
    pDistanceLLCR(i,:) = diag(pdist2([pLProj{1}(pSamp(i),:); pLProj{2}(pSamp(i),:); pLProj{3}(pSamp(i),:); pLProj{4}(pSamp(i),:); pLProj{5}(pSamp(i),:)]',[pLCrct{1}(pSamp(i),:); pLCrct{2}(pSamp(i),:); pLCrct{3}(pSamp(i),:); pLCrct{4}(pSamp(i),:); pLCrct{5}(pSamp(i),:)]',distanceMeasurement));

    wVelLLCR(i,:) = ([0 diff(wDistanceLLCR(i,:),[],2)]);
    pVelLLCR(i,:) = ([0 diff(pDistanceLLCR(i,:),[],2)]);

    wSamp(i) = datasample(1:min(size(wRProj{1},1), size(wRCrct{1},1)),1);
    pSamp(i) = datasample(1:min(size(pRProj{1},1), size(pRCrct{1},1)),1);
     
    wDistanceRRCR(i,:) = diag(pdist2([wRCrct{1}(wSamp(i),:); wRCrct{2}(wSamp(i),:); wRCrct{3}(wSamp(i),:); wRCrct{4}(wSamp(i),:); wRCrct{5}(wSamp(i),:)]',[wRProj{1}(wSamp(i),:); wRProj{2}(wSamp(i),:); wRProj{3}(wSamp(i),:); wRProj{4}(wSamp(i),:); wRProj{5}(wSamp(i),:)]',distanceMeasurement));
    pDistanceRRCR(i,:) = diag(pdist2([pRCrct{1}(pSamp(i),:); pRCrct{2}(pSamp(i),:); pRCrct{3}(pSamp(i),:); pRCrct{4}(pSamp(i),:); pRCrct{5}(pSamp(i),:)]',[pRProj{1}(pSamp(i),:); pRProj{2}(pSamp(i),:); pRProj{3}(pSamp(i),:); pRProj{4}(pSamp(i),:); pRProj{5}(pSamp(i),:)]',distanceMeasurement));

    wVelRRCR(i,:) = ([0 diff(wDistanceRRCR(i,:),[],2)]);
    pVelRRCR(i,:) = ([0 diff(pDistanceRRCR(i,:),[],2)]);
    
end

%% Create boxplots for distance
close all

figure('Units','normalized','Position',[0 0 1 1])
subplot(2,3,1)
boxplot([nanmean(wDistanceLR(:,1:25),2) nanmean(wDistanceLR(:,26:50),2) nanmean(wDistanceLR(:,51:75),2) nanmean(wDistanceLR(:,76:100),2)],'notch','on','symbol','')
ylim([-0.1 1.6])
title('Wistar | Left, Right')
xticklabels([{'Epoch1'},{'Epoch2'},{'Epoch3'}, {'Epoch4'}])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

subplot(2,3,4)
boxplot([nanmean(pDistanceLR(:,1:25),2) nanmean(pDistanceLR(:,26:50),2) nanmean(pDistanceLR(:,51:75),2) nanmean(pDistanceLR(:,76:100),2)],'notch','on','symbol','')
ylim([-0.1 1.6])
title('P rat | Left, Right')
xticklabels([{'Epoch1'},{'Epoch2'},{'Epoch3'}, {'Epoch4'}])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

subplot(2,3,2)
boxplot([nanmean(wDistanceLLCR(:,1:25),2) nanmean(wDistanceLLCR(:,26:50),2) nanmean(wDistanceLLCR(:,51:75),2) nanmean(wDistanceLLCR(:,76:100),2)],'notch','on','symbol','')
ylim([-0.1 1.6])
title('Wistar | Left, Redirect Left')
xticklabels([{'Epoch1'},{'Epoch2'},{'Epoch3'}, {'Epoch4'}])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

subplot(2,3,5)
boxplot([nanmean(pDistanceLLCR(:,1:25),2) nanmean(pDistanceLLCR(:,26:50),2) nanmean(pDistanceLLCR(:,51:75),2) nanmean(pDistanceLLCR(:,76:100),2)],'notch','on','symbol','')
ylim([-0.1 1.6])
title('P rat | Left, Redirect Left')
xticklabels([{'Epoch1'},{'Epoch2'},{'Epoch3'}, {'Epoch4'}])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

subplot(2,3,3)
boxplot([nanmean(wDistanceRRCR(:,1:25),2) nanmean(wDistanceRRCR(:,26:50),2) nanmean(wDistanceRRCR(:,51:75),2) nanmean(wDistanceRRCR(:,76:100),2)],'notch','on','symbol','')
ylim([-0.1 1.6])
title('Wistar | Right, Redirect Right')
xticklabels([{'Epoch1'},{'Epoch2'},{'Epoch3'}, {'Epoch4'}])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

subplot(2,3,6)
boxplot([nanmean(pDistanceRRCR(:,1:25),2) nanmean(pDistanceRRCR(:,26:50),2) nanmean(pDistanceRRCR(:,51:75),2) nanmean(pDistanceRRCR(:,76:100),2)],'notch','on','symbol','')
ylim([-0.1 1.6])
title('P rat | Right, Redirect Right')
xticklabels([{'Epoch1'},{'Epoch2'},{'Epoch3'}, {'Epoch4'}])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)


%% Create bar graphs 
% close all
% 
% figure('Units','normalized','Position',[0 0 1 1])
% plot([mean(mean(wDistanceLR(:,1:25),2)) mean(mean(wDistanceLR(:,26:50),2)) mean(mean(wDistanceLR(:,51:75),2)) mean(mean(wDistanceLR(:,76:100),2))],'ko-')
% ylim([0 2])
% figure('Units','normalized','Position',[0 0 1 1])
% plot([mean(mean(pDistanceLR(:,1:25),2)) mean(mean(pDistanceLR(:,26:50),2)) mean(mean(pDistanceLR(:,51:75),2)) mean(mean(pDistanceLR(:,76:100),2))],'ko-')
% ylim([0 2])
% 
% figure('Units','normalized','Position',[0 0 1 1])
% bar([mean(mean(wDistanceLLCR(:,1:25),2)) mean(mean(wDistanceLLCR(:,26:50),2)) mean(mean(wDistanceLLCR(:,51:75),2)) mean(mean(wDistanceLLCR(:,76:100),2))])
% ylim([0 2])
% figure('Units','normalized','Position',[0 0 1 1])
% bar([mean(mean(pDistanceLLCR(:,1:25),2)) mean(mean(pDistanceLLCR(:,26:50),2)) mean(mean(pDistanceLLCR(:,51:75),2)) mean(mean(pDistanceLLCR(:,76:100),2))])
% ylim([0 2])
% 
% figure('Units','normalized','Position',[0 0 1 1])
% bar([mean(mean(wDistanceRRCR(:,1:25),2)) mean(mean(wDistanceRRCR(:,26:50),2)) mean(mean(wDistanceRRCR(:,51:75),2)) mean(mean(wDistanceRRCR(:,76:100),2))])
% ylim([0 2])
% figure('Units','normalized','Position',[0 0 1 1])
% bar([mean(mean(pDistanceRRCR(:,1:25),2)) mean(mean(pDistanceRRCR(:,26:50),2)) mean(mean(pDistanceRRCR(:,51:75),2)) mean(mean(pDistanceRRCR(:,76:100),2))])
% ylim([0 2])

%% Kruksall-Wallis
close all
[p,h,stats] = kruskalwallis([mean(wDistanceLLCR(:,1:10),2) mean(wDistanceLLCR(:,11:50),2) mean(wDistanceLLCR(:,51:90),2) mean(wDistanceLLCR(:,91:100),2)]);
multcompare(stats)
[p,h,stats] = kruskalwallis([mean(pDistanceLLCR(:,1:10),2) mean(pDistanceLLCR(:,11:50),2) mean(pDistanceLLCR(:,51:90),2) mean(pDistanceLLCR(:,91:100),2)]);
multcompare(stats)
[p,h,stats] = kruskalwallis([mean(wDistanceRRCR(:,1:10),2) mean(wDistanceRRCR(:,11:50),2) mean(wDistanceRRCR(:,51:90),2) mean(wDistanceRRCR(:,91:100),2)]);
multcompare(stats)
[p,h,stats] = kruskalwallis([mean(pDistanceRRCR(:,1:10),2) mean(pDistanceRRCR(:,11:50),2) mean(pDistanceRRCR(:,51:90),2) mean(pDistanceRRCR(:,91:100),2)]);
multcompare(stats)

%% Create boxplots for velocity
%
close all

figure('Units','normalized','Position',[0 0 1 1])
subplot(2,3,1)
boxplot([max(wVelLR(:,1:25),[],2) max(wVelLR(:,26:50),[],2) max(wVelLR(:,51:75),[],2) max(wVelLR(:,76:100),[],2)],'notch','on','symbol','')
title('Wistar | Left, Right')
xticklabels([{'Epoch1'},{'Epoch2'},{'Epoch3'}, {'Epoch4'}])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

subplot(2,3,4)
boxplot([max(pVelLR(:,1:25),[],2) max(pVelLR(:,26:50),[],2) max(pVelLR(:,51:75),[],2) max(pVelLR(:,76:100),[],2)],'notch','on','symbol','')
title('P rat | Left, Right')
xticklabels([{'Epoch1'},{'Epoch2'},{'Epoch3'}, {'Epoch4'}])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

subplot(2,3,2)
boxplot([max(wVelLLCR(:,1:25),[],2) max(wVelLLCR(:,26:50),[],2) max(wVelLLCR(:,51:75),[],2) max(wVelLLCR(:,76:100),[],2)],'notch','on','symbol','')
title('Wistar | Left, Redirect Left')
xticklabels([{'Epoch1'},{'Epoch2'},{'Epoch3'}, {'Epoch4'}])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

subplot(2,3,5)
boxplot([max(pVelLLCR(:,1:25),[],2) max(pVelLLCR(:,26:50),[],2) max(pVelLLCR(:,51:75),[],2) max(pVelLLCR(:,76:100),[],2)],'notch','on','symbol','')
title('P rat | Left, Redirect Left')
xticklabels([{'Epoch1'},{'Epoch2'},{'Epoch3'}, {'Epoch4'}])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

subplot(2,3,3)
boxplot([max(wDistanceRRCR(:,1:25),[],2) max(wDistanceRRCR(:,26:50),[],2) max(wDistanceRRCR(:,51:75),[],2) max(wDistanceRRCR(:,76:100),[],2)],'notch','on','symbol','')
title('Wistar | Right, Redirect Right')
xticklabels([{'Epoch1'},{'Epoch2'},{'Epoch3'}, {'Epoch4'}])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

subplot(2,3,6)
boxplot([max(pDistanceRRCR(:,1:25),[],2) max(pDistanceRRCR(:,26:50),[],2) max(pDistanceRRCR(:,51:75),[],2) max(pDistanceRRCR(:,76:100),[],2)],'notch','on','symbol','')
title('P rat | Right, Redirect Right')
xticklabels([{'Epoch1'},{'Epoch2'},{'Epoch3'}, {'Epoch4'}])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
%%%%%%
%%
%%%%%%

close all

figure('Units','normalized','Position',[0 0 1 1])
subplot(2,3,1)
boxplot([min(wVelLR(:,1:25),[],2) min(wVelLR(:,26:50),[],2) min(wVelLR(:,51:75),[],2) min(wVelLR(:,76:100),[],2)],'notch','on','symbol','')
title('Wistar | Left, Right')
xticklabels([{'Epoch1'},{'Epoch2'},{'Epoch3'}, {'Epoch4'}])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

subplot(2,3,4)
boxplot([min(pVelLR(:,1:25),[],2) min(pVelLR(:,26:50),[],2) min(pVelLR(:,51:75),[],2) min(pVelLR(:,76:100),[],2)],'notch','on','symbol','')
title('P rat | Left, Right')
xticklabels([{'Epoch1'},{'Epoch2'},{'Epoch3'}, {'Epoch4'}])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

subplot(2,3,2)
boxplot([min(wVelLLCR(:,1:25),[],2) min(wVelLLCR(:,26:50),[],2) min(wVelLLCR(:,51:75),[],2) min(wVelLLCR(:,76:100),[],2)],'notch','on','symbol','')
title('Wistar | Left, Redirect Left')
xticklabels([{'Epoch1'},{'Epoch2'},{'Epoch3'}, {'Epoch4'}])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

subplot(2,3,5)
boxplot([min(pVelLLCR(:,1:25),[],2) min(pVelLLCR(:,26:50),[],2) min(pVelLLCR(:,51:75),[],2) min(pVelLLCR(:,76:100),[],2)],'notch','on','symbol','')
title('P rat | Left, Redirect Left')
xticklabels([{'Epoch1'},{'Epoch2'},{'Epoch3'}, {'Epoch4'}])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

subplot(2,3,3)
boxplot([min(wDistanceRRCR(:,1:25),[],2) min(wDistanceRRCR(:,26:50),[],2) min(wDistanceRRCR(:,51:75),[],2) min(wDistanceRRCR(:,76:100),[],2)],'notch','on','symbol','')
title('Wistar | Right, Redirect Right')
xticklabels([{'Epoch1'},{'Epoch2'},{'Epoch3'}, {'Epoch4'}])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

subplot(2,3,6)
boxplot([min(pDistanceRRCR(:,1:25),[],2) min(pDistanceRRCR(:,26:50),[],2) min(pDistanceRRCR(:,51:75),[],2) min(pDistanceRRCR(:,76:100),[],2)],'notch','on','symbol','')
title('P rat | Right, Redirect Right')
xticklabels([{'Epoch1'},{'Epoch2'},{'Epoch3'}, {'Epoch4'}])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

%%
%%%%%%

close all

figure('Units','normalized','Position',[0 0 1 1])
subplot(2,3,1)
boxplot([min(wVelLR(:,1:25),[],2) min(wVelLR(:,26:50),[],2) min(wVelLR(:,51:75),[],2) min(wVelLR(:,76:100),[],2)],'notch','on','symbol','')
title('Wistar | Left, Right')
xticklabels([{'Epoch1'},{'Epoch2'},{'Epoch3'}, {'Epoch4'}])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

subplot(2,3,4)
boxplot([min(pVelLR(:,1:25),[],2) min(pVelLR(:,26:50),[],2) min(pVelLR(:,51:75),[],2) min(pVelLR(:,76:100),[],2)],'notch','on','symbol','')
title('P rat | Left, Right')
xticklabels([{'Epoch1'},{'Epoch2'},{'Epoch3'}, {'Epoch4'}])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

subplot(2,3,2)
boxplot([min(wVelLLCR(:,1:25),[],2) min(wVelLLCR(:,26:50),[],2) min(wVelLLCR(:,51:75),[],2) min(wVelLLCR(:,76:100),[],2)],'notch','on','symbol','')
title('Wistar | Left, Redirect Left')
xticklabels([{'Epoch1'},{'Epoch2'},{'Epoch3'}, {'Epoch4'}])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

subplot(2,3,5)
boxplot([min(pVelLLCR(:,1:25),[],2) min(pVelLLCR(:,26:50),[],2) min(pVelLLCR(:,51:75),[],2) min(pVelLLCR(:,76:100),[],2)],'notch','on','symbol','')
title('P rat | Left, Redirect Left')
xticklabels([{'Epoch1'},{'Epoch2'},{'Epoch3'}, {'Epoch4'}])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

subplot(2,3,3)
boxplot([min(wDistanceRRCR(:,1:25),[],2) min(wDistanceRRCR(:,26:50),[],2) min(wDistanceRRCR(:,51:75),[],2) min(wDistanceRRCR(:,76:100),[],2)],'notch','on','symbol','')
title('Wistar | Right, Redirect Right')
xticklabels([{'Epoch1'},{'Epoch2'},{'Epoch3'}, {'Epoch4'}])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

subplot(2,3,6)
boxplot([min(pDistanceRRCR(:,1:25),[],2) min(pDistanceRRCR(:,26:50),[],2) min(pDistanceRRCR(:,51:75),[],2) min(pDistanceRRCR(:,76:100),[],2)],'notch','on','symbol','')
title('P rat | Right, Redirect Right')
xticklabels([{'Epoch1'},{'Epoch2'},{'Epoch3'}, {'Epoch4'}])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)


%% Track Raw Mean Firing Rates That Map Onto PCs that are important for left vs right choices
% 
% First, index coefficients based on loadings
% PCs 2-4 seem to have some directionality to them

% Variable allPCProjections contains coefficients for both Wistars (1) and
% P rats (2)
wPCIPOS = []; wPCINEG = [];

for i = 1:5
    wPCIPOS(i,:) = allPCCoeff{1}(:,i) > 0.01;
    wPCINEG(i,:) = allPCCoeff{1}(:,i) < -0.01;
    pPCIPOS(i,:) = allPCCoeff{2}(:,i) > 0.01;
    pPCINEG(i,:) = allPCCoeff{2}(:,i) < -0.01;

end
totalW = size(allPCCoeff{1},1);
totalP = size(allPCCoeff{2},1);

% Coefficient loadings per PC
figure('Units','normalized','Position',[0 0 1 1])
WNumPosPC = sum(wPCIPOS,2);
WNumNegPC = sum(wPCINEG,2);
WNeither = [totalW - (WNumPosPC + WNumNegPC)];

totalNumW = [WNumPosPC WNumNegPC WNeither];
bar(totalNumW,'stacked')
title('Wistars | Ratio of Positive Versus Negative Loaders')
xlabel('PC Number')
ylabel('Neurons')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

figure('Units','normalized','Position',[0 0 1 1])
PNumPosPC = sum(pPCIPOS,2);
PNumNegPC = sum(pPCINEG,2);
PNeither = [totalP - (PNumPosPC + PNumNegPC)];

totalNumP = [PNumPosPC PNumNegPC PNeither];
bar(totalNumP,'stacked')
title('P rats | Ratio of Postive Versus Negative Loaders')
xlabel('PC Number')
ylabel('Neurons')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

rawData = allPCCoeff{1} * allPCProjections{1}';
rawL = rawData(:,1:100);
rawR = rawData(:,101:end);

figure('Units','normalized','Position',[0 0 1 1])
    
plot(mean(rawL)); hold on; plot(mean(rawR));