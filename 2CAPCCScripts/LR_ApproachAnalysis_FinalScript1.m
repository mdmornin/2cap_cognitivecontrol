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
    parentPath = '/research/dissDat';
    figPath = '/research/dissDat/figs';
end
% LOAD
load([parentPath filesep 'ephysStruct.mat']);
load([parentPath filesep 'trlStruct.mat']);
load([parentPath filesep 'masterTable.mat']);

addpath(genpath([parentPath filesep 'analysisScripts']))
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
allMatI = [wMatLIC pMatLIC; wMatRIC pMatRIC];
% Below contain Congruent only sessino types
allMatL = [wMatLC pMatLC];
allMatR = [wMatRC pMatRC];
conditionIndexes = [length(wMatLC) length(pMatLC)];
conditionLabels = [{'Congruent W'}, {'Congruent P'}];
% Format matrix such that it expands 'in time / observations'
allMat = [allMatL; allMatR];
% allMat = zscore(allMat,[],2);
% allMat = [ephysStruct(WIdx).PSTH_ApproachMean ephysStruct(PIdx).PSTH_ApproachMean];

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
timeVec = [-2:0.1:2.1];


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
saveas(gca,[figPath filesep 'varianceExplained'     '_100ms'],'svg');     saveas(gca,[figPath filesep 'varianceExplained'    '_100ms'],'png')

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
pcSort = 2:4;
colors = [{'#0072BD'} {'#4DBEEE'} {'#77AC30'} {'#D95319'} {'#A2142F'} {'#0072BD'} {'#4DBEEE'} {'#77AC30'} {'#D95319'} {'#A2142F'} {'#0072BD'} {'#4DBEEE'} {'#77AC30'} {'#D95319'} {'#A2142F'} {'#0072BD'} {'#4DBEEE'} {'#77AC30'} {'#D95319'} {'#A2142F'} {'#0072BD'} {'#4DBEEE'} {'#77AC30'} {'#D95319'} {'#A2142F'} {'#0072BD'} {'#4DBEEE'} {'#77AC30'} {'#D95319'} {'#A2142F'} {'#0072BD'} {'#4DBEEE'} {'#77AC30'} {'#D95319'} {'#A2142F'} {'#0072BD'} {'#4DBEEE'} {'#77AC30'} {'#D95319'} {'#A2142F'} {'#0072BD'} {'#4DBEEE'} {'#77AC30'} {'#D95319'} {'#A2142F'} {'#0072BD'} {'#4DBEEE'} {'#77AC30'} {'#D95319'} {'#A2142F'}];
time = [-2:0.1:2];
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
        xline(cueOn,'k--','Approach','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
    end
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    if j == 1; legend([{'PC1'},{'PC2'},{'PC3'},{'PC4'},{'PC5'}],'AutoUpdate','off','location','best'); end

end

saveas(gca,[figPath filesep 'allNeurons_split_top5pcs_CongruentOnly'  num2str(pcSort)  '_100ms'],'svg');     saveas(gca,[figPath filesep 'allNeurons_split_top5pcs_CongruentOnly'  num2str(pcSort)  '_100ms'],'png')

%% Plot 2/3D Representations of PC Spaces
allPerms = nchoosek([1:5],2);

for j = 1:length(allPCProjections)
    for i = 1:length(allPerms)
        figure('Units','normalized','Position',[0 0 1 1])
        plot(allPCProjections{j}(1:timeIdx,allPerms(i,1)),allPCProjections{j}(1:timeIdx,allPerms(i,2)),'LineWidth',3,'Color',colors{i},'LineStyle','-')
        hold on
        plot(allPCProjections{j}(1,allPerms(i,1)),allPCProjections{j}(1,allPerms(i,2)),'ko','MarkerSize',10)

        plot(allPCProjections{j}(timeIdx+1:end,allPerms(i,1)),allPCProjections{j}(timeIdx+1:end,allPerms(i,2)),'LineWidth',3,'Color',colors{i},'LineStyle','--')
        plot(allPCProjections{j}(timeIdx+1:timeIdx+1,allPerms(i,1)),allPCProjections{j}(timeIdx+1:timeIdx+1,allPerms(i,2)),'ko','MarkerSize',10)

        xlabel(['PC' num2str(allPerms(i,1))]);
        ylabel(['PC' num2str(allPerms(i,2))]);

        title([titles{j}]);
        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
        saveas(gca,[figPath filesep '2DPROJNEW_CongruentOnly'  titles{j} num2str(allPerms(i,:))  '_100ms'],'svg');     saveas(gca,[figPath filesep '2DPROJNEW_CongruentOnly' titles{j} num2str(allPerms(i,:))  '_100ms'],'png')
    end
end

% %% Plot 3D Representations the Same Way
% allPerms = nchoosek([1:5],3);
% 
% for j = 1:length(allPCProjections)
%     for i = 1:length(allPerms)
%         figure('Units','normalized','Position',[0 0 1 1])
%         plot3(allPCProjections{j}(1:timeIdx,allPerms(i,1)),allPCProjections{j}(1:timeIdx,allPerms(i,2)),allPCProjections{j}(1:timeIdx,allPerms(i,3)),'LineWidth',3,'Color',colors{i},'LineStyle','-')
%         hold on
%         plot3(allPCProjections{j}(1,allPerms(i,1)),allPCProjections{j}(1,allPerms(i,2)),allPCProjections{j}(1,allPerms(i,3)),'ko','MarkerSize',15,'MarkerFaceColor','k')
% 
%         plot3(allPCProjections{j}(timeIdx+1:end,allPerms(i,1)),allPCProjections{j}(timeIdx+1:end,allPerms(i,2)),allPCProjections{j}(timeIdx+1:end,allPerms(i,3)),'LineWidth',3,'Color',colors{i},'LineStyle','--')
%         plot3(allPCProjections{j}(timeIdx+1:timeIdx+1,allPerms(i,1)),allPCProjections{j}(timeIdx+1:timeIdx+1,allPerms(i,2)),allPCProjections{j}(timeIdx+1:timeIdx+1,allPerms(i,3)),'ko','MarkerSize',15,'MarkerFaceColor','k')
% 
%         xlabel(['PC' num2str(allPerms(i,1))]);
%         ylabel(['PC' num2str(allPerms(i,2))]);
%         zlabel(['PC' num2str(allPerms(i,3))]);
% 
%         title([titles{j}]);
%         set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
%         saveas(gca,[figPath filesep '3DPROJNEW_CongruentOnly'  titles{j} num2str(allPerms(i,:))  '_100ms'],'svg');     saveas(gca,[figPath filesep '3DPROJNEW_CongruentOnly' titles{j} num2str(allPerms(i,:))  '_100ms'],'png')
%     end
% end

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
    numPC = 1:stickPoint;
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
for i = 1:stickPoint
    figure('Units','normalized','Position',[0 0 1 1])

    plot(time,mean(wLProj{i}),'LineWidth',3); hold on; plot(time,mean(wRProj{i}),'LineWidth',3,'LineStyle','--'); plot(time,mean(wLCrct{i}),'LineWidth',3); plot(time,mean(wRCrct{i}),'LineWidth',3,'LineStyle','--');
    title(['Ws | PC Number ' num2str(i)])
    legend([{'Left'},{'Right'},{'Correction to Left'},{'Correction to Right'}],'AutoUpdate','off','location','best');
    xlabel('Time (Bins)')
    ylabel('PC Score')
    xline(0,'k--','Approach','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figPath filesep  'WPCTime_CongruentOnly_PC_' num2str(i)  '_100ms'],'svg');     saveas(gca,[figPath filesep 'WPCTime_CongruentOnly_PC_' num2str(i) '_100ms'],'png')

    % P Rats
    figure('Units','normalized','Position',[0 0 1 1])

    plot(time,mean(pLProj{i}),'LineWidth',3); hold on; plot(time,mean(pRProj{i}),'LineWidth',3,'LineStyle','--'); plot(time,mean(pLCrct{i}),'LineWidth',3); plot(time,mean(pRCrct{i}),'LineWidth',3,'LineStyle','--');
    legend([{'Left'},{'Right'},{'Correction to Left'},{'Correction to Right'}],'AutoUpdate','off','location','best');
    title(['Ps | PC Number ' num2str(i)])
    xlabel('Time (Bins)')
    ylabel('PC Score')
    xline(0,'k--','Approach','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
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
numBootStraps = 500;
distanceMeasurement = 'euclidean';
wDistance = []; pDistance = [];
numDim = stickPoint;
WLeft = []; WRight = []; PLeft = []; PRight = [];
for i = 1:numBootStraps

    wSamp(i) = datasample(1:min(size(wLProj{1},1), size(wRProj{1},1)),1);
    pSamp(i) = datasample(1:min(size(pLProj{1},1), size(pRProj{1},1)),1);
    
    for k = 1:numDim
        WLeft = [WLeft; wLProj{k}(wSamp(i),:)];
        WRight = [WRight; wRProj{k}(wSamp(i),:)];
        PLeft = [PLeft; pLProj{k}(pSamp(i),:)];
        PRight = [PRight; pRProj{k}(pSamp(i),:)];
    end
     
    wDistance(i,:) = diag(pdist2(WLeft', WRight',distanceMeasurement));
    pDistance(i,:) = diag(pdist2(PLeft', PRight',distanceMeasurement));

end
figure('Units','normalized','Position',[0 0 1 1])
shadedErrorBar(time,mean(wDistance),std(wDistance)/sqrt(numBootStraps),'lineprops',{'b-','LineWidth',3});
hold on
shadedErrorBar(time,mean(pDistance),std(pDistance)/sqrt(numBootStraps),'lineprops',{'r-','LineWidth',3});
xlabel('Time (s)')
ylabel(['Distance ' distanceMeasurement])
title(['Distance of Top ' num2str(numDim) ' Principal Components'])
xline(0,'k--','Approach','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');

set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,[figPath filesep 'LR_Distance' '_100ms'],'svg');     saveas(gca,[figPath filesep 'LR_Distance'  '_100ms'],'fig')

figure
boxplot([mean(wDistance,2) mean(pDistance,2)],'notch','on')

%
wDistance = []; pDistance = [];
WLeft = []; WRight = []; PLeft = []; PRight = [];

for i = 1:numBootStraps

    wSamp(i) = datasample(1:min(size(wLProj{1},1), size(wLCrct{1},1)),1);
    pSamp(i) = datasample(1:min(size(pLProj{1},1), size(pLCrct{1},1)),1);
     for k = 1:numDim
        WLeft = [WLeft; wLProj{k}(wSamp(i),:)];
        WRight = [WRight; wLCrct{k}(wSamp(i),:)];
        PLeft = [PLeft; pLProj{k}(pSamp(i),:)];
        PRight = [PRight; pLCrct{k}(pSamp(i),:)];
     end

    wDistance(i,:) = diag(pdist2(WLeft', WRight',distanceMeasurement));
    pDistance(i,:) = diag(pdist2(PLeft', PRight',distanceMeasurement));


end
figure('Units','normalized','Position',[0 0 1 1])
shadedErrorBar(time,mean(wDistance),std(wDistance)/sqrt(numBootStraps),'lineprops',{'b-','LineWidth',3});
hold on
shadedErrorBar(time,mean(pDistance),std(pDistance)/sqrt(numBootStraps),'lineprops',{'r-','LineWidth',3});
xlabel('Time (s)')
ylabel(['Distance ' distanceMeasurement])
xline(0,'k--','Approach','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
title(['Distance of Top ' num2str(numDim) ' Principal Components'])

set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,[figPath filesep 'L_LCorrection_Distance' '_100ms'],'svg');     saveas(gca,[figPath filesep 'L_LCorrection_Distance'  '_100ms'],'fig')

figure
boxplot([mean(wDistance,2) mean(pDistance,2)],'notch','on')

%
wDistance = []; pDistance = [];
WLeft = []; WRight = []; PLeft = []; PRight = [];

for i = 1:numBootStraps

    wSamp(i) = datasample(1:min(size(wRProj{1},1), size(wRCrct{1},1)),1);
    pSamp(i) = datasample(1:min(size(pRProj{1},1), size(pRCrct{1},1)),1);
     
     for k = 1:numDim
        WLeft = [WLeft; wRProj{k}(wSamp(i),:)];
        WRight = [WRight; wRCrct{k}(wSamp(i),:)];
        PLeft = [PLeft; pRProj{k}(pSamp(i),:)];
        PRight = [PRight; pRCrct{k}(pSamp(i),:)];
     end

    wDistance(i,:) = diag(pdist2(WLeft', WRight',distanceMeasurement));
    pDistance(i,:) = diag(pdist2(PLeft', PRight',distanceMeasurement));

end
figure('Units','normalized','Position',[0 0 1 1])
shadedErrorBar(time,mean(wDistance),std(wDistance)/sqrt(numBootStraps),'lineprops',{'b-','LineWidth',3});
hold on
shadedErrorBar(time,mean(pDistance),std(pDistance)/sqrt(numBootStraps),'lineprops',{'r-','LineWidth',3});
xlabel('Time (s)')
ylabel(['Distance ' distanceMeasurement])
xline(0,'k--','Approach','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
% legend({'Wistar'},{'P rats'})
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,[figPath filesep 'R_RCorrection_Distance' '_100ms'],'svg');     saveas(gca,[figPath filesep 'R_RCorrection_Distance'  '_100ms'],'fig')

figure
boxplot([mean(wDistance,2) mean(pDistance,2)],'notch','on')


%% Find Dimensions That Explain Most Distance Between Left and Right Choices
% Of the N dimensions in numDim, iteratively sample each 2D combination
% until the top PCs are found for 
close all
numDim = stickPoint;
allPerms = nchoosek([1:numDim],3);
appTime = 1:20;

numBootStraps = 1000;
distanceMeasurement = 'euclidean';
wDistance = []; pDistance = [];
WLeft = []; WRight = []; PLeft = []; PRight = [];
meanPermValW = []; meanPermValP = [];
wDot = []; pDot = [];

wSamp = [randi(size(wLProj{1},1), 1, numBootStraps); randi(size(wRProj{1},1), 1, numBootStraps); randi(size(wLCrct{1},1), 1, numBootStraps); randi(size(wRCrct{1},1), 1, numBootStraps)];
pSamp = [randi(size(pLProj{1},1), 1, numBootStraps); randi(size(pRProj{1},1), 1, numBootStraps); randi(size(pLCrct{1},1), 1, numBootStraps); randi(size(pRCrct{1},1), 1, numBootStraps)];

for j = 1:length(allPerms)
    for i = 1:numBootStraps    

        WLeft = [wLProj{allPerms(j,1)}(wSamp(1,i),:); wLProj{allPerms(j,2)}(wSamp(1,i),:); wLProj{allPerms(j,3)}(wSamp(1,i),:)];
        WRight = [wRProj{allPerms(j,1)}(wSamp(2,i),:); wRProj{allPerms(j,2)}(wSamp(2,i),:); wRProj{allPerms(j,3)}(wSamp(2,i),:)];
        PLeft = [pLProj{allPerms(j,1)}(pSamp(1,i),:); pLProj{allPerms(j,2)}(pSamp(1,i),:); pLProj{allPerms(j,3)}(pSamp(1,i),:)];
        PRight = [pRProj{allPerms(j,1)}(pSamp(2,i),:); pRProj{allPerms(j,2)}(pSamp(2,i),:); pRProj{allPerms(j,3)}(pSamp(2,i),:)];
        
        WLC = [wLCrct{allPerms(j,1)}(wSamp(3,i),:); wLCrct{allPerms(j,2)}(wSamp(3,i),:); wLCrct{allPerms(j,3)}(wSamp(3,i),:)];
        WRC = [wRCrct{allPerms(j,1)}(wSamp(4,i),:); wRCrct{allPerms(j,2)}(wSamp(4,i),:); wRCrct{allPerms(j,3)}(wSamp(4,i),:)];
        PLC = [pLCrct{allPerms(j,1)}(wSamp(3,i),:); pLCrct{allPerms(j,2)}(wSamp(3,i),:); pLCrct{allPerms(j,3)}(wSamp(3,i),:)];
        PRC = [pRCrct{allPerms(j,1)}(wSamp(4,i),:); pRCrct{allPerms(j,2)}(wSamp(4,i),:); pRCrct{allPerms(j,3)}(wSamp(4,i),:)];
    

        wDistance{j}(i,:) = diag(pdist2(WLeft', WRight',distanceMeasurement));
        pDistance{j}(i,:) = diag(pdist2(PLeft', PRight',distanceMeasurement));

        wL_LC{j}(i,:) = diag(pdist2(WLeft', WLC',distanceMeasurement));
        wL_RC{j}(i,:) = diag(pdist2(WLeft', WRC',distanceMeasurement));
        wR_LC{j}(i,:) = diag(pdist2(WRight', WLC',distanceMeasurement));
        wR_RC{j}(i,:) = diag(pdist2(WRight', WRC',distanceMeasurement));
        wLC_RC{j}(i,:)  = diag(pdist2(WLC', WRC',distanceMeasurement));

        pL_LC{j}(i,:) = diag(pdist2(PLeft', PLC',distanceMeasurement));
        pL_RC{j}(i,:) = diag(pdist2(PLeft', PRC',distanceMeasurement));
        pR_LC{j}(i,:) = diag(pdist2(PRight', PLC',distanceMeasurement));
        pR_RC{j}(i,:) = diag(pdist2(PRight', PRC',distanceMeasurement));
        pLC_RC{j}(i,:)  = diag(pdist2(PLC', PRC',distanceMeasurement));

        wDot(j,i) = dot([WLeft(1,appTime) WLeft(2,appTime) WLeft(3,appTime)],[WRight(1,appTime) WRight(2,appTime) WRight(3,appTime)]);
        pDot(j,i) = dot([PLeft(1,appTime) PLeft(2,appTime) PLeft(3,appTime)],[PRight(1,appTime) PRight(2,appTime) PRight(3,appTime)]);

    end

    meanPermValW(j) = mean(mean(wDistance{j}(:,appTime)));
    meanPermValP(j) = mean(mean(pDistance{j}(:,appTime)));

    meanPermValW_L_LC(j) = mean(mean(wL_LC{j}(:,appTime)));
    meanPermValW_L_RC(j) = mean(mean(wL_RC{j}(:,appTime)));
    meanPermValW_R_LC(j) = mean(mean(wR_LC{j}(:,appTime)));
    meanPermValW_R_RC(j) = mean(mean(wR_RC{j}(:,appTime)));
    meanPermValW_LC_RC(j) = mean(mean(wLC_RC{j}(:,appTime)));

    meanPermValP_L_LC(j) = mean(mean(pL_LC{j}(:,appTime)));
    meanPermValP_L_RC(j) = mean(mean(pL_RC{j}(:,appTime)));
    meanPermValP_R_LC(j) = mean(mean(pR_LC{j}(:,appTime)));
    meanPermValP_R_RC(j) = mean(mean(pR_RC{j}(:,appTime)));
    meanPermValP_LC_RC(j) = mean(mean(pLC_RC{j}(:,appTime)));


end
wDotM = mean(wDot,2);
pDotM = mean(pDot,2);

idxW = []; idxP = [];
[~,idxW] = max(meanPermValW);
[~,idxP] = max(meanPermValP);

[~,idxWDot] = min(abs(meanPermValW));
[~,idxPDot] = min(abs(meanPermValP));

idxW = allPerms(idxW,:);
idxP = allPerms(idxP,:);

% idxW = allPerms(idxWDot,:);
% idxP = allPerms(idxPDot,:);
%%
figure
scatter(wDotM.^2,meanPermValW.^2,'bo'); hold on; scatter(pDotM.^2,meanPermValP.^2,'ro')
figure
scatter(wDotM,meanPermValW,'bo'); hold on; scatter(pDotM,meanPermValP,'ro')

%% Logistic regression
% Starting with just Wistars
LRClass = [zeros(1,41) ones(1,41)]; 
Strain = [zeros(size(wLProj{1},1),1); zeros(size(wRProj{1},1),1); ones(size(pLProj{1},1),1); ones(size(pRProj{1},1),1)]; 
% LRClass(LRClass == 1) = {'Right'};
% LRClass(LRClass == 0) = {'Left'};
% Strain(Strain == 0) = {'Wistar'};
% Strain(Strain == 1) = {'P rat'};

modelspec = 'y ~ x1-x41';
% Left is 0, Right is 1
mdl = fitglm(allMat,LRClass);

Mdl = fitclinear(allMat',LRClass,'ObservationsIn','columns','Solver','sparsa',...
    'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
    struct('AcquisitionFunctionName','expected-improvement-plus'))

betas = Mdl.Beta;

% hB = betas > 0.05;
% lB = betas < 0.05;
% 
% figure
% plot(mean(allMat(1:41,hB),2)); hold on; plot(mean(allMat(1:41,lB),2))
% figure
% plot(mean(allMat(42:end,hB),2)); hold on; plot(mean(allMat(42:end,lB),2))
% 
% rawDataW = allPCCoeff{1} * allPCProjections{1}';
% coeffWMax_Pos = rawDataW(allPCCoeff{1}(:,j) > 0, :);

wist_Betas = betas(1:size(wMatLC,2));
prat_Betas = betas(size(wMatLC,2)+1:end);

figure;
histogram(wist_Betas); hold on; histogram(prat_Betas)

%% Distance distribution
f = figure('Units','normalized','Position',[0 0 1 1]);
subplot(2,1,1)
h = histogram(meanPermValW,'FaceColor','blue','Normalization','pdf'); hold on; histogram(meanPermValP,'FaceColor','red','Normalization','pdf')
legend([{'Wistars'},{'P rats'}])
xlabel('Mean of Maximum Distance')
ylabel('Bin Counts')

set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
subplot(2,1,2)
b = boxchart([meanPermValW' meanPermValP'],'notch','on');
b.BoxWidth = 0.5;
% b.BoxLineColor = [0 0 0];
b.LineWidth = 3;
b.Orientation = 'horizontal';
yticklabels = {'Wistars', 'P rats'};
ylabel('Group')
xlabel('Mean of Maximum Distance')

set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

ax = findobj(f,'Type','Axes');
ax(1).XLim = ax(2).XLim;
ax(1).YTickLabel = {'Wistar','P rats'};

saveas(gca,[figPath filesep 'MaxDistHist'  '_100ms'],'svg');     saveas(gca,[figPath filesep 'MaxDistHist'  '_100ms'],'fig')

%% Distance distribution with corrections 
figure('Units','normalized','Position',[0 0 1 1]);
boxData = [meanPermValW; meanPermValW_L_LC; meanPermValW_L_RC; meanPermValW_LC_RC; meanPermValW_R_LC; meanPermValW_R_RC; ... 
    meanPermValP; meanPermValP_L_LC; meanPermValP_L_RC; meanPermValP_LC_RC; meanPermValP_R_LC; meanPermValP_R_RC];

boxchart(boxData','notch','on')
xlabel('Group')
ylabel(['Distance (' distanceMeasurement ')'])
xticklabels({'WLR','WL_LC','WL_RC','WLC_RC','WR_LC','WR_RC','PLR','PL_LC','PL_RC','PLC_RC','PR_LC','PR_RC'})

set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

%% Dot product distribution
f = figure('Units','normalized','Position',[0 0 1 1]);
subplot(2,1,1)
h = histogram(wDotM,'FaceColor','blue','Normalization','pdf'); hold on; histogram(pDotM,'FaceColor','red','Normalization','pdf')
legend([{'Wistars'},{'P rats'}])
xlabel('Mean of Maximum Distance')
ylabel('Bin Counts')

set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
subplot(2,1,2)
b = boxchart([wDotM pDotM],'notch','on');
b.BoxWidth = 0.5;
% b.BoxLineColor = [0 0 0];
b.LineWidth = 3;
b.Orientation = 'horizontal';
yticklabels = {'Wistars', 'P rats'};
ylabel('Group')
xlabel('Dot Product')

set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

ax = findobj(f,'Type','Axes');
ax(1).XLim = ax(2).XLim;
ax(1).YTickLabel = {'Wistar','P rats'};

saveas(gca,[figPath filesep 'DotProdHist'  '_100ms'],'svg');     saveas(gca,[figPath filesep 'DotProdHist'  '_100ms'],'fig')
%%
for i = 1
    figure('Units','normalized','Position',[0 0 1 1])
    plot3(mean(wLProj{idxW(i,1)}),mean(wLProj{idxW(i,2)}),mean(wLProj{idxW(i,3)}),'LineWidth',3,'Color',colors{i},'LineStyle','-')
    hold on
    plot3(mean(wRProj{idxW(i,1)}),mean(wRProj{idxW(i,2)}),mean(wRProj{idxW(i,3)}),'LineWidth',3,'Color',colors{i},'LineStyle','--')
    plot3(mean(wRProj{idxW(i,1)}(:,1)),mean(wRProj{idxW(i,2)}(:,1)),mean(wRProj{idxW(i,3)}(:,1)),'ko','MarkerSize',15,'MarkerFaceColor','k')
    plot3(mean(wLProj{idxW(i,1)}(:,1)),mean(wLProj{idxW(i,2)}(:,1)),mean(wLProj{idxW(i,3)}(:,1)),'ko','MarkerSize',15,'MarkerFaceColor','k')

    xlabel(['PC' num2str(idxW(i,1))]);
    ylabel(['PC' num2str(idxW(i,2))]);
    zlabel(['PC' num2str(idxW(i,3))]);

    title('Wistar Congruent')
    legend([{'Left'},{'Right'}])
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figPath filesep 'MaxDistW'  num2str(allPerms(i,:))  '_100ms'],'svg');     saveas(gca,[figPath filesep 'MaxDistW' num2str(allPerms(i,:))  '_100ms'],'fig')
end

%
for i = 1
    figure('Units','normalized','Position',[0 0 1 1])
    plot3(mean(pLProj{idxP(i,1)}),mean(pLProj{idxP(i,2)}),mean(pLProj{idxP(i,3)}),'LineWidth',3,'Color',colors{i},'LineStyle','-')
    hold on
    plot3(mean(pRProj{idxP(i,1)}),mean(pRProj{idxP(i,2)}),mean(pRProj{idxP(i,3)}),'LineWidth',3,'Color',colors{i},'LineStyle','--')
    plot3(mean(pRProj{idxP(i,1)}(:,1)),mean(pRProj{idxP(i,2)}(:,1)),mean(pRProj{idxP(i,3)}(:,1)),'ko','MarkerSize',15,'MarkerFaceColor','k')
    plot3(mean(pLProj{idxP(i,1)}(:,1)),mean(pLProj{idxP(i,2)}(:,1)),mean(pLProj{idxP(i,3)}(:,1)),'ko','MarkerSize',15,'MarkerFaceColor','k')

    xlabel(['PC' num2str(idxP(i,1))]);
    ylabel(['PC' num2str(idxP(i,2))]);
    zlabel(['PC' num2str(idxP(i,3))]);
    title('P rats Congruent')
    legend([{'Left'},{'Right'}])

    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figPath filesep 'MaxDistP'  num2str(allPerms(i,:))  '_100ms'],'svg');     saveas(gca,[figPath filesep 'MaxDistP' num2str(allPerms(i,:))  '_100ms'],'fig')
end

%%
i = 1;
figure('Units','normalized','Position',[0 0 1 1])
shadedErrorBar(time,mean(wLProj{idxW(i)}),std(wLProj{idxW(i)})/sqrt(length(wLProj{idxW(i)})),'lineprops',{'b-','LineWidth',3}); hold on;  shadedErrorBar(time,mean(wRProj{idxW(i)}),std(wRProj{idxW(i)})/sqrt(length(wRProj{idxW(i)})),'lineprops',{'r-','LineWidth',3});
shadedErrorBar(time,mean(wLCrct{idxW(i)}),std(wLCrct{idxW(i)})/sqrt(length(wLCrct{idxW(i)})),'lineprops',{'b--','LineWidth',3}); hold on;  shadedErrorBar(time,mean(wRCrct{idxW(i)}),std(wRCrct{idxW(i)})/sqrt(length(wRCrct{idxW(i)})),'lineprops',{'r--','LineWidth',3});
legend([{'Left'},{'Right'},{'Redirect to Left'},{'Redirect to Right'}])
xlabel('Time (s)')
ylabel('PC Score')
title('Congruent Wistar | Trial Specific Projected Scores')
xline(cueOn,'k--','Approach','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,[figPath filesep 'MaxDistW_Score'  num2str(i) '_' num2str(allPerms(i,:))  '_100ms'],'svg');     saveas(gca,[figPath filesep 'MaxDistW_Score' num2str(i) '_' num2str(allPerms(i,:))  '_100ms'],'fig')
%%
i = 1;
figure('Units','normalized','Position',[0 0 1 1])
shadedErrorBar(time,mean(pLProj{idxW(i,1)}),std(pLProj{idxW(i,1)})/sqrt(length(pLProj{idxW(i,1)})),'lineprops',{'b-','LineWidth',3}); hold on;  shadedErrorBar(time,mean(pRProj{idxW(i,1)}),std(pRProj{idxW(i,1)})/sqrt(length(pRProj{idxW(i,1)})),'lineprops',{'r-','LineWidth',3});
shadedErrorBar(time,mean(pLCrct{idxW(i,1)}),std(pLCrct{idxW(i,1)})/sqrt(length(pLCrct{idxW(i,1)})),'lineprops',{'b--','LineWidth',3}); hold on;  shadedErrorBar(time,mean(pRCrct{idxW(i,1)}),std(pRCrct{idxW(i,1)})/sqrt(length(pRCrct{idxW(i,1)})),'lineprops',{'r--','LineWidth',3});
legend([{'Left'},{'Right'},{'Redirect to Left'},{'Redirect to Right'}])
xlabel('Time (s)')
ylabel('PC Score')
title('Congruent P Rats | Trial Specific Projected Scores')
xline(cueOn,'k--','Approach','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,[figPath filesep 'MaxDistP_Score'  num2str(allPerms(i,:))  '_100ms'],'svg');     saveas(gca,[figPath filesep 'MaxDistP_Score' num2str(allPerms(i,:))  '_100ms'],'fig')


%% Now that dimenions that produce maximal distance are identified, we can now look at the behavior of neurons that attach to these dims 
% rawDataW = allPCCoeff{1} * allPCProjections{1}';
% rawDataP = allPCCoeff{2} * allPCProjections{2}';

rawDataW = allMat(:,wConIdx)'; 
rawDataP = allMat(:,pConIdx)';

% close all
for j = 1:stickPoint

        coeffWMax_Pos = rawDataW(allPCCoeff{1}(:,j) > 0, :);
        coeffWMax_Neg = rawDataW(allPCCoeff{1}(:,j) < 0, :);
        
        figure('Units','normalized','Position',[0 0 1 1])
        subplot(2,1,1)
        shadedErrorBar(time,mean(coeffWMax_Pos(:,1:41)),std(coeffWMax_Pos(:,1:41))/sqrt(length(coeffWMax_Pos(:,1:41))),'lineprops',{'b-','LineWidth',3}); hold on; shadedErrorBar(time,mean(coeffWMax_Pos(:,42:end)),std(coeffWMax_Pos(:,42:end))/sqrt(length(coeffWMax_Pos(:,42:end))),'lineprops',{'r-','LineWidth',3}); 
        xlabel('Time (s)')
        ylabel('Mean Firing Rate (Hz)')
        xline(cueOn,'k--','Approach','LineWidth',3,'FontSize',16,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        ylim([0.5 1])
        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
        title(['Wistar | Positive Coefficients for PC ' num2str(j)])
    
        subplot(2,1,2)
        shadedErrorBar(time,mean(coeffWMax_Neg(:,1:41)),std(coeffWMax_Neg(:,1:41))/sqrt(length(coeffWMax_Neg(:,1:41))),'lineprops',{'b-','LineWidth',3}); hold on; shadedErrorBar(time,mean(coeffWMax_Neg(:,42:end)),std(coeffWMax_Neg(:,42:end))/sqrt(length(coeffWMax_Neg(:,42:end))),'lineprops',{'r-','LineWidth',3}); 
        xlabel('Time (s)')
        ylabel('Mean Firing Rate (Hz)')
        title(['Wistar | Negative Coefficients for PC ' num2str(j)])
        xline(cueOn,'k--','Approach','LineWidth',3,'FontSize',16,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        ylim([0.5 1])

        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
        saveas(gca,[figPath filesep 'MaxDistW_FR_POS_PCNUM_'  num2str(j)  '_100ms'],'svg');     saveas(gca,[figPath filesep 'MaxDistW_FR_POS_PCNUM_' num2str(j) '_100ms'],'fig')
    
    
        coeffPMax_Pos = rawDataP(allPCCoeff{2}(:,j) > 0, :);
        coeffPMax_Neg = rawDataP(allPCCoeff{2}(:,j) < 0, :);
        
        figure('Units','normalized','Position',[0 0 1 1])
        subplot(2,1,1)
        shadedErrorBar(time,mean(coeffPMax_Pos(:,1:41)),std(coeffPMax_Pos(:,1:41))/sqrt(length(coeffPMax_Pos(:,1:41))),'lineprops',{'b-','LineWidth',3}); hold on; shadedErrorBar(time,mean(coeffPMax_Pos(:,42:end)),std(coeffPMax_Pos(:,42:end))/sqrt(length(coeffPMax_Pos(:,42:end))),'lineprops',{'r-','LineWidth',3}); 
        xlabel('Time (s)')
        ylabel('Mean Firing Rate (Hz)')
        title(['P rats | Positive Coefficients for PC ' num2str(j)])
        xline(cueOn,'k--','Approach','LineWidth',3,'FontSize',16,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        ylim([0.5 1])

        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
        subplot(2,1,2)
        shadedErrorBar(time,mean(coeffPMax_Neg(:,1:41)),std(coeffPMax_Neg(:,1:41))/sqrt(length(coeffPMax_Neg(:,1:41))),'lineprops',{'b-','LineWidth',3}); hold on; shadedErrorBar(time,mean(coeffPMax_Neg(:,42:end)),std(coeffPMax_Neg(:,42:end))/sqrt(length(coeffPMax_Neg(:,42:end))),'lineprops',{'r-','LineWidth',3}); 
        xlabel('Time (s)')
        ylabel('Mean Firing Rate (Hz)')
        title(['P rats | Negative Coefficients for PC ' num2str(j)])
        xline(cueOn,'k--','Approach','LineWidth',3,'FontSize',16,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        ylim([0.5 1])

        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
        saveas(gca,[figPath filesep 'MaxDistP_FR_NEG_PCNUM'  num2str(j)  '_100ms'],'svg');     saveas(gca,[figPath filesep 'MaxDistP_FR_NEG_PCNUM' num2str(j)  '_100ms'],'fig')
    
    
end