%%
%% Load and check data
clear all; close all % Refresh workspace
rng(1)
% PATHS
% For brain3 Linux Machine:
if ispc
    parentPath = 'F:/dissDat';
    figPath = 'F:/dissDat/newfigs';
else
    parentPath = '/research/dissDat';
    figPath = '/research/dissDat/figs';
end
% LOAD
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
% Pull ephys and trial data, create PSTH
    pBin = 0.1; % 33 ms bins for fr variable, match video framerate
    PSTH_Approach = [];
    if ~isfield(ephysStruct,'PSTH_Approach')
        
        for i = 1:length(trlStruct)
    
            % Load trial times from trlStruct
            trialTimes = trlStruct(i).trialTimes(1:48);
            [trialTimes, trialTimesIdx] = sort(trialTimes);


            approachVar = trlStruct(i).approachChangePoints;           
            sortedApproachVar = approachVar(trialTimesIdx);

%             approachVar = min(trlStruct(i).approach(:,1:2,2),[],2);
%             sortedApproachVar = approachVar(trialTimesIdx);

            approachIdx = trlStruct(i).approach(1:48,1,1) == 1 | trlStruct(i).approach(1:48,2,1) == 1;
            approachIdx = approachIdx(trialTimesIdx);
            sortedApproachVar(approachIdx == 0) = NaN;
            sortedApproachVar(approachIdx == 0) = randsample(sortedApproachVar(~isnan(sortedApproachVar)),1);

            if startsWith(masterTbl.SessionType(i),'Regular')            
                correctIndex = trlStruct(i).approach(1:48,1,1) == 1 & trlStruct(i).approach(1:48,2,1) == 0;
                incorrectIndex = trlStruct(i).approach(1:48,2,1) == 1;
                
                correctIndexL = trlStruct(i).approach(1:24,1,1) == 1 & trlStruct(i).approach(1:24,2,1) == 0;
                correctIndexR = trlStruct(i).approach(25:48,1,1) == 1 & trlStruct(i).approach(25:48,2,1) == 0;

             
            elseif startsWith(masterTbl.SessionType(i),'Reversal')
                correctIndex = trlStruct(i).approach(1:48,1,1) == 0 & trlStruct(i).approach(1:48,2,1) == 1;
                incorrectIndex = trlStruct(i).approach(1:48,1,1) == 1;
                
                correctIndexL = trlStruct(i).approach(25:48,1,1) == 0 & trlStruct(i).approach(25:48,2,1) == 1;
                correctIndexR = trlStruct(i).approach(1:24,1,1) == 0 & trlStruct(i).approach(1:24,2,1) == 1;

            end
            approachIndex = trlStruct(i).approach(1:48,1,1) == 1 | trlStruct(i).approach(1:48,2,1) == 1;
            noApproachIndex = approachIndex(approachIndex == 0);

            correctIndex = correctIndex(trialTimesIdx);
            incorrectIndex = incorrectIndex(trialTimesIdx);

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
            frTime = (1:length(fr))/10;   % Denominator corresponds with pBin
    
            dlcTimestamps = trlStruct(i).bodyCoords{1};
            dlcTimestamps = dlcTimestamps(dlcTimestamps >= 0);
            dlcTimestamps = (dlcTimestamps);
            trialTimes = trialTimes + sortedApproachVar;
            trialTimes = (trialTimes);
        
        for k = 1:length(trialTimes(1:48))
            [~,tMin] = min(abs((trialTimes(k) - 5) - dlcTimestamps));
            [~,tMax] = min(abs((trialTimes(k) + 2) - dlcTimestamps)); 

            [~,tMin_fr] = min(abs((dlcTimestamps(tMin)) - frTime));
            [~,tMax_fr] = min(abs((dlcTimestamps(tMax)) - frTime));
            
            trlSize(k) = size(fr(tMin_fr:tMax_fr,:),1);

            if size(fr(tMin_fr:tMax_fr,:),1) == 71
                PSTH(k,:,:) = fr(tMin_fr:tMax_fr,:);
            elseif size(fr(tMin_fr:tMax_fr,:),1) == 72
                PSTH(k,:,:) = fr(tMin_fr:tMax_fr-1,:);
            elseif size(fr(tMin_fr:tMax_fr,:),1) == 70
                PSTH(k,:,:) = fr(tMin_fr:tMax_fr+1,:);          
            end
            
                
               
        end

        

        % Remove low firing rate neurons. Not sure what the best value is
        % here.
        minFR = 0.1;
        if ~isempty(PSTH)            
            meanFR = squeeze(mean(mean(PSTH,1),2)); 
            lowFR = meanFR > minFR; % Find values greater than minFR
            PSTH = PSTH(:,:,lowFR);
            PSTH = zscore(PSTH,[],1);
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
        ephysStruct(i).PSTH_LeftMean= squeeze(mean(ephysStruct(i).PSTH_Left,1));
        ephysStruct(i).PSTH_RightMean = squeeze(mean(ephysStruct(i).PSTH_Right,1));
        PSTH_Approach = [];
        PSTH_NoApproach = [];
        PSTH = [];

        

    end

end

%%
ephysStruct(19).neuronNumApproach = 0;
ephysStruct(19).maxFR = [];
ephysStruct(19).PSTH_Correct = [];
ephysStruct(19).PSTH_Incorrect = [];
ephysStruct(19).PSTH_ApproachMean = [];
ephysStruct(19).PSTH_NoApproachMean = [];

ephysStruct(19).PSTH_CorrectMean = [];
ephysStruct(19).PSTH_IncorrectMean = [];
ephysStruct(19).PSTH_Approach = [];
ephysStruct(19).PSTH = [];

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
pMatLC = [ephysStruct(regPIdx).PSTH_ApproachMean];
wMatLC = [ephysStruct(regWIdx).PSTH_ApproachMean];

pMatRC = [ephysStruct(regPIdx).PSTH_NoApproachMean];
wMatRC = [ephysStruct(regWIdx).PSTH_NoApproachMean];
% Incongruent
pMatLIC = [ephysStruct(revPIdx).PSTH_ApproachMean];
wMatLIC = [ephysStruct(revWIdx).PSTH_ApproachMean];

pMatRIC = [ephysStruct(revPIdx).PSTH_NoApproachMean];
wMatRIC = [ephysStruct(revWIdx).PSTH_NoApproachMean];

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
%% Basic Figures First
% Plot mean of allMat
figure('Units','normalized','Position',[0 0 1 1])
plot(mean(allMat(1:71,:),2))
hold on
plot(mean(allMat(72:end,:),2))

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
xline(51,'--','LApproach')
xline(121,'--','RApproach')
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
%% Plot Separated PCs Over Time
pcSort = 1:5;
colors = [{'#0072BD'} {'#4DBEEE'} {'#77AC30'} {'#D95319'} {'#A2142F'} {'#0072BD'} {'#4DBEEE'} {'#77AC30'} {'#D95319'} {'#A2142F'} {'#0072BD'} {'#4DBEEE'} {'#77AC30'} {'#D95319'} {'#A2142F'} {'#0072BD'} {'#4DBEEE'} {'#77AC30'} {'#D95319'} {'#A2142F'} {'#0072BD'} {'#4DBEEE'} {'#77AC30'} {'#D95319'} {'#A2142F'} {'#0072BD'} {'#4DBEEE'} {'#77AC30'} {'#D95319'} {'#A2142F'} {'#0072BD'} {'#4DBEEE'} {'#77AC30'} {'#D95319'} {'#A2142F'} {'#0072BD'} {'#4DBEEE'} {'#77AC30'} {'#D95319'} {'#A2142F'} {'#0072BD'} {'#4DBEEE'} {'#77AC30'} {'#D95319'} {'#A2142F'} {'#0072BD'} {'#4DBEEE'} {'#77AC30'} {'#D95319'} {'#A2142F'}];
time = [-5:0.1:2];
timeIdx = length(wConScore)/2;
cueOn = 0;

figure('Units','normalized','Position',[0 0 1 1])

for j = 1:length(allPCProjections)
    for i=1:length(pcSort)
        
        subplot(1,2,j);
        plot(time,allPCProjections{j}(1:timeIdx,pcSort(i)),'LineWidth',3,'Color',colors{i},'LineStyle','-')

        hold on
        plot(time,allPCProjections{j}(timeIdx+1:end,pcSort(i)),'LineWidth',3,'Color',colors{i},'LineStyle','--')
        ylim([-3 3])
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

%% Plot 3D Representations the Same Way
allPerms = nchoosek([1:5],3);

for j = 1:length(allPCProjections)
    for i = 1:length(allPerms)
        figure('Units','normalized','Position',[0 0 1 1])
        plot3(allPCProjections{j}(1:timeIdx,allPerms(i,1)),allPCProjections{j}(1:timeIdx,allPerms(i,2)),allPCProjections{j}(1:timeIdx,allPerms(i,3)),'LineWidth',3,'Color',colors{i},'LineStyle','-')
        hold on
        plot3(allPCProjections{j}(1,allPerms(i,1)),allPCProjections{j}(1,allPerms(i,2)),allPCProjections{j}(1,allPerms(i,3)),'ko','MarkerSize',15,'MarkerFaceColor','k')

        plot3(allPCProjections{j}(timeIdx+1:end,allPerms(i,1)),allPCProjections{j}(timeIdx+1:end,allPerms(i,2)),allPCProjections{j}(timeIdx+1:end,allPerms(i,3)),'LineWidth',3,'Color',colors{i},'LineStyle','--')
        plot3(allPCProjections{j}(timeIdx+1:timeIdx+1,allPerms(i,1)),allPCProjections{j}(timeIdx+1:timeIdx+1,allPerms(i,2)),allPCProjections{j}(timeIdx+1:timeIdx+1,allPerms(i,3)),'ko','MarkerSize',15,'MarkerFaceColor','k')

        xlabel(['PC' num2str(allPerms(i,1))]);
        ylabel(['PC' num2str(allPerms(i,2))]);
        zlabel(['PC' num2str(allPerms(i,3))]);

        title([titles{j}]);
        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
        saveas(gca,[figPath filesep '3DPROJNEW_CongruentOnly'  titles{j} num2str(allPerms(i,:))  '_100ms'],'svg');     saveas(gca,[figPath filesep '3DPROJNEW_CongruentOnly' titles{j} num2str(allPerms(i,:))  '_100ms'],'png')
    end
end

%% Project PCs onto Trials 
% One way of doing this is to take the coefficients from each session and
% project those onto the full PSTH for that session. From there, we can
% take the average of all trial types for that sessions and concatenate
% those into one variable per genotype / session type
%
%
