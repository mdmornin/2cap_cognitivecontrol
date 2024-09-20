%% Basic script for analysis of 2CAP data
% 04/14/2022 -- MDM
% Pulls data from ephysStruct, trlStruct, masterTbl
% 
% - - - - - - - - 04/22/2022 - - - - - - - - 
% Good generic functionality, can handle and graph outputs from basic PCA
% well. There is some minor mis-alignment that is potentially caused by the
% gaussconv function that I utilized. This will need to be addressed at
% some point. 
%
% -- Investigated difference between gaussianconv and not, did not solve
% timestamps issue. Likely has to do with OE/DLC timestamp differences
% --
% 
% Next steps include... 
% 1. Separating behavior into correct and incorrect
% approach trials. From there, it should be a 2x2x2 analysis wherein there
% are correct and incorrect approaches for both Wistar and P rats in
% congruent and incongruent sessions.
%
% This should help the PCA explain more variance.
% 
% Additionally, split PCA into epochs? Investigate positive versus negative
% loaders?
%
% - - - - - - - - 04/25/2022 - - - - - - - - 
% Starting today by finishing behavior x ephys. Initially I was confused
% how to implement this since the matrix would have 2 unequal dimensions.
% This leaves the option to average over the neurons or average over the
% trials. Since we are interested in the neurons' behavior, I will take the
% average over trials for each session's correct and incorrect approaches
% that occur within the first 15 trials. 

% - - - - - - - - 05/03/2022 - - - - - - - - 
% Today I am going to lengthen the window of firing rates from 33ms to
% 100ms. This should do several things: increasing stability of PCs, mean
% values, and minimizing 0 values. This value also corresponds with what
% NMT used in his paper which will allow for better comparison. NMT also
% used a gaussian convlution, but I am unsure if the details are completely
% the same. 
%
% - - - - - - - - 05/06/2022 - - - - - - - - 
% Make an effort to replicate aspects of the way NMT did his PCA.
% Additionally, expand functionality of PSTH loop to include several
% timepoints. Eventually, will need 2 major functions: one to create or
% load PSTHs at a given timepoint input and one to facilitate replicable
% analysis of PCAs.
% For future analyses, flexiblity in data structures and rapidly and
% readily obtaining the information needed for dPCA, RNNs, and so forth
% will be critical. Seek ways to do this
%% Load and check data
clear all; close all % Refresh workspace
% PATHS
% For brain3 Linux Machine:
if ispc
    parentPath = 'F:/dissDat';
    figPath = 'F:/dissDat/figs';
else
    parentPath = '/research/dissDat';
    figPath = '/research/dissDat/figs';
end
% LOAD
load([parentPath filesep 'ephysStruct.mat']);
load([parentPath filesep 'trlStruct.mat']);
load([parentPath filesep 'masterTable.mat']);

addpath(genpath([parentPath filesep 'analysisScripts']))

    %% Pull ephys and trial data, create PSTH
    pBin = 0.1; % 33 ms bins for fr variable, match video framerate
    PSTH = [];
    if ~isfield(ephysStruct,'PSTH')
        
        for i = 1:length(trlStruct)
    
            % Load trial times from trlStruct
            trialTimes = trlStruct(i).trialTimes(1:48);
            [trialTimes, trialTimesIdx] = sort(trialTimes);
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
            trialTimes = (trialTimes);
        
        for k = 1:length(trialTimes(1:15))
            [~,tMin] = min(abs((trialTimes(k) - 1) - dlcTimestamps));
            [~,tMax] = min(abs((trialTimes(k) + 15) - dlcTimestamps)); 

            [~,tMin_fr] = min(abs((dlcTimestamps(tMin)) - frTime));
            [~,tMax_fr] = min(abs((dlcTimestamps(tMax)) - frTime));
            
            trlSize(k) = size(fr(tMin_fr:tMax_fr,:),1);

            if size(fr(tMin_fr:tMax_fr,:),1) == 161
                PSTH(k,:,:) = fr(tMin_fr:tMax_fr,:);
            elseif size(fr(tMin_fr:tMax_fr,:),1) == 162
                PSTH(k,:,:) = fr(tMin_fr:tMax_fr-1,:);
            elseif size(fr(tMin_fr:tMax_fr,:),1) == 160
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
        end
        
        ephysStruct(i).PSTH = PSTH;
        if ~isempty(PSTH)
            ephysStruct(i).neuronNum = size(PSTH,3);
        else
            ephysStruct(i).neuronNum = 0;
        end
        ephysStruct(i).maxFR = max(fr);
        
        PSTH = [];

    end

end


%% Set indexing based on Congruent versus Incongruent Sessions
regPIdx = startsWith(masterTbl.SessionType,'Regular') & startsWith(masterTbl.Strain,'P');
revPIdx = startsWith(masterTbl.SessionType,'Reversal') & startsWith(masterTbl.Strain,'P');
    
regWIdx = startsWith(masterTbl.SessionType,'Regular') & startsWith(masterTbl.Strain,'W');
revWIdx = startsWith(masterTbl.SessionType,'Reversal') & startsWith(masterTbl.Strain,'W');


%% Index 2x2 Wistar, P | Congruent, Incongruent 
pCon = cat(3,ephysStruct(regPIdx).PSTH);
pInc = cat(3,ephysStruct(revPIdx).PSTH);
wCon = cat(3,ephysStruct(regWIdx).PSTH);
wInc = cat(3,ephysStruct(revWIdx).PSTH);

%% Remove Neurons that Do Not Fire
%% Timepoints and other parameters
numTrials = 48;


sipDescent = 6;
sipAscent = 14;
cueOn = 1;

LeftTrials = 1:24;
RightTrials = 25:48;
colors = [{'#EDB120'}, {'#7E2F8E'}, {'#0072BD'},  {'#D95319'} ];
%% Quick Plot

plot(time,smooth(mean(mean(zscore(pCon(1:15,:,:),1,2),1),3),5))
hold on
plot(time,smooth(mean(mean(zscore(pInc(1:15,:,:),1,2),1),3),5))
plot(time,smooth(mean(mean(zscore(wCon(1:15,:,:),1,2),1),3),5))
plot(time,smooth(mean(mean(zscore(wInc(1:15,:,:),1,2),1),3),5))

xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');

%% PCA Preprocessing | Separating P's and Wistars, Congruent and Incongruent

pConPre = squeeze(mean(pCon(1:15,:,:),1));
pIncPre = squeeze(mean(pInc(1:15,:,:),1));
wConPre = squeeze(mean(wCon(1:15,:,:),1));
wIncPre = squeeze(mean(wInc(1:15,:,:),1));

[pConCoeff,pConScore,~,~,pConExplained] = pca(([pConPre]));
[pIncCoeff,pIncScore,~,~,pIncExplained] = pca(([pIncPre]));

[wConCoeff,wConScore,~,~,wConExplained] = pca(([wConPre]));
[wIncCoeff,wIncScore,~,~,wIncExplained] = pca(([wIncPre]));

% Centralize all output variables into a few cell arrays
allCoeff = [{pConCoeff} {pIncCoeff} {wConCoeff} {wIncCoeff}]; 
allScore = [{pConScore} {pIncScore} {wConScore} {wIncScore}];
allExplained = [{pConExplained} {pIncExplained} {wConExplained} {wIncExplained}];
titles = [{'Congruent P'} {'Incongruent P'} {'Congruent W'} {'Incongruent W'}];



%% Quick PCA Plot % 1st PC
pcSort = 1:3;
figure('Units','normalized','Position',[0 0 1 1])
time = (1:size(pCon,2))/10;

for j = 1:length(pcSort)
    for i=1:length(allCoeff)
        
        subplot(max(pcSort),1,j);
        plot(time,allScore{i}(:,j),'LineWidth',3,'Color',colors{i})
        hold on
        xlabel('Time (s)');
        title(['PC Num: ' num2str(pcSort(j))]);
        xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',10,'FontWeight','bold','LineWidth',4)
         
    end
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

end
saveas(gca,[figPath filesep 'pcScore_numPCs_' num2str(pcSort) '_100ms'],'svg');     saveas(gca,[figPath filesep 'pcScore_numPCs_' num2str(pcSort) '_100ms'],'png')

%% Sort Raw Values According to Coefficients

for j = 1:length(allCoeff)
    figure('Units','normalized','Position',[0 0 1 1])
    for i=1:length(pcSort)
        
        [b,k]=sort(allCoeff{j}(:,pcSort(i)));

        subplot(1,max(pcSort),i);
        centData = allCoeff{j} * allScore{j}';
        imagesc([0 max(time)], [1 length(b)],centData(k,:),[-1 1]);

        xlabel('Time (s)');
        if i==1; ylabel('Neuron # (Sorted by PC)');end
        xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        title(['PC Num: ' num2str(pcSort(i)) ',ExplVar=' num2str(sum(allExplained{j}(1:i)),'%.1f') '%']);
        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
         
    end
    sgtitle(titles{j})
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figPath filesep 'pcProj_' titles{j} '_100ms'],'svg');     saveas(gca,[figPath filesep 'pcProj_' titles{j} '_100ms'],'png')

end

%% PC Space Projections
figure('Units','normalized','Position',[0 0 1 1])
for j = 1:length(allCoeff)
    plot(allScore{j}(:,1),allScore{j}(:,2),'Color',colors{j},'LineWidth',3);
    hold on
    plot(allScore{j}(1,1),allScore{j}(1,2),'o','MarkerSize',10,'MarkerFaceColor',colors{j},'MarkerEdgeColor','k')
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)        
    

end
    saveas(gca,[figPath filesep 'pcSpaceProj_' '_100ms'],'svg');     saveas(gca,[figPath filesep 'pcSpaceProj_' '_100ms'],'png')
%%
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
%% Split electrophysiology data based on behaviors
for i = 1:length(trlStruct)
    if startsWith(masterTbl.SessionType{i},['Regular'])

            currentTrialTimes = trlStruct(i).trialTimes;
            [~,sortedIndex] = sort(currentTrialTimes(1:48));
            currentApproachVector = trlStruct(i).approach(1:48,:,1);
            sortedCorrectApproach{i} = double(currentApproachVector(sortedIndex,2,1) == 0 & currentApproachVector(sortedIndex,1,1) == 1);
            sortedIncrectApproach{i} = double(currentApproachVector(sortedIndex,2,1) == 1);
            
    elseif startsWith(masterTbl.SessionType{i},['Reversal'])
        
            currentTrialTimes = trlStruct(i).trialTimes;
            [~,sortedIndex] = sort(currentTrialTimes(1:48));
            currentApproachVector = trlStruct(i).approach(1:48,:,1);
            sortedCorrectApproach{i} = double(currentApproachVector(sortedIndex,2,1) == 1 & currentApproachVector(sortedIndex,1,1) == 0);
            sortedIncrectApproach{i} = double(currentApproachVector(sortedIndex,1,1) == 1);

    end
end

% Index out junk from approaches 
% P rat
pConCorrect = sortedCorrectApproach(regPIdx);
pConIncorrect = sortedIncrectApproach(regPIdx);
pIncCorrect = sortedCorrectApproach(revPIdx);
pIncIncorrect = sortedIncrectApproach(revPIdx);
% Wistar
wConCorrect = sortedCorrectApproach(regWIdx);
wConIncorrect = sortedIncrectApproach(regWIdx);
wIncCorrect = sortedCorrectApproach(revWIdx);
wIncIncorrect = sortedIncrectApproach(revWIdx);
%% Pull Mean Data
% Wistars
clear wConCorrectPSTH wConIncorrectPSTH wIncCorrectPSTH wIncIncorrectPSTH pConCorrectPSTH pConIncorrectPSTH pIncCorrectPSTH pIncIncorrectPSTH

timepoint = 'all';
for k = 1:length(wConCorrect)
    wConCorrectPSTH(k,:,:) = squeeze(mean(wCon(logical(wConCorrect{k}(1:15)),:,:),1));
    wConIncorrectPSTH(k,:,:) = squeeze(mean(wCon(logical(wConIncorrect{k}(1:15)),:,:),1));

    wIncCorrectPSTH(k,:,:) = squeeze(mean(wInc(logical(wIncCorrect{k}(1:15)),:,:),1));
    wIncIncorrectPSTH(k,:,:) = squeeze(mean(wInc(logical(wIncIncorrect{k}(1:15)),:,:),1));
end
% P rats
for k = 1:length(pConCorrect)
    pConCorrectPSTH(k,:,:) = squeeze(mean(pCon(logical(pConCorrect{k}(1:15)),:,:),1));
    pConIncorrectPSTH(k,:,:) = squeeze(mean(pCon(logical(pConIncorrect{k}(1:15)),:,:),1));

    pIncCorrectPSTH(k,:,:) = squeeze(mean(pInc(logical(pIncCorrect{k}(1:15)),:,:),1));
    pIncIncorrectPSTH(k,:,:) = squeeze(mean(pInc(logical(pIncIncorrect{k}(1:15)),:,:),1));
end

%% Preprocess data for PCA -- Take Mean over Session # (First Dimension)
wConCorrectPSTH = squeeze(mean(wConCorrectPSTH,1));
wConIncorrectPSTH = squeeze(mean(wConIncorrectPSTH,1));
wIncCorrectPSTH = squeeze(mean(wIncCorrectPSTH,1));
wIncIncorrectPSTH = squeeze(mean(wIncIncorrectPSTH,1));

pConCorrectPSTH = squeeze(mean(pConCorrectPSTH,1));
pConIncorrectPSTH = squeeze(mean(pConIncorrectPSTH,1));
pIncCorrectPSTH = squeeze(mean(pIncCorrectPSTH,1));
pIncIncorrectPSTH = squeeze(mean(pIncIncorrectPSTH,1));

%% Implement PCA
[wCCCoeff,wCCScore,~,~,wCCExplained] = pca((wConCorrectPSTH));
[wCICoeff,wCIScore,~,~,wCIExplained] = pca((wConIncorrectPSTH));
[wICCoeff,wICScore,~,~,wICExplained] = pca((wIncCorrectPSTH));
[wIICoeff,wIIScore,~,~,wIIExplained] = pca((wIncIncorrectPSTH));

[pCCCoeff,pCCScore,~,~,pCCExplained] = pca((pConCorrectPSTH));
[pCICoeff,pCIScore,~,~,pCIExplained] = pca((pConIncorrectPSTH));
[pICCoeff,pICScore,~,~,pICExplained] = pca((pIncCorrectPSTH));
[pIICoeff,pIIScore,~,~,pIIExplained] = pca((pIncIncorrectPSTH));

% Centralize all output variables into a few cell arrays
allCoeff = [{wCCCoeff} {wCICoeff} {wICCoeff} {wIICoeff} {pCCCoeff} {pCICoeff} {pICCoeff} {pIICoeff}]; 
allScore = [{wCCScore} {wCIScore} {wICScore} {wIIScore} {pCCScore} {pCIScore} {pICScore} {pIIScore}];
allExplained = [{wCCExplained} {wCIExplained} {wICExplained} {wIIExplained} {pCCExplained} ...
    {pCIExplained} {pICExplained} {pIIExplained}];


% 04/25/2022 --- Based on a quick look of the explained variance, splitting
% the PCA based on behaviors did not pan out. It may be more beneficial to
% project the sort the coefficients for each group and work backwards
%% Plot
% 8 elements, 8 colors
colors = [{[0.2 0.24 0.82]} {[0.65 0.67 0.81]} {[0.66 0.16 0.16]} {[0.71 0.55 0.55]} {[0.96 0.9 0.02]} {[0.95 0.91 0.54]} {[0.51 0.01 0.94]} {[0.81 0.69 0.92]}];
titles = [{'Congruent Correct (W)'} {'Congruent Incorrect (W)'} {'Incongruent Correct (W)'} {'Incongruent Incorrect (W)'} ...
    {'Congruent Correct (P)'} {'Congruent Incorrect (P)'} {'Incongruent Correct (P)'} {'Incongruent Incorrect (P)'}];
%% Quick PCA Plot % 1st PC
pcSort = 1:3;
figure('Units','normalized','Position',[0 0 1 1])
for j = 1:length(pcSort)
    for i=1:length(allCoeff)
        
        subplot(max(pcSort),1,j);
        plot(time,allScore{i}(:,j),'LineWidth',3,'Color',colors{i})
        hold on
        xlabel('Time (s)');
        title(['PC Num: ' num2str(pcSort(j))]);
        xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',10,'FontWeight','bold','LineWidth',4)
         
    end
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

end
saveas(gca,[figPath filesep 'pcScore_numPCs_' num2str(pcSort) timepoint '_100ms'],'svg');     saveas(gca,[figPath filesep 'pcScore_numPCs_' num2str(pcSort) timepoint '_100ms'],'png')

%% Sort Raw Values According to Coefficients

for j = 1:length(allCoeff)
    figure('Units','normalized','Position',[0 0 1 1])
    for i=1:length(pcSort)
        
        [b,k]=sort(allCoeff{j}(:,pcSort(i)));

        subplot(1,max(pcSort),i);
        centData = allCoeff{j} * allScore{j}';
        imagesc([0 max(time)], [1 length(b)],centData(k,:),[-1 1]);

        xlabel('Time (s)');
        if i==1; ylabel('Neuron # (Sorted by PC)');end
        xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        title(['PC Num: ' num2str(pcSort(i)) ',ExplVar=' num2str(sum(allExplained{j}(1:i)),'%.1f') '%']);
        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
         
    end
    sgtitle(titles{j})
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figPath filesep 'pcProj_' titles{j} '_' timepoint '_100ms'],'svg');     saveas(gca,[figPath filesep 'pcProj_' titles{j} '_' timepoint '_100ms'],'png')

end
%% Histograms of Positive, Negative Loaders

%%
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
%% Replicating Analyses Above on a Shorter Timescale
timepoint = 'sipDescent';

pConPre = squeeze(mean(pCon(1:15,1:sipDescent*10,:),1));
pIncPre = squeeze(mean(pInc(1:15,1:sipDescent*10,:),1));
wConPre = squeeze(mean(wCon(1:15,1:sipDescent*10,:),1));
wIncPre = squeeze(mean(wInc(1:15,1:sipDescent*10,:),1));

[pConCoeff,pConScore,~,~,pConExplained] = pca(zscore(pConPre));
[pIncCoeff,pIncScore,~,~,pIncExplained] = pca(zscore(pIncPre));

[wConCoeff,wConScore,~,~,wConExplained] = pca(zscore(wConPre));
[wIncCoeff,wIncScore,~,~,wIncExplained] = pca(zscore(wIncPre));

% Centralize all output variables into a few cell arrays
allCoeff = [{pConCoeff} {pIncCoeff} {wConCoeff} {wIncCoeff}]; 
allScore = [{pConScore} {pIncScore} {wConScore} {wIncScore}];
allExplained = [{pConExplained} {pIncExplained} {wConExplained} {wIncExplained}];
titles = [{'Congruent P'} {'Incongruent P'} {'Congruent W'} {'Incongruent W'}];

%% Quick PCA Plot % 1st PC
pcSort = 1:3;
time = (1:size(pConPre,1))/10;
colors = [{'#EDB120'}, {'#7E2F8E'}, {'#0072BD'},  {'#D95319'} ];

figure('Units','normalized','Position',[0 0 1 1])

for j = 1:length(pcSort)
    for i=1:length(allCoeff)
        
        subplot(max(pcSort),1,j);
        plot(time,allScore{i}(:,j),'LineWidth',3,'Color',colors{i})
        hold on
        xlabel('Time (s)');
        title(['PC Num: ' num2str(pcSort(j))]);
        xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',10,'FontWeight','bold','LineWidth',4)
         
    end
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

end
saveas(gca,[figPath filesep 'pcScore_numPCs_' timepoint num2str(pcSort) '_100ms'],'svg');     saveas(gca,[figPath filesep 'pcScore_numPCs_' timepoint num2str(pcSort) '_100ms'],'png')

%% Sort Raw Values According to Coefficients

for j = 1:length(allCoeff)
    figure('Units','normalized','Position',[0 0 1 1])
    for i=1:length(pcSort)
        
        [b,k]=sort(allCoeff{j}(:,pcSort(i)));

        subplot(1,max(pcSort),i);
        centData = allCoeff{j} * allScore{j}';
        imagesc([0 max(time)], [1 length(b)],centData(k,:),[-1 1]);

        xlabel('Time (s)');
        if i==1; ylabel('Neuron # (Sorted by PC)');end
        xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        title(['PC Num: ' num2str(pcSort(i)) ',ExplVar=' num2str(sum(allExplained{j}(1:i)),'%.1f') '%']);
        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
         
    end
    sgtitle(titles{j})
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figPath filesep 'pcProj_' timepoint titles{j} '_100ms'],'svg');     saveas(gca,[figPath filesep 'pcProj_' timepoint titles{j} '_100ms'],'png')

end

%% PC Space Projections
% 2D
figure('Units','normalized','Position',[0 0 1 1])
for j = 1:length(allCoeff)
    plot(allScore{j}(:,1),allScore{j}(:,2),'Color',colors{j},'LineWidth',3);
    hold on
    plot(allScore{j}(1,1),allScore{j}(1,2),'o','MarkerSize',10,'MarkerFaceColor',colors{j},'MarkerEdgeColor','k')
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)        
    

end
    saveas(gca,[figPath filesep 'pcSpaceProj_' timepoint '_100ms'],'svg');     saveas(gca,[figPath filesep 'pcSpaceProj_' timepoint '_100ms'],'png')
%% 3D
figure('Units','normalized','Position',[0 0 1 1])
for j = 1:length(allCoeff)
    plot3(allScore{j}(:,1),allScore{j}(:,2),allScore{j}(:,3),'Color',colors{j},'LineWidth',3);
    hold on
    plot3(allScore{j}(1,1),allScore{j}(1,2),allScore{j}(1,3),'o','MarkerSize',10,'MarkerFaceColor',colors{j},'MarkerEdgeColor','k')
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)        
    
end
    saveas(gca,[figPath filesep '3DpcSpaceProj_' timepoint '_100ms'],'svg');     saveas(gca,[figPath filesep '3DpcSpaceProj_' timepoint '_100ms'],'png')

%% Calculate and Plot the Distance Between PCs 
pCon3D = allScore{1}(:,1:3); wCon3D = allScore{3}(:,1:3);
pInc3D = allScore{2}(:,1:3); wInc3D = allScore{4}(:,1:3);

pDistance = diag(pdist2(pCon3D,pInc3D,'mahalanobis'));
wDistance = diag(pdist2(wCon3D,wInc3D,'mahalanobis'));

figure('Units','normalized','Position',[0 0 1 1])
plot(time,pDistance,'Color',colors{1},'LineWidth',3)
hold on
plot(time,wDistance,'Color',colors{3},'LineWidth',3)
xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
xlabel('Time (s)')
ylabel('Mahalanobis Distance')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,[figPath filesep '3DpcSpaceDistance_' timepoint '_100ms'],'svg');     saveas(gca,[figPath filesep '3DpcSpaceDistance_' timepoint '_100ms'],'png')

%% Determine a way to find how predictive each PC is of Correct and Incorrect Choices
% 05/04/2022
% Utilizing PCs from prior to the sipper descent, can we determine how
% predictive each PC is of later decisions? This activity then is likely
% critical for disambiguating correct/incorrect encoding & can be useful
% for subsequently analyzing congruent/incongruent choices. 

% One possible way to start --- find positive and negative loading
% coefficients for each PC and see how the mean of those neurons change
% over trials. We can also then determine how effective those neurons are
% at predicting later outcomes via the 'raw' data. 

% Find index of positive / negative loading coefficients, use that index on pCon,
% pInc, wCon, wInc variables. Take mean of those neurons to generate a heat
% map of mean activity over time and trials. 
time = (1:size(pCon,2))/10;

for j = 1:length(allCoeff)
    pcSort = 1:3;
    % Assign variable rawData depending on condition 
    if j == 1
        rawData = pCon;
    elseif j == 2
        rawData = pInc;
    elseif j == 3
        rawData = wCon;
    elseif j == 4
        rawData = wInc;
    end
    % Sort according to top 3 PCs, plot similar to other coefficient
    % figures

    figure('Units','normalized','Position',[0 0 1 1])

    for i = 1:length(pcSort)
        % Index Pos and Neg
        PosCoeff = rawData(:,:,allCoeff{j}(:,i) > 0.01);                  
        NegCoeff = rawData(:,:,allCoeff{j}(:,i) < -0.01);
        % Take Mean
        PosCoeff = zscore(mean(PosCoeff,3));
        NegCoeff = zscore(mean(NegCoeff,3));
        % Figure
        subplot(2,max(pcSort),i);          % Doubling subplot for pos+neg
        % Positive Plot
        imagesc([0 max(time)], [1 15],PosCoeff, [-3 3]);
        colormap jet
        colorbar
        xlabel('Time (s)');
        if i==1; ylabel('Trial #');end
        xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        
        title(['PC Num: ' num2str(pcSort(i)) ',ExplVar=' num2str(sum(allExplained{j}(1:i)),'%.1f') '%']);
        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
        % Negative Plot

        subplot(2,max(pcSort),i+3);          % Doubling subplot for pos+neg

        imagesc([0 max(time)], [1 15],NegCoeff,[-3 3]);
        colormap jet
        colorbar
        xlabel('Time (s)');
        if i==1; ylabel('Trial #');end
        xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        
        title(['PC Num: ' num2str(pcSort(i)) ',ExplVar=' num2str(sum(allExplained{j}(1:i)),'%.1f') '%']);
        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)        
    end

    sgtitle(titles{j})
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figPath filesep 'coeffNeurons_overTrials' titles{j} '_100ms'],'svg');     saveas(gca,[figPath filesep 'coeffNeurons_overTrials' titles{j} '_100ms'],'png')

end

% Of the 15 first trials that we have access to, can we see if there are
% meaningful differences in neural activity between correct and incorrect
% choices? 

% A strategy, find the index of correct and incorrect choices. Of those
% choices, take the mean of 


%% Focus on Congruent Sessions Only & The Differences Between Correct & 
% Incorrect Choices as Well as Corrections. We are essentially repeating
% the analysis from above except only on Congruent sessions. This is to
% determine differences between correct and incorrect choices at baseline.
% 04/27/2022
%% Split electrophysiology data based on behaviors
for i = 1:length(trlStruct)
    if startsWith(masterTbl.SessionType{i},['Regular'])

            currentTrialTimes = trlStruct(i).trialTimes;
            [~,sortedIndex] = sort(currentTrialTimes(1:48));
            currentApproachVector = trlStruct(i).approach(1:48,:,1);
            sortedCorrectApproach{i} = double(currentApproachVector(sortedIndex,2,1) == 0 & currentApproachVector(sortedIndex,1,1) == 1);
            sortedIncrectApproach{i} = double(currentApproachVector(sortedIndex,2,1) == 1 & currentApproachVector(sortedIndex,1,1) == 0);
            sortedCorrections{i} = double(currentApproachVector(sortedIndex,2,1) == 1 & currentApproachVector(sortedIndex,1,1) == 1);

    end
end

% Index out junk from approaches 
% P rat
pConCorrect = sortedCorrectApproach(regPIdx);
pConIncorrect = sortedIncrectApproach(regPIdx);
pConCorrections = sortedCorrections(regPIdx);
% Wistar
wConCorrect = sortedCorrectApproach(regWIdx);
wConIncorrect = sortedIncrectApproach(regWIdx);
wConCorrections = sortedCorrections(regWIdx);
% 
%% Pull Mean Data
% Wistars
clear wConCorrectPSTH wConIncorrectPSTH wIncCorrectPSTH wIncIncorrectPSTH pConCorrectPSTH pConIncorrectPSTH pIncCorrectPSTH pIncIncorrectPSTH wConCorrectionsPSTH pConCorrectionsPSTH

timepoint = 'all';
for k = 1:length(wConCorrect)
    wConCorrectPSTH(k,:,:) = squeeze(nanmean(wCon(logical(wConCorrect{k}(1:15)),:,:),1));
    wConIncorrectPSTH(k,:,:) = squeeze(nanmean(wCon(logical(wConIncorrect{k}(1:15)),:,:),1));
    wConCorrectionsPSTH(k,:,:) = squeeze(nanmean(wCon(logical(wConCorrections{k}(1:15)),:,:),1));
end
% P rats
for k = 1:length(pConCorrect)
    pConCorrectPSTH(k,:,:) = squeeze(nanmean(pCon(logical(pConCorrect{k}(1:15)),:,:),1));
    pConIncorrectPSTH(k,:,:) = squeeze(nanmean(pCon(logical(pConIncorrect{k}(1:15)),:,:),1));
    pConCorrectionsPSTH(k,:,:) = squeeze(nanmean(wCon(logical(pConCorrections{k}(1:15)),:,:),1));


end

%% Preprocess data for PCA -- Take Mean over Session # (First Dimension)
wConCorrectPSTH = squeeze(nanmean(wConCorrectPSTH,1));
wConIncorrectPSTH = squeeze(nanmean(wConIncorrectPSTH,1));
wConCorrectionsPSTH = squeeze(nanmean(wConCorrectionsPSTH,1));

pConCorrectPSTH = squeeze(nanmean(pConCorrectPSTH,1));
pConIncorrectPSTH = squeeze(nanmean(pConIncorrectPSTH,1));
pConCorrectionsPSTH = squeeze(nanmean(pConCorrectionsPSTH,1));

%% Implement PCA
% Code:
% CC = Congruent, Correct   IC = Incongruent, Correct
% CI = Congruent, Incorrect II = Incongruent, Incorrect
% CM = Congruent, Mulligan (Corrections) IM = Incongruent, Mulligan
% (Corrections)

[wCCCoeff,wCCScore,~,~,wCCExplained] = pca((wConCorrectPSTH));
[wCICoeff,wCIScore,~,~,wCIExplained] = pca((wConIncorrectPSTH));
[wCMCoeff,wCMScore,~,~,wCMExplained] = pca((wConCorrectionsPSTH));


[pCCCoeff,pCCScore,~,~,pCCExplained] = pca((pConCorrectPSTH));
[pCICoeff,pCIScore,~,~,pCIExplained] = pca((pConIncorrectPSTH));
[pCMCoeff,pCMScore,~,~,pCMExplained] = pca((pConCorrectionsPSTH));


% Centralize all output variables into a few cell arrays
allCoeff = [{wCCCoeff} {wCICoeff} {wCMCoeff} {pCCCoeff} {pCICoeff} {pCMCoeff}]; 
allScore = [{wCCScore} {wCIScore} {wCMScore} {pCCScore} {pCIScore} {pCMScore}];
allExplained = [{wCCExplained} {wCIExplained} {wCMExplained} {pCCExplained} ...
    {pCIExplained} {pCMExplained}];

%% Plot
% 8 elements, 8 colors
colors = [{[0.2 0.24 0.82]} {[0.2 0.24 0.82]} {[0.2 0.24 0.82]} {[0.96 0.9 0.02]} {[0.96 0.9 0.02]} {[0.96 0.9 0.02]}];
lineStyle = [{'-'} {':'} {'-.'} {'-'} {':'} {'-.'}];
titles = [{'Congruent Correct (W)'} {'Congruent Incorrect (W)'} {'Congruent Corrections (W)'} ...
    {'Congruent Correct (P)'} {'Congruent Incorrect (P)'} {'Congruent Corrections (P)'}];
%% Quick PCA Plot % 1st PC
pcSort = 1:3;
time = (1:size(pCon,2))/10;

figure('Units','normalized','Position',[0 0 1 1])

for j = 1:length(pcSort)
    for i=1:length(allCoeff)
        
        subplot(length(pcSort),1,j);
        plot(time,allScore{i}(:,pcSort(j)),'LineWidth',3,'Color',colors{i},'LineStyle',lineStyle{i})
        hold on
        xlabel('Time (s)');
        title(['PC Num: ' num2str(pcSort(j))]);
        xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',10,'FontWeight','bold','LineWidth',4)
         
    end
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

end
saveas(gca,[figPath filesep 'CongruentOnly_pcScore_numPCs_' num2str((pcSort)) timepoint '_100ms'],'svg');     saveas(gca,[figPath filesep 'CongruentOnly_pcScore_numPCs_' num2str((pcSort)) '_100ms'],'png')

%% Sort Raw Values According to Coefficients

for j = 1:length(allCoeff)
    figure('Units','normalized','Position',[0 0 1 1])
    for i=1:length(pcSort)
        
        [b,k]=sort(allCoeff{j}(:,pcSort(i)));

        subplot(1,max(pcSort),i);
        centData = allCoeff{j} * allScore{j}';
        imagesc([0 max(time)], [1 length(b)],centData(k,:),[-1 1]);

        xlabel('Time (s)');
        if i==1; ylabel('Neuron # (Sorted by PC)');end
        xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        title(['PC Num: ' num2str(pcSort(i)) ',ExplVar=' num2str(sum(allExplained{j}(1:i)),'%.1f') '%']);
        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
         
    end
    sgtitle(titles{j})
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figPath filesep 'CongruentOnly_pcProj_' titles{j} '_' timepoint '_100ms'],'svg');     saveas(gca,[figPath filesep 'CongruentOnly_pcProj_' titles{j} '_' timepoint '_100ms'],'png')

end


%% Repeat Above Analyses on Smaller Timescales 
% PreCue, CueOn, Sipper In
%% Pull Mean Data
% Wistars
clear wConCorrectPSTH wConIncorrectPSTH wIncCorrectPSTH wIncIncorrectPSTH pConCorrectPSTH pConIncorrectPSTH pIncCorrectPSTH pIncIncorrectPSTH wConCorrectionsPSTH pConCorrectionsPSTH

timepoint = 'PreCue_CueOn_SipperIn';
for k = 1:length(wConCorrect)
    wConCorrectPSTH(k,:,:) = squeeze(nanmean(wCon(logical(wConCorrect{k}(1:15)),1:sipDescent*10,:),1));
    wConIncorrectPSTH(k,:,:) = squeeze(nanmean(wCon(logical(wConIncorrect{k}(1:15)),1:sipDescent*10,:),1));
    wConCorrectionsPSTH(k,:,:) = squeeze(nanmean(wCon(logical(wConCorrections{k}(1:15)),1:sipDescent*10,:),1));
end
% P rats
for k = 1:length(pConCorrect)
    pConCorrectPSTH(k,:,:) = squeeze(nanmean(pCon(logical(pConCorrect{k}(1:15)),1:sipDescent*10,:),1));
    pConIncorrectPSTH(k,:,:) = squeeze(nanmean(pCon(logical(pConIncorrect{k}(1:15)),1:sipDescent*10,:),1));
    pConCorrectionsPSTH(k,:,:) = squeeze(nanmean(wCon(logical(pConCorrections{k}(1:15)),1:sipDescent*10,:),1));


end

%% Preprocess data for PCA -- Take Mean over Session # (First Dimension)
wConCorrectPSTH = squeeze(nanmean(wConCorrectPSTH,1));
wConIncorrectPSTH = squeeze(nanmean(wConIncorrectPSTH,1));
wConCorrectionsPSTH = squeeze(nanmean(wConCorrectionsPSTH,1));

pConCorrectPSTH = squeeze(nanmean(pConCorrectPSTH,1));
pConIncorrectPSTH = squeeze(nanmean(pConIncorrectPSTH,1));
pConCorrectionsPSTH = squeeze(nanmean(pConCorrectionsPSTH,1));

%% Implement PCA
% Code:
% CC = Congruent, Correct   IC = Incongruent, Correct
% CI = Congruent, Incorrect II = Incongruent, Incorrect
% CM = Congruent, Mulligan (Corrections) IM = Incongruent, Mulligan
% (Corrections)

[wCCCoeff,wCCScore,~,~,wCCExplained] = pca((wConCorrectPSTH));
[wCICoeff,wCIScore,~,~,wCIExplained] = pca((wConIncorrectPSTH));
[wCMCoeff,wCMScore,~,~,wCMExplained] = pca((wConCorrectionsPSTH));


[pCCCoeff,pCCScore,~,~,pCCExplained] = pca((pConCorrectPSTH));
[pCICoeff,pCIScore,~,~,pCIExplained] = pca((pConIncorrectPSTH));
[pCMCoeff,pCMScore,~,~,pCMExplained] = pca((pConCorrectionsPSTH));


% Centralize all output variables into a few cell arrays
allCoeff = [{wCCCoeff} {wCICoeff} {wCMCoeff} {pCCCoeff} {pCICoeff} {pCMCoeff}]; 
allScore = [{wCCScore} {wCIScore} {wCMScore} {pCCScore} {pCIScore} {pCMScore}];
allExplained = [{wCCExplained} {wCIExplained} {wCMExplained} {pCCExplained} ...
    {pCIExplained} {pCMExplained}];

%% Plot
% 8 elements, 8 colors
colors = [{[0.2 0.24 0.82]} {[0.2 0.24 0.82]} {[0.2 0.24 0.82]} {[0.96 0.9 0.02]} {[0.96 0.9 0.02]} {[0.96 0.9 0.02]}];
lineStyle = [{'-'} {':'} {'-.'} {'-'} {':'} {'-.'}];
titles = [{'Congruent Correct (W)'} {'Congruent Incorrect (W)'} {'Congruent Corrections (W)'} ...
    {'Congruent Correct (P)'} {'Congruent Incorrect (P)'} {'Congruent Corrections (P)'}];
%% Quick PCA Plot % 1st PC
pcSort = 1:3;
time = (1:size(pConCorrectPSTH,1))/10;

figure('Units','normalized','Position',[0 0 1 1])
for j = 1:length(pcSort)
    for i=1:length(allCoeff)
        
        subplot(length(pcSort),1,j);
        plot(time,allScore{i}(:,pcSort(j)),'LineWidth',3,'Color',colors{i},'LineStyle',lineStyle{i})
        hold on
        xlabel('Time (s)');
        title(['PC Num: ' num2str(pcSort(j))]);
        xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',10,'FontWeight','bold','LineWidth',4)
         
    end
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

end
saveas(gca,[figPath filesep 'CongruentOnly_pcScore_numPCs_' num2str((pcSort)) timepoint '_100ms'],'svg');     saveas(gca,[figPath filesep 'CongruentOnly_pcScore_numPCs_' num2str((pcSort)) timepoint '_100ms'],'png')

%% Sort Raw Values According to Coefficients

for j = 1:length(allCoeff)
    figure('Units','normalized','Position',[0 0 1 1])
    for i=1:length(pcSort)
        
        [b,k]=sort(allCoeff{j}(:,pcSort(i)));

        subplot(1,max(pcSort),i);
        centData = allCoeff{j} * allScore{j}';
        imagesc([0 max(time)], [1 length(b)],centData(k,:),[-1 1]);

        xlabel('Time (s)');
        if i==1; ylabel('Neuron # (Sorted by PC)');end
        xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        title(['PC Num: ' num2str(pcSort(i)) ',ExplVar=' num2str(sum(allExplained{j}(1:i)),'%.1f') '%']);
        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
         
    end
    sgtitle(titles{j})
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figPath filesep 'CongruentOnly_pcProj_' titles{j} '_' timepoint '_100ms'],'svg');     saveas(gca,[figPath filesep 'CongruentOnly_pcProj_' titles{j} '_' timepoint '_100ms'],'png')

end

%%
%% PC Space Projections
% 05/03/2022
% Find a systematic way of determining which PCs can be predictive of
% outcomes. 
%

figure('Units','normalized','Position',[0 0 1 1])
colors = [{[0.2 0.24 0.82]} {[0.65 0.67 0.81]} {[0.16 0.19 0.39]} {[0.96 0.9 0.02]} {[0.95 0.91 0.54]} {[0.7 0.65 0.33]}];
% allScoreMatrix = cat(3,allScore{:});
% allScoreMatrix = squeeze(allScoreMatrix(60,:,:));
% allScoreMatrix = reshape(allScoreMatrix,36,11);
% Y = tsne(allScoreMatrix);
% labels = [repmat({'WCC'},6,1); repmat({'WCI'},6,1); repmat({'WCCor'},6,1); repmat({'PCC'},6,1); repmat({'PCI'},6,1); repmat({'PCCor'},6,1)];
% gscatter(Y(:,1),Y(:,2))

for j = 1:length(allCoeff)
    scatter(allScore{j}(40:60,2),allScore{j}(40:60,8),'MarkerEdgeColor',colors{j},'LineWidth',3);
    hold on
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)        
end
    saveas(gca,[figPath filesep 'pcSeparations_' timepoint '_100ms'],'svg');     saveas(gca,[figPath filesep 'pcSeparations_' timepoint '_100ms'],'png')


%% Plot variance explained 

for j = 1:length(allCoeff)
    figure('Units','normalized','Position',[0 0 1 1])
        
    plot(allExplained{j},'ko');

    xlabel('PC');
    ylabel('% Explained')
    title(titles{j});         
   
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figPath filesep 'explainedVariance' titles{j} '_' timepoint '_100ms'],'svg');     saveas(gca,[figPath filesep 'explainedVariance' titles{j} '_' timepoint '_100ms'],'png')

end

%% Analysis Goals 05/09/2022
% Eliminate Low FR Neurons from Analysis 
% Create 2x2 Genotype x Session Type Figure with Top 5 PCs in each
% 
timepoint = 'sipDescent';

if strcmp(timepoint,'sipDescent')
pConPre = squeeze(mean(pCon(1:15,1:sipDescent*10,:),1));
pIncPre = squeeze(mean(pInc(1:15,1:sipDescent*10,:),1));
wConPre = squeeze(mean(wCon(1:15,1:sipDescent*10,:),1));
wIncPre = squeeze(mean(wInc(1:15,1:sipDescent*10,:),1));
time = (1:size(pConPre,1))/10;

[pConCoeff,pConScore,~,~,pConExplained] = pca(zscore(pConPre));
[pIncCoeff,pIncScore,~,~,pIncExplained] = pca(zscore(pIncPre));

[wConCoeff,wConScore,~,~,wConExplained] = pca(zscore(wConPre));
[wIncCoeff,wIncScore,~,~,wIncExplained] = pca(zscore(wIncPre));

% Centralize all output variables into a few cell arrays
allCoeff = [{pConCoeff} {pIncCoeff} {wConCoeff} {wIncCoeff}]; 
allScore = [{pConScore} {pIncScore} {wConScore} {wIncScore}];
allExplained = [{pConExplained} {pIncExplained} {wConExplained} {wIncExplained}];
titles = [{'Congruent P'} {'Incongruent P'} {'Congruent W'} {'Incongruent W'}];
%
elseif strcmp(timepoint,'all')
pConPre = squeeze(mean(pCon(1:15,:,:),1));
pIncPre = squeeze(mean(pInc(1:15,:,:),1));
wConPre = squeeze(mean(wCon(1:15,:,:),1));
wIncPre = squeeze(mean(wInc(1:15,:,:),1));
time = (1:size(pConPre,1))/10;

[pConCoeff,pConScore,~,~,pConExplained] = pca(zscore(pConPre));
[pIncCoeff,pIncScore,~,~,pIncExplained] = pca(zscore(pIncPre));

[wConCoeff,wConScore,~,~,wConExplained] = pca(zscore(wConPre));
[wIncCoeff,wIncScore,~,~,wIncExplained] = pca(zscore(wIncPre));

% Centralize all output variables into a few cell arrays
allCoeff = [{pConCoeff} {pIncCoeff} {wConCoeff} {wIncCoeff}]; 
allScore = [{pConScore} {pIncScore} {wConScore} {wIncScore}];
allExplained = [{pConExplained} {pIncExplained} {wConExplained} {wIncExplained}];
titles = [{'Congruent P'} {'Incongruent P'} {'Congruent W'} {'Incongruent W'}];
end
%% PCA Plot
% First 5 PCs per Each Condition

pcSort = 1:5;
colors = [{'#0072BD'} {'#4DBEEE'} {'#77AC30'} {'#D95319'} {'#A2142F'}];

figure('Units','normalized','Position',[0 0 1 1])

for j = 1:length(allCoeff)
    for i=1:length(pcSort)
        
        subplot(2,2,j);
        plot(time,allScore{j}(:,i),'LineWidth',3,'Color',colors{i})
        if j == 1; legend([{'PC1'},{'PC2'},{'PC3'},{'PC4'},{'PC5'}],'AutoUpdate','off','location','best'); end
        hold on
        xlabel('Time (s)');
        title([titles{j}]);
        ylim([-25 25])
        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',10,'FontWeight','bold','LineWidth',4)
        xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        if strcmp(timepoint,'all')
        xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        end
    end
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

end
saveas(gca,[figPath filesep 'top5pcs' timepoint num2str(pcSort) '_100ms'],'svg');     saveas(gca,[figPath filesep 'top5pcs' timepoint num2str(pcSort) '_100ms'],'png')


%% PCA: Force Neurons from Separate Conditions into the Same PCs
% 05/10/2022

%% PCA Preprocessing | No separation between P, W, Con, Inc groups
timepoint = 'all';

allConditions = cat(3,ephysStruct(:).PSTH);
if strcmp(timepoint,'all')
    allConditionsPre = squeeze(mean(allConditions(1:15,:,:),1));
    time = (1:size(allConditionsPre,1))/10;

elseif strcmp(timepoint,'sipDescent')
    allConditionsPre = squeeze(mean(allConditions(1:15,1:sipDescent*10,:),1));
    time = (1:size(allConditionsPre,1))/10;
end
[allConditionsCoeff,allConditionsScore,~,~,allConditionsExplained] = pca(allConditionsPre);

% Need to find an index of neurons per dataset, be able to sort that index
% into our conditions (W/P; Con/Inc). Should ideally be a 1x54 matrix
% containing N # of neurons.

neuronIndex = [ephysStruct(:).neuronNum];

% CLear structures
% ephysStruct = [];
% trlStruct = [];
% masterTbl = [];

%% Sanity check
if sum(neuronIndex) == size(allConditionsCoeff,1)
    disp('Matrix sizes match')
else
    disp('Matrix size mismatch, address')
end

%% Assign ID to Neurons
neuronID = [];
for i = 1:length(neuronIndex)
    neuronID = [neuronID ones(1,neuronIndex(i))*i];
end
%% Sanity check
if length(neuronID) == size(allConditionsCoeff,1)
    disp('Matrix sizes match')
else
    disp('Matrix size mismatch, address')
end

%% PC Scores Without Splits
pcSort = 1:5;
colors = [{'#0072BD'} {'#4DBEEE'} {'#77AC30'} {'#D95319'} {'#A2142F'}];
figure('Units','normalized','Position',[0 0 1 1])


 
for i=1:length(pcSort)
    plot(time,allConditionsScore(:,i),'LineWidth',3,'Color',colors{i})
    legend([{'PC1'},{'PC2'},{'PC3'},{'PC4'},{'PC5'}],'AutoUpdate','off','location','best');
    hold on
    xlabel('Time (s)');
    title('Top 5 PCs, All Neurons');
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',10,'FontWeight','bold','LineWidth',4)
    xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
    xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
    if strcmp(timepoint,'all')
    xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
    end

end
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)


saveas(gca,[figPath filesep 'allNeurons_top5pcs' timepoint num2str(pcSort) '_100ms'],'svg');     saveas(gca,[figPath filesep 'allNeurons_top5pcs' timepoint num2str(pcSort) '_100ms'],'png')

%% PC Coefficients Projections Without Splits
    figure('Units','normalized','Position',[0 0 1 1])
    for i=1:length(pcSort)
        
        [b,k]=sort(allConditionsCoeff(:,pcSort(i)));

        subplot(1,max(pcSort),i);
        centData = allConditionsCoeff * allConditionsScore';
        imagesc([0 max(time)], [1 length(b)],centData(k,:),[-1 1]);

        xlabel('Time (s)');
        if i==1; ylabel('Neuron # (Sorted by PC)');end
        xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        title(['PC Num: ' num2str(pcSort(i)) ',ExplVar=' num2str(sum(allConditionsExplained(1:i)),'%.1f') '%']);
        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
        if strcmp(timepoint,'all')
        xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        end
         
    end
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figPath filesep 'allNeurons_pcProj_' timepoint '_100ms'],'svg');     saveas(gca,[figPath filesep 'allNeurons_pcProj_' timepoint '_100ms'],'png')

%% Using neuronID variable, split data into 2x2 W/P Con/Inc
% Find session IDs associated with each 2x2 combination
conPIdx = find(regPIdx == 1);
incPIdx = find(revPIdx == 1); 
conWIdx = find(regWIdx == 1);
incWIdx = find(revWIdx == 1);

% Find raw, centered data
centData = allConditionsCoeff * allConditionsScore';

% Index 
conPCent = centData(any(neuronID == conPIdx),:);
incPCent = centData(any(neuronID == incPIdx),:);
conWCent = centData(any(neuronID == conWIdx),:);
incWCent = centData(any(neuronID == incWIdx),:);
% Coefficient
conPCoeff = allConditionsCoeff(any(neuronID == conPIdx),:);
incPCoeff = allConditionsCoeff(any(neuronID == incPIdx),:);
conWCoeff = allConditionsCoeff(any(neuronID == conWIdx),:);
incWCoeff = allConditionsCoeff(any(neuronID == incWIdx),:);
% Create singular variable
indexedCentData = {conPCent incPCent conWCent incWCent};
indexedCoeffData = {conPCoeff incPCoeff conWCoeff incWCoeff};
cellIndex = {conPIdx incPIdx conWIdx incWIdx};
%% PC Coefficients Projections With Splits
titles = [{'Congruent P'} {'Incongruent P'} {'Congruent W'} {'Incongruent W'}];
pcSort = 1:5;

for j = 1:length(indexedCentData)
    figure('Units','normalized','Position',[0 0 1 1])
    for i=1:length(pcSort)
        
        [b,k]=sort(indexedCoeffData{j}(:,pcSort(i)));      
        subplot(1,length(pcSort),i);
        centData = indexedCoeffData{j} * allConditionsScore';
        imagesc([0 max(time)], [1 length(b)],centData(k,:),[-1 1]);

        xlabel('Time (s)');
        if i==1; ylabel('Neuron # (Sorted by PC)');end
        xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        title(['PC Num: ' num2str(pcSort(i)) ',ExplVar=' num2str(sum(allConditionsExplained(1:pcSort(i))),'%.1f') '%']);
        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
        if strcmp(timepoint,'all')
        xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        end
         
    end
    sgtitle(titles{j})
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figPath filesep 'allNeurons_pcProj_' titles{j} timepoint num2str(pcSort) '_100ms'],'svg');     saveas(gca,[figPath filesep 'allNeurons_pcProj_' titles{j} timepoint num2str(pcSort) '_100ms'],'png')
end

%% Low-D PC Representation Based on Separated Matrices
% Parameters
allPCProjections = [];
pcSort = 1:5;
% Obtain data
for i = 1:length(cellIndex)
    for k = 1:length(pcSort)
        allPCProjections{i}(:,k) = allConditionsPre(:,any(neuronID == cellIndex{i})) * indexedCoeffData{i}(:,pcSort(k));
    end
end

% allPCProjections ends up as a 1x4 Cell containing a (time x pcSort)
% matrix 

%% Plot Top 5 Resulting Matrices for Each Condition
pcSort = 1:5;
colors = [{'#0072BD'} {'#4DBEEE'} {'#77AC30'} {'#D95319'} {'#A2142F'}];

figure('Units','normalized','Position',[0 0 1 1])

for j = 1:length(allPCProjections)
    for i=1:length(pcSort)
        
        subplot(2,2,j);
        plot(time,allPCProjections{j}(:,i),'LineWidth',3,'Color',colors{i})
        if j == 1; legend([{'PC1'},{'PC2'},{'PC3'},{'PC4'},{'PC5'}],'AutoUpdate','off','location','best'); end
        hold on
        xlabel('Time (s)');
        title([titles{j}]);
        ylim([-25 50])
        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',10,'FontWeight','bold','LineWidth',4)
        xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        if strcmp(timepoint,'all')
        xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        end
    end
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
end
saveas(gca,[figPath filesep 'allNeurons_split_top5pcs' timepoint num2str(pcSort) '_100ms'],'svg');     saveas(gca,[figPath filesep 'allNeurons_split_top5pcs' timepoint num2str(pcSort) '_100ms'],'png')
%% Plot resulting data in a 2D Space
dim = 2;
figure('Units','normalized','Position',[0 0 1 1])
colors = [{'#EDB120'}, {'#7E2F8E'}, {'#0072BD'},  {'#D95319'} ];

for i = 1:length(allPCProjections)
    plot3(allPCProjections{i}(:,1),allPCProjections{i}(:,2),allPCProjections{i}(:,3),'Color',colors{i},'LineWidth',3);
    hold on
    plot3(allPCProjections{i}(1,1),allPCProjections{i}(1,2),allPCProjections{i}(1,3),'o','MarkerSize',10,'MarkerFaceColor',colors{i},'MarkerEdgeColor','k')
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)            
end

saveas(gca,[figPath filesep 'allNeurons_split_pcSpaceProj_' timepoint '_100ms'],'svg');     saveas(gca,[figPath filesep 'allNeurons_split_pcSpaceProj_' timepoint '_100ms'],'png')

%% Calculate and Plot the Distance Between PCs 
pCon3D = allPCProjections{1}(:,1:5); wCon3D = allPCProjections{3}(:,1:5);
pInc3D = allPCProjections{2}(:,1:5); wInc3D = allPCProjections{4}(:,1:5);

pDistance = diag(pdist2(pCon3D,pInc3D,'mahalanobis'));
wDistance = diag(pdist2(wCon3D,wInc3D,'mahalanobis'));

figure('Units','normalized','Position',[0 0 1 1])
plot(time,pDistance,'Color',colors{1},'LineWidth',3)
hold on
plot(time,wDistance,'Color',colors{3},'LineWidth',3)
xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
if strcmp(timepoint,'all')
xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
end
xlabel('Time (s)')
ylabel('Mahalanobis Distance')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,[figPath filesep 'allNeurons_split_3DpcSpaceDistance_' timepoint '_100ms'],'svg');     saveas(gca,[figPath filesep 'allNeurons_split_3DpcSpaceDistance_' timepoint '_100ms'],'png')

%% Determine how PC Score Changes Over Trials
%% Trial-by-trial eventually 
% Project coefficients onto trial-by-trial data until
% % Parameters
allPCProjections = [];
pcSort = 1:5;
permAllConditions = permute(allConditions,[3 2 1]);
% Obtain data
for i = 1:length(cellIndex)
    for k = 1:length(pcSort)
        TRLallPCProjections{i}(:,:,:,k) = (permAllConditions(any(neuronID == cellIndex{i}),:,:) .* indexedCoeffData{i}(:,pcSort(k)));
    end
end

% for i = 1:length(allPCProjections)
%     allPCProjections{i} = squeeze(mean(squeeze(mean(allPCProjections{i},3))));
% end
% allPCProjections ends up as a 1x4 Cell containing a (time x pcSort)
% matrix 

%% Plot (Start with One Thing)

for j = 1:length(TRLallPCProjections)
    pcSort = 1:3;
    % Assign variable rawData depending on condition 
    % Sort according to top 3 PCs, plot similar to other coefficient
    % figures

    figure('Units','normalized','Position',[0 0 1 1])

    for i = 1:length(pcSort)
        % Index Pos and Neg
        PosCoeff = squeeze(mean(TRLallPCProjections{j}(:,:,:,i),1))';                  
        % Take Mean
        % Figure
        subplot(1,max(pcSort),i);          % Doubling subplot for pos+neg
        % Positive Plot
        imagesc([0:0.1:max(time)], [1 15],PosCoeff, [-0.03 0.03]);
        colormap jet
        colorbar
        xlabel('Time (s)');
        if i==1; ylabel('Trial #');end
        xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        
        title(['PC Num: ' num2str(pcSort(i)) ',ExplVar=' num2str(sum(allConditionsExplained(1:i)),'%.1f') '%']);
        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)    
        allPosCoeff{j}{i} = PosCoeff;
    end

    sgtitle(titles{j})
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figPath filesep 'allNeurons_split_scores_overTrials' titles{j} '_100ms'],'svg');     saveas(gca,[figPath filesep 'allNeurons_split_scores_overTrials' titles{j} '_100ms'],'png')

end

%% 05/23/2022
% Final stretch of analyses before shifting towards RSA Poster
% Things to Do:
% Plot loadings according to condition 
% Determine what each PC is encoding for
% Based on these, start finding ways to generate statistical samples that
% can be compared.
%
colors = [{'#EDB120'}, {'#7E2F8E'}, {'#0072BD'},  {'#D95319'} ];
pcSort = 1:5;
% Plot coefficient loadings
figure('Units','normalized','Position',[0 0 1 1])

    for i = 1:length(pcSort)
        subplot(1,max(pcSort),i);
        histogram(indexedCoeffData{1}(:,i),'Normalization','probability','DisplayStyle','stairs','LineWidth',2,'LineStyle','-','EdgeColor',colors{1})
        hold on
        histogram(indexedCoeffData{2}(:,i),'Normalization','probability','DisplayStyle','stairs','LineWidth',2,'LineStyle','-','EdgeColor',colors{2})
        histogram(indexedCoeffData{3}(:,i),'Normalization','probability','DisplayStyle','stairs','LineWidth',2,'LineStyle','-','EdgeColor',colors{3})
        histogram(indexedCoeffData{4}(:,i),'Normalization','probability','DisplayStyle','stairs','LineWidth',2,'LineStyle','-','EdgeColor',colors{4})

        xlabel(['PC ' num2str(i) ' Coefficients'])
        ylabel('Probability')
        xlim([-.2,.2])
        ylim([0 0.45])

        legend([{'PCon'},{'PInc'},{'WCon'},{'WInc'}])
        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)


      [h1(i),p1(i)] = kstest2(indexedCoeffData{1}(:,i),indexedCoeffData{2}(:,i));
      [h2(i),p2(i)] = kstest2(indexedCoeffData{1}(:,i),indexedCoeffData{3}(:,i));
      [h3(i),p3(i)] = kstest2(indexedCoeffData{1}(:,i),indexedCoeffData{4}(:,i));
      [h4(i),p4(i)] = kstest2(indexedCoeffData{2}(:,i),indexedCoeffData{3}(:,i));
      [h5(i),p5(i)] = kstest2(indexedCoeffData{2}(:,i),indexedCoeffData{4}(:,i));
      [h6(i),p6(i)] = kstest2(indexedCoeffData{3}(:,i),indexedCoeffData{4}(:,i));


    end
saveas(gca,[figPath filesep 'coefficienthistogram_split_scores' num2str(pcSort) '_100ms'],'svg');     saveas(gca,[figPath filesep 'coefficienthistogram_split_scores' num2str(pcSort) '_100ms'],'png')

%% Plot outputs of KSTEST2
allP = [p1; p2; p3; p4; p5; p6];
allH = [h1; h2; h3; h4; h5; h6];
ylabels = [{'PCon-PInc'}; {'PCon-WCon'}; {'PCon-WInc'}; {'PInc-WCon'}; {'PInc-WInc'}; {'WCon-WInc'}];
xlabels = [{'PC1'}; {'PC2'}; {'PC3'}; {'PC4'}; {'PC5'}];

figure('Units','normalized','Position',[0 0 1 1])
imagesc(allH,[0 1]); colorbar; colormap jet;
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4,'XTick',1:6,'XTickLabel', xlabels, 'YTick',1:6,'YTickLabel',ylabels)
saveas(gca,[figPath filesep 'coefficientKSTEST_HVALS' num2str(pcSort) '_100ms'],'svg');     saveas(gca,[figPath filesep 'coefficientKSTEST_HVALS' num2str(pcSort) '_100ms'],'png')

figure('Units','normalized','Position',[0 0 1 1])
imagesc(allP,[0 1]); colorbar; colormap jet;
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4,'XTick',1:6,'XTickLabel', xlabels, 'YTick',1:6,'YTickLabel',ylabels)
saveas(gca,[figPath filesep 'coefficientKSTEST_PVALS' num2str(pcSort) '_100ms'],'svg');     saveas(gca,[figPath filesep 'coefficientKSTEST_PVALS' num2str(pcSort) '_100ms'],'png')

%% Create environmental variable for regression
cueOnT = cueOn*10;
cueOffT = (cueOn+4)*10;

sipperInT = sipDescent*10;
sipperOutT = sipAscent*10;


numInput = 2; % Cue on, sipper out, sipper in 
envInput = zeros(numInput,size(pCon,2));
envInput(1,cueOnT:cueOffT) = 1;
envInput(2,sipperInT:sipperOutT) = 1;


%% Test Trial-by-Trial PC Loads
pcSort = 1:5;
addpath(genpath('F:\dissDat\restoredScripts'))
for i = 1:length(allPosCoeff)
    for j = 1:length(pcSort)
        meanTrialData{i}{j} = mean(allPosCoeff{i}{j});
        semTrialData{i}{j} = std(allPosCoeff{i}{j})/sqrt(15); % Hard coding 15 trials
    end
end

for i = 1:length(allPosCoeff)
    for j = 1:length(pcSort)
        if i == 1 || i == 3
            meanTrialDiff{i}{j} = mean(allPosCoeff{i+1}{j} - allPosCoeff{i}{j});
        
            semTrialDiff{i}{j} = std(allPosCoeff{i+1}{j} - allPosCoeff{i}{j})/sqrt(15); % Hard coding 15 trials
        end
    end
end
%%
for k = 1:length(pcSort)
    figure('Units','normalized','Position',[0 0 1 1])
    title(['P Congruent vs. P Incongruent: PC Number' num2str(pcSort(k))])
    shadedErrorBar(time,meanTrialData{1}{k},semTrialData{1}{k},'lineprops',{'y-','LineWidth',3})
    hold on
    shadedErrorBar(time,meanTrialData{2}{k},semTrialData{2}{k},'lineprops',{'p-','LineWidth',3})    
    xlabel('Time (s)')
    ylabel('Mean Neural Activity (a.u.)')
    xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
    xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
    xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figPath filesep 'PCon_vs_PInc' num2str(pcSort(k)) '_100ms'],'svg');     saveas(gca,[figPath filesep 'PCon_vs_PInc' num2str(pcSort(k)) '_100ms'],'png')


    figure('Units','normalized','Position',[0 0 1 1])
    title(['W Congruent vs. W Incongruent: PC Number' num2str(pcSort(k))])
    shadedErrorBar(time,meanTrialData{3}{k},semTrialData{3}{k},'lineprops',{'b-','LineWidth',3})
    hold on
    shadedErrorBar(time,meanTrialData{4}{k},semTrialData{4}{k},'lineprops',{'r-','LineWidth',3})    
    xlabel('Time (s)')
    ylabel('Mean Neural Activity (a.u.)')
    xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
    xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
    xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figPath filesep 'WCon_vs_WInc' num2str(pcSort(k)) '_100ms'],'svg');     saveas(gca,[figPath filesep 'WCon_vs_WInc' num2str(pcSort(k)) '_100ms'],'png')

     
%     figure('Units','normalized','Position',[0 0 1 1])
%     title(['P vs W: PC Number' num2str(pcSort(k))])
%     shadedErrorBar(time,meanTrialDiff{1}{k},semTrialDiff{1}{k},'lineprops',{'y-','LineWidth',3})
%     hold on
%     shadedErrorBar(time,meanTrialDiff{3}{k},semTrialDiff{3}{k},'lineprops',{'b-','LineWidth',3})    
%     xlabel('Time (s)')
%     ylabel('Mean Neural Activity (a.u.)')
%     xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
%     xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
%     xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
%     set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
%     saveas(gca,[figPath filesep 'P_vs_W_Difference_shadedError' num2str(pcSort(k)) '_100ms'],'svg');     saveas(gca,[figPath filesep 'P_vs_W_Difference_shadedError' num2str(pcSort(k)) '_100ms'],'png')

end


%% Based on data above, there are no real differences within Genotype (WCon,WInc OR PCon,PInc) coefficient loadings. However, there are consistent differences between genotype loadings
% Therefore, let us analyze P and Wistars separate from one another. 
% This will allow within-session type comparisons but not necessarily allow
% between genotype comparisons except in certain circumstances such as
% comparing the distance metrics between two PC subspaces. 

% A lot of this will repeat code from above. 
% Set indexing based on Genotype
PIdx = startsWith(masterTbl.Strain,'P');
WIdx = startsWith(masterTbl.Strain,'W');
% Index Wistar, P
pMat = cat(3,ephysStruct(PIdx).PSTH);
wMat = cat(3,ephysStruct(WIdx).PSTH);


%% PCA Preprocessing | Separation between P, W groups
timepoint = 'all';

if strcmp(timepoint,'all')
    WConditionsPre = squeeze(mean(wMat(1:15,:,:),1));
    PConditionsPre = squeeze(mean(pMat(1:15,:,:),1));
    time = (1:size(WConditionsPre,1))/10;

elseif strcmp(timepoint,'sipDescent')
    WConditionsPre = squeeze(mean(wMat(1:15,1:sipDescent*10,:),1));
    PConditionsPre = squeeze(mean(pMat(1:15,1:sipDescent*10,:),1)); 
    time = (1:size(WConditionsPre,1))/10;
end

[wConditionsCoeff,wConditionsScore,~,~,wConditionsExplained] = pca(WConditionsPre);
[pConditionsCoeff,pConditionsScore,~,~,pConditionsExplained] = pca(PConditionsPre);

% Need to find an index of neurons per dataset, be able to sort that index
% into our conditions (W/P; Con/Inc). Should ideally be a 1x54 matrix
% containing N # of neurons.

PneuronIndex = [ephysStruct(PIdx).neuronNum];
WneuronIndex = [ephysStruct(WIdx).neuronNum];

%% Sanity check
if sum(PneuronIndex) == size(pConditionsCoeff,1) && sum(WneuronIndex) == size(wConditionsCoeff,1)
    disp('Matrix sizes match')
else
    disp('Matrix size mismatch, address')
end

%% Assign ID to Neurons
PneuronID = []; WneuronID = [];
for i = 1:length(PneuronIndex)
    PneuronID = [PneuronID ones(1,PneuronIndex(i))*i];
end

for i = 1:length(WneuronIndex)
    WneuronID = [WneuronID ones(1,WneuronIndex(i))*i];
end


%% Sanity check
if length(PneuronID) == size(pConditionsCoeff,1) && length(WneuronID) == size(wConditionsCoeff,1)
    disp('Matrix sizes match')
else
    disp('Matrix size mismatch, address')
end

%% Pull data 
% A lot of extra indexing to get to the point where the neurons are broken
% up by session type within the genotype splits. This is due to the
% structure being made to be read 1:54 rather than splitting it beforehand.
% I think what I did here is adequate however and produces the required
% result. 
pSessionLabel = masterTbl.SessionType(PIdx);
wSessionLabel = masterTbl.SessionType(WIdx);

pConSessionIndex = startsWith(pSessionLabel,'Regular');
wConSessionIndex = startsWith(wSessionLabel,'Regular');

pIncSessionIndex = startsWith(pSessionLabel,'Reversal');
wIncSessionIndex = startsWith(wSessionLabel,'Reversal');

%% Find raw, centered data
conPIdx = find(pConSessionIndex == 1);
incPIdx = find(pIncSessionIndex == 1); 
conWIdx = find(wConSessionIndex == 1);
incWIdx = find(wIncSessionIndex == 1);


WcentData = wConditionsCoeff * wConditionsScore';
PcentData = pConditionsCoeff * pConditionsScore';

% Index 
conPCent = PcentData(any(PneuronID == find(pConSessionIndex == 1)),:);
incPCent = PcentData(any(PneuronID == find(pIncSessionIndex == 1)),:);
conWCent = WcentData(any(WneuronID == find(wConSessionIndex == 1)),:);
incWCent = WcentData(any(WneuronID == find(wIncSessionIndex == 1)),:);
% Coefficient
conPCoeff = pConditionsCoeff(any(PneuronID == find(pConSessionIndex == 1)),:);
incPCoeff = pConditionsCoeff(any(PneuronID == find(pIncSessionIndex == 1)),:);
conWCoeff = wConditionsCoeff(any(WneuronID == find(wConSessionIndex == 1)),:);
incWCoeff = wConditionsCoeff(any(WneuronID == find(wIncSessionIndex == 1)),:);
% Create singular variable
indexedCentData = {conPCent incPCent conWCent incWCent};
indexedCoeffData = {conPCoeff incPCoeff conWCoeff incWCoeff};
cellIndex = {conPIdx incPIdx conWIdx incWIdx};

%%
titles = [{'Congruent P'} {'Incongruent P'} {'Congruent W'} {'Incongruent W'}];
pcSort = 1:5;

for j = 1:length(indexedCentData)
    figure('Units','normalized','Position',[0 0 1 1])
    for i=1:length(pcSort)
        
        [b,k]=sort(indexedCoeffData{j}(:,pcSort(i)));      
        subplot(1,length(pcSort),i);
        if j == 1 || j == 2
        centData = indexedCoeffData{j} * pConditionsScore';
        elseif j == 3 || j == 4
        centData = indexedCoeffData{j} * wConditionsScore';
        end

        imagesc([0 max(time)], [1 length(b)],centData(k,:),[-1 1]);

        xlabel('Time (s)');
        if i==1; ylabel('Neuron # (Sorted by PC)');end
        xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        if j == 1 || j == 2
        title(['PC Num: ' num2str(pcSort(i)) ',ExplVar=' num2str(sum(pConditionsExplained(1:pcSort(i))),'%.1f') '%']);
        end
        if j == 3 || j == 4
        title(['PC Num: ' num2str(pcSort(i)) ',ExplVar=' num2str(sum(wConditionsExplained(1:pcSort(i))),'%.1f') '%']);
        end
        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
        if strcmp(timepoint,'all')
        xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        end
         
    end
    sgtitle(titles{j})
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figPath filesep 'allNeurons_pcProj_prePCAsplitGENO' titles{j} timepoint num2str(pcSort) '_100ms'],'svg');     saveas(gca,[figPath filesep 'allNeurons_pcProj_prePCAsplitGENO' titles{j} timepoint num2str(pcSort) '_100ms'],'png')
end
%%
colors = [{'#EDB120'}, {'#7E2F8E'}, {'#0072BD'},  {'#D95319'} ];
pcSort = 1:5;
% Plot coefficient loadings
figure('Units','normalized','Position',[0 0 1 1])

    for i = 1:length(pcSort)
        subplot(1,max(pcSort),i);
        histogram(indexedCoeffData{1}(:,i),'Normalization','probability','DisplayStyle','stairs','LineWidth',2,'LineStyle','-','EdgeColor',colors{1})
        hold on
        histogram(indexedCoeffData{2}(:,i),'Normalization','probability','DisplayStyle','stairs','LineWidth',2,'LineStyle','-','EdgeColor',colors{2})
        histogram(indexedCoeffData{3}(:,i),'Normalization','probability','DisplayStyle','stairs','LineWidth',2,'LineStyle','-','EdgeColor',colors{3})
        histogram(indexedCoeffData{4}(:,i),'Normalization','probability','DisplayStyle','stairs','LineWidth',2,'LineStyle','-','EdgeColor',colors{4})

        xlabel(['PC ' num2str(i) ' Coefficients'])
        ylabel('Probability')
        xlim([-.2,.2])
        ylim([0 0.45])

        legend([{'PCon'},{'PInc'},{'WCon'},{'WInc'}])
        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)


      [h1(i),p1(i)] = kstest2(indexedCoeffData{1}(:,i),indexedCoeffData{2}(:,i));
      [h2(i),p2(i)] = kstest2(indexedCoeffData{1}(:,i),indexedCoeffData{3}(:,i));
      [h3(i),p3(i)] = kstest2(indexedCoeffData{1}(:,i),indexedCoeffData{4}(:,i));
      [h4(i),p4(i)] = kstest2(indexedCoeffData{2}(:,i),indexedCoeffData{3}(:,i));
      [h5(i),p5(i)] = kstest2(indexedCoeffData{2}(:,i),indexedCoeffData{4}(:,i));
      [h6(i),p6(i)] = kstest2(indexedCoeffData{3}(:,i),indexedCoeffData{4}(:,i));


    end
    saveas(gca,[figPath filesep 'coefficienthistogram_split_scores_prePCAsplitGENO' timepoint num2str(pcSort) '_100ms'],'svg');     saveas(gca,[figPath filesep 'coefficienthistogram_split_scores_prePCAsplitGENO' timepoint num2str(pcSort) '_100ms'],'png')

%% Plot outputs of KSTEST2
allP = [p1; p2; p3; p4; p5; p6];
allH = [h1; h2; h3; h4; h5; h6];
ylabels = [{'PCon-PInc'}; {'PCon-WCon'}; {'PCon-WInc'}; {'PInc-WCon'}; {'PInc-WInc'}; {'WCon-WInc'}];
xlabels = [{'PC1'}; {'PC2'}; {'PC3'}; {'PC4'}; {'PC5'}];

figure('Units','normalized','Position',[0 0 1 1])
imagesc(allH,[0 1]); colorbar; colormap jet;
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4,'XTick',1:6,'XTickLabel', xlabels, 'YTick',1:6,'YTickLabel',ylabels)
saveas(gca,[figPath filesep 'coefficientKSTEST_HVALS_prePCAsplitGENO' timepoint num2str(pcSort) '_100ms'],'svg');     saveas(gca,[figPath filesep 'coefficientKSTEST_HVALS_prePCAsplitGENO' timepoint num2str(pcSort) '_100ms'],'png')

figure('Units','normalized','Position',[0 0 1 1])
imagesc(allP,[0 1]); colorbar; colormap jet;
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4,'XTick',1:6,'XTickLabel', xlabels, 'YTick',1:6,'YTickLabel',ylabels)
saveas(gca,[figPath filesep 'coefficientKSTEST_PVALS_prePCAsplitGENO' timepoint num2str(pcSort) '_100ms'],'svg');     saveas(gca,[figPath filesep 'coefficientKSTEST_PVALS_prePCAsplitGENO' timepoint num2str(pcSort) '_100ms'],'png')
%% Low-D PC Representation Based on Separated Matrices
% Parameters
allPCProjections = [];
pcSort = 1:5;
% Obtain data
for i = 1:length(cellIndex)

    if i == 1 || i == 2
        for k = 1:length(pcSort)        
            allPCProjections{i}(:,k) = PConditionsPre(:,any(PneuronID == cellIndex{i})) * indexedCoeffData{i}(:,pcSort(k));
        end
    elseif i == 3 || i == 4
        for k = 1:length(pcSort)
            allPCProjections{i}(:,k) = WConditionsPre(:,any(WneuronID == cellIndex{i})) * indexedCoeffData{i}(:,pcSort(k));
        end
    end

end
%% Plot Top 5 Resulting Matrices for Each Condition
pcSort = 1:5;
colors = [{'#0072BD'} {'#4DBEEE'} {'#77AC30'} {'#D95319'} {'#A2142F'}];

figure('Units','normalized','Position',[0 0 1 1])

for j = 1:length(allPCProjections)
    for i=1:length(pcSort)
        
        subplot(2,2,j);
        plot(time,allPCProjections{j}(:,i),'LineWidth',3,'Color',colors{i})
        if j == 1; legend([{'PC1'},{'PC2'},{'PC3'},{'PC4'},{'PC5'}],'AutoUpdate','off','location','best'); end
        hold on
        xlabel('Time (s)');
        title([titles{j}]);
        ylim([-25 50])
        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',10,'FontWeight','bold','LineWidth',4)
        xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        if strcmp(timepoint,'all')
        xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        end
    end
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
end
saveas(gca,[figPath filesep 'allNeurons_split_top5pcs_prePCAsplitGENO' timepoint num2str(pcSort) '_100ms'],'svg');     saveas(gca,[figPath filesep 'allNeurons_split_top5pcs_prePCAsplitGENO' timepoint num2str(pcSort) '_100ms'],'png')
%% Plot resulting data in a 2D Space
dim = 2;
figure('Units','normalized','Position',[0 0 1 1])
colors = [{'#EDB120'}, {'#7E2F8E'}, {'#0072BD'},  {'#D95319'} ];

for i = 1:length(allPCProjections)
    plot3(allPCProjections{i}(:,1),allPCProjections{i}(:,2),allPCProjections{i}(:,3),'Color',colors{i},'LineWidth',3);
    hold on
    plot3(allPCProjections{i}(1,1),allPCProjections{i}(1,2),allPCProjections{i}(1,3),'o','MarkerSize',10,'MarkerFaceColor',colors{i},'MarkerEdgeColor','k')
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)            
end

saveas(gca,[figPath filesep 'allNeurons_split_pcSpaceProj_prePCAsplitGENO' timepoint '_100ms'],'svg');     saveas(gca,[figPath filesep 'allNeurons_split_pcSpaceProj_prePCAsplitGENO' timepoint '_100ms'],'png')

%% Calculate and Plot the Distance Between PCs 
pCon3D = allPCProjections{1}(:,1:3); wCon3D = allPCProjections{3}(:,1:3);
pInc3D = allPCProjections{2}(:,1:3); wInc3D = allPCProjections{4}(:,1:3);

pDistance = diag(pdist2(pCon3D,pInc3D,'mahalanobis'));
wDistance = diag(pdist2(wCon3D,wInc3D,'mahalanobis'));

figure('Units','normalized','Position',[0 0 1 1])
plot(time,pDistance,'Color',colors{1},'LineWidth',3)
hold on
plot(time,wDistance,'Color',colors{3},'LineWidth',3)
xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
if strcmp(timepoint,'all')
xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
end
xlabel('Time (s)')
ylabel('Mahalanobis Distance')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,[figPath filesep 'allNeurons_split_3DpcSpaceDistance_prePCAsplitGENO' timepoint '_100ms'],'svg');     saveas(gca,[figPath filesep 'allNeurons_split_3DpcSpaceDistance_prePCAsplitGENO' timepoint '_100ms'],'png')


%% Plot Top 3 PCs for Wistars and P Rats Separately. Congruent and Incongruent Overlaid.
% From that, then plot the 3D Space, Distance, Velocity 
pcSort = 1:3;
colors = [{'#EDB120'}, {'#7E2F8E'}, {'#0072BD'},  {'#D95319'} ];
pcLabel = [{'PC1 Score'}, {'PC2 Score'}, {'PC3 Score'}];
%% Wistars
figure('Units','normalized','Position',[0 0 1 1])
sgtitle('Wistars')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',12,'FontWeight','bold','LineWidth',4)

for j = 3:4
    for i=1:length(pcSort)
        k = 1:2:5;
        subplot(7,2,k(i));
        plot(time,allPCProjections{j}(:,i),'LineWidth',3,'Color',colors{j})
        if i == 1; legend([{'Congruent'},{'Incongruent'}],'AutoUpdate','off','location','best'); end
        hold on
        xlabel('Time (s)');
        ylabel(pcLabel{i})
        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',10,'FontWeight','bold','LineWidth',4)
        xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        if strcmp(timepoint,'all')
        xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        end
        xlim([min(time) max(time)])

    end
    subplot(7,2,2:2:6)
    plot3(allPCProjections{j}(:,1),allPCProjections{j}(:,2),allPCProjections{j}(:,3),'LineWidth',3,'Color',colors{j});
    hold on
    plot3(allPCProjections{j}(1,1),allPCProjections{j}(1,2),allPCProjections{j}(1,3),'ko','MarkerSize',10,'MarkerFaceColor','k');
    % Projection
    plot3(allPCProjections{j}(:,1),allPCProjections{j}(:,2),ones(size(allPCProjections{j}(:,3)))*min(allPCProjections{j}(:,3)),'LineStyle','--','Color',colors{j},'LineWidth',3)
    xlabel('PC1')
    ylabel('PC2')
    zlabel('PC3')
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',12,'FontWeight','bold','LineWidth',4)
    
    % Distance Subplot
    subplot(7,2,7:8)
    plot(time,wDistance,'k-','LineWidth',3)
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',12,'FontWeight','bold','LineWidth',4)
        xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        if strcmp(timepoint,'all')
        xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        end
    xlabel('Time (s)')
    ylabel('Distance')
    xlim([min(time) max(time)])

    % Velocity Subplot
    subplot(7,2,9:10)
    plot(time(2:end),diff(allPCProjections{j}(:,1)),'Color',colors{j},'LineWidth',3)
    hold on
    ylabel('PC1 Velocity')
    xlabel('Time (s)')
        xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        if strcmp(timepoint,'all')
        xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        end
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',12,'FontWeight','bold','LineWidth',4)
    xlim([min(time) max(time)])

    subplot(7,2,11:12)
    plot(time(2:end),diff(allPCProjections{j}(:,2)),'Color',colors{j},'LineWidth',3)
    hold on
    ylabel('PC2 Velocity')
    xlabel('Time (s)')
        xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        if strcmp(timepoint,'all')
        xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        end
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',12,'FontWeight','bold','LineWidth',4)
    xlim([min(time) max(time)])

    subplot(7,2,13:14)
    plot(time(2:end),diff(allPCProjections{j}(:,3)),'Color',colors{j},'LineWidth',3)
    hold on
    ylabel('PC3 Velocity')
    xlabel('Time (s)')
        xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        if strcmp(timepoint,'all')
        xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        end
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',12,'FontWeight','bold','LineWidth',4)
    xlim([min(time) max(time)])



end
saveas(gca,[figPath filesep 'W_allPCInfo_allT' timepoint num2str(pcSort) '_100ms'],'svg');     saveas(gca,[figPath filesep 'W_allPCInfo_allT' timepoint num2str(pcSort) '_100ms'],'png')
print('-painters','-dsvg',[figPath filesep 'W_allPCInfo_allT' timepoint num2str(pcSort) '_100ms'])

%% P rats
figure('Units','normalized','Position',[0 0 1 1])
sgtitle('P Rats')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',12,'FontWeight','bold','LineWidth',4)

for j = 1:2
    for i=1:length(pcSort)
        k = 1:2:5;
        subplot(7,2,k(i));
        plot(time,allPCProjections{j}(:,i),'LineWidth',3,'Color',colors{j})
        if i == 1; legend([{'Congruent'},{'Incongruent'}],'AutoUpdate','off','location','best'); end
        hold on
        xlabel('Time (s)');
        ylabel(pcLabel{i})
        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',10,'FontWeight','bold','LineWidth',4)
        xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        if strcmp(timepoint,'all')
        xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        end
        xlim([min(time) max(time)])

    end
    subplot(7,2,2:2:6)
    plot3(allPCProjections{j}(:,1),allPCProjections{j}(:,2),allPCProjections{j}(:,3),'LineWidth',3,'Color',colors{j});
    hold on
    plot3(allPCProjections{j}(1,1),allPCProjections{j}(1,2),allPCProjections{j}(1,3),'ko','MarkerSize',10,'MarkerFaceColor','k');
    % Projection
    plot3(allPCProjections{j}(:,1),allPCProjections{j}(:,2),ones(size(allPCProjections{j}(:,3)))*min(allPCProjections{j}(:,3)),'LineStyle','--','Color',colors{j},'LineWidth',3)
    xlabel('PC1')
    ylabel('PC2')
    zlabel('PC3')
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',12,'FontWeight','bold','LineWidth',4)
    
    % Distance Subplot
    subplot(7,2,7:8)
    plot(time,pDistance,'k-','LineWidth',3)
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',12,'FontWeight','bold','LineWidth',4)
        xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        if strcmp(timepoint,'all')
        xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        end
    xlabel('Time (s)')
    ylabel('Distance')
    xlim([min(time) max(time)])


    % Velocity Subplot
    subplot(7,2,9:10)
    plot(time(2:end),diff(allPCProjections{j}(:,1)),'Color',colors{j},'LineWidth',3)
    hold on
    ylabel('PC1 Velocity')
    xlabel('Time (s)')
        xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        if strcmp(timepoint,'all')
        xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        end
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',12,'FontWeight','bold','LineWidth',4)
    xlim([min(time) max(time)])
   
    subplot(7,2,11:12)
    plot(time(2:end),diff(allPCProjections{j}(:,2)),'Color',colors{j},'LineWidth',3)
    hold on
    ylabel('PC2 Velocity')
    xlabel('Time (s)')
        xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        if strcmp(timepoint,'all')
        xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        end
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',12,'FontWeight','bold','LineWidth',4)
    xlim([min(time) max(time)])

    subplot(7,2,13:14)
    plot(time(2:end),diff(allPCProjections{j}(:,3)),'Color',colors{j},'LineWidth',3)
    hold on
    ylabel('PC3 Velocity')
    xlabel('Time (s)')
        xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        if strcmp(timepoint,'all')
        xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        end
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',12,'FontWeight','bold','LineWidth',4)
    xlim([min(time) max(time)])



end
saveas(gca,[figPath filesep 'P_allPCInfo_allT' timepoint num2str(pcSort) '_100ms'],'svg');     saveas(gca,[figPath filesep 'P_allPCInfo_allT' timepoint num2str(pcSort) '_100ms'],'png')
print('-painters','-dsvg',[figPath filesep 'P_allPCInfo_allT' timepoint num2str(pcSort) '_100ms'])
%% 
%% 05/24/2022
% Separating PCA based on Ps and Ws seems to work the best. Upcoming
% adventures prior to hard pivoting to RSA include quantifying how
% predictive each PC is of different epochs / events in the task.
% Additionally, reconstructing PSTHs to be centered on the first 'contact'
% with the sipper - correct or incorrect - will be the next essential thing
% to finish. Based on these two points, I will have an incredibly solid
% start on the Dissertation as a whole and a compelling story for the RSA
% poster. 

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
% Pull ephys and trial data, create PSTH
    pBin = 0.1; % 33 ms bins for fr variable, match video framerate
    PSTH_Approach = [];
    if ~isfield(ephysStruct,'PSTH_Approach')
        
        for i = 1:length(trlStruct)
    
            % Load trial times from trlStruct
            trialTimes = trlStruct(i).trialTimes(1:48);
            [trialTimes, trialTimesIdx] = sort(trialTimes);
            approachVar = min(trlStruct(i).approach(:,1:2,2),[],2);
            sortedApproachVar = approachVar(trialTimesIdx);

            if startsWith(masterTbl.SessionType(i),'Regular')            
                correctIndex = trlStruct(i).approach(1:48,1,1) == 1 & trlStruct(i).approach(1:48,2,1) == 0;
                incorrectIndex = trlStruct(i).approach(1:48,2,1) == 1;
            elseif startsWith(masterTbl.SessionType(i),'Reversal')
                correctIndex = trlStruct(i).approach(1:48,1,1) == 0 & trlStruct(i).approach(1:48,2,1) == 1;
                incorrectIndex = trlStruct(i).approach(1:48,1,1) == 1;
            end

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
            PSTH_Approach = nan(length(trialTimes),41,size(fr,2));
        for k = 1:length(trialTimes)
            [~,tMin] = min(abs((trialTimes(k) - 1) - dlcTimestamps));
            [~,tMax] = min(abs((trialTimes(k) + 3) - dlcTimestamps)); 

            [~,tMin_fr] = min(abs((dlcTimestamps(tMin)) - frTime));
            [~,tMax_fr] = min(abs((dlcTimestamps(tMax)) - frTime));
            
            trlSize(k) = size(fr(tMin_fr:tMax_fr,:),1);

            if size(fr(tMin_fr:tMax_fr,:),1) == 41
                PSTH_Approach(k,:,:) = fr(tMin_fr:tMax_fr,:);
            elseif size(fr(tMin_fr:tMax_fr,:),1) == 42
                PSTH_Approach(k,:,:) = fr(tMin_fr:tMax_fr-1,:);
            elseif size(fr(tMin_fr:tMax_fr,:),1) == 40
                PSTH_Approach(k,:,:) = fr(tMin_fr:tMax_fr+1,:);          
            end
            
                
               
        end

        

        % Remove low firing rate neurons. Not sure what the best value is
        % here.
        minFR = 0.1;
        if ~isempty(PSTH_Approach)            
            meanFR = squeeze(nanmean(nanmean(PSTH_Approach,1),2)); 
            lowFR = meanFR > minFR; % Find values greater than minFR
            PSTH_Approach = PSTH_Approach(:,:,lowFR);
        end
        
        if length(size(PSTH_Approach)) == 2
            PSTH_Approach = [];
        end
        
        if ~isempty(PSTH_Approach)
            ephysStruct(i).PSTH_Approach = PSTH_Approach(1:15,:,:);
            ephysStruct(i).neuronNumApproach = size(PSTH_Approach,3);
        else
            ephysStruct(i).neuronNumApproach = 0;
        end
        ephysStruct(i).maxFR = max(fr);

        if ~isempty(PSTH_Approach)
            ephysStruct(i).PSTH_Correct = nan(size(PSTH_Approach));
            ephysStruct(i).PSTH_Incorrect = nan(size(PSTH_Approach));

            ephysStruct(i).PSTH_Correct = PSTH_Approach(correctIndex(1:15),:,:);
            ephysStruct(i).PSTH_Incorrect = PSTH_Approach(incorrectIndex(1:15),:,:);                    
        else
            ephysStruct(i).PSTH_Correct = [];
            ephysStruct(i).PSTH_Incorrect = [];

        end
        ephysStruct(i).PSTH_ApproachMean = squeeze(nanmean(PSTH_Approach,1));
        ephysStruct(i).correctIndex = correctIndex;
        ephysStruct(i).incorrectIndex = incorrectIndex;
        ephysStruct(i).PSTH_CorrectMean = squeeze(nanmean(ephysStruct(i).PSTH_Correct,1));
        ephysStruct(i).PSTH_IncorrectMean = squeeze(nanmean(ephysStruct(i).PSTH_Incorrect,1));
        PSTH_Approach = [];

    end

end

%%
%% Based on data above, there are no real differences within Genotype (WCon,WInc OR PCon,PInc) coefficient loadings. However, there are consistent differences between genotype loadings
% Therefore, let us analyze P and Wistars separate from one another. 
% This will allow within-session type comparisons but not necessarily allow
% between genotype comparisons except in certain circumstances such as
% comparing the distance metrics between two PC subspaces. 

% A lot of this will repeat code from above. 
% Set indexing based on Genotype
pMat = [];
wMat = [];
%
PIdx = startsWith(masterTbl.Strain,'P');
WIdx = startsWith(masterTbl.Strain,'W');
% Index Wistar, P
psthFlag = 'incorrect';
if strcmp(psthFlag,'all')
    pMat = [ephysStruct(PIdx).PSTH_ApproachMean];
    wMat = [ephysStruct(WIdx).PSTH_ApproachMean];
elseif strcmp(psthFlag,'correct')
    pMat = [ephysStruct(PIdx).PSTH_CorrectMean];
    wMat = [ephysStruct(WIdx).PSTH_CorrectMean];
elseif strcmp(psthFlag,'incorrect')
    pMat = [ephysStruct(PIdx).PSTH_IncorrectMean];
    wMat = [ephysStruct(WIdx).PSTH_IncorrectMean];
end

%% PCA Preprocessing | Separation between P, W groups
timepoint = 'all';

if strcmp(timepoint,'all')
    WConditionsPre = wMat;
    PConditionsPre = pMat;
    time = (1:size(WConditionsPre,1))/10;

elseif strcmp(timepoint,'sipDescent')
    WConditionsPre = wMat(1:sipDescent*10,:);
    PConditionsPre = pMat(1:sipDescent*10,:); 
    time = (1:size(WConditionsPre,1))/10;
end

[wConditionsCoeff,wConditionsScore,~,~,wConditionsExplained] = pca(WConditionsPre);
[pConditionsCoeff,pConditionsScore,~,~,pConditionsExplained] = pca(PConditionsPre);

% Need to find an index of neurons per dataset, be able to sort that index
% into our conditions (W/P; Con/Inc). Should ideally be a 1x54 matrix
% containing N # of neurons.

PneuronIndex = [ephysStruct(PIdx).neuronNumApproach];
WneuronIndex = [ephysStruct(WIdx).neuronNumApproach];

%% Sanity check
if sum(PneuronIndex) == size(pConditionsCoeff,1) && sum(WneuronIndex) == size(wConditionsCoeff,1)
    disp('Matrix sizes match')
else
    disp('Matrix size mismatch, address')
end

%% Assign ID to Neurons
PneuronID = []; WneuronID = [];
for i = 1:length(PneuronIndex)
    PneuronID = [PneuronID ones(1,PneuronIndex(i))*i];
end

for i = 1:length(WneuronIndex)
    WneuronID = [WneuronID ones(1,WneuronIndex(i))*i];
end


%% Sanity check
if length(PneuronID) == size(pConditionsCoeff,1) && length(WneuronID) == size(wConditionsCoeff,1)
    disp('Matrix sizes match')
else
    disp('Matrix size mismatch, address')
end

%% Pull data 
% A lot of extra indexing to get to the point where the neurons are broken
% up by session type within the genotype splits. This is due to the
% structure being made to be read 1:54 rather than splitting it beforehand.
% I think what I did here is adequate however and produces the required
% result. 
pSessionLabel = masterTbl.SessionType(PIdx);
wSessionLabel = masterTbl.SessionType(WIdx);

pConSessionIndex = startsWith(pSessionLabel,'Regular');
wConSessionIndex = startsWith(wSessionLabel,'Regular');

pIncSessionIndex = startsWith(pSessionLabel,'Reversal');
wIncSessionIndex = startsWith(wSessionLabel,'Reversal');

%% Find raw, centered data
conPIdx = find(pConSessionIndex == 1);
incPIdx = find(pIncSessionIndex == 1); 
conWIdx = find(wConSessionIndex == 1);
incWIdx = find(wIncSessionIndex == 1);


WcentData = wConditionsCoeff * wConditionsScore';
PcentData = pConditionsCoeff * pConditionsScore';

% Index 
conPCent = PcentData(any(PneuronID == find(pConSessionIndex == 1)),:);
incPCent = PcentData(any(PneuronID == find(pIncSessionIndex == 1)),:);
conWCent = WcentData(any(WneuronID == find(wConSessionIndex == 1)),:);
incWCent = WcentData(any(WneuronID == find(wIncSessionIndex == 1)),:);
% Coefficient
conPCoeff = pConditionsCoeff(any(PneuronID == find(pConSessionIndex == 1)),:);
incPCoeff = pConditionsCoeff(any(PneuronID == find(pIncSessionIndex == 1)),:);
conWCoeff = wConditionsCoeff(any(WneuronID == find(wConSessionIndex == 1)),:);
incWCoeff = wConditionsCoeff(any(WneuronID == find(wIncSessionIndex == 1)),:);
% Create singular variable
indexedCentData = {conPCent incPCent conWCent incWCent};
indexedCoeffData = {conPCoeff incPCoeff conWCoeff incWCoeff};
cellIndex = {conPIdx incPIdx conWIdx incWIdx};

%%
titles = [{'Congruent P'} {'Incongruent P'} {'Congruent W'} {'Incongruent W'}];
pcSort = 1:5;

for j = 1:length(indexedCentData)
    figure('Units','normalized','Position',[0 0 1 1])
    for i=1:length(pcSort)
        
        [b,k]=sort(indexedCoeffData{j}(:,pcSort(i)));      
        subplot(1,length(pcSort),i);
        if j == 1 || j == 2
        centData = indexedCoeffData{j} * pConditionsScore';
        elseif j == 3 || j == 4
        centData = indexedCoeffData{j} * wConditionsScore';
        end

        imagesc([0 max(time)], [1 length(b)],centData(k,:),[-1 1]);

        xlabel('Time (s)');
        if i==1; ylabel('Neuron # (Sorted by PC)');end
        xline(cueOn,'k--','Approach','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        if j == 1 || j == 2
        title(['PC Num: ' num2str(pcSort(i)) ',ExplVar=' num2str(sum(pConditionsExplained(1:pcSort(i))),'%.1f') '%']);
        end
        if j == 3 || j == 4
        title(['PC Num: ' num2str(pcSort(i)) ',ExplVar=' num2str(sum(wConditionsExplained(1:pcSort(i))),'%.1f') '%']);
        end
        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
        if strcmp(timepoint,'all')
%         xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
        end
         
    end
    sgtitle(titles{j})
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figPath filesep 'allNeurons_pcProj_prePCAsplitGENO' psthFlag titles{j} timepoint num2str(pcSort) '_100ms'],'svg');     saveas(gca,[figPath filesep 'allNeurons_pcProj_prePCAsplitGENO' psthFlag titles{j} timepoint num2str(pcSort) '_100ms'],'png')
end
%%
colors = [{'#EDB120'}, {'#7E2F8E'}, {'#0072BD'},  {'#D95319'} ];
pcSort = 1:5;
% Plot coefficient loadings
figure('Units','normalized','Position',[0 0 1 1])

    for i = 1:length(pcSort)
        subplot(1,max(pcSort),i);
        histogram(indexedCoeffData{1}(:,i),'Normalization','probability','DisplayStyle','stairs','LineWidth',2,'LineStyle','-','EdgeColor',colors{1})
        hold on
        histogram(indexedCoeffData{2}(:,i),'Normalization','probability','DisplayStyle','stairs','LineWidth',2,'LineStyle','-','EdgeColor',colors{2})
        histogram(indexedCoeffData{3}(:,i),'Normalization','probability','DisplayStyle','stairs','LineWidth',2,'LineStyle','-','EdgeColor',colors{3})
        histogram(indexedCoeffData{4}(:,i),'Normalization','probability','DisplayStyle','stairs','LineWidth',2,'LineStyle','-','EdgeColor',colors{4})

        xlabel(['PC ' num2str(i) ' Coefficients'])
        ylabel('Probability')
        xlim([-.2,.2])
        ylim([0 0.45])

        legend([{'PCon'},{'PInc'},{'WCon'},{'WInc'}])
        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)


      [h1(i),p1(i)] = kstest2(indexedCoeffData{1}(:,i),indexedCoeffData{2}(:,i));
      [h2(i),p2(i)] = kstest2(indexedCoeffData{1}(:,i),indexedCoeffData{3}(:,i));
      [h3(i),p3(i)] = kstest2(indexedCoeffData{1}(:,i),indexedCoeffData{4}(:,i));
      [h4(i),p4(i)] = kstest2(indexedCoeffData{2}(:,i),indexedCoeffData{3}(:,i));
      [h5(i),p5(i)] = kstest2(indexedCoeffData{2}(:,i),indexedCoeffData{4}(:,i));
      [h6(i),p6(i)] = kstest2(indexedCoeffData{3}(:,i),indexedCoeffData{4}(:,i));


    end
    saveas(gca,[figPath filesep 'coefficienthistogram_split_scores_prePCAsplitGENO' psthFlag timepoint num2str(pcSort) '_100ms'],'svg');     saveas(gca,[figPath filesep 'coefficienthistogram_split_scores_prePCAsplitGENO' psthFlag timepoint num2str(pcSort) '_100ms'],'png')

%% Plot outputs of KSTEST2
allP = [p1; p2; p3; p4; p5; p6];
allH = [h1; h2; h3; h4; h5; h6];
ylabels = [{'PCon-PInc'}; {'PCon-WCon'}; {'PCon-WInc'}; {'PInc-WCon'}; {'PInc-WInc'}; {'WCon-WInc'}];
xlabels = [{'PC1'}; {'PC2'}; {'PC3'}; {'PC4'}; {'PC5'}];

figure('Units','normalized','Position',[0 0 1 1])
imagesc(allH,[0 1]); colorbar; colormap jet;
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4,'XTick',1:6,'XTickLabel', xlabels, 'YTick',1:6,'YTickLabel',ylabels)
saveas(gca,[figPath filesep 'coefficientKSTEST_HVALS_prePCAsplitGENO' psthFlag timepoint num2str(pcSort) '_100ms'],'svg');     saveas(gca,[figPath filesep 'coefficientKSTEST_HVALS_prePCAsplitGENO' psthFlag timepoint num2str(pcSort) '_100ms'],'png')

figure('Units','normalized','Position',[0 0 1 1])
imagesc(allP,[0 1]); colorbar; colormap jet;
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4,'XTick',1:6,'XTickLabel', xlabels, 'YTick',1:6,'YTickLabel',ylabels)
saveas(gca,[figPath filesep 'coefficientKSTEST_PVALS_prePCAsplitGENO' psthFlag timepoint num2str(pcSort) '_100ms'],'svg');     saveas(gca,[figPath filesep 'coefficientKSTEST_PVALS_prePCAsplitGENO' psthFlag timepoint num2str(pcSort) '_100ms'],'png')
%% Low-D PC Representation Based on Separated Matrices
% Parameters
allPCProjections = [];
pcSort = 1:5;
% Obtain data
for i = 1:length(cellIndex)

    if i == 1 || i == 2
        for k = 1:length(pcSort)        
            allPCProjections{i}(:,k) = PConditionsPre(:,any(PneuronID == cellIndex{i})) * indexedCoeffData{i}(:,pcSort(k));
        end
    elseif i == 3 || i == 4
        for k = 1:length(pcSort)
            allPCProjections{i}(:,k) = WConditionsPre(:,any(WneuronID == cellIndex{i})) * indexedCoeffData{i}(:,pcSort(k));
        end
    end

end
%% Plot Top 5 Resulting Matrices for Each Condition
pcSort = 1:5;
colors = [{'#0072BD'} {'#4DBEEE'} {'#77AC30'} {'#D95319'} {'#A2142F'}];

figure('Units','normalized','Position',[0 0 1 1])

for j = 1:length(allPCProjections)
    for i=1:length(pcSort)
        
        subplot(2,2,j);
        plot(time,allPCProjections{j}(:,i),'LineWidth',3,'Color',colors{i})
        if j == 1; legend([{'PC1'},{'PC2'},{'PC3'},{'PC4'},{'PC5'}],'AutoUpdate','off','location','best'); end
        hold on
        xlabel('Time (s)');
        title([titles{j}]);
        ylim([-25 50])
        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',10,'FontWeight','bold','LineWidth',4)
        xline(cueOn,'k--','Approach','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
    end
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
end
saveas(gca,[figPath filesep 'allNeurons_split_top5pcs_prePCAsplitGENO' psthFlag timepoint num2str(pcSort) '_100ms'],'svg');     saveas(gca,[figPath filesep 'allNeurons_split_top5pcs_prePCAsplitGENO' psthFlag timepoint num2str(pcSort) '_100ms'],'png')
%% Plot resulting data in a 2D Space
dim = 2;
figure('Units','normalized','Position',[0 0 1 1])
colors = [{'#EDB120'}, {'#7E2F8E'}, {'#0072BD'},  {'#D95319'} ];

for i = 1:length(allPCProjections)
    plot3(allPCProjections{i}(:,1),allPCProjections{i}(:,2),allPCProjections{i}(:,3),'Color',colors{i},'LineWidth',3);
    hold on
    plot3(allPCProjections{i}(1,1),allPCProjections{i}(1,2),allPCProjections{i}(1,3),'o','MarkerSize',10,'MarkerFaceColor',colors{i},'MarkerEdgeColor','k')
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)            
end

saveas(gca,[figPath filesep 'allNeurons_split_pcSpaceProj_prePCAsplitGENO' psthFlag timepoint '_100ms'],'svg');     saveas(gca,[figPath filesep 'allNeurons_split_pcSpaceProj_prePCAsplitGENO_COR' psthFlag timepoint '_100ms'],'png')

%% Calculate and Plot the Distance Between PCs 
pCon3D = allPCProjections{1}(:,1:3); wCon3D = allPCProjections{3}(:,1:3);
pInc3D = allPCProjections{2}(:,1:3); wInc3D = allPCProjections{4}(:,1:3);

pDistance = diag(pdist2(pCon3D,pInc3D,'mahalanobis'));
wDistance = diag(pdist2(wCon3D,wInc3D,'mahalanobis'));

figure('Units','normalized','Position',[0 0 1 1])
plot(time,pDistance,'Color',colors{1},'LineWidth',3)
hold on
plot(time,wDistance,'Color',colors{3},'LineWidth',3)
xline(cueOn,'k--','Approach','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
if strcmp(timepoint,'all')
end
xlabel('Time (s)')
ylabel('Mahalanobis Distance')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,[figPath filesep 'allNeurons_split_3DpcSpaceDistance_prePCAsplitGENO' psthFlag timepoint '_100ms'],'svg');     saveas(gca,[figPath filesep 'allNeurons_split_3DpcSpaceDistance_prePCAsplitGENO' psthFlag timepoint '_100ms'],'png')

%
%% Plot Top 3 PCs for Wistars and P Rats Separately. Congruent and Incongruent Overlaid.
% From that, then plot the 3D Space, Distance, Velocity 
pcSort = 1:3;
colors = [{'#EDB120'}, {'#7E2F8E'}, {'#0072BD'},  {'#D95319'} ];
pcLabel = [{'PC1 Score'}, {'PC2 Score'}, {'PC3 Score'}];
%% Wistars
figure('Units','normalized','Position',[0 0 1 1])
sgtitle('Wistars')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',12,'FontWeight','bold','LineWidth',4)

for j = 3:4
    for i=1:length(pcSort)
        k = 1:2:5;
        subplot(7,2,k(i));
        plot(time,allPCProjections{j}(:,i),'LineWidth',3,'Color',colors{j})
        if i == 1; legend([{'Congruent'},{'Incongruent'}],'AutoUpdate','off','location','best'); end
        hold on
        xlabel('Time (s)');
        ylabel(pcLabel{i})
        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',10,'FontWeight','bold','LineWidth',4)
        xline(cueOn,'k--','Approach','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');

        xlim([min(time) max(time)])

    end
    subplot(7,2,2:2:6)
    plot3(allPCProjections{j}(:,1),allPCProjections{j}(:,2),allPCProjections{j}(:,3),'LineWidth',3,'Color',colors{j});
    hold on
    plot3(allPCProjections{j}(1,1),allPCProjections{j}(1,2),allPCProjections{j}(1,3),'ko','MarkerSize',10,'MarkerFaceColor','k');
    % Projection
    plot3(allPCProjections{j}(:,1),allPCProjections{j}(:,2),ones(size(allPCProjections{j}(:,3)))*min(allPCProjections{j}(:,3)),'LineStyle','--','Color',colors{j},'LineWidth',3)
    xlabel('PC1')
    ylabel('PC2')
    zlabel('PC3')
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',12,'FontWeight','bold','LineWidth',4)
    
    % Distance Subplot
    subplot(7,2,7:8)
    plot(time,wDistance,'k-','LineWidth',3)
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',12,'FontWeight','bold','LineWidth',4)
    xline(cueOn,'k--','Approach','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');

    xlabel('Time (s)')
    ylabel('Distance')
    xlim([min(time) max(time)])

    % Velocity Subplot
    subplot(7,2,9:10)
    plot(time(2:end),diff(allPCProjections{j}(:,1)),'Color',colors{j},'LineWidth',3)
    hold on
    ylabel('PC1 Velocity')
    xlabel('Time (s)')
    xline(cueOn,'k--','Approach','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');

    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',12,'FontWeight','bold','LineWidth',4)
    xlim([min(time) max(time)])

    subplot(7,2,11:12)
    plot(time(2:end),diff(allPCProjections{j}(:,2)),'Color',colors{j},'LineWidth',3)
    hold on
    ylabel('PC2 Velocity')
    xlabel('Time (s)')
    xline(cueOn,'k--','Approach','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');

    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',12,'FontWeight','bold','LineWidth',4)
    xlim([min(time) max(time)])

    subplot(7,2,13:14)
    plot(time(2:end),diff(allPCProjections{j}(:,3)),'Color',colors{j},'LineWidth',3)
    hold on
    ylabel('PC3 Velocity')
    xlabel('Time (s)')
    xline(cueOn,'k--','Approach','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');

    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',12,'FontWeight','bold','LineWidth',4)
    xlim([min(time) max(time)])



end
saveas(gca,[figPath filesep 'W_allPCInfo_approachCentered' psthFlag num2str(pcSort) '_100ms'],'svg');     saveas(gca,[figPath filesep 'W_allPCInfo_approachCentered' psthFlag num2str(pcSort) '_100ms'],'png')

%% P rats
figure('Units','normalized','Position',[0 0 1 1])
sgtitle('P Rats')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',12,'FontWeight','bold','LineWidth',4)

for j = 1:2
    for i=1:length(pcSort)
        k = 1:2:5;
        subplot(7,2,k(i));
        plot(time,allPCProjections{j}(:,i),'LineWidth',3,'Color',colors{j})
        if i == 1; legend([{'Congruent'},{'Incongruent'}],'AutoUpdate','off','location','best'); end
        hold on
        xlabel('Time (s)');
        ylabel(pcLabel{i})
        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',10,'FontWeight','bold','LineWidth',4)
        xline(cueOn,'k--','Approach','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');

        xlim([min(time) max(time)])

    end
    subplot(7,2,2:2:6)
    plot3(allPCProjections{j}(:,1),allPCProjections{j}(:,2),allPCProjections{j}(:,3),'LineWidth',3,'Color',colors{j});
    hold on
    plot3(allPCProjections{j}(1,1),allPCProjections{j}(1,2),allPCProjections{j}(1,3),'ko','MarkerSize',10,'MarkerFaceColor','k');
    % Projection
    plot3(allPCProjections{j}(:,1),allPCProjections{j}(:,2),ones(size(allPCProjections{j}(:,3)))*min(allPCProjections{j}(:,3)),'LineStyle','--','Color',colors{j},'LineWidth',3)
    xlabel('PC1')
    ylabel('PC2')
    zlabel('PC3')
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',12,'FontWeight','bold','LineWidth',4)
    
    % Distance Subplot
    subplot(7,2,7:8)
    plot(time,pDistance,'k-','LineWidth',3)
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',12,'FontWeight','bold','LineWidth',4)
    xline(cueOn,'k--','Approach','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');

    xlabel('Time (s)')
    ylabel('Distance')
    xlim([min(time) max(time)])


    % Velocity Subplot
    subplot(7,2,9:10)
    plot(time(2:end),diff(allPCProjections{j}(:,1)),'Color',colors{j},'LineWidth',3)
    hold on
    ylabel('PC1 Velocity')
    xlabel('Time (s)')
    xline(cueOn,'k--','Approach','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');

    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',12,'FontWeight','bold','LineWidth',4)
    xlim([min(time) max(time)])
   
    subplot(7,2,11:12)
    plot(time(2:end),diff(allPCProjections{j}(:,2)),'Color',colors{j},'LineWidth',3)
    hold on
    ylabel('PC2 Velocity')
    xlabel('Time (s)')
    xline(cueOn,'k--','Approach','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');

    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',12,'FontWeight','bold','LineWidth',4)
    xlim([min(time) max(time)])

    subplot(7,2,13:14)
    plot(time(2:end),diff(allPCProjections{j}(:,3)),'Color',colors{j},'LineWidth',3)
    hold on
    ylabel('PC3 Velocity')
    xlabel('Time (s)')
    xline(cueOn,'k--','Approach','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');

    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',12,'FontWeight','bold','LineWidth',4)
    xlim([min(time) max(time)])



end
saveas(gca,[figPath filesep 'P_allPCInfo_approachCentered' psthFlag num2str(pcSort) '_100ms'],'svg');     saveas(gca,[figPath filesep 'P_allPCInfo_approachCentered' psthFlag num2str(pcSort) '_100ms'],'png')

%% 05/31/2022
% Back to analyzing approach differences between correct/incorrect choices
% This time, project data onto trial-by-trial analyses
% Then, average data and compare
%%
pCon = cat(3,ephysStruct(regPIdx).PSTH_Approach);
pInc = cat(3,ephysStruct(revPIdx).PSTH_Approach);
wCon = cat(3,ephysStruct(regWIdx).PSTH_Approach);
wInc = cat(3,ephysStruct(revWIdx).PSTH_Approach);
%%
for j = 1:length(indexedCoeffData)
    pcSort = 1:5;
    % Assign variable rawData depending on condition 
    if j == 1
        rawData = pCon;
    elseif j == 2
        rawData = pInc;
    elseif j == 3
        rawData = wCon;
    elseif j == 4
        rawData = wInc;
    end
    % Sort according to top 3 PCs, plot similar to other coefficient
    % figures

    figure('Units','normalized','Position',[0 0 1 1])

    for i = 1:length(pcSort)
        % Index Pos and Neg
        PosCoeff = rawData(:,:,indexedCoeffData{j}(:,i) > 0.01);                  
        NegCoeff = rawData(:,:,indexedCoeffData{j}(:,i) < -0.01);
        % Take Mean
        PosCoeff = (nanmean(PosCoeff,3));
        NegCoeff = (nanmean(NegCoeff,3));
        %
        mPosCoeff{j}(:,i) = nanmean(PosCoeff,1);
        mNegCoeff{j}(:,i) = nanmean(NegCoeff,1);
        stdPosCoeff{j}(:,i) = nanstd(PosCoeff,1)/sqrt(length(15));
        stdNegCoeff{j}(:,i) = nanstd(NegCoeff,1)/sqrt(length(15));
        % All
        allPosCoeff{j}(:,:,i) = PosCoeff; 
        allNegCoeff{j}(:,:,i) = NegCoeff; 

        % Figure
        subplot(2,max(pcSort),i);          % Doubling subplot for pos+neg
        % Positive Plot
        imagesc([0 max(time)], [1 15],PosCoeff, [0 10]);
        colormap jet
        colorbar
        xlabel('Time (s)');
        if i==1; ylabel('Trial #');end
        xline(cueOn,'k--','Approach','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        
        title(['PC Num: ' num2str(pcSort(i))]);
        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
        % Negative Plot

        subplot(2,max(pcSort),i+max(pcSort));          % Doubling subplot for pos+neg

        imagesc([0 max(time)], [1 15],NegCoeff,[0 10]);
        colormap jet
        colorbar
        xlabel('Time (s)');
        if i==1; ylabel('Trial #');end
        xline(cueOn,'k--','Approach','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        
        title(['PC Num: ' num2str(pcSort(i))]);
        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)        
    end

    sgtitle(titles{j})
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figPath filesep 'corincor_coeffNeurons_overTrials' num2str(pcSort) titles{j} '_100ms'],'svg');     saveas(gca,[figPath filesep 'corincor_coeffNeurons_overTrials' num2str(pcSort) titles{j} '_100ms'],'png')

end

% Of the 15 first trials that we have access to, can we see if there are
% meaningful differences in neural activity between correct and incorrect
% choices? 

% A strategy, find the index of correct and incorrect choices. Of those
% choices, take the mean of 

%% 
figure('Units','normalized','Position',[0 0 1 1])
for i = 1:length(mPosCoeff)
    plot3(mPosCoeff{i}(:,1),mPosCoeff{i}(:,2),mPosCoeff{i}(:,3),'LineWidth',3,'Color',colors{i})
    hold on

    plot3(mPosCoeff{i}(1,1),mPosCoeff{i}(1,2),mPosCoeff{i}(1,3),'o','LineWidth',3,'Color',colors{i})

    xlabel('Mean of PC1 Positive Loading Neurons')
    ylabel('Mean of PC2 Positive Loading Neurons')
    zlabel('Mean of PC3 Positive Loading Neurons')
end
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figPath filesep '3DProjection_Approach_POSCOEFF' num2str(pcSort) titles{j} '_100ms'],'svg');     saveas(gca,[figPath filesep '3DProjection_Approach' num2str(pcSort) titles{j} '_100ms'],'png')

figure('Units','normalized','Position',[0 0 1 1])
for i = 1:length(mPosCoeff)
    plot3(mNegCoeff{i}(:,1),mNegCoeff{i}(:,2),mNegCoeff{i}(:,3),'LineWidth',3,'Color',colors{i})
    hold on

    plot3(mNegCoeff{i}(1,1),mNegCoeff{i}(1,2),mNegCoeff{i}(1,3),'o','LineWidth',3,'Color',colors{i})

    xlabel('Mean of PC1 Negative Loading Neurons')
    ylabel('Mean of PC2 Negative Loading Neurons')
    zlabel('Mean of PC3 Negative Loading Neurons')
end
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figPath filesep '3DProjection_Approach_NEGCOEFF' num2str(pcSort) titles{j} '_100ms'],'svg');     saveas(gca,[figPath filesep '3DProjection_Approach' num2str(pcSort) titles{j} '_100ms'],'png')

%% Plot Positive / Negative Loaders. 
% Graph with raw values per each condition. Each graph, 1 PC
for i = 1:length(pcSort)

    figure('Units','normalized','Position',[0 0 1 1])
    shadedErrorBar(time,(mPosCoeff{1}(:,i)), stdPosCoeff{1}(:,i),'lineprops',{'Color',colors{1},'LineWidth',3})
    hold on
    shadedErrorBar(time,(mPosCoeff{2}(:,i)), stdPosCoeff{2}(:,i),'lineprops',{'Color',colors{2},'LineWidth',3})
    shadedErrorBar(time,(mPosCoeff{3}(:,i)), stdPosCoeff{3}(:,i),'lineprops',{'Color',colors{3},'LineWidth',3})
    shadedErrorBar(time,(mPosCoeff{4}(:,i)), stdPosCoeff{4}(:,i),'lineprops',{'Color',colors{4},'LineWidth',3})

end
%% Calculate distances
% 
clear posCoeffDistances negCoeffDistances
numTrials = 15;
method = 'mahalanobis';
    for i = 1:numTrials
        posCoeffDistances{1}(i,:) = diag(pdist2(squeeze(allPosCoeff{1}(i,:,:)),squeeze(allPosCoeff{2}(i,:,:)),method));
        negCoeffDistances{1}(i,:) = diag(pdist2(squeeze(allNegCoeff{1}(i,:,:)),squeeze(allNegCoeff{2}(i,:,:)),method));
        posCoeffDistances{2}(i,:) = diag(pdist2(squeeze(allPosCoeff{3}(i,:,:)),squeeze(allPosCoeff{4}(i,:,:)),method));
        negCoeffDistances{2}(i,:) = diag(pdist2(squeeze(allNegCoeff{3}(i,:,:)),squeeze(allNegCoeff{4}(i,:,:)),method));
    end

%% Plot distances
figure('Units','normalized','Position',[0 0 1 1])

shadedErrorBar(time,mean(posCoeffDistances{1}),std(posCoeffDistances{1}/sqrt(numTrials)),'lineprops',{'y-o','LineWidth',3})
hold on
shadedErrorBar(time,mean(posCoeffDistances{2}),std(posCoeffDistances{2}/sqrt(numTrials)),'lineprops',{'b-o','LineWidth',3})

xlabel('Time (s)')
ylabel(['Distance (' method ')'])
xline(cueOn,'k--','Approach','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

for i = 1:length(posCoeffDistances{1})
    [p(i),h(i)] = ranksum(posCoeffDistances{1}(:,i),posCoeffDistances{2}(:,i));
end

h_corrected = fdr_bh(p);
h_corrected(h_corrected == 0) = nan; 
h(h==0) = nan;
plot(time,h+0.2,'k*','MarkerSize',10)

saveas(gca,[figPath filesep 'posCoeff_corincor_shadedDistance_numTrials' num2str(numTrials) num2str(pcSort) method '_100ms'],'svg');     saveas(gca,[figPath filesep 'posCoeff_corincor_shadedDistance_numTrials' num2str(numTrials) num2str(pcSort) method '_100ms'],'png')

%% 05/25/2022
% Because I am missing the bootISI function on the new Linux, I may focus
% on the LFP Gamma Power to confirm that the approach timepoint analysis is
% broadly correct. This is a sanity check for future analyses using the
% approach timepoint. There is also some cleaning I can do of the script
% and perhaps creating functions to start cleaning up this giant script. 

%% Approach PSTH Creation for LFP Data
% Generally this should be similar to the spike data without the complexity
% of multiple neurons. Therefore it should simply be a TrialxSxT (Trial by Signal by Time)
% matrix. Signal is 1D, therefore we can collapse on it making a 2D matrix.
% 
% The trickiest part will be getting the timestamps to match.
%
% Pull ephys and trial data, create PSTH
    pBin = 0.1; % 33 ms bins for fr variable, match video framerate
    PSTH_Approach = [];
    fs = 1000;

    d = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
               'DesignMethod','butter','SampleRate',fs);


    if isfield(ephysStruct,'LFP_PSTH')
        
        for i = 1:length(trlStruct)
    
            % Load trial times from trlStruct
            trialTimes = trlStruct(i).trialTimes(1:48);
            [trialTimes, trialTimesIdx] = sort(trialTimes);
            approachVar = min(trlStruct(i).approach(:,1:2,2),[],2);
            sortedApproachVar = approachVar(trialTimesIdx);

            if startsWith(masterTbl.SessionType(i),'Regular')            
                correctIndex = trlStruct(i).approach(1:48,1,1) == 1 & trlStruct(i).approach(1:48,2,1) == 0;
                incorrectIndex = trlStruct(i).approach(1:48,2,1) == 1;
            elseif startsWith(masterTbl.SessionType(i),'Reversal')
                correctIndex = trlStruct(i).approach(1:48,1,1) == 0 & trlStruct(i).approach(1:48,2,1) == 1;
                incorrectIndex = trlStruct(i).approach(1:48,1,1) == 1;
            end

            correctIndex = correctIndex(trialTimesIdx);
            incorrectIndex = incorrectIndex(trialTimesIdx);

            % Load ephys related data
            LFP = ephysStruct(i).LFP;
            % LFP data has been previously downsampled to 1000Hz; allows
            % for analysis of data @ 500 Hz and below.

            % Find time, in seconds, that each fr index corresponds to.
            frTime = (1:length(LFP))/1000;   % Denominator corresponds with pBin
    
            dlcTimestamps = trlStruct(i).bodyCoords{1};
            dlcTimestamps = dlcTimestamps(dlcTimestamps >= 0);
            dlcTimestamps = round(dlcTimestamps,2);
            trialTimes = trialTimes + sortedApproachVar;
            trialTimes = round(trialTimes,2);
   %     
        for k = 1:length(trialTimes)
            [~,tMin] = min(abs((trialTimes(k) - 1) - dlcTimestamps));
            [~,tMax] = min(abs((trialTimes(k) + 1) - dlcTimestamps)); 
            [~,tMin_fr] = min(abs((dlcTimestamps(tMin)) - frTime));
            [~,tMax_fr] = min(abs((dlcTimestamps(tMax)) - frTime));            
            trlSize(k) = size(LFP(tMin_fr:tMax_fr,:),1);

            if size(LFP(tMin_fr:tMax_fr,:),1) == 2001
                LFP_PSTH(k,:,:) = LFP(tMin_fr:tMax_fr,:);
            elseif size(LFP(tMin_fr:tMax_fr,:),1) == 2002
                LFP_PSTH(k,:,:) = LFP(tMin_fr:tMax_fr-1,:);
            elseif size(LFP(tMin_fr:tMax_fr,:),1) == 2000
                LFP_PSTH(k,:,:) = LFP(tMin_fr:tMax_fr+1,:);          
            end

        end

        
        ephysStruct(i).LFP_PSTH = LFP_PSTH;
        ephysStruct(i).LFP_Correct = LFP_PSTH(correctIndex,:);
        ephysStruct(i).LFP_Incorrect = LFP_PSTH(incorrectIndex,:);
        
        ephysStruct(i).LFP_CorrectM = filtfilt(d,lowpass(mean(ephysStruct(i).LFP_Correct),100,fs));
        ephysStruct(i).LFP_IncorrectM = filtfilt(d,lowpass(mean(ephysStruct(i).LFP_Incorrect),100,fs));
        
        [s,w,t,f,ps] = spectrogram(ephysStruct(i).LFP_CorrectM,100,0,[25:100],fs,'psd'); 
        ephysStruct(i).Correct_SpectrogramRes = [{s} {w} {t} {f} {ps}];
        [s,w,t,f,ps] = spectrogram(ephysStruct(i).LFP_IncorrectM,100,0,[25:100],fs,'psd');
        ephysStruct(i).Incorrect_SpectrogramRes = [{s} {w} {t} {f} {ps}];

    end

    end

%% Given LFPs, create 2x2x2 analyses 
% Genotype x Session Type x Correct/Incorrect
% Take incorrect/correct trial mean per session, generate spectogram data per session,
% average spectrogram data, plot!
%% 
    PConCor = [ephysStruct(regPIdx).Correct_SpectrogramRes];
    PConCor = [PConCor(1:5:end)];
    PConCor = cat(3,PConCor{:});
    PIncCor = [ephysStruct(revPIdx).Correct_SpectrogramRes];
    PIncCor = [PIncCor(1:5:end)];
    PIncCor = cat(3,PIncCor{:});
    
    PConInc = [ephysStruct(regPIdx).Incorrect_SpectrogramRes];
    PConInc = [PConInc(1:5:end)];
    PConInc = cat(3,PConInc{:});
    PIncInc = [ephysStruct(revPIdx).Incorrect_SpectrogramRes];
    PIncInc = [PIncInc(1:5:end)];
    PIncInc = cat(3,PIncInc{:});
    
    
    WConCor = [ephysStruct(regWIdx).Correct_SpectrogramRes];
    WConCor = [WConCor(1:5:end)];
    WConCor = cat(3,WConCor{:});
    WIncCor = [ephysStruct(revWIdx).Correct_SpectrogramRes];
    WIncCor = [WIncCor(1:5:end)];
    WIncCor = cat(3,WIncCor{:});
    
    WConInc = [ephysStruct(regWIdx).Incorrect_SpectrogramRes];
    WConInc = [WConInc(1:5:end)];
    WConInc = cat(3,WConInc{:});
    WIncInc = [ephysStruct(revWIdx).Incorrect_SpectrogramRes];
    WIncInc = [WIncInc(1:5:end)];
    WIncInc = cat(3,WIncInc{:});
%% Find average per each condition, create 4x4 plot for each genotype
% Then create a difference/contrast plot with both genotypes on it
    mPConCor = mean(PConCor,3);
    mPIncCor = mean(PIncCor,3);
    mPConInc = mean(PConInc,3);
    mPIncInc = mean(PIncInc,3);

    mWConCor = mean(WConCor,3);
    mWIncCor = mean(WIncCor,3);
    mWConInc = mean(WConInc,3);
    mWIncInc = mean(WIncInc,3);
allMeans = [{mPConCor} {mPIncCor} {mPConInc} {mPIncInc} {mWConCor} {mWIncCor} {mWConInc} {mWIncInc}];
titles = [{'Congruent-Correct'} {'Incongruent-Correct'} {'Congruent-Incorrect'} {'Incongruent-Incorrect'} {'Congruent-Correct'} {'Incongruent-Correct'} {'Congruent-Incorrect'} {'Incongruent-Incorrect'}];
subMeans = [{abs(mPConInc) - abs(mPConCor)} {abs(mPIncInc) - abs(mPIncCor)} {abs(mWConInc) - abs(mWConCor)} {abs(mWIncInc) - abs(mWIncCor)}];
titles_sub = [ {'P Congruent'}, {'P Incongruent'}, {'W Congruent'}, {'W Incongruent'} ];

%% Create plots
    figure('Units','normalized','Position',[0 0 1 1]); sgtitle('P rats')
    for i = 1:length(allMeans)
        if i == 5; figure('Units','normalized','Position',[0 0 1 1]); sgtitle('Wistars'); end
        if i < 5
            subplot(2,2,i)
        elseif i >= 5
            subplot(2,2,i-4)
        end
    
        surf(t,(w),log10(abs(allMeans{i})),'EdgeColor','none')
        axis xy; axis tight; colormap(jet); view(0,90);
        xlabel('Time')
        ylabel('Frequency (Hz)')
        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
%         xline(1000,'k--','Approach','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        title(titles{i})
        colorbar 
    end

%% Subtracted plot
    figure('Units','normalized','Position',[0 0 1 1]); sgtitle('P rats')
    for i = 1:length(subMeans)
        subplot(2,2,i)        
        surf(t,(w),(zscore(subMeans{i})),'EdgeColor','none')
        axis xy; axis tight; colormap(jet); view(0,90);
        xlabel('Time')
        ylabel('Frequency (Hz)')
        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
%         xline(1000,'k--','Approach','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
        title(titles_sub{i})
        colorbar 
    end
