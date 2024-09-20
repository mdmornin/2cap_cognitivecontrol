%%
%% Load and check data
clear all; close all % Refresh workspace
% PATHS
% For brain3 Linux Machine:
if ispc
    parentPath = 'F:/dissDat';
    figPath = 'F:/dissDat/figs_sepPW';
else
    parentPath = '/research/dissDat';
    figPath = '/research/dissDat/figs';
end
% LOAD
load([parentPath filesep 'ephysStruct.mat']);
load([parentPath filesep 'trlStruct.mat']);
load([parentPath filesep 'masterTable.mat']);

addpath(genpath([parentPath filesep 'analysisScripts']))
addpath(genpath([parentPath filesep 'restoredScripts']))
%% Set indexing based on Congruent versus Incongruent Sessions
regPIdx = startsWith(masterTbl.SessionType,'Regular') & startsWith(masterTbl.Strain,'P');
revPIdx = startsWith(masterTbl.SessionType,'Reversal') & startsWith(masterTbl.Strain,'P');
    
regWIdx = startsWith(masterTbl.SessionType,'Regular') & startsWith(masterTbl.Strain,'W');
revWIdx = startsWith(masterTbl.SessionType,'Reversal') & startsWith(masterTbl.Strain,'W');

%% Timepoints and other parameters
numTrials = 48;


sipDescent = 6;
sipAscent = 14;
cueOn = 2;

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
numTrials = 15;
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
        
            ephysStruct(i).PSTH_redirectLeft = PSTH(correctionIndexL,:,:);
            ephysStruct(i).PSTH_redirectRight = PSTH(correctionIndexR,:,:);

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

        ephysStruct(i).PSTH_redirectLeft_mean = squeeze(nanmean(ephysStruct(i).PSTH_redirectLeft,1));
        ephysStruct(i).PSTH_redirectRight_mean = squeeze(nanmean(ephysStruct(i).PSTH_redirectRight,1));

        ephysStruct(i).correctionsIndex = correctionIndex;


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
ephysStruct(19).PSTH_ApproachSub = [];

ephysStruct(19).PSTH_CorrectMean = [];
ephysStruct(19).PSTH_IncorrectMean = [];
ephysStruct(19).PSTH_Approach = [];
ephysStruct(19).PSTH = [];
ephysStruct(19).PSTH_LeftMean = [];
ephysStruct(19).PSTH_RightMean = [];


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

pMat = cat(3,ephysStruct(PIdx).PSTH);
wMat = cat(3,ephysStruct(WIdx).PSTH);


% pMat = [ephysStruct(PIdx).PSTH_ApproachMean];
% wMat = [ephysStruct(WIdx).PSTH_ApproachMean];

meanType = 'All';

%% PCA Preprocessing | Separation between P, W groups
timepoint = 'all';

% WConditionsPre = squeeze(mean(wMat,1));
% PConditionsPre = squeeze(mean(pMat,1));

WConditionsPre = squeeze(mean(wMat(1:15,:,:),1));
PConditionsPre = squeeze(mean(pMat(1:15,:,:),1));

% WConditionsPre = wMat;
% PConditionsPre = pMat;

time = (1:size(WConditionsPre,1))/10;

[wConditionsCoeff,wConditionsScore,wConditionsLatent,~,wConditionsExplained] = pca(WConditionsPre);
[pConditionsCoeff,pConditionsScore,pConditionsLatent,~,pConditionsExplained] = pca(PConditionsPre);

% Need to find an index of neurons per dataset, be able to sort that index
% into our conditions (W/P; Con/Inc). Should ideally be a 1x54 matrix
% containing N # of neurons.

PneuronIndex = [ephysStruct(PIdx).neuronNumApproach];
WneuronIndex = [ephysStruct(WIdx).neuronNumApproach];

% PneuronIndex = [ephysStruct(PIdx).neuronNum];
% WneuronIndex = [ephysStruct(WIdx).neuronNum];
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
figure('Units','normalized','Position',[0 0 1 1])
subplot(2,2,1:2)
semWPre = std(WConditionsPre,[],2)/sqrt(length(WConditionsPre));
semPPre = std(PConditionsPre,[],2)/sqrt(length(PConditionsPre));

shadedErrorBar(time,mean(WConditionsPre,2),semWPre,'lineprops',{'b','LineWidth',3}); hold on; shadedErrorBar(time,mean(PConditionsPre,2),semPPre,'lineprops',{'r','LineWidth',3})
xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
if strcmp(timepoint,'all')
end
xlabel('Time (s)')
ylabel('Firing Rate (Hz)')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',25,'FontWeight','bold','LineWidth',5)
for j = 1:2

    if j == 1
        Explained = wConditionsExplained;
        titlevar = ['Explained Variance for Wistars'];
    elseif j == 2
        Explained = pConditionsExplained;
        titlevar = ['Explained Variance for P rats'];
    end

    p = length(Explained);
    pVec = 1:p;
    pExpected = zeros(length(pVec),1);
    
    for i = pVec
        pExpected(i) = sum( 1 ./ pVec(i:end)) / p;
    end
    pExpected = pExpected .* 100; % Convert to percentage
    stickPoint(j) = find(pExpected >= Explained,1);
    subplot(2,2,j+2)
    plot((Explained),'ko','MarkerSize',5);
    hold on
    plot(pExpected,'r--','LineWidth',2)
    xline(stickPoint(j),'k--','Broken Stick Point','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
    xlabel('PC Number')
    ylabel('Explained Variance')
    ylim([0 40])
    xlim([0 30])
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',25,'FontWeight','bold','LineWidth',5)
end
saveas(gca,[figPath filesep 'approach_varianceExplained_meanfiringrate'     '_100ms' ],'svg');     saveas(gca,[figPath filesep 'approach_varianceExplained_meanfiringrate'    '_100ms' ],'png')
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

         
    end
    sgtitle(titles{j})
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figPath filesep 'allNeurons_pcProj_prePCAsplitGENO_COR' titles{j}  num2str(pcSort)  '_100ms'],'svg');     saveas(gca,[figPath filesep 'allNeurons_pcProj_prePCAsplitGENO_COR' titles{j}  num2str(pcSort)  '_100ms'],'png')
end
%%
% colors = [{'#EDB120'}, {'#7E2F8E'}, {'#0072BD'},  {'#D95319'} ];
% pcSort = 1:5;
% % Plot coefficient loadings
% figure('Units','normalized','Position',[0 0 1 1])
% 
%     for i = 1:length(pcSort)
%         subplot(1,max(pcSort),i);
%         histogram(indexedCoeffData{1}(:,i),'Normalization','probability','DisplayStyle','stairs','LineWidth',2,'LineStyle','-','EdgeColor',colors{1})
%         hold on
%         histogram(indexedCoeffData{2}(:,i),'Normalization','probability','DisplayStyle','stairs','LineWidth',2,'LineStyle','-','EdgeColor',colors{2})
%         histogram(indexedCoeffData{3}(:,i),'Normalization','probability','DisplayStyle','stairs','LineWidth',2,'LineStyle','-','EdgeColor',colors{3})
%         histogram(indexedCoeffData{4}(:,i),'Normalization','probability','DisplayStyle','stairs','LineWidth',2,'LineStyle','-','EdgeColor',colors{4})
% 
%         xlabel(['PC ' num2str(i) ' Coefficients'])
%         ylabel('Probability')
%         xlim([-.2,.2])
%         ylim([0 0.45])
% 
%         legend([{'PCon'},{'PInc'},{'WCon'},{'WInc'}])
%         set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
% 
% 
%       [h1(i),p1(i)] = kstest2(indexedCoeffData{1}(:,i),indexedCoeffData{2}(:,i));
%       [h2(i),p2(i)] = kstest2(indexedCoeffData{1}(:,i),indexedCoeffData{3}(:,i));
%       [h3(i),p3(i)] = kstest2(indexedCoeffData{1}(:,i),indexedCoeffData{4}(:,i));
%       [h4(i),p4(i)] = kstest2(indexedCoeffData{2}(:,i),indexedCoeffData{3}(:,i));
%       [h5(i),p5(i)] = kstest2(indexedCoeffData{2}(:,i),indexedCoeffData{4}(:,i));
%       [h6(i),p6(i)] = kstest2(indexedCoeffData{3}(:,i),indexedCoeffData{4}(:,i));
% 
% 
%     end
%     saveas(gca,[figPath filesep 'coefficienthistogram_split_scores_prePCAsplitGENO_COR' num2str(pcSort)  '_100ms'],'svg');     saveas(gca,[figPath filesep 'coefficienthistogram_split_scores_prePCAsplitGENO_COR' num2str(pcSort)  '_100ms'],'png')
% 
% %% Plot outputs of KSTEST2
analysisType = 'app';
% allP = [p1; p2; p3; p4; p5; p6];
% allH = [h1; h2; h3; h4; h5; h6];
% ylabels = [{'PCon-PInc'}; {'PCon-WCon'}; {'PCon-WInc'}; {'PInc-WCon'}; {'PInc-WInc'}; {'WCon-WInc'}];
% xlabels = [{'PC1'}; {'PC2'}; {'PC3'}; {'PC4'}; {'PC5'}];
% 
% figure('Units','normalized','Position',[0 0 1 1])
% imagesc(allH,[0 1]); colorbar; colormap jet;
% set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4,'XTick',1:6,'XTickLabel', xlabels, 'YTick',1:6,'YTickLabel',ylabels)
% saveas(gca,[figPath filesep 'coefficientKSTEST_HVALS_prePCAsplitGENO_COR' analysisType num2str(pcSort)  '_100ms'],'svg');     saveas(gca,[figPath filesep 'coefficientKSTEST_HVALS_prePCAsplitGENO_COR' analysisType num2str(pcSort)  '_100ms'],'png')
% 
% figure('Units','normalized','Position',[0 0 1 1])
% imagesc(allP,[0 1]); colorbar; colormap jet;
% set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4,'XTick',1:6,'XTickLabel', xlabels, 'YTick',1:6,'YTickLabel',ylabels)
% saveas(gca,[figPath filesep 'coefficientKSTEST_PVALS_prePCAsplitGENO_COR' analysisType num2str(pcSort)  '_100ms'],'svg');     saveas(gca,[figPath filesep 'coefficientKSTEST_PVALS_prePCAsplitGENO_COR' analysisType num2str(pcSort)  '_100ms'],'png')
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
        ylim([-10 10])
        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',10,'FontWeight','bold','LineWidth',4)
        xline(cueOn,'k--','Approach','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
    end
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
end
saveas(gca,[figPath filesep 'allNeurons_split_top5pcs_prePCAsplitGENO_COR' analysisType num2str(pcSort)  '_100ms'],'svg');     saveas(gca,[figPath filesep 'allNeurons_split_top5pcs_prePCAsplitGENO_COR' analysisType num2str(pcSort)  '_100ms'],'png')


%% Plot resulting data in a 2D Space
dim = 2;
figure('Units','normalized','Position',[0 0 1 1])
colors = [{'#EDB120'}, {'#7E2F8E'}, {'#0072BD'},  {'#D95319'} ];

for i = 1:length(allPCProjections)
    plot3(allPCProjections{i}(1:end,1),allPCProjections{i}(1:end,2),allPCProjections{i}(1:end,3),'Color',colors{i},'LineWidth',3);
    hold on
    plot3(allPCProjections{i}(1,1),allPCProjections{i}(1,2),allPCProjections{i}(1,3),'o','MarkerSize',10,'MarkerFaceColor',colors{i},'MarkerEdgeColor','k')
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)            
end

saveas(gca,[figPath filesep 'allNeurons_split_pcSpaceProj_prePCAsplitGENO_COR 2' analysisType   '_100ms'],'svg');     saveas(gca,[figPath filesep 'allNeurons_split_pcSpaceProj_prePCAsplitGENO_COR 2'   analysisType  '_100ms'],'png')

%% Calculate and Plot the Distance Between PCs 
pCon3D = allPCProjections{1}(:,1:stickPoint(2)); wCon3D = allPCProjections{3}(:,1:stickPoint(1));
pInc3D = allPCProjections{2}(:,1:stickPoint(2)); wInc3D = allPCProjections{4}(:,1:stickPoint(1));

pDistance = diag(pdist2(pCon3D,pInc3D,'euclidean'));
wDistance = diag(pdist2(wCon3D,wInc3D,'euclidean'));

figure('Units','normalized','Position',[0 0 1 1])
plot(time,pDistance,'Color',colors{1},'LineWidth',3)
hold on
plot(time,wDistance,'Color',colors{3},'LineWidth',3)
xline(cueOn,'k--','Approach','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
if strcmp(timepoint,'all')
end
xlabel('Time (s)')
ylabel('Distance (euclidean)')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,[figPath filesep 'allNeurons_split_3DpcSpaceDistance_prePCAsplitGENO_COR' timepoint meanType '_100ms'],'svg');     saveas(gca,[figPath filesep 'allNeurons_split_3DpcSpaceDistance_prePCAsplitGENO_COR' timepoint meanType '_100ms'],'png')

%


%%
%% 

pCon = cat(3,ephysStruct(regPIdx).PSTH);
pInc = cat(3,ephysStruct(revPIdx).PSTH);
wCon = cat(3,ephysStruct(regWIdx).PSTH);
wInc = cat(3,ephysStruct(revWIdx).PSTH);



%%

% psthType = 'cueT';
% psthType = 'cueT';

numTrials = 15;
for j = 1:length(indexedCoeffData)
    pcSort = 1:5;
    % Assign variable rawData depending on condition 
    % Center rawData variable to then multiply by coeff to generate per
    % trial score.
    if j == 1
        rawData = pCon - mean(pCon,3);
        pcSort = 1:stickPoint(1);
    elseif j == 2
        rawData = pInc - mean(pInc,3);
        pcSort = 1:stickPoint(1);
    elseif j == 3
        rawData = wCon - mean(wCon,3);
        pcSort = 1:stickPoint(2);
    elseif j == 4
        rawData = wInc - mean(wInc,3);
        pcSort = 1:stickPoint(2);
    end
    % Sort according to top 3 PCs, plot similar to other coefficient
    % figures

    figure('Units','normalized','Position',[0 0 1 1])

    for i = 1:length(pcSort)
        % Index Pos and Neg
        PosCoeff = [];
        NegCoeff = [];
        for k = 1:numTrials
            PosCoeff(k,:,:) = indexedCoeffData{j}(indexedCoeffData{j}(:,i) >= 0.01,i)' .* squeeze(rawData(k,:,indexedCoeffData{j}(:,i) >= 0.01));                  
            NegCoeff(k,:,:) = indexedCoeffData{j}(indexedCoeffData{j}(:,i) <= -0.01,i)' .* squeeze(rawData(k,:,indexedCoeffData{j}(:,i) <= -0.01));
        end
        % Take Mean
        
        PosCoeff = (mean(PosCoeff,3));
        NegCoeff = (mean(NegCoeff,3));
        % Index extra mean information
        mPosCoeff{j,i} = mean(PosCoeff,1);
        mNegCoeff{j,i} = mean(NegCoeff,1);
        stdPosCoeff{j,i} = std(PosCoeff,1)/sqrt(15);
        stdNegCoeff{j,i} = std(PosCoeff,1)/sqrt(15);
        % Index Pos Neg Info
        indexedPosCoeff{j,i} = PosCoeff;
        indexedNegCoeff{j,i} = NegCoeff;


        % Figure
%         subplot(2,max(pcSort),i);          % Doubling subplot for pos+neg
%         % Positive Plot
%         imagesc([0 max(time)], [1 15], zscore(PosCoeff), [-3 3]);
%         colormap jet
%         colorbar
%         xlabel('Time (s)');
%         if i==1; ylabel('Trial #');end
%         xline(cueOn,'k--','Approach','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
%         xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
%         xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
% %         
%         if j < 3; title(['PC Num: ' num2str(pcSort(i)) ',ExplVar=' num2str(sum(pConditionsExplained(1:i)),'%.1f') '%']);
%         else
%             title(['PC Num: ' num2str(pcSort(i)) ',ExplVar=' num2str(sum(wConditionsExplained(1:i)),'%.1f') '%']);
%         end
%         set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
%         % Negative Plot
% 
%         subplot(2,max(pcSort),i+length(pcSort));          % Doubling subplot for pos+neg
% 
%         imagesc([0 max(time)], [1 15], zscore(NegCoeff), [-3 3]);
%         colormap jet
%         colorbar
%         xlabel('Time (s)');
%         if i==1; ylabel('Trial #');end
%         xline(cueOn,'k--','Approach','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
%         xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
%         xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
% %         
%         if j < 3; title(['PC Num: ' num2str(pcSort(i)) ',ExplVar=' num2str(sum(pConditionsExplained(1:i)),'%.1f') '%']);
%         else
%             title(['PC Num: ' num2str(pcSort(i)) ',ExplVar=' num2str(sum(wConditionsExplained(1:i)),'%.1f') '%']);
%         end
%         set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)        
    end
% 
%     sgtitle(titles{j})
%     set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
%     saveas(gca,[figPath filesep 'coeffNeurons_overTrials' titles{j} '_100ms'],'svg');     saveas(gca,[figPath filesep 'coeffNeurons_overTrials' titles{j}  '_100ms'],'png')

end


%% Test Trial-by-Trial PC Loads
addpath(genpath('/research/dissDat/restoredScripts'))
%% Calculate max val per epoch per trial similar to behavior
% Time periods to index epochs with
psthType = 'approachT';
pDist = [];
wDist = [];
if strcmp(psthType,'cueT')
    cOn = cueOn * 10;
    sipD = sipDescent * 10;
    sipA = sipAscent * 10;
    midPoint = (sipDescent + sipAscent) / 2 * 10;
    
    % Check that epoch variables align with time variable
    preCue = time(1:cOn);
    CueOn = time(cOn+1:sipD);
    EarlyEtOH = time(sipD+1:midPoint);
    LateEtOH = time(midPoint+1:sipA);
    PostEtOH = time(sipA+1:end);
    plot(time,[ones(size(preCue))* 1 ones(size(CueOn)) * 2 ones(size(EarlyEtOH)) * 3 ones(size(LateEtOH)) * 4 ones(size(PostEtOH)) * 5],'k');
    ylim([0 6])
    
    % Reassign epoch variables so that they are proper indices rather than
    % time values
    preCue = (1:cOn);
    CueOn = (cOn+1:sipD);
    EarlyEtOH = (sipD+1:midPoint);
    LateEtOH = (midPoint+1:sipA);
    PostEtOH = (sipA+1:length(time));
elseif strcmp(psthType,'approachT')
    cOn = cueOn * 10;
    preApproach = 1:cOn;
    postApproach = cOn+1:length(time);
end

%% Distance calculations for positive coefficients 
coeffValue = 'negCoeff';

if strcmp(coeffValue,'posCoeff')
    pCon3D = cat(3,indexedPosCoeff{1,:}); wCon3D = cat(3,indexedPosCoeff{3,:});
    pInc3D = cat(3,indexedPosCoeff{2,:}); wInc3D = cat(3,indexedPosCoeff{4,:});
elseif strcmp(coeffValue,'negCoeff')
    pCon3D = cat(3,indexedNegCoeff{1,:}); wCon3D = cat(3,indexedNegCoeff{3,:});
    pInc3D = cat(3,indexedNegCoeff{2,:}); wInc3D = cat(3,indexedNegCoeff{4,:});
end

method = 'mahalanobis';

for k = 1:numTrials
    pDist(k,:) = diag(pdist2(squeeze(pCon3D(k,:,:)),squeeze(pInc3D(k,:,:)),method)); 
    wDist(k,:) = diag(pdist2(squeeze(wCon3D(k,:,:)),squeeze(wInc3D(k,:,:)),method)); 
end

% Plot according to epochs (psthType == cueT)
figure('Units','normalized','Position',[0 0 1 1])


numEpochs = 5;
epochLabels = {'PreCue','CueOn','EarlyEtOH','LateEtOH','PostEtOH'};

if strcmp(psthType,'cueT')

    numEpochs = 5;
    epochLabels = {'PreCue','CueOn','EarlyEtOH','LateEtOH','PostEtOH'};

    for j = 1:numEpochs
            if j == 1
                tp = preCue;
            elseif j == 2
                tp = CueOn;
            elseif j == 3
                tp = EarlyEtOH;
            elseif j == 4
                tp = LateEtOH;
            elseif j == 5
                tp = PostEtOH;
            end
        
            pdistMaxMean(j) = mean(max(pDist(:,tp),[],2));
            wdistMaxMean(j) = mean(max(wDist(:,tp),[],2));

            pdm(:,j) = max(pDist(:,tp),[],2);
            wdm(:,j) = max(wDist(:,tp),[],2);
        
            pdistSTD(j) = std(max(pDist(:,tp),[],2))/sqrt(numTrials);
            wdistSTD(j) = std(max(wDist(:,tp),[],2))/sqrt(numTrials);            
    end

elseif strcmp(psthType,'approachT')

    numEpochs = 2;
    epochLabels = {'PreApproach','PostApproach'};
    for j = 1:numEpochs
        if j == 1
            tp = preApproach;
        elseif j == 2
            tp = postApproach;
        end
    
        pdistMaxMean(j) = mean(max(pDist(:,tp),[],2));
        wdistMaxMean(j) = mean(max(wDist(:,tp),[],2));

        pdm(:,j) = max(pDist(:,tp),[],2);
        wdm(:,j) = max(wDist(:,tp),[],2);
    
        pdistSTD(j) = std(max(pDist(:,tp),[],2))/sqrt(numTrials);
        wdistSTD(j) = std(max(wDist(:,tp),[],2))/sqrt(numTrials); 
    end
end

plot(pdistMaxMean,'k-o','LineWidth',3)
hold on
plot(wdistMaxMean,'b-o','LineWidth',3)
xticks(1:1:length(epochLabels))
xticklabels(epochLabels)
ylabel(['Distance (' method ')'])
xlim([0.9 numEpochs+0.1])
er = errorbar(pdistMaxMean, pdistSTD);
er.Color = [0 0 0];
er.LineStyle = 'none';
er.LineWidth = 3;

er = errorbar(wdistMaxMean, wdistSTD);
er.Color = [0 0 0];
er.LineStyle = 'none';
er.LineWidth = 3;
% ylim([0 10])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
saveas(gca,[figPath filesep 'Coeff_Distance_GenotypesComparison' coeffValue psthType '_' method timepoint '_100ms'],'svg');     saveas(gca,[figPath filesep 'Coeff_Distance_GenotypesComparison' psthType '_' method coeffValue timepoint '_100ms'],'png')

% Stats
trialRep = repmat([1:15],1,2); trialRep = trialRep';
tableData = [wdm; pdm];
maxTable = array2table(tableData,'VariableNames',{'PreApproach','PostApproach'});
maxTable.Strain = [repmat({'Wistar'},15,1); repmat({'P rat'},15,1)];
maxTable.Trial = trialRep;

% Create rm model
Epoch = table([1 2]','VariableNames',{'Epochs'});
rm = fitrm(maxTable,'PreApproach-PostApproach ~ Strain','WithinDesign',Epoch,'WithinModel','Epochs');
anovaTbl = ranova(rm,'WithinModel','Epochs')
multTbl1 = multcompare(rm,'Strain')
multTbl2 = multcompare(rm,'Strain','By','Epochs')

tblStr = formattedDisplayText(anovaTbl); 
% Write string to file
fid = fopen('anovaTable_maxDistance.txt', 'wt');
fileCleanup = onCleanup(@()fclose(fid));
formatSpec = '%s\n';
fprintf(fid, formatSpec, tblStr);
clear('fileCleanup')

tblStr = formattedDisplayText(multTbl2); 
% Write string to file
fid = fopen([coeffValue 'anovaTable_maxDistance_multcomp.txt'], 'wt');
fileCleanup = onCleanup(@()fclose(fid));
formatSpec = '%s\n';
fprintf(fid, formatSpec, tblStr);
clear('fileCleanup')

% Calculate distance per trial
% Consider each trial trace a PC, calculate distance between Trial N of
% Congruent and Trial N of Incongruent per each condition. 
%
% Given variables pDist and WDist -- create a shadedErrorBar graph that
% can be used on the poster.
numTrials = size(pDist,1);

mpDist = mean(pDist); sempDist = std(pDist)/sqrt(numTrials);
mwDist = mean(wDist); semwDist = std(wDist)/sqrt(numTrials);

figure('Units','normalized','Position',[0 0 1 1])
shadedErrorBar(time,mpDist,sempDist,'lineprops',{'r','LineWidth',3})
hold on
shadedErrorBar(time,mwDist,semwDist,'lineprops',{'b','LineWidth',3})    
xlabel('Time (s)')
ylabel(['Distance (' method ')'])
ylim([1 7])

if strcmp(psthType,'cueT')
    xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');    
    xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
    xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
elseif strcmp(psthType,'approachT')
    xline(cueOn,'k--','Approach','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');    
end

set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,[figPath filesep 'PostrlAvgDistances' psthType timepoint coeffValue '_100ms'],'svg');     saveas(gca,[figPath filesep 'PostrlAvgDistances' psthType timepoint coeffValue '_100ms'],'png')

% Omnibus stat test
totalT = size(pDist,2);
statTable = [wDist; pDist];
genotypeLabel = [repmat({'Wistar'},15,1); repmat({'P rat'},15,1)];

statTable = array2table(statTable);
statTable.Genotype = genotypeLabel;
Time = table([1:totalT]','VariableNames',{'Time'});

rm = fitrm(statTable,'statTable1-statTable41 ~ Genotype','WithinDesign',Time,'WithinModel','Time');
anovaTbl = ranova(rm,'WithinModel','Time')


tblStr = formattedDisplayText(anovaTbl); 
fid = fopen([ coeffValue 'appDistance_RANOVA.txt'], 'wt');
fileCleanup = onCleanup(@()fclose(fid));
formatSpec = '%s\n';
fprintf(fid, formatSpec, tblStr);
clear('fileCleanup')

% Run stats on distance per trial
figure('Units','normalized','Position',[0 0.5 1 .3])
for i = 1:totalT
    [h(i),p(i),stat{i}] = ttest2(wDist(:,i),pDist(:,i));
end
fdr_p = fdr_bh(p);
fdr_p_plot = double(fdr_p);
fdr_p_plot(fdr_p_plot == 0) = NaN;

% Combine
figure('Units','normalized','Position',[0 0 1 1])
shadedErrorBar(time,mpDist,sempDist,'lineprops',{'r','LineWidth',3})
hold on
shadedErrorBar(time,mwDist,semwDist,'lineprops',{'b','LineWidth',3})    
xlabel('Time (s)')
ylabel(['Distance (' method ')'])
ylim([1 7])
xlim([0 4.1])
if strcmp(psthType,'cueT')
    xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');    
    xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
    xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
elseif strcmp(psthType,'approachT')
    xline(cueOn,'k--','Approach','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');    
end
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

plot(time,fdr_p_plot+0.5, 'k*','MarkerSize',6)
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,[figPath filesep 'PostrlAvgDistances_withmarkers' psthType coeffValue timepoint '_100ms'],'svg');     saveas(gca,[figPath filesep 'PostrlAvgDistances_withmarkers' psthType coeffValue timepoint '_100ms'],'png')



% %% Reformatting data to contain stacks of L/R choices
% % L/R Choices will be stratified by genotype this time and stacked either
% % horiziontally or vertically depending on which explains the most
% % variance.
% 
% titles = [{'Congruent P'} {'Incongruent P'} {'Congruent W'} {'Incongruent W'}];
% 
% % A lot of this will repeat code from above. 
% % Set indexing based on Genotype
% pMat = [];
% wMat = [];
% %
% PIdx = startsWith(masterTbl.Strain,'P');
% WIdx = startsWith(masterTbl.Strain,'W');
% % Index Wistar, P
% pMat = [ephysStruct(PIdx).PSTH_LeftMean; ephysStruct(PIdx).PSTH_RightMean];
% wMat = [ephysStruct(WIdx).PSTH_LeftMean; ephysStruct(WIdx).PSTH_RightMean];
% 
% meanType = 'All';
% 
% %% PCA Preprocessing | Separation between P, W groups
% analysisType = 'stackedVertical';
% 
%     WConditionsPre = wMat;
%     PConditionsPre = pMat;
%     time = (1:size(WConditionsPre,1))/10;
% 
% 
% 
% [wConditionsCoeff,wConditionsScore,wConditionsLatent,~,wConditionsExplained] = pca(WConditionsPre);
% [pConditionsCoeff,pConditionsScore,pConditionsLatent,~,pConditionsExplained] = pca(PConditionsPre);
% 
% % Need to find an index of neurons per dataset, be able to sort that index
% % into our conditions (W/P; Con/Inc). Should ideally be a 1x54 matrix
% % containing N # of neurons.
% 
% PneuronIndex = [ephysStruct(PIdx).neuronNumApproach];
% WneuronIndex = [ephysStruct(WIdx).neuronNumApproach];
% 
% % PneuronIndex = [ephysStruct(PIdx).neuronNumApproach; ephysStruct(PIdx).neuronNumApproach];
% % WneuronIndex = [ephysStruct(WIdx).neuronNumApproach; ephysStruct(WIdx).neuronNumApproach];
% 
% %% Sanity check
% if sum(PneuronIndex) == size(pConditionsCoeff,1) && sum(WneuronIndex) == size(wConditionsCoeff,1)
%     disp('Matrix sizes match')
% else
%     disp('Matrix size mismatch, address')
% end
% 
% %% Assign ID to Neurons
% PneuronID = []; WneuronID = [];
% for i = 1:length(PneuronIndex)
%     PneuronID = [PneuronID ones(1,PneuronIndex(i))*i];
% end
% 
% for i = 1:length(WneuronIndex)
%     WneuronID = [WneuronID ones(1,WneuronIndex(i))*i];
% end
% 
% 
% %% Sanity check
% if length(PneuronID) == size(pConditionsCoeff,1) && length(WneuronID) == size(wConditionsCoeff,1)
%     disp('Matrix sizes match')
% else
%     disp('Matrix size mismatch, address')
% end
% 
% %% Pull data 
% % A lot of extra indexing to get to the point where the neurons are broken
% % up by session type within the genotype splits. This is due to the
% % structure being made to be read 1:54 rather than splitting it beforehand.
% % I think what I did here is adequate however and produces the required
% % result. 
% pSessionLabel = masterTbl.SessionType(PIdx);
% wSessionLabel = masterTbl.SessionType(WIdx);
% 
% pConSessionIndex = [startsWith(pSessionLabel,'Regular'); startsWith(pSessionLabel,'Regular')];
% wConSessionIndex = [startsWith(wSessionLabel,'Regular'); startsWith(wSessionLabel,'Regular')];
% 
% pIncSessionIndex = [startsWith(pSessionLabel,'Reversal'); startsWith(pSessionLabel,'Reversal')];
% wIncSessionIndex = [startsWith(wSessionLabel,'Reversal'); startsWith(wSessionLabel,'Reversal')];
% 
% %% Find raw, centered data
% conPIdx = find(pConSessionIndex == 1);
% incPIdx = find(pIncSessionIndex == 1); 
% conWIdx = find(wConSessionIndex == 1);
% incWIdx = find(wIncSessionIndex == 1);
% 
% 
% WcentData = wConditionsCoeff * wConditionsScore';
% PcentData = pConditionsCoeff * pConditionsScore';
% 
% % Index 
% conPCent = PcentData(any(PneuronID == find(pConSessionIndex == 1)),:);
% incPCent = PcentData(any(PneuronID == find(pIncSessionIndex == 1)),:);
% conWCent = WcentData(any(WneuronID == find(wConSessionIndex == 1)),:);
% incWCent = WcentData(any(WneuronID == find(wIncSessionIndex == 1)),:);
% % Coefficient
% conPCoeff = pConditionsCoeff(any(PneuronID == find(pConSessionIndex == 1)),:);
% incPCoeff = pConditionsCoeff(any(PneuronID == find(pIncSessionIndex == 1)),:);
% conWCoeff = wConditionsCoeff(any(WneuronID == find(wConSessionIndex == 1)),:);
% incWCoeff = wConditionsCoeff(any(WneuronID == find(wIncSessionIndex == 1)),:);
% % Create singular variable
% indexedCentData = {conPCent incPCent conWCent incWCent};
% indexedCoeffData = {conPCoeff incPCoeff conWCoeff incWCoeff};
% cellIndex = {conPIdx incPIdx conWIdx incWIdx};
% 
% %% Low-D PC Representation Based on Separated Matrices Stacked Horizontally
% % Parameters
% allPCProjections = [];
% pcSort = 1:5;
% % Obtain data
% for i = 1:length(cellIndex)
%     if i == 1 || i == 2
%         for k = 1:length(pcSort)
%             allPCProjections{i}(:,k) = PConditionsPre(:,any(PneuronID == cellIndex{i})) * indexedCoeffData{i}(:,pcSort(k));
%         end
%     elseif i == 3 || i == 4
%         for k = 1:length(pcSort)
%             allPCProjections{i}(:,k) = WConditionsPre(:,any(WneuronID == cellIndex{i})) * indexedCoeffData{i}(:,pcSort(k));
%         end
%     end
% end
% %%
% % %% Low-D PC Representation Based on Separated Matrices Stacked Vertically
% % % Parameters
% % allPCProjectionsL = [];
% % allPCProjectionsR = [];
% % 
% % pcSort = 1:5;
% % % Obtain data
% % for i = 1:length(cellIndex)
% % 
% %     if i == 1 || i == 2
% %             numN = size(indexedCoeffData{i},1);
% %             LIdx = 1:size(indexedCoeffData{i},1)/2;
% %             RIdx = LIdx(end)+1:size(indexedCoeffData{i},1);
% % 
% %         for k = 1:length(pcSort)
% % 
% %             allPCProjectionsL{i}(:,k) = PConditionsPre(:,any(PneuronID(1:numN) == cellIndex{i})) * indexedCoeffData{i}(LIdx,pcSort(k));
% %             allPCProjectionsR{i}(:,k) = PConditionsPre(:,any(PneuronID(numN+1:end) == cellIndex{i})) * indexedCoeffData{i}(RIdx,pcSort(k));
% % 
% %         end
% % 
% %     elseif i == 3 || i == 4
% % 
% %             LIdx = 1:size(indexedCoeffData{i},1)/2;
% %             RIdx = LIdx(end):size(indexedCoeffData{i},1);
% % 
% %         for k = 1:length(pcSort)
% % 
% %             allPCProjectionsL{i}(:,k) = WConditionsPre(:,any(WneuronID == cellIndex{i})) * indexedCoeffData{i}(:,pcSort(k));
% %             allPCProjectionsR{i}(:,k) = WConditionsPre(:,any(WneuronID == cellIndex{i})) * indexedCoeffData{i}(:,pcSort(k));
% % 
% %         end
% % 
% %     end
% % 
% % end
% %% Plot Top 5 Resulting Matrices for Each Condition
% pcSort = 1:5;
% colors = [{'#0072BD'} {'#4DBEEE'} {'#77AC30'} {'#D95319'} {'#A2142F'}];
% 
% figure('Units','normalized','Position',[0 0 1 1])
% 
% timePlot = -10:0.1:4.1;
% % Plotting Top 5 PCs As Normal, No Split Between L/R - just continuous 
% for j = 1:length(allPCProjections)
%     for i=1:length(pcSort)        
%         subplot(2,2,j);
%         plot(timePlot,allPCProjections{j}(:,i),'LineWidth',3,'Color',colors{i})
%         if j == 1; legend([{'PC1'},{'PC2'},{'PC3'},{'PC4'},{'PC5'}],'AutoUpdate','off','location','best'); end
%         hold on
%         xlabel('Time (s)');
%         title([titles{j}]);
%         set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',10,'FontWeight','bold','LineWidth',4)
%         xline(-5,'k--','ApproachL','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
%         xline(-3,'k--','LRSwap','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
%         xline(2,'k--','ApproachR','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
%     end
%     set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
% end
% saveas(gca,[figPath filesep 'allNeurons_split_top5pcs_prePCAsplitGENO_COR' analysisType num2str(pcSort)  '_100ms'],'svg');     saveas(gca,[figPath filesep 'allNeurons_split_top5pcs_prePCAsplitGENO_COR' analysisType num2str(pcSort)  '_100ms'],'png')
% 
% 
% %% Plotting TOp 5 PCs with a Split Between L/R (solid v dashed line)
% time = 1:size(wConditionsScore,1);
% timeL = time(1:71);
% timeR = time(72:end);
% timePlot = -5:0.1:2;
% colors = [{'#0072BD'} {'#4DBEEE'} {'#77AC30'} {'#D95319'} {'#A2142F'}];
% 
% figure('Units','normalized','Position',[0 0 1 1])
% for j = 1:length(allPCProjections)
%     for i=1:length(pcSort)        
%         subplot(2,2,j);
%         plot(timePlot,allPCProjections{j}(timeL,i),'LineWidth',3,'Color',colors{i})
%         legend([{'PC1'},{'PC2'},{'PC3'},{'PC4'},{'PC5'}],'AutoUpdate','off','location','best');
%         hold on
% 
%         plot(timePlot,allPCProjections{j}(timeR,i),'LineWidth',3,'Color',colors{i},'LineStyle','--')
%         xlabel('Time (s)');
%         title([titles{j}]);
%         set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',10,'FontWeight','bold','LineWidth',4)
%         xline(0,'k--','Approach','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold','HandleVisibility','off');
%     end
%     set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
% end
% saveas(gca,[figPath filesep 'allNeurons_split_LR' analysisType num2str(pcSort)  '_100ms'],'svg');     saveas(gca,[figPath filesep 'allNeurons_split_LR' analysisType num2str(pcSort)  '_100ms'],'png')
% %% Determine correlations between PCs -- see if anticorrelated PCs 
% % are meaningfully separating between L and R choices 


