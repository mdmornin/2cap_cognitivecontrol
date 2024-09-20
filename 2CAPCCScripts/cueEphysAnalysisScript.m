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
% 100ms. This should do several things: increasing stability of PCs, nanmean
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
    figPath = 'F:/dissDat/cueFigs';
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
            h=histc(stMtx,[min(stMtx(:)):pBin:max(stMtx(:))]);
            fr = gaussconv(stMtx,pBin);  % Gaussconv function is preferred
    
            % Find time, in seconds, that each fr index corresponds to.
            frTime = (1:length(fr))/10;   % Denominator corresponds with pBin
    
            dlcTimestamps = trlStruct(i).bodyCoords{1};
            dlcTimestamps = dlcTimestamps(dlcTimestamps >= 0);
            dlcTimestamps = (dlcTimestamps);
            trialTimes = (trialTimes);
        
        for k = 1:length(trialTimes(1:15))
            [~,tMin] = min(abs((trialTimes(k) - 4) - dlcTimestamps));
            [~,tMax] = min(abs((trialTimes(k) + 18) - dlcTimestamps)); 

            [~,tMin_fr] = min(abs((dlcTimestamps(tMin)) - frTime));
            [~,tMax_fr] = min(abs((dlcTimestamps(tMax)) - frTime));
            
            trlSize(k) = size(fr(tMin_fr:tMax_fr,:),1);
% 
            if size(fr(tMin_fr:tMax_fr,:),1) == 221
                PSTH(k,:,:) = fr(tMin_fr:tMax_fr,:);
            elseif size(fr(tMin_fr:tMax_fr,:),1) == 222
                PSTH(k,:,:) = fr(tMin_fr:tMax_fr-1,:);
            elseif size(fr(tMin_fr:tMax_fr,:),1) == 220
                PSTH(k,:,:) = fr(tMin_fr:tMax_fr+1,:);          
            end
            

               
        end

        

        % Remove low firing rate neurons. Not sure what the best value is
        % here.
        minFR = 0.01;
        if ~isempty(PSTH)            
            nanmeanFR = squeeze(nanmean(nanmean(PSTH,1),2)); 
            lowFR = nanmeanFR > minFR; % Find values greater than minFR
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

titles = [{'Congruent P'} {'Incongruent P'} {'Congruent W'} {'Incongruent W'}];

sipDescent = 9;
sipAscent = 17;
cueOn = 4;

LeftTrials = 1:24;
RightTrials = 25:48;
colors = [{'#EDB120'}, {'#7E2F8E'}, {'#0072BD'},  {'#D95319'} ];

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
    WConditionsPre = squeeze(nanmean(wMat(1:15,:,:),1));
    PConditionsPre = squeeze(nanmean(pMat(1:15,:,:),1));
    time = (1:size(WConditionsPre,1))/10;

elseif strcmp(timepoint,'sipDescent')
    WConditionsPre = squeeze(nanmean(wMat(1:15,1:sipDescent*10,:),1));
    PConditionsPre = squeeze(nanmean(pMat(1:15,1:sipDescent*10,:),1)); 
    time = (1:size(WConditionsPre,1))/10;
end


[wConditionsCoeff,wConditionsScore,wConditionsLatent,~,wConditionsExplained] = pca(WConditionsPre);
[pConditionsCoeff,pConditionsScore,pConditionsLatent,~,pConditionsExplained] = pca(PConditionsPre);

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
figure('Units','normalized','Position',[0 0 1 1])
subplot(2,2,1:2)
semWPre = std(WConditionsPre,[],2)/sqrt(length(WConditionsPre));
semPPre = std(PConditionsPre,[],2)/sqrt(length(PConditionsPre));

shadedErrorBar(time,nanmean(WConditionsPre,2),semWPre,'lineprops',{'b','LineWidth',3}); hold on; shadedErrorBar(time,nanmean(PConditionsPre,2),semPPre,'lineprops',{'r','LineWidth',3})
xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
if strcmp(timepoint,'all')
xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
end
xlabel('Time (s)')
ylabel('Firing Rate (Hz)')
axis tight
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
    ylim([0 20])
    xlim([0 30])
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',25,'FontWeight','bold','LineWidth',5)
end
saveas(gca,[figPath filesep 'varianceExplained_nanmeanfiringrate'     '_100ms' ],'svg');     saveas(gca,[figPath filesep 'varianceExplained_nanmeanfiringrate'    '_100ms' ],'png')

%% Reviewer Comment Analysis | Find number and strength of modulation of cue vs sipper access 
% For each neuron N in condition C, determine if Firing Rate (FR) is
% significantly different at time of cue or sipper 
% Use 2s of values before the cue and then -1s before to 1s after the cue.
pCon = squeeze(nanmean(cat(3,ephysStruct(regPIdx).PSTH)));
pInc = squeeze(nanmean(cat(3,ephysStruct(revPIdx).PSTH)));
wCon = squeeze(nanmean(cat(3,ephysStruct(regWIdx).PSTH)));
wInc = squeeze(nanmean(cat(3,ephysStruct(revWIdx).PSTH)));


wConStat = []; pConStat = []; pIncStat = []; wIncStat = [];
%%
cuePreWindow = (cueOn-2)*10:(cueOn)*10;
cueWindow = (cueOn)*10:(cueOn + 2)*10;
stat = [];
modIdx = [];
pIdx = [];

for i = 1:length(ephysStruct)
    modIdx = []; pIdx = []; stat = []; Stat = [];
    
    if i ~= 19
        psth = squeeze(nanmean(ephysStruct(i).PSTH,1));
        for k = 1:min(size(psth))
            [modIdx(k),pIdx(k),~,stat] = ttest(psth(cueWindow,k),psth(cuePreWindow,k));
            Stat(k) = stat.tstat;
            
        end
        modIdx = fdr_bh(pIdx);
        ephysStruct(i).modStatCue = Stat;
        ephysStruct(i).modSumCue = sum(modIdx == 1);
        ephysStruct(i).modPropCue = sum(modIdx == 1) / min(size(psth));
        ephysStruct(i).modIdxCue = modIdx;
        ephysStruct(i).modStatCueIdx = Stat(logical(modIdx));
        ephysStruct(i).posModCue = ephysStruct(i).modStatCueIdx(ephysStruct(i).modStatCueIdx > 0);
        ephysStruct(i).negModCue = ephysStruct(i).modStatCueIdx(ephysStruct(i).modStatCueIdx < 0);
        ephysStruct(i).posPropCue = length(ephysStruct(i).posModCue) / min(size(psth));
        ephysStruct(i).negPropCue = length(ephysStruct(i).negModCue) / min(size(psth));

    end
end

cuePreWindow = (sipDescent-2)*10:(sipDescent)*10;
cueWindow = (sipDescent)*10:(sipDescent + 2)*10;
stat = [];
modIdx = [];
pIdx = [];

for i = 1:length(ephysStruct)
    modIdx = []; pIdx = []; stat = []; Stat = [];
    
    if i ~= 19
        psth = squeeze(nanmean(ephysStruct(i).PSTH,1));
        for k = 1:min(size(psth))
            [modIdx(k),pIdx(k),~,stat] = ttest(psth(cueWindow,k),psth(cuePreWindow,k));
            Stat(k) = stat.tstat;
            
        end
        modIdx = fdr_bh(pIdx);
        ephysStruct(i).modStatSip = Stat;
        ephysStruct(i).modSumSip = sum(modIdx == 1);
        ephysStruct(i).modPropSip = sum(modIdx == 1) / min(size(psth));
        ephysStruct(i).modIdxSip = modIdx;
        ephysStruct(i).modStatSipIdx = Stat(logical(modIdx));
        ephysStruct(i).posModSip = ephysStruct(i).modStatSipIdx(ephysStruct(i).modStatSipIdx > 0);
        ephysStruct(i).negModSip = ephysStruct(i).modStatSipIdx(ephysStruct(i).modStatSipIdx < 0);
        ephysStruct(i).posPropSip = length(ephysStruct(i).posModSip) / min(size(psth));
        ephysStruct(i).negPropSip = length(ephysStruct(i).negModSip) / min(size(psth));

    end
end

%% 
wConPosMod = [ephysStruct(regWIdx).posModCue];  wIncPosMod = [ephysStruct(revWIdx).posModCue];
wConNegMod = [ephysStruct(regWIdx).negModCue];  wIncNegMod = [ephysStruct(revWIdx).negModCue];
        
pConPosMod = [ephysStruct(regPIdx).posModCue];  pIncPosMod = [ephysStruct(revPIdx).posModCue];
pConNegMod = [ephysStruct(regPIdx).negModCue];  pIncNegMod = [ephysStruct(revPIdx).negModCue];

testData = [wConPosMod abs(wConNegMod) wIncPosMod abs(wIncNegMod) pConPosMod abs(pConNegMod) pIncPosMod abs(pIncNegMod)];

g1 = [repmat({'Wistar'},1,length([wConPosMod wConNegMod wIncPosMod wIncNegMod])) repmat({'P rat'},1,length([pConPosMod pConNegMod pIncPosMod pIncNegMod]))];
g2 = [repmat({'Congruent'},1,length([wConPosMod wConNegMod])) repmat({'Incongruent'},1,length([wIncPosMod wIncNegMod])) repmat({'Congruent'},1,length([pConPosMod pConNegMod])) repmat({'Incongruent'},1,length([pIncPosMod pIncNegMod]))];
g3 = [repmat({'Positive'},1,length([wConPosMod])) repmat({'Negative'},1,length([wConNegMod])) repmat({'Positive'},1,length([wIncPosMod])) repmat({'Negative'},1,length([wIncNegMod])) repmat({'Positive'},1,length([pConPosMod])) ...
    repmat({'Negative'},1,length([pConNegMod])) repmat({'Positive'},1,length([pIncPosMod])) repmat({'Negative'},1,length([pIncNegMod]))];
if length(testData) == length(g1) && length(testData) == length(g2)
    disp('Data sizes are matched')
    [p,tbl,stat] = anovan(testData,{g1,g2,g3},'model','full');
    multcompare(stat,'Dimension',[1 2])
        % Write string to file
    tblStr = formattedDisplayText(tbl); 
    fid = fopen('anovaTable_3x2_cueModulation.txt', 'wt');
    fileCleanup = onCleanup(@()fclose(fid));
    formatSpec = '%s\n';
    fprintf(fid, formatSpec, tblStr);
    clear('fileCleanup')
end

%
wConPosProp = [ephysStruct(regWIdx).posPropCue];  wIncPosProp = [ephysStruct(revWIdx).posPropCue];
wConNegProp = [ephysStruct(regWIdx).negPropCue];  wIncNegProp = [ephysStruct(revWIdx).negPropCue];
        
pConPosProp = [ephysStruct(regPIdx).posPropCue];  pIncPosProp = [ephysStruct(revPIdx).posPropCue];
pConNegProp = [ephysStruct(regPIdx).negPropCue];  pIncNegProp = [ephysStruct(revPIdx).negPropCue];


testData = [wConPosProp wConNegProp wIncPosProp wIncNegProp pConPosProp pConNegProp pIncPosProp pIncNegProp];

g1 = [repmat({'Wistar'},1,length([wConPosProp wConNegProp wIncPosProp wIncNegProp])) repmat({'P rat'},1,length([pConPosProp pConNegProp pIncPosProp pIncNegProp]))];
g2 = [repmat({'Congruent'},1,length([wConPosProp wConNegProp])) repmat({'Incongruent'},1,length([wIncPosProp wIncNegProp])) repmat({'Congruent'},1,length([pConPosProp pConNegProp])) repmat({'Incongruent'},1,length([pIncPosProp pIncNegProp]))];
g3 = [repmat({'Positive'},1,length([wConPosProp])) repmat({'Negative'},1,length([wConNegProp])) repmat({'Positive'},1,length([wIncPosProp])) repmat({'Negative'},1,length([wIncNegProp])) repmat({'Positive'},1,length([pConPosProp])) ...
    repmat({'Negative'},1,length([pConNegProp])) repmat({'Positive'},1,length([pIncPosProp])) repmat({'Negative'},1,length([pIncNegProp]))];
if length(testData) == length(g1) && length(testData) == length(g2)
    disp('Data sizes are matched')
    [p,tbl,stat] = anovan(testData,{g1,g2,g3},'model','full');
    multcompare(stat,'Dimension',[1 2])
        % Write string to file
    tblStr = formattedDisplayText(tbl); 
    fid = fopen('anovaTable_3x2_cueModulationProportion.txt', 'wt');
    fileCleanup = onCleanup(@()fclose(fid));
    formatSpec = '%s\n';
    fprintf(fid, formatSpec, tblStr);
    clear('fileCleanup')
end


%% Attempt to plot this monstrosity 
% Plot figure
linSpaceOffset = 0.01;

figure('Units','normalized','Position',[0 0 1 1])
hb = bar([nanmean(wConPosMod) nanmean(abs(wConNegMod)) nanmean(wIncPosMod) nanmean(abs(wIncNegMod)); nanmean(pConPosMod) nanmean(abs(pConNegMod)) nanmean(pIncPosMod) nanmean(abs(pIncNegMod))],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);
hb(1).FaceColor = [0.1 0.1 0.1]; hb(2).FaceColor = [0.3 0.3 0.3]; hb(3).FaceColor = [0.7 0.7 0.7]; hb(4).FaceColor = [0.9 0.9 0.9];
offsetPos = [1+hb(1).XOffset 1+hb(2).XOffset 1+hb(3).XOffset 1+hb(4).XOffset 2+hb(1).XOffset 2+hb(2).XOffset 2+hb(3).XOffset 2+hb(4).XOffset];
scatterXData = {linspace(offsetPos(1)-linSpaceOffset,offsetPos(1)+linSpaceOffset,length(wConPosMod)), ... 
    linspace(offsetPos(2)-linSpaceOffset,offsetPos(2)+linSpaceOffset,length(wConNegMod)), ...
    linspace(offsetPos(3)-linSpaceOffset,offsetPos(3)+linSpaceOffset,length(wIncPosMod)), ...
    linspace(offsetPos(4)-linSpaceOffset,offsetPos(4)+linSpaceOffset,length(wIncNegMod)) ...
    linspace(offsetPos(5)-linSpaceOffset,offsetPos(5)+linSpaceOffset,length(pConPosMod)), ...
    linspace(offsetPos(6)-linSpaceOffset,offsetPos(6)+linSpaceOffset,length(pConNegMod)), ...
    linspace(offsetPos(7)-linSpaceOffset,offsetPos(7)+linSpaceOffset,length(pIncPosMod)), ...
    linspace(offsetPos(8)-linSpaceOffset,offsetPos(8)+linSpaceOffset,length(pIncNegMod))};
hold on
plot([scatterXData{1} scatterXData{2} scatterXData{3} scatterXData{4}]',[wConPosMod abs(wConNegMod) wIncPosMod abs(wIncNegMod)],'ko','LineWidth',2,'MarkerSize',2)
plot([scatterXData{5} scatterXData{6} scatterXData{7} scatterXData{8}]',[pConPosMod abs(pConNegMod) pIncPosMod abs(pIncNegMod)],'ko','LineWidth',2,'MarkerSize',2)

offsetPos = [1+hb(1).XOffset 1+hb(2).XOffset 1+hb(3).XOffset 1+hb(4).XOffset; 2+hb(1).XOffset 2+hb(2).XOffset 2+hb(3).XOffset 2+hb(4).XOffset];
er = errorbar(offsetPos,[nanmean(wConPosMod) nanmean(abs(wConNegMod)) nanmean(wIncPosMod) nanmean(abs(wIncNegMod)); nanmean(pConPosMod) nanmean(abs(pConNegMod)) nanmean(pIncPosMod) nanmean(abs(pIncNegMod))], ... 
     [nanstd(wConPosMod)/sqrt(length(wConPosMod)) nanstd(abs(wConNegMod))/sqrt(length(wConNegMod)) nanstd(wIncPosMod)/sqrt(length(wIncPosMod)) nanstd(abs(wIncNegMod))/sqrt(length(wIncNegMod)) ...
     ;nanstd(pConPosMod)/sqrt(length(pConPosMod)) nanstd(abs(pConNegMod))/sqrt(length(pConNegMod)) nanstd(pIncPosMod)/sqrt(length(pIncPosMod)) nanstd(abs(pIncNegMod))/sqrt(length(pIncNegMod))],'LineWidth',3);    

for i = 1:length(er)
    er(i).Color = [0 0 0];                        
    er(i).LineStyle = 'none';
end
legend([{'Congruent Positive'},{'Congruent Negative'},{'Incongruent Positive'},{'Incongruent Negative'}])
ylabel('Modulation Index')
ylim([0 40])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
saveas(gca,['modulationStrength_Cue'],'svg')

% Proportion
figure('Units','normalized','Position',[0 0 1 1])
hb = bar([nanmean(wConPosProp) nanmean(abs(wConNegProp)) nanmean(wIncPosProp) nanmean(abs(wIncNegProp)); nanmean(pConPosProp) nanmean(abs(pConNegProp)) nanmean(pIncPosProp) nanmean(abs(pIncNegProp))],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);
hb(1).FaceColor = [0.1 0.1 0.1]; hb(2).FaceColor = [0.3 0.3 0.3]; hb(3).FaceColor = [0.7 0.7 0.7]; hb(4).FaceColor = [0.9 0.9 0.9];
offsetPos = [1+hb(1).XOffset 1+hb(2).XOffset 1+hb(3).XOffset 1+hb(4).XOffset 2+hb(1).XOffset 2+hb(2).XOffset 2+hb(3).XOffset 2+hb(4).XOffset];
scatterXData = {linspace(offsetPos(1)-linSpaceOffset,offsetPos(1)+linSpaceOffset,length(wConPosProp)), ... 
    linspace(offsetPos(2)-linSpaceOffset,offsetPos(2)+linSpaceOffset,length(wConNegProp)), ...
    linspace(offsetPos(3)-linSpaceOffset,offsetPos(3)+linSpaceOffset,length(wIncPosProp)), ...
    linspace(offsetPos(4)-linSpaceOffset,offsetPos(4)+linSpaceOffset,length(wIncNegProp)) ...
    linspace(offsetPos(5)-linSpaceOffset,offsetPos(5)+linSpaceOffset,length(pConPosProp)), ...
    linspace(offsetPos(6)-linSpaceOffset,offsetPos(6)+linSpaceOffset,length(pConNegProp)), ...
    linspace(offsetPos(7)-linSpaceOffset,offsetPos(7)+linSpaceOffset,length(pIncPosProp)), ...
    linspace(offsetPos(8)-linSpaceOffset,offsetPos(8)+linSpaceOffset,length(pIncNegProp))};
hold on
plot([scatterXData{1} scatterXData{2} scatterXData{3} scatterXData{4}]',[wConPosProp abs(wConNegProp) wIncPosProp abs(wIncNegProp)],'ko','LineWidth',2,'MarkerSize',10)
plot([scatterXData{5} scatterXData{6} scatterXData{7} scatterXData{8}]',[pConPosProp abs(pConNegProp) pIncPosProp abs(pIncNegProp)],'ko','LineWidth',2,'MarkerSize',10)

offsetPos = [1+hb(1).XOffset 1+hb(2).XOffset 1+hb(3).XOffset 1+hb(4).XOffset; 2+hb(1).XOffset 2+hb(2).XOffset 2+hb(3).XOffset 2+hb(4).XOffset];
er = errorbar(offsetPos,[nanmean(wConPosProp) nanmean(abs(wConNegProp)) nanmean(wIncPosProp) nanmean(abs(wIncNegProp)); nanmean(pConPosProp) nanmean(abs(pConNegProp)) nanmean(pIncPosProp) nanmean(abs(pIncNegProp))], ... 
     [nanstd(wConPosProp)/sqrt(length(wConPosProp)) nanstd(abs(wConNegProp))/sqrt(length(wConNegProp)) nanstd(wIncPosProp)/sqrt(length(wIncPosProp)) nanstd(abs(wIncNegProp))/sqrt(length(wIncNegProp)) ...
     ;nanstd(pConPosProp)/sqrt(length(pConPosProp)) nanstd(abs(pConNegProp))/sqrt(length(pConNegProp)) nanstd(pIncPosProp)/sqrt(length(pIncPosProp)) nanstd(abs(pIncNegProp))/sqrt(length(pIncNegProp))],'LineWidth',3);    

for i = 1:length(er)
    er(i).Color = [0 0 0];                        
    er(i).LineStyle = 'none';
end
legend([{'Congruent Positive'},{'Congruent Negative'},{'Incongruent Positive'},{'Incongruent Negative'}])
ylabel('Proportion of Neurons')
ylim([0 1])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
saveas(gca,['modulationStrength_CueProp'],'svg')

%% 
wConPosModSip = [ephysStruct(regWIdx).posModSip];  wIncPosModSip = [ephysStruct(revWIdx).posModSip];
wConNegModSip = [ephysStruct(regWIdx).negModSip];  wIncNegModSip = [ephysStruct(revWIdx).negModSip];
        
pConPosModSip = [ephysStruct(regPIdx).posModSip];  pIncPosModSip = [ephysStruct(revPIdx).posModSip];
pConNegModSip = [ephysStruct(regPIdx).negModSip];  pIncNegModSip = [ephysStruct(revPIdx).negModSip];

testData = [wConPosModSip abs(wConNegModSip) wIncPosModSip abs(wIncNegModSip) pConPosModSip abs(pConNegModSip) pIncPosModSip abs(pIncNegModSip)];

g1 = [repmat({'Wistar'},1,length([wConPosModSip wConNegModSip wIncPosModSip wIncNegModSip])) repmat({'P rat'},1,length([pConPosModSip pConNegModSip pIncPosModSip pIncNegModSip]))];
g2 = [repmat({'Congruent'},1,length([wConPosModSip wConNegModSip])) repmat({'Incongruent'},1,length([wIncPosModSip wIncNegModSip])) repmat({'Congruent'},1,length([pConPosModSip pConNegModSip])) repmat({'Incongruent'},1,length([pIncPosModSip pIncNegModSip]))];
g3 = [repmat({'Positive'},1,length([wConPosModSip])) repmat({'Negative'},1,length([wConNegModSip])) repmat({'Positive'},1,length([wIncPosModSip])) repmat({'Negative'},1,length([wIncNegModSip])) repmat({'Positive'},1,length([pConPosModSip])) ...
    repmat({'Negative'},1,length([pConNegModSip])) repmat({'Positive'},1,length([pIncPosModSip])) repmat({'Negative'},1,length([pIncNegModSip]))];
if length(testData) == length(g1) && length(testData) == length(g2)
    disp('Data sizes are matched')
    [p,tbl,stat] = anovan(testData,{g1,g2,g3},'model','full');
    multcompare(stat,'Dimension',[1 3])
        % Write string to file
    tblStr = formattedDisplayText(tbl); 
    fid = fopen('anovaTable_3x2_sipperModulation.txt', 'wt');
    fileCleanup = onCleanup(@()fclose(fid));
    formatSpec = '%s\n';
    fprintf(fid, formatSpec, tblStr);
    clear('fileCleanup')
end


%
wConPosPropSip = [ephysStruct(regWIdx).posPropSip];  wIncPosPropSip = [ephysStruct(revWIdx).posPropSip];
wConNegPropSip = [ephysStruct(regWIdx).negPropSip];  wIncNegPropSip = [ephysStruct(revWIdx).negPropSip];
        
pConPosPropSip = [ephysStruct(regPIdx).posPropSip];  pIncPosPropSip = [ephysStruct(revPIdx).posPropSip];
pConNegPropSip = [ephysStruct(regPIdx).negPropSip];  pIncNegPropSip = [ephysStruct(revPIdx).negPropSip];

testData = [wConPosPropSip wConNegPropSip wIncPosPropSip wIncNegPropSip pConPosPropSip pConNegPropSip pIncPosPropSip pIncNegPropSip];

g1 = [repmat({'Wistar'},1,length([wConPosPropSip wConNegPropSip wIncPosPropSip wIncNegPropSip])) repmat({'P rat'},1,length([pConPosPropSip pConNegPropSip pIncPosPropSip pIncNegPropSip]))];
g2 = [repmat({'Congruent'},1,length([wConPosPropSip wConNegPropSip])) repmat({'Incongruent'},1,length([wIncPosPropSip wIncNegPropSip])) repmat({'Congruent'},1,length([pConPosPropSip pConNegPropSip])) repmat({'Incongruent'},1,length([pIncPosPropSip pIncNegPropSip]))];
g3 = [repmat({'Positive'},1,length([wConPosPropSip])) repmat({'Negative'},1,length([wConNegPropSip])) repmat({'Positive'},1,length([wIncPosPropSip])) repmat({'Negative'},1,length([wIncNegPropSip])) repmat({'Positive'},1,length([pConPosPropSip])) ...
    repmat({'Negative'},1,length([pConNegPropSip])) repmat({'Positive'},1,length([pIncPosPropSip])) repmat({'Negative'},1,length([pIncNegPropSip]))];
if length(testData) == length(g1) && length(testData) == length(g2)
    disp('Data sizes are matched')
    [p,tbl,stat] = anovan(testData,{g1,g2,g3},'model','full');
    multcompare(stat,'Dimension',[1 2])
        % Write string to file
    tblStr = formattedDisplayText(tbl); 
    fid = fopen('anovaTable_2x2_sipperModulationProportion.txt', 'wt');
    fileCleanup = onCleanup(@()fclose(fid));
    formatSpec = '%s\n';
    fprintf(fid, formatSpec, tblStr);
    clear('fileCleanup')
end

%% Attempt to plot this monstrosity 
% Plot figure
figure('Units','normalized','Position',[0 0 1 1])
hb = bar([nanmean(wConPosModSip) nanmean(abs(wConNegModSip)) nanmean(wIncPosModSip) nanmean(abs(wIncNegModSip)); nanmean(pConPosModSip) nanmean(abs(pConNegModSip)) nanmean(pIncPosModSip) nanmean(abs(pIncNegModSip))],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);
hb(1).FaceColor = [0.1 0.1 0.1]; hb(2).FaceColor = [0.3 0.3 0.3]; hb(3).FaceColor = [0.7 0.7 0.7]; hb(4).FaceColor = [0.9 0.9 0.9];
offsetPos = [1+hb(1).XOffset 1+hb(2).XOffset 1+hb(3).XOffset 1+hb(4).XOffset 2+hb(1).XOffset 2+hb(2).XOffset 2+hb(3).XOffset 2+hb(4).XOffset];
scatterXData = {linspace(offsetPos(1)-linSpaceOffset,offsetPos(1)+linSpaceOffset,length(wConPosModSip)), ... 
    linspace(offsetPos(2)-linSpaceOffset,offsetPos(2)+linSpaceOffset,length(wConNegModSip)), ...
    linspace(offsetPos(3)-linSpaceOffset,offsetPos(3)+linSpaceOffset,length(wIncPosModSip)), ...
    linspace(offsetPos(4)-linSpaceOffset,offsetPos(4)+linSpaceOffset,length(wIncNegModSip)) ...
    linspace(offsetPos(5)-linSpaceOffset,offsetPos(5)+linSpaceOffset,length(pConPosModSip)), ...
    linspace(offsetPos(6)-linSpaceOffset,offsetPos(6)+linSpaceOffset,length(pConNegModSip)), ...
    linspace(offsetPos(7)-linSpaceOffset,offsetPos(7)+linSpaceOffset,length(pIncPosModSip)), ...
    linspace(offsetPos(8)-linSpaceOffset,offsetPos(8)+linSpaceOffset,length(pIncNegModSip))};
hold on
plot([scatterXData{1} scatterXData{2} scatterXData{3} scatterXData{4}]',[wConPosModSip abs(wConNegModSip) wIncPosModSip abs(wIncNegModSip)],'ko','LineWidth',2,'MarkerSize',2)
plot([scatterXData{5} scatterXData{6} scatterXData{7} scatterXData{8}]',[pConPosModSip abs(pConNegModSip) pIncPosModSip abs(pIncNegModSip)],'ko','LineWidth',2,'MarkerSize',2)

offsetPos = [1+hb(1).XOffset 1+hb(2).XOffset 1+hb(3).XOffset 1+hb(4).XOffset; 2+hb(1).XOffset 2+hb(2).XOffset 2+hb(3).XOffset 2+hb(4).XOffset];
er = errorbar(offsetPos,[nanmean(wConPosModSip) nanmean(abs(wConNegModSip)) nanmean(wIncPosModSip) nanmean(abs(wIncNegModSip)); nanmean(pConPosModSip) nanmean(abs(pConNegModSip)) nanmean(pIncPosModSip) nanmean(abs(pIncNegModSip))], ... 
     [nanstd(wConPosModSip)/sqrt(length(wConPosModSip)) nanstd(abs(wConNegModSip))/sqrt(length(wConNegModSip)) nanstd(wIncPosModSip)/sqrt(length(wIncPosModSip)) nanstd(abs(wIncNegModSip))/sqrt(length(wIncNegModSip)) ...
     ;nanstd(pConPosModSip)/sqrt(length(pConPosModSip)) nanstd(abs(pConNegModSip))/sqrt(length(pConNegModSip)) nanstd(pIncPosModSip)/sqrt(length(pIncPosModSip)) nanstd(abs(pIncNegModSip))/sqrt(length(pIncNegModSip))],'LineWidth',3);    

for i = 1:length(er)
    er(i).Color = [0 0 0];                        
    er(i).LineStyle = 'none';
end
legend([{'Congruent Positive'},{'Congruent Negative'},{'Incongruent Positive'},{'Incongruent Negative'}])
ylabel('Modulation Index')
ylim([0 20])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
saveas(gca,['modulationStrength_Sip'],'svg')

% Proportion
figure('Units','normalized','Position',[0 0 1 1])
hb = bar([nanmean(wConPosPropSip) nanmean(abs(wConNegPropSip)) nanmean(wIncPosPropSip) nanmean(abs(wIncNegPropSip)); nanmean(pConPosPropSip) nanmean(abs(pConNegPropSip)) nanmean(pIncPosPropSip) nanmean(abs(pIncNegPropSip))],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);
hb(1).FaceColor = [0.1 0.1 0.1]; hb(2).FaceColor = [0.3 0.3 0.3]; hb(3).FaceColor = [0.7 0.7 0.7]; hb(4).FaceColor = [0.9 0.9 0.9];
offsetPos = [1+hb(1).XOffset 1+hb(2).XOffset 1+hb(3).XOffset 1+hb(4).XOffset 2+hb(1).XOffset 2+hb(2).XOffset 2+hb(3).XOffset 2+hb(4).XOffset];
scatterXData = {linspace(offsetPos(1)-linSpaceOffset,offsetPos(1)+linSpaceOffset,length(wConPosPropSip)), ... 
    linspace(offsetPos(2)-linSpaceOffset,offsetPos(2)+linSpaceOffset,length(wConNegPropSip)), ...
    linspace(offsetPos(3)-linSpaceOffset,offsetPos(3)+linSpaceOffset,length(wIncPosPropSip)), ...
    linspace(offsetPos(4)-linSpaceOffset,offsetPos(4)+linSpaceOffset,length(wIncNegPropSip)) ...
    linspace(offsetPos(5)-linSpaceOffset,offsetPos(5)+linSpaceOffset,length(pConPosPropSip)), ...
    linspace(offsetPos(6)-linSpaceOffset,offsetPos(6)+linSpaceOffset,length(pConNegPropSip)), ...
    linspace(offsetPos(7)-linSpaceOffset,offsetPos(7)+linSpaceOffset,length(pIncPosPropSip)), ...
    linspace(offsetPos(8)-linSpaceOffset,offsetPos(8)+linSpaceOffset,length(pIncNegPropSip))};
hold on
plot([scatterXData{1} scatterXData{2} scatterXData{3} scatterXData{4}]',[wConPosPropSip abs(wConNegPropSip) wIncPosPropSip abs(wIncNegPropSip)],'ko','LineWidth',2,'MarkerSize',10)
plot([scatterXData{5} scatterXData{6} scatterXData{7} scatterXData{8}]',[pConPosPropSip abs(pConNegPropSip) pIncPosPropSip abs(pIncNegPropSip)],'ko','LineWidth',2,'MarkerSize',10)

offsetPos = [1+hb(1).XOffset 1+hb(2).XOffset 1+hb(3).XOffset 1+hb(4).XOffset; 2+hb(1).XOffset 2+hb(2).XOffset 2+hb(3).XOffset 2+hb(4).XOffset];
er = errorbar(offsetPos,[nanmean(wConPosPropSip) nanmean(abs(wConNegPropSip)) nanmean(wIncPosPropSip) nanmean(abs(wIncNegPropSip)); nanmean(pConPosPropSip) nanmean(abs(pConNegPropSip)) nanmean(pIncPosPropSip) nanmean(abs(pIncNegPropSip))], ... 
     [nanstd(wConPosPropSip)/sqrt(length(wConPosPropSip)) nanstd(abs(wConNegPropSip))/sqrt(length(wConNegPropSip)) nanstd(wIncPosPropSip)/sqrt(length(wIncPosPropSip)) nanstd(abs(wIncNegPropSip))/sqrt(length(wIncNegPropSip)) ...
     ;nanstd(pConPosPropSip)/sqrt(length(pConPosPropSip)) nanstd(abs(pConNegPropSip))/sqrt(length(pConNegPropSip)) nanstd(pIncPosPropSip)/sqrt(length(pIncPosPropSip)) nanstd(abs(pIncNegPropSip))/sqrt(length(pIncNegPropSip))],'LineWidth',3);    

for i = 1:length(er)
    er(i).Color = [0 0 0];                        
    er(i).LineStyle = 'none';
end
legend([{'Congruent Positive'},{'Congruent Negative'},{'Incongruent Positive'},{'Incongruent Negative'}])
ylabel('Proportion of Neurons')
ylim([0 1])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
saveas(gca,['modulationStrength_SipProp'],'svg')

%%

for i = 1:length(pCon)
    [modIdx(i),pIdx(i),~,stat] = ttest(pCon(cuePreWindow,i),pCon(cueWindow,i)); 
    pConStat(i) = stat.tstat;
end
modIdx = fdr_bh(pIdx);
pConSum = sum(modIdx == 1);
pConModIdx = modIdx;

modIdx = [];
stat = [];
pIdx = [];

for i = 1:length(pInc)
    [modIdx(i),pIdx(i),~,stat] = ttest(pInc(cuePreWindow,i),pInc(cueWindow,i));
    
    pIncStat(i) = stat.tstat;
end
modIdx = fdr_bh(pIdx);
pIncSum = sum(modIdx == 1);
pIncModIdx = modIdx;

modIdx = [];
stat = [];
pIdx = [];

for i = 1:length(wCon)
    [modIdx(i),pIdx(i),~,stat] = ttest(wCon(cuePreWindow,i),wCon(cueWindow,i)); 
    
    wConStat(i) = stat.tstat;
end
modIdx = fdr_bh(pIdx);
wConSum = sum(modIdx == 1);
wConModIdx = modIdx;

modIdx = [];
stat = [];
pIdx = [];

for i = 1:length(wInc)
    [modIdx(i),pIdx(i),~,stat] = ttest(wInc(cuePreWindow,i),wInc(cueWindow,i)); 
    wIncStat(i) = stat.tstat;
end
modIdx = fdr_bh(pIdx);
wIncSum = sum(modIdx == 1);
wIncModIdx = modIdx;

bar([pConSum/max(size(pCon)) pIncSum/max(size(pInc)) wConSum/max(size(wCon)) wIncSum/max(size(wInc))])

%% Find positive vs. negative proportion
pConIdx = find(pConModIdx == 1);
pConMod = pConStat(pConIdx); pConModPosIdx = pConMod > 0; pConModNegIdx = pConMod < 0;
pConModPos = pConMod(pConModPosIdx); pConModNeg = pConMod(pConModNegIdx);

pIncIdx = find(pIncModIdx == 1);
pIncMod = pIncStat(pIncIdx); pIncModPosIdx = pIncMod > 0; pIncModNegIdx = pIncMod < 0;
pIncModPos = pIncMod(pIncModPosIdx); pIncModNeg = pIncMod(pIncModNegIdx);

%

wConIdx = find(wConModIdx == 1);
wConMod = wConStat(wConIdx); wConModPosIdx = wConMod > 0; wConModNegIdx = wConMod < 0;
wConModPos = wConMod(wConModPosIdx); wConModNeg = wConMod(wConModNegIdx);

wIncIdx = find(wIncModIdx == 1);
wIncMod = wIncStat(wIncIdx); wIncModPosIdx = wIncMod > 0; wIncModNegIdx = wIncMod < 0;
wIncModPos = wIncMod(wIncModPosIdx); wIncModNeg = wIncMod(wIncModNegIdx);

dataStackedBars = [length(wConModPos)/length(wCon) length(wConModNeg)/length(wCon); length(wIncModPos)/length(wInc) length(wIncModNeg)/length(wInc); length(pConModPos)/length(pCon) length(pConModNeg)/length(pCon); length(pIncModPos)/length(pInc) length(pIncModNeg)/length(pInc)];
bar(dataStackedBars)
ylim([0 1])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',25,'FontWeight','bold','LineWidth',5)



dataGroupedBars = [nanmean(wConModPos) nanmean(abs(wConModNeg)); nanmean(wIncModPos) nanmean(abs(wIncModNeg)); nanmean(pConModPos) nanmean(abs(pConModNeg)); nanmean(pIncModPos) nanmean(abs(pIncModNeg))];
bar(dataGroupedBars)

%% Attempt to make a table in order to visualize the amounts of positive / negative modulated neurons by session type 
dataSet = array2table([wConModPos; wIncModPos; pConModPos; pIncModPos])
%% Look at nanmean traces of neurons over trial blocks
idxN = [1:5; 6:10; 11:15];
pCon = cat(3,ephysStruct(regPIdx).PSTH);
pInc = cat(3,ephysStruct(revPIdx).PSTH);
wCon = cat(3,ephysStruct(regWIdx).PSTH);
wInc = cat(3,ephysStruct(revWIdx).PSTH);

colors = {'b','r','k'};
figure('Units','normalized','Position',[0 0 1 1])


for i = 1:min(size(idxN))
    WConditionsPre_Con = squeeze(nanmean(wCon(idxN(i),1:91,:),1));
    WConditionsPre_Inc = squeeze(nanmean(wInc(idxN(i),1:91,:),1));
    PConditionsPre_Con = squeeze(nanmean(pCon(idxN(i),1:91,:),1));
    PConditionsPre_Inc = squeeze(nanmean(pInc(idxN(i),1:91,:),1));

    WConData{i} = WConditionsPre_Con;
    trialBlockLabelWC = [repmat({'Block1'},1,max(size(WConditionsPre_Con))) repmat({'Block2'},1,max(size(WConditionsPre_Con))) repmat({'Block3'},1,max(size(WConditionsPre_Con)))];
    WIncData{i} = WConditionsPre_Inc;
    trialBlockLabelWI = [repmat({'Block1'},1,max(size(WConditionsPre_Inc))) repmat({'Block2'},1,max(size(WConditionsPre_Inc))) repmat({'Block3'},1,max(size(WConditionsPre_Inc)))];



end

trialBlockLabelWC = [repmat({'Block1'},1,max(size(WConditionsPre_Con))) repmat({'Block2'},1,max(size(WConditionsPre_Con))) repmat({'Block3'},1,max(size(WConditionsPre_Con)))];
trialBlockLabelWI = [repmat({'Block1'},1,max(size(WConditionsPre_Inc))) repmat({'Block2'},1,max(size(WConditionsPre_Inc))) repmat({'Block3'},1,max(size(WConditionsPre_Inc)))];
dataWC = [WConData{:}];
tableWC = array2table(dataWC');
tableWC.Block = trialBlockLabelWC';
TimeVar = table([1:91]','VariableNames',{'Time'});

rm = fitrm(tableWC,'Var1-Var91 ~ Block')%,'WithinDesign',TimeVar);
anovaTbl = ranova(rm,'WithinModel','Time')

dataWI = [WIncData{:}];
tableWI = array2table(dataWI');
tableWI.Block = trialBlockLabelWI';
rm = fitrm(tableWI,'Var1-Var91 ~ Block')%,'WithinDesign',TimeVar);
anovaTbl = ranova(rm,'WithinModel','Time')






for i = 1:min(size(idxN))
    WConditionsPre_Con = squeeze(nanmean(wCon(idxN(i),1:91,:),1));
    WConditionsPre_Inc = squeeze(nanmean(wInc(idxN(i),1:91,:),1));
    PConditionsPre_Con = squeeze(nanmean(pCon(idxN(i),1:91,:),1));
    PConditionsPre_Inc = squeeze(nanmean(pInc(idxN(i),1:91,:),1));

    % TTest Permutations 
    for k = 1:length(WConData)
        [h12wc,p12wc,stat12wc] = ttest(WConData{1}',WConData{2}');
        [h13wc,p13wc,stat13wc] = ttest(WConData{1}',WConData{3}');
        [h23wc,p23wc,stat23wc] = ttest(WConData{2}',WConData{3}');

    end

    h12wcfdr = fdr_bh(p12wc);
    h13wcfdr = fdr_bh(p13wc);
    h23wcfdr = fdr_bh(p23wc);

    %

    for k = 1:length(WIncData)
        [h12wi,p12wi,stat12wi] = ttest(WIncData{1}',WIncData{2}');
        [h13wi,p13wi,stat13wi] = ttest(WIncData{1}',WIncData{3}');
        [h23wi,p23wi,stat23wi] = ttest(WIncData{2}',WIncData{3}');
    end

    h12wifdr = fdr_bh(p12wi);
    h13wifdr = fdr_bh(p13wi);
    h23wifdr = fdr_bh(p23wi);

    %
    

    time = (1:size(WConditionsPre_Con,1))/10;

    
    semWPre_Con = std(WConditionsPre_Con,[],2)/sqrt(length(WConditionsPre_Con));
    semWPre_Inc = std(WConditionsPre_Inc,[],2)/sqrt(length(WConditionsPre_Inc));

    semPPre_Con = std(PConditionsPre_Con,[],2)/sqrt(length(PConditionsPre_Con));
    semPPre_Inc = std(PConditionsPre_Inc,[],2)/sqrt(length(PConditionsPre_Inc));

    subplot(1,2,1)
    shadedErrorBar(time,smooth(nanmean(WConditionsPre_Con,2),1),semWPre_Con,'lineprops',{colors{i},'LineWidth',3}); 
    hold on; 
    xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
    xlabel('Time (s)')
    ylabel('Firing Rate (Hz)')
    axis tight
    ylim([0.4 1.1])
    legend({'Trials 1-5','Trials 6-10','Trials 11-15'},'AutoUpdate','off')
    title('Congruent')
    if i == 1
        plot(time,double(h12wcfdr) - 0.55,'k*','MarkerSize',10)
        plot(time,double(h13wcfdr) - 0.53,'b*','MarkerSize',10)
        plot(time,double(h23wcfdr) - 0.51,'r*','MarkerSize',10)
        yline(0.45,'k--','Trials 1-5 v 6-10')
        yline(0.47,'k--','Trials 1-5 v 11-15')
        yline(0.49,'k--','Trials 6-10 v 11-15')

    end
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',25,'FontWeight','bold','LineWidth',5)

    subplot(1,2,2)
    shadedErrorBar(time,smooth(nanmean(WConditionsPre_Inc,2),1),semWPre_Con,'lineprops',{colors{i},'LineWidth',3}); 
    hold on; 
    xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
    
    xlabel('Time (s)')
    ylabel('Firing Rate (Hz)')
    axis tight
    ylim([0.4 1.1])
    legend({'Trials 1-5','Trials 6-10','Trials 11-15'})

    title('Incongruent')

    if i == 1
        plot(time,double(h12wifdr) - 0.55,'k*','MarkerSize',10)
        plot(time,double(h13wifdr) - 0.53,'b*','MarkerSize',10)
        plot(time,double(h23wifdr) - 0.51,'r*','MarkerSize',10)
        yline(0.45,'k--','Trials 1-5 v 6-10')
        yline(0.47,'k--','Trials 1-5 v 11-15')
        yline(0.49,'k--','Trials 6-10 v 11-15')

    end

    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',25,'FontWeight','bold','LineWidth',5)


end

saveas(gca,[figPath filesep 'trialnanmeanFiringRates_Wistars'     '_100ms' ],'svg');     saveas(gca,[figPath filesep 'trialnanmeanFiringRates_Wistars'    '_100ms' ],'png')


%% Look at nanmean traces of neurons over trial blocks ( P rats )
colors = {'b','r','k'};
figure('Units','normalized','Position',[0 0 1 1])


for i = 1:min(size(idxN))
    PConditionsPre_Con = squeeze(nanmean(pCon(idxN(i),1:91,:),1));
    PConditionsPre_Inc = squeeze(nanmean(pInc(idxN(i),1:91,:),1));

    PConData{i} = PConditionsPre_Con;
    PIncData{i} = PConditionsPre_Inc;

end

trialBlockLabelPC = [repmat({'Block1'},1,max(size(PConditionsPre_Con))) repmat({'Block2'},1,max(size(PConditionsPre_Con))) repmat({'Block3'},1,max(size(PConditionsPre_Con)))];
trialBlockLabelPI = [repmat({'Block1'},1,max(size(PConditionsPre_Inc))) repmat({'Block2'},1,max(size(PConditionsPre_Inc))) repmat({'Block3'},1,max(size(PConditionsPre_Inc)))];
dataPC = [PConData{:}];
tablePC = array2table(dataPC');
tablePC.Block = trialBlockLabelPC';
TimeVar = table([1:91]','VariableNames',{'Time'});

rm = fitrm(tablePC,'Var1-Var91 ~ Block')%,'WithinDesign',TimeVar);
anovaTbl = ranova(rm,'WithinModel','Time')

dataPI = [PIncData{:}];
tablePI = array2table(dataPI');
tablePI.Block = trialBlockLabelPI';
rm = fitrm(tablePI,'Var1-Var91 ~ Block')%,'WithinDesign',TimeVar);
anovaTbl = ranova(rm,'WithinModel','Time')

for i = 1:min(size(idxN))
    PConditionsPre_Con = squeeze(nanmean(pCon(idxN(i),1:91,:),1));
    PConditionsPre_Inc = squeeze(nanmean(pInc(idxN(i),1:91,:),1));

    % TTest Permutations 
    for k = 1:length(PConData)
        [h12pc,p12pc,stat12pc] = ttest(PConData{1}',PConData{2}');
        [h13pc,p13pc,stat13pc] = ttest(PConData{1}',PConData{3}');
        [h23pc,p23pc,stat23pc] = ttest(PConData{2}',PConData{3}');

    end

    h12pcfdr = fdr_bh(p12pc);
    h13pcfdr = fdr_bh(p13pc);
    h23pcfdr = fdr_bh(p23pc);

    %

    for k = 1:length(PIncData)
        [h12pi,p12pi,stat12pi] = ttest(PIncData{1}',PIncData{2}');
        [h13pi,p13pi,stat13pi] = ttest(PIncData{1}',PIncData{3}');
        [h23pi,p23pi,stat23pi] = ttest(PIncData{2}',PIncData{3}');
    end

    h12pifdr = fdr_bh(p12pi);
    h13pifdr = fdr_bh(p13pi);
    h23pifdr = fdr_bh(p23pi);

    %
    

    time = (1:size(PConditionsPre_Con,1))/10;

    
    semPPre_Con = std(PConditionsPre_Con,[],2)/sqrt(length(PConditionsPre_Con));
    semPPre_Inc = std(PConditionsPre_Inc,[],2)/sqrt(length(PConditionsPre_Inc));

    subplot(1,2,1)
    shadedErrorBar(time,smooth(nanmean(PConditionsPre_Con,2),1),semPPre_Con,'lineprops',{colors{i},'LineWidth',3}); 
    hold on; 
    xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
    xlabel('Time (s)')
    ylabel('Firing Rate (Hz)')
    axis tight
    ylim([0.4 1.1])
    legend({'Trials 1-5','Trials 6-10','Trials 11-15'},'AutoUpdate','off')
    title('Congruent')
    if i == 1
        plot(time,double(h12pcfdr) - 0.55,'k*','MarkerSize',10)
        plot(time,double(h13pcfdr) - 0.53,'b*','MarkerSize',10)
        plot(time,double(h23pcfdr) - 0.51,'r*','MarkerSize',10)
        yline(0.45,'k--','Trials 1-5 v 6-10')
        yline(0.47,'k--','Trials 1-5 v 11-15')
        yline(0.49,'k--','Trials 6-10 v 11-15')

    end
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',25,'FontWeight','bold','LineWidth',5)

    subplot(1,2,2)
    shadedErrorBar(time,smooth(nanmean(PConditionsPre_Inc,2),1),semPPre_Inc,'lineprops',{colors{i},'LineWidth',3}); 
    hold on; 
    xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
    
    xlabel('Time (s)')
    ylabel('Firing Rate (Hz)')
    axis tight
    ylim([0.4 1.1])
    legend({'Trials 1-5','Trials 6-10','Trials 11-15'})

    title('Incongruent')

    if i == 1
        plot(time,double(h12pifdr) - 0.55,'k*','MarkerSize',10)
        plot(time,double(h13pifdr) - 0.53,'b*','MarkerSize',10)
        plot(time,double(h23pifdr) - 0.51,'r*','MarkerSize',10)
        yline(0.45,'k--','Trials 1-5 v 6-10')
        yline(0.47,'k--','Trials 1-5 v 11-15')
        yline(0.49,'k--','Trials 6-10 v 11-15')

    end

    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',25,'FontWeight','bold','LineWidth',5)


end

saveas(gca,[figPath filesep 'trialnanmeanFiringRates_Prats'     '_100ms' ],'svg');     saveas(gca,[figPath filesep 'trialnanmeanFiringRates_Prats'    '_100ms' ],'png')


%%
idxN = [1:5; 6:10; 11:15];
pCon = cat(3,ephysStruct(regPIdx).PSTH);
pInc = cat(3,ephysStruct(revPIdx).PSTH);
wCon = cat(3,ephysStruct(regWIdx).PSTH);
wInc = cat(3,ephysStruct(revWIdx).PSTH);

colors = {'b','r','k'};
figure('Units','normalized','Position',[0 0 1 1])

clear WConditionsPre_Con WConditionsPre_Inc PConditionsPre_Con PConditionsPre_Inc
for i = 1:min(size(idxN))
    WConditionsPre_Con{i} = squeeze(nanmean(wCon(idxN(i),1:91,:),1));
    WConditionsPre_Inc{i} = squeeze(nanmean(wInc(idxN(i),1:91,:),1));
    PConditionsPre_Con{i} = squeeze(nanmean(pCon(idxN(i),1:91,:),1));
    PConditionsPre_Inc{i} = squeeze(nanmean(pInc(idxN(i),1:91,:),1));
end
allData = [WConditionsPre_Con; WConditionsPre_Inc; PConditionsPre_Con; PConditionsPre_Inc];

%
formatData = cell2mat(allData);
formatData = formatData';
datatable = array2table(formatData,'VariableNames',{'Block1','Block2','Block3'});
datatable.Strain = [repmat({'Wistar'},30,1); repmat({'P rat'},30,1)];
condition = [repmat({'Congruent'},15,2); repmat({'Incongruent'},15,2)];
condition = {(condition{:})};
datatable.Condition = condition';

% Create rm model
Block = table([1 2 3]','VariableNames',{'Blocks'});
rm = fitrm(datatable,'Block1-Block3 ~ Strain * Condition','WithinDesign',Block);
anovaTbl = ranova(rm,'WithinModel','Blocks')


%%
%%
%     time = (1:size(WConditionsPre_Con,1))/10;
    time = (80:200)/10;
    
    semWPre_Con = std(WConditionsPre_Con,[],2)/sqrt(length(WConditionsPre_Con));
    semWPre_Inc = std(WConditionsPre_Inc,[],2)/sqrt(length(WConditionsPre_Inc));

    semPPre_Con = std(PConditionsPre_Con,[],2)/sqrt(length(PConditionsPre_Con));
    semPPre_Inc = std(PConditionsPre_Inc,[],2)/sqrt(length(PConditionsPre_Inc));

    subplot(1,2,1)
    shadedErrorBar(time,smooth(nanmean(WConditionsPre_Con,2),7),semWPre_Con,'lineprops',{colors{i},'LineWidth',3}); 
    hold on; 
%     xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
    xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
%     if strcmp(timepoint,'all')
    xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
%     end
    xlabel('Time (s)')
    ylabel('Firing Rate (Hz)')
    axis tight
    ylim([0.5 1.1])
    legend({'Trials 1-5','Trials 6-10','Trials 11-15'})
    title('Congruent')
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',25,'FontWeight','bold','LineWidth',5)

    subplot(1,2,2)
    shadedErrorBar(time,smooth(nanmean(WConditionsPre_Inc,2),7),semWPre_Con,'lineprops',{colors{i},'LineWidth',3}); 
    hold on; 
%     xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
    xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
%     if strcmp(timepoint,'all')
    xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
%     end
    xlabel('Time (s)')
    ylabel('Firing Rate (Hz)')
    axis tight
    ylim([0.5 1.1])
    legend({'Trials 1-5','Trials 6-10','Trials 11-15'})
    title('Incongruent')
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',25,'FontWeight','bold','LineWidth',5)




%%
titles = [{'Congruent P'} {'Incongruent P'} {'Congruent W'} {'Incongruent W'}];
pcSort = 1:3;

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

        imagesc([0 max(time)], [1 length(b)],centData(k,:),[-.5 .5]);

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
%     saveas(gca,[figPath filesep 'coefficienthistogram_split_scores_prePCAsplitGENO' timepoint num2str(pcSort) '_100ms'],'svg');     saveas(gca,[figPath filesep 'coefficienthistogram_split_scores_prePCAsplitGENO' timepoint num2str(pcSort) '_100ms'],'png')
% 
% %% Plot outputs of KSTEST2
% allP = [p1; p2; p3; p4; p5; p6];
% allH = [h1; h2; h3; h4; h5; h6];
% ylabels = [{'PCon-PInc'}; {'PCon-WCon'}; {'PCon-WInc'}; {'PInc-WCon'}; {'PInc-WInc'}; {'WCon-WInc'}];
% xlabels = [{'PC1'}; {'PC2'}; {'PC3'}; {'PC4'}; {'PC5'}];
% 
% figure('Units','normalized','Position',[0 0 1 1])
% imagesc(allH,[0 1]); colorbar; colormap jet;
% set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4,'XTick',1:6,'XTickLabel', xlabels, 'YTick',1:6,'YTickLabel',ylabels)
% saveas(gca,[figPath filesep 'coefficientKSTEST_HVALS_prePCAsplitGENO' timepoint num2str(pcSort) '_100ms'],'svg');     saveas(gca,[figPath filesep 'coefficientKSTEST_HVALS_prePCAsplitGENO' timepoint num2str(pcSort) '_100ms'],'png')
% 
% figure('Units','normalized','Position',[0 0 1 1])
% imagesc(allP,[0 1]); colorbar; colormap jet;
% set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4,'XTick',1:6,'XTickLabel', xlabels, 'YTick',1:6,'YTickLabel',ylabels)
% saveas(gca,[figPath filesep 'coefficientKSTEST_PVALS_prePCAsplitGENO' timepoint num2str(pcSort) '_100ms'],'svg');     saveas(gca,[figPath filesep 'coefficientKSTEST_PVALS_prePCAsplitGENO' timepoint num2str(pcSort) '_100ms'],'png')
%% Low-D PC Representation Based on Separated Matrices
% Parameters
time = (1:size(WConditionsPre,1))/10;

allPCProjections = [];
pcSort = 1:25;
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
        plot(time,allPCProjections{j}(:,pcSort(i)),'LineWidth',3,'Color',colors{i})
        if j == 1; legend([{'PC1'},{'PC2'},{'PC3'},{'PC4'},{'PC5'}],'AutoUpdate','off','location','best'); end
        hold on
        xlabel('Time (s)');
        title([titles{j}]);
        ylim([-10 20])
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

for i = 1:4
    plot(allPCProjections{i}(:,1),allPCProjections{i}(:,2),'Color',colors{i},'LineWidth',3);
    hold on
    plot(allPCProjections{i}(1,1),allPCProjections{i}(1,2),'o','MarkerSize',10,'MarkerFaceColor',colors{i},'MarkerEdgeColor','k')
    plot(allPCProjections{i}(cueOn*10,1),allPCProjections{i}(cueOn*10,2),'s','MarkerSize',12,'MarkerFaceColor',colors{i},'MarkerEdgeColor','k')
    plot(allPCProjections{i}(sipDescent*10,1),allPCProjections{i}(sipDescent*10,2),'d','MarkerSize',12,'MarkerFaceColor',colors{i},'MarkerEdgeColor','k')
    plot(allPCProjections{i}(sipAscent*10,1),allPCProjections{i}(sipAscent*10,2),'h','MarkerSize',12,'MarkerFaceColor',colors{i},'MarkerEdgeColor','k')

    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)            
end

xlabel('PC1')
ylabel('PC2')
saveas(gca,[figPath filesep '2DWS_split_pcSpaceProj_prePCAsplitGENO' timepoint '_100ms'],'svg');     saveas(gca,[figPath filesep '2DWS_split_pcSpaceProj_prePCAsplitGENO' timepoint '_100ms'],'png')

%% Calculate and Plot the Distance Between PCs 
stickPoint = [3 3];
pCon3D = allPCProjections{1}(:,1:stickPoint(2)); wCon3D = allPCProjections{3}(:,1:stickPoint(1));
pInc3D = allPCProjections{2}(:,1:stickPoint(2)); wInc3D = allPCProjections{4}(:,1:stickPoint(1));

pDistance = [diag(pdist2(pCon3D,pInc3D,'mahalanobis'))];
wDistance = [diag(pdist2(wCon3D,wInc3D,'mahalanobis'))];

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
ylabel('Distance (euclidean)')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,[figPath filesep 'allNeurons_split_3DpcSpaceDistance_prePCAsplitGENO' timepoint '_100ms'],'svg');     saveas(gca,[figPath filesep 'allNeurons_split_3DpcSpaceDistance_prePCAsplitGENO' timepoint '_100ms'],'png')


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

%%
indexedPosCoeff = []; indexedNegCoeff = [];
psthType = 'cueT';

numTrials = 15;
for j = 1:length(indexedCoeffData)
    % Assign variable rawData depending on condition 
    % Center rawData variable to then multiply by coeff to generate per
    % trial score.
    if j == 1
        rawData = pCon - nanmean(pCon,3);
        pcSort = 1:stickPoint(1);
    elseif j == 2
        rawData = pInc - nanmean(pInc,3);
        pcSort = 1:stickPoint(1);
    elseif j == 3
        rawData = wCon - nanmean(wCon,3);
        pcSort = 1:stickPoint(2);
    elseif j == 4
        rawData = wInc - nanmean(wInc,3);
        pcSort = 1:stickPoint(2);
    end
    % Sort according to top 3 PCs, plot similar to other coefficient
    % figures
%
%     figure('Units','normalized','Position',[0 0 1 1])
%     pcSort = 1:5;
    for i = 1:length(pcSort)
        % Index Pos and Neg
        PosCoeff = [];
        NegCoeff = [];
        for k = 1:numTrials
            PosCoeff(k,:,:) = indexedCoeffData{j}(indexedCoeffData{j}(:,i) >= 0.01,i)' .* squeeze(rawData(k,:,indexedCoeffData{j}(:,i) >= 0.01));                  
            NegCoeff(k,:,:) = indexedCoeffData{j}(indexedCoeffData{j}(:,i) <= -0.01,i)' .* squeeze(rawData(k,:,indexedCoeffData{j}(:,i) <= -0.01));
        end
        % Take nanmean
        
        PosCoeff = (nanmean(PosCoeff,3));
        NegCoeff = (nanmean(NegCoeff,3));
        % Index extra nanmean information
        mPosCoeff{j,i} = nanmean(PosCoeff,1);
        mNegCoeff{j,i} = nanmean(NegCoeff,1);
        stdPosCoeff{j,i} = std(PosCoeff,1)/sqrt(15);
        stdNegCoeff{j,i} = std(PosCoeff,1)/sqrt(15);
        % Index Pos Neg Info
        indexedPosCoeff{j,i} = PosCoeff;
        indexedNegCoeff{j,i} = NegCoeff;
% 
%         % Figure
%         subplot(2,max(pcSort),i);          % Doubling subplot for pos+neg
%         % Positive Plot
%         imagesc([0 max(time)], [1 15],PosCoeff, [0 0.05]);
%         colormap jet
%         colorbar
%         xlabel('Time (s)');
%         if i==1; ylabel('Trial #');end
%         xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
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
%         imagesc([0 max(time)], [1 15],abs(NegCoeff),[0 0.05]);
%         colormap jet
%         colorbar
%         xlabel('Time (s)');
%         if i==1; ylabel('Trial #');end
%         xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
%         xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
%         xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
% %         
%         if j < 3; title(['PC Num: ' num2str(pcSort(i)) ',ExplVar=' num2str(sum(pConditionsExplained(1:i)),'%.1f') '%']);
%         else
%             title(['PC Num: ' num2str(pcSort(i)) ',ExplVar=' num2str(sum(wConditionsExplained(1:i)),'%.1f') '%']);
%         end
%         set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)        
%     end
% 
%     sgtitle(titles{j})
%     set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
%     saveas(gca,[figPath filesep 'coeffNeurons_overTrials' titles{j} psthType timepoint '_100ms'],'svg');     saveas(gca,[figPath filesep 'coeffNeurons_overTrials' titles{j} psthType timepoint '_100ms'],'png')
    end

end


%% Calculate how each PC activation changes over trials
% Each PC has 15 trials
% Reduce those 15 trials to 5 trials made up of the averages of each of the
% 15 
colors = {[0 0 0], [0.5 0.5 0.5], [0.9 0.9 0.9]};
binMat = zeros(length(PosCoeff),15);
binMat(:,1:5) = 1;
binMat(:,5:10) = 2;
binMat(:,11:15) = 3;
idxN = [1:5; 6:10; 11:15];
figure('Units','normalized','Position',[0 0 1 1])
for j = 1:4
    subplot(2,2,j)
    for k = 1:min(size(idxN))
        
        plot3(smooth(zscore(nanmean(indexedPosCoeff{j,1}(idxN(k,:),:))),11), smooth(zscore(nanmean(indexedPosCoeff{j,2}(idxN(k,:),:))),11), smooth(zscore(nanmean(indexedPosCoeff{j,3}(idxN(k,:),:))),11),'LineWidth',3,'Color',colors{k});
        hold on

        x = smooth(zscore(nanmean(indexedPosCoeff{j,1}(idxN(k,:),:))),11); x0 = x(1); x1 = x(cueOn*10); x2 = x(sipDescent*10); x3 = x(sipAscent*10);
        y = smooth(zscore(nanmean(indexedPosCoeff{j,2}(idxN(k,:),:))),11); y0 = y(1); y1 = y(cueOn*10); y2 = y(sipDescent*10); y3 = y(sipAscent*10);
        z = smooth(zscore(nanmean(indexedPosCoeff{j,3}(idxN(k,:),:))),11); z0 = z(1); z1 = z(cueOn*10); z2 = z(sipDescent*10); z3 = z(sipAscent*10);
        xAll = [x0 x1 x2 x3];
        yAll = [y0 y1 y2 y3];
        zAll = [z0 z1 z2 z3];
        
        plot3(xAll,yAll,zAll,'ko','MarkerSize',10,'MarkerFaceColor','k');
        hold on
        xlabel('PC 1')
        ylabel('PC 2')
        zlabel('PC 3')
        title(titles{j})
        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
    end
        text(xAll,yAll,zAll, {'Start','Cue','Sip In','Sip Out'},'VerticalAlignment','top','HorizontalAlignment','left');

end

set(gcf,'Renderer','Painter')
saveas(gcf,'2D.svg');
%% Find distance between trial N and every other trial in 3D PC Space
method = 'euclidean';
for i = 1:4
    data = reshape([indexedPosCoeff{i,:}],[15,221,3]);
    for k = 1:length(data)
        distance{i}(:,:,k) = squareform(pdist(squeeze(data(:,k,:)),method));
    end
end

for i = 1:length(distance)
    for k = 1:length(data)
        for j = 1:min(size(idxN))
            dissimilaritynanmeans{i}(j,k) = nanmean(squeeze(distance{i}(idxN(j),:,k)));
        end
    end
end

%% Representative figures of dissimilarity metrics
figure('Units','normalized','Position',[0 0 1 1])
subplot(1,2,1)
imagesc(squeeze(distance{3}(:,:,cueOn*10)),[0 0.05]);
ylabel('Trial #')
xlabel('Trial #')
title('Congruent Wistars')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
colormap jet

subplot(1,2,2)
imagesc(squeeze(distance{4}(:,:,cueOn*10)),[0 0.05]);
ylabel('Trial #')
xlabel('Trial #')
title('Incongruent Wistars')
colorbar
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
colormap jet

%%
figure('Units','normalized','Position',[0 0 1 1])
for i = 1:length(distance)
    subplot(2,2,i)
    plot(time,(dissimilaritynanmeans{i}))
    xlabel('Time (s)')
    ylabel('Disimilarity')
end
%%
dataHolding = [];
figure('Units','normalized','Position',[0 0 1 1])
tP = ((cueOn-1)*10:(cueOn+1)*10);
for i = 1:length(distance)
    subplot(2,2,i)
    for j = 1:3

        bar(j,nanmean(nanmean(nanmean(distance{i}(idxN(j,:),:,tP),1),3)));
        dataHolding{i}{j} = nanmean(nanmean(distance{i}(idxN(j,:),:,tP),1),3);
        hold on
        er = errorbar(j,nanmean(nanmean(nanmean(distance{i}(idxN(j,:),:,tP),1),3)),std(nanmean(nanmean(distance{i}(idxN(j,:),:,tP),1),3))/sqrt(15));
        er.Color = [0 0 0];
        er.LineStyle = 'none';
        er.LineWidth = 3;
    end
    xlabel('Trial Blocks')
    ylabel('Disimilarity')
    xticks([1 2 3])
    xticklabels({'Trials 1-5', 'Trials 6-10', 'Trials 11-15'})
    ylim([0 0.04])
    title(titles{i})
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)

end
%% Testing different ways to find distance 
dataHolding = [];
figure('Units','normalized','Position',[0 0 1 1])
tP = ((cueOn-1)*10:(cueOn+1)*10);
for i = 1:length(distance)
    subplot(2,2,i)
    for j = 1:3

        bar(j,nanmean(nanmean(nanmean(distance{i}(idxN(j,:),:,tP),1),3)));
        dataHolding{i}{j} = nanmean(nanmean(distance{i}(idxN(j,:),:,tP),1),3);
        hold on
        er = errorbar(j,nanmean(nanmean(nanmean(distance{i}(idxN(j,:),:,tP),1),3)),std(nanmean(nanmean(distance{i}(idxN(j,:),:,tP),1),3))/sqrt(15));
        er.Color = [0 0 0];
        er.LineStyle = 'none';
        er.LineWidth = 3;
    end
    xlabel('Trial Blocks')
    ylabel('Disimilarity')
    xticks([1 2 3])
    xticklabels({'Trials 1-5', 'Trials 6-10', 'Trials 11-15'})
    ylim([0 0.04])
    title(titles{i})
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)

end
%% Create single bar graph for all values above
figure('Units','normalized','Position',[0 0 1 1])
hb = bar([nanmean(dataHolding{3}{1}) nanmean(dataHolding{4}{1}) nanmean(dataHolding{1}{1}) nanmean(dataHolding{2}{1}); ... 
    nanmean(dataHolding{3}{2}) nanmean(dataHolding{4}{2}) nanmean(dataHolding{1}{2}) nanmean(dataHolding{2}{2}); ...
    nanmean(dataHolding{3}{3}) nanmean(dataHolding{4}{3}) nanmean(dataHolding{1}{3}) nanmean(dataHolding{2}{3})], 'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);

offsetPos = [];
scatterXData = [];

for i = 1:length(hb)
    offsetPos(i,:) = hb(i).XData + hb(i).XOffset;
end
for i = 1:length(offsetPos)
    if i == 1
        dataLength = 15;
    elseif i == 2
        dataLength = 15;
    elseif i == 3
        dataLength = 15;
    elseif i == 4
        dataLength = 15;
    end
    for j = 1:min(size(offsetPos))
        scatterXData{i,j} = linspace(offsetPos(i,j)-0.01,offsetPos(i,j)+0.01,dataLength);
    end
end
WConData = [dataHolding{3}{1}; dataHolding{3}{2}; dataHolding{3}{3}];
WIncData = [dataHolding{4}{1}; dataHolding{4}{2}; dataHolding{4}{3}];

PConData = [dataHolding{1}{1}; dataHolding{1}{2}; dataHolding{1}{3}];
PIncData = [dataHolding{2}{1}; dataHolding{2}{2}; dataHolding{2}{3}];
hold on
allData = {WConData, WIncData, PConData, PIncData};
for i = 1:length(offsetPos)
    plot([scatterXData{i,:}],allData{i}(:),'ko','LineWidth',2,'MarkerSize',8)

    for j = 1:min(size(offsetPos))
        er = errorbar(offsetPos(i,j),nanmean(allData{i}(j,:)),std(allData{i}(j,:))/sqrt(length(allData{i}(j,:))));
        er.Color = [0 0 0];
        er.LineStyle = 'none';
    end
end
ylim([0 0.05])
ylabel('Dissimilarity')
xticklabels({'Trials 1-5', 'Trials 6-10','Trials 11-15'})

set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
figPath = 'F:\dissDat\revisionfigs';
saveas(gca,[figPath filesep 'dissimilarity_cueon'     '_100ms' ],'svg');     saveas(gca,[figPath filesep 'dissimilarity_cueon '    '_100ms' ],'png')
%%
% Stats
formatData = cell2mat(allData);
formatData = formatData';
datatable = array2table(formatData,'VariableNames',{'Block1','Block2','Block3'});
datatable.Strain = [repmat({'Wistar'},30,1); repmat({'P rat'},30,1)];
condition = [repmat({'Congruent'},15,2); repmat({'Incongruent'},15,2)];
condition = {(condition{:})};
datatable.Condition = condition';

% Create rm model
Block = table([1 2 3]','VariableNames',{'Blocks'});
rm = fitrm(datatable,'Block1-Block3 ~ Strain * Condition','WithinDesign',Block);
anovaTbl = ranova(rm,'WithinModel','Blocks')
multTbl1 = multcompare(rm,'Condition','By','Strain','By','Blocks')
multTbl1 = multcompare(rm,'Condition','By','Blocks')
multTbl1 = multcompare(rm,'Blocks','By','Strain','By','Condition')
multTbl2 = multcompare(rm,'Blocks','By','Strain','By','Condition')

tblStr = formattedDisplayText(anovaTbl); 
% Write string to file
fid = fopen('dissimilarity_cueon_ranova.txt', 'wt');
fileCleanup = onCleanup(@()fclose(fid));
formatSpec = '%s\n';
fprintf(fid, formatSpec, tblStr);
clear('fileCleanup')

tblStr = formattedDisplayText(multTbl2); 
% Write string to file
fid = fopen([coeffValue 'dissimilarity_cueon.txt'], 'wt');
fileCleanup = onCleanup(@()fclose(fid));
formatSpec = '%s\n';
fprintf(fid, formatSpec, tblStr);
clear('fileCleanup')

filename = 'F:\dissDat\csvs\datatable_dissimilarity_cueon';
writetable(datatable,filename);
%%

bin_nanmeans = [];
for j = 1:4
    for i = 1:stickPoint
        data = indexedPosCoeff{j,i}';
        for n = 1:3
            for k = 1:length(data)
                bin_nanmeans{j,i}(k,n) = nanmean(data(k,binMat(k,:)==n));
            end
        end
    end
end

%%
titles = [{'Congruent P'} {'Incongruent P'} {'Congruent W'} {'Incongruent W'}];

figure('Units','normalized','Position',[0 0 1 1])

for i = 1:length(bin_nanmeans)
    subplot(2,2,i)
    for k = 1:min(size(bin_nanmeans))
        plot(time,zscore(bin_nanmeans{i,k}))
    end
    title(titles{i})
    xlabel('Time (s)')
    ylabel('PC Score (Z-Scored)')
    
end


%% Calculate max val per epoch per trial similar to behavior
% Time periods to index epochs with
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
coeffValue = 'posCoeff';

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
        
            pdistMaxnanmean(j) = nanmean(max(pDist(:,tp),[],2));
            wdistMaxnanmean(j) = nanmean(max(wDist(:,tp),[],2));

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
    
        pdistMaxnanmean(j) = nanmean(max(pDist(:,tp),[],2));
        wdistMaxnanmean(j) = nanmean(max(wDist(:,tp),[],2));

        pdm(:,j) = max(pDist(:,tp),[],2);
        wdm(:,j) = max(wDist(:,tp),[],2);
    
        pdistSTD(j) = std(max(pDist(:,tp),[],2))/sqrt(numTrials);
        wdistSTD(j) = std(max(wDist(:,tp),[],2))/sqrt(numTrials); 
    end
end

plot(pdistMaxnanmean,'k-o','LineWidth',3)
hold on
plot(wdistMaxnanmean,'b-o','LineWidth',3)
xticks(1:1:length(epochLabels))
xticklabels(epochLabels)
ylabel(['Distance (' method ')'])
xlim([0.9 numEpochs+0.1])
er = errorbar(pdistMaxnanmean, pdistSTD);
er.Color = [0 0 0];
er.LineStyle = 'none';
er.LineWidth = 3;

er = errorbar(wdistMaxnanmean, wdistSTD);
er.Color = [0 0 0];
er.LineStyle = 'none';
er.LineWidth = 3;
% ylim([0 10])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
saveas(gca,[figPath filesep 'Coeff_Distance_GenotypesComparison' coeffValue psthType '_' method timepoint '_100ms'],'svg');     saveas(gca,[figPath filesep 'Coeff_Distance_GenotypesComparison' psthType '_' method coeffValue timepoint '_100ms'],'png')

% Stats
trialRep = repmat([1:15],1,2); trialRep = trialRep';
tableData = [wdm; pdm];
maxTable = array2table(tableData,'VariableNames',{'PreCue','CueOn','EarlyEtOHAccess', 'LateEtOHAccess','PostAccess'});
maxTable.Strain = [repmat({'Wistar'},15,1); repmat({'P rat'},15,1)];
maxTable.Trial = trialRep;

% Create rm model
Epoch = table([1 2 3 4 5]','VariableNames',{'Epochs'});
rm = fitrm(maxTable,'PreCue-PostAccess ~ Strain','WithinDesign',Epoch);
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

mpDist = nanmean(pDist); sempDist = std(pDist)/sqrt(numTrials);
mwDist = nanmean(wDist); semwDist = std(wDist)/sqrt(numTrials);

figure('Units','normalized','Position',[0 0 1 1])
shadedErrorBar(time,mpDist,sempDist,'lineprops',{'r','LineWidth',3})
hold on
shadedErrorBar(time,mwDist,semwDist,'lineprops',{'b','LineWidth',3})    
xlabel('Time (s)')
ylabel(['Distance (' method ')'])
ylim([5 20])

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

rm = fitrm(statTable,'statTable1-statTable251 ~ Genotype','WithinDesign',Time);
anovaTbl = ranova(rm,'WithinModel','Time')


tblStr = formattedDisplayText(anovaTbl); 
fid = fopen([ coeffValue 'cueDistance_RANOVA.txt'], 'wt');
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
ylim([5 20])
xlim([0 22.1])
if strcmp(psthType,'cueT')
    xline(cueOn,'k--','Cue On','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');    
    xline(sipDescent,'k--','Sipper In','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
    xline(sipAscent,'k--','Sipper Out','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');
elseif strcmp(psthType,'approachT')
    xline(cueOn,'k--','Approach','LineWidth',3,'FontSize',10,'FontName','Arial','FontWeight','bold');    
end
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

plot(time,fdr_p_plot+4.5, 'k*','MarkerSize',6)
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,[figPath filesep 'PostrlAvgDistances_withmarkers' psthType coeffValue timepoint '_100ms'],'svg');     saveas(gca,[figPath filesep 'PostrlAvgDistances_withmarkers' psthType coeffValue timepoint '_100ms'],'png')


%% Assess neural correlations following task events over trials
% Create PSTHs for W's P's
% Within wMat and pMat there are several neurons but they are not
% delineated by the condition they are (congruent or incongruent). Adding
% an additional dimension to each neuron to signify their allegiance will
% help. From there, we can assess how each of these neurons align to PCs
% associated with session event like the cue and take their nanmean, first,
% and then correlation coefficent over trials.
    
for i = 1:length(trlStruct)
    if ismember(i,find(regPIdx == 1)) || ismember(i,find(regWIdx == 1))
        ephysStruct(i).sessionLabel = ones(1,ephysStruct(i).neuronNum);
    else
        ephysStruct(i).sessionLabel = zeros(1,ephysStruct(i).neuronNum);
    end
end
%         for i = 1:5
%             firstPC_P{i} = pConditionsCoeff(:,i) > 0.01;
%             firstPC_W{i} = wConditionsCoeff(:,i) > 0.01;
%         end

firstPC_P = pConditionsCoeff(:,i) > 0.01;
firstPC_W = wConditionsCoeff(:,i) > 0.01;

pMatCon = squeeze(pMat(:,10,:));% - pMat(:,8,:)); % Take first time point as a control 
pMatCue = squeeze(pMat(:,cOn,:));% - pMat(:,cOn-2,:));
pMatDes = squeeze(pMat(:,sipD,:));% - pMat(:,sipD-2,:));
pMatAes = squeeze(pMat(:,sipA,:));% - pMat(:,sipA-2,:));

pMatCueSession = [ephysStruct(PIdx).sessionLabel];

wMatCon = squeeze(wMat(:,10,:));% - wMat(:,8,:)); 
wMatCue = squeeze(wMat(:,cOn,:));%- wMat(:,cOn-2,:));% - wMat(:,cueOn*10+5));
wMatDes = squeeze(wMat(:,sipD,:));% - wMat(:,sipD-2,:));
wMatAes = squeeze(wMat(:,sipA,:));% - wMat(:,sipA-2,:));

wMatCueSession = [ephysStruct(WIdx).sessionLabel];



WconCon = wMatCon(:,wMatCueSession' == 1 & firstPC_W == 1);
WcueCon = wMatCue(:,wMatCueSession' == 1 & firstPC_W == 1);
WDesCon = wMatDes(:,wMatCueSession' == 1 & firstPC_W == 1);
WAesCon = wMatAes(:,wMatCueSession' == 1 & firstPC_W == 1); 

WconInc = wMatCon(:,wMatCueSession' == 0 & firstPC_W == 1);
WcueInc = wMatCue(:,wMatCueSession' == 0 & firstPC_W == 1);
WDesInc = wMatDes(:,wMatCueSession' == 0 & firstPC_W == 1);
WAesInc = wMatAes(:,wMatCueSession' == 0 & firstPC_W == 1);

PconCon = pMatCon(:,pMatCueSession' == 1 & firstPC_P == 1);
PcueCon = pMatCue(:,pMatCueSession' == 1 & firstPC_P == 1);
PDesCon = pMatDes(:,pMatCueSession' == 1 & firstPC_P == 1);
PAesCon = pMatAes(:,pMatCueSession' == 1 & firstPC_P == 1);

PconInc = pMatCon(:,pMatCueSession' == 0 & firstPC_P == 1);
PcueInc = pMatCue(:,pMatCueSession' == 0 & firstPC_P == 1);
PDesInc = pMatDes(:,pMatCueSession' == 0 & firstPC_P == 1);
PAesInc = pMatAes(:,pMatCueSession' == 0 & firstPC_P == 1);

%% Calculate mutual info
for i = 1:15
    WconCueMut(i) = mutualinfo(WcueCon(i,:),WconCon(i,:));
    WconDesMut(i) = mutualinfo(WDesCon(i,:),WconCon(i,:));
    WconAesMut(i) = mutualinfo(WAesCon(i,:),WconCon(i,:));

    WincCueMut(i) = mutualinfo(WcueInc(i,:),WconInc(i,:));
    WincDesMut(i) = mutualinfo(WDesInc(i,:),WconInc(i,:));
    WincAesMut(i) = mutualinfo(WAesInc(i,:),WconInc(i,:));

    PconCueMut(i) = mutualinfo(PcueCon(i,:),PconCon(i,:));
    PconDesMut(i) = mutualinfo(PDesCon(i,:),PconCon(i,:));
    PconAesMut(i) = mutualinfo(PAesCon(i,:),PconCon(i,:));

    PincCueMut(i) = mutualinfo(PcueInc(i,:),PconInc(i,:));
    PincDesMut(i) = mutualinfo(PDesInc(i,:),PconInc(i,:));
    PincAesMut(i) = mutualinfo(PAesInc(i,:),PconInc(i,:));
end


%%
WCorrCon_con = [];
WCorrInc_con = [];
PCorrCon_con = [];
PCorrInc_con = [];
WCorrCon_cue = [];
WCorrInc_cue = [];
PCorrCon_cue = [];
PCorrInc_cue = [];
WCorrCon_des = [];
WCorrInc_des = [];
PCorrCon_des = [];
PCorrInc_des = [];
WCorrCon_aes = [];
WCorrInc_aes = [];
PCorrCon_aes = [];
PCorrInc_aes = [];
k = 1;
for i = 1:15
    WCorrCon_con(:,:,k) = corrcoef(WconCon(i,:), WconCon(i,:));
    WCorrInc_con(:,:,k) = corrcoef(WconInc(i,:), WconInc(i,:));
    PCorrCon_con(:,:,k) = corrcoef(PconCon(i,:), PconCon(i,:));
    PCorrInc_con(:,:,k) = corrcoef(PconInc(i,:), PconInc(i,:));

    WCorrCon_cue(:,:,k) = corrcoef(WconCon(i,:), WcueCon(i,:));
    WCorrInc_cue(:,:,k) = corrcoef(WconInc(i,:), WcueInc(i,:));
    PCorrCon_cue(:,:,k) = corrcoef(PconCon(i,:), PcueCon(i,:));
    PCorrInc_cue(:,:,k) = corrcoef(PconInc(i,:), PcueInc(i,:));

    WCorrCon_des(:,:,k) = corrcoef(WconCon(i,:), WDesCon(i,:));
    WCorrInc_des(:,:,k) = corrcoef(WconInc(i,:), WDesInc(i,:));
    PCorrCon_des(:,:,k) = corrcoef(PconCon(i,:), PDesCon(i,:));
    PCorrInc_des(:,:,k) = corrcoef(PconInc(i,:), PDesInc(i,:));

    WCorrCon_aes(:,:,k) = corrcoef(WconCon(i,:), WAesCon(i,:));
    WCorrInc_aes(:,:,k) = corrcoef(WconInc(i,:), WAesInc(i,:));
    PCorrCon_aes(:,:,k) = corrcoef(PconCon(i,:), PAesCon(i,:));
    PCorrInc_aes(:,:,k) = corrcoef(PconInc(i,:), PAesInc(i,:));
    k = 1+k;
end
%%
WCorrCon_con = squeeze(WCorrCon_con(1,2,:));
WCorrInc_con = squeeze(WCorrInc_con(1,2,:));
PCorrCon_con = squeeze(PCorrCon_con(1,2,:));
PCorrInc_con = squeeze(PCorrInc_con(1,2,:));

WCorrCon_cue = squeeze(WCorrCon_cue(1,2,:));
WCorrInc_cue = squeeze(WCorrInc_cue(1,2,:));
PCorrCon_cue = squeeze(PCorrCon_cue(1,2,:));
PCorrInc_cue = squeeze(PCorrInc_cue(1,2,:));

WCorrCon_des = squeeze(WCorrCon_des(1,2,:));
WCorrInc_des = squeeze(WCorrInc_des(1,2,:));
PCorrCon_des = squeeze(PCorrCon_des(1,2,:));
PCorrInc_des = squeeze(PCorrInc_des(1,2,:));

WCorrCon_aes = squeeze(WCorrCon_aes(1,2,:));
WCorrInc_aes = squeeze(WCorrInc_aes(1,2,:));
PCorrCon_aes = squeeze(PCorrCon_aes(1,2,:));
PCorrInc_aes = squeeze(PCorrInc_aes(1,2,:));

ind = [1 1 1 2 2 2 3 3 3 4 4 4 5 5 5]; % Grouping index
ind = ind';

%%
figure

subplot(2,2,1)
plot(groupsummary(WCorrCon_con,ind,'nanmean'),'-o')
hold on
plot(groupsummary(WCorrInc_con,ind,'nanmean'),'-o')
plot(groupsummary(PCorrCon_con,ind,'nanmean'),'-o')
plot(groupsummary(PCorrInc_con,ind,'nanmean'),'-o')
xlabel('Trials')
ylabel('nanmean Correlation Coefficient')
title('Control (Pre-Cue)')
ylim([0 1])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

subplot(2,2,2)
plot(groupsummary(WCorrCon_cue,ind,'nanmean'),'-o')
hold on
plot(groupsummary(WCorrInc_cue,ind,'nanmean'),'-o')
plot(groupsummary(PCorrCon_cue,ind,'nanmean'),'-o')
plot(groupsummary(PCorrInc_cue,ind,'nanmean'),'-o')
xlabel('Trials')
ylabel('nanmean Correlation Coefficient')
title('Cue On')
ylim([0 1])

set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

subplot(2,2,3)
plot(groupsummary(WCorrCon_des,ind,'nanmean'),'-o')
hold on
plot(groupsummary(WCorrInc_des,ind,'nanmean'),'-o')
plot(groupsummary(PCorrCon_des,ind,'nanmean'),'-o')
plot(groupsummary(PCorrInc_des,ind,'nanmean'),'-o')
xlabel('Trials')
ylabel('nanmean Correlation Coefficient')
title('Sipper In')
ylim([0 1])

set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

subplot(2,2,4)
plot(groupsummary(WCorrCon_aes,ind,'nanmean'),'-o')
hold on
plot(groupsummary(WCorrInc_aes,ind,'nanmean'),'-o')
plot(groupsummary(PCorrCon_aes,ind,'nanmean'),'-o')
plot(groupsummary(PCorrInc_aes,ind,'nanmean'),'-o')
xlabel('Trials')
ylabel('nanmean Correlation Coefficient')
title('Sipper Out')
ylim([0 1])

set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
% 
% WCorrnanmeanCor = nannanmean(WCorrCon);
% WcorrnanmeanInc
% 
% PCorrnanmeanCor
% PCorrnanmeanInc