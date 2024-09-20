%% 2CAP Cognitive Control Study
% 2 Important Revisions
% Calculate behavior in the early vs. late sessions (argument against
% intoxication).
% Calculate distance from sipper of each strain at different epochs
% pre-cue, cue-on, early-etoh, late-etoh, post-etoh.
%% Load and check data
clear all; close all % Refresh workspace
% PATHS
% For brain3 Linux Machine:
if ispc
    parentPath = 'F:/dissDat';
    figPath = 'F:/dissDat/revisionfigs';
else
    parentPath = '/research/dissDat';
    figPath = '/research/dissDat/figs';
end
% LOAD
load([parentPath filesep 'ephysStruct.mat']);
load([parentPath filesep 'trlStruct.mat']);
load([parentPath filesep 'masterTable.mat']);
addpath(genpath('F:\dissDat\restoredScripts'))
addpath(genpath([parentPath filesep 'analysisScripts']))
%% 
% First, calculate time from 756
time = (1:751) / 30; % Behavior is at 30 FPS
binVal = [1:5; 6:10; 11:15];
sipDescent = 10*30;
sipAscent = 18*30;
cueOn = 5*30;

% Set indexing 
regPIdx = startsWith(masterTbl.SessionType,'Regular') & startsWith(masterTbl.Strain,'P');
revPIdx = startsWith(masterTbl.SessionType,'Reversal') & startsWith(masterTbl.Strain,'P');
    
regWIdx = startsWith(masterTbl.SessionType,'Regular') & startsWith(masterTbl.Strain,'W');
revWIdx = startsWith(masterTbl.SessionType,'Reversal') & startsWith(masterTbl.Strain,'W');
%%
for i = 1:length(trlStruct)
    trialTimes = trlStruct(i).trialTimes(1:48);
    [trialTimes,trialTimesIdx] = sort(trialTimes);
    LDist = trlStruct(i).LcorrectApproachTrls;
    RDist = trlStruct(i).RcorrectApproachTrls;
    allCorrectAppDistance = [LDist; RDist];
    allCorrectAppDistance = allCorrectAppDistance(trialTimesIdx,:);


    correctChoices = [trlStruct(i).LcorrectIdx; trlStruct(i).RcorrectIdx];
    correctChoices = correctChoices(trialTimesIdx);

    allCorrectAppDistance = allCorrectAppDistance(correctChoices,:);

    trlStruct(i).allCorrectAppDistance = allCorrectAppDistance;
    

end

%%
for i = 1:length(trlStruct)
    trialTimes = trlStruct(i).trialTimes(1:48);
    [trialTimes,trialTimesIdx] = sort(trialTimes);
%   PIdx = startsWith(masterTbl.Strain,'P');
%   WIdx = startsWith(masterTbl.Strain,'W');

    if regPIdx(i) == 1 || regWIdx(i) == 1
        RDist = trlStruct(i).trlRSipDist(1:24,:);
        LDist = trlStruct(i).trlLSipDist(25:48,:);
    elseif revPIdx(i) == 1 || revWIdx(i) == 1
        RDist = trlStruct(i).trlRSipDist(25:48,:);
        LDist = trlStruct(i).trlLSipDist(1:24,:);        
    end


    allCorrectAppDistance = [LDist; RDist];
    allCorrectAppDistance = allCorrectAppDistance(trialTimesIdx,:);


    incorrectChoices = [trlStruct(i).LincorrectIdx; trlStruct(i).RincorrectIdx];
    incorrectChoices = incorrectChoices(trialTimesIdx);

    allInCorrectAppDistance = allCorrectAppDistance(incorrectChoices,:);

    trlStruct(i).allInCorrectAppDistance = allInCorrectAppDistance;
    
end
%% Pull data, save trial structure if possible

WCon = {trlStruct(regWIdx).allCorrectAppDistance};
WInc = {trlStruct(revWIdx).allCorrectAppDistance};

PCon = {trlStruct(regPIdx).allCorrectAppDistance};
PInc = {trlStruct(revPIdx).allCorrectAppDistance};

% NanPad 
% Not all sessions have 15 trials. This is rough! Therefore, we are padding
% the trials with NaNs to circumvent the headache this causes.

endPad = 15;
for i = 1:length(WCon)
    sizeIdxWC = size(WCon{i});
    WCon{i}(sizeIdxWC(1):endPad,:) = NaN;

    sizeIdxWIC = size(WInc{i});
    WInc{i}(sizeIdxWIC(1):endPad,:) = NaN;
end
for i = 1:length(PCon)
    sizeIdxPC = size(PCon{i});
    PCon{i}(sizeIdxPC(1):endPad,:) = NaN;

    sizeIdxPIC = size(PInc{i});
    PInc{i}(sizeIdxPIC(1):endPad,:) = NaN;
end
WCon_All = []; WInc_All = []; PCon_All = []; PInc_All = [];
for i = 1:length(WCon)
    for j = 1:min(size(binVal))
        WCon_All{i}(:,:,j) = cat(1,WCon{i}(binVal(j,:),:));
        WInc_All{i}(:,:,j) = cat(1,WInc{i}(binVal(j,:),:));
    end
end

for i = 1:length(PCon)
    for j = 1:min(size(binVal))
        PCon_All{i}(:,:,j) = PCon{i}(binVal(j,:),:);
        PInc_All{i}(:,:,j) = PInc{i}(binVal(j,:),:);
    end
end

WCon_All = cat(1,WCon_All{:});
WInc_All = cat(1,WInc_All{:});

PCon_All = cat(1,PCon_All{:});
PInc_All = cat(1,PInc_All{:});

%% Plot traces of each category vis a vis subplot 
% First group data by subject
% Grouping Var for Wistars
groupingVar = ones(1,min(size(WCon_All)));
idxVar = 1:5:65;
for i = 1:length(idxVar)
    if i < 13
        groupingVar2(idxVar(i):idxVar(i+1) - 1,:) = i;
    elseif i == 13
        groupingVar2(idxVar(i):idxVar(i) + 4,:) = i;
    end
end

for i = 1:3
    sqz_WCon = squeeze(WCon_All(:,:,i));
    WCon_All_GRP(:,:,i) = grpstats(sqz_WCon,groupingVar2);
    sqz_WInc = squeeze(WInc_All(:,:,i));
    WInc_All_GRP(:,:,i) = grpstats(sqz_WInc,groupingVar2);
end
% P rats
groupingVar = ones(1,min(size(PCon_All)));
idxVar = 1:5:70;
for i = 1:length(idxVar)
    if i < 14
        groupingVar2(idxVar(i):idxVar(i+1) - 1,:) = i;
    elseif i == 14
        groupingVar2(idxVar(i):idxVar(i) + 4,:) = i;
    end
end
%
for i = 1:3
    sqz_PCon = squeeze(PCon_All(:,:,i));
    PCon_All_GRP(:,:,i) = grpstats(sqz_PCon,groupingVar2);
    sqz_PInc = squeeze(PInc_All(:,:,i));
    PInc_All_GRP(:,:,i) = grpstats(sqz_PInc,groupingVar2);
end
%%
for i = 1:3
    WCon_Trl(i,:) = nanmean(squeeze(WCon_All_GRP(:,:,i)));
    WInc_Trl(i,:) = nanmean(squeeze(WInc_All_GRP(:,:,i)));
    PCon_Trl(i,:) = nanmean(squeeze(PCon_All_GRP(:,:,i)));
    PInc_Trl(i,:) = nanmean(squeeze(PInc_All_GRP(:,:,i)));

    WCon_TrlSEM(i,:) = nanstd(squeeze(WCon_All_GRP(:,:,i))/sqrt(13));
    WInc_TrlSEM(i,:) = nanstd(squeeze(WInc_All_GRP(:,:,i))/sqrt(13));
    PCon_TrlSEM(i,:) = nanstd(squeeze(PCon_All_GRP(:,:,i))/sqrt(14));
    PInc_TrlSEM(i,:) = nanstd(squeeze(PInc_All_GRP(:,:,i))/sqrt(14));
end
LineStyles = {'b-','r-','g-'};
figure('Units','normalized','Position',[0 0 1 1])
subplot(2,2,1)
for i = 1:3
    shadedErrorBar(time,WCon_Trl(i,:),WCon_TrlSEM(i,:),'lineprops',{LineStyles{i},'LineWidth',3}); 
end
ylim([0 300])
xline(cueOn/30,'k--','Cue On','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(sipAscent/30,'k--','SipOut','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(sipDescent/30,'k--','SipIn','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
title('Wistar Congruent')
xlabel('Time (s)')
ylabel('Distance (pixels)')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',25,'FontWeight','bold','LineWidth',5)

subplot(2,2,2)
for i = 1:3
    shadedErrorBar(time,WInc_Trl(i,:),WInc_TrlSEM(i,:),'lineprops',{LineStyles{i},'LineWidth',3}); 
end
ylim([0 300])
xline(cueOn/30,'k--','Cue On','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(sipAscent/30,'k--','SipOut','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(sipDescent/30,'k--','SipIn','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
title('Wistar Incongruent')
xlabel('Time (s)')
ylabel('Distance (pixels)')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',25,'FontWeight','bold','LineWidth',5)

subplot(2,2,3)
for i = 1:3
    shadedErrorBar(time,PCon_Trl(i,:),PCon_TrlSEM(i,:),'lineprops',{LineStyles{i},'LineWidth',3}); 
end
ylim([0 300])
xline(cueOn/30,'k--','Cue On','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(sipAscent/30,'k--','SipOut','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(sipDescent/30,'k--','SipIn','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
title('P rat Congruent')
xlabel('Time (s)')
ylabel('Distance (pixels)')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',25,'FontWeight','bold','LineWidth',5)

subplot(2,2,4)
for i = 1:3
    shadedErrorBar(time,PInc_Trl(i,:),PInc_TrlSEM(i,:),'lineprops',{LineStyles{i},'LineWidth',3}); 
end
ylim([0 300])
xline(cueOn/30,'k--','Cue On','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(sipAscent/30,'k--','SipOut','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(sipDescent/30,'k--','SipIn','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
title('P rat Incongruent')
xlabel('Time (s)')
ylabel('Distance (pixels)')
legend({'Trials 1-5', 'Trials 6-10', 'Trials 11-15'})
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',25,'FontWeight','bold','LineWidth',5)

saveas(gca,[figPath filesep 'distancetocorrectsipper_correctTrials_4x4_binned15'     '_100ms' ],'svg');     saveas(gca,[figPath filesep 'distancetocorrectsipper_correctTrials_4x4_binned15'    '_100ms' ],'png')
%%

%% Measure of modulation of each epoch (incorrect)

% Find the difference between each time point, divide by the lenght of time
% between those timepoints to obtain pixels per 1/30s. Multipy by 30 to
% obtain velocity in pixel /s. 
for i = 1:3
    preCueWCon3(:,:,i) = abs((WCon_All_GRP(:,1,i) - WCon_All_GRP(:,cueOn,i)) / length(1:cueOn)) * 30;
    preCueWInc3(:,:,i) = abs((WInc_All_GRP(:,1,i) - WInc_All_GRP(:,cueOn,i)) / length(1:cueOn)) * 30;
    preCuePCon3(:,:,i) = abs((PCon_All_GRP(:,1,i) - PCon_All_GRP(:,cueOn,i)) / length(1:cueOn)) * 30;
    preCuePInc3(:,:,i) = abs((PInc_All_GRP(:,1,i) - PInc_All_GRP(:,cueOn,i)) / length(1:cueOn)) * 30;
    
    cueOnWCon3(:,:,i) = abs((WCon_All_GRP(:,cueOn,i) - WCon_All_GRP(:,sipDescent,i)) / length(cueOn:sipDescent)) * 30;
    cueOnWInc3(:,:,i) = abs((WInc_All_GRP(:,cueOn,i) - WInc_All_GRP(:,sipDescent,i)) / length(cueOn:sipDescent)) * 30;
    cueOnPCon3(:,:,i) = abs((PCon_All_GRP(:,cueOn,i) - PCon_All_GRP(:,sipDescent,i)) / length(cueOn:sipDescent)) * 30;
    cueOnPInc3(:,:,i) = abs((PInc_All_GRP(:,cueOn,i) - PInc_All_GRP(:,sipDescent,i)) / length(cueOn:sipDescent)) * 30;
    
    sipInWCon3(:,:,i) = abs((WCon_All_GRP(:,sipDescent,i) - WCon_All_GRP(:,(sipAscent),i)) / length(sipDescent:sipAscent)) * 30;
    sipInWInc3(:,:,i) = abs((WInc_All_GRP(:,sipDescent,i) - WInc_All_GRP(:,(sipAscent),i)) / length(sipDescent:sipAscent)) * 30;
    sipInPCon3(:,:,i) = abs((PCon_All_GRP(:,sipDescent,i) - PCon_All_GRP(:,(sipAscent),i)) / length(sipDescent:sipAscent)) * 30;
    sipInPInc3(:,:,i) = abs((PInc_All_GRP(:,sipDescent,i) - PInc_All_GRP(:,(sipAscent),i)) / length(sipDescent:sipAscent)) * 30;
end
%% Grab data and pull into single variable for group stats 
allCueData = {squeeze(cueOnWCon3),squeeze(cueOnWInc3),squeeze(cueOnPCon3),squeeze(cueOnPInc3)};

figure('Units','normalized','Position',[0 0 1 1])
hb = bar([nanmean(allCueData{1}(:,1)) nanmean(allCueData{2}(:,1)) nanmean(allCueData{3}(:,1)) nanmean(allCueData{4}(:,1)); nanmean(allCueData{1}(:,2)) nanmean(allCueData{2}(:,2)) nanmean(allCueData{3}(:,2)) nanmean(allCueData{4}(:,2)); nanmean(allCueData{1}(:,3)) nanmean(allCueData{2}(:,3)) nanmean(allCueData{3}(:,3)) nanmean(allCueData{4}(:,3))], 'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);

offsetPos = [];
scatterXData = [];


for i = 1:length(hb)
    offsetPos(i,:) = hb(i).XData + hb(i).XOffset;
end
for i = 1:length(offsetPos)
    if i == 1
        dataLength = size(WCon,1);
    elseif i == 2
        dataLength = size(WInc,1);
    elseif i == 3
        dataLength = size(PCon,1);
    elseif i == 4
        dataLength = size(PInc,1);
    end
    for j = 1:min(size(offsetPos))
        scatterXData{i,j} = linspace(offsetPos(i,j)-0.01,offsetPos(i,j)+0.01,dataLength);
    end
end
hold on

for i = 1:length(offsetPos)
    plot([scatterXData{i,:}],allCueData{i},'ko','LineWidth',2,'MarkerSize',8)

    for j = 1:min(size(offsetPos))
        er = errorbar(offsetPos(i,j),nanmean(allCueData{i}(:,j)),nanstd(allCueData{i}(:,j))/sqrt(length(allCueData{i}(:,j))));
        er.Color = [0 0 0];
        er.LineStyle = 'none';
        er.LineWidth = 3;
    end
end

ylabel('Speed (Pixels / Second)')
ylim([0 100])
xticklabels({'Trials 1-5', 'Trials 6-10','Trials 11-15'})

set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)

saveas(gca,[figPath filesep 'speedtocorrectsipper_correctTrials_15Bin'     '_100ms' ],'svg');     saveas(gca,[figPath filesep 'speedtocorrectsipper_correctTrials_15Bin'    '_100ms' ],'png')


%%
% figure('Units','normalized','Position',[0 0 1 1])
% hb = bar([mean(preCueWCon) mean(preCueWInc) mean(preCuePCon) mean(preCuePInc); mean(cueOnWCon) mean(cueOnWInc) mean(cueOnPCon) mean(cueOnPInc); mean(sipInWCon) mean(sipInWInc) mean(sipInPCon) mean(sipInPInc)], 'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);
% 
% offsetPos = [];
% scatterXData = [];
% 
% for i = 1:length(hb)
%     offsetPos(i,:) = hb(i).XData + hb(i).XOffset;
% end
% for i = 1:length(offsetPos)
%     if i == 1
%         dataLength = size(WCon,1);
%     elseif i == 2
%         dataLength = size(WInc,1);
%     elseif i == 3
%         dataLength = size(PCon,1);
%     elseif i == 4
%         dataLength = size(PInc,1);
%     end
%     for j = 1:min(size(offsetPos))
%         scatterXData{i,j} = linspace(offsetPos(i,j)-0.01,offsetPos(i,j)+0.01,dataLength);
%     end
% end
% WConData = [preCueWCon cueOnWCon sipInWCon];
% WIncData = [preCueWInc cueOnWInc sipInWInc];
% 
% PConData = [preCuePCon cueOnPCon sipInPCon];
% PIncData = [preCuePInc cueOnPInc sipInPInc];
% hold on
% allData = {WConData, WIncData, PConData, PIncData};
% 
% for i = 1:length(offsetPos)
%     plot([scatterXData{i,:}],allData{i}(:),'ko','LineWidth',2,'MarkerSize',8)
% 
%     for j = 1:min(size(offsetPos))
%         er = errorbar(offsetPos(i,j),mean(allData{i}(:,j)),std(allData{i}(:,j))/sqrt(length(allData{i}(:,j))));
%         er.Color = [0 0 0];
%         er.LineStyle = 'none';
%     end
% end
% 
% legend([{'Congruent'},{'Incongruent'}])
% ylabel('Speed (Pixels / Second)')
% ylim([0 100])
% 
% xticklabels({'PreCue', 'CueOn','SipIn'})
% 
% set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
% 
% saveas(gca,[figPath filesep 'speedtocorrectsipper_correctTrials'     '_100ms' ],'svg');     saveas(gca,[figPath filesep 'speedtocorrectsipper_correctTrials'    '_100ms' ],'png')

%% INCORRECT
%% Pull data, save trial structure if possible

WCon = {trlStruct(regWIdx).allInCorrectAppDistance};
WInc = {trlStruct(revWIdx).allInCorrectAppDistance};

PCon = {trlStruct(regPIdx).allInCorrectAppDistance};
PInc = {trlStruct(revPIdx).allInCorrectAppDistance};

% NanPad 
% Not all sessions have 15 trials. This is rough! Therefore, we are padding
% the trials with NaNs to circumvent the headache this causes.

endPad = 15;
for i = 1:length(WCon)
    sizeIdxWC = size(WCon{i});
    WCon{i}(sizeIdxWC(1):endPad,:) = NaN;

    sizeIdxWIC = size(WInc{i});
    WInc{i}(sizeIdxWIC(1):endPad,:) = NaN;
end
for i = 1:length(PCon)
    sizeIdxPC = size(PCon{i});
    PCon{i}(sizeIdxPC(1):endPad,:) = NaN;

    sizeIdxPIC = size(PInc{i});
    PInc{i}(sizeIdxPIC(1):endPad,:) = NaN;
end
WCon_All = []; WInc_All = []; PCon_All = []; PInc_All = [];
binVal = [1:5; 6:10; 11:15];
for i = 1:length(WCon)
    for j = 1:min(size(binVal))
        WCon_All{i}(:,:,j) = cat(1,WCon{i}(binVal(j,:),:));
        WInc_All{i}(:,:,j) = cat(1,WInc{i}(binVal(j,:),:));
    end
end

for i = 1:length(PCon)
    for j = 1:min(size(binVal))
        PCon_All{i}(:,:,j) = PCon{i}(binVal(j,:),:);
        PInc_All{i}(:,:,j) = PInc{i}(binVal(j,:),:);
    end
end

WCon_All = cat(1,WCon_All{:});
WInc_All = cat(1,WInc_All{:});

PCon_All = cat(1,PCon_All{:});
PInc_All = cat(1,PInc_All{:});

%% Plot traces of each category vis a vis subplot 
% First group data by subject
% Grouping Var for Wistars
groupingVar2 = [];
groupingVar = ones(1,min(size(WCon_All)));
idxVar = 1:5:65;
for i = 1:length(idxVar)
    if i < 13
        groupingVar2(idxVar(i):idxVar(i+1) - 1,:) = i;
    elseif i == 13
        groupingVar2(idxVar(i):idxVar(i) + 4,:) = i;
    end
end

for i = 1:3
    sqz_WCon = squeeze(WCon_All(:,:,i));
    WCon_All_GRP(:,:,i) = grpstats(sqz_WCon,groupingVar2);
    sqz_WInc = squeeze(WInc_All(:,:,i));
    WInc_All_GRP(:,:,i) = grpstats(sqz_WInc,groupingVar2);
end
% P rats
groupingVar = ones(1,min(size(PCon_All)));
idxVar = 1:5:70;
for i = 1:length(idxVar)
    if i < 14
        groupingVar2(idxVar(i):idxVar(i+1) - 1,:) = i;
    elseif i == 14
        groupingVar2(idxVar(i):idxVar(i) + 4,:) = i;
    end
end
%
for i = 1:3
    sqz_PCon = squeeze(PCon_All(:,:,i));
    PCon_All_GRP(:,:,i) = grpstats(sqz_PCon,groupingVar2);
    sqz_PInc = squeeze(PInc_All(:,:,i));
    PInc_All_GRP(:,:,i) = grpstats(sqz_PInc,groupingVar2);
end

%% Plot traces of each category vis a vis subplot 

for i = 1:3
    WCon_Trl(i,:) = nanmean(squeeze(WCon_All_GRP(:,:,i)));
    WInc_Trl(i,:) = nanmean(squeeze(WInc_All_GRP(:,:,i)));
    PCon_Trl(i,:) = nanmean(squeeze(PCon_All_GRP(:,:,i)));
    PInc_Trl(i,:) = nanmean(squeeze(PInc_All_GRP(:,:,i)));

    WCon_TrlSEM(i,:) = nanstd(squeeze(WCon_All_GRP(:,:,i))/sqrt(13));
    WInc_TrlSEM(i,:) = nanstd(squeeze(WInc_All_GRP(:,:,i))/sqrt(13));
    PCon_TrlSEM(i,:) = nanstd(squeeze(PCon_All_GRP(:,:,i))/sqrt(14));
    PInc_TrlSEM(i,:) = nanstd(squeeze(PInc_All_GRP(:,:,i))/sqrt(14));
end
LineStyles = {'b-','r-','g-'};
figure('Units','normalized','Position',[0 0 1 1])
subplot(2,2,1)
for i = 1:3
    shadedErrorBar(time,WCon_Trl(i,:),WCon_TrlSEM(i,:),'lineprops',{LineStyles{i},'LineWidth',3}); 
end
ylim([0 500])
xline(cueOn/30,'k--','Cue On','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(sipAscent/30,'k--','SipOut','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(sipDescent/30,'k--','SipIn','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
title('Wistar Congruent')
xlabel('Time (s)')
ylabel('Distance (pixels)')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',25,'FontWeight','bold','LineWidth',5)

subplot(2,2,2)
for i = 1:3
    shadedErrorBar(time,WInc_Trl(i,:),WInc_TrlSEM(i,:),'lineprops',{LineStyles{i},'LineWidth',3}); 
end
ylim([0 500])
xline(cueOn/30,'k--','Cue On','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(sipAscent/30,'k--','SipOut','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(sipDescent/30,'k--','SipIn','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
title('Wistar Incongruent')
xlabel('Time (s)')
ylabel('Distance (pixels)')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',25,'FontWeight','bold','LineWidth',5)

subplot(2,2,3)
for i = 1:3
    shadedErrorBar(time,PCon_Trl(i,:),PCon_TrlSEM(i,:),'lineprops',{LineStyles{i},'LineWidth',3}); 
end
ylim([0 500])
xline(cueOn/30,'k--','Cue On','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(sipAscent/30,'k--','SipOut','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(sipDescent/30,'k--','SipIn','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
title('P rat Congruent')
xlabel('Time (s)')
ylabel('Distance (pixels)')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',25,'FontWeight','bold','LineWidth',5)

subplot(2,2,4)
for i = 1:3
    shadedErrorBar(time,PInc_Trl(i,:),PInc_TrlSEM(i,:),'lineprops',{LineStyles{i},'LineWidth',3}); 
end
ylim([0 500])
xline(cueOn/30,'k--','Cue On','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(sipAscent/30,'k--','SipOut','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(sipDescent/30,'k--','SipIn','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
title('P rat Incongruent')
xlabel('Time (s)')
ylabel('Distance (pixels)')
legend({'Trials 1-5', 'Trials 6-10', 'Trials 11-15'})
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',25,'FontWeight','bold','LineWidth',5)

saveas(gca,[figPath filesep 'distancetoincorrectsipper_incorrectTrials_4x4_binned15'     '_100ms' ],'svg');     saveas(gca,[figPath filesep 'distancetoincorrectsipper_incorrectTrials_4x4_binned15'    '_100ms' ],'png')
%%
%% Measure of modulation of each epoch (correct)

% Find the difference between each time point, divide by the lenght of time
% between those timepoints to obtain pixels per 1/30s. Multipy by 30 to
% obtain velocity in pixel /s. 
for i = 1:3
    preCueWCon3(:,:,i) = abs((WCon_All_GRP(:,1,i) - WCon_All_GRP(:,cueOn,i)) / length(1:cueOn)) * 30;
    preCueWInc3(:,:,i) = abs((WInc_All_GRP(:,1,i) - WInc_All_GRP(:,cueOn,i)) / length(1:cueOn)) * 30;
    preCuePCon3(:,:,i) = abs((PCon_All_GRP(:,1,i) - PCon_All_GRP(:,cueOn,i)) / length(1:cueOn)) * 30;
    preCuePInc3(:,:,i) = abs((PInc_All_GRP(:,1,i) - PInc_All_GRP(:,cueOn,i)) / length(1:cueOn)) * 30;
    
    cueOnWCon3(:,:,i) = abs((WCon_All_GRP(:,cueOn,i) - WCon_All_GRP(:,sipDescent,i)) / length(cueOn:sipDescent)) * 30;
    cueOnWInc3(:,:,i) = abs((WInc_All_GRP(:,cueOn,i) - WInc_All_GRP(:,sipDescent,i)) / length(cueOn:sipDescent)) * 30;
    cueOnPCon3(:,:,i) = abs((PCon_All_GRP(:,cueOn,i) - PCon_All_GRP(:,sipDescent,i)) / length(cueOn:sipDescent)) * 30;
    cueOnPInc3(:,:,i) = abs((PInc_All_GRP(:,cueOn,i) - PInc_All_GRP(:,sipDescent,i)) / length(cueOn:sipDescent)) * 30;
    
    sipInWCon3(:,:,i) = abs((WCon_All_GRP(:,sipDescent,i) - WCon_All_GRP(:,(sipAscent),i)) / length(sipDescent:sipAscent)) * 30;
    sipInWInc3(:,:,i) = abs((WInc_All_GRP(:,sipDescent,i) - WInc_All_GRP(:,(sipAscent),i)) / length(sipDescent:sipAscent)) * 30;
    sipInPCon3(:,:,i) = abs((PCon_All_GRP(:,sipDescent,i) - PCon_All_GRP(:,(sipAscent),i)) / length(sipDescent:sipAscent)) * 30;
    sipInPInc3(:,:,i) = abs((PInc_All_GRP(:,sipDescent,i) - PInc_All_GRP(:,(sipAscent),i)) / length(sipDescent:sipAscent)) * 30;
end
%% Grab data and pull into single variable for group stats 
allCueData = {squeeze(cueOnWCon3),squeeze(cueOnWInc3),squeeze(cueOnPCon3),squeeze(cueOnPInc3)};

figure('Units','normalized','Position',[0 0 1 1])
hb = bar([nanmean(allCueData{1}(:,1)) nanmean(allCueData{2}(:,1)) nanmean(allCueData{3}(:,1)) nanmean(allCueData{4}(:,1)); nanmean(allCueData{1}(:,2)) nanmean(allCueData{2}(:,2)) nanmean(allCueData{3}(:,2)) nanmean(allCueData{4}(:,2)); nanmean(allCueData{1}(:,3)) nanmean(allCueData{2}(:,3)) nanmean(allCueData{3}(:,3)) nanmean(allCueData{4}(:,3))], 'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);

offsetPos = [];
scatterXData = [];


for i = 1:length(hb)
    offsetPos(i,:) = hb(i).XData + hb(i).XOffset;
end
for i = 1:length(offsetPos)
    if i == 1
        dataLength = size(WCon,1);
    elseif i == 2
        dataLength = size(WInc,1);
    elseif i == 3
        dataLength = size(PCon,1);
    elseif i == 4
        dataLength = size(PInc,1);
    end
    for j = 1:min(size(offsetPos))
        scatterXData{i,j} = linspace(offsetPos(i,j)-0.01,offsetPos(i,j)+0.01,dataLength);
    end
end
hold on

for i = 1:length(offsetPos)
    plot([scatterXData{i,:}],allCueData{i},'ko','LineWidth',2,'MarkerSize',8)

    for j = 1:min(size(offsetPos))
        er = errorbar(offsetPos(i,j),nanmean(allCueData{i}(:,j)),nanstd(allCueData{i}(:,j))/sqrt(length(allCueData{i}(:,j))));
        er.Color = [0 0 0];
        er.LineStyle = 'none';
        er.LineWidth = 3;
    end
end

ylabel('Speed (Pixels / Second)')
xticklabels({'Trials 1-5', 'Trials 6-10','Trials 11-15'})
ylim([0 100])

set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)

saveas(gca,[figPath filesep 'speedtoincorrectsipper_incorrectTrials_15Bin'     '_100ms' ],'svg');     saveas(gca,[figPath filesep 'speedtoincorrectsipper_incorrectTrials_15Bin'    '_100ms' ],'png')


%%

% formatData = cell2mat(allCueData);
% formatData = formatData';
% datatable = array2table(formatData,'VariableNames',{'Block1','Block2','Block3'});
% datatable.Strain = [repmat({'Wistar'},30,1); repmat({'P rat'},30,1)];
% condition = [repmat({'Congruent'},15,2); repmat({'Incongruent'},15,2)];
% condition = {(condition{:})};
% datatable.Condition = condition';
% 
% % Create rm model
% Block = table([1 2 3]','VariableNames',{'Blocks'});
% rm = fitrm(datatable,'Block1-Block3 ~ Strain * Condition','WithinDesign',Block);
% anovaTbl = ranova(rm,'WithinModel','Blocks')


%% Pull data based on categories (correct)

WCon = cellfun(@mean,{trlStruct(regWIdx).allCorrectAppDistance},'UniformOutput',false);
WCon = cat(1,WCon{:});
WInc = cellfun(@mean,{trlStruct(revWIdx).allCorrectAppDistance},'UniformOutput',false);
WInc = cat(1,WInc{:});

PCon = cellfun(@mean,{trlStruct(regPIdx).allCorrectAppDistance},'UniformOutput',false);
PCon = cat(1,PCon{:});
PInc = cellfun(@mean,{trlStruct(revPIdx).allCorrectAppDistance},'UniformOutput',false);
PInc = cat(1,PInc{:});

%% Make table to peform RANOVA on distance data. (Correct)

tableStarter = [WCon; WInc; PCon; PInc];
strainVar = [repmat({'Wistar'},1,min(size(WCon))+min(size(WInc))) repmat({'P rat'},1,min(size(PCon))+min(size(PInc)))];
strainVar = strainVar';
sessionVar = [repmat({'Congruent'},1,min(size(WCon))) repmat({'Incongruent'},1,min(size(WInc))) repmat({'Congruent'},1,min(size(PCon))) repmat({'Incongruent'},1,min(size(PInc)))];
sessionVar = sessionVar';
TimeVar = table([1:751]','VariableNames',{'Time'});
% TimeVarTable = mat2cell([1:751],1,751);
dataTable = array2table(tableStarter);
dataTable.Strain = strainVar;
dataTable.Session = sesssionVar;


% Create rm model
rm = fitrm(dataTable,'tableStarter1-tableStarter751 ~ Strain * Session','WithinDesign',TimeVar);
anovaTbl = ranova(rm,'WithinModel','Time')
multTbl1_Incor = multcompare(rm,'Strain','By','Session')
multTbl2_Incor = multcompare(rm,'Strain')

%% First pass plot of distance over time
sipDescent = 10*30;
sipAscent = 18*30;
cueOn = 5*30;
% SEM for ShadedErrorBar
wConError = std(WCon)/sqrt(size(WCon,1));
wIncError = std(WInc)/sqrt(size(WInc,1));

PConError = std(PCon)/sqrt(size(PCon,1));
PIncError = std(PInc)/sqrt(size(PInc,1));

figure('Units','normalized','Position',[0 0 1 1])

shadedErrorBar(time,mean(WCon),wConError,'lineprops',{'b','LineWidth',3});
hold on
shadedErrorBar(time,mean(WInc),wIncError,'lineprops',{'b--','LineWidth',3});
shadedErrorBar(time,mean(PCon),PConError,'lineprops',{'r','LineWidth',3});
shadedErrorBar(time,mean(PInc),PIncError,'lineprops',{'r--','LineWidth',3});

xline(cueOn/30,'k--','Cue On','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(sipAscent/30,'k--','SipOut','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(sipDescent/30,'k--','SipIn','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
ylim([0 450])

xlabel('Time (s)')
ylabel('Distance from Correct Sipper (pixels)')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',25,'FontWeight','bold','LineWidth',5)

saveas(gca,[figPath filesep 'distancetocorrectsipper_correctTrials'     '_100ms' ],'svg');     saveas(gca,[figPath filesep 'distancetocorrectsipper_correctTrials'    '_100ms' ],'png')

%% Measure of modulation of each epoch
% Find the difference between each time point, divide by the lenght of time
% between those timepoints to obtain pixels per 1/30s. Multipy by 30 to
% obtain velocity in pixel /s. 

preCueWCon = abs((WCon(:,1) - WCon(:,cueOn)) / length(1:cueOn)) * 30;
preCueWInc = abs((WInc(:,1) - WInc(:,cueOn)) / length(1:cueOn)) * 30;
preCuePCon = abs((PCon(:,1) - PCon(:,cueOn)) / length(1:cueOn)) * 30;
preCuePInc = abs((PInc(:,1) - PInc(:,cueOn)) / length(1:cueOn)) * 30;

cueOnWCon = abs((WCon(:,cueOn) - WCon(:,sipDescent)) / length(cueOn:sipDescent)) * 30;
cueOnWInc = abs((WInc(:,cueOn) - WInc(:,sipDescent)) / length(cueOn:sipDescent)) * 30;
cueOnPCon = abs((PCon(:,cueOn) - PCon(:,sipDescent)) / length(cueOn:sipDescent)) * 30;
cueOnPInc = abs((PInc(:,cueOn) - PInc(:,sipDescent)) / length(cueOn:sipDescent)) * 30;

sipInWCon = abs((WCon(:,sipDescent) - WCon(:,(sipAscent))) / length(sipDescent:sipAscent)) * 30;
sipInWInc = abs((WInc(:,sipDescent) - WInc(:,(sipAscent))) / length(sipDescent:sipAscent)) * 30;
sipInPCon = abs((PCon(:,sipDescent) - PCon(:,(sipAscent))) / length(sipDescent:sipAscent)) * 30;
sipInPInc = abs((PInc(:,sipDescent) - PInc(:,(sipAscent))) / length(sipDescent:sipAscent)) * 30;

%% Plot velocity
figure('Units','normalized','Position',[0 0 1 1])
subplot(1,3,1)
bar([mean(preCueWCon) mean(preCueWInc); mean(preCuePCon) mean(preCuePInc)])
title('PreCue')
xticklabels({'Wistars', 'P rats'})
ylabel('Velocity (pixels / s)')
ylim([0 30])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',25,'FontWeight','bold','LineWidth',5)


subplot(1,3,2)

bar([mean(cueOnWCon) mean(cueOnWInc); mean(cueOnPCon) mean(cueOnPInc)])
title('CueOn')
xticklabels({'Wistars', 'P rats'})
ylabel('Velocity (pixels / s)')
ylim([0 30])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',25,'FontWeight','bold','LineWidth',5)

subplot(1,3,3)
bar([mean(sipInWCon) mean(sipInWInc); mean(sipInPCon) mean(sipInPInc)])
title('SipperIn')
xticklabels({'Wistars', 'P rats'})
ylabel('Velocity (pixels / s)')
ylim([0 30])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',25,'FontWeight','bold','LineWidth',5)

figure('Units','normalized','Position',[0 0 1 1])
bar([mean(preCueWCon) mean(preCueWInc) mean(preCuePCon) mean(preCuePInc); mean(cueOnWCon) mean(cueOnWInc) mean(cueOnPCon) mean(cueOnPInc); mean(sipInWCon) mean(sipInWInc) mean(sipInPCon) mean(sipInPInc)])
xticklabels({'PreCue', 'CueOn','SipIn'})
ylabel('Velocity (pixels / s)')
ylim([0 30])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',25,'FontWeight','bold','LineWidth',5)

%% Fancy bar graph
figure('Units','normalized','Position',[0 0 1 1])
hb = bar([mean(preCueWCon) mean(preCueWInc) mean(preCuePCon) mean(preCuePInc); mean(cueOnWCon) mean(cueOnWInc) mean(cueOnPCon) mean(cueOnPInc); mean(sipInWCon) mean(sipInWInc) mean(sipInPCon) mean(sipInPInc)], 'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);

offsetPos = [];
scatterXData = [];

for i = 1:length(hb)
    offsetPos(i,:) = hb(i).XData + hb(i).XOffset;
end
for i = 1:length(offsetPos)
    if i == 1
        dataLength = size(WCon,1);
    elseif i == 2
        dataLength = size(WInc,1);
    elseif i == 3
        dataLength = size(PCon,1);
    elseif i == 4
        dataLength = size(PInc,1);
    end
    for j = 1:min(size(offsetPos))
        scatterXData{i,j} = linspace(offsetPos(i,j)-0.01,offsetPos(i,j)+0.01,dataLength);
    end
end
WConData = [preCueWCon cueOnWCon sipInWCon];
WIncData = [preCueWInc cueOnWInc sipInWInc];

PConData = [preCuePCon cueOnPCon sipInPCon];
PIncData = [preCuePInc cueOnPInc sipInPInc];
hold on
allData = {WConData, WIncData, PConData, PIncData};

cueOnData = {cueOnWCon cueOnWInc cueOnPCon cueOnPInc};
for i = 1:length(offsetPos)
    plot([scatterXData{i,:}],allData{i}(:),'ko','LineWidth',2,'MarkerSize',8)

    for j = 1:min(size(offsetPos))
        er = errorbar(offsetPos(i,j),mean(allData{i}(:,j)),std(allData{i}(:,j))/sqrt(length(allData{i}(:,j))));
        er.Color = [0 0 0];
        er.LineStyle = 'none';
    end
end

legend([{'Congruent'},{'Incongruent'}])
ylabel('Speed (Pixels / Second)')
xticklabels({'PreCue', 'CueOn','SipIn'})
ylim([0 50])

set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)

saveas(gca,[figPath filesep 'speedtocorrectsipper_correctTrials'     '_100ms' ],'svg');     saveas(gca,[figPath filesep 'speedtocorrectsipper_correctTrials'    '_100ms' ],'png')


%%
%% Pull data based on categories (incorrect)

WCon = cellfun(@mean,{trlStruct(regWIdx).allInCorrectAppDistance},'UniformOutput',false);
WCon = cat(1,WCon{:});
WInc = cellfun(@mean,{trlStruct(revWIdx).allInCorrectAppDistance},'UniformOutput',false);
WInc = cat(1,WInc{:});

PCon = cellfun(@mean,{trlStruct(regPIdx).allInCorrectAppDistance},'UniformOutput',false);
PCon = cat(1,PCon{:});
PInc = cellfun(@mean,{trlStruct(revPIdx).allInCorrectAppDistance},'UniformOutput',false);
PInc = cat(1,PInc{:});


%% Make table to peform RANOVA on distance data.

tableStarter = [WCon; WInc; PCon; PInc];
strainVar = [repmat({'Wistar'},1,min(size(WCon))+min(size(WInc))) repmat({'P rat'},1,min(size(PCon))+min(size(PInc)))];
strainVar = strainVar';
sessionVar = [repmat({'Congruent'},1,min(size(WCon))) repmat({'Incongruent'},1,min(size(WInc))) repmat({'Congruent'},1,min(size(PCon))) repmat({'Incongruent'},1,min(size(PInc)))];
sessionVar = sessionVar';
TimeVar = table([1:751]','VariableNames',{'Time'});
% TimeVarTable = mat2cell([1:751],1,751);
dataTable = array2table(tableStarter);
dataTable.Strain = strainVar;
dataTable.Session = sesssionVar;


% Create rm model
rm = fitrm(dataTable,'tableStarter1-tableStarter751 ~ Strain * Session')%,'WithinDesign',TimeVar);
anovaTbl = ranova(rm,'WithinModel','Time')
multTbl1_Incor = multcompare(rm,'Strain','By','Session')
multTbl2_Incor = multcompare(rm,'Strain')
multTbl2_Incor = multcompare(rm,'Session')


%% First pass plot of distance over time
sipDescent = 10*30;
sipAscent = 18*30;
cueOn = 5*30;
% SEM for ShadedErrorBar
wConError = std(WCon)/sqrt(size(WCon,1));
wIncError = std(WInc)/sqrt(size(WInc,1));

PConError = std(PCon)/sqrt(size(PCon,1));
PIncError = std(PInc)/sqrt(size(PInc,1));

figure('Units','normalized','Position',[0 0 1 1])

shadedErrorBar(time,mean(WCon),wConError,'lineprops',{'b','LineWidth',3});
hold on
shadedErrorBar(time,mean(WInc),wIncError,'lineprops',{'b--','LineWidth',3});
shadedErrorBar(time,mean(PCon),PConError,'lineprops',{'r','LineWidth',3});
shadedErrorBar(time,mean(PInc),PIncError,'lineprops',{'r--','LineWidth',3});

xline(cueOn/30,'k--','Cue On','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(sipAscent/30,'k--','SipOut','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(sipDescent/30,'k--','SipIn','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');

xlabel('Time (s)')
ylabel('Distance from Correct Sipper (pixels)')
ylim([0 450])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',25,'FontWeight','bold','LineWidth',5)

saveas(gca,[figPath filesep 'distancetocorrectsipper_incorrectTrials'     '_100ms' ],'svg');     saveas(gca,[figPath filesep 'distancetocorrectsipper_incorrectTrials'    '_100ms' ],'png')

%% Measure of modulation of each epoch
% Find the difference between each time point, divide by the lenght of time
% between those timepoints to obtain pixels per 1/30s. Multipy by 30 to
% obtain velocity in pixel /s. 

preCueWCon = abs((WCon(:,1) - WCon(:,cueOn)) / length(1:cueOn)) * 30;
preCueWInc = abs((WInc(:,1) - WInc(:,cueOn)) / length(1:cueOn)) * 30;
preCuePCon = abs((PCon(:,1) - PCon(:,cueOn)) / length(1:cueOn)) * 30;
preCuePInc = abs((PInc(:,1) - PInc(:,cueOn)) / length(1:cueOn)) * 30;

cueOnWCon = abs((WCon(:,cueOn) - WCon(:,sipDescent)) / length(cueOn:sipDescent)) * 30;
cueOnWInc = abs((WInc(:,cueOn) - WInc(:,sipDescent)) / length(cueOn:sipDescent)) * 30;
cueOnPCon = abs((PCon(:,cueOn) - PCon(:,sipDescent)) / length(cueOn:sipDescent)) * 30;
cueOnPInc = abs((PInc(:,cueOn) - PInc(:,sipDescent)) / length(cueOn:sipDescent)) * 30;

sipInWCon = abs((WCon(:,sipDescent) - WCon(:,(sipAscent))) / length(sipDescent:sipAscent)) * 30;
sipInWInc = abs((WInc(:,sipDescent) - WInc(:,(sipAscent))) / length(sipDescent:sipAscent)) * 30;
sipInPCon = abs((PCon(:,sipDescent) - PCon(:,(sipAscent))) / length(sipDescent:sipAscent)) * 30;
sipInPInc = abs((PInc(:,sipDescent) - PInc(:,(sipAscent))) / length(sipDescent:sipAscent)) * 30;

%% Plot velocity


%% Fancy bar graph
figure('Units','normalized','Position',[0 0 1 1])
hb = bar([mean(preCueWCon) mean(preCueWInc) mean(preCuePCon) mean(preCuePInc); mean(cueOnWCon) mean(cueOnWInc) mean(cueOnPCon) mean(cueOnPInc); mean(sipInWCon) mean(sipInWInc) mean(sipInPCon) mean(sipInPInc)], 'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);

offsetPos = [];
scatterXData = [];

for i = 1:length(hb)
    offsetPos(i,:) = hb(i).XData + hb(i).XOffset;
end
for i = 1:length(offsetPos)
    if i == 1
        dataLength = size(WCon,1);
    elseif i == 2
        dataLength = size(WInc,1);
    elseif i == 3
        dataLength = size(PCon,1);
    elseif i == 4
        dataLength = size(PInc,1);
    end
    for j = 1:min(size(offsetPos))
        scatterXData{i,j} = linspace(offsetPos(i,j)-0.01,offsetPos(i,j)+0.01,dataLength);
    end
end
WConData = [preCueWCon cueOnWCon sipInWCon];
WIncData = [preCueWInc cueOnWInc sipInWInc];

PConData = [preCuePCon cueOnPCon sipInPCon];
PIncData = [preCuePInc cueOnPInc sipInPInc];
hold on
allData = {WConData, WIncData, PConData, PIncData};

cueOnData = {cueOnWCon cueOnWInc cueOnPCon cueOnPInc};
for i = 1:length(offsetPos)
    plot([scatterXData{i,:}],allData{i}(:),'ko','LineWidth',2,'MarkerSize',8)

    for j = 1:min(size(offsetPos))
        er = errorbar(offsetPos(i,j),mean(allData{i}(:,j)),std(allData{i}(:,j))/sqrt(length(allData{i}(:,j))));
        er.Color = [0 0 0];
        er.LineStyle = 'none';
    end
end

legend([{'Congruent'},{'Incongruent'}])
ylabel('Speed (Pixels / Second)')
xticklabels({'PreCue', 'CueOn','SipIn'})
ylim([0 50])

set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)

saveas(gca,[figPath filesep 'speedtocorrectsipper_incorrectTrials'     '_100ms' ],'svg');     saveas(gca,[figPath filesep 'speedtocorrectsipper_incorrectTrials'    '_100ms' ],'png')
