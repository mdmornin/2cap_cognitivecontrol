%% Combine Nick's Code with my own to create PSTHs of Distances
% Variable information

%% approach : contains nTrials x approachType x approachLength
% % nTrials should equal 96 and is organized by type from variable
% trialTypes

% % approachTypes are either correct approach during access, incorrect
% approach during access, correct approach during cue, incorrect approach
% during cue

% % 3rd dimension either has binary {1,0} signifying approach types or a
% value in seconds specifying amount of time at cue or descended sipper

%% boxCoords : contains coordinates {X,Y} and Likelihood for each box pos
% boxCoords = [sipLDLC;sipRDLC;ULC;LLC;URC;LRC]

%% bodyCoords : contains coordinates {X,Y} and likelhiood for each body part
% bodyCoords = {frameT,SnoutTip,HeadCap,MidBack,BackTail};

%% trialTypes & trialTimes
% trialTypes : 1 = CS Plus Left, 2 = CS Plus Right, 3 = CS Negative
% trialTimes : time of each trial organized by type and not continuous

%% Goal -- Create PSTHs of distance from sipper then organize by approach variable
% Calculate distance from left and right sipper at every time point and
% then create cue-centered PSTH for both the left and right sipper at that
% time point. 
%% Loads & Paths
addpath(genpath('F:\dissDat\analysis-tools-master'))
path = 'F:\dissDat';
load(['F:\dissDat\restoredScripts\masterTable.mat'])
addpath(genpath('F:\dissDat\rawData'))
fpath = 'F:/dissDat'
NewDataDir = [fpath,filesep,'res'];
DLCDir = [fpath,filesep,'dlcRes'];
DLCsuf = '_videoDLC_resnet50_2CAPFullRun1Mar28shuffle1_1030000.csv';
dirList = dir([fpath,filesep,'rawData']);
folderNames = {dirList([dirList.isdir]).name};
folderNames = folderNames(~ismember(folderNames ,{'.','..'}));
dataSetNames = folderNames;
dataDirs = dataSetNames;
for i = 1:length(dataSetNames)
    splitName = strsplit(dataSetNames{i},'-');
    DLCNames{i} = [splitName{1} '-' splitName{2} '-' splitName{3} splitName{4}];
end
dataSetIDs = DLCNames;
stMtx = []; % For now

for i = 1:length(dataSetIDs)
    load([path,filesep,'res',filesep,dataSetIDs{i} '_AllData.mat'])
    dlcTimestamps = bodyCoords{1};
    nosePos = bodyCoords{2};
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
        end
    end
    
    % Update structure

    trlStruct(i).Session = dataSetIDs{i};
    trlStruct(i).trlLSipDist = trlLSipDist;
    trlStruct(i).trlRSipDist = trlRSipDist;
    trlStruct(i).trlnosePos = trlnosePos;

    % Add Session Type Conditional

    if startsWith(masterTbl.SessionType{i},'Regular')
        LcorrectIdx = (approach(LeftTrials,1,1));
        LincorrectIdx = (approach(LeftTrials,2,1));
        RcorrectIdx = (approach(RightTrials,1,1));
        RincorrectIdx = (approach(RightTrials,2,1));
        LcorrectApproachTrls = trlLSipDist(LeftTrials,:);
        RcorrectApproachTrls = trlRSipDist(RightTrials,:);
    elseif startsWith(masterTbl.SessionType{i},'Reversal')
        LcorrectIdx = (approach(LeftTrials,2,1));
        LincorrectIdx = (approach(LeftTrials,1,1));

        RcorrectIdx = (approach(RightTrials,2,1));
        RincorrectIdx = (approach(RightTrials,1,1));

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

    masterTbl.NeuronYield(i) = min(size(stMtx));

end

%% Find some summary information on Reg and Rev Trials
close all
sessionNum = num2str([]);
strain = 'P';

regLDist = []; regRDist = []; revLDist = []; revRDist = [];


for i = 1:length(dataSetIDs)
    if strcmp(masterTbl.SessionType{i},['Regular' sessionNum]) && strcmp(masterTbl.Strain{i},strain)
        regLDist(:,i) = nanmedian(trlStruct(i).LDistCorrectAppr,1);
        regRDist(:,i) = nanmedian(trlStruct(i).RDistCorrectAppr,1);
    elseif strcmp(masterTbl.SessionType{i},['Reversal' sessionNum]) && strcmp(masterTbl.Strain{i},strain)
        revLDist(:,i) = nanmedian(trlStruct(i).LDistCorrectAppr,1);
        revRDist(:,i) = nanmedian(trlStruct(i).RDistCorrectAppr,1);
    end

end

% Remove 0's 
regLDist=regLDist(:,any(regLDist)); revLDist=revLDist(:,any(revLDist));
regRDist=regRDist(:,any(regRDist)); revRDist=revRDist(:,any(revRDist));

% Test plots
timeAdj = 1:length(regLDist);
timeAdj = timeAdj./30;
timeAdj = timeAdj - 5;
figure('Units','normalized','Position',[0 0 1 1])
plot(timeAdj,nanmedian([regLDist regRDist],2),'r-','LineWidth',3);
hold on
plot(timeAdj,nanmedian([revLDist revRDist],2),'b-','LineWidth',3);
xlabel('Time from Cue Onset (s)')
ylabel('Distance from Correct Sipper')
xline(0,'k--','Cue On','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(5,'k--','Sipper In','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(13,'k--','Sipper Out','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
title('Correct Approach Trials')
legend([{'Congruent'},{'Incongruent'}])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,['correctApproachDist' sessionNum strain],'png')

% Determine Correct Approaches per Session Type
regCorrect = []; revCorrect = [];

for i = 1:length(dataSetIDs)
    if strcmp(masterTbl.SessionType{i},['Regular' sessionNum]) && strcmp(masterTbl.Strain{i},strain)
        regCorrect = [regCorrect size(trlStruct(i).LDistCorrectAppr,1) + size(trlStruct(i).RDistCorrectAppr,1)];
    elseif strcmp(masterTbl.SessionType{i},['Reversal' sessionNum]) && strcmp(masterTbl.Strain{i},strain)
        revCorrect = [revCorrect size(trlStruct(i).LDistCorrectAppr,1) + size(trlStruct(i).RDistCorrectAppr,1)];
    end

end
semregCR = std(regCorrect)/sqrt(length(regCorrect));
semrevCR = std(revCorrect)/sqrt(length(revCorrect));

figure('Units','normalized','Position',[0 0 1 1])
bar(categorical({'Congruent Correct Approach','Incongruent Correct Approach'}),[mean(regCorrect) mean(revCorrect)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3)
hold on
plot([linspace(0.95,1.05,length(regCorrect));linspace(1.95,2.05,length(revCorrect))],[regCorrect; revCorrect],'k--o','LineWidth',3,'MarkerSize',10)
er = errorbar(categorical({'Congruent Correct Approach','Incongruent Correct Approach'}),[mean(regCorrect) mean(revCorrect)],[semregCR semrevCR],'LineWidth',3);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
ylabel('Trials')
ylim([0 48])

set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,['correctApproachNum' sessionNum strain],'png')



% Line incorrect
regLDist = []; regRDist = []; revLDist = []; revRDist = [];


for i = 1:length(dataSetIDs)
    if strcmp(masterTbl.SessionType{i},['Regular' sessionNum]) && strcmp(masterTbl.Strain{i},strain)
        regLDist(:,i) = nanmedian(trlStruct(i).LDistInCorrectAppr,1);
        regRDist(:,i) = nanmedian(trlStruct(i).RDistInCorrectAppr,1);
    elseif strcmp(masterTbl.SessionType{i},['Reversal' sessionNum]) && strcmp(masterTbl.Strain{i},strain)
        revLDist(:,i) = nanmedian(trlStruct(i).LDistInCorrectAppr,1);
        revRDist(:,i) = nanmedian(trlStruct(i).RDistInCorrectAppr,1);
    end

end

% Remove 0's 
regLDist=regLDist(:,any(regLDist)); revLDist=revLDist(:,any(revLDist));
regRDist=regRDist(:,any(regRDist)); revRDist=revRDist(:,any(revRDist));

% Test plots
timeAdj = 1:length(regLDist);
timeAdj = timeAdj./30;
timeAdj = timeAdj - 5;
figure('Units','normalized','Position',[0 0 1 1])
plot(timeAdj,nanmedian([regLDist regRDist],2),'r-','LineWidth',3);
hold on
plot(timeAdj,nanmedian([revLDist revRDist],2),'b--','LineWidth',3);
xlabel('Time from Cue Onset (s)')
ylabel('Distance from Correct Sipper')
xline(0,'k--','Cue On','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(5,'k--','Sipper In','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(13,'k--','Sipper Out','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
title('Incorrect Approach Trials')
legend([{'Congruent'},{'Incongruent'}])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,['incorrectApproachDist' sessionNum strain],'png')

% Determine incorrect Approaches per Session Type
regInCorrect = []; revInCorrect = [];

for i = 1:length(dataSetIDs)
    if strcmp(masterTbl.SessionType{i},['Regular' sessionNum]) && strcmp(masterTbl.Strain{i},strain)
        regInCorrect = [regInCorrect size(trlStruct(i).LDistInCorrectAppr,1) + size(trlStruct(i).RDistInCorrectAppr,1)];
    elseif strcmp(masterTbl.SessionType{i},['Reversal' sessionNum]) && strcmp(masterTbl.Strain{i},strain)
        revInCorrect = [revInCorrect size(trlStruct(i).LDistInCorrectAppr,1) + size(trlStruct(i).RDistInCorrectAppr,1)];
    end

end

semregIN = std(regInCorrect)/sqrt(length(regInCorrect));
semrevIN = std(revInCorrect)/sqrt(length(revInCorrect));

figure('Units','normalized','Position',[0 0 1 1])
bar(categorical({'Congruent Incorrect Approach','Incongruent Incorrect Approach'}),[mean(regInCorrect) mean(revInCorrect)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3)
hold on
plot([linspace(0.95,1.05,length(regInCorrect));linspace(1.95,2.05,length(revInCorrect))],[regInCorrect; revInCorrect],'k--o','LineWidth',3,'MarkerSize',10)
er = errorbar(categorical({'Congruent Incorrect Approach','Incongruent Incorrect Approach'}),[mean(regInCorrect) mean(revInCorrect)],[semregIN semrevIN],'LineWidth',3);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
ylabel('Trials')
ylim([0 48])

set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,['incorrectApproachNum' sessionNum strain],'png')

% Find omissions
regOmit = []; revOmit = [];

for i = 1:length(dataSetIDs)
    if strcmp(masterTbl.SessionType{i},['Regular' sessionNum]) && strcmp(masterTbl.Strain{i},strain)
        regOmit = [regOmit sum(sum([trlStruct(i).approach(1:48,1:2,1) == 1],2)==0)];
    elseif strcmp(masterTbl.SessionType{i},['Reversal' sessionNum]) && strcmp(masterTbl.Strain{i},strain)
        revOmit = [revOmit sum(sum([trlStruct(i).approach(1:48,1:2,1) == 1],2)==0)];
    end

end
semregOmit = std(regOmit)/sqrt(length(regOmit));
semrevcorLat = std(revOmit)/sqrt(length(revOmit));

figure('Units','normalized','Position',[0 0 1 1])
bar(categorical({'Congruent Omissions','Incongruent Omissions'}),[mean(regOmit) mean(revOmit)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3)
hold on
plot([linspace(0.95,1.05,length(regOmit));linspace(1.95,2.05,length(revOmit))],[regOmit; revOmit],'k--o','LineWidth',3,'MarkerSize',10)
er = errorbar(categorical({'Congruent Omissions','Incongruent Omissions'}),[mean(regOmit) mean(revOmit)],[semregOmit semrevcorLat],'LineWidth',3);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
ylabel('Trials')
ylim([0 48])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,['omitNum' sessionNum strain],'png')

[h,p,stat] = ttest(regOmit,revOmit);

% Determine secondary statistics such as latency
% Find first moment following sipper descent where animal makes contact
% with sipper pixels. Use threshold set by Nick's code (< 9 pixels) and
% threshold for timebins (> 3 timebins). First start by indexing all
% moments where distance to correct sipper is less than 9 pixels, then
% ensure that the # of pixels is greater than 3. The first index where this
% occurs will be the latency to approach the sipper. 
disThresh = 9;
disRun = 3;
sipDescent = 10*30;
sipAscent = 18*30;
cueOn = 5*30;
ReglatCor = []; ReglatinCor = []; RevlatCor = []; RevlatinCor = []; % Correct and Incorrect latencies

for i = 1:length(dataSetIDs)
    if strcmp(masterTbl.SessionType{i},['Regular' sessionNum]) && strcmp(masterTbl.Strain{i},strain)
        latCorIdx = [trlStruct(i).RDistCorrectAppr <= disThresh; trlStruct(i).LDistCorrectAppr <= disThresh];
        latInCorIdx = [trlStruct(i).RDistInCorrectAppr <= disThresh; trlStruct(i).LDistInCorrectAppr <= disThresh];
        [~,sessCorLats] = max(latCorIdx(:,cueOn:sipAscent),[],2);
        [~,sessInCorLats] = max(latInCorIdx(:,cueOn:sipAscent),[],2);
        ReglatCor = [ReglatCor; sessCorLats/30;];
        ReglatinCor = [ReglatinCor; sessInCorLats/30;];
    elseif strcmp(masterTbl.SessionType{i},['Reversal' sessionNum]) && strcmp(masterTbl.Strain{i},strain)
        latCorIdx = [trlStruct(i).RDistCorrectAppr <= disThresh; trlStruct(i).LDistCorrectAppr <= disThresh];
        latInCorIdx = [trlStruct(i).RDistInCorrectAppr <= disThresh; trlStruct(i).LDistInCorrectAppr <= disThresh];
        [~,sessCorLats] = max(latCorIdx(:,cueOn:sipAscent),[],2);
        [~,sessInCorLats] = max(latInCorIdx(:,cueOn:sipAscent),[],2);
        RevlatCor = [RevlatCor; sessCorLats/30;];
        RevlatinCor = [RevlatinCor; sessInCorLats/30;];
    end
end

% Create latency figure
semregcorLat = std(ReglatCor)/sqrt(length(ReglatCor));
semrevcorLat = std(RevlatCor)/sqrt(length(RevlatCor));
semregIncorLat = std(ReglatinCor)/sqrt(length(ReglatinCor));
semrevIncorLat = std(RevlatinCor)/sqrt(length(RevlatinCor));

figure('Units','normalized','Position',[0 0 1 1])
hb = bar(categorical({'Correct Latencies','Incorrect Latencies'}),[mean(ReglatCor) mean(RevlatCor); mean(ReglatinCor) mean(RevlatinCor)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);
hb(1).FaceColor = [0.5 0.5 0.5]; hb(2).FaceColor = [0.9 0.9 0.9];
offsetPos = [1+hb(1).XOffset 1+hb(2).XOffset; 2+hb(1).XOffset 2+hb(2).XOffset];
scatterXData = [linspace(offsetPos(1)-0.1,offsetPos(1)+0.1,length(ReglatCor)),linspace(offsetPos(2)-0.1,offsetPos(2)+0.1,length(ReglatinCor)),linspace(offsetPos(3)-0.1,offsetPos(3)+0.1,length(RevlatCor)),linspace(offsetPos(4)-0.1,offsetPos(4)+0.1,length(RevlatinCor))];
hold on
scatter(scatterXData',[ReglatCor; RevlatCor; ReglatinCor; RevlatinCor],50,'ko','LineWidth',2)
er = errorbar(offsetPos,[mean(ReglatCor) mean(RevlatCor); mean(ReglatinCor) mean(RevlatinCor)],[semregcorLat semrevcorLat; semregIncorLat semrevIncorLat],'LineWidth',3);    
er(1).Color = [0 0 0]; er(2).Color = [0 0 0];                          
er(1).LineStyle = 'none'; er(2).LineStyle = 'none'; 
ylabel('Time (s)')
legend([{'Congruent'},{'Incongruent'}])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
saveas(gca,['Latencies' sessionNum strain],'png')
% xData = hb(1).XData+hb(1).XOffset;

%% Approach Likelihood, Metrics over Trials, Correct or Incorrect
%sessionNum = []; 
close all

sessionNum = num2str([]);

strain = 'P';

RegcorApproachVector = []; RevcorApproachVector = [];
RegincApproachVector = []; RevincApproachVector = [];
for i = 1:length(dataSetIDs)
    if strcmp(masterTbl.SessionType{i},['Regular' sessionNum]) && strcmp(masterTbl.Strain{i},strain)
        currentTrialTimes = trlStruct(i).trialTimes;
        [sortedTimes,sortedIndex] = sort(currentTrialTimes(1:48));
        currentApproachVector = trlStruct(i).approach(1:48,:,1);
        sortedCorrectApproach = currentApproachVector(sortedIndex,1);
        sortedIncrectApproach = currentApproachVector(sortedIndex,2);
        RegcorApproachVector = [RegcorApproachVector sortedCorrectApproach];
        RegincApproachVector = [RegincApproachVector sortedIncrectApproach];
    elseif strcmp(masterTbl.SessionType{i},['Reversal' sessionNum]) && strcmp(masterTbl.Strain{i},strain)
        currentTrialTimes = trlStruct(i).trialTimes;
        [sortedTimes,sortedIndex] = sort(currentTrialTimes(1:48));
        currentApproachVector = trlStruct(i).approach(1:48,:,1);
        sortedCorrectApproach = currentApproachVector(sortedIndex,1);
        sortedIncrectApproach = currentApproachVector(sortedIndex,2);
        RevcorApproachVector = [RevcorApproachVector sortedCorrectApproach];
        RevincApproachVector = [RevincApproachVector sortedIncrectApproach];
    end
end



RegCorM = movmean(RegcorApproachVector,3); RegincorM = movmean(RegincApproachVector,3);
RevCorM = movmean(RevcorApproachVector,3); RevincorM = movmean(RevincApproachVector,3);

RegCorMSEM = std(RegCorM')/sqrt(size(RegCorM,2)); RegincorSEM = std(RegincorM')/sqrt(size(RegincorM,2));
RevCorMSEM = std(RevCorM')/sqrt(size(RevCorM,2)); RevincorSEM = std(RevincorM')/sqrt(size(RevincorM,2));

figure('Units','normalized','Position',[0 0 1 1])
plot(mean(RegCorM,2),'r-o','LineWidth',3); hold on; plot(mean(RevCorM,2),'b-o','LineWidth',3);
er = errorbar([mean(RegCorM,2) mean(RevCorM,2)],[RegCorMSEM' RevCorMSEM'],'LineWidth',3); 
er(1).Color = 'r'; er(2).Color = 'b';                          
er(1).LineStyle = 'none'; er(2).LineStyle = 'none'; 
ylim([0 1.1])
xlabel('Trial')
ylabel('Approach Likelihood')
title(['Correct Approach Likelihood for Session ' sessionNum ' of ' strain])
legend([{'Congruent'},{'Incongruent'}])

set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
saveas(gca,['approachLikelihood_Correct' sessionNum strain],'png')

figure('Units','normalized','Position',[0 0 1 1])
plot(mean(RegincorM,2),'r-o','LineWidth',3); hold on; plot(mean(RevincorM,2),'b-o','LineWidth',3);
er = errorbar([mean(RegincorM,2) mean(RevincorM,2)],[RegincorSEM' RevincorSEM'],'LineWidth',3);  
er(1).Color = 'r'; er(2).Color = 'b';                          
er(1).LineStyle = 'none'; er(2).LineStyle = 'none'; 
ylim([0 1.1])
xlabel('Trial')
ylabel('Approach Likelihood')
title(['Incorrect Approach Likelihood for Session ' sessionNum ' of ' strain])
legend([{'Congruent'},{'Incongruent'}])

set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
saveas(gca,['approachLikelihood_Incorrect' sessionNum strain],'png')

% Ratio Between Correct and Incorrect


RegRatioIncorrect = RegCorM - RegincorM;
RevRatioIncorrect = RevCorM - RevincorM;

RegRatioSEM = std(RegRatioIncorrect')/sqrt(size(RegRatioIncorrect,2));
RevRatioSEM = std(RevRatioIncorrect')/sqrt(size(RevRatioIncorrect,2));

figure('Units','normalized','Position',[0 0 1 1])
plot(mean(RegRatioIncorrect,2),'r-o','LineWidth',3); hold on; plot(mean(RevRatioIncorrect,2),'b-o','LineWidth',3);
er = errorbar([mean(RegRatioIncorrect,2) mean(RevRatioIncorrect,2)],[RegRatioSEM' RevRatioSEM'],'LineWidth',3); 
er(1).Color = 'r'; er(2).Color = 'b';                          
er(1).LineStyle = 'none'; er(2).LineStyle = 'none'; 
ylim([-1 1])
xlabel('Trial')
ylabel('\Delta Likelihood')
title(['Difference of Correct & Incorrect Approaches for Session ' sessionNum ' of ' strain])
legend([{'Congruent'},{'Incongruent'}])

set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
saveas(gca,['llDifferenceCorrect' sessionNum strain],'png')

%% Correct and Incorrect Approaches Combined

close all

sessionNum = num2str([]);

strain = 'W';

RegcorApproachVector = []; RevcorApproachVector = [];
RegincApproachVector = []; RevincApproachVector = [];

for i = 1:length(dataSetIDs)
    if strcmp(masterTbl.SessionType{i},['Regular' sessionNum]) && strcmp(masterTbl.Strain{i},strain)
        currentTrialTimes = trlStruct(i).trialTimes;
        [sortedTimes,sortedIndex] = sort(currentTrialTimes(1:48));
        currentApproachVector = trlStruct(i).approach(1:48,:,1);
        sortedCorrectApproach = currentApproachVector(sortedIndex,1);
        sortedIncrectApproach = currentApproachVector(sortedIndex,2);
        RegcorApproachVector = [RegcorApproachVector sortedCorrectApproach];
        RegincApproachVector = [RegincApproachVector sortedIncrectApproach];
    elseif strcmp(masterTbl.SessionType{i},['Reversal' sessionNum]) && strcmp(masterTbl.Strain{i},strain)
        currentTrialTimes = trlStruct(i).trialTimes;
        [sortedTimes,sortedIndex] = sort(currentTrialTimes(1:48));
        currentApproachVector = trlStruct(i).approach(1:48,:,1);
        sortedCorrectApproach = currentApproachVector(sortedIndex,1);
        sortedIncrectApproach = currentApproachVector(sortedIndex,2);
        RevcorApproachVector = [RevcorApproachVector sortedCorrectApproach];
        RevincApproachVector = [RevincApproachVector sortedIncrectApproach];
    end
end

RegApproachVector = RegcorApproachVector > 0 | RegincApproachVector > 0;
RevApproachVecotr = RevcorApproachVector > 0 | RevincApproachVector > 0;


RegApp = movmean(RegcorApproachVector,3); 
RevApp = movmean(RevcorApproachVector,3); 

RegAppSEM = std(RegApp')/sqrt(size(RegApp,2)); 
RevAppSEM = std(RevApp')/sqrt(size(RevApp,2)); 

figure('Units','normalized','Position',[0 0 1 1])
plot(mean(RegApp,2),'r-o','LineWidth',3); hold on; plot(mean(RevApp,2),'b-o','LineWidth',3);
er = errorbar([mean(RegApp,2) mean(RevApp,2)],[RegAppSEM' RevAppSEM'],'LineWidth',3); 
er(1).Color = 'r'; er(2).Color = 'b';                          
er(1).LineStyle = 'none'; er(2).LineStyle = 'none'; 
ylim([0 1.1])
xlabel('Trial')
ylabel('Approach Likelihood')
title(['Approach Likelihood, All Trials, for Session ' sessionNum ' of ' strain])
legend([{'Congruent'},{'Incongruent'}])

set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
saveas(gca,['approachLikelihood_ALLAPP' sessionNum strain],'png')


%% Cue Discrimination per Session
close all

sessionNum = num2str([]);

strain = 'W';

RegcorCue = []; RevcorCue = [];  RegCSMinCue = [];
RegincCue = []; RevincCue = [];  RevCSMinCue = [];  

for i = 1:length(dataSetIDs)
    if strcmp(masterTbl.SessionType{i},['Regular' sessionNum]) && strcmp(masterTbl.Strain{i},strain)
        RegcorCue = [RegcorCue sum(trlStruct(i).approach(1:48,3,1))];
        RegincCue = [RegincCue sum(trlStruct(i).approach(1:48,4,1))];
        RegCSMinCue = [RegCSMinCue sum(trlStruct(i).approach(49:end,4,1))];
    elseif strcmp(masterTbl.SessionType{i},['Reversal' sessionNum]) && strcmp(masterTbl.Strain{i},strain)
        RevcorCue = [RevcorCue sum(trlStruct(i).approach(1:48,4,1))];
        RevincCue = [RevincCue sum(trlStruct(i).approach(1:48,3,1))];
        RevCSMinCue = [RevCSMinCue sum(trlStruct(i).approach(49:end,4,1))];

    end
end

regCorMin = RegcorCue./RegCSMinCue;
regCorInc = RegcorCue./RegincCue;

revCorMin = RevcorCue./RevCSMinCue;
revCorInc = RevcorCue./RevincCue;

semregCorMin = std(regCorMin)/sqrt(length(regCorMin));
semregCorInc = std(regCorInc)/sqrt(length(regCorInc));

semrevCorMin = std(revCorMin)/sqrt(length(revCorMin));
semrevCorInc = std(revCorInc)/sqrt(length(revCorInc));

figure('Units','normalized','Position',[0 0 1 1])
hb = bar(categorical({'Congruent Sessions','Incongruent Sessions'}),[mean(regCorMin) mean(regCorInc); mean(revCorMin) mean(revCorInc)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);
hb(1).FaceColor = [0.5 0.5 0.5]; hb(2).FaceColor = [0.9 0.9 0.9];
offsetPos = [1+hb(1).XOffset 1+hb(2).XOffset; 2+hb(1).XOffset 2+hb(2).XOffset];
scatterXData = [linspace(offsetPos(1)-0.1,offsetPos(1)+0.1,length(regCorMin)),linspace(offsetPos(2)-0.1,offsetPos(2)+0.1,length(revCorMin)),linspace(offsetPos(3)-0.1,offsetPos(3)+0.1,length(revCorMin)),linspace(offsetPos(4)-0.1,offsetPos(4)+0.1,length(revCorInc))];
hold on
scatter(scatterXData',[regCorMin regCorInc revCorMin revCorInc],50,'ko','LineWidth',2)
er = errorbar(offsetPos,[mean(regCorMin) mean(regCorInc); mean(revCorMin) mean(revCorInc)],[semregCorMin semregCorInc; semrevCorMin semrevCorInc],'LineWidth',3);    
er(1).Color = [0 0 0]; er(2).Color = [0 0 0];                          
er(1).LineStyle = 'none'; er(2).LineStyle = 'none'; 
ylabel('Ratio (a.u.)')
legend([{'CS+/CS-'},{'Correct Cue / Incorrect Cue'}])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
% saveas(gca,['Latencies' sessionNum strain],'png')

%% Determine # Corrections per Session Type
close all

sessionNum = num2str([]);

strain = 'W';

regCorrections = []; regIncorrect = [];
revCorrections = []; revIncorrect = [];

for i = 1:length(dataSetIDs)
    vecCorrections = [];
    if strcmp(masterTbl.SessionType{i},['Regular' sessionNum]) && strcmp(masterTbl.Strain{i},strain)
        vecCorrections = trlStruct(i).approach(1:48,1:2,1);
        regCorrections = [regCorrections nansum(vecCorrections(:,1) == 1 & vecCorrections(:,2) == 1)];
        regIncorrect = [regIncorrect nansum(vecCorrections(:,2))];
    elseif strcmp(masterTbl.SessionType{i},['Reversal' sessionNum]) && strcmp(masterTbl.Strain{i},strain)
        vecCorrections = trlStruct(i).approach(1:48,1:2,1);
        revCorrections = [revCorrections nansum(vecCorrections(:,1) == 1 & vecCorrections(:,2) == 1)];
        revIncorrect = [revIncorrect nansum(vecCorrections(:,2))];

    end
end

% Number of Corrections
if isempty(sessionNum)
    sessionNum = num2str(1);
end
semregCorrections = std(regCorrections)/sqrt(length(regCorrections));
semrevCorrections = std(revCorrections)/sqrt(length(revCorrections));

figure('Units','normalized','Position',[0 0 1 1])
hb = bar(categorical({'Congruent Sessions','Incongruent Sessions'}),[mean(regCorrections); mean(revCorrections)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);
hb.FaceColor = 'flat'; 
colorData = hb.CData;
hb.CData(1,:) = [0.5 0.5 0.5]; hb.CData(2,:) = [0.9 0.9 0.9];

hold on
offsetPos = [1+hb(1).XOffset ; 2+hb(1).XOffset ];
plot([linspace(0.95,1.05,length(regCorrections));linspace(1.95,2.05,length(revCorrections))],[regCorrections; revCorrections],'k--o','LineWidth',3,'MarkerSize',10)
er = errorbar(offsetPos,[mean(regCorrections); mean(revCorrections)],[semregCorrections; semrevCorrections],'LineWidth',3);    
er.Color = [0 0 0];                           
er.LineStyle = 'none';

ylabel('Number of Corrections')
title(['Number of Corrections for Session ' sessionNum ' of ' strain])

set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
saveas(gca,['Corrections' sessionNum strain],'png')

% Number of Corrections Standardized to Incorrect Approaches
regStCorrections = regCorrections./regIncorrect;
revStCorrections = revCorrections./revIncorrect;


semregStCorrections = std(regStCorrections)/sqrt(length(regStCorrections));
semrevStCorrections = std(revStCorrections)/sqrt(length(revStCorrections));

figure('Units','normalized','Position',[0 0 1 1])
hb = bar(categorical({'Congruent Sessions','Incongruent Sessions'}),[mean(regStCorrections); mean(revStCorrections)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);
hb.FaceColor = 'flat'; 
colorData = hb.CData;
hb.CData(1,:) = [0.5 0.5 0.5]; hb.CData(2,:) = [0.9 0.9 0.9];

hold on

plot([linspace(0.95,1.05,length(regStCorrections));linspace(1.95,2.05,length(revStCorrections))],[regStCorrections; revStCorrections],'k--o','LineWidth',3,'MarkerSize',10)
er = errorbar(offsetPos,[mean(regStCorrections); mean(revStCorrections)],[semregStCorrections; semrevStCorrections],'LineWidth',3);    
er.Color = [0 0 0];                           
er.LineStyle = 'none';
ylim([-0.01 1.1])
ylabel('Proportion of Corrections')
title(['Proportion of Corrections for Session ' sessionNum ' of ' strain])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
saveas(gca,['PropCorrections' sessionNum strain],'png')
%% Combine All Sessions & Genotypes, Redo All Prior Analyses
%% Corrections (All Sessions & Genotypes)
close all

% sessionNum = num2str([]);
% 
% strain = 'P';

regCorrections = []; regIncorrect = [];
revCorrections = []; revIncorrect = [];

for i = 1:length(dataSetIDs)
    vecCorrections = [];
    if strcmp(masterTbl.SessionType{i},['Regular'])
        vecCorrections = trlStruct(i).approach(1:48,1:2,1);
        regCorrections = [regCorrections nansum(vecCorrections(:,1) == 1 & vecCorrections(:,2) == 1)];
        regIncorrect = [regIncorrect nansum(vecCorrections(:,2))];
    elseif strcmp(masterTbl.SessionType{i},['Reversal'])
        vecCorrections = trlStruct(i).approach(1:48,1:2,1);
        revCorrections = [revCorrections nansum(vecCorrections(:,1) == 1 & vecCorrections(:,2) == 1)];
        revIncorrect = [revIncorrect nansum(vecCorrections(:,2))];

    end
end

% Number of Corrections
% if isempty(sessionNum)
%     sessionNum = num2str(1);
% end
semregCorrections = std(regCorrections)/sqrt(length(regCorrections));
semrevCorrections = std(revCorrections)/sqrt(length(revCorrections));

figure('Units','normalized','Position',[0 0 1 1])
hb = bar(categorical({'Congruent Sessions','Incongruent Sessions'}),[mean(regCorrections); mean(revCorrections)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);
hb.FaceColor = 'flat'; 
colorData = hb.CData;
hb.CData(1,:) = [0.5 0.5 0.5]; hb.CData(2,:) = [0.9 0.9 0.9];

hold on
offsetPos = [1+hb(1).XOffset; 2+hb(1).XOffset];

plot([linspace(0.95,1.05,length(regCorrections));linspace(1.95,2.05,length(revCorrections))],[regCorrections; revCorrections],'k--o','LineWidth',3,'MarkerSize',10)
er = errorbar(offsetPos,[mean(regCorrections); mean(revCorrections)],[semregCorrections; semrevCorrections],'LineWidth',3);    
er.Color = [0 0 0];                           
er.LineStyle = 'none';

ylabel('Number of Corrections')
title(['Number of Corrections for All Sessions and Genotypes'])

set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
saveas(gca,['ALLCorrections'],'png')

% Number of Corrections Standardized to Incorrect Approaches
regStCorrections = regCorrections./regIncorrect;
revStCorrections = revCorrections./revIncorrect;


semregStCorrections = std(regStCorrections)/sqrt(length(regStCorrections));
semrevStCorrections = std(revStCorrections)/sqrt(length(revStCorrections));

figure('Units','normalized','Position',[0 0 1 1])
hb = bar(categorical({'Congruent Sessions','Incongruent Sessions'}),[mean(regStCorrections); mean(revStCorrections)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);
hb.FaceColor = 'flat'; 
colorData = hb.CData;
hb.CData(1,:) = [0.5 0.5 0.5]; hb.CData(2,:) = [0.9 0.9 0.9];

hold on

plot([linspace(0.95,1.05,length(regStCorrections));linspace(1.95,2.05,length(revStCorrections))],[regStCorrections; revStCorrections],'k--o','LineWidth',3,'MarkerSize',10)
er = errorbar(offsetPos,[mean(regStCorrections); mean(revStCorrections)],[semregStCorrections; semrevStCorrections],'LineWidth',3);    
er.Color = [0 0 0];                           
er.LineStyle = 'none';
ylim([-0.01 1.1])
ylabel('Proportion of Corrections')
title(['Proportion of Corrections for All Sessions and Genotypes'])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
saveas(gca,['ALLPropCorrections'],'png')

%% Approach Likelihood Combined Correct and Incorrect Approaches (All Sessions and Genotypes)

close all

sessionNum = num2str([]);

strain = 'W';

RegcorApproachVector = []; RevcorApproachVector = [];
RegincApproachVector = []; RevincApproachVector = [];

for i = 1:length(dataSetIDs)
    if strcmp(masterTbl.SessionType{i},['Regular'])
        currentTrialTimes = trlStruct(i).trialTimes;
        [sortedTimes,sortedIndex] = sort(currentTrialTimes(1:48));
        currentApproachVector = trlStruct(i).approach(1:48,:,1);
        sortedCorrectApproach = currentApproachVector(sortedIndex,1);
        sortedIncrectApproach = currentApproachVector(sortedIndex,2);
        RegcorApproachVector = [RegcorApproachVector sortedCorrectApproach];
        RegincApproachVector = [RegincApproachVector sortedIncrectApproach];
    elseif strcmp(masterTbl.SessionType{i},['Reversal'])
        currentTrialTimes = trlStruct(i).trialTimes;
        [sortedTimes,sortedIndex] = sort(currentTrialTimes(1:48));
        currentApproachVector = trlStruct(i).approach(1:48,:,1);
        sortedCorrectApproach = currentApproachVector(sortedIndex,1);
        sortedIncrectApproach = currentApproachVector(sortedIndex,2);
        RevcorApproachVector = [RevcorApproachVector sortedCorrectApproach];
        RevincApproachVector = [RevincApproachVector sortedIncrectApproach];
    end
end

RegApproachVector = RegcorApproachVector > 0 | RegincApproachVector > 0;
RevApproachVecotr = RevcorApproachVector > 0 | RevincApproachVector > 0;


RegApp = movmean(RegApproachVector,3); 
RevApp = movmean(RevApproachVecotr,3); 

RegAppSEM = std(RegApp')/sqrt(size(RegApp,2)); 
RevAppSEM = std(RevApp')/sqrt(size(RevApp,2)); 

figure('Units','normalized','Position',[0 0 1 1])
plot(mean(RegApp,2),'r-o','LineWidth',3); hold on; plot(mean(RevApp,2),'b-o','LineWidth',3);
er = errorbar([mean(RegApp,2) mean(RevApp,2)],[RegAppSEM' RevAppSEM'],'LineWidth',3); 
er(1).Color = 'r'; er(2).Color = 'b';                          
er(1).LineStyle = 'none'; er(2).LineStyle = 'none'; 
ylim([0 1.1])
xlabel('Trial')
ylabel('Approach Likelihood')
title(['Approach Likelihood, All Trials, All Sessions & Genotypes'])
legend([{'Congruent'},{'Incongruent'}])

set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
saveas(gca,['ALLapproachLikelihood_ALLAPP'],'png')

%% Find some summary information on Reg and Rev Trials (All Sessions and Genotypes)
addpath(genpath('F:/dissDat'))
% sessionNum = []; %num2str(3);
% strain = 'W';
close all 

regLDist = []; regRDist = []; revLDist = []; revRDist = [];


for i = 1:length(dataSetIDs)
    if strcmp(masterTbl.SessionType{i},['Regular' ]) 
        regLDist(:,i) = nanmedian(trlStruct(i).LDistCorrectAppr,1);
        regRDist(:,i) = nanmedian(trlStruct(i).RDistCorrectAppr,1);
    elseif strcmp(masterTbl.SessionType{i},['Reversal' ]) 
        revLDist(:,i) = nanmedian(trlStruct(i).LDistCorrectAppr,1);
        revRDist(:,i) = nanmedian(trlStruct(i).RDistCorrectAppr,1);
    end

end

% Remove 0's 
regLDist=regLDist(:,any(regLDist)); revLDist=revLDist(:,any(revLDist));
regRDist=regRDist(:,any(regRDist)); revRDist=revRDist(:,any(revRDist));
% SEM
regSEM = std([regLDist regRDist]')/sqrt(min(size(regLDist)));
revSEM = std([revLDist revRDist]')/sqrt(min(size(revLDist)));
% Test plots
timeAdj = 1:length(regLDist);
timeAdj = timeAdj./30;
timeAdj = timeAdj - 5;
figure('Units','normalized','Position',[0 0 1 1])
sfh1 = subplot(2,1,1);
shadedErrorBar(timeAdj,nanmean([regLDist regRDist],2),regSEM,'lineprops',{'r-','LineWidth',3});
hold on
shadedErrorBar(timeAdj,nanmean([revLDist revRDist],2),revSEM,'lineprops',{'b-','LineWidth',3});
xlabel('Time from Cue Onset (s)')
ylabel('Distance from Correct Sipper')
xline(0,'k--','Cue On','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(5,'k--','Sipper In','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(13,'k--','Sipper Out','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
title('Correct Approach Trials (All Sessions and Genotypes)')
legend([{'Congruent'},{'Incongruent'}])
xlim([-5 20])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

sfh2 = subplot(2,1,2);
KLVal = KLDiv([regLDist regRDist],[revLDist revRDist]);
plot(timeAdj,KLVal,'k-','LineWidth',3)
ylabel('KL Divergence')
xlabel('Time (s)')
xline(0,'k--','Cue On','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(5,'k--','Sipper In','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(13,'k--','Sipper Out','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
ylim([0 2])
xlim([-5 20])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

sfh2.Position = [sfh2.Position(1) sfh2.Position(2) sfh2.Position(3) sfh2.Position(4)/2];
sfh1.Position = [sfh1.Position(1) sfh1.Position(2)-sfh2.Position(4) sfh1.Position(3) sfh1.Position(4)*1.5];

saveas(gca,['ALLcorrectApproachDist_withFRINGE'],'png')

% Determine Correct Approaches per Session Type
regCorrect = []; revCorrect = [];

for i = 1:length(dataSetIDs)
    if strcmp(masterTbl.SessionType{i},['Regular']) 
        regCorrect = [regCorrect size(trlStruct(i).LDistCorrectAppr,1) + size(trlStruct(i).RDistCorrectAppr,1)];
    elseif strcmp(masterTbl.SessionType{i},['Reversal']) 
        revCorrect = [revCorrect size(trlStruct(i).LDistCorrectAppr,1) + size(trlStruct(i).RDistCorrectAppr,1)];
    end

end
semregCR = std(regCorrect)/sqrt(length(regCorrect));
semrevCR = std(revCorrect)/sqrt(length(revCorrect));

figure('Units','normalized','Position',[0 0 1 1])
bar(categorical({'Congruent Correct Approach','Incongruent Correct Approach'}),[mean(regCorrect) mean(revCorrect)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3)
hold on
plot([linspace(0.95,1.05,length(regCorrect));linspace(1.95,2.05,length(revCorrect))],[regCorrect; revCorrect],'k--o','LineWidth',3,'MarkerSize',10)
er = errorbar(categorical({'Congruent Correct Approach','Incongruent Correct Approach'}),[mean(regCorrect) mean(revCorrect)],[semregCR semrevCR],'LineWidth',3);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
ylabel('Trials')
ylim([0 48])

set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,['ALLcorrectApproachNum'],'png')



% Line incorrect
regLDist = []; regRDist = []; revLDist = []; revRDist = [];


for i = 1:length(dataSetIDs)
    if strcmp(masterTbl.SessionType{i},['Regular'])
        regLDist(:,i) = nanmedian(trlStruct(i).LDistInCorrectAppr,1);
        regRDist(:,i) = nanmedian(trlStruct(i).RDistInCorrectAppr,1);
    elseif strcmp(masterTbl.SessionType{i},['Reversal'])
        revLDist(:,i) = nanmedian(trlStruct(i).LDistInCorrectAppr,1);
        revRDist(:,i) = nanmedian(trlStruct(i).RDistInCorrectAppr,1);
    end

end

% Remove 0's 
regLDist=regLDist(:,any(regLDist)); revLDist=revLDist(:,any(revLDist));
regRDist=regRDist(:,any(regRDist)); revRDist=revRDist(:,any(revRDist));

% SEM
regSEM = std([regLDist regRDist]')/sqrt(min(size(regLDist)));
revSEM = std([revLDist revRDist]')/sqrt(min(size(revLDist)));
% Test plots
timeAdj = 1:length(regLDist);
timeAdj = timeAdj./30;
timeAdj = timeAdj - 5;
figure('Units','normalized','Position',[0 0 1 1])
sfh1 = subplot(2,1,1);
shadedErrorBar(timeAdj,nanmean([regLDist regRDist],2),regSEM,'lineprops',{'r-','LineWidth',3});
hold on
shadedErrorBar(timeAdj,nanmean([revLDist revRDist],2),revSEM,'lineprops',{'b-','LineWidth',3});
xlabel('Time from Cue Onset (s)')
ylabel('Distance from Correct Sipper')
xline(0,'k--','Cue On','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(5,'k--','Sipper In','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(13,'k--','Sipper Out','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
title('Incorrect Approach Trials (All Sessions and Genotypes)')
legend([{'Congruent'},{'Incongruent'}])
xlim([-5 20])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)


sfh2 = subplot(2,1,2);
KLVal = KLDiv([regLDist regRDist],[revLDist revRDist]);
plot(timeAdj,KLVal,'k-','LineWidth',3)
ylabel('KL Divergence')
xlabel('Time (s)')
xline(0,'k--','Cue On','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(5,'k--','Sipper In','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(13,'k--','Sipper Out','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xlim([-5 20])
ylim([0 2])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

sfh2.Position = [sfh2.Position(1) sfh2.Position(2) sfh2.Position(3) sfh2.Position(4)/2];
sfh1.Position = [sfh1.Position(1) sfh1.Position(2)-sfh2.Position(4) sfh1.Position(3) sfh1.Position(4)*1.5];

saveas(gca,['ALLincorrectApproachDist_withFRINGE'],'png')

% Determine incorrect Approaches per Session Type
regInCorrect = []; revInCorrect = [];

for i = 1:length(dataSetIDs)
    if strcmp(masterTbl.SessionType{i},['Regular'])
        regInCorrect = [regInCorrect size(trlStruct(i).LDistInCorrectAppr,1) + size(trlStruct(i).RDistInCorrectAppr,1)];
    elseif strcmp(masterTbl.SessionType{i},['Reversal'])
        revInCorrect = [revInCorrect size(trlStruct(i).LDistInCorrectAppr,1) + size(trlStruct(i).RDistInCorrectAppr,1)];
    end

end

semregIN = std(regInCorrect)/sqrt(length(regInCorrect));
semrevIN = std(revInCorrect)/sqrt(length(revInCorrect));

figure('Units','normalized','Position',[0 0 1 1])
bar(categorical({'Congruent Incorrect Approach','Incongruent Incorrect Approach'}),[mean(regInCorrect) mean(revInCorrect)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3)
hold on
plot([linspace(0.95,1.05,length(regInCorrect));linspace(1.95,2.05,length(revInCorrect))],[regInCorrect; revInCorrect],'k--o','LineWidth',3,'MarkerSize',10)
er = errorbar(categorical({'Congruent Incorrect Approach','Incongruent Incorrect Approach'}),[mean(regInCorrect) mean(revInCorrect)],[semregIN semrevIN],'LineWidth',3);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
ylabel('Trials')
ylim([0 48])

set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,['ALLincorrectApproachNum'],'png')

% Find omissions
regOmit = []; revOmit = [];

for i = 1:length(dataSetIDs)
    if strcmp(masterTbl.SessionType{i},['Regular'])
        regOmit = [regOmit sum(sum([trlStruct(i).approach(1:48,1:2,1) == 1],2)==0)];
    elseif strcmp(masterTbl.SessionType{i},['Reversal'])
        revOmit = [revOmit sum(sum([trlStruct(i).approach(1:48,1:2,1) == 1],2)==0)];
    end

end
semregOmit = std(regOmit)/sqrt(length(regOmit));
semrevcorLat = std(revOmit)/sqrt(length(revOmit));

figure('Units','normalized','Position',[0 0 1 1])
bar(categorical({'Congruent Omissions','Incongruent Omissions'}),[mean(regOmit) mean(revOmit)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3)
hold on
plot([linspace(0.95,1.05,length(regOmit));linspace(1.95,2.05,length(revOmit))],[regOmit; revOmit],'k--o','LineWidth',3,'MarkerSize',10)
er = errorbar(categorical({'Congruent Omissions','Incongruent Omissions'}),[mean(regOmit) mean(revOmit)],[semregOmit semrevcorLat],'LineWidth',3);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
ylabel('Trials')
ylim([0 48])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,['ALLomitNum'],'png')

[h,p,stat] = ttest(regOmit,revOmit);

%% Determine secondary statistics such as latency
% Find first moment following sipper descent where animal makes contact
% with sipper pixels. Use threshold set by Nick's code (< 9 pixels) and
% threshold for timebins (> 3 timebins). First start by indexing all
% moments where distance to correct sipper is less than 9 pixels, then
% ensure that the # of pixels is greater than 3. The first index where this
% occurs will be the latency to approach the sipper. 
disThresh = 9;
disRun = 3;
sipDescent = 10*30;
sipAscent = 18*30;
cueOn = 5*30;
ReglatCor = []; ReglatinCor = []; RevlatCor = []; RevlatinCor = []; % Correct and Incorrect latencies

for i = 1:length(dataSetIDs)
    if regexp(masterTbl.SessionType{i},regexptranslate('wildcard','Regular*'))
        latCorIdx = [trlStruct(i).RDistCorrectAppr <= disThresh; trlStruct(i).LDistCorrectAppr <= disThresh];
        latInCorIdx = [trlStruct(i).RDistInCorrectAppr <= disThresh; trlStruct(i).LDistInCorrectAppr <= disThresh];
        [~,sessCorLats] = max(latCorIdx(:,sipDescent:sipAscent),[],2);
        [~,sessInCorLats] = max(latInCorIdx(:,sipDescent:sipAscent),[],2);
        ReglatCor = [ReglatCor; sessCorLats/30;];
        ReglatinCor = [ReglatinCor; sessInCorLats/30;];
    elseif regexp(masterTbl.SessionType{i},regexptranslate('wildcard','Reversal*'))
        latCorIdx = [trlStruct(i).RDistCorrectAppr <= disThresh; trlStruct(i).LDistCorrectAppr <= disThresh];
        latInCorIdx = [trlStruct(i).RDistInCorrectAppr <= disThresh; trlStruct(i).LDistInCorrectAppr <= disThresh];
        [~,sessCorLats] = max(latCorIdx(:,sipDescent:sipAscent),[],2);
        [~,sessInCorLats] = max(latInCorIdx(:,sipDescent:sipAscent),[],2);
        RevlatCor = [RevlatCor; sessCorLats/30;];
        RevlatinCor = [RevlatinCor; sessInCorLats/30;];
    end
end

% Create latency figure
semregcorLat = std(ReglatCor)/sqrt(length(ReglatCor));
semrevcorLat = std(RevlatCor)/sqrt(length(RevlatCor));
semregIncorLat = std(ReglatinCor)/sqrt(length(ReglatinCor));
semrevIncorLat = std(RevlatinCor)/sqrt(length(RevlatinCor));

figure('Units','normalized','Position',[0 0 1 1])
hb = bar(categorical({'Correct Latencies','Incorrect Latencies'}),[mean(ReglatCor) mean(RevlatCor); mean(ReglatinCor) mean(RevlatinCor)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);
hb(1).FaceColor = [0.5 0.5 0.5]; hb(2).FaceColor = [0.9 0.9 0.9];
offsetPos = [1+hb(1).XOffset 1+hb(2).XOffset 2+hb(1).XOffset 2+hb(2).XOffset];
scatterXData = [linspace(offsetPos(1)-0.1,offsetPos(1)+0.1,length(ReglatCor)),linspace(offsetPos(2)-0.1,offsetPos(2)+0.1,length(RevlatCor)),linspace(offsetPos(3)-0.1,offsetPos(3)+0.1,length(ReglatinCor)),linspace(offsetPos(4)-0.1,offsetPos(4)+0.1,length(RevlatinCor))];
hold on
scatter(scatterXData,[ReglatCor' RevlatCor' ReglatinCor' RevlatinCor'],50,'ko','LineWidth',2)
er = errorbar(offsetPos,[mean(ReglatCor) mean(RevlatCor) mean(ReglatinCor) mean(RevlatinCor)],[semregcorLat semrevcorLat semregIncorLat semrevIncorLat],'LineWidth',3);    
er(1).Color = [0 0 0];                        
er(1).LineStyle = 'none'; 
ylabel('Time (s)')
legend([{'Congruent'},{'Incongruent'}])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
saveas(gca,['ALLLatencies'],'png')
% xData = hb(1).XData+hb(1).XOffset;

%% Histogram of Latency Data
figure('Units','normalized','Position',[0 0 1 1])
histogram(ReglatCor,14,'FaceColor','r','Normalization','pdf'); hold on; histogram(RevlatCor,14,'FaceColor','b','Normalization','pdf')
title('Histogram of Correct Approach Trials | Congruent & Incongruent')
legend([{'Congruent'},{'Incongruent'}])
xlabel('Latency (s)')
ylabel('Counts')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)

figure('Units','normalized','Position',[0 0 1 1])
histogram(ReglatinCor,14,'FaceColor','r'); hold on; histogram(RevlatinCor,14,'FaceColor','b')
title('Histogram of Incorrect Approach Trials | Congruent & Incongruent')
legend([{'Congruent'},{'Incongruent'}])
xlabel('Latency (s)')
ylabel('Counts')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)

%% Fit distributions to latency data
xwin = [0:1/30:max(ReglatCor)];
regLatFit = fitdist(ReglatCor,'Kernel');   regLatFitInc = fitdist(ReglatinCor,'Kernel');  
revLatFit = fitdist(RevlatCor,'Kernel');  revLatFitInc = fitdist(RevlatinCor,'Kernel');

regLatDist = pdf(regLatFit,xwin);    regLatDistInc = pdf(regLatFitInc,xwin); 
revLatDist = pdf(revLatFit,xwin);   revLatDistInc = pdf(revLatFitInc,xwin); 

figure('Units','normalized','Position',[0 0 1 1])
plot(xwin,regLatDist,'r-','LineWidth',3); hold on; plot(xwin,revLatDist,'b-','LineWidth',3);
title('Probability Distribution of Latencies in Correct Approach Trials')
xlabel('Latency to Sipper (s)')
ylabel('Probability')
legend([{'Congruent'},{'Incongruent'}])
%ylim([0 1])
% xlim([min(xwin) max(xwin)])
% plot(sort(ReglatCor),1:length(ReglatCor),'ro')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)

figure('Units','normalized','Position',[0 0 1 1])
plot(xwin,regLatDistInc,'r-','LineWidth',3); hold on; plot(xwin,revLatDistInc,'b-','LineWidth',3);
title('Probability Distribution of Latencies in Incorrect Approach Trials')
xlabel('Latency to Sipper (s)')
ylabel('Probability')
%ylim([0 1])
% xlim([min(xwin) max(xwin)])
legend([{'Congruent'},{'Incongruent'}])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)

%% Approach Likelihood, Metrics over Trials, Correct or Incorrect
%sessionNum = []; 
close all

sessionNum = num2str([]);

strain = 'P';

RegcorApproachVector = []; RevcorApproachVector = [];
RegincApproachVector = []; RevincApproachVector = [];
for i = 1:length(dataSetIDs)
    if strcmp(masterTbl.SessionType{i},['Regular'])
        currentTrialTimes = trlStruct(i).trialTimes;
        [sortedTimes,sortedIndex] = sort(currentTrialTimes(1:48));
        currentApproachVector = trlStruct(i).approach(1:48,:,1);
        sortedCorrectApproach = currentApproachVector(sortedIndex,1);
        sortedIncrectApproach = currentApproachVector(sortedIndex,2);
        RegcorApproachVector = [RegcorApproachVector sortedCorrectApproach];
        RegincApproachVector = [RegincApproachVector sortedIncrectApproach];
    elseif strcmp(masterTbl.SessionType{i},['Reversal'])
        currentTrialTimes = trlStruct(i).trialTimes;
        [sortedTimes,sortedIndex] = sort(currentTrialTimes(1:48));
        currentApproachVector = trlStruct(i).approach(1:48,:,1);
        sortedCorrectApproach = currentApproachVector(sortedIndex,1);
        sortedIncrectApproach = currentApproachVector(sortedIndex,2);
        RevcorApproachVector = [RevcorApproachVector sortedCorrectApproach];
        RevincApproachVector = [RevincApproachVector sortedIncrectApproach];
    end
end



RegCorM = movmean(RegcorApproachVector,3); RegincorM = movmean(RegincApproachVector,3);
RevCorM = movmean(RevcorApproachVector,3); RevincorM = movmean(RevincApproachVector,3);

RegCorMSEM = std(RegCorM')/sqrt(size(RegCorM,2)); RegincorSEM = std(RegincorM')/sqrt(size(RegincorM,2));
RevCorMSEM = std(RevCorM')/sqrt(size(RevCorM,2)); RevincorSEM = std(RevincorM')/sqrt(size(RevincorM,2));

figure('Units','normalized','Position',[0 0 1 1])
plot(mean(RegCorM,2),'r-o','LineWidth',3); hold on; plot(mean(RevCorM,2),'b-o','LineWidth',3);
er = errorbar([mean(RegCorM,2) mean(RevCorM,2)],[RegCorMSEM' RevCorMSEM'],'LineWidth',3); 
er(1).Color = 'r'; er(2).Color = 'b';                          
er(1).LineStyle = 'none'; er(2).LineStyle = 'none'; 
ylim([0 1.1])
xlabel('Trial')
ylabel('Approach Likelihood')
title(['Correct Approach Likelihood for for All Sessions & Genotypes'])
legend([{'Congruent'},{'Incongruent'}])

set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
saveas(gca,['ALLapproachLikelihood_Correct'],'png')

figure('Units','normalized','Position',[0 0 1 1])
plot(mean(RegincorM,2),'r-o','LineWidth',3); hold on; plot(mean(RevincorM,2),'b-o','LineWidth',3);
er = errorbar([mean(RegincorM,2) mean(RevincorM,2)],[RegincorSEM' RevincorSEM'],'LineWidth',3);  
er(1).Color = 'r'; er(2).Color = 'b';                          
er(1).LineStyle = 'none'; er(2).LineStyle = 'none'; 
ylim([0 1.1])
xlabel('Trial')
ylabel('Approach Likelihood')
title(['Incorrect Approach Likelihood for All Sessions & Genotypes'])
legend([{'Congruent'},{'Incongruent'}])

set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
saveas(gca,['ALLapproachLikelihood_Incorrect'],'png')

% Ratio Between Correct and Incorrect


RegRatioIncorrect = RegCorM - RegincorM;
RevRatioIncorrect = RevCorM - RevincorM;

RegRatioSEM = std(RegRatioIncorrect')/sqrt(size(RegRatioIncorrect,2));
RevRatioSEM = std(RevRatioIncorrect')/sqrt(size(RevRatioIncorrect,2));

figure('Units','normalized','Position',[0 0 1 1])
plot(mean(RegRatioIncorrect,2),'r-o','LineWidth',3); hold on; plot(mean(RevRatioIncorrect,2),'b-o','LineWidth',3);
er = errorbar([mean(RegRatioIncorrect,2) mean(RevRatioIncorrect,2)],[RegRatioSEM' RevRatioSEM'],'LineWidth',3); 
er(1).Color = 'r'; er(2).Color = 'b';                          
er(1).LineStyle = 'none'; er(2).LineStyle = 'none'; 
ylim([-1 1])
xlabel('Trial')
ylabel('\Delta Likelihood')
title(['Difference of Correct & Incorrect Approaches for All Sessions & Genotypes'])
legend([{'Congruent'},{'Incongruent'}])

set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
saveas(gca,['ALLllDifferenceCorrect'],'png')

%% Likelihood to Occupy Sipper
disThresh = 9;
disRun = 3;
sipDescent = 10*30;
sipAscent = 18*30;
cueOn = 5*30;
llReg = []; llRev = []; % Correct and Incorrect latencies

for i = 1:length(dataSetIDs)
    if regexp(masterTbl.SessionType{i},regexptranslate('wildcard','Regular*'))
        CorIdx = [trlStruct(i).RDistCorrectAppr <= disThresh; trlStruct(i).LDistCorrectAppr <= disThresh];
        CorDbl = double(CorIdx);
        regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
        llOccupy = movmean(CorDbl',3);
        trlStruct(i).llOccupy = llOccupy;
        trlStruct(i).CorCorOcc = regOcc;
    elseif regexp(masterTbl.SessionType{i},regexptranslate('wildcard','Reversal*'))
        CorIdx = [trlStruct(i).RDistCorrectAppr <= disThresh; trlStruct(i).LDistCorrectAppr <= disThresh];
        CorDbl = double(CorIdx);
        regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
        llOccupy = movmean(CorDbl',3);
        trlStruct(i).llOccupy = llOccupy;
        trlStruct(i).CorCorOcc = regOcc;
    end
end

for i = 1:length(dataSetIDs)
    if regexp(masterTbl.SessionType{i},regexptranslate('wildcard','Regular*'))
        llReg = [llReg trlStruct(i).llOccupy];
    elseif regexp(masterTbl.SessionType{i},regexptranslate('wildcard','Reversal*'))
        llRev = [llRev trlStruct(i).llOccupy];
    end
end

% Plot 
regSEM = std(llReg')/sqrt(min(size(llReg)));
revSEM = std(llRev')/sqrt(min(size(llRev)));
% Test plots
timeAdj = 1:length(llRev);
timeAdj = timeAdj./30;
timeAdj = timeAdj - 5;
figure('Units','normalized','Position',[0 0 1 1])
shadedErrorBar(timeAdj,nanmean(llReg,2),regSEM,'lineprops',{'r-','LineWidth',3});
hold on
shadedErrorBar(timeAdj,nanmean(llRev,2),revSEM,'lineprops',{'b-','LineWidth',3});
xlabel('Time from Cue Onset (s)')
ylabel('Likelihood to Occupy Correct Sipper')
xline(0,'k--','Cue On','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(5,'k--','Sipper In','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(13,'k--','Sipper Out','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
title('Likelihood to Occupy Correct Sipper on Correct Approach Trials (All Sessions and Genotypes)')
legend([{'Congruent'},{'Incongruent'}])
xlim([-5 20])
ylim([0 1])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

%% Likelihood to Occupy Correct Sipper on Incorrect Approach Trials
disThresh = 9;
disRun = 3;
sipDescent = 10*30;
sipAscent = 18*30;
cueOn = 5*30;
llReg = []; llRev = []; % Correct and Incorrect latencies

for i = 1:length(dataSetIDs)
    if regexp(masterTbl.SessionType{i},regexptranslate('wildcard','Regular*'))
        CorIdx = [trlStruct(i).RDistInCorrectAppr <= disThresh; trlStruct(i).LDistInCorrectAppr <= disThresh];
        CorDbl = double(CorIdx);
        regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
        llOccupy = movmean(CorDbl',3);
        trlStruct(i).CorIncOcc = regOcc;
        trlStruct(i).llOccupy = llOccupy;
    elseif regexp(masterTbl.SessionType{i},regexptranslate('wildcard','Reversal*'))
        CorIdx = [trlStruct(i).RDistInCorrectAppr <= disThresh; trlStruct(i).LDistInCorrectAppr <= disThresh];
        CorDbl = double(CorIdx);
        regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
        llOccupy = movmean(CorDbl',3);
        trlStruct(i).CorIncOcc = regOcc;
        trlStruct(i).llOccupy = llOccupy;
    end
end

for i = 1:length(dataSetIDs)
    if regexp(masterTbl.SessionType{i},regexptranslate('wildcard','Regular*'))
        llReg = [llReg trlStruct(i).llOccupy];
    elseif regexp(masterTbl.SessionType{i},regexptranslate('wildcard','Reversal*'))
        llRev = [llRev trlStruct(i).llOccupy];
    end
end

% Plot 
regSEM = std(llReg')/sqrt(min(size(llReg)));
revSEM = std(llRev')/sqrt(min(size(llRev)));
% Test plots
timeAdj = 1:length(llRev);
timeAdj = timeAdj./30;
timeAdj = timeAdj - 5;
figure('Units','normalized','Position',[0 0 1 1])
shadedErrorBar(timeAdj,nanmean(llReg,2),regSEM,'lineprops',{'r-','LineWidth',3});
hold on
shadedErrorBar(timeAdj,nanmean(llRev,2),revSEM,'lineprops',{'b-','LineWidth',3});
xlabel('Time from Cue Onset (s)')
ylabel('Likelihood to Occupy Correct Sipper')
xline(0,'k--','Cue On','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(5,'k--','Sipper In','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(13,'k--','Sipper Out','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
title('Likelihood to Occupy Correct Sipper on Incorrect Approach Trials (All Sessions and Genotypes)')
legend([{'Congruent'},{'Incongruent'}])
xlim([-5 20])
ylim([0 1])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

%% LL Incorrect Sipper Occupancy on Correct Approach Trials
disThresh = 9;
disRun = 3;
sipDescent = 10*30;
sipAscent = 18*30;
cueOn = 5*30;
llReg = []; llRev = []; % Correct and Incorrect latencies

for i = 1:length(dataSetIDs)
    if regexp(masterTbl.SessionType{i},regexptranslate('wildcard','Regular*'))
        CorIdx = [trlStruct(i).trlRSipDist(LeftTrials(trlStruct(i).LcorrectIdx),:) <= disThresh; trlStruct(i).trlLSipDist(RightTrials(trlStruct(i).RcorrectIdx),:) <= disThresh];
        CorDbl = double(CorIdx);
        llOccupy = movmean(CorDbl',3);
        regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
        trlStruct(i).IncCorOcc = regOcc;
        trlStruct(i).llOccupy = llOccupy;
    elseif regexp(masterTbl.SessionType{i},regexptranslate('wildcard','Reversal*'))
        CorIdx = [trlStruct(i).trlRSipDist(RightTrials(trlStruct(i).RcorrectIdx),:) <= disThresh; trlStruct(i).trlLSipDist(LeftTrials(trlStruct(i).LcorrectIdx),:) <= disThresh];
        CorDbl = double(CorIdx);
        llOccupy = movmean(CorDbl',3);
        regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
        trlStruct(i).IncCorOcc = regOcc;
        trlStruct(i).llOccupy = llOccupy;
    end
end

for i = 1:length(dataSetIDs)
    if regexp(masterTbl.SessionType{i},regexptranslate('wildcard','Regular*'))
        llReg = [llReg trlStruct(i).llOccupy];
    elseif regexp(masterTbl.SessionType{i},regexptranslate('wildcard','Reversal*'))
        llRev = [llRev trlStruct(i).llOccupy];
    end
end

% Plot 
regSEM = std(llReg')/sqrt(min(size(llReg)));
revSEM = std(llRev')/sqrt(min(size(llRev)));
% Test plots
timeAdj = 1:length(regSEM);
timeAdj = timeAdj./30;
timeAdj = timeAdj - 5;
figure('Units','normalized','Position',[0 0 1 1])
shadedErrorBar(timeAdj,nanmean(llReg,2),regSEM,'lineprops',{'r-','LineWidth',3});
hold on
shadedErrorBar(timeAdj,nanmean(llRev,2),revSEM,'lineprops',{'b-','LineWidth',3});
xlabel('Time from Cue Onset (s)')
ylabel('Likelihood to Occupy Incorrect Sipper')
xline(0,'k--','Cue On','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(5,'k--','Sipper In','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(13,'k--','Sipper Out','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
title('Likelihood to Occupy Incorrect Sipper on Correct Approach Trials (All Sessions and Genotypes)')
legend([{'Congruent'},{'Incongruent'}])
xlim([-5 20])
ylim([0 1])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
%% LL Incorrect Sipper on Incorrect Approach Trials
disThresh = 9;
disRun = 3;
sipDescent = 10*30;
sipAscent = 18*30;
cueOn = 5*30;
llReg = []; llRev = []; % Correct and Incorrect latencies

for i = 1:length(dataSetIDs)
    if regexp(masterTbl.SessionType{i},regexptranslate('wildcard','Regular*'))
        CorIdx = [trlStruct(i).trlRSipDist(LeftTrials(trlStruct(i).LincorrectIdx),:) <= disThresh; trlStruct(i).trlLSipDist(RightTrials(trlStruct(i).RincorrectIdx),:) <= disThresh];
        CorDbl = double(CorIdx);
        llOccupy = movmean(CorDbl',3);
        regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
        trlStruct(i).IncIncOcc = regOcc;
        trlStruct(i).llOccupy = llOccupy;
    elseif regexp(masterTbl.SessionType{i},regexptranslate('wildcard','Reversal*'))
        CorIdx = [trlStruct(i).trlRSipDist(RightTrials(trlStruct(i).RincorrectIdx),:) <= disThresh; trlStruct(i).trlLSipDist(LeftTrials(trlStruct(i).LincorrectIdx),:) <= disThresh];
        CorDbl = double(CorIdx);
        llOccupy = movmean(CorDbl',3);
        regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
        trlStruct(i).IncIncOcc = regOcc;
        trlStruct(i).llOccupy = llOccupy;
    end
end

for i = 1:length(dataSetIDs)
    if regexp(masterTbl.SessionType{i},regexptranslate('wildcard','Regular*'))
        llReg = [llReg trlStruct(i).llOccupy];
    elseif regexp(masterTbl.SessionType{i},regexptranslate('wildcard','Reversal*'))
        llRev = [llRev trlStruct(i).llOccupy];
    end
end

% Plot 
regSEM = std(llReg')/sqrt(min(size(llReg)));
revSEM = std(llRev')/sqrt(min(size(llRev)));
% Test plots
timeAdj = 1:length(regSEM);
timeAdj = timeAdj./30;
timeAdj = timeAdj - 5;
figure('Units','normalized','Position',[0 0 1 1])
shadedErrorBar(timeAdj,nanmean(llReg,2),regSEM,'lineprops',{'r-','LineWidth',3});
hold on
shadedErrorBar(timeAdj,nanmean(llRev,2),revSEM,'lineprops',{'b-','LineWidth',3});
xlabel('Time from Cue Onset (s)')
ylabel('Likelihood to Occupy Incorrect Sipper')
xline(0,'k--','Cue On','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(5,'k--','Sipper In','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(13,'k--','Sipper Out','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
title('Likelihood to Occupy Incorrect Sipper on Incorrect Approach Trials (All Sessions and Genotypes)')
legend([{'Congruent'},{'Incongruent'}])
xlim([-5 20])
ylim([0 1])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

%% Graph Alcohol Intake
regPIdx = startsWith(masterTbl.SessionType,'Regular') & strcmp(masterTbl.Strain,'P');
revPIdx = startsWith(masterTbl.SessionType,'Reversal') & strcmp(masterTbl.Strain,'P');

regWIdx = startsWith(masterTbl.SessionType,'Regular') & strcmp(masterTbl.Strain,'W');
revWIdx = startsWith(masterTbl.SessionType,'Reversal') & strcmp(masterTbl.Strain,'W');

intakePReg = masterTbl.Intake(regPIdx);
intakePRev = masterTbl.Intake(revPIdx);

intakeWReg = masterTbl.Intake(regWIdx);
intakeWRev = masterTbl.Intake(revWIdx);

g1 = repmat({'Regular P'},length(intakePReg),1);
g2 = repmat({'Reversal P'},length(intakePRev),1);
g3 = repmat({'Regular W'},length(intakeWReg),1);
g4 = repmat({'Reversal W'},length(intakeWRev),1);

g = [g1; g2; g3; g4];

semregcorLat = std(intakePReg)/sqrt(length(intakePReg));
semrevcorLat = std(intakePRev)/sqrt(length(intakePRev));
semregIncorLat = std(intakeWReg)/sqrt(length(intakeWReg));
semrevIncorLat = std(intakeWRev)/sqrt(length(intakeWRev));

figure('Units','normalized','Position',[0 0 1 1])
hb = bar(categorical({'P Rats','Wistars'}),[mean(intakePReg) mean(intakePRev); mean(intakeWReg) mean(intakeWRev)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);
hb(1).FaceColor = [0.5 0.5 0.5]; hb(2).FaceColor = [0.9 0.9 0.9];
offsetPos = [1+hb(1).XOffset 1+hb(2).XOffset 2+hb(1).XOffset 2+hb(2).XOffset];
scatterXData = [linspace(offsetPos(1)-0.1,offsetPos(1)+0.1,length(intakePReg)),linspace(offsetPos(2)-0.1,offsetPos(2)+0.1,length(intakePRev)),linspace(offsetPos(3)-0.1,offsetPos(3)+0.1,length(intakeWReg)),linspace(offsetPos(4)-0.1,offsetPos(4)+0.1,length(intakeWRev))];
hold on
plot([linspace(offsetPos(1)-0.05,offsetPos(1)+0.05,length(intakePReg)); linspace(offsetPos(2)-0.05,offsetPos(2)+0.05,length(intakePRev))],[intakePReg intakePRev]','k--o','LineWidth',2,'MarkerSize',8)
plot([linspace(offsetPos(3)-0.05,offsetPos(3)+0.05,length(intakeWReg)); linspace(offsetPos(4)-0.05,offsetPos(4)+0.05,length(intakeWRev))],[intakeWReg intakeWRev]','k--o','LineWidth',2,'MarkerSize',8)

er = errorbar(offsetPos,[mean(intakePReg) mean(intakePRev) mean(intakeWReg) mean(intakeWRev)],[semregcorLat semrevcorLat semregIncorLat semrevIncorLat],'LineWidth',3);    
er(1).Color = [0 0 0];                        
er(1).LineStyle = 'none'; 
legend([{'Congruent'},{'Incongruent'}])
ylabel('Intake (g/kg)')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)

%% Gross Location Data
% Standardize positional data based on box coordinates 
% boxCoords = [sipLDLC;sipRDLC;ULC;LLC;URC;LRC]; Left Sipper, Right Sipper,
% Upper Left Corner (ULC), Lower Left Corner (LLC), Upper Right Corner
% (URC), Lower Right Corner (LRC)
posRegCSPos = []; posRegCSMin = [];
posRevCSPos = []; posRevCSMin = [];
for i = 1:length(trlStruct)

    [nmBoxCoords,nmCenter,nmScale] = normalize(trlStruct(i).boxCoords);
    for k = 1:length(trialTimes)
        nmNosePos(k,:,:) = normalize(squeeze(trlStruct(i).trlnosePos(k,:,:)),'center',nmCenter,'scale',nmScale);
    end

    trlStruct(i).nmNosePos = nmNosePos;
    trlStruct(i).nmBoxCoords = nmBoxCoords;
end

for i = 1:length(dataSetIDs)
    if regexp(masterTbl.SessionType{i},regexptranslate('wildcard','Regular*'))
        posRegCSPos = [posRegCSPos; trlStruct(i).nmNosePos(1:48,:,:)];
        posRegCSMin = [posRegCSMin; trlStruct(i).nmNosePos(49:end,:,:)];
    elseif regexp(masterTbl.SessionType{i},regexptranslate('wildcard','Reversal*'))
        posRevCSPos = [posRevCSPos; trlStruct(i).nmNosePos(1:48,:,:)];
        posRevCSMin = [posRevCSMin; trlStruct(i).nmNosePos(49:end,:,:)];
    end
end

CSPosTrls = nmNosePos(1:48,:,:);
CSMinTrls = nmNosePos(49:end,:,:);

figure
h = histogram2(posRegCSPos(:,:,1),posRegCSPos(:,:,2),[20 20],'DisplayStyle','tile','Normalization','pdf');
hVal = h.Values;
figure
h2 = histogram2(posRegCSMin(:,:,1),posRegCSMin(:,:,2),[20 20],'DisplayStyle','tile','Normalization','pdf');
h2Val = h2.Values;

XLim = h.XBinEdges; YLim = h.YBinEdges;
diffValReg = normalize(hVal);

figure
h = histogram2(posRevCSPos(:,:,1),posRevCSPos(:,:,2),[20 20],'DisplayStyle','tile','Normalization','pdf');
hVal = h.Values;
figure
h2 = histogram2(posRevCSMin(:,:,1),posRevCSMin(:,:,2),[20 20],'DisplayStyle','tile','Normalization','pdf');
h2Val = h2.Values;
diffValRev = normalize(hVal);


figure
s = pcolor(XLim(2:end),YLim(2:end),normalize([diffValReg]));  s.EdgeColor = 'none';
colormap jet; grid off
hold on
scatter(nmBoxCoords(:,1),nmBoxCoords(:,2),50,'ko')
title('Congruent Session Position Estimates')
xlabel('Standardized X Coordinate')
ylabel('Standardized Y Coordinate')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)

figure
s = pcolor(XLim(2:end),YLim(2:end),normalize([diffValRev]));  s.EdgeColor = 'none';
colormap jet; grid off
hold on
scatter(nmBoxCoords(:,1),nmBoxCoords(:,2),50,'ko')
title('Incongruent Session Position Estimates')
xlabel('Standardized X Coordinate')
ylabel('Standardized Y Coordinate')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)


figure
s = pcolor(XLim(2:end),YLim(2:end),normalize([diffValRev - diffValReg])); s.EdgeColor = 'none';
colormap jet; grid off
hold on
scatter(nmBoxCoords(:,1),nmBoxCoords(:,2),50,'ko')
title('Incongruent - Congruent Session Position Estimates')
xlabel('Standardized X Coordinate')
ylabel('Standardized Y Coordinate')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)

%% Probability Testing on Location Data
% Find histcounts for every trial

%% Determine perserverative errors 
cueOn = 5*30;

for k = 1:length(trlStruct)

    test = [trlStruct(k).trlLSipDist <= 9; trlStruct(k).trlRSipDist <= 9];
    for i = 1:length(trialTimes)
        [~,firstSipSide] = min([find(test(i,cueOn:end),1,'first') find(test(i+length(trialTimes),cueOn:end),1,'first')]);
        if isempty(find(test(i,cueOn:end),1,'first')) && ~isempty(find(test(i+length(trialTimes),cueOn:end),1,'first'))
            firstSipTrl(i) = firstSipSide + 1;
        elseif isempty(find(test(i,cueOn:end),1,'first')) && isempty(find(test(i+length(trialTimes),cueOn:end),1,'first'))
            firstSipTrl(i) = 3;
        else
            firstSipTrl(i) = firstSipSide;
        end

    end

    trlStruct(k).firstSipSide = firstSipTrl;

end

% Key: 1 = Left Sipper First Visited, 2 = Right Sipper First Visited, 3 =
% Omission

figure
plot(trialTimes,firstSipTrl,'ko','MarkerSize',12)
hold on
plot(trialTimes,trialTypes,'r+','MarkerSize',12)
xlabel('Time (s)')
ylabel('Trial Type')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
ytick = [1 2 3];
yticklabels = {'Left','Right','CS-/Omission'};
ylim([0 4])

%% Sipper Occupancy Times
regCorOccu = [];   regIncOccu = [];
revCorOccu = [];   revIncOccu = [];                                                                         

for i = 1:length(trlStruct)

    if regexp(masterTbl.SessionType{i},regexptranslate('wildcard','Regular*'))
        regCorOccu = [regCorOccu; trlStruct(i).CorCorOcc];
        regIncOccu = [regIncOccu; trlStruct(i).CorIncOcc];
    elseif regexp(masterTbl.SessionType{i},regexptranslate('wildcard','Reversal*'))
        revCorOccu = [revCorOccu; trlStruct(i).CorCorOcc];
        revIncOccu = [revIncOccu; trlStruct(i).CorIncOcc];
    end

end

regCorOccu = regCorOccu(~isnan(regCorOccu))./30;
regIncOccu = regIncOccu(~isnan(regIncOccu))./30;
revCorOccu = revCorOccu(~isnan(revCorOccu))./30;
revIncOccu = revIncOccu(~isnan(revIncOccu))./30;

% Plot
semregcorOcc = std(regCorOccu)/sqrt(length(regCorOccu));
semrevcorOcc = std(regIncOccu)/sqrt(length(regIncOccu));
semregIncorOcc = std(revCorOccu)/sqrt(length(revCorOccu));
semrevIncorOcc = std(revIncOccu)/sqrt(length(revIncOccu));

figure('Units','normalized','Position',[0 0 1 1])
hb = bar(categorical({'Correct Approach Trials','Incorrect Approach Trials'}),[mean(regCorOccu) mean(revCorOccu); mean(regIncOccu) mean(revIncOccu)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);
hb(1).FaceColor = [0.5 0.5 0.5]; hb(2).FaceColor = [0.9 0.9 0.9];
offsetPos = [1+hb(1).XOffset 1+hb(2).XOffset 2+hb(1).XOffset 2+hb(2).XOffset];
scatterXData = [linspace(offsetPos(1)-0.1,offsetPos(1)+0.1,length(regCorOccu)),linspace(offsetPos(2)-0.1,offsetPos(2)+0.1,length(revCorOccu)),linspace(offsetPos(3)-0.1,offsetPos(3)+0.1,length(regIncOccu)),linspace(offsetPos(4)-0.1,offsetPos(4)+0.1,length(revIncOccu))];
hold on
scatter(scatterXData,[regCorOccu' revCorOccu' regIncOccu' revIncOccu'],50,'ko','LineWidth',2)
er = errorbar(offsetPos,[mean(regCorOccu) mean(revCorOccu) mean(regIncOccu) mean(revIncOccu)],[semregcorOcc semrevcorOcc semregIncorOcc semrevIncorOcc],'LineWidth',3);    
er(1).Color = [0 0 0];                        
er(1).LineStyle = 'none'; 
legend([{'Congruent'},{'Incongruent'}])
ylabel('Sipper Occupancy Time (s)')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)

%% Create Distance per Trial Metrics
for i = 1:length(trlStruct)
    for k = 1:length(trialTimes)
        trlStruct(i).trlDist(k) = sum(diag(squareform(pdist(squeeze(trlStruct(i).trlnosePos(k,cueOn:sipAscent,:)))),1)); % Sum of movement per trial. Taken from offdiagonal of pdist squareform per trial.
    end
end
%% Find distance moved per trial
distRegCSPos = []; distRegCSMin = [];
distRevCSPos = []; distRevCSMin = [];


for i = 1:length(trlStruct)

    if regexp(masterTbl.SessionType{i},regexptranslate('wildcard','Regular*'))
        currentTrialTimes = trlStruct(i).trialTimes;
        [sortedTimes,sortedIndex] = sort(currentTrialTimes(1:48));
        distRegCSPos = [distRegCSPos; trlStruct(i).trlDist(sortedIndex)];
        
    elseif regexp(masterTbl.SessionType{i},regexptranslate('wildcard','Reversal*'))
        currentTrialTimes = trlStruct(i).trialTimes;
        [sortedTimes,sortedIndex] = sort(currentTrialTimes(1:48));
        distRevCSPos = [distRevCSPos; trlStruct(i).trlDist(sortedIndex)];
        
    end

end

distRegSEM_CSPos = std(distRegCSPos)/sqrt(27);
distRevSEM_CSPos = std(distRevCSPos)/sqrt(27);

figure('Units','normalized','Position',[0 0 1 1])
plot(mean(distRegCSPos,1),'r-o','LineWidth',3); hold on; plot(mean(distRevCSPos,1),'b-o','LineWidth',3);
% er = errorbar([mean(distRegCSPos,1) mean(distRevCSPos,1)],[distRegSEM_CSPos distRevSEM_CSPos],'LineWidth',3);  
% er(1).Color = 'r'; er(2).Color = 'b';                          
% er(1).LineStyle = 'none'; er(2).LineStyle = 'none'; 
xlabel('Trial')
ylabel('Total Distance Traveled (pixels)')
title(['Total Distance Traveled per Trial'])
legend([{'Congruent'},{'Incongruent'}])

set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)

%% Plot neuron yield
neuronYield = [];
for i = 1:length(masterTbl.NeuronYield)

    neuronYield(i) = masterTbl.NeuronYield(i);
    
end

figure('Units','normalized','Position',[0 0 1 1])
bar(sort(neuronYield))
xlabel('Sorted Data Sets')
ylabel('Neuron Yield')
txt = 'Total Yield: 2830';
text(10,250,txt,'FontSize',30,'FontName','Arial','FontWeight','bold')
y1 = yline(mean(neuronYield),'--','Mean: 52.41','LineWidth',3,'FontSize',30,'FontName','Arial','FontWeight','bold');
y1.LabelHorizontalAlignment = 'left';
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)

%% Sipper Occupancy Per Trial
% Find distances below threshold for each CS+ trial, sort trials based on
% trial time & then pull and plot data. 
disThresh = 9;
disRun = 3;
sipDescent = 10*30;
sipAscent = 18*30;
cueOn = 5*30;
occReg = []; occRev = []; % Correct and Incorrect latencies

for i = 1:length(dataSetIDs)
    if regexp(masterTbl.SessionType{i},regexptranslate('wildcard','Regular*'))
        CorIdx = [trlStruct(i).LcorrectApproachTrls <= disThresh; trlStruct(i).RcorrectApproachTrls <= disThresh];
        currentTrialTimes = trlStruct(i).trialTimes;
        [sortedTimes,sortedIndex] = sort(currentTrialTimes(1:48));
        CorDbl = double(CorIdx);
        regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
        regOcc = regOcc(sortedIndex);
        trlStruct(i).trlOcc = regOcc;
    elseif regexp(masterTbl.SessionType{i},regexptranslate('wildcard','Reversal*'))
        CorIdx = [trlStruct(i).LcorrectApproachTrls <= disThresh; trlStruct(i).RcorrectApproachTrls <= disThresh];
        currentTrialTimes = trlStruct(i).trialTimes;
        [sortedTimes,sortedIndex] = sort(currentTrialTimes(1:48));
        CorDbl = double(CorIdx);
        regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
        regOcc = regOcc(sortedIndex);
        trlStruct(i).trlOcc = regOcc;
    end
end
% Pull data
for i = 1:length(dataSetIDs)
    if regexp(masterTbl.SessionType{i},regexptranslate('wildcard','Regular*'))
        occReg = [occReg trlStruct(i).trlOcc];
    elseif regexp(masterTbl.SessionType{i},regexptranslate('wildcard','Reversal*'))
        occRev = [occRev trlStruct(i).trlOcc];
    end
end
occReg = occReg ./ 30; % Convert time to seconds
occRev = occRev ./ 30; % Convert time to seconds
% Plot 
regSEM = std(occReg')/sqrt(min(size(occReg)));
revSEM = std(occRev')/sqrt(min(size(occRev)));
% Test plots
figure('Units','normalized','Position',[0 0 1 1])
shadedErrorBar(1:length(sortedIndex),nanmean(occReg,2),regSEM,'lineprops',{'r-o','LineWidth',3});
hold on
shadedErrorBar(1:length(sortedIndex),nanmean(occRev,2),revSEM,'lineprops',{'b-o','LineWidth',3});
xlabel('Trials')
ylabel('Sipper Occupancy (s)')

title('Time Spent at Sipper per Trial')
legend([{'Congruent'},{'Incongruent'}])
xlim([1 48])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

%% Calculate Inferred Rate of Intake Based on Total Sipper Occupancy per Session and Intake
% First find total sipper time
for i = 1:length(trlStruct)
    trlStruct(i).totalCorrectSipperTime = sum(trlStruct(i).CorCorOcc)./30;
end

% Divide total intake by total sipper time to find rate (g/kg/s)
for i = 1:length(trlStruct)
  trlStruct(i).intakeRate = masterTbl.Intake(i) / trlStruct(i).totalCorrectSipperTime;
end
    
regIntakeRate = []; revIntakeRate = [];
regIntake = [];     revIntake = [];

for i = 1:length(dataSetIDs)
    if regexp(masterTbl.SessionType{i},regexptranslate('wildcard','Regular*'))
        regIntakeRate = [regIntakeRate trlStruct(i).intakeRate];
        regIntake = [regIntake masterTbl.Intake(i)];
    elseif regexp(masterTbl.SessionType{i},regexptranslate('wildcard','Reversal*'))
        revIntakeRate = [revIntakeRate trlStruct(i).intakeRate];
        revIntake = [revIntake masterTbl.Intake(i)];

    end
end

semregIntakeRate  = std(regIntakeRate )/sqrt(length(regIntakeRate ));
semrevIntakeRate = std(revIntakeRate)/sqrt(length(revIntakeRate));

figure('Units','normalized','Position',[0 0 1 1])
bar(categorical({'Congruent Intake Rate','Incongruent Intake Rate'}),[mean(regIntakeRate ) mean(revIntakeRate)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3)
hold on
plot([linspace(0.95,1.05,length(regIntakeRate ));linspace(1.95,2.05,length(revIntakeRate))],[regIntakeRate ; revIntakeRate],'k--o','LineWidth',3,'MarkerSize',10)
er = errorbar(categorical({'Congruent Intake Rate','Incongruent Intake Rate'}),[mean(regIntakeRate ) mean(revIntakeRate)],[semregIntakeRate semrevIntakeRate],'LineWidth',3);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
ylabel('Intake Rate (g/kg/s)')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4) 

% Scatter plot, intake vs. rate of intake
figure('Units','normalized','Position',[0 0 1 1])
scatter(regIntakeRate,regIntake,100,'ro','LineWidth',3); hold on; scatter(revIntakeRate,revIntake,100,'bo','LineWidth',3);
xlabel("Intake Rate (g/kg/s)")
ylabel('Intake (g/kg)')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4) 

hold on

% Fit line
order = 1; % Linear fit

regCoefficients = polyfit(regIntakeRate,regIntake,order);
xReg = linspace(min(regIntake),max(regIntake),100);
yReg = polyval(regCoefficients,xReg);

revCoefficients = polyfit(revIntakeRate,revIntake,order);
xRev = linspace(min(revIntake),max(revIntake),100);
yRev = polyval(revCoefficients,xRev);

plot(xReg,yReg,'r-','LineWidth',3)
plot(xRev,yRev,'b-','LineWidth',3)

