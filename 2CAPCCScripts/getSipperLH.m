function [sipperOccupancyStatistics,trlStruct] = getSipperLH(trlStruct,figSavePath,masterTbl,flags)
% Parameters for every 
disThresh = 9;
disRun = 3;
sipDescent = 10*30;
sipAscent = 18*30;
cueOn = 5*30;
LeftTrials = 1:24;
RightTrials = 25:48;
% Variable init
llReg = []; llRev = []; % Correct and Incorrect latencies
%% Correct sipper on correct sipper trials
% Pulling data
if strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'all') 

    for i = 1:length(trlStruct)
        if startsWith(masterTbl.SessionType{i},['Regular'])
            CorIdx = [trlStruct(i).RDistCorrectAppr <= disThresh; trlStruct(i).LDistCorrectAppr <= disThresh];
            CorDbl = double(CorIdx);
            regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
            llOccupy = movmean(CorDbl',3);
            trlStruct(i).llOccupy_CC = llOccupy;
            trlStruct(i).CorCorOcc = regOcc;
            trlStruct(i).llOccupy = llOccupy;
        elseif startsWith(masterTbl.SessionType{i},['Reversal'])
            CorIdx = [trlStruct(i).RDistCorrectAppr <= disThresh; trlStruct(i).LDistCorrectAppr <= disThresh];
            CorDbl = double(CorIdx);
            regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
            llOccupy = movmean(CorDbl',3);
            trlStruct(i).llOccupy_CC = llOccupy;
            trlStruct(i).CorCorOcc = regOcc;
            trlStruct(i).llOccupy = llOccupy;
        end
    end

elseif strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'P')

    for i = 1:length(trlStruct)
        if startsWith(masterTbl.SessionType{i},['Regular'])  && strcmp(masterTbl.Strain{i},'P')
            CorIdx = [trlStruct(i).RDistCorrectAppr <= disThresh; trlStruct(i).LDistCorrectAppr <= disThresh];
            CorDbl = double(CorIdx);
            regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
            llOccupy = movmean(CorDbl',3);
            trlStruct(i).llOccupy = llOccupy;
            trlStruct(i).CorCorOcc = regOcc;
        elseif startsWith(masterTbl.SessionType{i},['Reversal'])  && strcmp(masterTbl.Strain{i},'P')
            CorIdx = [trlStruct(i).RDistCorrectAppr <= disThresh; trlStruct(i).LDistCorrectAppr <= disThresh];
            CorDbl = double(CorIdx);
            regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
            llOccupy = movmean(CorDbl',3);
            trlStruct(i).llOccupy = llOccupy;
            trlStruct(i).CorCorOcc = regOcc;
        end
    end

elseif strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'W')

    for i = 1:length(trlStruct)
        if startsWith(masterTbl.SessionType{i},['Regular'])  && strcmp(masterTbl.Strain{i},'W')
            CorIdx = [trlStruct(i).RDistCorrectAppr <= disThresh; trlStruct(i).LDistCorrectAppr <= disThresh];
            CorDbl = double(CorIdx);
            regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
            llOccupy = movmean(CorDbl',3);
            trlStruct(i).llOccupy = llOccupy;
            trlStruct(i).CorCorOcc = regOcc;
        elseif startsWith(masterTbl.SessionType{i},['Reversal'])  && strcmp(masterTbl.Strain{i},'W')
            CorIdx = [trlStruct(i).RDistCorrectAppr <= disThresh; trlStruct(i).LDistCorrectAppr <= disThresh];
            CorDbl = double(CorIdx);
            regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
            llOccupy = movmean(CorDbl',3);
            trlStruct(i).llOccupy = llOccupy;
            trlStruct(i).CorCorOcc = regOcc;
        end
    end

end
% Further pulling data into vector for plotting
for i = 1:length(trlStruct)
    if regexp(masterTbl.SessionType{i},regexptranslate('wildcard','Regular*'))
        llReg = [llReg trlStruct(i).llOccupy];
    elseif regexp(masterTbl.SessionType{i},regexptranslate('wildcard','Reversal*'))
        llRev = [llRev trlStruct(i).llOccupy];
    end
end
% Plotting data
% Find SEM 
regSEM = std(llReg')/sqrt(min(size(llReg)));
revSEM = std(llRev')/sqrt(min(size(llRev)));
% Adjust time to be Seconds and Center on Cue Presentation
timeAdj = 1:length(llRev);
timeAdj = timeAdj./30;
timeAdj = timeAdj - 5;
% Generate plot
figure('Units','normalized','Position',[0 0 1 1])
shadedErrorBar(timeAdj,nanmean(llReg,2),regSEM,'lineprops',{'r-','LineWidth',3});
hold on
shadedErrorBar(timeAdj,nanmean(llRev,2),revSEM,'lineprops',{'b-','LineWidth',3});
xlabel('Time from Cue Onset (s)')
ylabel('Likelihood to Occupy Correct Sipper')
xline(0,'k--','Cue On','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(5,'k--','Sipper In','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(13,'k--','Sipper Out','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
title('Likelihood to Occupy Correct Sipper on Correct Approach Trials')
legend([{'Congruent'},{'Incongruent'}])
xlim([-5 20])
ylim([0 1])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,[figSavePath filesep 'sipperOccupyLH_CorrectSipper_CorrectApproach_' flags.SessionN '_Strain_' flags.Genotype],'png')

sipperOccupancyStatistics.CorCorReg = llReg;
sipperOccupancyStatistics.CorCorRev = llRev;

%% Correct Sipper, Incorrect Approach Trials
% Reset variables
llReg = []; llRev = []; % Correct and Incorrect latencies
trlStruct = rmfield(trlStruct,"llOccupy");
if strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'all') 

    for i = 1:length(trlStruct)
        if startsWith(masterTbl.SessionType{i},['Regular'])
            CorIdx = [trlStruct(i).RDistInCorrectAppr <= disThresh; trlStruct(i).LDistInCorrectAppr <= disThresh];
            CorDbl = double(CorIdx);
            regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
            llOccupy = movmean(CorDbl',3);
            trlStruct(i).CorIncOcc = regOcc;
            trlStruct(i).llOccupy_CI = llOccupy;
            trlStruct(i).llOccupy = llOccupy;
        elseif startsWith(masterTbl.SessionType{i},['Reversal'])
            CorIdx = [trlStruct(i).RDistInCorrectAppr <= disThresh; trlStruct(i).LDistInCorrectAppr <= disThresh];
            CorDbl = double(CorIdx);
            regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
            llOccupy = movmean(CorDbl',3);
            trlStruct(i).CorIncOcc = regOcc;
            trlStruct(i).llOccupy_CI = llOccupy;
            trlStruct(i).llOccupy = llOccupy;

        end
    end

elseif strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'P')

    for i = 1:length(trlStruct)
        if startsWith(masterTbl.SessionType{i},['Regular'])  && strcmp(masterTbl.Strain{i},'P')
            CorIdx = [trlStruct(i).RDistInCorrectAppr <= disThresh; trlStruct(i).LDistInCorrectAppr <= disThresh];
            CorDbl = double(CorIdx);
            regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
            llOccupy = movmean(CorDbl',3);
            trlStruct(i).CorIncOcc = regOcc;
            trlStruct(i).llOccupy = llOccupy;
        elseif startsWith(masterTbl.SessionType{i},['Reversal'])  && strcmp(masterTbl.Strain{i},'P')
            CorIdx = [trlStruct(i).RDistInCorrectAppr <= disThresh; trlStruct(i).LDistInCorrectAppr <= disThresh];
            CorDbl = double(CorIdx);
            regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
            llOccupy = movmean(CorDbl',3);
            trlStruct(i).CorIncOcc = regOcc;
            trlStruct(i).llOccupy = llOccupy;
        end
    end

elseif strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'W')

    for i = 1:length(trlStruct)
        if startsWith(masterTbl.SessionType{i},['Regular'])  && strcmp(masterTbl.Strain{i},'W')
            CorIdx = [trlStruct(i).RDistInCorrectAppr <= disThresh; trlStruct(i).LDistInCorrectAppr <= disThresh];
            CorDbl = double(CorIdx);
            regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
            llOccupy = movmean(CorDbl',3);
            trlStruct(i).CorIncOcc = regOcc;
            trlStruct(i).llOccupy = llOccupy;
        elseif startsWith(masterTbl.SessionType{i},['Reversal'])  && strcmp(masterTbl.Strain{i},'W')
            CorIdx = [trlStruct(i).RDistInCorrectAppr <= disThresh; trlStruct(i).LDistInCorrectAppr <= disThresh];
            CorDbl = double(CorIdx);
            regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
            llOccupy = movmean(CorDbl',3);
            trlStruct(i).CorIncOcc = regOcc;
            trlStruct(i).llOccupy = llOccupy;
        end
    end

end
% Further pulling data into vector for plotting

for i = 1:length(trlStruct)
    if startsWith(masterTbl.SessionType{i},['Regular'])
        llReg = [llReg trlStruct(i).llOccupy];
    elseif startsWith(masterTbl.SessionType{i},['Reversal'])
        llRev = [llRev trlStruct(i).llOccupy];
    end
end

% Plotting data
% Find SEM 
regSEM = std(llReg')/sqrt(min(size(llReg)));
revSEM = std(llRev')/sqrt(min(size(llRev)));
% Adjust time to be Seconds and Center on Cue Presentation
timeAdj = 1:length(llRev);
timeAdj = timeAdj./30;
timeAdj = timeAdj - 5;
% Generate plot
figure('Units','normalized','Position',[0 0 1 1])
shadedErrorBar(timeAdj,nanmean(llReg,2),regSEM,'lineprops',{'r-','LineWidth',3});
hold on
shadedErrorBar(timeAdj,nanmean(llRev,2),revSEM,'lineprops',{'b-','LineWidth',3});
xlabel('Time from Cue Onset (s)')
ylabel('Likelihood to Occupy Correct Sipper')
xline(0,'k--','Cue On','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(5,'k--','Sipper In','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(13,'k--','Sipper Out','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
title('Likelihood to Occupy Correct Sipper on Incorrect Approach Trials')
legend([{'Congruent'},{'Incongruent'}])
xlim([-5 20])
ylim([0 1])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,[figSavePath filesep 'sipperOccupyLH_CorrectSipper_IncorrectApproach_' flags.SessionN '_Strain_' flags.Genotype],'png')

sipperOccupancyStatistics.CorIncorrectCong = llReg;
sipperOccupancyStatistics.CorIncorrectIncong = llRev;
%% Incorrect Sipper, Correct Approach Trials
llReg = []; llRev = []; % Correct and Incorrect latencies
trlStruct = rmfield(trlStruct,"llOccupy");
if strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'all') 

    for i = 1:length(trlStruct)
        if startsWith(masterTbl.SessionType{i},['Regular'])
            CorIdx = [trlStruct(i).trlRSipDist(LeftTrials(trlStruct(i).LcorrectIdx),:) <= disThresh; trlStruct(i).trlLSipDist(RightTrials(trlStruct(i).RcorrectIdx),:) <= disThresh];
            CorDbl = double(CorIdx);
            llOccupy = movmean(CorDbl',3);
            regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
            trlStruct(i).IncCorOcc = regOcc;
            trlStruct(i).llOccupy_IC = llOccupy;
            trlStruct(i).llOccupy = llOccupy;
        elseif startsWith(masterTbl.SessionType{i},['Reversal'])
            CorIdx = [trlStruct(i).trlRSipDist(RightTrials(trlStruct(i).RcorrectIdx),:) <= disThresh; trlStruct(i).trlLSipDist(LeftTrials(trlStruct(i).LcorrectIdx),:) <= disThresh];
            CorDbl = double(CorIdx);
            llOccupy = movmean(CorDbl',3);
            regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
            trlStruct(i).IncCorOcc = regOcc;
            trlStruct(i).llOccupy_IC = llOccupy;
            trlStruct(i).llOccupy = llOccupy;
        end
    end

elseif strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'P')

    for i = 1:length(trlStruct)
        if startsWith(masterTbl.SessionType{i},['Regular'])  && strcmp(masterTbl.Strain{i},'P')
            CorIdx = [trlStruct(i).trlRSipDist(LeftTrials(trlStruct(i).LcorrectIdx),:) <= disThresh; trlStruct(i).trlLSipDist(RightTrials(trlStruct(i).RcorrectIdx),:) <= disThresh];
            CorDbl = double(CorIdx);
            llOccupy = movmean(CorDbl',3);
            regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
            trlStruct(i).IncCorOcc = regOcc;
            trlStruct(i).llOccupy = llOccupy;
        elseif startsWith(masterTbl.SessionType{i},['Reversal'])  && strcmp(masterTbl.Strain{i},'P')
            CorIdx = [trlStruct(i).trlRSipDist(RightTrials(trlStruct(i).RcorrectIdx),:) <= disThresh; trlStruct(i).trlLSipDist(LeftTrials(trlStruct(i).LcorrectIdx),:) <= disThresh];
            CorDbl = double(CorIdx);
            llOccupy = movmean(CorDbl',3);
            regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
            trlStruct(i).IncCorOcc = regOcc;
            trlStruct(i).llOccupy = llOccupy;
        end
    end

elseif strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'W')

    for i = 1:length(trlStruct)
        if startsWith(masterTbl.SessionType{i},['Regular'])  && strcmp(masterTbl.Strain{i},'W')
            CorIdx = [trlStruct(i).trlRSipDist(LeftTrials(trlStruct(i).LcorrectIdx),:) <= disThresh; trlStruct(i).trlLSipDist(RightTrials(trlStruct(i).RcorrectIdx),:) <= disThresh];
            CorDbl = double(CorIdx);
            llOccupy = movmean(CorDbl',3);
            regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
            trlStruct(i).IncCorOcc = regOcc;
            trlStruct(i).llOccupy = llOccupy;
        elseif startsWith(masterTbl.SessionType{i},['Reversal'])  && strcmp(masterTbl.Strain{i},'W')
            CorIdx = [trlStruct(i).trlRSipDist(RightTrials(trlStruct(i).RcorrectIdx),:) <= disThresh; trlStruct(i).trlLSipDist(LeftTrials(trlStruct(i).LcorrectIdx),:) <= disThresh];
            CorDbl = double(CorIdx);
            llOccupy = movmean(CorDbl',3);
            regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
            trlStruct(i).IncCorOcc = regOcc;
            trlStruct(i).llOccupy = llOccupy;
        end
    end

end
% Further pulling data into vector for plotting

for i = 1:length(trlStruct)
    if startsWith(masterTbl.SessionType{i},['Regular'])
        llReg = [llReg trlStruct(i).llOccupy];
    elseif startsWith(masterTbl.SessionType{i},['Reversal'])
        llRev = [llRev trlStruct(i).llOccupy];
    end
end

% Plotting data
% Find SEM 
regSEM = std(llReg')/sqrt(min(size(llReg)));
revSEM = std(llRev')/sqrt(min(size(llRev)));
% Adjust time to be Seconds and Center on Cue Presentation
timeAdj = 1:length(llRev);
timeAdj = timeAdj./30;
timeAdj = timeAdj - 5;
% Generate plot
figure('Units','normalized','Position',[0 0 1 1])
shadedErrorBar(timeAdj,nanmean(llReg,2),regSEM,'lineprops',{'r-','LineWidth',3});
hold on
shadedErrorBar(timeAdj,nanmean(llRev,2),revSEM,'lineprops',{'b-','LineWidth',3});
xlabel('Time from Cue Onset (s)')
ylabel('Likelihood to Occupy Correct Sipper')
xline(0,'k--','Cue On','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(5,'k--','Sipper In','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(13,'k--','Sipper Out','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
title('Likelihood to Occupy Incorrect Sipper on Correct Approach Trials')
legend([{'Congruent'},{'Incongruent'}])
xlim([-5 20])
ylim([0 1])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,[figSavePath filesep 'sipperOccupyLH_IncorrectSipper_CorrectApproach_' flags.SessionN '_Strain_' flags.Genotype],'png')

sipperOccupancyStatistics.IncCorrectCong = llReg;
sipperOccupancyStatistics.IncCorrectIncong = llRev;

%% Incorrect Sipper on Incorrect Approach Trials
% Reset variables
llReg = []; llRev = []; % Correct and Incorrect latencies
trlStruct = rmfield(trlStruct,"llOccupy");

if strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'all') 

    for i = 1:length(trlStruct)
        if startsWith(masterTbl.SessionType{i},['Regular'])
            CorIdx = [trlStruct(i).trlRSipDist(LeftTrials(trlStruct(i).LincorrectIdx),:) <= disThresh; trlStruct(i).trlLSipDist(RightTrials(trlStruct(i).RincorrectIdx),:) <= disThresh];
            CorDbl = double(CorIdx);
            llOccupy = movmean(CorDbl',3);
            regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
            trlStruct(i).IncIncOcc = regOcc;
            trlStruct(i).llOccupy_II = llOccupy;
            trlStruct(i).llOccupy = llOccupy;
        elseif startsWith(masterTbl.SessionType{i},['Reversal'])
            CorIdx = [trlStruct(i).trlRSipDist(RightTrials(trlStruct(i).RincorrectIdx),:) <= disThresh; trlStruct(i).trlLSipDist(LeftTrials(trlStruct(i).LincorrectIdx),:) <= disThresh];
            CorDbl = double(CorIdx);
            llOccupy = movmean(CorDbl',3);
            regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
            trlStruct(i).IncIncOcc = regOcc;
            trlStruct(i).llOccupy_II = llOccupy;
            trlStruct(i).llOccupy = llOccupy;
        end
    end

elseif strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'P')

    for i = 1:length(trlStruct)
        if startsWith(masterTbl.SessionType{i},['Regular'])  && strcmp(masterTbl.Strain{i},'P')
            CorIdx = [trlStruct(i).trlRSipDist(LeftTrials(trlStruct(i).LincorrectIdx),:) <= disThresh; trlStruct(i).trlLSipDist(RightTrials(trlStruct(i).RincorrectIdx),:) <= disThresh];
            CorDbl = double(CorIdx);
            llOccupy = movmean(CorDbl',3);
            regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
            trlStruct(i).IncIncOcc = regOcc;
            trlStruct(i).llOccupy = llOccupy;
        elseif startsWith(masterTbl.SessionType{i},['Reversal'])  && strcmp(masterTbl.Strain{i},'P')
            CorIdx = [trlStruct(i).trlRSipDist(RightTrials(trlStruct(i).RincorrectIdx),:) <= disThresh; trlStruct(i).trlLSipDist(LeftTrials(trlStruct(i).LincorrectIdx),:) <= disThresh];
            CorDbl = double(CorIdx);
            llOccupy = movmean(CorDbl',3);
            regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
            trlStruct(i).IncIncOcc = regOcc;
            trlStruct(i).llOccupy = llOccupy;
        end
    end

elseif strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'W')

    for i = 1:length(trlStruct)
        if startsWith(masterTbl.SessionType{i},['Regular'])  && strcmp(masterTbl.Strain{i},'W')
            CorIdx = [trlStruct(i).trlRSipDist(LeftTrials(trlStruct(i).LincorrectIdx),:) <= disThresh; trlStruct(i).trlLSipDist(RightTrials(trlStruct(i).RincorrectIdx),:) <= disThresh];
            CorDbl = double(CorIdx);
            llOccupy = movmean(CorDbl',3);
            regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
            trlStruct(i).IncIncOcc = regOcc;
            trlStruct(i).llOccupy = llOccupy;
        elseif startsWith(masterTbl.SessionType{i},['Reversal'])  && strcmp(masterTbl.Strain{i},'W')
            CorIdx = [trlStruct(i).trlRSipDist(RightTrials(trlStruct(i).RincorrectIdx),:) <= disThresh; trlStruct(i).trlLSipDist(LeftTrials(trlStruct(i).LincorrectIdx),:) <= disThresh];
            CorDbl = double(CorIdx);
            llOccupy = movmean(CorDbl',3);
            regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
            trlStruct(i).IncIncOcc = regOcc;
            trlStruct(i).llOccupy = llOccupy;
        end
    end

end
% Further pulling data into vector for plotting

for i = 1:length(trlStruct)
    if startsWith(masterTbl.SessionType{i},['Regular'])
        llReg = [llReg trlStruct(i).llOccupy];
    elseif startsWith(masterTbl.SessionType{i},['Reversal'])
        llRev = [llRev trlStruct(i).llOccupy];
    end
end

% Plotting data
% Find SEM 
regSEM = std(llReg')/sqrt(min(size(llReg)));
revSEM = std(llRev')/sqrt(min(size(llRev)));
% Adjust time to be Seconds and Center on Cue Presentation
timeAdj = 1:length(llRev);
timeAdj = timeAdj./30;
timeAdj = timeAdj - 5;
% Generate plot
figure('Units','normalized','Position',[0 0 1 1])
shadedErrorBar(timeAdj,nanmean(llReg,2),regSEM,'lineprops',{'r-','LineWidth',3});
hold on
shadedErrorBar(timeAdj,nanmean(llRev,2),revSEM,'lineprops',{'b-','LineWidth',3});
xlabel('Time from Cue Onset (s)')
ylabel('Likelihood to Occupy Correct Sipper')
xline(0,'k--','Cue On','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(5,'k--','Sipper In','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(13,'k--','Sipper Out','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
title('Likelihood to Occupy Incorrect Sipper on Incorrect Approach Trials')
legend([{'Congruent'},{'Incongruent'}])
xlim([-5 20])
ylim([0 1])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,[figSavePath filesep 'sipperOccupyLH_IncorrectSipper_IncorrectApproach_' flags.SessionN '_Strain_' flags.Genotype],'png')

sipperOccupancyStatistics.IncIncorrectCong = llReg;
sipperOccupancyStatistics.IncIncorrectIncong = llRev;

%% Output trlStruct because it contains new fields relevant for analysis

end