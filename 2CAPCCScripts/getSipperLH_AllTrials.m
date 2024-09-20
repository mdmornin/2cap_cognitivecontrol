function [sipperOccupancyStatistics,trlStruct] = getSipperLH_AllTrials(trlStruct,figSavePath,masterTbl,flags)
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
            [~,sortIdx] = sort(trlStruct(i).trialTimes(1:48));
            CorIdx = [trlStruct(i).trlLSipDist(LeftTrials,:) <= disThresh; trlStruct(i).trlRSipDist(RightTrials,:) <= disThresh];
            CorIdx = CorIdx(sortIdx,:);
            CorDbl = double(CorIdx);
            trlStruct(i).noCorrectionsIdx = trlStruct(i).approach(1:48,1,1) == 1 & trlStruct(i).approach(1:48,2,1) == 0;
            regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
            llOccupy = movmean(CorDbl',3);
            trlStruct(i).llOccupy_CC = llOccupy;
            trlStruct(i).CorCorOcc = regOcc;
            trlStruct(i).llOccupy = llOccupy;
        elseif startsWith(masterTbl.SessionType{i},['Reversal'])
            [~,sortIdx] = sort(trlStruct(i).trialTimes(1:48));
            CorIdx = [trlStruct(i).trlLSipDist(RightTrials,:) <= disThresh; trlStruct(i).trlRSipDist(LeftTrials,:) <= disThresh];
            CorIdx = CorIdx(sortIdx,:);
            CorDbl = double(CorIdx);
            trlStruct(i).noCorrectionsIdx = trlStruct(i).approach(1:48,1,1) == 0 & trlStruct(i).approach(1:48,2,1) == 1;
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
            [~,sortIdx] = sort(trlStruct(i).trialTimes(1:48));
            CorIdx = [trlStruct(i).trlLSipDist(LeftTrials,:) <= disThresh; trlStruct(i).trlRSipDist(RightTrials,:) <= disThresh];
            CorIdx = CorIdx(sortIdx,:);
            CorDbl = double(CorIdx);
            regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
            llOccupy = movmean(CorDbl',3);
            trlStruct(i).llOccupy_CC = llOccupy;
            trlStruct(i).CorCorOcc = regOcc;
            trlStruct(i).llOccupy = llOccupy;
        elseif startsWith(masterTbl.SessionType{i},['Reversal'])  && strcmp(masterTbl.Strain{i},'P')
            [~,sortIdx] = sort(trlStruct(i).trialTimes(1:48));
            CorIdx = [trlStruct(i).trlLSipDist(RightTrials,:) <= disThresh; trlStruct(i).trlRSipDist(LeftTrials,:) <= disThresh];
            CorIdx = CorIdx(sortIdx,:);
            CorDbl = double(CorIdx);
            regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
            llOccupy = movmean(CorDbl',3);
            trlStruct(i).llOccupy_CC = llOccupy;
            trlStruct(i).CorCorOcc = regOcc;
            trlStruct(i).llOccupy = llOccupy;
        end
    end

elseif strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'W')

    for i = 1:length(trlStruct)
        if startsWith(masterTbl.SessionType{i},['Regular'])  && strcmp(masterTbl.Strain{i},'W')
            [~,sortIdx] = sort(trlStruct(i).trialTimes(1:48));
            CorIdx = [trlStruct(i).trlLSipDist(LeftTrials,:) <= disThresh; trlStruct(i).trlRSipDist(RightTrials,:) <= disThresh];
            CorIdx = CorIdx(sortIdx,:);
            CorDbl = double(CorIdx);
            regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
            llOccupy = movmean(CorDbl',3);
            trlStruct(i).llOccupy_CC = llOccupy;
            trlStruct(i).CorCorOcc = regOcc;
            trlStruct(i).llOccupy = llOccupy;
        elseif startsWith(masterTbl.SessionType{i},['Reversal'])  && strcmp(masterTbl.Strain{i},'W')
            [~,sortIdx] = sort(trlStruct(i).trialTimes(1:48));
            CorIdx = [trlStruct(i).trlLSipDist(RightTrials,:) <= disThresh; trlStruct(i).trlRSipDist(LeftTrials,:) <= disThresh];
            CorIdx = CorIdx(sortIdx,:);
            CorDbl = double(CorIdx);
            regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
            llOccupy = movmean(CorDbl',3);
            trlStruct(i).llOccupy_CC = llOccupy;
            trlStruct(i).CorCorOcc = regOcc;
            trlStruct(i).llOccupy = llOccupy;
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
timeAdj = 1:min(size(llRev));
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
title('Likelihood to Occupy Correct Sipper, All Trials')
legend([{'Congruent'},{'Incongruent'}])
xlim([-5 20])
ylim([0 1])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,[figSavePath filesep 'allTrials_sipperOccupyLH_CorrectSipper_CorrectApproach_' flags.SessionN '_Strain_' flags.Genotype],'png')

sipperOccupancyStatistics.CorCorReg = llReg;
sipperOccupancyStatistics.CorCorRev = llRev;

end