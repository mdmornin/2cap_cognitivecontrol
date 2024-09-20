%% Notes
% DLC Left Sipper is OE Right Sipper
% DLC Right Sipper is OE Left Sipper


%% Analyze Trials DLC
leftLight = 0;
leftSipper = 2;
rightLight = 1;
rightSipper = 3;

addpath(genpath('H:\dissDat\analysis-tools-master'))
path = 'H:\dissDat\';
load([path '/trialStructureDLC.mat'])
load([path '/masterTable.mat'])

for i = 1:length(trialStructure)
    trials = size(trialStructure(i).snoutCSPlusPSTH,3); 
    timeBins = size(trialStructure(i).snoutCSPlusPSTH,1);
    for k = 1:trials
            for j = 1:timeBins
                if trialStructure(i).CorrectSipper(k) == 3  % Correct DLC Sipper is Right 
                    distanceFromCorrectSipper(j,k) = pdist([trialStructure(i).snoutCSPlusPSTH(j,:,k);trialStructure(i).rightSipperLocation]);
                elseif trialStructure(i).CorrectSipper(k) == 2 % Correct DLC Sipper is Left
                    distanceFromCorrectSipper(j,k) = pdist([trialStructure(i).snoutCSPlusPSTH(j,:,k);trialStructure(i).leftSipperLocation]); 
                end
            end
    end
    trialStructure(i).distanceFromCorrectSipper = distanceFromCorrectSipper;
end
allDistReg = [];
allDistRev = [];

for i = 1:length(trialStructure)
    for k = 1:27
    if strcmp(masterTbl.SessionType(i),'Regular')
        allDistReg(:,:,k) = trialStructure(i).distanceFromCorrectSipper;
    elseif strcmp(masterTbl.SessionType(i),'Reversal')
        allDistRev(:,:,k) = trialStructure(i).distanceFromCorrectSipper;
    end
    end
end

timeStamps = (1:size(allDistReg,1))/30;
timeStamps = timeStamps - 5;

% Errors
meanErrorReg = std(mean(allDistReg,3)')/sqrt(size(allDistReg,2));
meanErrorRev = std(mean(allDistRev,3)')/sqrt(size(allDistRev,2));

meanErrorReg124 = std(mean(allDistReg(:,1:24,:),3)')/sqrt(size(allDistReg,2));
meanErrorRev124 = std(mean(allDistRev(:,1:24,:),3)')/sqrt(size(allDistRev,2));

meanErrorReg2548 = std(mean(allDistReg(:,25:48,:),3)')/sqrt(size(allDistReg,2));
meanErrorRev2548 = std(mean(allDistRev(:,25:48,:),3)')/sqrt(size(allDistRev,2));

figure

s1 = shadedErrorBar(timeStamps,mean(mean(allDistReg,3),2),meanErrorReg,'lineprops','b-');
hold on
s2 = shadedErrorBar(timeStamps,mean(mean(allDistRev,3),2),meanErrorRev,'lineprops','r-');


xline(5,'-','CS On')
xline(10,'-','Sipper In')
xline(18, '-', 'Sipper Out')

xlabel('Time (s)')
ylabel('Distance from Correct Sipper (pixels)')
title('Distance from Sipper, All Trials')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold')
maxDist = pdist([trialStructure(1).rightSipperLocation;trialStructure(1).leftSipperLocation]);
yline(maxDist,'-','Distance Between Sippers')
ylim([120 350])

figure
shadedErrorBar(timeStamps,mean(mean(allDistReg(:,1:24,:),3),2),meanErrorReg124,'lineprops','b-');
hold on
shadedErrorBar(timeStamps,mean(mean(allDistRev(:,1:24,:),3),2),meanErrorRev124,'lineprops','r-');
xline(5,'-','CS On')
xline(10,'-','Sipper In')
xline(18, '-', 'Sipper Out')

xlabel('Time (s)')
ylabel('Distance from Correct Sipper (pixels)')
title('Distance from Sipper, First 24')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold')
ylim([120 350])
yline(maxDist,'-','Distance Between Sippers')


figure
shadedErrorBar(timeStamps,mean(mean(allDistReg(:,25:48,:),3),2),meanErrorReg124,'lineprops','b-');
hold on
shadedErrorBar(timeStamps,mean(mean(allDistRev(:,25:48,:),3),2),meanErrorRev124,'lineprops','r-');
xline(5,'-','CS On')
xline(10,'-','Sipper In')
xline(18, '-', 'Sipper Out')

xlabel('Time (s)')
ylabel('Distance from Correct Sipper (pixels)')
title('Distance from Sipper, Last 24')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold')
ylim([120 350])
yline(maxDist,'-','Distance Between Sippers')

