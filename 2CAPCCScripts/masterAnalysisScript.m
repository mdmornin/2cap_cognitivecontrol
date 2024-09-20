%% Master script for running 2CAP Congruent/Incongruent Session Analyses
% Run functions based on prior script, output enough data for future
% analyses and statistics. Output & save figures to defined figure path. 
% Figure save names are dynamic. 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% - - - - - - - - - - - - - - - Change Log - - - - - - - - - - - - - - - 
% 
% 02/02/2022 - 02/03/2022 
% Reimplemented previous analyses from 'analyzeApproach' script into
% functions. All functions currently have the ability to separate analyses
% by genotype BUT not by session #. I am hoping to be able to find more
% efficient ways of implementing flagging systems into the functions. Right
% now it is a series of if, elseif statements that repeat based on the
% flags. Sub-functions within each function may aid in this, but it is a
% question of whether it is worthwhile.
%
% Other future goals include supression of figures and dynamic inputs to
% configure figure extensions, i.e. .fig, .png, .svg.
% 
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% - - - - - - - - - - - - - - - General notes - - - - - - - - - - - - - - 
% MDM - IUPUI 02/03/2022
% Analyses for 2CAP Behaviors
% Behavioral only as of now.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% - - - - - - - - - - - - - - - Variable names - - - - - - - - - - - - - 

%% Add paths, path variables
addpath(genpath('F:\dissDat\analysis-tools-master'))
addpath(genpath('F:\dissDat\res'))

dataPath = 'F:/dissDat/res'; % Path to outputs from pWistDataOrganizerVer3. These are .mat files for each dataset.
checkPath = 'F:/dissDat/trlStruct.mat'; % Path to check existence of trlStruct variable. Preexisting variables will be loaded.
savePath = 'F:/dissDat/'; % Path to save output trlStruct
tablePath = 'F:/dissDat/masterTable.mat'; % Path to masterTable.mat which contains data set information and names
figSavePath = 'F:/dissDat/figs';
addpath(genpath('F:/dissDat/restoredScripts'))
addpath(genpath('F:/dissDat/analysisScripts'))
load(tablePath);
%% Generate or load trlStruct
if ~exist('trlStruct','var')
    trlStruct = trlStructLoadorMake(dataPath,checkPath,savePath,tablePath);
else
    disp('trlStruct.mat already exists in the workspace')
end
%% Analysis 1. Session statistics - Correct approaches, incorrect approaches, omissions, corrections
%
flags.Genotype = 'all'; % Genotype must be 'all' 'P' or 'W'
flags.SessionN = 'all'; % SessionN must be 'all' '1' '2' or '3'

[trlStruct,sessionStats] = getSessionStats(trlStruct,figSavePath,masterTbl,flags);

%% Analysis 2. Latencies to Approach Correct or Incorrect Sipper
flags.Genotype = 'all'; % Genotype must be 'all' 'P' or 'W'
flags.SessionN = 'all'; % SessionN must be 'all' '1' '2' or '3'

[trlStruct,latencyStats] = getLatency(trlStruct,figSavePath,masterTbl,flags);
%% Analysis 3. Likelihood to Approach, Correct or Incorrect and combined
flags.Genotype = 'all'; % Genotype must be 'all' 'P' or 'W'
flags.SessionN = 'all'; % SessionN must be 'all' '1' '2' or '3'
flags.combined = 0;     % 1 means combine correct and incorrect approaches, 0 means separate 

[approachStatistics,trlStruct] = getApproachLH(trlStruct,masterTbl,figSavePath,flags);

%% Analysis 4. Likelihood to Occupy Correct, Incorrect Sipper on Correct or Incorrect Approach Trials
flags.Genotype = 'all'; % Genotype must be 'all' 'P' or 'W'
flags.SessionN = 'all'; % SessionN must be 'all' '1' '2' or '3'
[sipperOccupancyStatistics,trlStruct] = getSipperLH(trlStruct,figSavePath,masterTbl,flags);

%% Analysis 5. Sipper Occupancy Time per Session Type
% Utilizes trlStruct output from Analysis 4. 
flags.Genotype = 'all'; % Genotype must be 'all' 'P' or 'W'
flags.SessionN = 'all'; % SessionN must be 'all' '1' '2' or '3'
[sipperTimeStats,trlStruct] = getSipperTimePerSession(trlStruct,figSavePath,masterTbl,flags);

%% Analysis 6. Sipper Occupancy Time per Trial 
flags.Genotype = 'all'; % Genotype must be 'all' 'P' or 'W'
flags.SessionN = 'all'; % SessionN must be 'all' '1' '2' or '3'
[trlStruct,trialSipperTime] = getTrialSipperTime(trlStruct,figSavePath,masterTbl,flags);

%% Analysis 7. Correction LH
flags.Genotype = 'all'; % Genotype must be 'all' 'P' or 'W'
flags.SessionN = 'all'; % SessionN must be 'all' '1' '2' or '3'
[correctionLH,trlStruct] = getCorrectionsLH(trlStruct,figSavePath,masterTbl,flags);

%% Analysis 8. Drinking behaviors!
% These will initially be kept 2x2 (P Rats: Congruent v Incongruent &
% Wistars: Congruent v Incongruent since overall intakes are very different
% between P rats and Wistars. Therefore, currently no reason to throw flags
% into this function. However, am keeping it an option for later. 

[drinkingStats,trlStruct] = getDrinkingBehaviors(trlStruct,figSavePath,masterTbl,flags);

%% Analysis 9. Change points (Approach & Corrections).
% Utilizes change points saved from Analysis 3 and 7. Plots and analyzes
% these change points. 

flags.Genotype = 'all'; % Genotype must be 'all' 'P' or 'W'
flags.SessionN = 'all'; % SessionN must be 'all' '1' '2' or '3'
flags.figextension = 'png';

changePoints = getChangePoints(trlStruct,figSavePath,masterTbl,flags);

%% Analysis 10. Change points (over sessions).
% changePoints variable exists from Analysis 9. 
% flags.Genotype = 'all'; % Genotype must be 'all' 'P' or 'W'
% flags.SessionN = '3'; % SessionN must be 'all' '1' '2' or '3'
% flags.figextension = 'png';
% 
% changebySession = sessionChangePoints(changebySession,trlStruct,figSavePath,masterTbl,flags);

%% Analysis 11. Find Correction Latencies & Times.
flags.Genotype = 'all'; % Genotype must be 'all' 'P' or 'W'
flags.SessionN = 'all'; % SessionN must be 'all' '1' '2' or '3'
flags.figextension = 'png';

[trlStruct,Lat2Correct] = getLat2Correct(trlStruct,figSavePath,masterTbl,flags);

%% Analysis 12. Find characteristics of movement of nose during corrections
% Function now has nothing to do with angles. The animal did not exhibit
% that drastic of a change that I noticed. Will look again some day. 
% Function here, currently makes a line between the two sippers and
% calculates the distance between that line (thought to be the most
% efficient path available) and the animal's trajectory. 

flags.Genotype = 'all'; % Genotype must be 'all' 'P' or 'W'
flags.SessionN = 'all'; % SessionN must be 'all' '1' '2' or '3'
flags.figextension = 'png';

[trlStruct,angStruct] = getAng(trlStruct,figSavePath,masterTbl,flags);

%% Analysis 13. Find nature of incorrect approaches.
flags.Genotype = 'all'; % Genotype must be 'all' 'P' or 'W'
flags.SessionN = 'all'; % SessionN must be 'all' '1' '2' or '3'
flags.figextension = 'png';

[trlStruct,errorInfoSt] = getErrorInfo(trlStruct,masterTbl,figSavePath,flags);

%% Analysis 14. Approach Likelihood, revisited & corrected.
flags.Genotype = 'P'; % Genotype must be 'all' 'P' or 'W'
flags.SessionN = 'all'; % SessionN must be 'all' '1' '2' or '3'
flags.combined = 2;     % 1 means combine correct and incorrect approaches, 0 means separate 

[approachStatistics,trlStruct] = getApproachLH_V2(trlStruct,masterTbl,figSavePath,flags);
statTable = [];
statTable = [approachStatistics.correctCongruent approachStatistics.incorrectCongruent approachStatistics.correctIncongruent approachStatistics.incorrectIncongruent];
% statTable = [approachStatistics.correctIncongruent approachStatistics.incorrectIncongruent];

statTable = statTable';
statTable = array2table(statTable);

sessionLabel = [repmat({'Congruent'},size(approachStatistics.correctCongruent,2),2); repmat({'Incongruent'},size(approachStatistics.correctIncongruent,2),2)]; sessionLabel = sessionLabel'; sessionLabel = sessionLabel(:);
trialLabel = [repmat({'Correct'},size(approachStatistics.correctCongruent,2),2); repmat({'Incorrect'},size(approachStatistics.incorrectCongruent,2),2)]; trialLabel = trialLabel(:);
% trialLabel = [repmat({'Correct'},size(approachStatistics.correctCongruent,2),1); repmat({'Incorrect'},size(approachStatistics.incorrectCongruent,2),1)]; trialLabel = trialLabel(:);

IDMat = repmat(1:size(approachStatistics.correctIncongruent,2),1,4); IDMat = IDMat(:);
% IDMat = repmat(1:size(approachStatistics.correctIncongruent,2),1,2); IDMat = IDMat(:);

statTable.SessionType = sessionLabel;
statTable.trialType = trialLabel;
statTable.ID = IDMat;
Trial = table([1:48]','VariableNames',{'Trials'});

rm = fitrm(statTable,'statTable1-statTable48 ~ trialType + SessionType + SessionType * trialType','WithinDesign',Trial);
anovaTbl = ranova(rm,'WithinModel','Trials')
multcompare(rm,'trialType')


tblStr = formattedDisplayText(anovaTbl); 
fid = fopen([flags.Genotype flags.SessionN 'anovaTable_ranova_approach.txt'], 'wt');
fileCleanup = onCleanup(@()fclose(fid));
formatSpec = '%s\n';
fprintf(fid, formatSpec, tblStr);
clear('fileCleanup')

%%
for i = 1:min(size(approachStatistics.correctCongruent))
    correctCong_ipt(i) = findchangepts(approachStatistics.correctCongruent(:,i));
    correctIncong_ipt(i) = findchangepts(approachStatistics.correctIncongruent(:,i));

    incorrectCong_ipt(i) = findchangepts(approachStatistics.incorrectCongruent(:,i));
    incorrectIncong_ipt(i) = findchangepts(approachStatistics.incorrectIncongruent(:,i));
end
%% Plot Stacked Bar Plot where the proportions are Correct, Incorrect, and
% Omissions for Congruent and Incongruent sessions for first 15 trials

congruent = [mean(approachStatistics.correctCongruent(1:15,:),2), mean(approachStatistics.incorrectCongruent(1:15,:),2), (1 - (mean(approachStatistics.correctCongruent(1:15,:),2) + mean(approachStatistics.incorrectCongruent(1:15,:),2)))];
incongruent = [mean(approachStatistics.correctIncongruent(1:15,:),2), mean(approachStatistics.incorrectIncongruent(1:15,:),2), (1 - (mean(approachStatistics.correctIncongruent(1:15,:),2) + mean(approachStatistics.incorrectIncongruent(1:15,:),2)))];

% Plot type 1
% data = cat(3,congruent,incongruent);
% data = permute(data,[1 3 2]);
% plotBarStackGroups(data,{1:15});
% ylim([0 1.1])
% xlabel('Trials')
% ylabel('Proportion of Trials')
% legend({'Correct','Incorrect','Omissions'})

% Plot type 2
figure('Units','normalized','Position',[0 0 1 1])
subplot(1,2,1)
bar(congruent,'Stacked')
title('Congruent Sessions')
xlabel('Trials')
ylabel('Proportion of Trials')
legend({'Correct','Incorrect','Omissions'})
ylim([0 1.1])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
subplot(1,2,2)
bar(incongruent,'Stacked')
title('Incongruent Sessions')
xlabel('Trials')
ylabel('Proportion of Trials')
legend({'Correct','Incorrect','Omissions'})
ylim([0 1.1])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

saveas(gca,[figSavePath filesep 'StackedBarLikelihoods_v2_' flags.Genotype '_100ms'],'svg');     saveas(gca,[figSavePath filesep 'StackedBarLikelihoods_v2_' flags.Genotype '_100ms'],'fig')
%% Attempt to do stats on stacked bar graph

dataConNumP = floor(congruent(1,:)*48);
dataIncNumP = floor(incongruent(1,:)*48);
X = [dataConNumP(:),dataIncNumP(:)];
[h,p,XSquared] = chi2cont(X)


%% Analysis 15. Sipper Occupancy Likelihood. Not all trials
flags.Genotype = 'all'; % Genotype must be 'all' 'P' or 'W'
flags.SessionN = 'all'; % SessionN must be 'all' '1' '2' or '3'

[sipperOccupancyStatistics,trlStruct] = getSipperLH_AllTrials(trlStruct,figSavePath,masterTbl,flags);


%% Analysis 16. Same as 15 but incorrect trials.
flags.Genotype = 'all'; % Genotype must be 'all' 'P' or 'W'
flags.SessionN = 'all'; % SessionN must be 'all' '1' '2' or '3'

[sipperOccupancyStatistics,trlStruct] = getSipperLH_AllTrials_IncorrectSipper(trlStruct,figSavePath,masterTbl,flags);

%% Analysis 17. Velocity information of nose per trial 
flags.Genotype = 'all'; % Genotype must be 'all' 'P' or 'W'
flags.SessionN = 'all'; % SessionN must be 'all' '1' '2' or '3'

[trlStruct,~] = getVelocity(trlStruct,figSavePath,masterTbl,flags);

%% Analysis 18. CSDiscrim information per trial 
flags.Genotype = 'all'; % Genotype must be 'all' 'P' or 'W'
flags.SessionN = 'all'; % SessionN must be 'all' '1' '2' or '3'

[~,trlStruct] = getCSdiscrim(trlStruct,masterTbl,flags);
% %% Ephys Structure
% checkPath = 'F:/dissDat/dataStructures/ephysStruct.mat'; % Path to check existence of trlStruct variable. Preexisting variables will be loaded.
% 
% if ~exist('ephysStruct','var')
%     ephysStruct = ephysStructLoadorMake(dataPath,checkPath,savePath,tablePath);
% else
%     disp('ephysStruct.mat already exists in the workspace')
% end

%% Run script twofactoranalyses.m for graphing and statistics
run('twofactoranalyses.m')

%% Behavioral analyses ends here
% Run scripts associated with electrophysiology next
% Inside the script there are paths that need changed.
% Additionally, on line 21 change: sessionType = 'Incongruent' to 
% sessionType = 'Congruent' and vis versa depending on which session type
% you would like to look at.
% 
% run("LR_ApproachAnalysis_FinalScript.m")
% run("plotDistLat.m")
% run("cueEphysAnalysisScript.m")