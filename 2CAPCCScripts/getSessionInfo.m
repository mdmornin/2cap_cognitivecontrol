masterTbl = readtable('DataLog3.xlsx','Sheet','RecordingLog');
reversalSessions = strcmp(masterTbl.SessionType, 'Reversal');
numReversalSessions = reversalSessions(reversalSessions == 1);

%% Check existence of local datasets
dataDir = dir;
folderNames = {dataDir([dataDir.isdir]).name};
folderNames = folderNames(~ismember(folderNames ,{'.','..'}));

dataSetIDs = masterTbl.DataSetName;

matchingDatasets = strcmp(dataSetIDs,folderNames);

% %% Input Bottle In/Out and Animal Weight into masterTbl
% path = 'H:\dissDat\Data\';
% varInterest = [{'Initial Left Bottle Weight (g)'},{'Initial Right Bottle Weight (g)'},{'Animal Weight (g)'},{'Final Left Bottle Weight (g)'},{'Final Right Bottle Weight (g)'}];
% for j = 1:length(dataSetIDs)
%     folderName = dataSetIDs{j};
%     filePath = [path folderName '\sessionInfo.xlsx'];
%     sessionInfo = readtable(filePath,'Sheet','Sheet1');
%     for k = 1:length(varInterest)
%         idx = find(strcmp(sessionInfo.Var1,varInterest{k}));
%         columnName = varInterest{k};
%         masterTbl.(columnName)(j) = str2double(sessionInfo.Var2(idx));
%     end
% end 
%% Check existence of DLC Results
% Adds 1 if .csv exists in file directory versus does not exist (0)
affix = 'DLC_resnet50_2CAPFullRun1Mar28shuffle1_1030000.csv';
for j = 1:length(dataSetIDs)
    strSplit = split(dataSetIDs{j});
    prefix = [strSplit{1} strSplit{3} '_video'];
    dlcFileName = [path dataSetIDs{j} '\' prefix affix];
    if exist(dlcFileName,'file') == 2
        masterTbl.RGBVideoInDeepLabCut_{j} = 1;
    else
        masterTbl.RGBVideoInDeepLabCut_{j} = 0;
    end
end
%% Calculate Alcohol Intake (g/kg)
masterTbl.RightBottleDiff = masterTbl.("Initial Right Bottle Weight (g)") - masterTbl.("Final Right Bottle Weight (g)");
masterTbl.LeftBottleDiff = masterTbl.("Initial Left Bottle Weight (g)") - masterTbl.("Final Left Bottle Weight (g)");
masterTbl.Intake = (masterTbl.RightBottleDiff + masterTbl.LeftBottleDiff) ./ (masterTbl.("Animal Weight (g)")/1000) * 0.1;

%% Intake Per Session Type

regIdx = strcmp(masterTbl.SessionType,'Regular');
    regIntake = masterTbl.Intake(regIdx);
revIdx = strcmp(masterTbl.SessionType,'Reversal');
    revIntake = masterTbl.Intake(revIdx);

figure
boxplot([regIntake,revIntake],[{'Regular'},{'Reversal'}],'notch','on')
xlabel('Session Type')
ylabel('Intake (g/kg)')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold')
ylim([0 3])


%% Intake Per Session Type & Strain
regPIdx = strcmp(masterTbl.SessionType,'Regular') & strcmp(masterTbl.Strain,'P');
revPIdx = strcmp(masterTbl.SessionType,'Reversal') & strcmp(masterTbl.Strain,'P');

regWIdx = strcmp(masterTbl.SessionType,'Regular') & strcmp(masterTbl.Strain,'W');
revWIdx = strcmp(masterTbl.SessionType,'Reversal') & strcmp(masterTbl.Strain,'W');

intakePReg = masterTbl.Intake(regPIdx);
intakePRev = masterTbl.Intake(revPIdx);

intakeWReg = masterTbl.Intake(regWIdx);
intakeWRev = masterTbl.Intake(revWIdx);

g1 = repmat({'Regular P'},length(intakePReg),1);
g2 = repmat({'Reversal P'},length(intakePRev),1);
g3 = repmat({'Regular W'},length(intakeWReg),1);
g4 = repmat({'Reversal W'},length(intakeWRev),1);

g = [g1; g2; g3; g4];

figure
boxplot([intakePReg;intakePRev;intakeWReg;intakeWRev],g,'notch','on')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold')
xlabel('Strain & Session Type')
ylabel('Intake (g/kg)')
ylim([0 3])


