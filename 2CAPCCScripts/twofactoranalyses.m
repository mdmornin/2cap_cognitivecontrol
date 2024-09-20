close all;
%% Perserverations on Incorrect Approaches
% Indices for congruent, incongruent, P, Wistars
regPIdx = startsWith(masterTbl.SessionType,'Regular') & strcmp(masterTbl.Strain,'P');
revPIdx = startsWith(masterTbl.SessionType,'Reversal') & strcmp(masterTbl.Strain,'P');
    
regWIdx = startsWith(masterTbl.SessionType,'Regular') & strcmp(masterTbl.Strain,'W');
revWIdx = startsWith(masterTbl.SessionType,'Reversal') & strcmp(masterTbl.Strain,'W');

errorPersPCon = [trlStruct(regPIdx).Pers];
errorPersPInc = [trlStruct(revPIdx).Pers];

errorPersWCon = [trlStruct(regWIdx).Pers];
errorPersWInc = [trlStruct(revWIdx).Pers];
% Plotting

% Sem Intakes
semerrorPersPCon = std(errorPersPCon)/sqrt(length(errorPersPCon));
semerrorPersPInc = std(errorPersPInc)/sqrt(length(errorPersPInc));
semerrorPersWCon= std(errorPersWCon)/sqrt(length(errorPersWCon));
semerrorPersWInc = std(errorPersWInc)/sqrt(length(errorPersWInc));

% Plot figure
figure('Units','normalized','Position',[0 0 1 1])
hb = bar(categorical({'P Rats','Wistars'}),[mean(errorPersPCon) mean(errorPersPInc); mean(errorPersWCon) mean(errorPersWInc)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);
hb(1).FaceColor = [0.5 0.5 0.5]; hb(2).FaceColor = [0.9 0.9 0.9];
offsetPos = [1+hb(1).XOffset 1+hb(2).XOffset 2+hb(1).XOffset 2+hb(2).XOffset];
scatterXData = [linspace(offsetPos(1)-0.1,offsetPos(1)+0.1,length(errorPersPCon)),linspace(offsetPos(2)-0.1,offsetPos(2)+0.1,length(errorPersPInc)),linspace(offsetPos(3)-0.1,offsetPos(3)+0.1,length(errorPersWCon)),linspace(offsetPos(4)-0.1,offsetPos(4)+0.1,length(errorPersWInc))];
hold on
plot([linspace(offsetPos(1)-0.05,offsetPos(1)+0.05,length(errorPersPCon)) linspace(offsetPos(2)-0.05,offsetPos(2)+0.05,length(errorPersPInc))]',[errorPersPCon errorPersPInc],'ko','LineWidth',2,'MarkerSize',8)
plot([linspace(offsetPos(3)-0.05,offsetPos(3)+0.05,length(errorPersWCon)) linspace(offsetPos(4)-0.05,offsetPos(4)+0.05,length(errorPersWInc))]',[errorPersWCon errorPersWInc],'ko','LineWidth',2,'MarkerSize',8)

er = errorbar(offsetPos,[mean(errorPersPCon) mean(errorPersPInc) mean(errorPersWCon) mean(errorPersWInc)],[semerrorPersPCon semerrorPersPInc semerrorPersWCon semerrorPersWInc],'LineWidth',3);    
er(1).Color = [0 0 0];                        
er(1).LineStyle = 'none'; 
legend([{'Congruent'},{'Incongruent'}])
ylabel('Errors in a Row')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
% Save figure
saveas(gca,[figSavePath filesep 'twobyPers_'],'svg')
%% Measure number of perserveration bouts, defined as >=3 in a row.
for i = 1:length(trlStruct)
    trlStruct(i).PersBouts = sum(trlStruct(i).Pers >= 3);
end
persBoutPCon = [trlStruct(regPIdx).PersBouts];
persBoutPInc = [trlStruct(revPIdx).PersBouts];
persBoutWCon = [trlStruct(regWIdx).PersBouts];
persBoutWInc = [trlStruct(revWIdx).PersBouts];

% persBoutPCon = persBoutPCon(persBoutPCon > 0);
% persBoutPInc = persBoutPInc(persBoutPInc > 0);
% 
% persBoutWCon = persBoutWCon(persBoutWCon > 0);
% persBoutWInc = persBoutWInc(persBoutWInc > 0);

% Sem 
sempersBoutPCon= std(persBoutPCon)/sqrt(length(persBoutPCon));
sempersBoutPInc = std(persBoutPInc)/sqrt(length(persBoutPInc));
sempersBoutWCon = std(persBoutWCon)/sqrt(length(persBoutWCon));
sempersBoutWInc = std(persBoutWInc)/sqrt(length(persBoutWInc));

% Plot figure
figure('Units','normalized','Position',[0 0 1 1])
hb = bar(categorical({'P Rats','Wistars'}),[mean(persBoutPCon) mean(persBoutPInc); mean(persBoutWCon) mean(persBoutWInc)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);
hb(1).FaceColor = [0.5 0.5 0.5]; hb(2).FaceColor = [0.9 0.9 0.9];
offsetPos = [1+hb(1).XOffset 1+hb(2).XOffset 2+hb(1).XOffset 2+hb(2).XOffset];
scatterXData = [linspace(offsetPos(1)-0.1,offsetPos(1)+0.1,length(persBoutPCon)),linspace(offsetPos(2)-0.1,offsetPos(2)+0.1,length(persBoutPInc)),linspace(offsetPos(3)-0.1,offsetPos(3)+0.1,length(persBoutWCon)),linspace(offsetPos(4)-0.1,offsetPos(4)+0.1,length(persBoutWInc))];
hold on
plot([linspace(offsetPos(1)-0.05,offsetPos(1)+0.05,length(persBoutPCon)) linspace(offsetPos(2)-0.05,offsetPos(2)+0.05,length(persBoutPInc))]',[persBoutPCon persBoutPInc],'ko','LineWidth',2,'MarkerSize',8)
plot([linspace(offsetPos(3)-0.05,offsetPos(3)+0.05,length(persBoutWCon)) linspace(offsetPos(4)-0.05,offsetPos(4)+0.05,length(persBoutWInc))]',[persBoutWCon persBoutWInc],'ko','LineWidth',2,'MarkerSize',8)

er = errorbar(offsetPos,[mean(persBoutPCon) mean(persBoutPInc) mean(persBoutWCon) mean(persBoutWInc)],[sempersBoutPCon sempersBoutPInc sempersBoutWCon sempersBoutWInc],'LineWidth',3);    
er(1).Color = [0 0 0];                        
er(1).LineStyle = 'none'; 
legend([{'Congruent'},{'Incongruent'}])
ylabel('Perseveration Bouts')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
saveas(gca,[figSavePath filesep 'twobyPersBouts_'],'svg')

%% Difference in Side where Errors Occur
errorDPCon = abs([trlStruct(regPIdx).errorD]);
errorDPInc = abs([trlStruct(revPIdx).errorD]);

errorDWCon = abs([trlStruct(regWIdx).errorD]);
errorDWInc = abs([trlStruct(revWIdx).errorD]);

% SEM values
semerrorDPCon = std(errorDPCon)/sqrt(length(errorDPCon));
semerrorDPInc = std(errorDPInc)/sqrt(length(errorDPInc));
semerrorDWCon = std(errorDWCon)/sqrt(length(errorDWCon));
semerrorDWInc = std(errorDWInc)/sqrt(length(errorDWInc));

% Plot figure
figure('Units','normalized','Position',[0 0 1 1])
hb = bar(categorical({'P Rats','Wistars'}),[mean(errorDPCon) mean(errorDPInc); mean(errorDWCon) mean(errorDWInc)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);
hb(1).FaceColor = [0.5 0.5 0.5]; hb(2).FaceColor = [0.9 0.9 0.9];
offsetPos = [1+hb(1).XOffset 1+hb(2).XOffset 2+hb(1).XOffset 2+hb(2).XOffset];
scatterXData = [linspace(offsetPos(1)-0.1,offsetPos(1)+0.1,length(errorDPCon)),linspace(offsetPos(2)-0.1,offsetPos(2)+0.1,length(errorDPInc)),linspace(offsetPos(3)-0.1,offsetPos(3)+0.1,length(errorDWCon)),linspace(offsetPos(4)-0.1,offsetPos(4)+0.1,length(errorDWInc))];
hold on
plot([linspace(offsetPos(1)-0.05,offsetPos(1)+0.05,length(errorDPCon)) linspace(offsetPos(2)-0.05,offsetPos(2)+0.05,length(errorDPInc))]',[errorDPCon errorDPInc],'ko','LineWidth',2,'MarkerSize',8)
plot([linspace(offsetPos(3)-0.05,offsetPos(3)+0.05,length(errorDWCon)) linspace(offsetPos(4)-0.05,offsetPos(4)+0.05,length(errorDWInc))]',[errorDWCon errorDWInc],'ko','LineWidth',2,'MarkerSize',8)

er = errorbar(offsetPos,[mean(errorDPCon) mean(errorDPInc) mean(errorDWCon) mean(errorDWInc)],[semerrorDPCon semerrorDPInc semerrorDWCon semerrorDWInc],'LineWidth',3);    
er(1).Color = [0 0 0];                        
er(1).LineStyle = 'none'; 
legend([{'Congruent'},{'Incongruent'}])
ylabel('Side Bias')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
saveas(gca,[figSavePath filesep 'twobyErrorD_'],'svg')

%% Latency to Correct Data
corLatPCon = [trlStruct(regPIdx).correctionLatsD]./30;
corLatPInc = [trlStruct(revPIdx).correctionLatsD]./30;
corLatWCon = [trlStruct(regWIdx).correctionLatsD]./30;
corLatWInc = [trlStruct(revWIdx).correctionLatsD]./30;

% SEM values
semcorLatPCon = std(corLatPCon)/sqrt(length(corLatPCon));
semcorLatPInc = std(corLatPInc)/sqrt(length(corLatPInc));
semcorLatWCon = std(corLatWCon)/sqrt(length(corLatWCon));
semcorLatWInc = std(corLatWInc)/sqrt(length(corLatWInc));

% Plot figure
figure('Units','normalized','Position',[0 0 1 1])
hb = bar(categorical({'P Rats','Wistars'}),[mean(corLatPCon) mean(corLatPInc); mean(corLatWCon) mean(corLatWInc)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);
hb(1).FaceColor = [0.5 0.5 0.5]; hb(2).FaceColor = [0.9 0.9 0.9];
offsetPos = [1+hb(1).XOffset 1+hb(2).XOffset 2+hb(1).XOffset 2+hb(2).XOffset];
scatterXData = [linspace(offsetPos(1)-0.1,offsetPos(1)+0.1,length(corLatPCon)),linspace(offsetPos(2)-0.1,offsetPos(2)+0.1,length(corLatPInc)),linspace(offsetPos(3)-0.1,offsetPos(3)+0.1,length(corLatWCon)),linspace(offsetPos(4)-0.1,offsetPos(4)+0.1,length(corLatWInc))];
hold on
plot([linspace(offsetPos(1)-0.05,offsetPos(1)+0.05,length(corLatPCon)) linspace(offsetPos(2)-0.05,offsetPos(2)+0.05,length(corLatPInc))]',[corLatPCon corLatPInc],'ko','LineWidth',2,'MarkerSize',8)
plot([linspace(offsetPos(3)-0.05,offsetPos(3)+0.05,length(corLatWCon)) linspace(offsetPos(4)-0.05,offsetPos(4)+0.05,length(corLatWInc))]',[corLatWCon corLatWInc],'ko','LineWidth',2,'MarkerSize',8)

er = errorbar(offsetPos,[mean(corLatPCon) mean(corLatPInc) mean(corLatWCon) mean(corLatWInc)],[semcorLatPCon semcorLatPInc semcorLatWCon semcorLatWInc],'LineWidth',3);    
er(1).Color = [0 0 0];                        
er(1).LineStyle = 'none'; 
legend([{'Congruent'},{'Incongruent'}])
ylabel('Latency to Correct (s)')
ylim([0 10])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
saveas(gca,[figSavePath filesep 'twobycorrectionLatenciesD_Session1'],'svg')

%% Correct, Incorrect, Corrections Statistics
regPIdx = startsWith(masterTbl.SessionType,'Regular') & strcmp(masterTbl.Strain,'P');
revPIdx = startsWith(masterTbl.SessionType,'Reversal') & strcmp(masterTbl.Strain,'P');
    
regWIdx = startsWith(masterTbl.SessionType,'Regular') & strcmp(masterTbl.Strain,'W');
revWIdx = startsWith(masterTbl.SessionType,'Reversal') & strcmp(masterTbl.Strain,'W');

% Pull all data desired
corPCon = [trlStruct(regPIdx).CorrectNoCorrections];     incorPCon = [trlStruct(regPIdx).Incorrect];       omitPCon = [trlStruct(regPIdx).Omissions];       correctionsPCon = [trlStruct(regPIdx).correctionRaw];
corPInc = [trlStruct(revPIdx).CorrectNoCorrections];     incorPInc = [trlStruct(revPIdx).Incorrect];       omitPInc = [trlStruct(revPIdx).Omissions];       correctionsPInc = [trlStruct(revPIdx).correctionRaw];
corWCon = [trlStruct(regWIdx).CorrectNoCorrections];     incorWCon = [trlStruct(regWIdx).Incorrect];       omitWCon = [trlStruct(regWIdx).Omissions];       correctionsWCon = [trlStruct(regWIdx).correctionRaw];
corWInc = [trlStruct(revWIdx).CorrectNoCorrections];     incorWInc = [trlStruct(revWIdx).Incorrect];       omitWInc = [trlStruct(revWIdx).Omissions];       correctionsWInc = [trlStruct(revWIdx).correctionRaw];

correctionsPCon = correctionsPCon./(48 - omitPCon); correctionsPInc = correctionsPInc./(48 - omitPInc);
correctionsWCon = correctionsWCon./(48 - omitWCon); correctionsWInc = correctionsWInc./(48 - omitWInc);
intakePReg = masterTbl.Intake(regPIdx); intakeWReg = masterTbl.Intake(regWIdx); intakeWRev = masterTbl.Intake(revWIdx); intakePRev = masterTbl.Intake(revPIdx);





% Calculate SEM values

semcorPCon = std(corPCon)/sqrt(length(corPCon));    semincorPCon = std(incorPCon)/sqrt(length(incorPCon));  semomitPCon = std(omitPCon)/sqrt(length(omitPCon));  semcorrectionsPCon = std(correctionsPCon)/sqrt(length(correctionsPCon));  
semcorPInc = std(corPInc)/sqrt(length(corPInc));    semincorPInc = std(incorPInc)/sqrt(length(incorPInc));  semomitPInc = std(omitPInc)/sqrt(length(omitPInc));  semcorrectionsPInc = std(correctionsPInc)/sqrt(length(correctionsPInc));
semcorWCon = std(corWCon)/sqrt(length(corWCon));    semincorWCon = std(incorWCon)/sqrt(length(incorWCon));  semomitWCon = std(omitWCon)/sqrt(length(omitWCon));  semcorrectionsWCon = std(correctionsWCon)/sqrt(length(correctionsWCon));
semcorWInc = std(corWInc)/sqrt(length(corWInc));    semincorWInc = std(incorWInc)/sqrt(length(incorWInc));  semomitWInc = std(omitWInc)/sqrt(length(omitWInc));  semcorrectionsWInc = std(correctionsWInc)/sqrt(length(correctionsWInc));
% Statistics
% Correct approach ANOVA

testData = [corPCon corPInc corWCon corWInc];
g1 = [repmat({'P rat'},1,length([corPCon corPInc])) repmat({'Wistar'},1,length([corWCon corWInc]))];
g2 = [repmat({'Congruent'},1,length(corPCon)) repmat({'Incongruent'},1,length(corPInc)) repmat({'Congruent'},1,length(corWCon)) repmat({'Incongruent'},1,length(corWInc))];
if length(testData) == length(g1) && length(testData) == length(g2)
    disp('Data sizes are matched')
    [p,tbl,stat] = anovan(testData,{g1,g2},'model','full');
    multcompare(stat,'Dimension',[1 2]);

    tblStr = formattedDisplayText(tbl); 
    % Write string to file
    fid = fopen('anovaTable_2x2_correct.txt', 'wt');
    fileCleanup = onCleanup(@()fclose(fid));
    formatSpec = '%s\n';
    fprintf(fid, formatSpec, tblStr);
    clear('fileCleanup')
end



% Incorrect approach ANOVA
testData = [incorPCon incorPInc incorWCon incorWInc];
g1 = [repmat({'P rat'},1,length([incorPCon incorPInc])) repmat({'Wistar'},1,length([incorWCon incorWInc]))];
g2 = [repmat({'Congruent'},1,length(incorPCon)) repmat({'Incongruent'},1,length(incorPInc)) repmat({'Congruent'},1,length(incorWCon)) repmat({'Incongruent'},1,length(incorWInc))];
if length(testData) == length(g1) && length(testData) == length(g2)
    disp('Data sizes are matched')
    [p,tbl,stat] = anovan(testData,{g1,g2},'model','full');
    multcompare(stat,'Dimension',[1 2])
        % Write string to file
    tblStr = formattedDisplayText(tbl); 
    fid = fopen('anovaTable_2x2_incorrect.txt', 'wt');
    fileCleanup = onCleanup(@()fclose(fid));
    formatSpec = '%s\n';
    fprintf(fid, formatSpec, tblStr);
    clear('fileCleanup')
end

% Omissions approach ANOVA
testData = [omitPCon omitPInc omitWCon omitWInc];
g1 = [repmat({'P rat'},1,length([incorPCon incorPInc])) repmat({'Wistar'},1,length([incorWCon incorWInc]))];
g2 = [repmat({'Congruent'},1,length(incorPCon)) repmat({'Incongruent'},1,length(incorPInc)) repmat({'Congruent'},1,length(incorWCon)) repmat({'Incongruent'},1,length(incorWInc))];
if length(testData) == length(g1) && length(testData) == length(g2)
    disp('Data sizes are matched')
    [p,tbl,stat] = anovan(testData,{g1,g2},'model','full');
    multcompare(stat,'Dimension',[1 2])
        % Write string to file
    tblStr = formattedDisplayText(tbl); 
    fid = fopen('anovaTable_2x2_omissions.txt', 'wt');
    fileCleanup = onCleanup(@()fclose(fid));
    formatSpec = '%s\n';
    fprintf(fid, formatSpec, tblStr);
    clear('fileCleanup')
end

% Corrections approach ANOVA
testData = [correctionsPCon correctionsPInc correctionsWCon correctionsWInc];
g1 = [repmat({'P rat'},1,length([incorPCon incorPInc])) repmat({'Wistar'},1,length([incorWCon incorWInc]))];
g2 = [repmat({'Congruent'},1,length(incorPCon)) repmat({'Incongruent'},1,length(incorPInc)) repmat({'Congruent'},1,length(incorWCon)) repmat({'Incongruent'},1,length(incorWInc))];
if length(testData) == length(g1) && length(testData) == length(g2)
    disp('Data sizes are matched')
    [p,tbl,stat] = anovan(testData,{g1,g2},'model','full');
    multcompare(stat,'Dimension',[1 2])
        % Write string to file
    tblStr = formattedDisplayText(tbl); 
    fid = fopen('anovaTable_2x2_corrections.txt', 'wt');
    fileCleanup = onCleanup(@()fclose(fid));
    formatSpec = '%s\n';
    fprintf(fid, formatSpec, tblStr);
    clear('fileCleanup')
end

% Intake ANOVA
testData = [intakePReg; intakePRev; intakeWReg; intakeWRev];
g1 = [repmat({'P rat'},1,length([incorPCon incorPInc])) repmat({'Wistar'},1,length([incorWCon incorWInc]))];
g2 = [repmat({'Congruent'},1,length(incorPCon)) repmat({'Incongruent'},1,length(incorPInc)) repmat({'Congruent'},1,length(incorWCon)) repmat({'Incongruent'},1,length(incorWInc))];
if length(testData) == length(g1) && length(testData) == length(g2)
    disp('Data sizes are matched')
    [p,tbl,stat] = anovan(testData,{g1,g2},'model','full');
    multcompare(stat,'Dimension',[1 2])
        % Write string to file
    tblStr = formattedDisplayText(tbl); 
    fid = fopen('anovaTable_2x2_intake.txt', 'wt');
    fileCleanup = onCleanup(@()fclose(fid));
    formatSpec = '%s\n';
    fprintf(fid, formatSpec, tblStr);
    clear('fileCleanup')
end

% Plot 2x2 Correct Approaches
% Plot figure
figure('Units','normalized','Position',[0 0 1 1])
subplot(1,2,1)
hb = bar(categorical({'P Rats','Wistars'}),[mean(corPCon) mean(corPInc); mean(corWCon) mean(corWInc)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);
hb(1).FaceColor = [0.5 0.5 0.5]; hb(2).FaceColor = [0.9 0.9 0.9];
offsetPos = [1+hb(1).XOffset 1+hb(2).XOffset 2+hb(1).XOffset 2+hb(2).XOffset];
scatterXData = [linspace(offsetPos(1)-0.1,offsetPos(1)+0.1,length(corPCon)),linspace(offsetPos(2)-0.1,offsetPos(2)+0.1,length(corPInc)),linspace(offsetPos(3)-0.1,offsetPos(3)+0.1,length(corWCon)),linspace(offsetPos(4)-0.1,offsetPos(4)+0.1,length(corWInc))];
hold on
plot([linspace(offsetPos(1)-0.05,offsetPos(1)+0.05,length(corPCon)) linspace(offsetPos(2)-0.05,offsetPos(2)+0.05,length(corPInc))]',[corPCon corPInc],'ko','LineWidth',2,'MarkerSize',8)
plot([linspace(offsetPos(3)-0.05,offsetPos(3)+0.05,length(corWCon)) linspace(offsetPos(4)-0.05,offsetPos(4)+0.05,length(corWInc))]',[corWCon corWInc],'ko','LineWidth',2,'MarkerSize',8)

er = errorbar(offsetPos,[mean(corPCon) mean(corPInc) mean(corWCon) mean(corWInc)],[semcorPCon semcorPInc semcorWCon semcorWInc],'LineWidth',3);    
er(1).Color = [0 0 0];                        
er(1).LineStyle = 'none'; 
legend([{'Congruent'},{'Incongruent'}])
ylabel('Correct Trials')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)

% Plot 2x2 Incorrect Approaches % % % % % % % % % % % %
subplot(1,2,2)
hb = bar(categorical({'P Rats','Wistars'}),[mean(incorPCon) mean(incorPInc); mean(incorWCon) mean(incorWInc)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);
hb(1).FaceColor = [0.5 0.5 0.5]; hb(2).FaceColor = [0.9 0.9 0.9];
offsetPos = [1+hb(1).XOffset 1+hb(2).XOffset 2+hb(1).XOffset 2+hb(2).XOffset];
scatterXData = [linspace(offsetPos(1)-0.1,offsetPos(1)+0.1,length(incorPCon)),linspace(offsetPos(2)-0.1,offsetPos(2)+0.1,length(incorPInc)),linspace(offsetPos(3)-0.1,offsetPos(3)+0.1,length(incorWCon)),linspace(offsetPos(4)-0.1,offsetPos(4)+0.1,length(incorWInc))];
hold on
plot([linspace(offsetPos(1)-0.05,offsetPos(1)+0.05,length(incorPCon)) linspace(offsetPos(2)-0.05,offsetPos(2)+0.05,length(incorPInc))]',[incorPCon incorPInc],'ko','LineWidth',2,'MarkerSize',8)
plot([linspace(offsetPos(3)-0.05,offsetPos(3)+0.05,length(incorWCon)) linspace(offsetPos(4)-0.05,offsetPos(4)+0.05,length(incorWInc))]',[incorWCon incorWInc],'ko','LineWidth',2,'MarkerSize',8)

er = errorbar(offsetPos,[mean(incorPCon) mean(incorPInc) mean(incorWCon) mean(incorWInc)],[semincorPCon semincorPInc semincorWCon semincorWInc],'LineWidth',3);    
er(1).Color = [0 0 0];                        
er(1).LineStyle = 'none'; 
legend([{'Congruent'},{'Incongruent'}])
ylabel('Incorrect Trials')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
saveas(gca,[figSavePath filesep 'correct_incorrectApproach_2x2'],'svg')

% Plot 2x2 Omissions
figure('Units','normalized','Position',[0 0 1 1])
subplot(1,2,1)
hb = bar(categorical({'P Rats','Wistars'}),[mean(omitPCon) mean(omitPInc); mean(omitWCon) mean(omitWInc)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);
hb(1).FaceColor = [0.5 0.5 0.5]; hb(2).FaceColor = [0.9 0.9 0.9];
offsetPos = [1+hb(1).XOffset 1+hb(2).XOffset 2+hb(1).XOffset 2+hb(2).XOffset];
scatterXData = [linspace(offsetPos(1)-0.1,offsetPos(1)+0.1,length(omitPCon)),linspace(offsetPos(2)-0.1,offsetPos(2)+0.1,length(omitPInc)),linspace(offsetPos(3)-0.1,offsetPos(3)+0.1,length(omitWCon)),linspace(offsetPos(4)-0.1,offsetPos(4)+0.1,length(omitWInc))];
hold on
plot([linspace(offsetPos(1)-0.05,offsetPos(1)+0.05,length(omitPCon)) linspace(offsetPos(2)-0.05,offsetPos(2)+0.05,length(omitPInc))]',[omitPCon omitPInc],'ko','LineWidth',2,'MarkerSize',8)
plot([linspace(offsetPos(3)-0.05,offsetPos(3)+0.05,length(omitWCon)) linspace(offsetPos(4)-0.05,offsetPos(4)+0.05,length(omitWInc))]',[omitWCon omitWInc],'ko','LineWidth',2,'MarkerSize',8)

er = errorbar(offsetPos,[mean(omitPCon) mean(omitPInc) mean(omitWCon) mean(omitWInc)],[semomitPCon semomitPInc semomitWCon semomitWInc],'LineWidth',3);    
er(1).Color = [0 0 0];                        
er(1).LineStyle = 'none'; 
legend([{'Congruent'},{'Incongruent'}])
ylabel('Omissions')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
saveas(gca,[figSavePath filesep 'omissions_2x2'],'svg')

% Plot 2x2 Corrections
% correctionsPCon = correctionsPCon./(48); correctionsPInc = correctionsPInc./(48);
% correctionsWCon = correctionsWCon./(48); correctionsWInc = correctionsWInc./(48);

subplot(1,2,2)
hb = bar(categorical({'P Rats','Wistars'}),[mean(correctionsPCon) mean(correctionsPInc); mean(correctionsWCon) mean(correctionsWInc)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);
hb(1).FaceColor = [0.5 0.5 0.5]; hb(2).FaceColor = [0.9 0.9 0.9];
offsetPos = [1+hb(1).XOffset 1+hb(2).XOffset 2+hb(1).XOffset 2+hb(2).XOffset];
scatterXData = [linspace(offsetPos(1)-0.1,offsetPos(1)+0.1,length(correctionsPCon)),linspace(offsetPos(2)-0.1,offsetPos(2)+0.1,length(correctionsPInc)),linspace(offsetPos(3)-0.1,offsetPos(3)+0.1,length(correctionsWCon)),linspace(offsetPos(4)-0.1,offsetPos(4)+0.1,length(correctionsWInc))];
hold on
plot([linspace(offsetPos(1)-0.05,offsetPos(1)+0.05,length(correctionsPCon)) linspace(offsetPos(2)-0.05,offsetPos(2)+0.05,length(correctionsPInc))]',[correctionsPCon correctionsPInc],'ko','LineWidth',2,'MarkerSize',8)
plot([linspace(offsetPos(3)-0.05,offsetPos(3)+0.05,length(correctionsWCon)) linspace(offsetPos(4)-0.05,offsetPos(4)+0.05,length(correctionsWInc))]',[correctionsWCon correctionsWInc],'ko','LineWidth',2,'MarkerSize',8)

er = errorbar(offsetPos,[mean(correctionsPCon) mean(correctionsPInc) mean(correctionsWCon) mean(correctionsWInc)],[semcorrectionsPCon semcorrectionsPInc semcorrectionsWCon semcorrectionsWInc],'LineWidth',3);    
er(1).Color = [0 0 0];                        
er(1).LineStyle = 'none'; 
legend([{'Congruent'},{'Incongruent'}])
ylabel('Likelihood to Correct')
ylim([0 1])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
saveas(gca,[figSavePath filesep 'omissions_propcorrections_2x2'],'svg')

%% Latency to Make Correct or Incorrect Approach
regPIdx = startsWith(masterTbl.SessionType,'Regular') & strcmp(masterTbl.Strain,'P');
revPIdx = startsWith(masterTbl.SessionType,'Reversal') & strcmp(masterTbl.Strain,'P');
    
regWIdx = startsWith(masterTbl.SessionType,'Regular') & strcmp(masterTbl.Strain,'W');
revWIdx = startsWith(masterTbl.SessionType,'Reversal') & strcmp(masterTbl.Strain,'W');

for i = 1:length(trlStruct)
    trlStruct(i).CorrectLatenciesM = mean(trlStruct(i).CorrectLatencies);
    trlStruct(i).IncorrectLatenciesM = mean(trlStruct(i).IncorrectLatencies);
end
% Pull all data desired
corLatPCon = [trlStruct(regPIdx).CorrectLatenciesM];     incorLatPCon = [trlStruct(regPIdx).IncorrectLatenciesM];       
corLatPInc = [trlStruct(revPIdx).CorrectLatenciesM];     incorLatPInc = [trlStruct(revPIdx).IncorrectLatenciesM];       
corLatWCon = [trlStruct(regWIdx).CorrectLatenciesM];     incorLatWCon = [trlStruct(regWIdx).IncorrectLatenciesM];      
corLatWInc = [trlStruct(revWIdx).CorrectLatenciesM];     incorLatWInc = [trlStruct(revWIdx).IncorrectLatenciesM];       

% Calculate SEM values

semcorLatPCon = std(corLatPCon)/sqrt(length(corLatPCon));    semincorLatPCon = std(incorLatPCon)/sqrt(length(incorLatPCon));
semcorLatPInc = std(corLatPInc)/sqrt(length(corLatPInc));    semincorLatPInc = std(incorLatPInc)/sqrt(length(incorLatPInc));
semcorLatWCon = std(corLatWCon)/sqrt(length(corLatWCon));    semincorLatWCon = std(incorLatWCon)/sqrt(length(incorLatWCon));
semcorLatWInc = std(corLatWInc)/sqrt(length(corLatWInc));    semincorLatWInc = std(incorLatWInc)/sqrt(length(incorLatWInc)); 


% ANOVA ZONE
testData = [corLatPCon corLatPInc corLatWCon corLatWInc];
g1 = [repmat({'P rat'},1,length([corLatPCon corLatPInc])) repmat({'Wistar'},1,length([corLatWCon corLatWInc]))];
g2 = [repmat({'Congruent'},1,length(corLatPCon)) repmat({'Incongruent'},1,length(corLatPInc)) repmat({'Congruent'},1,length(corLatWCon)) repmat({'Incongruent'},1,length(corLatWInc))];

if length(testData) == length(g1) && length(testData) == length(g2)
    disp('Data sizes are matched')
    [p,tbl,stat] = anovan(testData,{g1,g2},'model','full');
    multcompare(stat,'Dimension',[2])
        % Write string to file
    tblStr = formattedDisplayText(tbl); 
    fid = fopen('anovaTable_2x2_correctlatencies.txt', 'wt');
    fileCleanup = onCleanup(@()fclose(fid));
    formatSpec = '%s\n';
    fprintf(fid, formatSpec, tblStr);
    clear('fileCleanup')
end

% ANOVA ZONE INCORRECT LATENCIES
testData = [incorLatPCon incorLatPInc incorLatWCon incorLatWInc];
g1 = [repmat({'P rat'},1,length([incorLatPCon incorLatPInc])) repmat({'Wistar'},1,length([incorLatWCon incorLatWInc]))];
g2 = [repmat({'Congruent'},1,length(incorLatPCon)) repmat({'Incongruent'},1,length(incorLatPInc)) repmat({'Congruent'},1,length(incorLatWCon)) repmat({'Incongruent'},1,length(incorLatWInc))];

if length(testData) == length(g1) && length(testData) == length(g2)
    disp('Data sizes are matched')
    [p,tbl,stat] = anovan(testData,{g1,g2},'model','full');
    multcompare(stat,'Dimension',[1 2])
        % Write string to file
    tblStr = formattedDisplayText(tbl); 
    fid = fopen('anovaTable_2x2_incorrectlatencies.txt', 'wt');
    fileCleanup = onCleanup(@()fclose(fid));
    formatSpec = '%s\n';
    fprintf(fid, formatSpec, tblStr);
    clear('fileCleanup')
end

% Plot 2x2 Latency to Correct Approach
figure('Units','normalized','Position',[0 0 1 1])
subplot(1,2,1)
hb = bar(categorical({'P Rats','Wistars'}),[mean(corLatPCon) mean(corLatPInc); mean(corLatWCon) mean(corLatWInc)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);
hb(1).FaceColor = [0.5 0.5 0.5]; hb(2).FaceColor = [0.9 0.9 0.9];
offsetPos = [1+hb(1).XOffset 1+hb(2).XOffset 2+hb(1).XOffset 2+hb(2).XOffset];
scatterXData = [linspace(offsetPos(1)-0.1,offsetPos(1)+0.1,length(corLatPCon)),linspace(offsetPos(2)-0.1,offsetPos(2)+0.1,length(corLatPInc)),linspace(offsetPos(3)-0.1,offsetPos(3)+0.1,length(corLatWCon)),linspace(offsetPos(4)-0.1,offsetPos(4)+0.1,length(corLatWInc))];
hold on
plot([linspace(offsetPos(1)-0.05,offsetPos(1)+0.05,length(corLatPCon)) linspace(offsetPos(2)-0.05,offsetPos(2)+0.05,length(corLatPInc))]',[corLatPCon corLatPInc],'ko','LineWidth',2,'MarkerSize',8)
plot([linspace(offsetPos(3)-0.05,offsetPos(3)+0.05,length(corLatWCon)) linspace(offsetPos(4)-0.05,offsetPos(4)+0.05,length(corLatWInc))]',[corLatWCon corLatWInc],'ko','LineWidth',2,'MarkerSize',8)

er = errorbar(offsetPos,[mean(corLatPCon) mean(corLatPInc) mean(corLatWCon) mean(corLatWInc)],[semcorLatPCon semcorLatPInc semcorLatWCon semcorLatWInc],'LineWidth',3);    
er(1).Color = [0 0 0];                        
er(1).LineStyle = 'none'; 
legend([{'Congruent'},{'Incongruent'}])
ylim([0 3])

ylabel('Latency to Approach (Correct Sipper) (s)')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)

% Plot 2x2 Latency to Incorrect Approach
subplot(1,2,2)
hb = bar(categorical({'P Rats','Wistars'}),[mean(incorLatPCon) mean(incorLatPInc); mean(incorLatWCon) mean(incorLatWInc)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);
hb(1).FaceColor = [0.5 0.5 0.5]; hb(2).FaceColor = [0.9 0.9 0.9];
offsetPos = [1+hb(1).XOffset 1+hb(2).XOffset 2+hb(1).XOffset 2+hb(2).XOffset];
scatterXData = [linspace(offsetPos(1)-0.1,offsetPos(1)+0.1,length(incorLatPCon)),linspace(offsetPos(2)-0.1,offsetPos(2)+0.1,length(incorLatPInc)),linspace(offsetPos(3)-0.1,offsetPos(3)+0.1,length(incorLatWCon)),linspace(offsetPos(4)-0.1,offsetPos(4)+0.1,length(incorLatWInc))];
hold on
plot([linspace(offsetPos(1)-0.05,offsetPos(1)+0.05,length(incorLatPCon)) linspace(offsetPos(2)-0.05,offsetPos(2)+0.05,length(incorLatPInc))]',[incorLatPCon incorLatPInc],'ko','LineWidth',2,'MarkerSize',8)
plot([linspace(offsetPos(3)-0.05,offsetPos(3)+0.05,length(incorLatWCon)) linspace(offsetPos(4)-0.05,offsetPos(4)+0.05,length(incorLatWInc))]',[incorLatWCon incorLatWInc],'ko','LineWidth',2,'MarkerSize',8)

er = errorbar(offsetPos,[mean(incorLatPCon) mean(incorLatPInc) mean(incorLatWCon) mean(incorLatWInc)],[semincorLatPCon semincorLatPInc semincorLatWCon semincorLatWInc],'LineWidth',3);    
er(1).Color = [0 0 0];                        
er(1).LineStyle = 'none'; 
legend([{'Congruent'},{'Incongruent'}])
ylim([0 3])
ylabel('Latency to Approach (Incorrect Sipper) (s)')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
saveas(gca,[figSavePath filesep 'correct_incorrectApproachLatency_2x2_allSessions'],'svg')

% %% 2x2 Maximum Distance from Straight Line Between Sippers During Correction
% regPIdx = startsWith(masterTbl.SessionType,'Regular') & strcmp(masterTbl.Strain,'P');
% revPIdx = startsWith(masterTbl.SessionType,'Reversal') & strcmp(masterTbl.Strain,'P');
%     
% regWIdx = startsWith(masterTbl.SessionType,'Regular') & strcmp(masterTbl.Strain,'W');
% revWIdx = startsWith(masterTbl.SessionType,'Reversal') & strcmp(masterTbl.Strain,'W');
% 
% % Pull all data desired
% maxDistPCon = [trlStruct(regPIdx).aucDistBestLine];         
% maxDistPInc = [trlStruct(revPIdx).aucDistBestLine];           
% maxDistWCon = [trlStruct(regWIdx).aucDistBestLine];           
% maxDistWInc = [trlStruct(revWIdx).aucDistBestLine];           
% 
% % Calculate SEM values
% 
% semmaxDistPCon = std(maxDistPCon)/sqrt(length(maxDistPCon));    
% semmaxDistPInc = std(maxDistPInc)/sqrt(length(maxDistPInc));    
% semmaxDistWCon = std(maxDistWCon)/sqrt(length(maxDistWCon));    
% semmaxDistWInc = std(maxDistWInc)/sqrt(length(maxDistWInc));   
% 
% % 2x2 Plot
% figure('Units','normalized','Position',[0 0 1 1])
% hb = bar(categorical({'P Rats','Wistars'}),[mean(maxDistPCon) mean(maxDistPInc); mean(maxDistWCon) mean(maxDistWInc)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);
% hb(1).FaceColor = [0.5 0.5 0.5]; hb(2).FaceColor = [0.9 0.9 0.9];
% offsetPos = [1+hb(1).XOffset 1+hb(2).XOffset 2+hb(1).XOffset 2+hb(2).XOffset];
% scatterXData = [linspace(offsetPos(1)-0.1,offsetPos(1)+0.1,length(maxDistPCon)),linspace(offsetPos(2)-0.1,offsetPos(2)+0.1,length(maxDistPInc)),linspace(offsetPos(3)-0.1,offsetPos(3)+0.1,length(maxDistWCon)),linspace(offsetPos(4)-0.1,offsetPos(4)+0.1,length(maxDistWInc))];
% hold on
% plot([linspace(offsetPos(1)-0.05,offsetPos(1)+0.05,length(maxDistPCon)) linspace(offsetPos(2)-0.05,offsetPos(2)+0.05,length(maxDistPInc))]',[maxDistPCon maxDistPInc],'ko','LineWidth',2,'MarkerSize',8)
% plot([linspace(offsetPos(3)-0.05,offsetPos(3)+0.05,length(maxDistWCon)) linspace(offsetPos(4)-0.05,offsetPos(4)+0.05,length(maxDistWInc))]',[maxDistWCon maxDistWInc],'ko','LineWidth',2,'MarkerSize',8)
% 
% er = errorbar(offsetPos,[mean(maxDistPCon) mean(maxDistPInc) mean(maxDistWCon) mean(maxDistWInc)],[semmaxDistPCon semmaxDistPInc semmaxDistWCon semmaxDistWInc],'LineWidth',3);    
% er(1).Color = [0 0 0];                        
% er(1).LineStyle = 'none'; 
% legend([{'Congruent'},{'Incongruent'}])
% ylabel('Deviation from Ideal Correction Path (pixels)')
% set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
% saveas(gca,[figSavePath filesep 'maxDistIdealLine_2x2_allSessions'],'svg')

%% 2x2 Sipper Occupancy Time (s)
% Total sipper occupancy time 
% Indexing
regPIdx = startsWith(masterTbl.SessionType,'Regular') & strcmp(masterTbl.Strain,'P');
revPIdx = startsWith(masterTbl.SessionType,'Reversal') & strcmp(masterTbl.Strain,'P');
    
regWIdx = startsWith(masterTbl.SessionType,'Regular') & strcmp(masterTbl.Strain,'W');
revWIdx = startsWith(masterTbl.SessionType,'Reversal') & strcmp(masterTbl.Strain,'W');
% Pulling data from structure
occPCon = [trlStruct(regPIdx).totalCorrectSipperTime];
occPInc = [trlStruct(revPIdx).totalCorrectSipperTime];
occWCon = [trlStruct(regWIdx).totalCorrectSipperTime];
occWInc = [trlStruct(revWIdx).totalCorrectSipperTime];

% Calculate SEM
semoccPCon = std(occPCon)/sqrt(length(occPCon));    
semoccPInc = std(occPInc)/sqrt(length(occPInc));    
semoccWCon = std(occWCon)/sqrt(length(occWCon));    
semoccWInc = std(occWInc)/sqrt(length(occWInc));  

% Plot 2x2
figure('Units','normalized','Position',[0 0 1 1])
subplot(1,2,1)
hb = bar(categorical({'P Rats','Wistars'}),[mean(occPCon) mean(occPInc); mean(occWCon) mean(occWInc)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);
hb(1).FaceColor = [0.5 0.5 0.5]; hb(2).FaceColor = [0.9 0.9 0.9];
offsetPos = [1+hb(1).XOffset 1+hb(2).XOffset 2+hb(1).XOffset 2+hb(2).XOffset];
scatterXData = [linspace(offsetPos(1)-0.1,offsetPos(1)+0.1,length(occPCon)),linspace(offsetPos(2)-0.1,offsetPos(2)+0.1,length(occPInc)),linspace(offsetPos(3)-0.1,offsetPos(3)+0.1,length(occWCon)),linspace(offsetPos(4)-0.1,offsetPos(4)+0.1,length(occWInc))];
hold on
plot([linspace(offsetPos(1)-0.05,offsetPos(1)+0.05,length(occPCon)) linspace(offsetPos(2)-0.05,offsetPos(2)+0.05,length(occPInc))]',[occPCon occPInc],'ko','LineWidth',2,'MarkerSize',8)
plot([linspace(offsetPos(3)-0.05,offsetPos(3)+0.05,length(occWCon)) linspace(offsetPos(4)-0.05,offsetPos(4)+0.05,length(occWInc))]',[occWCon occWInc],'ko','LineWidth',2,'MarkerSize',8)

er = errorbar(offsetPos,[mean(occPCon) mean(occPInc) mean(occWCon) mean(occWInc)],[semoccPCon semoccPInc semoccWCon semoccWInc],'LineWidth',3);    
er(1).Color = [0 0 0];                        
er(1).LineStyle = 'none'; 
legend([{'Congruent'},{'Incongruent'}])
ylabel('Time at Sipper (s)')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
saveas(gca,[figSavePath filesep 'totalSipperOccupancy_2x2_allSessions'],'svg')

% ANOVA ZONE SIPPER OCCUPANCY
testData = [occPCon occPInc occWCon occWInc];
g1 = [repmat({'P rat'},1,length([occPCon occPInc])) repmat({'Wistar'},1,length([occWCon occWInc]))];
g2 = [repmat({'Congruent'},1,length(occPCon)) repmat({'Incongruent'},1,length(occPInc)) repmat({'Congruent'},1,length(occWCon)) repmat({'Incongruent'},1,length(occWInc))];

if length(testData) == length(g1) && length(testData) == length(g2)
    disp('Data sizes are matched')
    [p,tbl,stat] = anovan(testData,{g1,g2},'model','full');
    multcompare(stat,'Dimension',[2])
        % Write string to file
    tblStr = formattedDisplayText(tbl); 
    fid = fopen('anovaTable_2x2_sipperoccupancy.txt', 'wt');
    fileCleanup = onCleanup(@()fclose(fid));
    formatSpec = '%s\n';
    fprintf(fid, formatSpec, tblStr);
    clear('fileCleanup')
end

% ANOVA ZONE INTAKE RATE
occPCon = [trlStruct(regPIdx).intakeRate];
occPInc = [trlStruct(revPIdx).intakeRate];
occWCon = [trlStruct(regWIdx).intakeRate];
occWInc = [trlStruct(revWIdx).intakeRate];

testData = [occPCon occPInc occWCon occWInc];
g1 = [repmat({'P rat'},1,length([occPCon occPInc])) repmat({'Wistar'},1,length([occWCon occWInc]))];
g2 = [repmat({'Congruent'},1,length(occPCon)) repmat({'Incongruent'},1,length(occPInc)) repmat({'Congruent'},1,length(occWCon)) repmat({'Incongruent'},1,length(occWInc))];

if length(testData) == length(g1) && length(testData) == length(g2)
    disp('Data sizes are matched')
    [p,tbl,stat] = anovan(testData,{g1,g2},'model','full');
    multcompare(stat,'Dimension',[2])
        % Write string to file
    tblStr = formattedDisplayText(tbl); 
    fid = fopen('anovaTable_2x2_intakerate.txt', 'wt');
    fileCleanup = onCleanup(@()fclose(fid));
    formatSpec = '%s\n';
    fprintf(fid, formatSpec, tblStr);
    clear('fileCleanup')
end


%% 2x2 Likelihood to Occupy Sipper 
corcorPCon = [trlStruct(regPIdx).llOccupy_CC];
corcorPInc = [trlStruct(revPIdx).llOccupy_CC];
corcorWCon = [trlStruct(regWIdx).llOccupy_CC];
corcorWInc = [trlStruct(revWIdx).llOccupy_CC];

% Calculate SEM
semcorcorPCon = std(corcorPCon')/sqrt(min(size(corcorPCon)));
semcorcorPInc = std(corcorPInc')/sqrt(min(size(corcorPInc)));
semcorcorWCon = std(corcorWCon')/sqrt(min(size(corcorWCon)));
semcorcorWInc = std(corcorWInc')/sqrt(min(size(corcorWInc)));

% Plot 2x2 Line Graph
timeAdj = 1:length(corcorWInc);
timeAdj = timeAdj./30;
timeAdj = timeAdj - 5;

figure('Units','normalized','Position',[0 0 1 1])
shadedErrorBar(timeAdj,nanmean(corcorPCon,2),semcorcorPCon,'lineprops',{'r-','LineWidth',3});
hold on
shadedErrorBar(timeAdj,nanmean(corcorPInc,2),semcorcorPInc,'lineprops',{'b-','LineWidth',3});
shadedErrorBar(timeAdj,nanmean(corcorWCon,2),semcorcorWCon,'lineprops',{'g-','LineWidth',3});
shadedErrorBar(timeAdj,nanmean(corcorWInc,2),semcorcorWInc,'lineprops',{'y-','LineWidth',3});
xlabel('Time from Cue Onset (s)')
ylabel('Likelihood to Occupy Correct Sipper')
xline(0,'k--','Cue On','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(5,'k--','Sipper In','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(13,'k--','Sipper Out','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
title('Likelihood to Occupy Correct Sipper on Correct Approach Trials')
legend([{'Congruent-P'},{'Incongruent-P'},{'Congruent-W'},{'Incongruent-W'}])
xlim([-5 20])
ylim([0 1])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,[figSavePath filesep 'sipperOccupyLH_CorrectSipper_CorrectApproach_2x2'],'svg')


%% Correct Sipper, Incorrect Approach. Keeping variable names, changing save file name.
% Draw data
corcorPCon = [trlStruct(regPIdx).llOccupy_CI];
corcorPInc = [trlStruct(revPIdx).llOccupy_CI];
corcorWCon = [trlStruct(regWIdx).llOccupy_CI];
corcorWInc = [trlStruct(revWIdx).llOccupy_CI];

% Calculate SEM
semcorcorPCon = std(corcorPCon')/sqrt(min(size(corcorPCon)));
semcorcorPInc = std(corcorPInc')/sqrt(min(size(corcorPInc)));
semcorcorWCon = std(corcorWCon')/sqrt(min(size(corcorWCon)));
semcorcorWInc = std(corcorWInc')/sqrt(min(size(corcorWInc)));

% Plot 2x2 Line Graph
timeAdj = 1:length(corcorWInc);
timeAdj = timeAdj./30;
timeAdj = timeAdj - 5;

figure('Units','normalized','Position',[0 0 1 1])
shadedErrorBar(timeAdj,nanmean(corcorPCon,2),semcorcorPCon,'lineprops',{'r-','LineWidth',3});
hold on
shadedErrorBar(timeAdj,nanmean(corcorPInc,2),semcorcorPInc,'lineprops',{'b-','LineWidth',3});
shadedErrorBar(timeAdj,nanmean(corcorWCon,2),semcorcorWCon,'lineprops',{'g-','LineWidth',3});
shadedErrorBar(timeAdj,nanmean(corcorWInc,2),semcorcorWInc,'lineprops',{'y-','LineWidth',3});
xlabel('Time from Cue Onset (s)')
ylabel('Likelihood to Occupy Correct Sipper')
xline(0,'k--','Cue On','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(5,'k--','Sipper In','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(13,'k--','Sipper Out','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
title('Likelihood to Occupy Correct Sipper on Incorrect Approach Trials')
legend([{'Congruent-P'},{'Incongruent-P'},{'Congruent-W'},{'Incongruent-W'}])
xlim([-5 20])
ylim([0 1])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,[figSavePath filesep 'sipperOccupyLH_CorrectSipper_IncorrectApproach_2x2'],'svg')

%% Correct Sipper, Incorrect Approach. Keeping variable names, changing save file name.
% Draw data
corcorPCon = [trlStruct(regPIdx).llOccupy_IC];
corcorPInc = [trlStruct(revPIdx).llOccupy_IC];
corcorWCon = [trlStruct(regWIdx).llOccupy_IC];
corcorWInc = [trlStruct(revWIdx).llOccupy_IC];

% Calculate SEM
semcorcorPCon = std(corcorPCon')/sqrt(min(size(corcorPCon)));
semcorcorPInc = std(corcorPInc')/sqrt(min(size(corcorPInc)));
semcorcorWCon = std(corcorWCon')/sqrt(min(size(corcorWCon)));
semcorcorWInc = std(corcorWInc')/sqrt(min(size(corcorWInc)));

% Plot 2x2 Line Graph
timeAdj = 1:length(corcorWInc);
timeAdj = timeAdj./30;
timeAdj = timeAdj - 5;

figure('Units','normalized','Position',[0 0 1 1])
shadedErrorBar(timeAdj,nanmean(corcorPCon,2),semcorcorPCon,'lineprops',{'r-','LineWidth',3});
hold on
shadedErrorBar(timeAdj,nanmean(corcorPInc,2),semcorcorPInc,'lineprops',{'b-','LineWidth',3});
shadedErrorBar(timeAdj,nanmean(corcorWCon,2),semcorcorWCon,'lineprops',{'g-','LineWidth',3});
shadedErrorBar(timeAdj,nanmean(corcorWInc,2),semcorcorWInc,'lineprops',{'y-','LineWidth',3});
xlabel('Time from Cue Onset (s)')
ylabel('Likelihood to Occupy Correct Sipper')
xline(0,'k--','Cue On','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(5,'k--','Sipper In','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(13,'k--','Sipper Out','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
title('Likelihood to Occupy Incorrect Sipper on Correct Approach Trials')
legend([{'Congruent-P'},{'Incongruent-P'},{'Congruent-W'},{'Incongruent-W'}])
xlim([-5 20])
ylim([0 1])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,[figSavePath filesep 'sipperOccupyLH_IncorrectSipper_CorrectApproach_2x2'],'svg')

%% Incorrect Sipper, Incorrect Approach. Keeping variable names, changing save file name.
% Draw data
corcorPCon = [trlStruct(regPIdx).llOccupy_II];
corcorPInc = [trlStruct(revPIdx).llOccupy_II];
corcorWCon = [trlStruct(regWIdx).llOccupy_II];
corcorWInc = [trlStruct(revWIdx).llOccupy_II];

% Calculate SEM
semcorcorPCon = std(corcorPCon')/sqrt(min(size(corcorPCon)));
semcorcorPInc = std(corcorPInc')/sqrt(min(size(corcorPInc)));
semcorcorWCon = std(corcorWCon')/sqrt(min(size(corcorWCon)));
semcorcorWInc = std(corcorWInc')/sqrt(min(size(corcorWInc)));

% Plot 2x2 Line Graph
timeAdj = 1:length(corcorWInc);
timeAdj = timeAdj./30;
timeAdj = timeAdj - 5;

figure('Units','normalized','Position',[0 0 1 1])
shadedErrorBar(timeAdj,nanmean(corcorPCon,2),semcorcorPCon,'lineprops',{'r-','LineWidth',3});
hold on
shadedErrorBar(timeAdj,nanmean(corcorPInc,2),semcorcorPInc,'lineprops',{'b-','LineWidth',3});
shadedErrorBar(timeAdj,nanmean(corcorWCon,2),semcorcorWCon,'lineprops',{'g-','LineWidth',3});
shadedErrorBar(timeAdj,nanmean(corcorWInc,2),semcorcorWInc,'lineprops',{'y-','LineWidth',3});
xlabel('Time from Cue Onset (s)')
ylabel('Likelihood to Occupy Correct Sipper')
xline(0,'k--','Cue On','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(5,'k--','Sipper In','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(13,'k--','Sipper Out','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
title('Likelihood to Occupy Incorrect Sipper on Incorrect Approach Trials')
legend([{'Congruent-P'},{'Incongruent-P'},{'Congruent-W'},{'Incongruent-W'}])
xlim([-5 20])
ylim([0 1])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,[figSavePath filesep 'sipperOccupyLH_IncorrectSipper_IncorrectApproach_2x2'],'svg')
%%
%% 2x2 Likelihood to Occupy Sipper (Correct and Incorrect Sippers)
trialRange = 1:15;
numTrials = 48;
numP = 14;
numW = 13;

corcorPCon = [trlStruct(regPIdx).llOccupy];
corcorPInc = [trlStruct(revPIdx).llOccupy];
corcorWCon = [trlStruct(regWIdx).llOccupy];
corcorWInc = [trlStruct(revWIdx).llOccupy];


%
corcorPCon = reshape(corcorPCon,751,numTrials,672/numTrials);
corcorPInc = reshape(corcorPInc,751,numTrials,672/numTrials);
corcorWCon = reshape(corcorWCon,751,numTrials,624/numTrials);
corcorWInc = reshape(corcorWInc,751,numTrials,624/numTrials);

corcorPCon = mean(corcorPCon,3);
corcorPInc = mean(corcorPInc,3);
corcorWCon = mean(corcorWCon,3);
corcorWInc = mean(corcorWInc,3);

corcorPCon = corcorPCon(:,trialRange);
corcorPInc = corcorPInc(:,trialRange);
corcorWCon = corcorWCon(:,trialRange);
corcorWInc = corcorWInc(:,trialRange);
% Calculate SEM
semcorcorPCon = std(corcorPCon')/sqrt(min(size(corcorPCon)));
semcorcorPInc = std(corcorPInc')/sqrt(min(size(corcorPInc)));
semcorcorWCon = std(corcorWCon')/sqrt(min(size(corcorWCon)));
semcorcorWInc = std(corcorWInc')/sqrt(min(size(corcorWInc)));

% Plot 2x2 Line Graph
timeAdj = 1:length(corcorWInc);
timeAdj = timeAdj./30;
timeAdj = timeAdj - 5;

figure('Units','normalized','Position',[0 0 1 1])
shadedErrorBar(timeAdj,nanmean(corcorPCon,2),semcorcorPCon,'lineprops',{'r-','LineWidth',3});
hold on
shadedErrorBar(timeAdj,nanmean(corcorPInc,2),semcorcorPInc,'lineprops',{'b-','LineWidth',3});
shadedErrorBar(timeAdj,nanmean(corcorWCon,2),semcorcorWCon,'lineprops',{'g-','LineWidth',3});
shadedErrorBar(timeAdj,nanmean(corcorWInc,2),semcorcorWInc,'lineprops',{'m-','LineWidth',3});
xlabel('Time from Cue Onset (s)')
ylabel('Likelihood to Occupy Correct Sipper')
xline(0,'k--','Cue On','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(5,'k--','Sipper In','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(13,'k--','Sipper Out','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
title('Likelihood to Occupy Correct Sipper')
legend([{'Congruent-P'},{'Incongruent-P'},{'Congruent-W'},{'Incongruent-W'}])
xlim([-5 20])
ylim([0 1])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,[figSavePath filesep 'sipperOccupyLH_CorrectSipper_2x2_' num2str(length(trialRange))],'svg')
%% 2x2 Likelihood to Occupy Sipper (Incorrect Sippers)
trialRange = 1:15;
numTrials = 48;
numP = 14;
numW = 13;

corcorPCon = [trlStruct(regPIdx).llOccupy_IC];
corcorPInc = [trlStruct(revPIdx).llOccupy_IC];
corcorWCon = [trlStruct(regWIdx).llOccupy_IC];
corcorWInc = [trlStruct(revWIdx).llOccupy_IC];


%
corcorPCon = reshape(corcorPCon,751,numTrials,672/numTrials);
corcorPInc = reshape(corcorPInc,751,numTrials,672/numTrials);
corcorWCon = reshape(corcorWCon,751,numTrials,624/numTrials);
corcorWInc = reshape(corcorWInc,751,numTrials,624/numTrials);

corcorPCon = mean(corcorPCon,3);
corcorPInc = mean(corcorPInc,3);
corcorWCon = mean(corcorWCon,3);
corcorWInc = mean(corcorWInc,3);

corcorPCon = corcorPCon(:,trialRange);
corcorPInc = corcorPInc(:,trialRange);
corcorWCon = corcorWCon(:,trialRange);
corcorWInc = corcorWInc(:,trialRange);
% Calculate SEM
semcorcorPCon = std(corcorPCon')/sqrt(min(size(corcorPCon)));
semcorcorPInc = std(corcorPInc')/sqrt(min(size(corcorPInc)));
semcorcorWCon = std(corcorWCon')/sqrt(min(size(corcorWCon)));
semcorcorWInc = std(corcorWInc')/sqrt(min(size(corcorWInc)));

% Plot 2x2 Line Graph
timeAdj = 1:length(corcorWInc);
timeAdj = timeAdj./30;
timeAdj = timeAdj - 5;

figure('Units','normalized','Position',[0 0 1 1])
shadedErrorBar(timeAdj,nanmean(corcorPCon,2),semcorcorPCon,'lineprops',{'r-','LineWidth',3});
hold on
shadedErrorBar(timeAdj,nanmean(corcorPInc,2),semcorcorPInc,'lineprops',{'b-','LineWidth',3});
shadedErrorBar(timeAdj,nanmean(corcorWCon,2),semcorcorWCon,'lineprops',{'g-','LineWidth',3});
shadedErrorBar(timeAdj,nanmean(corcorWInc,2),semcorcorWInc,'lineprops',{'m-','LineWidth',3});
xlabel('Time from Cue Onset (s)')
ylabel('Likelihood to Occupy Correct Sipper')
xline(0,'k--','Cue On','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(5,'k--','Sipper In','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(13,'k--','Sipper Out','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
title('Likelihood to Occupy Incorrect Sipper')
legend([{'Congruent-P'},{'Incongruent-P'},{'Congruent-W'},{'Incongruent-W'}])
xlim([-5 20])
ylim([0 1])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
saveas(gca,[figSavePath filesep 'sipperOccupyLH_InCorrectSipper_2x2_' num2str(length(trialRange))],'svg')

%% Approach and Correction Change Points
regPIdx = startsWith(masterTbl.SessionType,'Regular') & strcmp(masterTbl.Strain,'P');
revPIdx = startsWith(masterTbl.SessionType,'Reversal') & strcmp(masterTbl.Strain,'P');
    
regWIdx = startsWith(masterTbl.SessionType,'Regular') & strcmp(masterTbl.Strain,'W');
revWIdx = startsWith(masterTbl.SessionType,'Reversal') & strcmp(masterTbl.Strain,'W');

% Pulling data from structure
appCP_PCon = [trlStruct(regPIdx).changePoint];
appCP_PInc = [trlStruct(revPIdx).changePoint];
appCP_Con = [trlStruct(regWIdx).changePoint];
appCP_Inc = [trlStruct(revWIdx).changePoint];

% Calculate SEM
semappCP_PCon = std(appCP_PCon)/sqrt(length(appCP_PCon));    
semappCP_PInc = std(appCP_PInc)/sqrt(length(appCP_PInc));    
semappCP_Con = std(appCP_Con)/sqrt(length(appCP_Con));    
semappCP_Inc = std(appCP_Inc)/sqrt(length(appCP_Inc));  

% ANOVA ZONE SIPPER OCCUPANCY
testData = [appCP_PCon appCP_PInc appCP_Con appCP_Inc];
g1 = [repmat({'P rat'},1,length([appCP_PCon appCP_PInc])) repmat({'Wistar'},1,length([appCP_Con appCP_Inc]))];
g2 = [repmat({'Congruent'},1,length(appCP_PCon)) repmat({'Incongruent'},1,length(appCP_PInc)) repmat({'Congruent'},1,length(appCP_Con)) repmat({'Incongruent'},1,length(appCP_Inc))];

if length(testData) == length(g1) && length(testData) == length(g2)
    disp('Data sizes are matched')
    [p,tbl,stat] = anovan(testData,{g1,g2},'model','full');
    multcompare(stat,'Dimension',[2])
        % Write string to file
    tblStr = formattedDisplayText(tbl); 
    fid = fopen('anovaTable_2x2_changepoint.txt', 'wt');
    fileCleanup = onCleanup(@()fclose(fid));
    formatSpec = '%s\n';
    fprintf(fid, formatSpec, tblStr);
    clear('fileCleanup')
end

figure('Units','normalized','Position',[0 0 1 1])
subplot(1,2,1)
hb = bar(categorical({'P Rats','Wistars'}),[nanmean(appCP_PCon) nanmean(appCP_PInc); nanmean(appCP_Con) nanmean(appCP_Inc)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);
hb(1).FaceColor = [0.5 0.5 0.5]; hb(2).FaceColor = [0.9 0.9 0.9];
offsetPos = [1+hb(1).XOffset 1+hb(2).XOffset 2+hb(1).XOffset 2+hb(2).XOffset];
scatterXData = [linspace(offsetPos(1)-0.1,offsetPos(1)+0.1,length(appCP_PCon)),linspace(offsetPos(2)-0.1,offsetPos(2)+0.1,length(appCP_PInc)),linspace(offsetPos(3)-0.1,offsetPos(3)+0.1,length(appCP_Con)),linspace(offsetPos(4)-0.1,offsetPos(4)+0.1,length(appCP_Inc))];
hold on
plot([linspace(offsetPos(1)-0.05,offsetPos(1)+0.05,length(appCP_PCon)) linspace(offsetPos(2)-0.05,offsetPos(2)+0.05,length(appCP_PInc))]',[appCP_PCon appCP_PInc],'ko','LineWidth',2,'MarkerSize',8)
plot([linspace(offsetPos(3)-0.05,offsetPos(3)+0.05,length(appCP_Con)) linspace(offsetPos(4)-0.05,offsetPos(4)+0.05,length(appCP_Inc))]',[appCP_Con appCP_Inc],'ko','LineWidth',2,'MarkerSize',8)

er = errorbar(offsetPos,[mean(appCP_PCon) mean(appCP_PInc) mean(appCP_Con) mean(appCP_Inc)],[semappCP_PCon semappCP_PInc semappCP_Con semappCP_Inc],'LineWidth',3);    
er(1).Color = [0 0 0];                        
er(1).LineStyle = 'none'; 
legend([{'Congruent'},{'Incongruent'}])
ylabel('Session Change Point (Approaches)')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
saveas(gca,[figSavePath filesep 'approachChangePoint_2x2_allSessions'],'svg')

%%
% Pulling data from structure
appCP_PCon = [trlStruct(regPIdx).changePointCorrections];
appCP_PInc = [trlStruct(revPIdx).changePointCorrections];
appCP_Con = [trlStruct(regWIdx).changePointCorrections];
appCP_Inc = [trlStruct(revWIdx).changePointCorrections];

% Calculate SEM
semappCP_PCon = std(appCP_PCon)/sqrt(length(appCP_PCon));    
semappCP_PInc = std(appCP_PInc)/sqrt(length(appCP_PInc));    
semappCP_Con = std(appCP_Con)/sqrt(length(appCP_Con));    
semappCP_Inc = std(appCP_Inc)/sqrt(length(appCP_Inc));  

figure('Units','normalized','Position',[0 0 1 1])
hb = bar(categorical({'P Rats','Wistars'}),[mean(appCP_PCon) mean(appCP_PInc); mean(appCP_Con) mean(appCP_Inc)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);
hb(1).FaceColor = [0.5 0.5 0.5]; hb(2).FaceColor = [0.9 0.9 0.9];
offsetPos = [1+hb(1).XOffset 1+hb(2).XOffset 2+hb(1).XOffset 2+hb(2).XOffset];
scatterXData = [linspace(offsetPos(1)-0.1,offsetPos(1)+0.1,length(appCP_PCon)),linspace(offsetPos(2)-0.1,offsetPos(2)+0.1,length(appCP_PInc)),linspace(offsetPos(3)-0.1,offsetPos(3)+0.1,length(appCP_Con)),linspace(offsetPos(4)-0.1,offsetPos(4)+0.1,length(appCP_Inc))];
hold on
plot([linspace(offsetPos(1)-0.05,offsetPos(1)+0.05,length(appCP_PCon)) linspace(offsetPos(2)-0.05,offsetPos(2)+0.05,length(appCP_PInc))]',[appCP_PCon appCP_PInc],'ko','LineWidth',2,'MarkerSize',8)
plot([linspace(offsetPos(3)-0.05,offsetPos(3)+0.05,length(appCP_Con)) linspace(offsetPos(4)-0.05,offsetPos(4)+0.05,length(appCP_Inc))]',[appCP_Con appCP_Inc],'ko','LineWidth',2,'MarkerSize',8)

er = errorbar(offsetPos,[mean(appCP_PCon) mean(appCP_PInc) mean(appCP_Con) mean(appCP_Inc)],[semappCP_PCon semappCP_PInc semappCP_Con semappCP_Inc],'LineWidth',3);    
er(1).Color = [0 0 0];                        
er(1).LineStyle = 'none'; 
legend([{'Congruent'},{'Incongruent'}])
ylabel('Session Change Point (Corrections)')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
saveas(gca,[figSavePath filesep 'correctionChangePoint_2x2_allSessions'],'svg')

%% Surface Plot? 
numTrials = 48;
disThresh = 9;
disRun = 3;
sipDescent = 10*30;
sipAscent = 18*30;
cueOn = 5*30;
LeftTrials = 1:24;
RightTrials = 25:48;

numTrials = 48;

regPIdx = startsWith(masterTbl.SessionType,'Regular') & strcmp(masterTbl.Strain,'P');
revPIdx = startsWith(masterTbl.SessionType,'Reversal') & strcmp(masterTbl.Strain,'P');
    
regWIdx = startsWith(masterTbl.SessionType,'Regular') & strcmp(masterTbl.Strain,'W');
revWIdx = startsWith(masterTbl.SessionType,'Reversal') & strcmp(masterTbl.Strain,'W');

sipLL_PCon = [trlStruct(regPIdx).llOccupy];
sipLL_PInc = [trlStruct(revPIdx).llOccupy];
sipLL_WCon = [trlStruct(regWIdx).llOccupy];
sipLL_WInc = [trlStruct(revWIdx).llOccupy];

% PConIdx = [trlStruct(regPIdx).noCorrectionsIdx]; 
% PIncIdx = [trlStruct(revPIdx).noCorrectionsIdx];
% WConIdx = [trlStruct(regWIdx).noCorrectionsIdx];
% WIncIdx = [trlStruct(revWIdx).noCorrectionsIdx];

sipLL_PCon3 = reshape(sipLL_PCon,751,numTrials,672/numTrials);
sipLL_PInc3 = reshape(sipLL_PInc,751,numTrials,672/numTrials);
sipLL_WCon3 = reshape(sipLL_WCon,751,numTrials,624/numTrials);
sipLL_WInc3 = reshape(sipLL_WInc,751,numTrials,624/numTrials);

% sipLL_PCon3(:,PConIdx == 0) = NaN;
% sipLL_PInc3(:,PIncIdx == 0) = NaN;
% sipLL_WCon3(:,WConIdx == 0) = NaN;
% sipLL_WInc3(:,WIncIdx == 0) = NaN;

mean_sipLL_PCon = mean(sipLL_PCon3,3,'omitnan');
mean_sipLL_PInc = mean(sipLL_PInc3,3,'omitnan');
mean_sipLL_WCon = mean(sipLL_WCon3,3,'omitnan');
mean_sipLL_WInc = mean(sipLL_WInc3,3,'omitnan');

mean_sipLL_WDiff = (mean_sipLL_WInc - mean_sipLL_WCon);
mean_sipLL_PDiff = (mean_sipLL_PInc - mean_sipLL_PCon);

%% Plot Wistars
figure('Units','normalized','Position',[0 0 1 1])
time = [1:751]/30;
s1 = imagesc(time,1:15,(mean_sipLL_WDiff(:,1:15))',[-0.5 0.5]); colormap jet; colorbar
xline(cueOn/30,'w--','Cue On','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(sipDescent/30,'w--','Sipper In','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(sipAscent/30,'w--','Sipper Out','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xlabel('Time (s)')
ylabel('Trials')
title('Wistar Correct Sipper Probability')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
saveas(gca,[figSavePath filesep 'diff_congruent_incongruent_SipLL_perTrial_Ws_115'],'svg')
    

%% Plot p rats
figure('Units','normalized','Position',[0 0 1 1])
time = [1:751]/30;
s1 = imagesc(time,1:15,(mean_sipLL_PDiff(:,1:15))',[-0.5 0.5]); colormap jet; colorbar
xline(cueOn/30,'w--','Cue On','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(sipDescent/30,'w--','Sipper In','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(sipAscent/30,'w--','Sipper Out','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xlabel('Time (s)')
ylabel('Trials')
title('P rat Correct Sipper Probability')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
saveas(gca,[figSavePath filesep 'diff_congruent_incongruent_SipLL_perTrial_Ps_115'],'svg')

%%
%% Plot difference per epoch 
trialRange = 1:15;
sipTimeVal = abs(sipDescent - sipAscent)/2;
sipTime1 = sipDescent:sipDescent+sipTimeVal;
sipTime2 = sipDescent+sipTimeVal:sipAscent;
trialRange = 1:15;
sipTimeVal = abs(sipDescent - sipAscent)/2;
sipTime1 = sipDescent:sipDescent+sipTimeVal;
sipTime2 = sipDescent+sipTimeVal:sipAscent;


CorPreCue_PCon = nanmean(max(mean_sipLL_PDiff(1:cueOn,trialRange)));
CorCue_PCon = nanmean(max(mean_sipLL_PDiff(cueOn:sipDescent-30,trialRange)));
CorSip_PCon1 = nanmean(max(mean_sipLL_PDiff(sipTime1,trialRange)));
CorSip_PCon2 = nanmean(max(mean_sipLL_PDiff(sipTime2,trialRange)));
CorPostSip_PCon = nanmean(max(mean_sipLL_PDiff(sipAscent:end,trialRange)));


CorPreCue_WCon = nanmean(max(mean_sipLL_WDiff(1:cueOn,trialRange)));
CorCue_WCon = nanmean(max(mean_sipLL_WDiff(cueOn:sipDescent-30,trialRange)));
CorSip_WCon1 = nanmean(max(mean_sipLL_WDiff(sipTime1,trialRange)));
CorSip_WCon2 = nanmean(max(mean_sipLL_WDiff(sipTime2,trialRange)));
CorPostSip_WCon = nanmean(max(mean_sipLL_WDiff(sipAscent:end,trialRange)));


plotMat = [CorPreCue_WCon CorCue_WCon CorSip_WCon1 CorSip_WCon2 CorPostSip_WCon;... 
     CorPreCue_PCon CorCue_PCon CorSip_PCon1 CorSip_PCon2 CorPostSip_PCon];
plotMat = plotMat';

semMat = [max(mean_sipLL_WDiff(1:cueOn,trialRange)) max(mean_sipLL_WDiff(cueOn:sipDescent-30,trialRange)) max(mean_sipLL_WDiff(sipTime1,trialRange)) max(mean_sipLL_WDiff(sipTime2,trialRange)) max(mean_sipLL_WDiff(sipAscent:end,trialRange)); ...
    max(mean_sipLL_PDiff(1:cueOn,trialRange)) max(mean_sipLL_PDiff(cueOn:sipDescent-30,trialRange)) max(mean_sipLL_PDiff(sipTime1,trialRange)) max(mean_sipLL_PDiff(sipTime2,trialRange)) max(mean_sipLL_PDiff(sipAscent:end,trialRange))];

semMat = semMat';
semMat = reshape(semMat,length(trialRange),5,2);
calc_semMat = squeeze(std(semMat)/sqrt(length(trialRange)));
%%

figure('Units','normalized','Position',[0 0 1 1])

plot(plotMat(:,1),'-o','MarkerSize',7,'LineWidth',3,'Color',[0, 0.4470, 0.7410])
hold on
plot(plotMat(:,2),'--o','MarkerSize',7,'LineWidth',3,'Color',[0.8500, 0.3250, 0.0980])

% plot(plotMat(:,2:2:4),'--o','MarkerSize',10,'LineWidth',3)
er = errorbar(plotMat, calc_semMat);
for k = 1:2
    er(k).Color = [0 0 0];
    er(k).LineStyle = 'none';
    er(k).LineWidth = 3;
end
ylabel('Difference in Likelihood to Occupy Incorrect Sipper (Incongruent - Congruent)')

xlim([0.9 5.1])
xticks([1 2 3 4 5])
xticklabels({'Pre-Cue','Cue On','Early EtOH Access', 'Late EtOH Access','Post Access'})
ylim([0 0.5])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
legend([{'Wistars'},{'P Rats'}])
saveas(gca,[figSavePath filesep 'noCorrections_epochanalysis_xyz_DIFF_CORRECT' num2str(length(trialRange))],'svg')

%% Stats
tableLength = size(semMat,1)*size(semMat,3);
trialRep = repmat([1:length(trialRange)],1,size(semMat,3)); trialRep = trialRep';
tableData = [];
for i = 1:size(semMat,3)
    tableData = [tableData; squeeze(semMat(:,:,i))];
end
maxTable = array2table(tableData,'VariableNames',{'PreCue','CueOn','EarlyEtOHAccess', 'LateEtOHAccess','PostAccess'});
maxTable.Strain = [repmat({'Wistar'},length(trialRange),1); repmat({'P rat'},length(trialRange),1)];
% maxTable.Session = [repmat({'Congruent'},15,1); repmat({'Incongruent'},15,1); repmat({'Congruent'},15,1); repmat({'Incongruent'},15,1)];
maxTable.Trial = trialRep;

% Create rm model
Epoch = table([1 2 3 4 5]','VariableNames',{'Epochs'});
rm = fitrm(maxTable,'PreCue-PostAccess ~ Strain','WithinDesign',Epoch);
anovaTbl = ranova(rm,'WithinModel','Epochs')
multTbl1_Incor = multcompare(rm,'Strain','By','Epochs')
multTbl2_Incor = multcompare(rm,'Strain')

tblStr = formattedDisplayText(anovaTbl); 
% Write string to file
fid = fopen('anovaTable_correct15.txt', 'wt');
fileCleanup = onCleanup(@()fclose(fid));
formatSpec = '%s\n';
fprintf(fid, formatSpec, tblStr);
clear('fileCleanup')



%% Incorrect sipper likelihood 
numTrials = 48;
disThresh = 9;
disRun = 3;
sipDescent = 10*30;
sipAscent = 18*30;
cueOn = 5*30;
LeftTrials = 1:24;
RightTrials = 25:48;

numTrials = 48;

regPIdx = startsWith(masterTbl.SessionType,'Regular') & strcmp(masterTbl.Strain,'P');
revPIdx = startsWith(masterTbl.SessionType,'Reversal') & strcmp(masterTbl.Strain,'P');
    
regWIdx = startsWith(masterTbl.SessionType,'Regular') & strcmp(masterTbl.Strain,'W');
revWIdx = startsWith(masterTbl.SessionType,'Reversal') & strcmp(masterTbl.Strain,'W');

sipLL_PCon = [trlStruct(regPIdx).llOccupy_IC];
sipLL_PInc = [trlStruct(revPIdx).llOccupy_IC];
sipLL_WCon = [trlStruct(regWIdx).llOccupy_IC];
sipLL_WInc = [trlStruct(revWIdx).llOccupy_IC];

% PConIdx = [trlStruct(regPIdx).noCorrectionsIdx]; 
% PIncIdx = [trlStruct(revPIdx).noCorrectionsIdx];
% WConIdx = [trlStruct(regWIdx).noCorrectionsIdx];
% WIncIdx = [trlStruct(revWIdx).noCorrectionsIdx];

sipLL_PCon3 = reshape(sipLL_PCon,751,numTrials,672/numTrials);
sipLL_PInc3 = reshape(sipLL_PInc,751,numTrials,672/numTrials);
sipLL_WCon3 = reshape(sipLL_WCon,751,numTrials,624/numTrials);
sipLL_WInc3 = reshape(sipLL_WInc,751,numTrials,624/numTrials);

% sipLL_PCon3(:,PConIdx == 0) = NaN;
% sipLL_PInc3(:,PIncIdx == 0) = NaN;
% sipLL_WCon3(:,WConIdx == 0) = NaN;
% sipLL_WInc3(:,WIncIdx == 0) = NaN;

mean_sipLL_PCon = mean(sipLL_PCon3,3,'omitnan');
mean_sipLL_PInc = mean(sipLL_PInc3,3,'omitnan');
mean_sipLL_WCon = mean(sipLL_WCon3,3,'omitnan');
mean_sipLL_WInc = mean(sipLL_WInc3,3,'omitnan');

mean_sipLL_WDiff = (mean_sipLL_WInc - mean_sipLL_WCon);
mean_sipLL_PDiff = (mean_sipLL_PInc - mean_sipLL_PCon);

%% Plot testing
figure('Units','normalized','Position',[0 0 1 1])
s1 = imagesc(time,1:15,(mean_sipLL_WDiff(:,1:15))',[-0.5 0.5]); colormap jet; colorbar
xline(cueOn/30,'w--','Cue On','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(sipDescent/30,'w--','Sipper In','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(sipAscent/30,'w--','Sipper Out','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xlabel('Time (s)')
ylabel('Trials')
title('Wistars Incorrect Sipper Probability')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
saveas(gca,[figSavePath filesep 'incorrectsipper_diff_congruent_incongruent_SipLL_perTrial_Ws_115'],'svg')

%% Plot p rats
figure('Units','normalized','Position',[0 0 1 1])
s1 = imagesc(time,1:15,(mean_sipLL_PDiff(:,1:15))',[-0.5 0.5]); colormap jet; colorbar
xline(cueOn/30,'w--','Cue On','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(sipDescent/30,'w--','Sipper In','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xline(sipAscent/30,'w--','Sipper Out','LineWidth',3,'FontSize',20,'FontName','Arial','FontWeight','bold');
xlabel('Time (s)')
ylabel('Trials')
title('P rats Incorrect Sipper Probability')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
saveas(gca,[figSavePath filesep 'incorrectsipper_diff_congruent_incongruent_SipLL_perTrial_Ps_115'],'svg')


%% Plot difference per epoch 
trialRange = 1:15;
sipTimeVal = abs(sipDescent - sipAscent)/2;
sipTime1 = sipDescent:sipDescent+sipTimeVal;
sipTime2 = sipDescent+sipTimeVal:sipAscent;
trialRange = 1:15;
sipTimeVal = abs(sipDescent - sipAscent)/2;
sipTime1 = sipDescent:sipDescent+sipTimeVal;
sipTime2 = sipDescent+sipTimeVal:sipAscent;


CorPreCue_PCon = nanmean(max(mean_sipLL_PDiff(1:cueOn,trialRange)));
CorCue_PCon = nanmean(max(mean_sipLL_PDiff(cueOn:sipDescent-30,trialRange)));
CorSip_PCon1 = nanmean(max(mean_sipLL_PDiff(sipTime1,trialRange)));
CorSip_PCon2 = nanmean(max(mean_sipLL_PDiff(sipTime2,trialRange)));
CorPostSip_PCon = nanmean(max(mean_sipLL_PDiff(sipAscent:end,trialRange)));


CorPreCue_WCon = nanmean(max(mean_sipLL_WDiff(1:cueOn,trialRange)));
CorCue_WCon = nanmean(max(mean_sipLL_WDiff(cueOn:sipDescent-30,trialRange)));
CorSip_WCon1 = nanmean(max(mean_sipLL_WDiff(sipTime1,trialRange)));
CorSip_WCon2 = nanmean(max(mean_sipLL_WDiff(sipTime2,trialRange)));
CorPostSip_WCon = nanmean(max(mean_sipLL_WDiff(sipAscent:end,trialRange)));


plotMat = [CorPreCue_WCon CorCue_WCon CorSip_WCon1 CorSip_WCon2 CorPostSip_WCon;... 
     CorPreCue_PCon CorCue_PCon CorSip_PCon1 CorSip_PCon2 CorPostSip_PCon];
plotMat = plotMat';

semMat = [max(mean_sipLL_WDiff(1:cueOn,trialRange)) max(mean_sipLL_WDiff(cueOn:sipDescent-30,trialRange)) max(mean_sipLL_WDiff(sipTime1,trialRange)) max(mean_sipLL_WDiff(sipTime2,trialRange)) max(mean_sipLL_WDiff(sipAscent:end,trialRange)); ...
    max(mean_sipLL_PDiff(1:cueOn,trialRange)) max(mean_sipLL_PDiff(cueOn:sipDescent-30,trialRange)) max(mean_sipLL_PDiff(sipTime1,trialRange)) max(mean_sipLL_PDiff(sipTime2,trialRange)) max(mean_sipLL_PDiff(sipAscent:end,trialRange))];

semMat = semMat';
semMat = reshape(semMat,length(trialRange),5,2);
calc_semMat = squeeze(std(semMat)/sqrt(length(trialRange)));

%% Plot

figure('Units','normalized','Position',[0 0 1 1])

plot(plotMat(:,1),'-o','MarkerSize',7,'LineWidth',3,'Color',[0, 0.4470, 0.7410])
hold on
plot(plotMat(:,2),'--o','MarkerSize',7,'LineWidth',3,'Color',[0.8500, 0.3250, 0.0980])

% plot(plotMat(:,2:2:4),'--o','MarkerSize',10,'LineWidth',3)
er = errorbar(plotMat, calc_semMat);
for k = 1:2
    er(k).Color = [0 0 0];
    er(k).LineStyle = 'none';
    er(k).LineWidth = 3;
end
ylabel('Difference in Likelihood to Occupy Incorrect Sipper (Incongruent - Congruent)')

xlim([0.9 5.1])
xticks([1 2 3 4 5])
xticklabels({'Pre-Cue','Cue On','Early EtOH Access', 'Late EtOH Access','Post Access'})
ylim([0 0.5])
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
legend([{'Wistars'},{'P Rats'}])
saveas(gca,[figSavePath filesep 'noCorrections_epochanalysis_xyz_DIFF_INCORRECT' num2str(length(trialRange))],'svg')

%% Stats
tableLength = size(semMat,1)*size(semMat,3);
trialRep = repmat([1:length(trialRange)],1,size(semMat,3)); trialRep = trialRep';
tableData = [];
for i = 1:size(semMat,3)
    tableData = [tableData; squeeze(semMat(:,:,i))];
end
maxTable = array2table(tableData,'VariableNames',{'PreCue','CueOn','EarlyEtOHAccess', 'LateEtOHAccess','PostAccess'});
maxTable.Strain = [repmat({'Wistar'},length(trialRange),1); repmat({'P rat'},length(trialRange),1)];
% maxTable.Session = [repmat({'Congruent'},15,1); repmat({'Incongruent'},15,1); repmat({'Congruent'},15,1); repmat({'Incongruent'},15,1)];
maxTable.Trial = trialRep;

% Create rm model
Epoch = table([1 2 3 4 5]','VariableNames',{'Epochs'});
rm = fitrm(maxTable,'PreCue-PostAccess ~ Strain','WithinDesign',Epoch);
anovaTbl = ranova(rm,'WithinModel','Epochs')
multTbl1_Incor = multcompare(rm,'Strain','By','Epochs')
multTbl2_Incor = multcompare(rm,'Strain')

tblStr = formattedDisplayText(anovaTbl); 
% Write string to file
fid = fopen('anovaTable_incorrect15.txt', 'wt');
fileCleanup = onCleanup(@()fclose(fid));
formatSpec = '%s\n';
fprintf(fid, formatSpec, tblStr);
clear('fileCleanup')


%%
numTrials = 48;
disThresh = 9;
disRun = 3;
sipDescent = 10*30;
sipAscent = 18*30;
cueOn = 5*30;
LeftTrials = 1:24;
RightTrials = 25:48;

numTrials = 48;

regPIdx = startsWith(masterTbl.SessionType,'Regular') & strcmp(masterTbl.Strain,'P');
revPIdx = startsWith(masterTbl.SessionType,'Reversal') & strcmp(masterTbl.Strain,'P');
    
regWIdx = startsWith(masterTbl.SessionType,'Regular') & strcmp(masterTbl.Strain,'W');
revWIdx = startsWith(masterTbl.SessionType,'Reversal') & strcmp(masterTbl.Strain,'W');

regPApp = [trlStruct(regPIdx).CorrectApproach_NoCorrections];
revPApp = [trlStruct(revPIdx).CorrectApproach_NoCorrections];
regWApp = [trlStruct(regWIdx).CorrectApproach_NoCorrections];
revWApp = [trlStruct(revWIdx).CorrectApproach_NoCorrections];

regPSem = std(regPApp)/sqrt(length(regPApp));
revPSem = std(revPApp)/sqrt(length(revPApp));
regWSem = std(regWApp)/sqrt(length(regWApp));
revWSem = std(revWApp)/sqrt(length(revWApp));


% Plot corrections
figure('Units','normalized','Position',[0 0 1 1])
hb = bar(categorical({'P Rats','Wistars'}),[mean(regPApp) mean(revPApp); mean(regWApp) mean(revWApp)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);
hb(1).FaceColor = [0.5 0.5 0.5]; hb(2).FaceColor = [0.9 0.9 0.9];
offsetPos = [1+hb(1).XOffset 1+hb(2).XOffset 2+hb(1).XOffset 2+hb(2).XOffset];
scatterXData = [linspace(offsetPos(1)-0.1,offsetPos(1)+0.1,length(regPApp)),linspace(offsetPos(2)-0.1,offsetPos(2)+0.1,length(revPApp)),linspace(offsetPos(3)-0.1,offsetPos(3)+0.1,length(regWApp)),linspace(offsetPos(4)-0.1,offsetPos(4)+0.1,length(revWApp))];
hold on
plot([linspace(offsetPos(1)-0.05,offsetPos(1)+0.05,length(regPApp)) linspace(offsetPos(2)-0.05,offsetPos(2)+0.05,length(revPApp))]',[regPApp revPApp],'ko','LineWidth',2,'MarkerSize',8)
plot([linspace(offsetPos(3)-0.05,offsetPos(3)+0.05,length(regWApp)) linspace(offsetPos(4)-0.05,offsetPos(4)+0.05,length(revWApp))]',[regWApp revWApp],'ko','LineWidth',2,'MarkerSize',8)

er = errorbar(offsetPos,[mean(regPApp) mean(revPApp) mean(regWApp) mean(revWApp)],[regPSem revPSem regWSem revWSem],'LineWidth',3);    
er(1).Color = [0 0 0];                        
er(1).LineStyle = 'none'; 
legend([{'Congruent'},{'Incongruent'}])
ylabel('Correct Trials')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
saveas(gca,[figSavePath filesep 'correctApproach_2x2_allSessions_noCorrections'],'svg')

%%
regPIdx = startsWith(masterTbl.SessionType,'Regular') & strcmp(masterTbl.Strain,'P');
revPIdx = startsWith(masterTbl.SessionType,'Reversal') & strcmp(masterTbl.Strain,'P');
    
regWIdx = startsWith(masterTbl.SessionType,'Regular') & strcmp(masterTbl.Strain,'W');
revWIdx = startsWith(masterTbl.SessionType,'Reversal') & strcmp(masterTbl.Strain,'W');

regPNumN = [trlStruct(regPIdx).NeuronYield];
revPNumN = [trlStruct(revPIdx).NeuronYield];
regWNumN = [trlStruct(regWIdx).NeuronYield];
revWNumN = [trlStruct(revWIdx).NeuronYield];

%% CS Min Approaches
numTrials = 48;
disThresh = 9;
disRun = 3;
sipDescent = 10*30;
sipAscent = 18*30;
cueOn = 5*30;
LeftTrials = 1:24;
RightTrials = 25:48;

numTrials = 48;

regPIdx = startsWith(masterTbl.SessionType,'Regular') & strcmp(masterTbl.Strain,'P');
revPIdx = startsWith(masterTbl.SessionType,'Reversal') & strcmp(masterTbl.Strain,'P');
    
regWIdx = startsWith(masterTbl.SessionType,'Regular') & strcmp(masterTbl.Strain,'W');
revWIdx = startsWith(masterTbl.SessionType,'Reversal') & strcmp(masterTbl.Strain,'W');

regPApp = [trlStruct(regPIdx).CSRatio]; regPApp = 1 - regPApp; 
revPApp = [trlStruct(revPIdx).CSRatio]; revPApp = 1 - revPApp;
regWApp = [trlStruct(regWIdx).CSRatio]; regWApp = 1 - regWApp;
revWApp = [trlStruct(revWIdx).CSRatio]; revWApp = 1 - revWApp;

regPSem = std(regPApp)/sqrt(length(regPApp));
revPSem = std(revPApp)/sqrt(length(revPApp));
regWSem = std(regWApp)/sqrt(length(regWApp));
revWSem = std(revWApp)/sqrt(length(revWApp));

[h1,p1,~,stat1] = ttest(regPApp,0.5)
[h2,p2,~,stat2] = ttest(revPApp,0.5)
[h3,p3,~,stat3] = ttest(regWApp,0.5)
[h4,p4,~,stat4] = ttest(revWApp,0.5)

% ANOVA ZONE
testData = [regPApp revPApp regWApp revWApp];
g1 = [repmat({'P rat'},1,length([regPApp revPApp])) repmat({'Wistar'},1,length([regWApp revWApp]))];
g2 = [repmat({'Congruent'},1,length(regPApp)) repmat({'Incongruent'},1,length(revPApp)) repmat({'Congruent'},1,length(regWApp)) repmat({'Incongruent'},1,length(revWApp))];

if length(testData) == length(g1) && length(testData) == length(g2)
    disp('Data sizes are matched')
    [p,tbl,stat] = anovan(testData,{g1,g2},'model','full');
    multcompare(stat,'Dimension',[2])
        % Write string to file
    tblStr = formattedDisplayText(tbl); 
    fid = fopen('anovaTable_2x2_csminapproaches.txt', 'wt');
    fileCleanup = onCleanup(@()fclose(fid));
    formatSpec = '%s\n';
    fprintf(fid, formatSpec, tblStr);
    clear('fileCleanup')
end


% Plot corrections
figure('Units','normalized','Position',[0 0 1 1])
subplot(1,2,1)
hb = bar(categorical({'P Rats','Wistars'}),[mean(regPApp) mean(revPApp); mean(regWApp) mean(revWApp)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);
hb(1).FaceColor = [0.5 0.5 0.5]; hb(2).FaceColor = [0.9 0.9 0.9];
offsetPos = [1+hb(1).XOffset 1+hb(2).XOffset 2+hb(1).XOffset 2+hb(2).XOffset];
scatterXData = [linspace(offsetPos(1)-0.1,offsetPos(1)+0.1,length(regPApp)),linspace(offsetPos(2)-0.1,offsetPos(2)+0.1,length(revPApp)),linspace(offsetPos(3)-0.1,offsetPos(3)+0.1,length(regWApp)),linspace(offsetPos(4)-0.1,offsetPos(4)+0.1,length(revWApp))];
hold on
plot([linspace(offsetPos(1)-0.05,offsetPos(1)+0.05,length(regPApp)) linspace(offsetPos(2)-0.05,offsetPos(2)+0.05,length(revPApp))]',[regPApp revPApp],'ko','LineWidth',2,'MarkerSize',8)
plot([linspace(offsetPos(3)-0.05,offsetPos(3)+0.05,length(regWApp)) linspace(offsetPos(4)-0.05,offsetPos(4)+0.05,length(revWApp))]',[regWApp revWApp],'ko','LineWidth',2,'MarkerSize',8)
ylim([0 1])
er = errorbar(offsetPos,[mean(regPApp) mean(revPApp) mean(regWApp) mean(revWApp)],[regPSem revPSem regWSem revWSem],'LineWidth',3);    
er(1).Color = [0 0 0];                        
er(1).LineStyle = 'none'; 
legend([{'Congruent'},{'Incongruent'}])
ylabel('CS Minus Approaches')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
saveas(gca,[figSavePath filesep 'CSRatio_versusCorrect_2x2_allSession_skinny'],'svg')


%% Correction LH First 15 Trials Mean
regPApp = [trlStruct(regPIdx).correctionLH];
revPApp = [trlStruct(revPIdx).correctionLH];
regWApp = [trlStruct(regWIdx).correctionLH];
revWApp = [trlStruct(revWIdx).correctionLH];

regPApp = movmean(regPApp,3);
revPApp = movmean(revPApp,3);
regWApp = movmean(regWApp,3);
revWApp = movmean(revWApp,3);

regPApp = mean(regPApp(1:15,:));
revPApp = mean(revPApp(1:15,:));
regWApp = mean(regWApp(1:15,:));
revWApp = mean(revWApp(1:15,:));

regPSem = std(regPApp')/sqrt(length(regPApp));
revPSem = std(revPApp')/sqrt(length(revPApp));
regWSem = std(regWApp')/sqrt(length(regWApp));
revWSem = std(revWApp')/sqrt(length(revWApp));

testData = [regPApp revPApp regWApp revWApp];
g1 = [repmat({'P rat'},1,length([regPApp revPApp])) repmat({'Wistar'},1,length([regWApp revWApp]))];
g2 = [repmat({'Congruent'},1,length(regPApp)) repmat({'Incongruent'},1,length(revPApp)) repmat({'Congruent'},1,length(regWApp)) repmat({'Incongruent'},1,length(revWApp))];

if length(testData) == length(g1) && length(testData) == length(g2)
    disp('Data sizes are matched')
    [p,tbl,stat] = anovan(testData,{g1,g2},'model','full');
    multcompare(stat,'Dimension',[1 2])
        % Write string to file
    tblStr = formattedDisplayText(tbl); 
    fid = fopen('anovaTable_2x2_corrections_1_15.txt', 'wt');
    fileCleanup = onCleanup(@()fclose(fid));
    formatSpec = '%s\n';
    fprintf(fid, formatSpec, tblStr);
    clear('fileCleanup')
end


figure('Units','normalized','Position',[0 0 1 1])
subplot(1,2,1)
hb = bar(categorical({'P Rats','Wistars'}),[mean(regPApp) mean(revPApp); mean(regWApp) mean(revWApp)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);
hb(1).FaceColor = [0.5 0.5 0.5]; hb(2).FaceColor = [0.9 0.9 0.9];
offsetPos = [1+hb(1).XOffset 1+hb(2).XOffset 2+hb(1).XOffset 2+hb(2).XOffset];
scatterXData = [linspace(offsetPos(1)-0.1,offsetPos(1)+0.1,length(regPApp)),linspace(offsetPos(2)-0.1,offsetPos(2)+0.1,length(revPApp)),linspace(offsetPos(3)-0.1,offsetPos(3)+0.1,length(regWApp)),linspace(offsetPos(4)-0.1,offsetPos(4)+0.1,length(revWApp))];
hold on
plot([linspace(offsetPos(1)-0.05,offsetPos(1)+0.05,length(regPApp)) linspace(offsetPos(2)-0.05,offsetPos(2)+0.05,length(revPApp))]',[regPApp revPApp],'ko','LineWidth',2,'MarkerSize',8)
plot([linspace(offsetPos(3)-0.05,offsetPos(3)+0.05,length(regWApp)) linspace(offsetPos(4)-0.05,offsetPos(4)+0.05,length(revWApp))]',[regWApp revWApp],'ko','LineWidth',2,'MarkerSize',8)
ylim([0 1])
er = errorbar(offsetPos,[mean(regPApp) mean(revPApp) mean(regWApp) mean(revWApp)],[regPSem revPSem regWSem revWSem],'LineWidth',3);    
er(1).Color = [0 0 0];                        
er(1).LineStyle = 'none'; 
legend([{'Congruent'},{'Incongruent'}])
ylabel('Likelihood to Correct (First 15 Trials)')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
saveas(gca,[figSavePath filesep 'LHCorrect_1_15__2x2_allSession_skinny'],'svg')



