
function [trlStruct,latencyStats] = getLatency(trlStruct,figSavePath,masterTbl,flags)
%% Determine secondary statistics such as latency
% Find first moment following sipper descent where animal makes contact
% with sipper pixels. Use threshold set by Nick's code (< 9 pixels) and
% threshold for timebins (> 3 timebins). First start by indexing all
% moments where distance to correct sipper is less than 9 pixels, then
% ensure that the # of pixels is greater than 3. The first index where this
% occurs will be the latency to approach the sipper. 

%% Params
LeftTrials = 1:24;
RightTrials = 25:48;
disThresh = 9;
disRun = 3;
sipDescent = 10*30;
sipAscent = 18*30;
cueOn = 5*30;
ReglatCor = []; ReglatinCor = []; RevlatCor = []; RevlatinCor = []; % Correct and Incorrect latencies
%% Extract data according to flags, 

if strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'all')    

    for i = 1:length(trlStruct)
        if startsWith(masterTbl.SessionType{i},['Regular'])
            latCorIdx = [trlStruct(i).RDistCorrectAppr <= disThresh; trlStruct(i).LDistCorrectAppr <= disThresh];
            latInCorIdx = [trlStruct(i).trlRSipDist(LeftTrials(trlStruct(i).LincorrectIdx),:) <= disThresh; trlStruct(i).trlLSipDist(RightTrials(trlStruct(i).RincorrectIdx),:) <= disThresh];
            [~,sessCorLats] = max(latCorIdx(:,sipDescent:sipAscent),[],2);
            [~,sessInCorLats] = max(latInCorIdx(:,sipDescent:sipAscent),[],2);
            trlStruct(i).CorrectLatencies = sessCorLats'./30;
            trlStruct(i).IncorrectLatencies = sessInCorLats'./30;
            excludeCorrections = trlStruct(i).approach(1:48,1,1) == 1 & trlStruct(i).approach(1:48,2,1) == 0;
            ReglatCor = [ReglatCor; sessCorLats/30;];
            ReglatinCor = [ReglatinCor; sessInCorLats/30;];
        elseif startsWith(masterTbl.SessionType{i},['Reversal'])
            latCorIdx = [trlStruct(i).RDistCorrectAppr <= disThresh; trlStruct(i).LDistCorrectAppr <= disThresh];
            latInCorIdx = [trlStruct(i).trlRSipDist(RightTrials(trlStruct(i).RincorrectIdx),:) <= disThresh; trlStruct(i).trlLSipDist(LeftTrials(trlStruct(i).LincorrectIdx),:) <= disThresh];
            [~,sessCorLats] = max(latCorIdx(:,sipDescent:sipAscent),[],2);
            [~,sessInCorLats] = max(latInCorIdx(:,sipDescent:sipAscent),[],2);
            trlStruct(i).CorrectLatencies = sessCorLats'./30;
            trlStruct(i).IncorrectLatencies = sessInCorLats'./30;
            excludeCorrections = trlStruct(i).approach(1:48,1,1) == 0 & trlStruct(i).approach(1:48,2,1) == 1;
            RevlatCor = [RevlatCor; sessCorLats/30;];
            RevlatinCor = [RevlatinCor; sessInCorLats/30;];
        end
    end

elseif strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'P')    

    for i = 1:length(trlStruct)
        if startsWith(masterTbl.SessionType{i},['Regular'])  && strcmp(masterTbl.Strain{i},'P') 
            latCorIdx = [trlStruct(i).RDistCorrectAppr <= disThresh; trlStruct(i).LDistCorrectAppr <= disThresh];
            latInCorIdx = [trlStruct(i).trlRSipDist(LeftTrials(trlStruct(i).LincorrectIdx),:) <= disThresh; trlStruct(i).trlLSipDist(RightTrials(trlStruct(i).RincorrectIdx),:) <= disThresh];
            [~,sessCorLats] = max(latCorIdx(:,cueOn:sipAscent),[],2);
            [~,sessInCorLats] = max(latInCorIdx(:,cueOn:sipAscent),[],2);
            trlStruct(i).CorrectLatencies = sessCorLats./30;
            trlStruct(i).IncorrectLatencies = sesInCorLats./30;
            ReglatCor = [ReglatCor; sessCorLats/30;];
            ReglatinCor = [ReglatinCor; sessInCorLats/30;];
        elseif startsWith(masterTbl.SessionType{i},['Reversal']) && strcmp(masterTbl.Strain{i},'P') 
            latCorIdx = [trlStruct(i).RDistCorrectAppr <= disThresh; trlStruct(i).LDistCorrectAppr <= disThresh];
            latInCorIdx = [trlStruct(i).trlRSipDist(RightTrials(trlStruct(i).RincorrectIdx),:) <= disThresh; trlStruct(i).trlLSipDist(LeftTrials(trlStruct(i).LincorrectIdx),:) <= disThresh];
            [~,sessCorLats] = max(latCorIdx(:,cueOn:sipAscent),[],2);
            [~,sessInCorLats] = max(latInCorIdx(:,cueOn:sipAscent),[],2);
            trlStruct(i).CorrectLatencies = sessCorLats./30;
            trlStruct(i).IncorrectLatencies = sesInCorLats./30;
            RevlatCor = [RevlatCor; sessCorLats/30;];
            RevlatinCor = [RevlatinCor; sessInCorLats/30;];
        end
    end

elseif strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'W')    

    for i = 1:length(trlStruct)
        if startsWith(masterTbl.SessionType{i},['Regular'])  && strcmp(masterTbl.Strain{i},'W') 
            latCorIdx = [trlStruct(i).RDistCorrectAppr <= disThresh; trlStruct(i).LDistCorrectAppr <= disThresh];
            latInCorIdx = [trlStruct(i).trlRSipDist(LeftTrials(trlStruct(i).LincorrectIdx),:) <= disThresh; trlStruct(i).trlLSipDist(RightTrials(trlStruct(i).RincorrectIdx),:) <= disThresh];
            [~,sessCorLats] = max(latCorIdx(:,cueOn:sipAscent),[],2);
            [~,sessInCorLats] = max(latInCorIdx(:,cueOn:sipAscent),[],2);
            trlStruct(i).CorrectLatencies = sessCorLats./30;
            trlStruct(i).IncorrectLatencies = sesInCorLats./30;
            ReglatCor = [ReglatCor; sessCorLats/30;];
            ReglatinCor = [ReglatinCor; sessInCorLats/30;];
        elseif startsWith(masterTbl.SessionType{i},['Reversal']) && strcmp(masterTbl.Strain{i},'W') 
            latCorIdx = [trlStruct(i).RDistCorrectAppr <= disThresh; trlStruct(i).LDistCorrectAppr <= disThresh];
            latInCorIdx = [trlStruct(i).trlRSipDist(RightTrials(trlStruct(i).RincorrectIdx),:) <= disThresh; trlStruct(i).trlLSipDist(LeftTrials(trlStruct(i).LincorrectIdx),:) <= disThresh];
            [~,sessCorLats] = max(latCorIdx(:,cueOn:sipAscent),[],2);
            [~,sessInCorLats] = max(latInCorIdx(:,cueOn:sipAscent),[],2);
            trlStruct(i).CorrectLatencies = sessCorLats./30;
            trlStruct(i).IncorrectLatencies = sesInCorLats./30;
            RevlatCor = [RevlatCor; sessCorLats/30;];
            RevlatinCor = [RevlatinCor; sessInCorLats/30;];
        end
    end

end

%% Create latency figure
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
saveas(gca,[figSavePath filesep 'Latencies_Sessions_' flags.SessionN '_Strain_' flags.Genotype],'png')

%% Output structure
latencyStats.regCorrect = ReglatCor;
latencyStats.regIncorrect = ReglatinCor;
latencyStats.revCorrect = RevlatCor;
latencyStats.revIncorrect = RevlatinCor;
