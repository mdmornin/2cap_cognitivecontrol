function [trlStruct,sessionStatistics] = getSessionStats(trlStruct,figSavePath,masterTbl,flags)
%% 2 Components: Generate data, plot data.
% Generation of data will be hard coded based on flags utilized, plotting will be
% flexible to flags

%% Find some summary information on Reg and Rev Trials (All Sessions and Genotypes)
% Determine Correct Approaches per Session Type
% Currently only capable of handling genotype differences
regCorrect = []; revCorrect = [];
if strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'all')    
    for i = 1:length(trlStruct)
        if startsWith(masterTbl.SessionType{i},['Regular']) 
            regCorrect = [regCorrect sum([min(size(trlStruct(i).LDistCorrectAppr)) min(size(trlStruct(i).RDistCorrectAppr))])];
            trlStruct(i).Correct = sum([min(size(trlStruct(i).LDistCorrectAppr)) min(size(trlStruct(i).RDistCorrectAppr))]);
        elseif startsWith(masterTbl.SessionType{i},['Reversal']) 
            revCorrect = [revCorrect sum([min(size(trlStruct(i).LDistCorrectAppr)) min(size(trlStruct(i).RDistCorrectAppr))])];
            trlStruct(i).Correct = sum([min(size(trlStruct(i).LDistCorrectAppr)) min(size(trlStruct(i).RDistCorrectAppr))]);
        end
    end
elseif strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'P') 
    for i = 1:length(trlStruct)
        if startsWith(masterTbl.SessionType{i},['Regular']) && strcmp(masterTbl.Strain{i},'P')
            regCorrect = [regCorrect sum([min(size(trlStruct(i).LDistCorrectAppr)) min(size(trlStruct(i).RDistCorrectAppr))])];
        elseif startsWith(masterTbl.SessionType{i},['Reversal']) && strcmp(masterTbl.Strain{i},'P')
            revCorrect = [revCorrect sum([min(size(trlStruct(i).LDistCorrectAppr)) min(size(trlStruct(i).RDistCorrectAppr))])];
        end
    end
elseif strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'W') 
    for i = 1:length(trlStruct)
        if startsWith(masterTbl.SessionType{i},['Regular']) && strcmp(masterTbl.Strain{i},'W')
            regCorrect = [regCorrect sum([min(size(trlStruct(i).LDistCorrectAppr)) min(size(trlStruct(i).RDistCorrectAppr))])];
        elseif startsWith(masterTbl.SessionType{i},['Reversal']) && strcmp(masterTbl.Strain{i},'W')
            revCorrect = [revCorrect sum([min(size(trlStruct(i).LDistCorrectAppr)) min(size(trlStruct(i).RDistCorrectAppr))])];
        end
    end
end

%% Create Figures for Correct Approaches
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
saveas(gca,[figSavePath filesep 'correctApproachNum_Session_' flags.SessionN '_Strain_' flags.Genotype],'svg')
    
%% Determine incorrect Approaches per Session Type
    regInCorrect = []; revInCorrect = [];

    if strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'all')    
        for i = 1:length(trlStruct)
            if startsWith(masterTbl.SessionType{i},['Regular']) 
                regInCorrect = [regInCorrect sum([min(size(trlStruct(i).LDistInCorrectAppr)) min(size(trlStruct(i).RDistInCorrectAppr))])];
                trlStruct(i).Incorrect = sum([min(size(trlStruct(i).LDistInCorrectAppr)) min(size(trlStruct(i).RDistInCorrectAppr))]);
            elseif startsWith(masterTbl.SessionType{i},['Reversal']) 
                revInCorrect = [revInCorrect sum([min(size(trlStruct(i).LDistInCorrectAppr)) min(size(trlStruct(i).RDistInCorrectAppr))])];
                trlStruct(i).Incorrect = sum([min(size(trlStruct(i).LDistInCorrectAppr)) min(size(trlStruct(i).RDistInCorrectAppr))]);
            end
        end
    elseif strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'P') 
        for i = 1:length(trlStruct)
            if startsWith(masterTbl.SessionType{i},['Regular']) && strcmp(masterTbl.Strain{i},'P')
                regInCorrect = [regInCorrect sum([min(size(trlStruct(i).LDistInCorrectAppr)) min(size(trlStruct(i).RDistInCorrectAppr))])];
            elseif startsWith(masterTbl.SessionType{i},['Reversal']) && strcmp(masterTbl.Strain{i},'P')
                revInCorrect = [revInCorrect sum([min(size(trlStruct(i).LDistInCorrectAppr)) min(size(trlStruct(i).RDistInCorrectAppr))])];
            end
        end
    elseif strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'W') 
        for i = 1:length(trlStruct)
            if startsWith(masterTbl.SessionType{i},['Regular']) && strcmp(masterTbl.Strain{i},'W')
                regInCorrect = [regInCorrect sum([min(size(trlStruct(i).LDistInCorrectAppr)) min(size(trlStruct(i).RDistInCorrectAppr))])];
            elseif startsWith(masterTbl.SessionType{i},['Reversal']) && strcmp(masterTbl.Strain{i},'W')
                revInCorrect = [revInCorrect sum([min(size(trlStruct(i).LDistInCorrectAppr)) min(size(trlStruct(i).RDistInCorrectAppr))])];
            end
        end
    end

 %% Create figure for incorrect approaches    
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
    saveas(gca,[figSavePath filesep 'incorrectApproachNum_Session_' flags.SessionN '_Strain_' flags.Genotype],'svg')
    
%% Determine omissions per session type
% Init variables
regOmit = []; revOmit = [];

% Pull data
if strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'all')  

    for i = 1:length(trlStruct)

        if startsWith(masterTbl.SessionType{i},['Regular'])
            regOmit = [regOmit sum(sum([trlStruct(i).approach(1:48,1:2,1) == 1],2)==0)];
            trlStruct(i).Omissions = sum(sum([trlStruct(i).approach(1:48,1:2,1) == 1],2)==0);
        elseif startsWith(masterTbl.SessionType{i},['Reversal'])
            revOmit = [revOmit sum(sum([trlStruct(i).approach(1:48,1:2,1) == 1],2)==0)];
            trlStruct(i).Omissions = sum(sum([trlStruct(i).approach(1:48,1:2,1) == 1],2)==0);
        end
    
    end

elseif strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'P')

    for i = 1:length(trlStruct)

        if startsWith(masterTbl.SessionType{i},['Regular']) && strcmp(masterTbl.Strain{i},'P')
            regOmit = [regOmit sum(sum([trlStruct(i).approach(1:48,1:2,1) == 1],2)==0)];
        elseif startsWith(masterTbl.SessionType{i},['Reversal']) && strcmp(masterTbl.Strain{i},'P')
            revOmit = [revOmit sum(sum([trlStruct(i).approach(1:48,1:2,1) == 1],2)==0)];
        end
    
    end

elseif strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'W')

    for i = 1:length(trlStruct)

        if startsWith(masterTbl.SessionType{i},['Regular']) && strcmp(masterTbl.Strain{i},'W')
            regOmit = [regOmit sum(sum([trlStruct(i).approach(1:48,1:2,1) == 1],2)==0)];
        elseif startsWith(masterTbl.SessionType{i},['Reversal']) && strcmp(masterTbl.Strain{i},'W')
            revOmit = [revOmit sum(sum([trlStruct(i).approach(1:48,1:2,1) == 1],2)==0)];
        end
        
    end

end
% Calculate SEM for plots
semregOmit = std(regOmit)/sqrt(length(regOmit));
semrevcorLat = std(revOmit)/sqrt(length(revOmit));


% Plotting!
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
% Saving!
saveas(gca,[figSavePath filesep 'ommissionNum_Session_' flags.SessionN '_Strain_' flags.Genotype],'svg')

%% Determine corrections per session type
regCorrections = []; regIncorrect = [];
revCorrections = []; revIncorrect = [];

% Pull data
if strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'all')  

    for i = 1:length(trlStruct)
        vecCorrections = [];
        if startsWith(masterTbl.SessionType{i},['Regular'])
            vecCorrections = trlStruct(i).approach(1:48,1:2,1);
            regCorrections = [regCorrections nansum(vecCorrections(:,1) == 1 & vecCorrections(:,2) == 1)];
            regIncorrect = [regIncorrect nansum(vecCorrections(:,2))];
            trlStruct(i).Corrections = nansum(vecCorrections(:,1) == 1 & vecCorrections(:,2) == 1)/nansum(vecCorrections(:,2));
            trlStruct(i).correctionRaw = nansum(vecCorrections(:,1) == 1 & vecCorrections(:,2) == 1);
        elseif startsWith(masterTbl.SessionType{i},['Reversal'])
            vecCorrections = trlStruct(i).approach(1:48,1:2,1);
            revCorrections = [revCorrections nansum(vecCorrections(:,1) == 1 & vecCorrections(:,2) == 1)];
            revIncorrect = [revIncorrect nansum(vecCorrections(:,1))];
            trlStruct(i).Corrections = nansum(vecCorrections(:,1) == 1 & vecCorrections(:,2) == 1)/nansum(vecCorrections(:,1));
            trlStruct(i).correctionRaw = nansum(vecCorrections(:,1) == 1 & vecCorrections(:,2) == 1);
    
        end
    end

elseif strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'P')

    for i = 1:length(trlStruct)
        vecCorrections = [];
        if startsWith(masterTbl.SessionType{i},['Regular']) && strcmp(masterTbl.Strain{i},'P')
            vecCorrections = trlStruct(i).approach(1:48,1:2,1);
            regCorrections = [regCorrections nansum(vecCorrections(:,1) == 1 & vecCorrections(:,2) == 1)];
            regIncorrect = [regIncorrect nansum(vecCorrections(:,2))];
            trlStruct(i).correctionRaw = nansum(vecCorrections(:,1) == 1 & vecCorrections(:,2) == 1);
        elseif startsWith(masterTbl.SessionType{i},['Reversal']) && strcmp(masterTbl.Strain{i},'P')
            vecCorrections = trlStruct(i).approach(1:48,1:2,1);
            revCorrections = [revCorrections nansum(vecCorrections(:,1) == 1 & vecCorrections(:,2) == 1)];
            revIncorrect = [revIncorrect nansum(vecCorrections(:,1))];
            trlStruct(i).correctionRaw = nansum(vecCorrections(:,1) == 1 & vecCorrections(:,2) == 1);
    
        end
    end

elseif strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'W')

    for i = 1:length(trlStruct)
        vecCorrections = [];
        if startsWith(masterTbl.SessionType{i},['Regular']) && strcmp(masterTbl.Strain{i},'W')
            vecCorrections = trlStruct(i).approach(1:48,1:2,1);
            regCorrections = [regCorrections nansum(vecCorrections(:,1) == 1 & vecCorrections(:,2) == 1)];
            regIncorrect = [regIncorrect nansum(vecCorrections(:,2))];
            trlStruct(i).correctionRaw = nansum(vecCorrections(:,1) == 1 & vecCorrections(:,2) == 1);
        elseif startsWith(masterTbl.SessionType{i},['Reversal']) && strcmp(masterTbl.Strain{i},'W')
            vecCorrections = trlStruct(i).approach(1:48,1:2,1);
            revCorrections = [revCorrections nansum(vecCorrections(:,1) == 1 & vecCorrections(:,2) == 1)];
            revIncorrect = [revIncorrect nansum(vecCorrections(:,1))];
            trlStruct(i).correctionRaw = nansum(vecCorrections(:,1) == 1 & vecCorrections(:,2) == 1);
    
        end
    end

end
% Plot Corrections, Raw Value
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

set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
saveas(gca,[figSavePath filesep 'corrections_Session_' flags.SessionN '_Strain_' flags.Genotype],'svg')

% Plot number of Corrections Standardized to Total Incorrect Approaches
regStCorrections = regCorrections./regIncorrect;
revStCorrections = revCorrections./revIncorrect;


semregStCorrections = std(regStCorrections)/sqrt(length(regStCorrections));
semrevStCorrections = std(revStCorrections)/sqrt(length(revStCorrections));

figure('Units','normalized','Position',[0 0 1 1])
hb = bar(categorical({'Congruent Sessions','Incongruent Sessions'}),[mean(regStCorrections); mean(revStCorrections)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);
hb.FaceColor = 'flat';
offsetPos = [1+hb(1).XOffset; 2+hb(1).XOffset];
colorData = hb.CData;
hb.CData(1,:) = [0.5 0.5 0.5]; hb.CData(2,:) = [0.9 0.9 0.9];

hold on

plot([linspace(0.95,1.05,length(regStCorrections));linspace(1.95,2.05,length(revStCorrections))],[regStCorrections; revStCorrections],'k--o','LineWidth',3,'MarkerSize',10)
er = errorbar(offsetPos,[mean(regStCorrections); mean(revStCorrections)],[semregStCorrections; semrevStCorrections],'LineWidth',3);    
er.Color = [0 0 0];                           
er.LineStyle = 'none';
ylim([-0.01 1.1])
ylabel('Proportion of Corrections')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
saveas(gca,[figSavePath filesep 'correctionsProportion_Session_' flags.SessionN '_Strain_' flags.Genotype],'svg')
%% Find # of Correct Approaches without Corrections
regCorrect = []; revCorrect = [];
if strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'all')    
    for i = 1:length(trlStruct)
        if startsWith(masterTbl.SessionType{i},['Regular']) 
            regVec = trlStruct(i).approach(1:48,1,1) == 1 & trlStruct(i).approach(1:48,2,1) == 0;
            trlStruct(i).CorrectNoCorrections = sum(regVec);
            regCorrect = [regCorrect sum(regVec)];

        elseif startsWith(masterTbl.SessionType{i},['Reversal']) 
            revVec = trlStruct(i).approach(1:48,1,1) == 0 & trlStruct(i).approach(1:48,2,1) == 1;
            trlStruct(i).CorrectNoCorrections = sum(revVec);
            revCorrect = [revCorrect sum(revVec)];

        end
    end
elseif strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'P') 
    for i = 1:length(trlStruct)
        if startsWith(masterTbl.SessionType{i},['Regular']) && strcmp(masterTbl.Strain{i},'P')
            regVec = trlStruct(i).approach(1:48,1,1) == 1 & trlStruct(i).approach(1:48,2,1) == 0;
            trlStruct(i).CorrectNoCorrections = sum(regVec);
            regCorrect = [regCorrect sum(regVec)];

        elseif startsWith(masterTbl.SessionType{i},['Reversal']) && strcmp(masterTbl.Strain{i},'P')
            revVec = trlStruct(i).approach(1:48,1,1) == 0 & trlStruct(i).approach(1:48,2,1) == 1;
            trlStruct(i).CorrectNoCorrections = sum(revVec);    
            revCorrect = [revCorrect sum(revVec)];

        end
    end
elseif strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'W') 
    for i = 1:length(trlStruct)
        if startsWith(masterTbl.SessionType{i},['Regular']) && strcmp(masterTbl.Strain{i},'W')
            regVec= trlStruct(i).approach(1:48,1,1) == 1 & trlStruct(i).approach(1:48,2,1) == 0;
            trlStruct(i).CorrectNoCorrections = sum(regVec);
            regCorrect = [regCorrect sum(regVec)];

        elseif startsWith(masterTbl.SessionType{i},['Reversal']) && strcmp(masterTbl.Strain{i},'W')
            revVec = trlStruct(i).approach(1:48,1,1) == 0 & trlStruct(i).approach(1:48,2,1) == 1;
            trlStruct(i).CorrectNoCorrections = sum(revVec);
            revCorrect = [revCorrect sum(revVec)];

        end
    end
end

% Plot
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
saveas(gca,[figSavePath filesep 'NOCORRECTIONS_correctApproachNum_Session_' flags.SessionN '_Strain_' flags.Genotype],'svg')

%% Output statistics structure
sessionStatistics.regCorrect = regCorrect;
sessionStatistics.revCorrect = revCorrect;

sessionStatistics.regInCorrect = regInCorrect;
sessionStatistics.revInCorrect = revInCorrect;

sessionStatistics.regOmission = regOmit;
sessionStatistics.revOmission = revOmit;

sessionStatistics.regCorrection = regCorrections; % Put raw values in since the proprotion of corrections is readily and easily derivable 
sessionStatistics.revCorrection = revCorrections;
end
