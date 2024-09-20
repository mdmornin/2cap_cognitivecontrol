function [approachStatistics,trlStruct] = getApproachLH(trlStruct,masterTbl,figSavePath,flags)

% Global params

    % Apporach Likelihood. Correct or Incorrect Approaches (Not combined).
    % First gather data, then plot

RegcorApproachVector = []; RevcorApproachVector = [];
RegincApproachVector = []; RevincApproachVector = [];

if strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'all')    

    for i = 1:length(trlStruct)
        if startsWith(masterTbl.SessionType{i},['Regular'])
            currentTrialTimes = trlStruct(i).trialTimes;
            [~,sortedIndex] = sort(currentTrialTimes(1:48));
            currentApproachVector = trlStruct(i).approach(1:48,:,1);
            sortedCorrectApproach = currentApproachVector(sortedIndex,1);
            sortedIncrectApproach = currentApproachVector(sortedIndex,2);
            RegcorApproachVector = [RegcorApproachVector sortedCorrectApproach];
            RegincApproachVector = [RegincApproachVector sortedIncrectApproach];
            trlStruct(i).changePoint = findchangepts(movmean(sortedCorrectApproach > 0 | sortedIncrectApproach > 0,3));

        elseif startsWith(masterTbl.SessionType{i},['Reversal'])
            currentTrialTimes = trlStruct(i).trialTimes;
            [~,sortedIndex] = sort(currentTrialTimes(1:48));
            currentApproachVector = trlStruct(i).approach(1:48,:,1);
            sortedCorrectApproach = currentApproachVector(sortedIndex,2);
            sortedIncrectApproach = currentApproachVector(sortedIndex,1);
            RevcorApproachVector = [RevcorApproachVector sortedCorrectApproach];
            RevincApproachVector = [RevincApproachVector sortedIncrectApproach];
            trlStruct(i).changePoint = findchangepts(movmean(sortedCorrectApproach > 0 | sortedIncrectApproach > 0,3));

        end
    end

elseif strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'P')

    for i = 1:length(trlStruct)
        if startsWith(masterTbl.SessionType{i},['Regular'])  && strcmp(masterTbl.Strain{i},'P')
            currentTrialTimes = trlStruct(i).trialTimes;
            [~,sortedIndex] = sort(currentTrialTimes(1:48));
            currentApproachVector = trlStruct(i).approach(1:48,:,1);
            sortedCorrectApproach = currentApproachVector(sortedIndex,1);
            sortedIncrectApproach = currentApproachVector(sortedIndex,2);
            RegcorApproachVector = [RegcorApproachVector sortedCorrectApproach];
            RegincApproachVector = [RegincApproachVector sortedIncrectApproach];
            trlStruct(i).changePoint = findchangepts(movmean(sortedCorrectApproach > 0 | sortedIncrectApproach > 0,3));
        elseif startsWith(masterTbl.SessionType{i},['Reversal'])  && strcmp(masterTbl.Strain{i},'P')
            currentTrialTimes = trlStruct(i).trialTimes;
            [~,sortedIndex] = sort(currentTrialTimes(1:48));
            currentApproachVector = trlStruct(i).approach(1:48,:,1);
            sortedCorrectApproach = currentApproachVector(sortedIndex,2);
            sortedIncrectApproach = currentApproachVector(sortedIndex,1);
            RevcorApproachVector = [RevcorApproachVector sortedCorrectApproach];
            RevincApproachVector = [RevincApproachVector sortedIncrectApproach];
            trlStruct(i).changePoint = findchangepts(movmean(sortedCorrectApproach > 0 | sortedIncrectApproach > 0,3));
        end
    end

elseif strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'W')

    for i = 1:length(trlStruct)
        if startsWith(masterTbl.SessionType{i},['Regular'])  && strcmp(masterTbl.Strain{i},'W') 
            currentTrialTimes = trlStruct(i).trialTimes;
            [~,sortedIndex] = sort(currentTrialTimes(1:48));
            currentApproachVector = trlStruct(i).approach(1:48,:,1);
            sortedCorrectApproach = currentApproachVector(sortedIndex,1);
            sortedIncrectApproach = currentApproachVector(sortedIndex,2);
            RegcorApproachVector = [RegcorApproachVector sortedCorrectApproach];
            RegincApproachVector = [RegincApproachVector sortedIncrectApproach];
            trlStruct(i).changePoint = findchangepts(movmean(sortedCorrectApproach > 0 | sortedIncrectApproach > 0,3));
        elseif startsWith(masterTbl.SessionType{i},['Reversal'])  && strcmp(masterTbl.Strain{i},'W')
            currentTrialTimes = trlStruct(i).trialTimes;
            [~,sortedIndex] = sort(currentTrialTimes(1:48));
            currentApproachVector = trlStruct(i).approach(1:48,:,1);
            sortedCorrectApproach = currentApproachVector(sortedIndex,2);
            sortedIncrectApproach = currentApproachVector(sortedIndex,1);
            RevcorApproachVector = [RevcorApproachVector sortedCorrectApproach];
            RevincApproachVector = [RevincApproachVector sortedIncrectApproach];
            trlStruct(i).changePoint = findchangepts(movmean(sortedCorrectApproach > 0 | sortedIncrectApproach > 0,3));
        end
    end

end


% Plotting and Pulling Data

if flags.combined == 0
    % Approach likelihood. Correct OR incorrect approaches, not combined. 
    RegCorM = movmean(RegcorApproachVector,3); RegincorM = movmean(RegincApproachVector,3);
    RevCorM = movmean(RevcorApproachVector,3); RevincorM = movmean(RevincApproachVector,3);
    RegCorMSEM = std(RegCorM')/sqrt(size(RegCorM,2)); RegincorSEM = std(RegincorM')/sqrt(size(RegincorM,2));
    RevCorMSEM = std(RevCorM')/sqrt(size(RevCorM,2)); RevincorSEM = std(RevincorM')/sqrt(size(RevincorM,2));
    % Plotting data, 2 figures
    figure('Units','normalized','Position',[0 0 1 1])
    plot(mean(RegCorM,2),'r-o','LineWidth',3); hold on; plot(mean(RevCorM,2),'b-o','LineWidth',3);
    er = errorbar([mean(RegCorM,2) mean(RevCorM,2)],[RegCorMSEM' RevCorMSEM'],'LineWidth',3); 
    er(1).Color = 'r'; er(2).Color = 'b';                          
    er(1).LineStyle = 'none'; er(2).LineStyle = 'none'; 
    ylim([0 1.1])
    xlabel('Trial')
    ylabel('Approach Likelihood')
    title(['Correct Approach Likelihood'])
    legend([{'Congruent'},{'Incongruent'}])
    
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figSavePath filesep 'approachLikelihood_Correct_Sessions_' flags.SessionN '_Strain_' flags.Genotype],'svg')
    
    figure('Units','normalized','Position',[0 0 1 1])
    plot(mean(RegincorM,2),'r-o','LineWidth',3); hold on; plot(mean(RevincorM,2),'b-o','LineWidth',3);
    er = errorbar([mean(RegincorM,2) mean(RevincorM,2)],[RegincorSEM' RevincorSEM'],'LineWidth',3);  
    er(1).Color = 'r'; er(2).Color = 'b';                          
    er(1).LineStyle = 'none'; er(2).LineStyle = 'none'; 
    ylim([0 1.1])
    xlabel('Trial')
    ylabel('Approach Likelihood')
    title(['Incorrect Approach Likelihood'])
    legend([{'Congruent'},{'Incongruent'}])
    
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figSavePath filesep 'approachLikelihood_Incorrect_Sessions_' flags.SessionN '_Strain_' flags.Genotype],'svg')

    % Put data into structure

    approachStatistics.correctCongruent = RegCorM;     approachStatistics.incorrectCongruent = RegincorM;
    approachStatistics.correctIncongruent = RevCorM;    approachStatistics.incorrectIncongruent = RevincorM;

elseif flags.combined == 1
    % Approach Likelihood. Correct and Incorrect Approaches Combined 
    % Pull data
    RegApproachVector = RegcorApproachVector > 0 | RegincApproachVector > 0;
    RevApproachVector = RevcorApproachVector > 0 | RevincApproachVector > 0;
    RegApp = movmean(RegApproachVector,3); 
    RevApp = movmean(RevApproachVector,3); 
    
    % Plotting data (1 figure since combined)
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
    legend([{'Congruent'},{'Incongruent'}])   
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)

    saveas(gca,[figSavePath filesep 'approachLikelihood_combined_Sessions_' flags.SessionN '_Strain_' flags.Genotype],'svg')

    % Put data into structure
    approachStatistics.combinedCongruent = RegApp;
    approachStatistics.combinedIncongruent = RevApp;

end


end