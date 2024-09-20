function [correctionLH,trlStruct] = getCorrectionsLH(trlStruct,figSavePath,masterTbl,flags)
%% Likelihood to Correct over Trials
% Similar analysis to likelihood to approach. Changed logical operator on
% the data from OR to AND when looking at the approach variable. The
% approach variable looks like this each trial: [1 0] [0 1] [1 1] [0 0].
% Currently ASSUMING [1 1] indicates a correction. Will need a finer
% grained analysis on this eventually (02/03/2022).


    % Initialize blank variables
    RegCorrectionApproachVector = []; RevCorrectionApproachVector = [];
    RegCorrectionApproachVectorInc = []; RevCorrectionApproachVectorInc = [];
    
    % Pull data from structure
    if strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'all') 
        
        for i = 1:length(trlStruct)

            if startsWith(masterTbl.SessionType{i},['Regular']) %&& strcmp(masterTbl.Strain{i},'P')
                currentTrialTimes = trlStruct(i).trialTimes;
                [~,sortedIndex] = sort(currentTrialTimes(1:48));
                currentApproachVector = trlStruct(i).approach(1:48,:,1);
                sortedCorrectApproach = currentApproachVector(sortedIndex,1);
                sortedIncrectApproach = currentApproachVector(sortedIndex,2);
                RegCorrectionApproachVector = [RegCorrectionApproachVector sortedCorrectApproach];
                RegCorrectionApproachVectorInc = [RegCorrectionApproachVectorInc sortedIncrectApproach];
                trlStruct(i).changePointCorrections = findchangepts(movmean(sortedCorrectApproach > 0 & sortedIncrectApproach > 0,3));
                trlStruct(i).correctionLH = sortedCorrectApproach > 0 & sortedIncrectApproach > 0;

        
            elseif startsWith(masterTbl.SessionType{i},['Reversal']) %&& strcmp(masterTbl.Strain{i},'P')
                currentTrialTimes = trlStruct(i).trialTimes;
                [~,sortedIndex] = sort(currentTrialTimes(1:48));
                currentApproachVector = trlStruct(i).approach(1:48,:,1);
                sortedCorrectApproach = currentApproachVector(sortedIndex,1);
                sortedIncrectApproach = currentApproachVector(sortedIndex,2);
                RevCorrectionApproachVector = [RevCorrectionApproachVector sortedCorrectApproach];
                RevCorrectionApproachVectorInc = [RevCorrectionApproachVectorInc sortedIncrectApproach];
                trlStruct(i).changePointCorrections = findchangepts(movmean(sortedCorrectApproach > 0 & sortedIncrectApproach > 0,3));
                trlStruct(i).correctionLH = sortedCorrectApproach > 0 & sortedIncrectApproach > 0;

            end

        end

    elseif strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'P')

        for i = 1:length(trlStruct)

            if startsWith(masterTbl.SessionType{i},['Regular']) && strcmp(masterTbl.Strain{i},'P') 
                currentTrialTimes = trlStruct(i).trialTimes;
                [~,sortedIndex] = sort(currentTrialTimes(1:48));
                currentApproachVector = trlStruct(i).approach(1:48,:,1);
                sortedCorrectApproach = currentApproachVector(sortedIndex,1);
                sortedIncrectApproach = currentApproachVector(sortedIndex,2);
                RegCorrectionApproachVector = [RegCorrectionApproachVector sortedCorrectApproach];
                RegCorrectionApproachVectorInc = [RegCorrectionApproachVectorInc sortedIncrectApproach];
                trlStruct(i).changePointCorrections = findchangepts(movmean(sortedCorrectApproach > 0 & sortedIncrectApproach > 0,3));
                trlStruct(i).correctionLH = sortedCorrectApproach > 0 & sortedIncrectApproach > 0;

        
            elseif startsWith(masterTbl.SessionType{i},['Reversal']) && strcmp(masterTbl.Strain{i},'P') 
                currentTrialTimes = trlStruct(i).trialTimes;
                [~,sortedIndex] = sort(currentTrialTimes(1:48));
                currentApproachVector = trlStruct(i).approach(1:48,:,1);
                sortedCorrectApproach = currentApproachVector(sortedIndex,1);
                sortedIncrectApproach = currentApproachVector(sortedIndex,2);
                RevCorrectionApproachVector = [RevCorrectionApproachVector sortedCorrectApproach];
                RevCorrectionApproachVectorInc = [RevCorrectionApproachVectorInc sortedIncrectApproach];
                trlStruct(i).changePointCorrections = findchangepts(movmean(sortedCorrectApproach > 0 & sortedIncrectApproach > 0,3));
                trlStruct(i).correctionLH = sortedCorrectApproach > 0 & sortedIncrectApproach > 0;

            end

        end

    elseif strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'W')

        for i = 1:length(trlStruct)

            if startsWith(masterTbl.SessionType{i},['Regular']) && strcmp(masterTbl.Strain{i},'W') 
                currentTrialTimes = trlStruct(i).trialTimes;
                [~,sortedIndex] = sort(currentTrialTimes(1:48));
                currentApproachVector = trlStruct(i).approach(1:48,:,1);
                sortedCorrectApproach = currentApproachVector(sortedIndex,1);
                sortedIncrectApproach = currentApproachVector(sortedIndex,2);
                RegCorrectionApproachVector = [RegCorrectionApproachVector sortedCorrectApproach];
                RegCorrectionApproachVectorInc = [RegCorrectionApproachVectorInc sortedIncrectApproach];
                trlStruct(i).changePointCorrections = findchangepts(movmean(sortedCorrectApproach > 0 & sortedIncrectApproach > 0,3));
                trlStruct(i).correctionLH = sortedCorrectApproach > 0 & sortedIncrectApproach > 0;
        
            elseif startsWith(masterTbl.SessionType{i},['Reversal']) && strcmp(masterTbl.Strain{i},'W') 
                currentTrialTimes = trlStruct(i).trialTimes;
                [~,sortedIndex] = sort(currentTrialTimes(1:48));
                currentApproachVector = trlStruct(i).approach(1:48,:,1);
                sortedCorrectApproach = currentApproachVector(sortedIndex,1);
                sortedIncrectApproach = currentApproachVector(sortedIndex,2);
                RevCorrectionApproachVector = [RevCorrectionApproachVector sortedCorrectApproach];
                RevCorrectionApproachVectorInc = [RevCorrectionApproachVectorInc sortedIncrectApproach];
                trlStruct(i).changePointCorrections = findchangepts(movmean(sortedCorrectApproach > 0 & sortedIncrectApproach > 0,3));
                trlStruct(i).correctionLH = sortedCorrectApproach > 0 & sortedIncrectApproach > 0;
            end

        end

    end
    
    
    % Collect data where there is both an incorrect AND correct approach
    RegApproachVector = RegCorrectionApproachVector > 0 & RegCorrectionApproachVectorInc > 0;
    RevApproachVector = RevCorrectionApproachVector > 0 & RevCorrectionApproachVectorInc > 0;
    
    
    % Take a 3-trial moving average to 'smooth' data
    RegApp = movmean(RegApproachVector,3); 
    RevApp = movmean(RevApproachVector,3); 
    
    % Calculate SEM for plots
    RegAppSEM = std(RegApp')/sqrt(size(RegApp,2)); 
    RevAppSEM = std(RevApp')/sqrt(size(RevApp,2)); 
    
    % Plotting!
    figure('Units','normalized','Position',[0 0 1 1])
    plot(mean(RegApp,2),'r-o','LineWidth',3); hold on; plot(mean(RevApp,2),'b-o','LineWidth',3);
    er = errorbar([mean(RegApp,2) mean(RevApp,2)],[RegAppSEM' RevAppSEM'],'LineWidth',3); 
    er(1).Color = 'r'; er(2).Color = 'b';                          
    er(1).LineStyle = 'none'; er(2).LineStyle = 'none'; 
    ylim([0 1.1])
    xlabel('Trial')
    ylabel('Correction Likelihood')
    title(['Correction Likelihood, All Sessions, All Genotypes'])
    legend([{'Congruent'},{'Incongruent'}])
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
    % Saving!
    saveas(gca,[figSavePath filesep 'correctionsLikelihood_Sessions_' flags.SessionN '_Strain_' flags.Genotype],'png')


    % Outputting!
    correctionLH.CongruentCorrectionLH = RegApp;
    correctionLH.IncongruentCorrectionLH = RevApp;

    % Done!
end