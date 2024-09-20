function [trlStruct,trialSipperTime] = getTrialSipperTime(trlStruct,figSavePath,masterTbl,flags)
    % Parameters
    disThresh = 9;
    disRun = 3;
    sipDescent = 10*30;
    sipAscent = 18*30;
    cueOn = 5*30;
    % Initialize variables
    occReg = []; occRev = []; % Correct and Incorrect latencies
   
    if strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'all')    

        for i = 1:length(trlStruct)
            if startsWith(masterTbl.SessionType{i},['Regular'])
                CorIdx = [trlStruct(i).LcorrectApproachTrls <= disThresh; trlStruct(i).RcorrectApproachTrls <= disThresh];
                currentTrialTimes = trlStruct(i).trialTimes;
                [~,sortedIndex] = sort(currentTrialTimes(1:48));
                CorDbl = double(CorIdx);
                regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
                regOcc = regOcc(sortedIndex);
                trlStruct(i).trlOcc = regOcc;
            elseif startsWith(masterTbl.SessionType{i},['Reversal'])
                CorIdx = [trlStruct(i).LcorrectApproachTrls <= disThresh; trlStruct(i).RcorrectApproachTrls <= disThresh];
                currentTrialTimes = trlStruct(i).trialTimes;
                [~,sortedIndex] = sort(currentTrialTimes(1:48));
                CorDbl = double(CorIdx);
                regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
                regOcc = regOcc(sortedIndex);
                trlStruct(i).trlOcc = regOcc;
            end
        end

    elseif strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'P') 

        for i = 1:length(trlStruct)
            if startsWith(masterTbl.SessionType{i},['Regular'])  && strcmp(masterTbl.Strain{i},'P')
                CorIdx = [trlStruct(i).LcorrectApproachTrls <= disThresh; trlStruct(i).RcorrectApproachTrls <= disThresh];
                currentTrialTimes = trlStruct(i).trialTimes;
                [~,sortedIndex] = sort(currentTrialTimes(1:48));
                CorDbl = double(CorIdx);
                regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
                regOcc = regOcc(sortedIndex);
                trlStruct(i).trlOcc = regOcc;
            elseif startsWith(masterTbl.SessionType{i},['Reversal'])  && strcmp(masterTbl.Strain{i},'P')
                CorIdx = [trlStruct(i).LcorrectApproachTrls <= disThresh; trlStruct(i).RcorrectApproachTrls <= disThresh];
                currentTrialTimes = trlStruct(i).trialTimes;
                [~,sortedIndex] = sort(currentTrialTimes(1:48));
                CorDbl = double(CorIdx);
                regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
                regOcc = regOcc(sortedIndex);
                trlStruct(i).trlOcc = regOcc;
            end
        end

    elseif strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'W') 

        for i = 1:length(trlStruct)
            if startsWith(masterTbl.SessionType{i},['Regular'])  && strcmp(masterTbl.Strain{i},'W')
                CorIdx = [trlStruct(i).LcorrectApproachTrls <= disThresh; trlStruct(i).RcorrectApproachTrls <= disThresh];
                currentTrialTimes = trlStruct(i).trialTimes;
                [~,sortedIndex] = sort(currentTrialTimes(1:48));
                CorDbl = double(CorIdx);
                regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
                regOcc = regOcc(sortedIndex);
                trlStruct(i).trlOcc = regOcc;
            elseif startsWith(masterTbl.SessionType{i},['Reversal'])  && strcmp(masterTbl.Strain{i},'W')
                CorIdx = [trlStruct(i).LcorrectApproachTrls <= disThresh; trlStruct(i).RcorrectApproachTrls <= disThresh];
                currentTrialTimes = trlStruct(i).trialTimes;
                [~,sortedIndex] = sort(currentTrialTimes(1:48));
                CorDbl = double(CorIdx);
                regOcc = sum(CorDbl(:,sipDescent:sipAscent),2);
                regOcc = regOcc(sortedIndex);
                trlStruct(i).trlOcc = regOcc;
            end
        end

    end
    % Pull data
    for i = 1:length(trlStruct)
        if startsWith(masterTbl.SessionType{i},['Regular'])
            occReg = [occReg trlStruct(i).trlOcc];
        elseif startsWith(masterTbl.SessionType{i},['Reversal'])
            occRev = [occRev trlStruct(i).trlOcc];
        end
    end

    % Convert to seconds
    occReg = occReg ./ 30; % Convert time to seconds
    occRev = occRev ./ 30; % Convert time to seconds
    % Find SEM for plots
    regSEM = std(occReg')/sqrt(min(size(occReg)));
    revSEM = std(occRev')/sqrt(min(size(occRev)));
    
    % Plot data
    figure('Units','normalized','Position',[0 0 1 1])
    shadedErrorBar(1:length(sortedIndex),nanmean(occReg,2),regSEM,'lineprops',{'r-o','LineWidth',3});
    hold on
    shadedErrorBar(1:length(sortedIndex),nanmean(occRev,2),revSEM,'lineprops',{'b-o','LineWidth',3});
    xlabel('Trials')
    ylabel('Sipper Occupancy (s)')
    
    title('Time Spent at Sipper per Trial')
    legend([{'Congruent'},{'Incongruent'}])
    xlim([1 48])
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figSavePath filesep 'sipperOccupancy_perTrial_' flags.SessionN '_Strain_' flags.Genotype],'png')

    trialSipperTime.occReg = occReg;
    trialSipperTime.occRev = occRev;
    
end