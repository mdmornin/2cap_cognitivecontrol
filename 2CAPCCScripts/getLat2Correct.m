function [trlStruct,Lat2Correct] = getLat2Correct(trlStruct,figSavePath,masterTbl,flags)
% Find latency to correct
% Also find other metrics related to animals' posture, VTE

% First, find each corrections corresponding trial times 
% Based on the trial time, first find when the animal approaches the wrong
% sipper and then find how long it takes the animal to approach the correct
% sipper
% Once a basic latency metric is obtained, we can look at the angle between
% the animal's snout and head and see when and how quickly it rotates (in
% radians) once the animal realizes they are at the wrong sipper port. 

% Params
    regLatCorrections = [];
    revLatCorrections = [];
    LeftTrials = 1:24;
    RightTrials = 25:48;
    disThresh = 9;
    disRun = 3;
    sipDescent = 10*30;
    sipAscent = 18*30;
    cueOn = 5*30;
    b4cueOn = 3*30; % Some of Nick's analyses focus on moment 2 seconds before cue onset
    if strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'all')    
        for i = 1:length(trlStruct)
            if startsWith(masterTbl.SessionType{i},['Regular'])
                correctionsVector = trlStruct(i).approach(1:48,1:2,1);
                correctionsIdx = correctionsVector(:,1) == 1 & correctionsVector(:,2) == 1;
                correctionsDist_helper = [trlStruct(i).trlLSipDist; trlStruct(i).trlRSipDist];
                correctionsDist_helper = reshape(correctionsDist_helper',[751,96,2]);
                correctionsDist = correctionsDist_helper(:,correctionsIdx,:);
                correctionsDist = correctionsDist(sipDescent:end,:,:);
                [~,CorrectionLats] = max(correctionsDist < disThresh,[],1);
                correctionsTimes = (trlStruct(i).trialTimes(correctionsIdx == 1));
                trlStruct(i).correctionLats = CorrectionLats;
                trlStruct(i).correctionTimes = correctionsTimes;
                regLatCorrections = [regLatCorrections abs(diff(CorrectionLats,1,3))];
                trlStruct(i).correctionLatsD = abs(diff(CorrectionLats,1,3));

            elseif startsWith(masterTbl.SessionType{i},['Reversal'])
                correctionsVector = trlStruct(i).approach(1:48,1:2,1);
                correctionsIdx = correctionsVector(:,1) == 1 & correctionsVector(:,2) == 1;
                correctionsDist_helper = [trlStruct(i).trlLSipDist; trlStruct(i).trlRSipDist];
                correctionsDist_helper = reshape(correctionsDist_helper',[751,96,2]);
                correctionsDist = correctionsDist_helper(:,correctionsIdx,:);
                correctionsDist = correctionsDist(sipDescent:end,:,:);
                [~,CorrectionLats] = max(correctionsDist < disThresh,[],1);
                correctionsTimes = (trlStruct(i).trialTimes(correctionsIdx == 1));
                trlStruct(i).correctionLats = CorrectionLats;
                trlStruct(i).correctionTimes = correctionsTimes;
                revLatCorrections = [revLatCorrections abs(diff(CorrectionLats,1,3))];
                trlStruct(i).correctionLatsD = abs(diff(CorrectionLats,1,3));

            end
        end
    elseif strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'P') 
        for i = 1:length(trlStruct)
            if startsWith(masterTbl.SessionType{i},['Regular']) && strcmp(masterTbl.Strain{i},'P') 
                correctionsVector = trlStruct(i).approach(1:48,1:2,1);
                correctionsIdx = correctionsVector(:,1) == 1 & correctionsVector(:,2) == 1;
                correctionsDist_helper = [trlStruct(i).trlLSipDist; trlStruct(i).trlRSipDist];
                correctionsDist_helper = reshape(correctionsDist_helper',[751,96,2]);
                correctionsDist = correctionsDist_helper(:,correctionsIdx,:);
                correctionsDist = correctionsDist(sipDescent:end,:,:);
                [~,CorrectionLats] = max(correctionsDist < disThresh,[],1);
                correctionsTimes = (trlStruct(i).trialTimes(correctionsIdx == 1));
                trlStruct(i).correctionLats = CorrectionLats;
                trlStruct(i).correctionTimes = correctionsTimes;
                regLatCorrections = [regLatCorrections abs(diff(CorrectionLats,1,3))];
                trlStruct(i).correctionLatsD = abs(diff(CorrectionLats,1,3));

            elseif startsWith(masterTbl.SessionType{i},['Reversal']) && strcmp(masterTbl.Strain{i},'P') 
                correctionsVector = trlStruct(i).approach(1:48,1:2,1);
                correctionsIdx = correctionsVector(:,1) == 1 & correctionsVector(:,2) == 1;
                correctionsDist_helper = [trlStruct(i).trlLSipDist; trlStruct(i).trlRSipDist];
                correctionsDist_helper = reshape(correctionsDist_helper',[751,96,2]);
                correctionsDist = correctionsDist_helper(:,correctionsIdx,:);
                correctionsDist = correctionsDist(sipDescent:end,:,:);
                [~,CorrectionLats] = max(correctionsDist < disThresh,[],1);
                correctionsTimes = (trlStruct(i).trialTimes(correctionsIdx == 1));
                trlStruct(i).correctionLats = CorrectionLats;
                trlStruct(i).correctionTimes = correctionsTimes;
                revLatCorrections = [revLatCorrections abs(diff(CorrectionLats,1,3))];
                trlStruct(i).correctionLatsD = abs(diff(CorrectionLats,1,3));

            end    
        end
    elseif strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'W') 
        for i = 1:length(trlStruct)
            if startsWith(masterTbl.SessionType{i},['Regular']) && strcmp(masterTbl.Strain{i},'W') 
                correctionsVector = trlStruct(i).approach(1:48,1:2,1);
                correctionsIdx = correctionsVector(:,1) == 1 & correctionsVector(:,2) == 1;
                correctionsDist_helper = [trlStruct(i).trlLSipDist; trlStruct(i).trlRSipDist];
                correctionsDist_helper = reshape(correctionsDist_helper',[751,96,2]);
                correctionsDist = correctionsDist_helper(:,correctionsIdx,:);
                correctionsDist = correctionsDist(sipDescent:end,:,:);
                [~,CorrectionLats] = max(correctionsDist < disThresh,[],1);
                correctionsTimes = (trlStruct(i).trialTimes(correctionsIdx == 1));
                trlStruct(i).correctionLats = CorrectionLats;
                trlStruct(i).correctionTimes = correctionsTimes;
                regLatCorrections = [regLatCorrections abs(diff(CorrectionLats,1,3))];
                trlStruct(i).correctionLatsD = abs(diff(CorrectionLats,1,3));

            elseif startsWith(masterTbl.SessionType{i},['Reversal']) && strcmp(masterTbl.Strain{i},'W') 
                correctionsVector = trlStruct(i).approach(1:48,1:2,1);
                correctionsIdx = correctionsVector(:,1) == 1 & correctionsVector(:,2) == 1;
                correctionsDist_helper = [trlStruct(i).trlLSipDist; trlStruct(i).trlRSipDist];
                correctionsDist_helper = reshape(correctionsDist_helper',[751,96,2]);
                correctionsDist = correctionsDist_helper(:,correctionsIdx,:);
                correctionsDist = correctionsDist(sipDescent:end,:,:);
                [~,CorrectionLats] = max(correctionsDist < disThresh,[],1);
                correctionsTimes = (trlStruct(i).trialTimes(correctionsIdx == 1));
                trlStruct(i).correctionLats = CorrectionLats;
                trlStruct(i).correctionTimes = correctionsTimes;
                revLatCorrections = [revLatCorrections abs(diff(CorrectionLats,1,3))];
                trlStruct(i).correctionLatsD = abs(diff(CorrectionLats,1,3));


            end    
        end
    end

    trlStruct = trlStruct;
    Lat2Correct.correctionLats = CorrectionLats;
    Lat2Correct.correctionTimes = correctionsTimes;

    % Plot Latency to Correct
    % First, divide by 30 to put into seconds
    regLatCorrections = regLatCorrections/30;
    revLatCorrections = revLatCorrections/30;

    semregcorLat = std(regLatCorrections)/sqrt(length(regLatCorrections));
    semrevcorLat = std(revLatCorrections)/sqrt(length(revLatCorrections));
    

    figure('Units','normalized','Position',[0 0 1 1])
    bar(categorical({'Congruent Correction Latencies','Incongruent Correction Latencies'}),[mean(regLatCorrections) mean(revLatCorrections)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);
%     scatterXData = [linspace(offsetPos(1)-0.1,offsetPos(1)+0.1,length(ReglatCor)),linspace(offsetPos(2)-0.1,offsetPos(2)+0.1,length(RevlatCor))];

    hold on
    scatter([linspace(0.95,1.05,length(regLatCorrections)) linspace(1.95,2.05,length(revLatCorrections))],[regLatCorrections revLatCorrections],50,'ko','LineWidth',2)
    er = errorbar(categorical({'Congruent Correction Latencies','Incongruent Correction Latencies'}),[mean(regLatCorrections) mean(revLatCorrections)],[semregcorLat semrevcorLat],'LineWidth',3);    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
    ylabel('Latency to Correct (s)')
    
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figSavePath filesep 'correctionLatencies_Sessions_' flags.SessionN '_Strain_' flags.Genotype],flags.figextension)
end

