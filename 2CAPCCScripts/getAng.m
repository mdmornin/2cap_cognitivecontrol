function [trlStruct,angStruct] = getAng(trlStruct,figSavePath,masterTbl,flags)

% Params
    LeftTrials = 1:24;
    RightTrials = 25:48;
    disThresh = 9;
    disRun = 3;
    sipDescent = 10*30;
    sipAscent = 18*30;
    cueOn = 5*30;
% trlStruct contains information relevant to finding the angular velocity
% as the animal corrects. 
% Given the head and nose points, find the line between them at the first
% time point where the animal is at the wrong port and make that the
% reference. Then track the change in angle from that reference until it is
% opposite. Measure the time it took to make this change. 
% Animal does not actually make a drastic change in nose angle. Changing
% approach to finding the difference between a straight, optimal line and
% the path the animal actually takes from incorrect to the correct sipper.
% Then need to check if this is correlated with latency. If not, it could
% be an important measure. 

% Blank vars
regAUC = [];
revAUC = [];

    if strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'all')    

        for i = 1:length(trlStruct)            
            if ~isempty(trlStruct(i).correctionLatsD)
            if startsWith(masterTbl.SessionType{i},['Regular'])

                correctionsVector = trlStruct(i).approach(1:48,1:2,1);
                correctionsIdx = correctionsVector(:,1) == 1 & correctionsVector(:,2) == 1;    
   
                headIdx = trlStruct(i).trlheadPos(correctionsIdx,sipDescent:end,:);
                noseIdx = trlStruct(i).trlnosePos(correctionsIdx,sipDescent:end,:);
                
                correctionTrlTimes = squeeze(trlStruct(i).correctionLats);
    
                RightSipperPos = trlStruct(i).boxCoords(1,:);
                LeftSipperPos = trlStruct(i).boxCoords(2,:);
                
                xsip = linspace(LeftSipperPos(1),RightSipperPos(1),100);
                ysip = linspace(LeftSipperPos(2),RightSipperPos(2),100);
                
                sipLine = [xsip; ysip];

                for k = 1:length(trlStruct(i).correctionLatsD)

                    if length(trlStruct(i).correctionLatsD) > 1
                        trialTraj = squeeze(noseIdx(k,min(correctionTrlTimes(k,:)):max(correctionTrlTimes(k,:)),:))';
                    else
                        trialTraj = squeeze(noseIdx(k,min(correctionTrlTimes):max(correctionTrlTimes),:))';
                    end

                    distBestLine = min(pdist2(sipLine',trialTraj'));
                    aucDistBestLine = max(distBestLine);
    
                    trlStruct(i).trialTraj{k} = trialTraj;
                    trlStruct(i).distBestLine{k} = distBestLine;
                    trlStruct(i).aucDistBestLine(k) = aucDistBestLine;
                    regAUC = [regAUC aucDistBestLine];
                end

            elseif startsWith(masterTbl.SessionType{i},['Reversal'])
                correctionsVector = trlStruct(i).approach(1:48,1:2,1);
                correctionsIdx = correctionsVector(:,1) == 1 & correctionsVector(:,2) == 1;    
    
                headIdx = trlStruct(i).trlheadPos(correctionsIdx,sipDescent:end,:);
                noseIdx = trlStruct(i).trlnosePos(correctionsIdx,sipDescent:end,:);
                
                correctionTrlTimes = squeeze(trlStruct(i).correctionLats);
    
                RightSipperPos = trlStruct(i).boxCoords(1,:);
                LeftSipperPos = trlStruct(i).boxCoords(2,:);
                
                xsip = linspace(LeftSipperPos(1),RightSipperPos(1),100);
                ysip = linspace(LeftSipperPos(2),RightSipperPos(2),100);
                
                sipLine = [xsip; ysip];

                for k = 1:length(trlStruct(i).correctionLatsD)
                    if length(trlStruct(i).correctionLatsD) > 1
                        trialTraj = squeeze(noseIdx(k,min(correctionTrlTimes(k,:)):max(correctionTrlTimes(k,:)),:))';
                    else
                        trialTraj = squeeze(noseIdx(k,min(correctionTrlTimes):max(correctionTrlTimes),:))';
                    end
                    distBestLine = min(pdist2(sipLine',trialTraj'));
                    aucDistBestLine = max(distBestLine);
    
                    trlStruct(i).trialTraj{k} = trialTraj;
                    trlStruct(i).distBestLine{k} = distBestLine; 
                    trlStruct(i).aucDistBestLine(k) = aucDistBestLine;
                    revAUC = [revAUC aucDistBestLine];
                end

            end
            end
        end
        

        
    elseif strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'P') 
        for i = 1:length(trlStruct)
            disp(i)
            
            if ~isempty(trlStruct(i).correctionLatsD)
            
            if startsWith(masterTbl.SessionType{i},['Regular'])  && strcmp(masterTbl.Strain{i},'P') 
            
                correctionsVector = trlStruct(i).approach(1:48,1:2,1);
                correctionsIdx = correctionsVector(:,1) == 1 & correctionsVector(:,2) == 1;    
    
                headIdx = trlStruct(i).trlheadPos(correctionsIdx,sipDescent:end,:);
                noseIdx = trlStruct(i).trlnosePos(correctionsIdx,sipDescent:end,:);
                
                correctionTrlTimes = squeeze(trlStruct(i).correctionLats);
    
                RightSipperPos = trlStruct(i).boxCoords(1,:);
                LeftSipperPos = trlStruct(i).boxCoords(2,:);
                
                xsip = linspace(LeftSipperPos(1),RightSipperPos(1),100);
                ysip = linspace(LeftSipperPos(2),RightSipperPos(2),100);
                
                sipLine = [xsip; ysip];
    
                for k = 1:length(trlStruct(i).correctionLatsD)
                    if length(trlStruct(i).correctionLatsD) > 1
                        trialTraj = squeeze(noseIdx(k,min(correctionTrlTimes(k,:)):max(correctionTrlTimes(k,:)),:))';
                    else
                        trialTraj = squeeze(noseIdx(k,min(correctionTrlTimes):max(correctionTrlTimes),:))';
                    end
                    distBestLine = min(pdist2(sipLine',trialTraj'));
                    aucDistBestLine = max(distBestLine);
    
                    trlStruct(i).trialTraj{k} = trialTraj;
                    trlStruct(i).distBestLine{k} = distBestLine;
                    trlStruct(i).aucDistBestLine(k) = aucDistBestLine;
                    regAUC = [regAUC aucDistBestLine];
                end

            
            elseif startsWith(masterTbl.SessionType{i},['Reversal'])  && strcmp(masterTbl.Strain{i},'P') 
                correctionsVector = trlStruct(i).approach(1:48,1:2,1);
                correctionsIdx = correctionsVector(:,1) == 1 & correctionsVector(:,2) == 1;    
    
                headIdx = trlStruct(i).trlheadPos(correctionsIdx,sipDescent:end,:);
                noseIdx = trlStruct(i).trlnosePos(correctionsIdx,sipDescent:end,:);
                
                correctionTrlTimes = squeeze(trlStruct(i).correctionLats);
    
                RightSipperPos = trlStruct(i).boxCoords(1,:);
                LeftSipperPos = trlStruct(i).boxCoords(2,:);
                
                xsip = linspace(LeftSipperPos(1),RightSipperPos(1),100);
                ysip = linspace(LeftSipperPos(2),RightSipperPos(2),100);
                
                sipLine = [xsip; ysip];
    
                for k = 1:length(trlStruct(i).correctionLatsD)
                    if length(trlStruct(i).correctionLatsD) > 1
                        trialTraj = squeeze(noseIdx(k,min(correctionTrlTimes(k,:)):max(correctionTrlTimes(k,:)),:))';
                    else
                        trialTraj = squeeze(noseIdx(k,min(correctionTrlTimes):max(correctionTrlTimes),:))';
                    end
                    distBestLine = min(pdist2(sipLine',trialTraj'));
                    aucDistBestLine = max(distBestLine);
    
                    trlStruct(i).trialTraj{k} = trialTraj;
                    trlStruct(i).distBestLine{k} = distBestLine; 
                    trlStruct(i).aucDistBestLine(k) = aucDistBestLine;
                    revAUC = [revAUC aucDistBestLine];
                end

            end   
            end
        end
    elseif strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'W') 

        for i = 1:length(trlStruct)
            disp(i)
            if ~isempty(trlStruct(i).correctionLatsD)
            
            if startsWith(masterTbl.SessionType{i},['Regular']) && strcmp(masterTbl.Strain{i},'W') 
            
                correctionsVector = trlStruct(i).approach(1:48,1:2,1);
                correctionsIdx = correctionsVector(:,1) == 1 & correctionsVector(:,2) == 1;    
    
                headIdx = trlStruct(i).trlheadPos(correctionsIdx,sipDescent:end,:);
                noseIdx = trlStruct(i).trlnosePos(correctionsIdx,sipDescent:end,:);
                
                correctionTrlTimes = squeeze(trlStruct(i).correctionLats);
    
                RightSipperPos = trlStruct(i).boxCoords(1,:);
                LeftSipperPos = trlStruct(i).boxCoords(2,:);
                
                xsip = linspace(LeftSipperPos(1),RightSipperPos(1),100);
                ysip = linspace(LeftSipperPos(2),RightSipperPos(2),100);
                
                sipLine = [xsip; ysip];
    
                for k = 1:length(trlStruct(i).correctionLatsD)
                    if length(trlStruct(i).correctionLatsD) > 1
                        trialTraj = squeeze(noseIdx(k,min(correctionTrlTimes(k,:)):max(correctionTrlTimes(k,:)),:))';
                    else
                        trialTraj = squeeze(noseIdx(k,min(correctionTrlTimes):max(correctionTrlTimes),:))';
                    end
                    distBestLine = min(pdist2(sipLine',trialTraj'));
                    aucDistBestLine = max(distBestLine);
    
                    trlStruct(i).trialTraj{k} = trialTraj;
                    trlStruct(i).distBestLine{k} = distBestLine;
                    trlStruct(i).aucDistBestLine(k) = aucDistBestLine;
                    regAUC = [regAUC aucDistBestLine];
                end

            
            elseif startsWith(masterTbl.SessionType{i},['Reversal']) && strcmp(masterTbl.Strain{i},'W') 
                correctionsVector = trlStruct(i).approach(1:48,1:2,1);
                correctionsIdx = correctionsVector(:,1) == 1 & correctionsVector(:,2) == 1;    
    
                headIdx = trlStruct(i).trlheadPos(correctionsIdx,sipDescent:end,:);
                noseIdx = trlStruct(i).trlnosePos(correctionsIdx,sipDescent:end,:);
                
                correctionTrlTimes = squeeze(trlStruct(i).correctionLats);
    
                RightSipperPos = trlStruct(i).boxCoords(1,:);
                LeftSipperPos = trlStruct(i).boxCoords(2,:);
                
                xsip = linspace(LeftSipperPos(1),RightSipperPos(1),100);
                ysip = linspace(LeftSipperPos(2),RightSipperPos(2),100);
                
                sipLine = [xsip; ysip];
    
                for k = 1:length(trlStruct(i).correctionLatsD)
                    if length(trlStruct(i).correctionLatsD) > 1
                        trialTraj = squeeze(noseIdx(k,min(correctionTrlTimes(k,:)):max(correctionTrlTimes(k,:)),:))';
                    else
                        trialTraj = squeeze(noseIdx(k,min(correctionTrlTimes):max(correctionTrlTimes),:))';
                    end
                    distBestLine = min(pdist2(sipLine',trialTraj'));
                    aucDistBestLine = max(distBestLine);
    
                    trlStruct(i).trialTraj{k} = trialTraj;
                    trlStruct(i).distBestLine{k} = distBestLine;
                    trlStruct(i).aucDistBestLine(k) = aucDistBestLine;
                    revAUC = [revAUC aucDistBestLine];
                end

            end        
            end
        end

    end

    %% Plot data

    semregAUC = std(regAUC)/sqrt(length(regAUC));
    semrevAUC = std(revAUC)/sqrt(length(revAUC));

    figure('Units','normalized','Position',[0 0 1 1])
    bar(categorical({'Congruent MaxVal','Incongruent MaxVal'}),[mean(regAUC) mean(revAUC)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);
%     scatterXData = [linspace(offsetPos(1)-0.1,offsetPos(1)+0.1,length(ReglatCor)),linspace(offsetPos(2)-0.1,offsetPos(2)+0.1,length(RevlatCor))];

    hold on
    scatter([linspace(0.95,1.05,length(regAUC)) linspace(1.95,2.05,length(revAUC))],[regAUC revAUC],50,'ko','LineWidth',2)
    er = errorbar(categorical({'Congruent MaxVal','Incongruent MaxVal'}),[mean(regAUC) mean(revAUC)],[semregAUC semrevAUC],'LineWidth',3);    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
    ylabel('Max Distance from Most Efficient Line')
    
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figSavePath filesep 'MaxValBestLine_Sessions_' flags.SessionN '_Strain_' flags.Genotype],flags.figextension)
    %% Output data
    angStruct.regAUC = regAUC;
    angStruct.revAUC = revAUC;
end