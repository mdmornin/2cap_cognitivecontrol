function [trlStruct,errorInfoSt] = getErrorInfo(trlStruct,masterTbl,figSavePath,flags)
    LeftTrials = 1:24;
    RightTrials = 25:48;
    
    regErrorsD = [];    regPers = [];
    revErrorsD = [];    revPers = [];
    if strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'all') 
            for i = 1:length(trlStruct)
                               
                if startsWith(masterTbl.SessionType{i},['Regular'])
                    vecErrors = trlStruct(i).approach(1:48,2,1);
                    trlTimes = trlStruct(i).trialTimes(1:48);
                    [~,sortedIndex] = sort(trlTimes);
                    sortedErrors = vecErrors(sortedIndex);
                    labeledMatrix = bwlabel(sortedErrors);
                    measurements = regionprops(labeledMatrix, 'Area');
                    rightErrors = sum(vecErrors(LeftTrials));
                    leftErrors = sum(vecErrors(RightTrials));
                    pers = [measurements.Area];
                    regPers = [regPers pers];
                    regErrorsD = [regErrorsD (rightErrors - leftErrors)];
                    trlStruct(i).Pers = pers;
                    trlStruct(i).errorD = (rightErrors - leftErrors);
                elseif startsWith(masterTbl.SessionType{i},['Reversal'])
                    vecErrors = trlStruct(i).approach(1:48,1,1);
                    trlTimes = trlStruct(i).trialTimes(1:48);
                    [~,sortedIndex] = sort(trlTimes);
                    sortedErrors = vecErrors(sortedIndex);
                    labeledMatrix = bwlabel(sortedErrors);
                    measurements = regionprops(labeledMatrix, 'Area');                    
                    leftErrors = sum(vecErrors(LeftTrials));
                    rightErrors = sum(vecErrors(RightTrials));
                    pers = [measurements.Area];
                    revPers = [revPers pers];
                    revErrorsD = [revErrorsD (rightErrors - leftErrors)];
                    trlStruct(i).Pers = pers;
                    trlStruct(i).errorD = (rightErrors - leftErrors);
                end
            end
    elseif strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'P')
            for i = 1:length(trlStruct)
                               
                if startsWith(masterTbl.SessionType{i},['Regular']) && strcmp(masterTbl.Strain{i},'P')
                    vecErrors = trlStruct(i).approach(1:48,2,1);
                    trlTimes = trlStruct(i).trialTimes(1:48);
                    [~,sortedIndex] = sort(trlTimes);
                    sortedErrors = vecErrors(sortedIndex);
                    labeledMatrix = bwlabel(sortedErrors);
                    measurements = regionprops(labeledMatrix, 'Area');                                        
                    rightErrors = sum(vecErrors(LeftTrials));
                    leftErrors = sum(vecErrors(RightTrials));
                    pers = [measurements.Area];
                    regPers = [regPers pers];
                    regErrorsD = [regErrorsD (rightErrors - leftErrors)];   
                    trlStruct(i).Pers = pers;
                    trlStruct(i).errorD = (rightErrors - leftErrors);
                elseif startsWith(masterTbl.SessionType{i},['Reversal']) && strcmp(masterTbl.Strain{i},'P')
                    vecErrors = trlStruct(i).approach(1:48,1,1);
                    trlTimes = trlStruct(i).trialTimes(1:48);
                    [~,sortedIndex] = sort(trlTimes);
                    sortedErrors = vecErrors(sortedIndex);
                    labeledMatrix = bwlabel(sortedErrors);
                    measurements = regionprops(labeledMatrix, 'Area');                   
                    leftErrors = sum(vecErrors(LeftTrials));
                    rightErrors = sum(vecErrors(RightTrials));
                    pers = [measurements.Area];
                    revPers = [revPers pers];
                    revErrorsD = [revErrorsD (rightErrors - leftErrors)];  
                    trlStruct(i).Pers = pers;
                    trlStruct(i).errorD = (rightErrors - leftErrors);
                end
            end 
    elseif strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'W')
            for i = 1:length(trlStruct)
                
                if startsWith(masterTbl.SessionType{i},['Regular']) && strcmp(masterTbl.Strain{i},'W')
                    vecErrors = trlStruct(i).approach(1:48,2,1);
                    trlTimes = trlStruct(i).trialTimes(1:48);
                    [~,sortedIndex] = sort(trlTimes);
                    sortedErrors = vecErrors(sortedIndex);
                    labeledMatrix = bwlabel(sortedErrors);
                    measurements = regionprops(labeledMatrix, 'Area');                                       
                    rightErrors = sum(vecErrors(LeftTrials));
                    leftErrors = sum(vecErrors(RightTrials));
                    pers = [measurements.Area];
                    regPers = [regPers pers];
                    regErrorsD = [regErrorsD (rightErrors - leftErrors)];
                    trlStruct(i).Pers = pers;
                    trlStruct(i).errorD = (rightErrors - leftErrors);
                elseif startsWith(masterTbl.SessionType{i},['Reversal']) && strcmp(masterTbl.Strain{i},'W')
                    vecErrors = trlStruct(i).approach(1:48,1,1);
                    trlTimes = trlStruct(i).trialTimes(1:48);
                    [~,sortedIndex] = sort(trlTimes);
                    sortedErrors = vecErrors(sortedIndex);
                    labeledMatrix = bwlabel(sortedErrors);
                    measurements = regionprops(labeledMatrix, 'Area');                                        
                    leftErrors = sum(vecErrors(LeftTrials));
                    rightErrors = sum(vecErrors(RightTrials));
                    pers = [measurements.Area];
                    revPers = [revPers pers];
                    revErrorsD = [revErrorsD (rightErrors - leftErrors)];
                    trlStruct(i).Pers = pers;
                    trlStruct(i).errorD = (rightErrors - leftErrors);
                end
            end
    end

% Plot
    semreg = std(regErrorsD)/sqrt(length(regErrorsD));
    semrev = std(revErrorsD)/sqrt(length(revErrorsD));

    figure('Units','normalized','Position',[0 0 1 1])
    bar(categorical({'Congruent ErrorD','Incongruent ErrorD'}),[mean(regErrorsD) mean(revErrorsD)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);
    hold on
    scatter([linspace(0.95,1.05,length(regErrorsD)) linspace(1.95,2.05,length(revErrorsD))],[regErrorsD revErrorsD],50,'ko','LineWidth',2)
    er = errorbar(categorical({'Congruent ErrorD','Incongruent ErrorD'}),[mean(regErrorsD) mean(revErrorsD)],[semreg semrev],'LineWidth',3);    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
    ylabel('Difference in Error Sides')
    
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figSavePath filesep 'ErrorD_Sessions_' flags.SessionN '_Strain_' flags.Genotype],flags.figextension)

% Plot perserverations
    semreg = std(regPers)/sqrt(length(regPers));
    semrev = std(revPers)/sqrt(length(revPers));

    figure('Units','normalized','Position',[0 0 1 1])
    bar(categorical({'Congruent Perserveration','Incongruent Perserveration'}),[mean(regPers) mean(revPers)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);
    hold on
    scatter([linspace(0.95,1.05,length(regPers)) linspace(1.95,2.05,length(revPers))],[regPers revPers],50,'ko','LineWidth',2)
    er = errorbar(categorical({'Congruent Perserveration','Incongruent Perserveration'}),[mean(regPers) mean(revPers)],[semreg semrev],'LineWidth',3);    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
    ylabel('Errors in a Row')
    ylim([0 10])
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figSavePath filesep 'Pers_Sessions_' flags.SessionN '_Strain_' flags.Genotype],flags.figextension)

% Output
    errorInfoSt.regErrorD = regErrorsD;
    errorInfoSt.revErrorD = revErrorsD;
    errorInfoSt.regPers = regPers;
    errorInfoSt.revPers = revPers;

end