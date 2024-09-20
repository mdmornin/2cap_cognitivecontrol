function [sipperTimeStats,trlStruct] = getSipperTimePerSession(trlStruct,figSavePath,masterTbl,flags)
    
    % Set blank variables
    regCorOccu = [];   regIncOccu = [];
    revCorOccu = [];   revIncOccu = [];                                                                         
    
    if strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'all')    

        for i = 1:length(trlStruct)   
            if startsWith(masterTbl.SessionType{i},['Regular'])
                regCorOccu = [regCorOccu; trlStruct(i).CorCorOcc];
                regIncOccu = [regIncOccu; trlStruct(i).CorIncOcc];
            elseif startsWith(masterTbl.SessionType{i},['Reversal'])
                revCorOccu = [revCorOccu; trlStruct(i).CorCorOcc];
                revIncOccu = [revIncOccu; trlStruct(i).CorIncOcc];
            end  
        end

    elseif strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'P')

        for i = 1:length(trlStruct)   
            if startsWith(masterTbl.SessionType{i},['Regular']) && strcmp(flags.Genotype,'P')
                regCorOccu = [regCorOccu; trlStruct(i).CorCorOcc];
                regIncOccu = [regIncOccu; trlStruct(i).CorIncOcc];
            elseif startsWith(masterTbl.SessionType{i},['Reversal']) && strcmp(flags.Genotype,'P')
                revCorOccu = [revCorOccu; trlStruct(i).CorCorOcc];
                revIncOccu = [revIncOccu; trlStruct(i).CorIncOcc];
            end  
        end

    elseif strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'W')

        for i = 1:length(trlStruct)   
            if startsWith(masterTbl.SessionType{i},['Regular']) && strcmp(flags.Genotype,'W')
                regCorOccu = [regCorOccu; trlStruct(i).CorCorOcc];
                regIncOccu = [regIncOccu; trlStruct(i).CorIncOcc];
            elseif startsWith(masterTbl.SessionType{i},['Reversal']) && strcmp(flags.Genotype,'W')
                revCorOccu = [revCorOccu; trlStruct(i).CorCorOcc];
                revIncOccu = [revIncOccu; trlStruct(i).CorIncOcc];
            end  
        end

    end
    
    regCorOccu = regCorOccu(~isnan(regCorOccu))./30;
    regIncOccu = regIncOccu(~isnan(regIncOccu))./30;
    revCorOccu = revCorOccu(~isnan(revCorOccu))./30;
    revIncOccu = revIncOccu(~isnan(revIncOccu))./30;
    
    % Plot
    semregcorOcc = std(regCorOccu)/sqrt(length(regCorOccu));
    semrevcorOcc = std(regIncOccu)/sqrt(length(regIncOccu));
    semregIncorOcc = std(revCorOccu)/sqrt(length(revCorOccu));
    semrevIncorOcc = std(revIncOccu)/sqrt(length(revIncOccu));
    
    figure('Units','normalized','Position',[0 0 1 1])
    hb = bar(categorical({'Correct Approach Trials','Incorrect Approach Trials'}),[mean(regCorOccu) mean(revCorOccu); mean(regIncOccu) mean(revIncOccu)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);
    hb(1).FaceColor = [0.5 0.5 0.5]; hb(2).FaceColor = [0.9 0.9 0.9];
    offsetPos = [1+hb(1).XOffset 1+hb(2).XOffset 2+hb(1).XOffset 2+hb(2).XOffset];
    scatterXData = [linspace(offsetPos(1)-0.1,offsetPos(1)+0.1,length(regCorOccu)),linspace(offsetPos(2)-0.1,offsetPos(2)+0.1,length(revCorOccu)),linspace(offsetPos(3)-0.1,offsetPos(3)+0.1,length(regIncOccu)),linspace(offsetPos(4)-0.1,offsetPos(4)+0.1,length(revIncOccu))];
    hold on
    scatter(scatterXData,[regCorOccu' revCorOccu' regIncOccu' revIncOccu'],50,'ko','LineWidth',2)
    er = errorbar(offsetPos,[mean(regCorOccu) mean(revCorOccu) mean(regIncOccu) mean(revIncOccu)],[semregcorOcc semrevcorOcc semregIncorOcc semrevIncorOcc],'LineWidth',3);    
    er(1).Color = [0 0 0];                        
    er(1).LineStyle = 'none'; 
    legend([{'Congruent'},{'Incongruent'}])
    ylabel('Sipper Occupancy Time (s)')
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figSavePath filesep 'sipperOccupancy_perSession_' flags.SessionN '_Strain_' flags.Genotype],'png')

    sipperTimeStats.regCorOccu = regCorOccu;
    sipperTimeStats.regIncOccu = regIncOccu;
    sipperTimeStats.revCorOccu = revCorOccu;
    sipperTimeStats.revIncOccu = revIncOccu;

end