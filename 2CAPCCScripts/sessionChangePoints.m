function changebySession = sessionChangePoints(changebySession,trlStruct,figSavePath,masterTbl,flags)


    if strcmp(flags.SessionN,'1') && strcmp(flags.Genotype,'all')    
        regIdx = strcmp(masterTbl.SessionType,'Regular');
        revIdx = strcmp(masterTbl.SessionType,'Reversal');

        regAppCP = [trlStruct(regIdx).changePoint];
        revAppCP = [trlStruct(revIdx).changePoint];

        regCorrectionCP = [trlStruct(regIdx).changePointCorrections];
        revCorrectionCP = [trlStruct(revIdx).changePointCorrections];

    elseif strcmp(flags.SessionN,'2') && strcmp(flags.Genotype,'all')
        regIdx = strcmp(masterTbl.SessionType,'Regular2'); 
        revIdx = strcmp(masterTbl.SessionType,'Reversal2');

        regAppCP = [trlStruct(regIdx).changePoint];
        revAppCP = [trlStruct(revIdx).changePoint];
        regCorrectionCP = [trlStruct(regIdx).changePointCorrections];
        revCorrectionCP = [trlStruct(revIdx).changePointCorrections];

    elseif strcmp(flags.SessionN,'3') && strcmp(flags.Genotype,'all')
        regIdx = strcmp(masterTbl.SessionType,'Regular3');
        revIdx = strcmp(masterTbl.SessionType,'Reversal3');

        regAppCP = [trlStruct(regIdx).changePoint];
        revAppCP = [trlStruct(revIdx).changePoint];
        regCorrectionCP = [trlStruct(regIdx).changePointCorrections];
        revCorrectionCP = [trlStruct(revIdx).changePointCorrections];
    end
    
    % Plot Approach Change Points
    semregChPt = std(regAppCP)/sqrt(length(regAppCP));
    semrevChPt = std(revAppCP)/sqrt(length(revAppCP));
    
    figure('Units','normalized','Position',[0 0 1 1])
    bar(categorical({'Congruent Change Point','Incongruent Change Point'}),[mean(regAppCP) mean(revAppCP)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3)
    hold on
    plot([linspace(0.95,1.05,length(regAppCP));linspace(1.95,2.05,length(revAppCP))],[regAppCP; revAppCP],'k--o','LineWidth',3,'MarkerSize',10)
    er = errorbar(categorical({'Congruent Change Point','Incongruent Change Point'}),[mean(regAppCP) mean(revAppCP)],[semregChPt semrevChPt],'LineWidth',3);    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
    ylabel('Trials')
    ylim([0 48])
    
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figSavePath filesep 'approachChangePoints_Sessions_' flags.SessionN '_Strain_' flags.Genotype],flags.figextension)

    % Plot Correction Change Points
    semregChPtcorrections = std(regCorrectionCP)/sqrt(length(regCorrectionCP));
    semrevChPtcorrections = std(revCorrectionCP)/sqrt(length(revCorrectionCP));
    
    figure('Units','normalized','Position',[0 0 1 1])
    bar(categorical({'Congruent Change Point','Incongruent Change Point'}),[mean(regCorrectionCP) mean(revCorrectionCP)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3)
    hold on
    plot([linspace(0.95,1.05,length(regCorrectionCP));linspace(1.95,2.05,length(revCorrectionCP))],[regCorrectionCP; revCorrectionCP],'k--o','LineWidth',3,'MarkerSize',10)
    er = errorbar(categorical({'Congruent Change Point','Incongruent Change Point'}),[mean(regCorrectionCP) mean(revCorrectionCP)],[semregChPtcorrections semrevChPtcorrections],'LineWidth',3);    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
    ylabel('Trials')
    ylim([0 48])
    
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figSavePath filesep 'correctionChangePoints_Sessions_' flags.SessionN '_Strain_' flags.Genotype],flags.figextension)

    % Output data

    if strcmp(flags.SessionN,'1') && strcmp(flags.Genotype,'all')
        changebySession.session1_regappcor = regAppCP;
        changebySession.session1_revappcor = revAppCP;
    elseif strcmp(flags.SessionN,'2') && strcmp(flags.Genotype,'all')
        changebySession.session2_regappcor = regAppCP;
        changebySession.session2_revappcor = revAppCP;
    elseif strcmp(flags.SessionN,'3') && strcmp(flags.Genotype,'all')
        changebySession.session3_regappcor = regAppCP;
        changebySession.session3_revappcor = revAppCP;
    end

end