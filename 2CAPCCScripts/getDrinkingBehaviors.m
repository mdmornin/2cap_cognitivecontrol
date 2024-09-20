function [drinkingStats,trlStruct] = getDrinkingBehaviors(trlStruct,figSavePath,masterTbl,flags)
    % Function extracts and compiles drinking behavior.
    % Overall intake per session in plotted. Overall rate of intake per
    % session is plotted. We infer rate of intake from sipper occupancy
    % time. The sipper occupancy time is computed in
    % getSipperTimePerSession. Overall intake is computed from the
    % difference in bottle weights pre-and-post session divided by animals
    % weight (g/kg).

    %% Intake & Rate of Intake. Generate data according to flags. 
    % Pull data
    % Intake
    regPIdx = startsWith(masterTbl.SessionType,'Regular') & strcmp(masterTbl.Strain,'P');
    revPIdx = startsWith(masterTbl.SessionType,'Reversal') & strcmp(masterTbl.Strain,'P');
    
    regWIdx = startsWith(masterTbl.SessionType,'Regular') & strcmp(masterTbl.Strain,'W');
    revWIdx = startsWith(masterTbl.SessionType,'Reversal') & strcmp(masterTbl.Strain,'W');
    
    intakePReg = masterTbl.Intake(regPIdx);
    intakePRev = masterTbl.Intake(revPIdx);
    
    intakeWReg = masterTbl.Intake(regWIdx);
    intakeWRev = masterTbl.Intake(revWIdx);
    

    % Rate of Intake
    for i = 1:length(trlStruct)
        trlStruct(i).totalCorrectSipperTime = sum(trlStruct(i).CorCorOcc)./30;
    end
    
    % Divide total intake by total sipper time to find rate (g/kg/s)
    for i = 1:length(trlStruct)
      trlStruct(i).intakeRate = masterTbl.Intake(i) / trlStruct(i).totalCorrectSipperTime;
    end
    

    intakeRatePReg = [trlStruct(regPIdx).intakeRate];
    intakeRatePRev = [trlStruct(revPIdx).intakeRate];
    
    intakeRateWReg = [trlStruct(regWIdx).intakeRate];
    intakeRateWRev = [trlStruct(revWIdx).intakeRate];

    % Plotting intakes
    % Sem Intakes
    semregcorLat = std(intakePReg)/sqrt(length(intakePReg));
    semrevcorLat = std(intakePRev)/sqrt(length(intakePRev));
    semregIncorLat = std(intakeWReg)/sqrt(length(intakeWReg));
    semrevIncorLat = std(intakeWRev)/sqrt(length(intakeWRev));

    % Plot figure
    figure('Units','normalized','Position',[0 0 1 1])
    subplot(1,2,1)
    hb = bar(categorical({'P Rats','Wistars'}),[mean(intakePReg) mean(intakePRev); mean(intakeWReg) mean(intakeWRev)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);
    hb(1).FaceColor = [0.5 0.5 0.5]; hb(2).FaceColor = [0.9 0.9 0.9];
    offsetPos = [1+hb(1).XOffset 1+hb(2).XOffset 2+hb(1).XOffset 2+hb(2).XOffset];
    scatterXData = [linspace(offsetPos(1)-0.1,offsetPos(1)+0.1,length(intakePReg)),linspace(offsetPos(2)-0.1,offsetPos(2)+0.1,length(intakePRev)),linspace(offsetPos(3)-0.1,offsetPos(3)+0.1,length(intakeWReg)),linspace(offsetPos(4)-0.1,offsetPos(4)+0.1,length(intakeWRev))];
    hold on
    plot([linspace(offsetPos(1)-0.05,offsetPos(1)+0.05,length(intakePReg)); linspace(offsetPos(2)-0.05,offsetPos(2)+0.05,length(intakePRev))],[intakePReg intakePRev]','ko','LineWidth',2,'MarkerSize',8)
    plot([linspace(offsetPos(3)-0.05,offsetPos(3)+0.05,length(intakeWReg)); linspace(offsetPos(4)-0.05,offsetPos(4)+0.05,length(intakeWRev))],[intakeWReg intakeWRev]','ko','LineWidth',2,'MarkerSize',8)
    
    er = errorbar(offsetPos,[mean(intakePReg) mean(intakePRev) mean(intakeWReg) mean(intakeWRev)],[semregcorLat semrevcorLat semregIncorLat semrevIncorLat],'LineWidth',3);    
    er(1).Color = [0 0 0];                        
    er(1).LineStyle = 'none'; 
    legend([{'Congruent'},{'Incongruent'}])
    ylabel('Intake (g/kg)')
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
    % Save figure
    saveas(gca,[figSavePath filesep 'intake_'],'svg')

    % Plotting rate of intakes
    semintakeregcorLat = std(intakeRatePReg)/sqrt(length(intakeRatePReg));
    semintakerevcorLat = std(intakeRatePRev)/sqrt(length(intakeRatePRev));
    semintakeregIncorLat = std(intakeRateWReg)/sqrt(length(intakeRateWReg));
    semintakerevIncorLat = std(intakeRateWRev)/sqrt(length(intakeRateWRev));

    % Plot figure
    subplot(1,2,2)
    hb2 = bar(categorical({'P Rats','Wistars'}),[mean(intakeRatePReg) mean(intakeRatePRev); mean(intakeRateWReg) mean(intakeRateWRev)],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);
    hb2(1).FaceColor = [0.5 0.5 0.5]; hb2(2).FaceColor = [0.9 0.9 0.9];
    offsetPos = [1+hb2(1).XOffset 1+hb2(2).XOffset 2+hb2(1).XOffset 2+hb2(2).XOffset];
    hold on
    plot([linspace(offsetPos(1)-0.05,offsetPos(1)+0.05,length(intakeRatePReg)); linspace(offsetPos(2)-0.05,offsetPos(2)+0.05,length(intakeRatePRev))],[intakeRatePReg; intakeRatePRev],'ko','LineWidth',2,'MarkerSize',8)
    plot([linspace(offsetPos(3)-0.05,offsetPos(3)+0.05,length(intakeRateWReg)); linspace(offsetPos(4)-0.05,offsetPos(4)+0.05,length(intakeRateWRev))],[intakeRateWReg; intakeRateWRev],'ko','LineWidth',2,'MarkerSize',8)
    
    er = errorbar(offsetPos,[mean(intakeRatePReg) mean(intakeRatePRev) mean(intakeRateWReg) mean(intakeRateWRev)],[semintakeregcorLat semintakerevcorLat semintakeregIncorLat semintakerevIncorLat],'LineWidth',3);    
    er(1).Color = [0 0 0];                        
    er(1).LineStyle = 'none'; 
    legend([{'Congruent'},{'Incongruent'}])
    ylabel('Rate of Intake (g/kg/s)')
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
    % Save figure
    saveas(gca,[figSavePath filesep 'intake_intakeRate_'],'svg')

    drinkingStats.regPintake = intakePReg; drinkingStats.regPintakeRate = intakeRatePReg;
    drinkingStats.revPintake = intakePRev; drinkingStats.revPintakeRate = intakeRatePRev;
    drinkingStats.regWintake = intakeWReg; drinkingStats.regWintakeRate = intakeRateWReg;
    drinkingStats.revWintake = intakeWRev; drinkingStats.revWintakeRate = intakeRateWRev;

end