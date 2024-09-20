%% Script for correlations of behavioral data in 2CAP
% Necessary information from the behavior is in the masterAnalysis script,
% run that here
% Add path if needed
% Find a way to skip if data exists
if ~isfield(trlStruct,'CSRatio')
    run('analysisMasterScript.m')
else
    disp('Data is already loaded')
end

%% Make indexing variables
% Arrange data similar to how wData and pData are 
regPIdx = startsWith(masterTbl.SessionType,'Regular') & strcmp(masterTbl.Strain,'P');
revPIdx = startsWith(masterTbl.SessionType,'Reversal') & strcmp(masterTbl.Strain,'P');
    
regWIdx = startsWith(masterTbl.SessionType,'Regular') & strcmp(masterTbl.Strain,'W');
revWIdx = startsWith(masterTbl.SessionType,'Reversal') & strcmp(masterTbl.Strain,'W');

% Pull all data desired
corPCon = [trlStruct(regPIdx).CorrectNoCorrections];     incorPCon = [trlStruct(regPIdx).Incorrect];       omitPCon = [trlStruct(regPIdx).Omissions];       correctionsPCon = [trlStruct(regPIdx).Corrections];
corPInc = [trlStruct(revPIdx).CorrectNoCorrections];     incorPInc = [trlStruct(revPIdx).Incorrect];       omitPInc = [trlStruct(revPIdx).Omissions];       correctionsPInc = [trlStruct(revPIdx).Corrections];
corWCon = [trlStruct(regWIdx).CorrectNoCorrections];     incorWCon = [trlStruct(regWIdx).Incorrect];       omitWCon = [trlStruct(regWIdx).Omissions];       correctionsWCon = [trlStruct(regWIdx).Corrections];
corWInc = [trlStruct(revWIdx).CorrectNoCorrections];     incorWInc = [trlStruct(revWIdx).Incorrect];       omitWInc = [trlStruct(revWIdx).Omissions];       correctionsWInc = [trlStruct(revWIdx).Corrections];

for i = 1:length(trlStruct)
    trlStruct(i).meanLat = mean(trlStruct(i).CorrectLatencies);
end
csRatPCon = [trlStruct(regPIdx).CSRatio];                latPCon = ([trlStruct(regPIdx).meanLat]);
csRatPInc = [trlStruct(revPIdx).CSRatio];                latPInc = ([trlStruct(revPIdx).meanLat]);
csRatWCon = [trlStruct(regWIdx).CSRatio];                latWCon = ([trlStruct(regWIdx).meanLat]);
csRatWInc = [trlStruct(revWIdx).CSRatio];                latWInc = ([trlStruct(revWIdx).meanLat]);

intakePReg = masterTbl.Intake(regPIdx); intakeWReg = masterTbl.Intake(regWIdx); intakeWRev = masterTbl.Intake(revWIdx); intakePRev = masterTbl.Intake(revPIdx);



% Calculate SEM values

semcorPCon = std(corPCon)/sqrt(length(corPCon));    semincorPCon = std(incorPCon)/sqrt(length(incorPCon));  semomitPCon = std(omitPCon)/sqrt(length(omitPCon));  semcorrectionsPCon = std(correctionsPCon)/sqrt(length(correctionsPCon));  
semcorPInc = std(corPInc)/sqrt(length(corPInc));    semincorPInc = std(incorPInc)/sqrt(length(incorPInc));  semomitPInc = std(omitPInc)/sqrt(length(omitPInc));  semcorrectionsPInc = std(correctionsPInc)/sqrt(length(correctionsPInc));
semcorWCon = std(corWCon)/sqrt(length(corWCon));    semincorWCon = std(incorWCon)/sqrt(length(incorWCon));  semomitWCon = std(omitWCon)/sqrt(length(omitWCon));  semcorrectionsWCon = std(correctionsWCon)/sqrt(length(correctionsWCon));
semcorWInc = std(corWInc)/sqrt(length(corWInc));    semincorWInc = std(incorWInc)/sqrt(length(incorWInc));  semomitWInc = std(omitWInc)/sqrt(length(omitWInc));  semcorrectionsWInc = std(correctionsWInc)/sqrt(length(correctionsWInc));

%%
figPath = 'F:/dissDat/behavioralCorrelationFigs';
close all
    figure('Units','normalized','Position',[0 0 1 1]);

for i = 1:2

    if i == 1
        wData = [latWCon; corWCon];
        pData = [latPCon; corPCon];
        sessionType = 'Congruent';
    elseif i == 2
        wData = [latWInc; corWInc];
        pData = [latPInc; corPInc];
        sessionType = 'Incongruent';
    end
    % Remove NaN columns
    nanidxW = find(isnan(wData(1,:)));
    nanidxP = find(isnan(pData(1,:)));
    wData(:,nanidxW) = [];
    pData(:,nanidxP) = [];
    % work
    [wcorrcoef_r,wcorrcoef_p] = corr(wData(1,:)',wData(2,:)');
    [pcorrcoef_r,pcorrcoef_p] = corr(pData(1,:)',pData(2,:)');
    subplot(1,2,i)
    plot((wData(2,:)),(wData(1,:)),'bo','MarkerSize',14,'MarkerFaceColor','b');
    hold on
    plot((pData(2,:)),(pData(1,:)),'ro','MarkerSize',14,'MarkerFaceColor','r');
    xlim([0 25])
    ylim([0 4])
    xlabel('Correct Approaches')
    ylabel('Latency')
    title(['Correlation Between Correct Approaches and Latency (s)']);
    
    str=sprintf('  Wistar correlation r = %1.2f',wcorrcoef_r);
    str2=sprintf('  Wistar p-value = %1.2f',wcorrcoef_p);
    str3=sprintf('  P rat correlation r = %1.2f',pcorrcoef_r);
    str4=sprintf('  P rat p-value = %1.2f',pcorrcoef_p);
    
    T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
    T2 = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')-1), str2); 
    T3 = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')-2), str3); 
    T4 = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')-3), str4); 
    
    set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    set(T2, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    set(T3, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    set(T4, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figPath filesep 'correct_latencycorrelation' '_100ms' sessionType],'svg');     saveas(gca,[figPath filesep 'correct_latencycorrelation' '_100ms' sessionType],'fig')

end

%%
close all
    figure('Units','normalized','Position',[0 0 1 1]);

for i = 1:2

    if i == 1
        wData = [latWCon; csRatWCon];
        pData = [latPCon; csRatPCon];
        sessionType = 'Congruent';
    elseif i == 2
        wData = [latWInc; csRatWInc];
        pData = [latPInc; csRatPInc];
        sessionType = 'Incongruent';
    end
    % Remove NaN columns
    nanidxW = find(isnan(wData(1,:)));
    nanidxP = find(isnan(pData(1,:)));
    wData(:,nanidxW) = [];
    pData(:,nanidxP) = [];
    % work
    [wcorrcoef_r,wcorrcoef_p] = corr(wData(1,:)',wData(2,:)');
    [pcorrcoef_r,pcorrcoef_p] = corr(pData(1,:)',pData(2,:)');
    subplot(1,2,i)
    plot((wData(1,:)),(wData(2,:)),'bo','MarkerSize',14,'MarkerFaceColor','b');
    hold on
    plot((pData(1,:)),(pData(2,:)),'ro','MarkerSize',14,'MarkerFaceColor','r');
    xlim([0 4])
    ylim([0 1])
    ylabel('CS Ratio')
    xlabel('Latency')
    title(['Correlation Between CS Ratio and Latency (s)']);
    
    str=sprintf('  Wistar correlation r = %1.2f',wcorrcoef_r);
    str2=sprintf('  Wistar p-value = %1.2f',wcorrcoef_p);
    str3=sprintf('  P rat correlation r = %1.2f',pcorrcoef_r);
    str4=sprintf('  P rat p-value = %1.2f',pcorrcoef_p);
    
    T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
    T2 = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')-.11), str2); 
    T3 = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')-.2), str3); 
    T4 = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')-.3), str4); 
    
    set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    set(T2, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    set(T3, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    set(T4, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figPath filesep 'csratio_latencycorrelation' '_100ms' sessionType],'svg');     saveas(gca,[figPath filesep 'csratio_latencycorrelation' '_100ms' sessionType],'fig')

end

%%
figPath = 'F:/dissDat/behavioralCorrelationFigs';
close all
    figure('Units','normalized','Position',[0 0 1 1]);

for i = 1:2

    if i == 1
        wData = [latWCon; incorWCon];
        pData = [latPCon; incorPCon];
        sessionType = 'Congruent';
    elseif i == 2
        wData = [latWInc; incorWInc];
        pData = [latPInc; incorPInc];
        sessionType = 'Incongruent';
    end
    % Remove NaN columns
    nanidxW = find(isnan(wData(1,:)));
    nanidxP = find(isnan(pData(1,:)));
    wData(:,nanidxW) = [];
    pData(:,nanidxP) = [];
    % work
    [wcorrcoef_r,wcorrcoef_p] = corr(wData(1,:)',wData(2,:)');
    [pcorrcoef_r,pcorrcoef_p] = corr(pData(1,:)',pData(2,:)');
    subplot(1,2,i)
    plot((wData(2,:)),(wData(1,:)),'bo','MarkerSize',14,'MarkerFaceColor','b');
    hold on
    plot((pData(2,:)),(pData(1,:)),'ro','MarkerSize',14,'MarkerFaceColor','r');
    xlim([0 30])
    ylim([0 4])
    xlabel('Incorrect Approaches')
    ylabel('Latency')
    title(['Correlation Between Incorrect Approaches and Latency (s)']);
    
    str=sprintf('  Wistar correlation r = %1.2f',wcorrcoef_r);
    str2=sprintf('  Wistar p-value = %1.2f',wcorrcoef_p);
    str3=sprintf('  P rat correlation r = %1.2f',pcorrcoef_r);
    str4=sprintf('  P rat p-value = %1.2f',pcorrcoef_p);
    
    T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
    T2 = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')-1), str2); 
    T3 = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')-2), str3); 
    T4 = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')-3), str4); 
    
    set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    set(T2, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    set(T3, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    set(T4, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figPath filesep 'incorrect_latencycorrelation' '_100ms' sessionType],'svg');     saveas(gca,[figPath filesep 'incorrect_latencycorrelation' '_100ms' sessionType],'fig')

end

%%
figPath = 'F:/dissDat/behavioralCorrelationFigs';
close all
    figure('Units','normalized','Position',[0 0 1 1]);

for i = 1:2

    if i == 1
        wData = [csRatWCon  ; incorWCon];
        pData = [csRatPCon; incorPCon];
        sessionType = 'Congruent';
    elseif i == 2
        wData = [csRatWInc; incorWInc];
        pData = [csRatPInc; incorPInc];
        sessionType = 'Incongruent';
    end
    % Remove NaN columns
    nanidxW = find(isnan(wData(1,:)));
    nanidxP = find(isnan(pData(1,:)));
    wData(:,nanidxW) = [];
    pData(:,nanidxP) = [];
    % work
    [wcorrcoef_r,wcorrcoef_p] = corr(wData(1,:)',wData(2,:)');
    [pcorrcoef_r,pcorrcoef_p] = corr(pData(1,:)',pData(2,:)');
    subplot(1,2,i)
    plot((wData(2,:)),(wData(1,:)),'bo','MarkerSize',14,'MarkerFaceColor','b');
    hold on
    plot((pData(2,:)),(pData(1,:)),'ro','MarkerSize',14,'MarkerFaceColor','r');
    xlim([0 30])
    ylim([0 1])
    xlabel('Incorrect Approaches')
    ylabel('CS Ratio')
    title(['Correlation Between Incorrect Approaches and CSRatio']);
    
    str=sprintf('  Wistar correlation r = %1.2f',wcorrcoef_r);
    str2=sprintf('  Wistar p-value = %1.2f',wcorrcoef_p);
    str3=sprintf('  P rat correlation r = %1.2f',pcorrcoef_r);
    str4=sprintf('  P rat p-value = %1.2f',pcorrcoef_p);
    
    T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
    T2 = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')-.1), str2); 
    T3 = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')-.2), str3); 
    T4 = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')-.3), str4); 
    
    set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    set(T2, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    set(T3, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    set(T4, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figPath filesep 'incorrect_csratiocorrelation' '_100ms' sessionType],'svg');     saveas(gca,[figPath filesep 'incorrect_csratiocorrelation' '_100ms' sessionType],'fig')

end

%%
%%
figPath = 'F:/dissDat/behavioralCorrelationFigs';
close all
    figure('Units','normalized','Position',[0 0 1 1]);

for i = 1:2

    if i == 1
        wData = [csRatWCon  ; corWCon];
        pData = [csRatPCon; corPCon];
        sessionType = 'Congruent';
    elseif i == 2
        wData = [csRatWInc; corWInc];
        pData = [csRatPInc; corPInc];
        sessionType = 'Incongruent';
    end
    % Remove NaN columns
    nanidxW = find(isnan(wData(1,:)));
    nanidxP = find(isnan(pData(1,:)));
    wData(:,nanidxW) = [];
    pData(:,nanidxP) = [];
    % work
    [wcorrcoef_r,wcorrcoef_p] = corr(wData(1,:)',wData(2,:)');
    [pcorrcoef_r,pcorrcoef_p] = corr(pData(1,:)',pData(2,:)');
    subplot(1,2,i)
    plot((wData(2,:)),(wData(1,:)),'bo','MarkerSize',14,'MarkerFaceColor','b');
    hold on
    plot((pData(2,:)),(pData(1,:)),'ro','MarkerSize',14,'MarkerFaceColor','r');
    xlim([0 25])
    ylim([0 1])
    xlabel('Correct Approaches')
    ylabel('CS Ratio')
    title(['Correlation Between Correct Approaches and CSRatio']);
    
    str=sprintf('  Wistar correlation r = %1.2f',wcorrcoef_r);
    str2=sprintf('  Wistar p-value = %1.2f',wcorrcoef_p);
    str3=sprintf('  P rat correlation r = %1.2f',pcorrcoef_r);
    str4=sprintf('  P rat p-value = %1.2f',pcorrcoef_p);
    
    T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
    T2 = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')-.1), str2); 
    T3 = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')-.2), str3); 
    T4 = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')-.3), str4); 
    
    set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    set(T2, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    set(T3, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    set(T4, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figPath filesep 'correct_csratiocorrelation' '_100ms' sessionType],'svg');     saveas(gca,[figPath filesep 'correct_csratiocorrelation' '_100ms' sessionType],'fig')

end

%% Intake v Correct
figPath = 'F:/dissDat/behavioralCorrelationFigs';
close all
    figure('Units','normalized','Position',[0 0 1 1]);

for i = 1:2

    if i == 1
        wData = [intakeWReg'  ; corWCon];
        pData = [intakePReg'; corPCon];
        sessionType = 'Congruent';
    elseif i == 2
        wData = [intakeWRev'; corWInc];
        pData = [intakePRev'; corPInc];
        sessionType = 'Incongruent';
    end
    % Remove NaN columns
    nanidxW = find(isnan(wData(1,:)));
    nanidxP = find(isnan(pData(1,:)));
    wData(:,nanidxW) = [];
    pData(:,nanidxP) = [];
    % work
    [wcorrcoef_r,wcorrcoef_p] = corr(wData(1,:)',wData(2,:)');
    [pcorrcoef_r,pcorrcoef_p] = corr(pData(1,:)',pData(2,:)');
    subplot(1,2,i)
    plot((wData(2,:)),(wData(1,:)),'bo','MarkerSize',14,'MarkerFaceColor','b');
    hold on
    plot((pData(2,:)),(pData(1,:)),'ro','MarkerSize',14,'MarkerFaceColor','r');
    xlim([0 25])
    ylim([0 3])
    xlabel('Correct Approaches')
    ylabel('Intake (g/kg)')
    title(['Correlation Between Intake and Correct Approaches']);
    
    str=sprintf('  Wistar correlation r = %1.2f',wcorrcoef_r);
    str2=sprintf('  Wistar p-value = %1.2f',wcorrcoef_p);
    str3=sprintf('  P rat correlation r = %1.2f',pcorrcoef_r);
    str4=sprintf('  P rat p-value = %1.2f',pcorrcoef_p);
    
    T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
    T2 = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')-.1), str2); 
    T3 = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')-.2), str3); 
    T4 = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')-.3), str4); 
    
    set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    set(T2, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    set(T3, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    set(T4, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figPath filesep 'correct_intakecorrelation' '_100ms' sessionType],'svg');     saveas(gca,[figPath filesep 'correct_intakecorrelation' '_100ms' sessionType],'fig')

end
%% Intake v CSRatio
figPath = 'F:/dissDat/behavioralCorrelationFigs';
close all
    figure('Units','normalized','Position',[0 0 1 1]);

for i = 1:2

    if i == 1
        wData = [intakeWReg'  ; csRatWCon];
        pData = [intakePReg'; csRatPCon];
        sessionType = 'Congruent';
    elseif i == 2
        wData = [intakeWRev'; csRatWInc];
        pData = [intakePRev'; csRatPInc];
        sessionType = 'Incongruent';
    end
    % Remove NaN columns
    nanidxW = find(isnan(wData(1,:)));
    nanidxP = find(isnan(pData(1,:)));
    wData(:,nanidxW) = [];
    pData(:,nanidxP) = [];
    % work
    [wcorrcoef_r,wcorrcoef_p] = corr(wData(1,:)',wData(2,:)');
    [pcorrcoef_r,pcorrcoef_p] = corr(pData(1,:)',pData(2,:)');
    subplot(1,2,i)
    plot((wData(1,:)),(wData(2,:)),'bo','MarkerSize',14,'MarkerFaceColor','b');
    hold on
    plot((pData(1,:)),(pData(2,:)),'ro','MarkerSize',14,'MarkerFaceColor','r');
    xlim([0 3])
    ylim([0 1])
    xlabel('Intake (g/kg)')
    ylabel('CS Ratio')
    title(['Correlation Between Intake and CS Ratio Approaches']);
    
    str=sprintf('  Wistar correlation r = %1.2f',wcorrcoef_r);
    str2=sprintf('  Wistar p-value = %1.2f',wcorrcoef_p);
    str3=sprintf('  P rat correlation r = %1.2f',pcorrcoef_r);
    str4=sprintf('  P rat p-value = %1.2f',pcorrcoef_p);
    
    T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
    T2 = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')-.1), str2); 
    T3 = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')-.2), str3); 
    T4 = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')-.3), str4); 
    
    set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    set(T2, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    set(T3, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    set(T4, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figPath filesep 'csratio_intakecorrelation' '_100ms' sessionType],'svg');     saveas(gca,[figPath filesep 'csratio_intakecorrelation' '_100ms' sessionType],'fig')

end

%% Intake v Latency
figPath = 'F:/dissDat/behavioralCorrelationFigs';
close all
    figure('Units','normalized','Position',[0 0 1 1]);

for i = 1:2

    if i == 1
        wData = [intakeWReg'  ; latWCon];
        pData = [intakePReg'; latPCon];
        sessionType = 'Congruent';
    elseif i == 2
        wData = [intakeWRev'; latWInc];
        pData = [intakePRev'; latPInc];
        sessionType = 'Incongruent';
    end
    % Remove NaN columns
    nanidxW = find(isnan(wData(1,:)));
    nanidxP = find(isnan(pData(1,:)));
    wData(:,nanidxW) = [];
    pData(:,nanidxP) = [];
    % work
    [wcorrcoef_r,wcorrcoef_p] = corr(wData(1,:)',wData(2,:)');
    [pcorrcoef_r,pcorrcoef_p] = corr(pData(1,:)',pData(2,:)');
    subplot(1,2,i)
    plot((wData(1,:)),(wData(2,:)),'bo','MarkerSize',14,'MarkerFaceColor','b');
    hold on
    plot((pData(1,:)),(pData(2,:)),'ro','MarkerSize',14,'MarkerFaceColor','r');
    xlim([0 3])
    ylim([0 5])
    xlabel('Intake (g/kg)')
    ylabel('Latency')
    title(['Correlation Between Intake and Latency Approaches']);
    
    str=sprintf('  Wistar correlation r = %1.2f',wcorrcoef_r);
    str2=sprintf('  Wistar p-value = %1.2f',wcorrcoef_p);
    str3=sprintf('  P rat correlation r = %1.2f',pcorrcoef_r);
    str4=sprintf('  P rat p-value = %1.2f',pcorrcoef_p);
    
    T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
    T2 = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')-.1), str2); 
    T3 = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')-.2), str3); 
    T4 = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')-.3), str4); 
    
    set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    set(T2, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    set(T3, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    set(T4, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figPath filesep 'lat_intakecorrelation' '_100ms' sessionType],'svg');     saveas(gca,[figPath filesep 'lat_intakecorrelation' '_100ms' sessionType],'fig')

end