%% Plotting variables from both congruent and incongruent sessions together
% Scratch pad to do little plots in little plot world
close all; clear all
%%
filepath = 'F:/dissDat/csvs/lat_distance_data';
cong = load([filepath 'Congruent.mat']);
incong = load([filepath 'Incongruent.mat']);
figPath = 'F:/dissDat/lrfigs';

% Clean data by excluding columns with NaN
nanidxW = find(isnan(incong.wData(2,:)));
nanidxP = find(isnan(incong.pData(2,:)));
incong.wData(:,nanidxW) = [];
incong.pData(:,nanidxP) = [];

semWDist = [std(cong.wData(2,:))/sqrt(length(cong.wData)) std(incong.wData(2,:))/sqrt(length(incong.wData))];
semWLate = [std(cong.wData(1,:))/sqrt(length(cong.wData)) std(incong.wData(1,:))/sqrt(length(incong.wData))];

semPDist = [std(cong.pData(2,:))/sqrt(length(cong.pData)) std(incong.pData(2,:))/sqrt(length(incong.pData))];
semPLate = [std(cong.pData(1,:))/sqrt(length(cong.pData)) std(incong.pData(1,:))/sqrt(length(incong.pData))];

figure('Units','normalized','Position',[0 0 1 1])
subplot(1,2,1)
hb = bar([mean(cong.wData(2,:)) mean(incong.wData(2,:)); mean(cong.pData(2,:)) mean(incong.pData(2,:))],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);

offsetPos = [1+hb(1).XOffset 1+hb(2).XOffset; 2+hb(1).XOffset 2+hb(2).XOffset];
scatterXData = [linspace(offsetPos(1)-0.1,offsetPos(1)+0.1,length(cong.wData)) linspace(offsetPos(3)-0.1,offsetPos(3)+0.1,length(incong.wData)) linspace(offsetPos(2)-0.1,offsetPos(2)+0.1,length(cong.pData)), linspace(offsetPos(4)-0.1,offsetPos(4)+0.1,length(incong.pData))];
hold on
scatter(scatterXData,[(cong.wData(2,:)) (incong.wData(2,:)) (cong.pData(2,:)) (incong.pData(2,:))],50,'ko','LineWidth',2)
er = errorbar(offsetPos,[mean(cong.wData(2,:)) mean(incong.wData(2,:)); mean(cong.pData(2,:)) mean(incong.pData(2,:))],[semWDist(1) semWDist(2); semPDist(1) semPDist(2)],'LineWidth',3);    
[er.Color] = deal([0 0 0],[0 0 0]);
[er.LineStyle] = deal('none','none');
ylabel('Distance (a.u.)')
xticklabels({'Wistar','P rat'})
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)

subplot(1,2,2)
hb = bar([mean(cong.wData(1,:)) mean(incong.wData(1,:)); mean(cong.pData(1,:)) mean(incong.pData(1,:))],'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',3);

offsetPos = [1+hb(1).XOffset 1+hb(2).XOffset; 2+hb(1).XOffset 2+hb(2).XOffset];
scatterXData = [linspace(offsetPos(1)-0.1,offsetPos(1)+0.1,length(cong.wData)) linspace(offsetPos(3)-0.1,offsetPos(3)+0.1,length(incong.wData)) linspace(offsetPos(2)-0.1,offsetPos(2)+0.1,length(cong.pData)), linspace(offsetPos(4)-0.1,offsetPos(4)+0.1,length(incong.pData))];
hold on
scatter(scatterXData,[(cong.wData(1,:)) (incong.wData(1,:)) (cong.pData(1,:)) (incong.pData(1,:))],50,'ko','LineWidth',2)
er = errorbar(offsetPos,[mean(cong.wData(1,:)) mean(incong.wData(1,:)); mean(cong.pData(1,:)) mean(incong.pData(1,:))],[semWLate(1) semWLate(2); semPLate(1) semPLate(2)],'LineWidth',3);    
[er.Color] = deal([0 0 0],[0 0 0]);
[er.LineStyle] = deal('none','none');
ylabel('Latency (s)')
xticklabels({'Wistar','P rat'})
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)

saveas(gca,[figPath filesep 'DistanceLatency_CongruentIncongruent_100ms'],'svg');     saveas(gca,[figPath filesep 'DistanceLatency_CongruentIncongruent_100ms'],'fig')
%% Put Wistar Congruent/Incongruent and P rat Congruent/Incongruent 
% on 1x2 correlation 
close all
    figure('Units','normalized','Position',[0 0 1 1]);

for i = 1:2

    if i == 1
        wDataCon = [cong.wData(2,:); cong.wData(1,:)];
        wDataInc = [incong.wData(2,:); incong.wData(1,:)];
        sessionType = 'Wistars';
        nanidxW = find(isnan(wDataCon(1,:)));
        nanidxP = find(isnan(wDataInc(1,:)));
        wDataCon(:,nanidxW) = [];
        wDataInc(:,nanidxP) = [];
        % work
        [wcorrcoef_r,wcorrcoef_p] = corr(wDataCon(1,:)',wDataCon(2,:)');
        [pcorrcoef_r,pcorrcoef_p] = corr(wDataInc(1,:)',wDataInc(2,:)');
        subplot(1,2,i)
        fitcon = polyfit(wDataCon(1,:),wDataCon(2,:),1);
        fitinc = polyfit(wDataInc(1,:),wDataInc(2,:),1);
        x1 = linspace(min(wDataCon(1,:)), max(wDataCon(1,:)), 100);
        x2 = linspace(min(wDataInc(1,:)), max(wDataInc(1,:)), 100);
        plot((wDataCon(1,:)),(wDataCon(2,:)),'ko','MarkerSize',14,'MarkerFaceColor',[0.4 0.4 1]); alpha 0.2157;
        hold on
        plot(x1,polyval(fitcon,x1),'LineWidth',5,'color',[0.4 0.4 1]); alpha 0.2157;
        plot((wDataInc(1,:)),(wDataInc(2,:)),'ko','MarkerSize',14,'MarkerFaceColor',[0.8706 0.9176 1]); alpha 0.3333;
        plot(x2,polyval(fitinc,x2),'LineWidth',5,'color',[0.8706 0.9176 1]); alpha 0.3333;
        ylim([0 4])
        xlim([0 1.6])
        xlabel('Distance (a.u.)')
        ylabel('Latency (s)')
        title(['Correlation Between Latency and Distance (' sessionType ')']);
        
        str=sprintf('Con correlation r = %1.3f',wcorrcoef_r);
        str2=sprintf('Con p-value = %1.3f',wcorrcoef_p);
        str3=sprintf('Inc correlation r = %1.3f',pcorrcoef_r);
        str4=sprintf('Inc p-value = %1.3f',pcorrcoef_p);
        
        T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
        T2 = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')-.1), str2); 
        T3 = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')-.2), str3); 
        T4 = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')-.3), str4); 
        
        set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
        set(T2, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
        set(T3, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
        set(T4, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
        
        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

    elseif i == 2

        pDataCon = [cong.pData(2,:); cong.pData(1,:)];
        pDataInc = [incong.pData(2,:); incong.pData(1,:)];
        sessionType = 'P rats';
        nanidxW = find(isnan(pDataCon(1,:)));
        nanidxP = find(isnan(pDataInc(1,:)));
        pDataCon(:,nanidxW) = [];
        pDataInc(:,nanidxP) = [];
        % work
        [wcorrcoef_r,wcorrcoef_p] = corr(pDataCon(1,:)',pDataCon(2,:)');
        [pcorrcoef_r,pcorrcoef_p] = corr(pDataInc(1,:)',pDataInc(2,:)');
        fitcon = polyfit(pDataCon(1,:),pDataCon(2,:),1);
        fitinc = polyfit(pDataInc(1,:),pDataInc(2,:),1);
        x1 = linspace(min(pDataCon(1,:)), max(pDataCon(1,:)), 100);
        x2 = linspace(min(pDataInc(1,:)), max(pDataInc(1,:)), 100);
        subplot(1,2,i)
        plot((pDataCon(1,:)),(pDataCon(2,:)),'ko','MarkerSize',14,'MarkerFaceColor',[1 0.8706 0.8706]);        
        hold on
        plot(x1,polyval(fitcon,x1),'LineWidth',5,'color',[1 0.8706 0.8706]); 
        plot((pDataInc(1,:)),(pDataInc(2,:)),'ko','MarkerSize',14,'MarkerFaceColor',[0.9 0.8745 0.8745]);
        plot(x2,polyval(fitinc,x2),'LineWidth',5,'color',[0.9 0.8745 0.8745]);
        ylim([0 4])
        xlim([0 1.6])
        xlabel('Distance (a.u.)')
        ylabel('Latency (s)')
        title(['Correlation Between Latency and Distance (' sessionType ')']);
        
        str=sprintf('Con correlation r = %1.3f',wcorrcoef_r);
        str2=sprintf('Con p-value = %1.3f',wcorrcoef_p);
        str3=sprintf('Inc correlation r = %1.3f',pcorrcoef_r);
        str4=sprintf('Inc p-value = %1.3f',pcorrcoef_p);
        
        T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
        T2 = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')-.1), str2); 
        T3 = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')-.2), str3); 
        T4 = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')-.3), str4); 
        
        set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
        set(T2, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
        set(T3, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
        set(T4, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
        
        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

    end
    % Remove NaN columns

    saveas(gca,[figPath filesep 'latency_distance_wp' '_100ms' sessionType],'svg');     saveas(gca,[figPath filesep 'latency_distance_wp' '_100ms' sessionType],'fig')

end

%% Put Wistar Congruent/Incongruent and P rat Congruent/Incongruent 
% on 1x2 correlation 
close all
filepath = 'F:/dissDat/csvs/lat_distance_data';
filepath2 = 'F:/dissDat/csvs/';

cong = load([filepath 'Congruent.mat']);
incong = load([filepath 'Incongruent.mat']);
figPath = 'F:/dissDat/lrfigs';
load([filepath2 'cueonspeeddata.mat'])
figure('Units','normalized','Position',[0 0 1 1]);

for i = 1:2

    if i == 1
        wDataCon = [cong.wData(2,:); cueOnData{1}'];
        wDataInc = [incong.wData(2,:); cueOnData{2}'];

        sessionType = 'Wistars';
        nanidxW = find(isnan(wDataCon(1,:)));
        nanidxP = find(isnan(wDataInc(1,:)));
        wDataCon(:,nanidxW) = [];
        wDataInc(:,nanidxP) = [];
%         wDataCon = zscore(wDataCon,[],2);
%         wDataInc = zscore(wDataInc,[],2); 
        % work
        [wcorrcoef_r,wcorrcoef_p] = corr(wDataCon(1,:)',wDataCon(2,:)');
        [pcorrcoef_r,pcorrcoef_p] = corr(wDataInc(1,:)',wDataInc(2,:)');
        subplot(1,2,i)
        fitcon = polyfit(wDataCon(1,:),wDataCon(2,:),1);
        fitinc = polyfit(wDataInc(1,:),wDataInc(2,:),1);
        x1 = linspace(min(wDataCon(1,:)), max(wDataCon(1,:)), 100);
        x2 = linspace(min(wDataInc(1,:)), max(wDataInc(1,:)), 100);
        plot((wDataCon(1,:)),(wDataCon(2,:)),'ko','MarkerSize',14,'MarkerFaceColor',[0.4 0.4 1]); alpha 0.2157;
        hold on
        plot(x1,polyval(fitcon,x1),'LineWidth',5,'color',[0.4 0.4 1]); alpha 0.2157;
        plot((wDataInc(1,:)),(wDataInc(2,:)),'ko','MarkerSize',14,'MarkerFaceColor',[0.8706 0.9176 1]); alpha 0.3333;
        plot(x2,polyval(fitinc,x2),'LineWidth',5,'color',[0.8706 0.9176 1]); alpha 0.3333;
%         ylim([0 4])
%         xlim([0 1.6])
        xlabel('Distance (a.u.)')
        ylabel('Speed (px/s)')
        title(['Correlation Between Speed and Distance (' sessionType ')']);
        
%         str=sprintf('Con correlation r = %1.3f',wcorrcoef_r);
%         str2=sprintf('Con p-value = %1.3f',wcorrcoef_p);
%         str3=sprintf('Inc correlation r = %1.3f',pcorrcoef_r);
%         str4=sprintf('Inc p-value = %1.3f',pcorrcoef_p);
%         
%         T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
%         T2 = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')-.1), str2); 
%         T3 = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')-.2), str3); 
%         T4 = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')-.3), str4); 
%         
%         set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
%         set(T2, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
%         set(T3, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
%         set(T4, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
        
        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

    elseif i == 2

        pDataCon = [cong.pData(2,:); cueOnData{3}'];
        pDataInc = [incong.pData(2,:); cueOnData{4}'];

        sessionType = 'P rats';
        nanidxW = find(isnan(pDataCon(1,:)));
        nanidxP = find(isnan(pDataInc(1,:)));
        pDataCon(:,nanidxW) = [];
        pDataInc(:,nanidxP) = [];
%         pDataCon = zscore(pDataCon,[],2);
%         pDataInc = zscore(pDataInc,[],2);
        % work
        [wcorrcoef_r,wcorrcoef_p] = corr(pDataCon(1,:)',pDataCon(2,:)');
        [pcorrcoef_r,pcorrcoef_p] = corr(pDataInc(1,:)',pDataInc(2,:)');
        fitcon = polyfit(pDataCon(1,:),pDataCon(2,:),1);
        fitinc = polyfit(pDataInc(1,:),pDataInc(2,:),1);
        x1 = linspace(min(pDataCon(1,:)), max(pDataCon(1,:)), 100);
        x2 = linspace(min(pDataInc(1,:)), max(pDataInc(1,:)), 100);
        subplot(1,2,i)
        plot((pDataCon(1,:)),(pDataCon(2,:)),'ko','MarkerSize',14,'MarkerFaceColor',[1 0.8706 0.8706]);        
        hold on
        plot(x1,polyval(fitcon,x1),'LineWidth',5,'color',[1 0.8706 0.8706]); 
        plot((pDataInc(1,:)),(pDataInc(2,:)),'ko','MarkerSize',14,'MarkerFaceColor',[0.9 0.8745 0.8745]);
        plot(x2,polyval(fitinc,x2),'LineWidth',5,'color',[0.9 0.8745 0.8745]);
%         ylim([0 4])
%         xlim([0 1.6])
        xlabel('Distance (a.u.)')
        ylabel('Speed (px/s)')
        title(['Correlation Between Speed and Distance (' sessionType ')']);
        
%         str=sprintf('Con correlation r = %1.3f',wcorrcoef_r);
%         str2=sprintf('Con p-value = %1.3f',wcorrcoef_p);
%         str3=sprintf('Inc correlation r = %1.3f',pcorrcoef_r);
%         str4=sprintf('Inc p-value = %1.3f',pcorrcoef_p);
%         
%         T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
%         T2 = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')-.1), str2); 
%         T3 = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')-.2), str3); 
%         T4 = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')-.3), str4); 
%         
%         set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
%         set(T2, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
%         set(T3, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
%         set(T4, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
        
        set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)

    end
    % Remove NaN columns

    saveas(gca,[figPath filesep 'speed_distance_wp' '_100ms' sessionType],'svg');     saveas(gca,[figPath filesep 'speed_distance_wp' '_100ms' sessionType],'fig')

end
%% 
tbl = readtable('speed_distance_wonly.csv');

lm = fitlm(tbl,'Speed~Distance*Condition');


%% Try other correlations like correct or incorrect approaches
% Reload data
cong = load([filepath 'Congruent.mat']);
incong = load([filepath 'Incongruent.mat']);
% Get correct and incorrect values
for i = 1:length(trlStruct)
    trlStruct(i).NumCorrect = trlStruct(i).CorrectNoCorrections;
    trlStruct(i).NumIncorrect = size(trlStruct(i).LDistInCorrectAppr,1) + size(trlStruct(i).RDistInCorrectAppr,1);
end
% Arrange data similar to how wData and pData are 
PIdxCon = startsWith(masterTbl.Strain,'P') & startsWith(masterTbl.SessionType,'Regular');
WIdxCon = startsWith(masterTbl.Strain,'W') & startsWith(masterTbl.SessionType,'Regular');
conditionLabelsCon = [{'Congruent W'}, {'Congruent W'}, {'Congruent P'}, {'Congruent P'}];


PIdxInc = startsWith(masterTbl.Strain,'P') & startsWith(masterTbl.SessionType,'Reversal');
WIdxInc = startsWith(masterTbl.Strain,'W') & startsWith(masterTbl.SessionType,'Reversal');
conditionLabelsInc = [{'Incongruent W'}, {'Incongruent W'}, {'Incongruent P'}, {'Incongruent P'}];

% Indexes
dataPosCon = [find(WIdxCon == 1); find(PIdxCon == 1)]; dataPosInc = [find(WIdxInc == 1); find(PIdxInc == 1)];

% Congruent
cong.Correct = [trlStruct(dataPosCon).NumCorrect]; cong.wCorrect = cong.Correct(1:13);

cong.pCorrect = cong.Correct(14:end); 

cong.Incorrect = [trlStruct(dataPosCon).NumIncorrect]; cong.wIncorrect = cong.Incorrect(1:13);
cong.pIncorrect= cong.Incorrect(14:end);

cong.ChngPt = [trlStruct(dataPosCon).changePoint]; cong.wChngPt = cong.ChngPt(1:13);
cong.pChngPt = cong.ChngPt(14:end);

cong.csratio = [trlStruct(dataPosCon).CSRatio]; cong.wcsratio = cong.csratio(1:13);
cong.pcsratio = cong.csratio(14:end);

% Incongruent
incong.Correct = [trlStruct(dataPosInc).NumCorrect]; incong.wCorrect = incong.Correct(1:13);
incong.pCorrect = incong.Correct(14:end);

incong.Incorrect = [trlStruct(dataPosInc).NumIncorrect]; incong.wIncorrect = incong.Incorrect(1:13);
incong.pIncorrect= incong.Incorrect(14:end);

incong.ChngPt = [trlStruct(dataPosInc).changePoint]; incong.wChngPt = incong.ChngPt(1:13);
incong.pChngPt = incong.ChngPt(14:end);

incong.csratio = [trlStruct(dataPosInc).CSRatio]; incong.wcsratio = incong.csratio(1:13);
incong.pcsratio = incong.csratio(14:end);

%% Plot and calculate Correct v. Distance 
close all
    figure('Units','normalized','Position',[0 0 1 1]);

for i = 1:2

    if i == 1
        wData = [cong.wData(2,:); cong.wCorrect];
        pData = [cong.pData(2,:); cong.pCorrect];
        sessionType = 'Congruent';
    elseif i == 2
        wData = [incong.wData(2,:); incong.wCorrect];
        pData = [incong.pData(2,:); incong.pCorrect];
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
    ylim([0 30])
    xlim([0 1.6])
    xlabel('Distance (a.u.)')
    ylabel('Correct Approaches')
    title(['Correlation Between Correct Approaches and Distance (Session Type ' sessionType ')']);
    
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
    saveas(gca,[figPath filesep 'correct_distancecorrelation' '_100ms' sessionType],'svg');     saveas(gca,[figPath filesep 'correct_distancecorrelation' '_100ms' sessionType],'fig')

end
%% Incorrect versus distance
 figure('Units','normalized','Position',[0 0 1 1]);

for i = 1:2

    if i == 1
        wData = [cong.wData(2,:); cong.wIncorrect];
        pData = [cong.pData(2,:); cong.pIncorrect];
        sessionType = 'Congruent';
    elseif i == 2
        wData = [incong.wData(2,:); incong.wIncorrect];
        pData = [incong.pData(2,:); incong.pIncorrect];
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
    xlabel('Distance (a.u.)')
    ylabel('Incorrect Approaches')
    ylim([0 30])
    xlim([0 1.6])

    title(['Correlation Between Incorrect Approaches and Distance (Session Type ' sessionType ')']);
    
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
    saveas(gca,[figPath filesep 'incorrect_distancecorrelation' '_100ms' sessionType],'svg');     saveas(gca,[figPath filesep 'incorrect_distancecorrelation' '_100ms' sessionType],'fig')

end

%% Other interesting correlations between distance: CS Discrimination, Change Point
%  After running masterAnalysisScript, trlStruct is populated with lots of
%  scalar values. Can run 
    figure('Units','normalized','Position',[0 0 1 1]);

% Incorrect versus distance
for i = 1:2

    if i == 1
        wData = [cong.wData(2,:); cong.wChngPt];
        pData = [cong.pData(2,:); cong.pChngPt];
        sessionType = 'Congruent';
    elseif i == 2
        wData = [incong.wData(2,:); incong.wChngPt];
        pData = [incong.pData(2,:); incong.pChngPt];
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
    ylim([0 30])
    xlim([0 1.6])

    xlabel('Distance (a.u.)')
    ylabel('Change Point')
    title(['Correlation Between Session Change Point and Distance (Session Type ' sessionType ')']);
    
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
    saveas(gca,[figPath filesep 'changepoint_distancecorrelation' '_100ms' sessionType],'svg');     saveas(gca,[figPath filesep 'changepoint_distancecorrelation' '_100ms' sessionType],'fig')

end
    figure('Units','normalized','Position',[0 0 1 1]);

% CS Ratio
for i = 1:2

    if i == 1
        wData = [cong.wData(2,:); cong.wcsratio];
        pData = [cong.pData(2,:); cong.pcsratio];
        sessionType = 'Congruent';
    elseif i == 2
        wData = [incong.wData(2,:); incong.wcsratio];
        pData = [incong.pData(2,:); incong.pcsratio];
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
    ylim([0 1])
    xlim([0 1.6])

    xlabel('Distance (a.u.)')
    ylabel('CS Discrimination')
    title(['Correlation Between CS Discrimination and Distance (Session Type ' sessionType ')']);
    
    str=sprintf('  Wistar correlation r = %1.2f',wcorrcoef_r);
    str2=sprintf('  Wistar p-value = %1.2f',wcorrcoef_p);
    str3=sprintf('  P rat correlation r = %1.2f',pcorrcoef_r);
    str4=sprintf('  P rat p-value = %1.2f',pcorrcoef_p);
    
    T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
    T2 = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')-0.1), str2); 
    T3 = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')-0.2), str3); 
    T4 = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')-0.3), str4); 
    
    set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    set(T2, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    set(T3, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    set(T4, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
    saveas(gca,[figPath filesep 'csratio_distancecorrelation' '_100ms' sessionType],'svg');     saveas(gca,[figPath filesep 'csratio_distancecorrelation' '_100ms' sessionType],'fig')

end


%% Plotting variables from both congruent and incongruent sessions together
% In this instance, plotting the values of the difference score between L/R
% raw firing rates indexed by their corresponding pc coefficients 

% close all; clear all
filepath = 'F:/dissDat/csvs';
cong = load([filepath filesep 'Congruentindexed_diffScores.mat']);
incong = load([filepath filesep 'Incongruentindexed_diffScores.mat']);
figPath = 'F:/dissDat/lrfigs';

addpath(genpath('F:/dissDat/restoredScripts'))

% Quick plot to see what the time course looks like 
% Organize data into 2x2 for Postive Coefficients, 2x2 for Negative
% Coefficients
% Means
appTime1 = 1:20;
appTime2 = 22:41;
dataPos = [mean(mean(cong.dWPos{1}(:,appTime1),2)) mean(mean(incong.dWPos{2}(:,appTime1),2)) mean(mean(cong.dPPos{1}(:,appTime1),2)) mean(mean(incong.dPPos{2}(:,appTime1),2)); mean(mean(cong.dWPos{1}(:,appTime2),2)) ...
    mean(mean(incong.dWPos{2}(:,appTime2),2)) mean(mean(cong.dPPos{1}(:,appTime2),2)) mean(mean(incong.dPPos{2}(:,appTime2),2))];

dataNeg = [mean(mean(cong.dWNeg{1}(:,appTime1),2)) mean(mean(incong.dWNeg{2}(:,appTime1),2)) mean(mean(cong.dPNeg{1}(:,appTime1),2)) mean(mean(incong.dPNeg{2}(:,appTime1),2)); mean(mean(cong.dWNeg{1}(:,appTime2),2)) ...
    mean(mean(incong.dWNeg{2}(:,appTime2),2)) mean(mean(cong.dPNeg{1}(:,appTime2),2)) mean(mean(incong.dPNeg{2}(:,appTime2),2))];

% SEM

dataSemPos = [std(mean(cong.dWPos{1}(:,appTime1),2))/sqrt(length(cong.dWPos{1})) std(mean(incong.dWPos{2}(:,appTime1),2))/sqrt(length(incong.dWPos{2})) ...
    std(mean(cong.dPPos{1}(:,appTime1),2))/sqrt(length(cong.dPPos{1})) std(mean(incong.dPPos{2}(:,appTime1),2))/sqrt(length(incong.dPPos{2})); std(mean(cong.dWPos{1}(:,appTime2),2))/sqrt(length(cong.dWPos{1})) ...
    std(mean(incong.dWPos{2}(:,appTime2),2))/sqrt(length(incong.dWPos{2})) std(mean(cong.dPPos{1}(:,appTime2),2))/sqrt(length(cong.dPPos{1})) std(mean(incong.dPPos{2}(:,appTime2),2))/sqrt(length(incong.dPPos{2}))];

dataSemNeg = [std(mean(cong.dWNeg{1}(:,appTime1),2))/sqrt(length(cong.dWNeg{1})) std(mean(incong.dWNeg{2}(:,appTime1),2))/sqrt(length(incong.dWNeg{2})) std(mean(cong.dPNeg{1}(:,appTime1),2))/sqrt(length(cong.dPNeg{1})) ...
    std(mean(incong.dPNeg{2}(:,appTime1),2))/sqrt(length(incong.dPNeg{2})); std(mean(cong.dWNeg{1}(:,appTime2),2))/sqrt(length(cong.dWNeg{1})) ...
    std(mean(incong.dWNeg{2}(:,appTime2),2))/sqrt(length(incong.dWNeg{2})) std(mean(cong.dPNeg{1}(:,appTime2),2))/sqrt(length(cong.dPNeg{1})) std(mean(incong.dPNeg{2}(:,appTime2),2))/sqrt(length(incong.dPNeg{2}))];



figure('Units','normalized','Position',[0 0 1 1]);
subplot(1,2,1)
hb = bar(dataPos);
hold on
offsetPos = [1+hb(1).XOffset 1+hb(2).XOffset 1+hb(3).XOffset 1+hb(4).XOffset; 2+hb(1).XOffset 2+hb(2).XOffset 2+hb(3).XOffset 2+hb(4).XOffset];
er = errorbar(offsetPos,dataPos,dataSemPos,'LineWidth',3);    
[er.Color] = deal([0 0 0],[0 0 0],[0 0 0],[0 0 0]);
[er.LineStyle] = deal('none','none','none','none');
ylim([-0.3 0.3])
xticklabels({'Pre-Approach','Post-Approach'})
legend({'Congruent Wistar', 'Incongruent Wistar', 'Congruent P rat', 'Incongruent P rat'})
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)

subplot(1,2,2)
hb = bar(dataNeg);
hold on
offsetPos = [1+hb(1).XOffset 1+hb(2).XOffset 1+hb(3).XOffset 1+hb(4).XOffset; 2+hb(1).XOffset 2+hb(2).XOffset 2+hb(3).XOffset 2+hb(4).XOffset];
er = errorbar(offsetPos,dataNeg,dataSemNeg,'LineWidth',3);    
[er.Color] = deal([0 0 0],[0 0 0],[0 0 0],[0 0 0]);
[er.LineStyle] = deal('none','none','none','none');
ylim([-0.3 0.3])
xticklabels({'Pre-Approach','Post-Approach'})
legend({'Congruent Wistar', 'Incongruent Wistar', 'Congruent P rat', 'Incongruent P rat'})
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',4)
saveas(gca,[figPath filesep 'prepostbar' '_100ms'],'svg');     saveas(gca,[figPath filesep 'prepostbar' '_100ms'],'fig')

% Need to now run stats

dataPos = [(mean(cong.dWPos{1}(:,1:20))); (mean(incong.dWPos{2}(:,1:20))); (mean(cong.dPPos{1}(:,1:20))); (mean(incong.dPPos{2}(:,1:20))); (mean(cong.dWPos{1}(:,22:end))); ...
    (mean(incong.dWPos{2}(:,22:end))); (mean(cong.dPPos{1}(:,22:end))); (mean(incong.dPPos{2}(:,22:end)))];

dataNeg = [(mean(cong.dWNeg{1}(:,1:20))); (mean(incong.dWNeg{2}(:,1:20))); (mean(cong.dPNeg{1}(:,1:20))); (mean(incong.dPNeg{2}(:,1:20))); (mean(cong.dWNeg{1}(:,22:end))); ...
    (mean(incong.dWNeg{2}(:,22:end))); (mean(cong.dPNeg{1}(:,22:end))); (mean(incong.dPNeg{2}(:,22:end)))];

% First test again zero
for i = 1:size(dataPos,1)
    [h(i),p(i),stat] = ttest(dataPos(i,:),0);
end
p_pos = fdr_bh(p);

for i = 1:size(dataNeg,1)
    [h(i),p(i),stat] = ttest(dataNeg(i,:),0);

end
p_neg = fdr_bh(p);

% 2x2 RANOVA?



