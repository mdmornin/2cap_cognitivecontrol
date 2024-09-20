%% Create dot plot for each animal's session-by-session history 
AnimalID = {'1-15P','1-21P','1-3P','1-13P','2-17W','2-1W','2-14W','2-5W','3-3P','3-1P','3-11P','3-2P','4-22W','4-10W','4-1W'};
sessionLabels = {'Congruent','Quinine','Incongruent','Free Access','Free Access + Quinine','High Quinine','Extinction','CS+ Reversal','CNO Injection','None'};
path = 'F:/dissDat/csvs/';
data = xlsread([path 'dotData.xlsx']);
data(isnan(data)) = 9;
figPath = 'F:/dissDat/figs';


for k = 1:min(size(data))
    for i = 1:length(data)
        s1 = plot(i,data(k,i),'ko');
        hold on
    end
end

figure('Units','normalized','Position',[0 0 1 1])
cmap = [1 0 0; 1 0.6 0.3333; 1 0.8353 0.8353;
    0.1647 0.4980 1; 0.5294 0.8039 0.8706; 0.8353 0.9647 1;
    0.2667 0.6667 0; 0.6667 0.8706 0.5294; 0.6667 1 0.6667; 0 0 0];
cmap = ["fd7f6f", "7eb0d5", "b2e061", "bd7ebe", "ffb55a", "ffee65", "beb9db", "fdcce5", "8bd3c7","000000"];
cmap = cellstr(cmap);  
cmap = hex2rgb(cmap);
cmap = cmap/255;
colormap(cmap)
s1 = imagesc(1:length(data),1:min(size(data)),data);
yticks(1:min(size(data)))
ax = gca;
yt = ax.YTick;
set(gca,'YTick',yt,'YTickLabel',AnimalID)
xlabel('Session Number')
ylabel('Animal ID')

cb = colorbar;
cb.TickLabels = sessionLabels;
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',8)
saveas(gca,[figPath filesep 'dotPlotSessionInfo'],'svg');     saveas(gca,[figPath filesep 'dotPlotSessionInfo'],'fig')
