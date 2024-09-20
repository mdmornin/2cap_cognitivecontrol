filepath = 'F:/dissDat/csvs/';
cong = load([filepath 'Congruent_distancedata.mat']);
incong = load([filepath 'Incongruent_distancedata.mat']);
figPath = 'F:/dissDat/lrfigs';

wCon = cong.meanPermValW; 
wInc = incong.meanPermValW;
pCon = cong.meanPermValP;
pInc = incong.meanPermValP;

labels = [repmat({'WistarCongruent'},1,length(wCon)) repmat({'WistarIncongruent'},1,length(wInc)) repmat({'PCongruent'},1,length(pCon)) repmat({'PIncongruent'},1,length(pInc))];
data = [wCon wInc pCon pInc];
[p, table, stats] = kruskalwallis(data,labels)
multcompare(stats)

figure('Units','normalized','Position',[0 0 1 1])
boxplot(data,labels,'notch','on')
ylabel('Distance (a.u.)')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',30,'FontWeight','bold','LineWidth',8)
saveas(gca,[figPath filesep 'DistanceEachSessionGenotype'],'svg');     saveas(gca,[figPath filesep 'DistanceEachSessionGenotype'],'fig')
