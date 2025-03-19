%%
% Fig.2(c)
% creator: WxyZ
% Date: 20250314

%%
ft_defaults;

load('.\mat\SensorLevel\GFP_stat_ICA.mat')
load('.\mat\SensorLevel\AVG_gfp.mat')
timeidx = find(GMFP_stat.mask ==1);  % GMFP: Global mean field power - GFP

%% plot
h = figure('Units','centimeters','Position',[10 10 5.9 3.5], 'IntegerHandle', 'off');
tiledlayout(1, 1, 'TileSpacing','tight','Padding','tight');
nexttile
rectangle('Position',[GMFP_stat.time(timeidx(1)) 0.1 (GMFP_stat.time(timeidx(623))-GMFP_stat.time(timeidx(1))) 49.8],'FaceColor',[0.9 0.9 0.9],...
    'EdgeColor','none')

rectangle('Position',[GMFP_stat.time(timeidx(624)) 0.1 (GMFP_stat.time(timeidx(end))-GMFP_stat.time(timeidx(624))) 49.8],'FaceColor',[0.9 0.9 0.9],...
    'EdgeColor','none')
hold on
plot(AVG_gmfp.time, 1e15*AVG_gmfp.avg,'k','LineWidth',1.2);

plot([0 0], [0 120],'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',0.5)
plot([-2.0 2.0], [1e15*base_gmfp 1e15*base_gmfp],'Color',[0.5 0.5 0.5],'LineStyle','--')

xlim([-1.5 2.0])
ylim([0 50])

set(gca,'XTick',[-1.5 0 2.0],'FontSize', 10)
set(gca,'YTick',[0 50],'FontSize', 10)
box on
ylabel('GFP (fT)','FontSize', 10, 'VerticalAlignment','bottom');
xlabel('Time (s)','FontSize', 10, 'VerticalAlignment','cap');
set(gca,'FontName', 'Calibri','FontWeight', 'bold');  % Times New Roman, Calibri



