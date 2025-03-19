%%
% Fig.2(b)
% creator: WxyZ
% Date: 20250314

%%
ft_defaults;

baseline = [-1.5 -1.0];
load('.\mat\SensorLevel\ITPC2plot.mat')
load('.\mat\SensorLevel\ITPC_stat_ICA.mat')

%%
h = figure('Units','centimeters','Position',[10 10 8.1 3.6], 'IntegerHandle', 'off');
drawnow;
tiledlayout(1,1,'TileSpacing','tight','Padding','tight');
cfg                 = [];
cfg.parameter       = 'itpc';
cfg.xlim            = [-1.5 2.0];
cfg.ylim            = [0.5 30];
cfg.zlim            = [0 0.5];
cfg.layout          = 'CTF275_helmet.mat';
cfg.colorbar        = 'no';
cfg.masknans        = 'no';
cfg.showlabels      = 'no';
cfg.showoutline     = 'yes'; % show head/helmet outline
cfg.showscale       = 'no'; % remove bottom right corner graph
cfg.figure          = gcf;
cfg.style              = 'contour';

nexttile
cfg.channel         = 'MLF67';   
ft_singleplotTFR(cfg, ITPC2plot)

hold on
contour(ITPC2plot.time(26:201),0.5:0.5:30,squeeze(ITPC_stat.mask(55,1:60,1:176)),[1 1],'LineWidth',1,'EdgeColor','k')
% 55 - MLF67; 1:60 - 0.5-30 Hz

plot([0 0], [0 120],'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',0.5)
ylim([0.5 30])
set(gca,'YTick',[8 30],'FontSize', 10)
box off

colormap(wxyz_colormap(6))
c = colorbar;
set(c,'Ticks',[0 0.5],'YTickLabel',{'0','0.5'},'FontSize', 10)
c.Label.String = 'ITPC';
c.Label.VerticalAlignment = 'bottom';
xlabel('Time/s'); 
ylabel('Frequency (Hz)','FontSize', 10, 'VerticalAlignment','bottom');
title([]);
box on;
set(gca,'FontName', 'Calibri', 'FontWeight', 'bold');  % Times New Roman



















