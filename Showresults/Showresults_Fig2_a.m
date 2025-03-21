%%
% Finger extension movement
% Fig.2(a)
% creator: WxyZ
% Date: 20250314

%%
ft_defaults;
ft_version;

load('.\mat\SensorLevel\Timelock_AVG_Thumb.mat')
load('.\mat\SensorLevel\Timelock_AVG_Index.mat')
load('.\mat\SensorLevel\Timelock_AVG_Middle.mat')
load('.\mat\SensorLevel\Timelock_AVG_Little.mat')
load('.\mat\SensorLevel\TFWAVE_AVG_alltaskavg.mat')

cm = [59 135 199; 242 120 115; 255 211 115; 54 151 88]/255;

%% 
h = figure('Units','centimeters','Position',[10 10 8.8 5.7], 'IntegerHandle', 'off');
drawnow;
tiledlayout(2,1,'TileSpacing','tight','Padding','tight');

%% AVG - MLF67
cfg                 = [];
cfg.baseline        = [-1.5 -1.0];
cfg.xlim            = [-1.5 2.0];
cfg.linewidth       = 1.5;
cfg.linecolor       = cm;
cfg.figure          = gcf;

nexttile([1 1]);
cfg.channel         =  'MLF67';  
ft_singleplotER(cfg, AVG_Thumb, AVG_Index, AVG_Middle, AVG_Little);

hold on
plot([-2 2], [0 0], 'LineStyle','--', 'Color',[0 0 0 0.5],'LineWidth',0.8)
plot([0 0], [-30 30], 'LineStyle','--','Color',[0 0 0 0.5], 'LineWidth',0.8);

ylim([-1.5e-13 2e-13])
title([]);
box on;
set(gca,'XTick',[ ],'XColor','none')
set(gca,'YTick',[-1.5e-13 0 2e-13],'FontSize', 10,'YTickLabel',{'-150','0','200'}) 
ylabel('Amplitude (fT)','FontSize', 10, 'VerticalAlignment','middle');

legend({'Thumb','Index','Middle','Little'},'FontSize',10,'Box','off','Location','eastoutside') 
set(gca,'FontName', 'Calibri', 'FontWeight', 'bold');  

%% TFR - MLF67
cfg                 = [];
cfg.baseline        = [-1.5 -1.0];
cfg.baselinetype    = 'relative';
cfg.xlim            = [-1.5 2.0];
cfg.ylim            = [8 90];
cfg.zlim            = [0.5 1.5];
cfg.layout          = 'CTF275_helmet.mat';
cfg.colorbar        = 'no';
cfg.masknans        = 'no';
cfg.showlabels      = 'no';
cfg.limittext       = ['Time:[' num2str(cfg.xlim(1)) ' ' num2str(cfg.xlim(2)) ']s' newline ...
    'Freq:[' num2str(cfg.ylim(1)) ' ' num2str(cfg.ylim(2)) ']Hz' newline ...
    'Baseline:[' num2str(cfg.baseline(1)) ' ' num2str(cfg.baseline(2)) ']s'];
cfg.showoutline     = 'yes'; % show head/helmet outline
cfg.showscale       = 'no'; % remove bottom right corner graph
cfg.figure          = gcf;

nexttile([1 1]);
cfg.channel         = 'MLF67';   
ft_singleplotTFR(cfg, TFWAVE_AVG_all)

hold on
plot([0 0], [8 90], 'LineStyle','--','Color',[0 0 0 0.5], 'LineWidth',0.8);

set(gca,'XTick',[-1.5 0 2.0],'FontSize', 10)
set(gca,'YTick',[8 90],'FontSize', 10)
ylabel('Frequency (Hz)','FontSize', 10,'VerticalAlignment','bottom');   %, 'FontWeight', 'bold'
colormap(wxyz_colormap(1)) % RdBu4
c= colorbar;

set(c,'Box','off','FontSize',10,'Ticks',[0.5 1 1.5],'TickLabels',{'-50', '0','50'})  % ,'Position',c_pos
c.Label.String = 'Power (%)';
c.Label.VerticalAlignment = 'baseline';
box on;
title([]);
set(gca, 'FontName', 'Calibri','FontWeight', 'bold');  % , 'FontWeight', 'bold','FontSize', 12, Arial


