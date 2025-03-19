%%
% Finger extension movement
% Fig.3(a),(b)
% creator: WxyZ - Yu Zheng
% Date: 20250314

%%
ft_defaults;

cm = [0 0.447 0.741; 0.85 0.325 0.098; 0.929 0.694 0.125; 0.466 0.674 0.188];

%%  Low-freqband
opt.freqname = 'LowFreqBand';
load(['.\mat\SourceLevel\' opt.freqname '\Source_Group.mat'])
load(['.\mat\SourceLevel\' opt.freqname '\Source_time_avg_pertime.mat'])

%% plot pseudo-statistics - Lowfreqband
% plot Time
h = figure('Units','centimeters','Position',[10 10 7 7], 'IntegerHandle', 'off');
tiledlayout(2, 1, 'TileSpacing','compact','Padding','tight');
drawnow;

% boxplot
nexttile([1 1])
time2plot = reshape(Sourcetime_Group(:,:,:),16,[])+0.01;
boxplot(time2plot,'Widths',0.65 ,'Symbol','+k','PlotStyle','compact','MedianStyle','line','BoxStyle','outline','Orientation','horizontal'...
    ,'Colors',cm) 
lineObj = findobj(gca,'Type','Line');
for i = 1:length(lineObj)
    lineObj(i).LineWidth = 0.75;
end

hold on
yline(2.5, ':', 'LineWidth', 0.5, 'Color',[0.3 0.3 0.3],'HandleVisibility', 'off')
yline(6.5, ':', 'LineWidth', 0.5, 'Color',[0.3 0.3 0.3],'HandleVisibility', 'off')
yline(10.5, ':', 'LineWidth', 0.5, 'Color',[0.3 0.3 0.3],'HandleVisibility', 'off')
yline(14.5, ':', 'LineWidth', 0.5, 'Color',[0.3 0.3 0.3],'HandleVisibility', 'off')

plot([0 0], [0 20],'Color',[0.5 0.5 0.5],'LineStyle','--')
ylim([0 17])
xlim([-0.5 1.0])

set(gca,'XTick',[ ], 'XColor','none')
set(gca,'YTick',2.5:4:16,'YTickLabel',{'t1','t2','t3','t4'})
box off


set(gca,'FontName', 'Calibri','FontSize', 10, 'FontWeight', 'bold');  % Times New Roman

%% plot norm pZ - joy plot
mean2plot = [];
for t = 1:4
     for ss = 1:16
        Pz_norm(ss,t,:) = squeeze(pZ_Group(ss,t,:))/max(pZ_Group(ss,t,:));
     end
    mean2plot(t,:) = squeeze(mean(Pz_norm(:,t,:),1));
end

data2plot = mean2plot';

lgLable= {'Thumb','Index','Middle','Little'};
colors=[0 0.447 0.741; 0.85 0.325 0.098; 0.929 0.694 0.125; 0.466 0.674 0.188];

nexttile([1 1])

p = 0.2;
x = (-0.5:0.02:0.98) + 0.01;
for i = size(data2plot, 2):-1:1
    f = data2plot(:,i)';
    fShifted = f + i * p;
    pHandle = plot(x, fShifted,'color',colors(i,:),'LineStyle','-','LineWidth', 1.2,'HandleVisibility', 'off');
    hold on;
    yline(i * p, '-', 'LineWidth', 0.5, 'Color',[0.7 0.7 0.7],'HandleVisibility', 'off')
    Xfill = [x, fliplr(x)];
    Yfill = [fShifted, ones(1, length(x)) * i * p];
    fill(Xfill, Yfill, pHandle.Color, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
end
% yTick = (1:size(data2plot, 2)) * p;
box off

plot([0 0], [0 10],'Color',[0.5 0.5 0.5],'LineStyle','--')

xlim([-0.5 1.0])
ylim([0.15 1.55])
% set(gca, 'YTick', yTick, 'YTickLabel', lgLable, 'TickDir', 'out','YColor','k');
set(gca, 'YTick', []);
set(gca,'XTick',[-0.5 0 1.0])
xlabel('Time (s)','VerticalAlignment','Middle');
ylabel(['Normalized ' '\itpseudo-z'],'VerticalAlignment','bottom');

for t = 1:4
lin(t) = line(nan,nan,'Color',cm(t,:),'Linestyle','-','LineWidth',1.5);
end
legend([lin(1) lin(2) lin(3) lin(4)],{'Thumb','Index','Middle','Little'},'FontSize',10,'Box','off','Location','northeast')

sgtitle('Low-Frequency Band','FontName', 'Calibri','FontSize', 12, 'FontWeight', 'bold');
set(gca,'FontName', 'Calibri','FontSize', 10, 'FontWeight', 'bold');  % Times New Roman



%%  plot peak time - alpha & beta & gamma
h = figure('Units','centimeters','Position',[10 10 7 7], 'IntegerHandle', 'off');
tiledlayout(8, 1, 'TileSpacing','tight','Padding','tight');
drawnow;

% alpha
opt.freqname = 'Alpha';
load(['.\mat\SourceLevel\' opt.freqname '\Source_Group.mat'])
load(['.\mat\SourceLevel\' opt.freqname '\Source_time_avg_pertime.mat'])

MEANpeaktime = squeeze(mean(mean(Sourcetime_Group(:,:,:),1)));
for pk = 1:2
    stdpeak(pk) = std(Sourcetime_Group(:,:,pk),0,'all');
end

nexttile([3 1])
time2plot = reshape(Sourcetime_Group(:,:,:),16,[])+0.25;

boxplot(time2plot,'Widths',0.6 ,'Symbol','+k','PlotStyle','compact','MedianStyle','line','BoxStyle','outline','Orientation','horizontal'...
    ,'Colors',cm) 
lineObj = findobj(gca,'Type','Line');
for i = 1:length(lineObj)
    lineObj(i).LineWidth = 0.75;
end

hold on
yline(2.5, ':', 'LineWidth', 0.5, 'Color',[0.3 0.3 0.3],'HandleVisibility', 'off')
yline(6.5, ':', 'LineWidth', 0.5, 'Color',[0.3 0.3 0.3],'HandleVisibility', 'off')

plot([0 0], [0 20],'Color',[0.5 0.5 0.5],'LineStyle','--')
ylim([0 9])
xlim([-0.5 2.0])

set(gca,'XTick',[ ] ,'XColor','none')
set(gca,'YTick',2.5:4:10,'YTickLabel',{'ERD','ERS'})
box off
set(gca,'FontName', 'Calibri','FontSize', 10, 'FontWeight', 'bold');  % Times New Roman



% beta
opt.freqname = 'Beta';
load(['.\mat\SourceLevel\' opt.freqname '\Source_Group.mat'])
load(['.\mat\SourceLevel\' opt.freqname '\Source_time_avg_pertime.mat'])

MEANpeaktime = squeeze(mean(mean(Sourcetime_Group(:,:,:),1)));
for pk = 1:2
    stdpeak(pk) = std(Sourcetime_Group(:,:,pk),0,'all');
end

nexttile([3 1])
time2plot = reshape(Sourcetime_Group(:,:,:),16,[])+0.25;
boxplot(time2plot,'Widths',0.6 ,'Symbol','+k','PlotStyle','compact','MedianStyle','line','BoxStyle','outline','Orientation','horizontal'...
    ,'Colors',cm)
lineObj = findobj(gca,'Type','Line');
for i = 1:length(lineObj)
    lineObj(i).LineWidth = 0.75;
end

hold on
yline(2.5, ':', 'LineWidth', 0.5, 'Color',[0.3 0.3 0.3],'HandleVisibility', 'off')
yline(6.5, ':', 'LineWidth', 0.5, 'Color',[0.3 0.3 0.3],'HandleVisibility', 'off')
plot([0 0], [0 20],'Color',[0.5 0.5 0.5],'LineStyle','--')
ylim([0 9])
xlim([-0.5 2.0])

set(gca,'XTick',[ ] ,'XColor','none')
set(gca,'YTick',2.5:4:10,'YTickLabel',{'ERD','ERS'})
box off
set(gca,'FontName', 'Calibri','FontSize', 10, 'FontWeight', 'bold');  % Times New Roman



% gamma
opt.freqname = 'Gamma';
load(['.\mat\SourceLevel\' opt.freqname '\Source_Group.mat'])
load(['.\mat\SourceLevel\' opt.freqname '\Source_time_avg_pertime.mat'])

% boxplot
nexttile([2 1])
time2plot = reshape(Sourcetime_Group(:,:,:),16,[])+0.25;
boxplot(time2plot,'Widths',0.6 ,'Symbol','+k','PlotStyle','compact','MedianStyle','line','BoxStyle','outline','Orientation','horizontal'...
    ,'Colors',cm) 
lineObj = findobj(gca,'Type','Line');
for i = 1:length(lineObj)
    lineObj(i).LineWidth = 1;
end

hold on
yline(2.5, ':', 'LineWidth', 0.5, 'Color',[0.3 0.3 0.3],'HandleVisibility', 'off')
plot([0 0], [0 20],'Color',[0.5 0.5 0.5],'LineStyle','--')
ylim([0 5])
xlim([-0.5 2.0])

set(gca,'XTick',[-0.5 0 1.0 2.0])
xlabel('Time (s)','VerticalAlignment','cap');
set(gca,'YTick',2.5,'YTickLabel','ERS')
box off
set(gca,'FontName', 'Calibri','FontSize', 10, 'FontWeight', 'bold');  % Times New Roman


