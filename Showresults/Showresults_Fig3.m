%%
% Finger extension movement
% Fig.3(a),(b)
% creator: WxyZ
% Date: 20250314

%%
ft_defaults;

cm = [59 135 199; 242 120 115; 255 211 115; 54 151 88]/255;

%%  Low-freqband
opt.freqname = 'LowFreqBand';
load(['.\mat\SourceLevel\' opt.freqname '\Source_Group.mat'])
load(['.\mat\SourceLevel\' opt.freqname '\Source_time_avg_pertime.mat'])

% plot pseudo-statistics - Lowfreqband
% plot Time
h = figure('Units','centimeters','Position',[10 10 7 9], 'IntegerHandle', 'off');
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

% plot norm pZ - joy plot
mean2plot = [];
for t = 1:4
     for ss = 1:16
        Pz_norm(ss,t,:) = squeeze(pZ_Group(ss,t,:))/max(pZ_Group(ss,t,:));
     end
    mean2plot(t,:) = squeeze(mean(Pz_norm(:,t,:),1));
end

data2plot = mean2plot';

nexttile([1 1])
p = 0.2;
x = (-0.5:0.02:0.98) + 0.01;
for i = size(data2plot, 2):-1:1
    f = data2plot(:,i)';
    fShifted = f + i * p;
    pHandle = plot(x, fShifted,'color',cm(i,:),'LineStyle','-','LineWidth', 1.2,'HandleVisibility', 'off');
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

sgtitle('Low-Frequency Band','FontName', 'Calibri','FontSize', 12, 'FontWeight', 'bold');
set(gca,'FontName', 'Calibri','FontSize', 10, 'FontWeight', 'bold');  % Times New Roman



%%  plot peak time - alpha & beta & gamma
%% alpha
h = figure('Units','centimeters','Position',[10 10 7 9], 'IntegerHandle', 'off');
tiledlayout(3, 1, 'TileSpacing','tight','Padding','tight');
drawnow;

opt.freqname = 'Alpha';
load(['.\mat\SourceLevel\' opt.freqname '\Source_Group.mat'])
load(['.\mat\SourceLevel\' opt.freqname '\Source_time_avg_pertime.mat'])

MEANpeaktime = squeeze(mean(mean(Sourcetime_Group(:,:,:),1)));
for pk = 1:2
    stdpeak(pk) = std(Sourcetime_Group(:,:,pk),0,'all');
end

nexttile
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
xlim([-0.25 1.75])

set(gca,'XTick',[ ] ,'XColor','none')
set(gca,'YTick',2.5:4:10,'YTickLabel',{'ERD','ERS'})
box off
set(gca,'FontName', 'Calibri','FontSize', 10, 'FontWeight', 'bold');  % Times New Roman


% ERS
mean2plot = [];
Pt_norm = [];
for t = 1:4
     for ss = 1:16
        Pt_norm(ss,t,:) = squeeze(pTmax_Group(ss,t,:))/max(pTmax_Group(ss,t,:));
     end
    mean2plot(t,:) = squeeze(mean(Pt_norm(:,t,:),1));
end

data2plot = mean2plot';
lgLable= {'Thumb','Index','Middle','Little'};

nexttile
p = 0.25;
x = (-0.5:0.1:1.5) + 0.25;
for i = size(data2plot, 2):-1:1
    f = data2plot(:,i)';
    fShifted = f + i * p;
    pHandle = plot(x, fShifted,'color',cm(i,:),'LineStyle','-','LineWidth', 1.2,'HandleVisibility', 'off');
    hold on;
    yline(i * p, '-', 'LineWidth', 0.5, 'Color',[0.7 0.7 0.7],'HandleVisibility', 'off')
    Xfill = [x, fliplr(x)];
    Yfill = [fShifted, ones(1, length(x)) * i * p];
    fill(Xfill, Yfill, pHandle.Color, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
end
% yTick = (1:size(data2plot, 2)) * p;
box off

plot([0 0], [0 10],'Color',[0.5 0.5 0.5],'LineStyle','--')

xlim([-0.25 1.75])
ylim([0.15 1.6])

set(gca,'XTick',[ ] ,'XColor','none')
set(gca, 'YTick', []);
ylabel(['Normalized ' '\itpseudo-t'],'VerticalAlignment','bottom');
title('ERS','FontName', 'Calibri','FontSize', 12, 'FontWeight', 'bold');
set(gca,'FontName', 'Calibri','FontSize', 10, 'FontWeight', 'bold');


% ERD
mean2plot = [];
Pt_norm = [];
for t = 1:4
     for ss = 1:16
        Pt_norm(ss,t,:) = squeeze(pTmin_Group(ss,t,:))/max(pTmax_Group(ss,t,:));
     end
    mean2plot(t,:) = squeeze(mean(Pt_norm(:,t,:),1));
end

data2plot = mean2plot';

lgLable= {'Thumb','Index','Middle','Little'};

nexttile([1 1])
p = 0.5;
for i = 1:1:size(data2plot, 2)  %:-1:1
    f =  data2plot(:,i)';
    fShifted = f + i * p;
    pHandle = plot(x, fShifted,'color',cm(i,:),'LineStyle','-','LineWidth', 1.2,'HandleVisibility', 'off');
    hold on;
    yline(i * p, '-', 'LineWidth', 0.5, 'Color',[0.7 0.7 0.7],'HandleVisibility', 'off')
    Xfill = [x, fliplr(x)];
    Yfill = [fShifted, ones(1, length(x)) * i * p];
    fill(Xfill, Yfill, pHandle.Color, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
end
% yTick = (1:size(data2plot, 2)) * p;
box off

plot([0 0], [-10 10],'Color',[0.5 0.5 0.5],'LineStyle','--')

xlim([-0.25 1.75])
ylim([-0.5 2.25])

set(gca,'XTick',[-0.25 0 1.75]);
set(gca, 'YTick', []);
title('ERD','FontName', 'Calibri','FontSize', 12, 'FontWeight', 'bold');
xlabel('Time (s)','VerticalAlignment','Middle');
set(gca,'FontName', 'Calibri','FontSize', 10, 'FontWeight', 'bold');  


%% beta
h = figure('Units','centimeters','Position',[10 10 7 9], 'IntegerHandle', 'off');
tiledlayout(3, 1, 'TileSpacing','tight','Padding','tight');
drawnow;

opt.freqname = 'Beta';
load(['.\mat\SourceLevel\' opt.freqname '\Source_Group.mat'])
load(['.\mat\SourceLevel\' opt.freqname '\Source_time_avg_pertime.mat'])

MEANpeaktime = squeeze(mean(mean(Sourcetime_Group(:,:,:),1)));
for pk = 1:2
    stdpeak(pk) = std(Sourcetime_Group(:,:,pk),0,'all');
end

nexttile
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


% ERS
mean2plot = [];
Pt_norm = [];
for t = 1:4
     for ss = 1:16
        Pt_norm(ss,t,:) = squeeze(pTmax_Group(ss,t,:))/max(pTmax_Group(ss,t,:));
     end
    mean2plot(t,:) = squeeze(mean(Pt_norm(:,t,:),1));
end

data2plot = mean2plot';
lgLable= {'Thumb','Index','Middle','Little'};

nexttile
p = 0.25;
x = (-0.5:0.1:1.5) + 0.25;
for i = size(data2plot, 2):-1:1
    f = data2plot(:,i)';
    fShifted = f + i * p;
    pHandle = plot(x, fShifted,'color',cm(i,:),'LineStyle','-','LineWidth', 1.2,'HandleVisibility', 'off');
    hold on;
    yline(i * p, '-', 'LineWidth', 0.5, 'Color',[0.7 0.7 0.7],'HandleVisibility', 'off')
    Xfill = [x, fliplr(x)];
    Yfill = [fShifted, ones(1, length(x)) * i * p];
    fill(Xfill, Yfill, pHandle.Color, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
end
% yTick = (1:size(data2plot, 2)) * p;
box off

plot([0 0], [0 10],'Color',[0.5 0.5 0.5],'LineStyle','--')

xlim([-0.25 1.75])
ylim([0.15 1.8])

set(gca,'XTick',[ ] ,'XColor','none')
set(gca, 'YTick', []);
ylabel(['Normalized ' '\itpseudo-t'],'VerticalAlignment','bottom');
title('ERS','FontName', 'Calibri','FontSize', 12, 'FontWeight', 'bold');
set(gca,'FontName', 'Calibri','FontSize', 10, 'FontWeight', 'bold');


% ERD
mean2plot = [];
Pt_norm = [];
for t = 1:4
     for ss = 1:16
        Pt_norm(ss,t,:) = squeeze(pTmin_Group(ss,t,:))/max(pTmax_Group(ss,t,:));
     end
    mean2plot(t,:) = squeeze(mean(Pt_norm(:,t,:),1));
end

data2plot = mean2plot';

lgLable= {'Thumb','Index','Middle','Little'};

nexttile
p = 0.5;
for i = 1:1:size(data2plot, 2)  %:-1:1
    f =  data2plot(:,i)';
    fShifted = f + i * p;
    pHandle = plot(x, fShifted,'color',cm(i,:),'LineStyle','-','LineWidth', 1.2,'HandleVisibility', 'off');
    hold on;
    yline(i * p, '-', 'LineWidth', 0.5, 'Color',[0.7 0.7 0.7],'HandleVisibility', 'off')
    Xfill = [x, fliplr(x)];
    Yfill = [fShifted, ones(1, length(x)) * i * p];
    fill(Xfill, Yfill, pHandle.Color, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
end
% yTick = (1:size(data2plot, 2)) * p;
box off

plot([0 0], [-10 10],'Color',[0.5 0.5 0.5],'LineStyle','--')

xlim([-0.25 1.75])
ylim([-0.5 2.25])

set(gca,'XTick',[-0.25 0 1.75]);
set(gca, 'YTick', []);
title('ERD','FontName', 'Calibri','FontSize', 12, 'FontWeight', 'bold');
xlabel('Time (s)','VerticalAlignment','Middle');
set(gca,'FontName', 'Calibri','FontSize', 10, 'FontWeight', 'bold');  


%% gamma
h = figure('Units','centimeters','Position',[10 10 7 8], 'IntegerHandle', 'off');
tiledlayout(5, 1, 'TileSpacing','tight','Padding','tight');
drawnow;

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
xlim([-0.25 1.25])
set(gca,'XTick',[ ] ,'XColor','none')
set(gca,'YTick',2.5,'YTickLabel','ERS')
box off
set(gca,'FontName', 'Calibri','FontSize', 10, 'FontWeight', 'bold');  % Times New Roman


% ERS
mean2plot = [];
Pt_norm = [];
for t = 1:4
     for ss = 1:16
        Pt_norm(ss,t,:) = squeeze(pTmax_Group(ss,t,:))/max(pTmax_Group(ss,t,:));
     end
    mean2plot(t,:) = squeeze(mean(Pt_norm(:,t,:),1));
end

data2plot = mean2plot';
lgLable= {'Thumb','Index','Middle','Little'};

nexttile([3 1])
p = 0.3;
x = (-0.5:0.1:1) + 0.25;
for i = size(data2plot, 2):-1:1
    f = data2plot(:,i)';
    fShifted = f + i * p;
    pHandle = plot(x, fShifted,'color',cm(i,:),'LineStyle','-','LineWidth', 1.2,'HandleVisibility', 'off');
    hold on;
    yline(i * p, '-', 'LineWidth', 0.5, 'Color',[0.7 0.7 0.7],'HandleVisibility', 'off')
    Xfill = [x, fliplr(x)];
    Yfill = [fShifted, ones(1, length(x)) * i * p];
    fill(Xfill, Yfill, pHandle.Color, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
end
% yTick = (1:size(data2plot, 2)) * p;
box off

plot([0 0], [0 10],'Color',[0.5 0.5 0.5],'LineStyle','--')

xlim([-0.25 1.25])
ylim([0.2 2.25])
set(gca,'XTick',[-0.5 0 1.0 2.0])
xlabel('Time (s)','VerticalAlignment','cap');
set(gca, 'YTick', []);
ylabel(['Normalized ' '\itpseudo-t'],'VerticalAlignment','bottom');
title('ERS','FontName', 'Calibri','FontSize', 12, 'FontWeight', 'bold');

xlabel('Time (s)','VerticalAlignment','Middle');
set(gca,'FontName', 'Calibri','FontSize', 10, 'FontWeight', 'bold');  

for t = 1:4
lin(t) = line(nan,nan,'Color',cm(t,:),'Linestyle','-','LineWidth',1.5);
end
legend([lin(1) lin(2) lin(3) lin(4)],{'Thumb','Index','Middle','Little'},'FontSize',10,'Box','off','Location','southoutside')


