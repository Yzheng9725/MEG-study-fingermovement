%% This figure is shown in Figure 2.
%%
clc; clear; close all;
ft_defaults;
ft_version;
[ftv, ftpath] = ft_version;

%% add toolbox path
addpath(fullfile(pwd, 'toolbox', 'Violinplot'));

%%
opt = [];
opt.filepath = fileparts(pwd);
opt.sub = 'S00';
opt.paradigm = 'Extension_4tasks';
opt.motion = 'Extension';
opt.mripath = fullfile(opt.filepath, 'mat', 'other', opt.sub, 'forMEG');
opt.task = {'Thumb','Index','Middle','Little'};
opt.ana = 'EMG';
opt.groupmatdir = fullfile(opt.filepath, 'mat', opt.ana);
opt.groupresultdir = fullfile(opt.filepath, 'results', opt.ana);

if ~exist(opt.groupresultdir, 'dir')
    mkdir(opt.groupresultdir);
end

cm = [59 135 199; 242 120 115; 255 211 115; 54 151 88]/255;
Subidx = 1:16; % 16 subjects

%% Load EMG parameter
load(fullfile(opt.groupmatdir, 'showresults_EMG_parameter.mat'));

%% Plot EMG parameter. Using violinplot toolbox
% Reaction time
fig = figure('Units','centimeters','Position',[10 10 5 4], 'IntegerHandle', 'off');
tiledlayout(1, 1, 'TileSpacing','compact','Padding','tight'); hold on;
violinobj = violinplot(Reactime_show.mean, opt.task, 'QuartileStyle', 'boxplot', 'HalfViolin', 'full', 'DataStyle', 'scatter', 'ViolinColor', cm, 'ShowData', true, 'ViolinAlpha', 0.5, 'Width', 0.3, 'ShowMean', false, 'ShowMedian', true, 'MarkerSize', 10);
for i = 1:length(violinobj)
    violinobj(i).BoxPlot.LineWidth = 3;
    violinobj(i).EdgeColor = 'none';
    violinobj(i).MedianPlot.LineWidth = 1;
    violinobj(i).MedianPlot.SizeData = 30;
    violinobj(i).ScatterPlot.MarkerEdgeColor = cm(i,:);
    violinobj(i).ScatterPlot.MarkerFaceAlpha = 0.3;
    violinobj(i).ScatterPlot.SizeData = 10;
    violinobj(i).ViolinPlot.FaceColor = 'none';
    violinobj(i).ViolinPlot.EdgeColor = 'none';
    violinobj(i).MeanPlot.LineWidth = 2;
    violinobj(i).WhiskerPlot.LineWidth = 1;
end
xlim([0.5 4.5]); ylim([450 950]);
ylabel('Time (ms)', 'FontWeight','bold','VerticalAlignment','middle');
box off;
set(gca, 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'none');
set(gca, 'YTickLabelRotation', 90, 'YAxisLocation', 'right');
set(gcf, 'Name', 'Show-EMG-parameter-16sub-Reactiontime-violinplot');
print(fig, fullfile(opt.groupresultdir, strcat(fig.Name, '.png')), '-dpng', '-r600', '-vector');
% print(fig, fullfile(opt.groupresultdir, strcat(fig.Name, '.eps')), '-depsc2', '-r600', '-vector');

% Movement time
fig = figure('Units','centimeters','Position',[10 10 5 4], 'IntegerHandle', 'off');
tiledlayout(1, 1, 'TileSpacing','compact','Padding','tight'); hold on;
violinobj = violinplot(Motortime_show.mean, opt.task, 'QuartileStyle', 'boxplot', 'HalfViolin', 'full', 'DataStyle', 'scatter', 'ViolinColor', cm, 'ShowData', true, 'ViolinAlpha', 0.5, 'Width', 0.3, 'ShowMean', false, 'ShowMedian', true, 'MarkerSize', 10);
for i = 1:length(violinobj)
    violinobj(i).BoxPlot.LineWidth = 3;
    violinobj(i).EdgeColor = 'none';
    violinobj(i).MedianPlot.LineWidth = 1;
    violinobj(i).MedianPlot.SizeData = 30;
    violinobj(i).ScatterPlot.MarkerEdgeColor = cm(i,:);
    violinobj(i).ScatterPlot.MarkerFaceAlpha = 0.3;
    violinobj(i).ScatterPlot.SizeData = 10;
    violinobj(i).ViolinPlot.FaceColor = 'none';
    violinobj(i).ViolinPlot.EdgeColor = 'none';
    violinobj(i).MeanPlot.LineWidth = 2;
    violinobj(i).WhiskerPlot.LineWidth = 1;
end
xlim([0.5 4.5]); ylim([50 650]);
ylabel('Time (ms)', 'FontWeight','bold','VerticalAlignment','middle');
box off;
set(gca, 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'none');
set(gca, 'YTickLabelRotation', 90, 'YAxisLocation', 'right');
set(gcf, 'Name', 'Show-EMG-parameter-16sub-Movementtime-violinplot');
print(fig, fullfile(opt.groupresultdir, strcat(fig.Name, '.png')), '-dpng', '-r600', '-vector');
% print(fig, fullfile(opt.groupresultdir, strcat(fig.Name, '.eps')), '-depsc2', '-r600', '-vector');

%% Load classification results
load(fullfile(opt.groupmatdir, 'showresults_confusionmartix.mat'));
load(fullfile(opt.groupmatdir, 'showresults_ACC_EMG_-500_2000ms.mat'));

%% Plot Classification. Using barplot+scatter
subname = arrayfun(@(x) strcat('S', num2str(x, '%02d')), 1:numel(Subidx), 'UniformOutput', false); subname = cat(2, subname, {'AVG'});
acc = Acc_show.mean*100; accstd = std(Acc_show.mean*100, 0, 1);
fig = figure('Units','centimeters','Position',[10 10 2 5], 'IntegerHandle', 'off');
tiledlayout(1, 1, 'TileSpacing','compact','Padding','tight'); hold on;
nexttile; hold on;
cdata = [84 39 143]/255;
hbar = bar(mean(acc), 'FaceColor', 'none', 'EdgeColor', cdata, 'LineWidth', 1, 'BarWidth', 1);
scatter(hbar.XEndPoints.'*ones(1, numel(acc)), acc, 15, 'filled', 'o', 'MarkerEdgeColor', cdata, 'MarkerEdgeAlpha', 1, 'MarkerFaceAlpha', 0.3,...
    'CData', cdata, 'XJitter', 'rand', 'XJitterWidth', 0.6); % draw data scatter
errorbar(hbar.XEndPoints, mean(acc), accstd, 'vertical', 'LineStyle', 'none', 'LineWidth', 1.2, 'Color', 'k', 'CapSize', 6);
set(gca, 'xtick', 1, 'xticklabel', {'EMG'}, 'Color', 'none');
ylim([50 100]);
ylabel('Accuracy (%)');
set(gcf, 'Name', 'Show-EMG-Classification-Accuracy');
print(fig, fullfile(opt.groupresultdir, strcat(fig.Name, '.png')), '-dpng', '-r600', '-vector');
% print(fig, fullfile(opt.groupresultdir, strcat(fig.Name, '.eps')), '-depsc2', '-r600', '-vector');

%% Plot confusion matrix
[confusionMatrix, order] = confusionmat(truelabels, predictlabels);
confusionMatrix = confusionMatrix/(sum(confusionMatrix(:))/4)*100;
xname = {'Thumb','Index','Middle','Little'};
yname = {'Thumb','Index','Middle','Little'};
fig = figure('Units','centimeters','Position',[10 10 6 5], 'IntegerHandle', 'off');
tiledlayout(1, 1, 'TileSpacing', 'tight', 'Padding', 'tight');
nexttile;
f1 = heatmap(xname, yname, confusionMatrix);
set(f1, 'Title', 'Confusion matrix', 'XLabel', 'Predict class', 'YLabel', 'True class', 'FontName', 'Arial', 'FontSize', 10, 'GridVisible', 'off', 'CellLabelFormat', '%.2f'); % , 'CellLabelFormat', '%.4f'
colormap(ft_colormap('Purples'));
clim([0 90]);
set(gca, 'FontName', 'Arial', 'FontSize', 10);
set(gcf, 'Name', 'Show-EMG-ConfusionMatrix');
print(fig, fullfile(opt.groupresultdir, strcat(fig.Name, '.png')), '-dpng', '-r600', '-vector');
% print(fig, fullfile(opt.groupresultdir, strcat(fig.Name, '.eps')), '-depsc2', '-r600', '-vector');

%% remove toolbox path
rmpath(fullfile(pwd, 'toolbox', 'Violinplot'));