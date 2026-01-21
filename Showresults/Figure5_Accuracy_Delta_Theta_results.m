%% This figure is shown in Figure 5.
%%
clc; clear; close all;
ft_defaults;

%%
opt = [];
opt.filepath = fileparts(pwd);
opt.sub = 'S00';
opt.paradigm = 'Extension_4tasks';
opt.motion = 'Extension';
opt.mripath = fullfile(opt.filepath, 'mat', 'other', opt.sub, 'forMEG');
opt.task = {'Thumb','Index','Middle','Little'};
opt.ana = 'Source_Space';
opt.subana = 'Classification';
opt.groupmatdir = fullfile(opt.filepath, 'mat', opt.ana, opt.subana);
opt.groupresultdir = fullfile(opt.filepath, 'results', opt.ana, opt.subana);

if ~exist(opt.groupresultdir, 'dir')
    mkdir(opt.groupresultdir);
end

cm = [0 0.447 0.741; 0.85 0.325 0.098; 0.929 0.694 0.125; 0.466 0.674 0.188];

FN = {'Allband','Delta','Theta','Alpha','Beta','Gamma','High-gamma','LowFreqBand'};
FB = {[0.5 90],[0.5 4],[4 8],[8 13],[13 30],[30 60],[60 90],[0.5 8]};

%% Color map
cmap_grad = genColormap();
cmap_allvalue = [0 0 0];

%% Delta
%% Load data
opt.freqname = 'Delta';
load(fullfile(opt.groupmatdir, strcat(opt.freqname, '_meanACC_allwin_group_16sub.mat')));
load(fullfile(opt.groupmatdir, strcat(opt.freqname, '_meanACC_perwin_group_16sub.mat')));

%% plot ACC
Acc_show = cat(2, meanACC_maxAI_perwindow.mean, meanACC_maxAI_allwin.mean'); % win01-19 + allwin

% plot ACC
acc2show = mean(Acc_show*100)';
accstd = std(Acc_show*100)'; 

colors = cat(1, cmap_grad, cmap_allvalue);
fig = figure('Units','centimeters','Position',[10 10 10 5], 'IntegerHandle', 'off');
tiledlayout(1, 1, 'TileSpacing','tight','Padding','tight');
nexttile; hold on;
plot([size(acc2show, 1)-0.5 size(acc2show, 1)-0.5], get(gca, 'YLim'), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.6);
plot([0 size(acc2show, 1)+0.5], [25 25], '--', 'Color', 'k', 'LineWidth', 0.6);
for i = 1:size(Acc_show, 2)
    facealpha = 0.6;
    edgealpha = 1;
    scatter(repmat(i, size(Acc_show, 1),1), Acc_show(:,i)*100, 15, 'filled', 'o', 'MarkerFaceColor', colors(i,:), 'MarkerFaceAlpha', facealpha, 'MarkerEdgeColor', colors(i,:), 'MarkerEdgeAlpha', edgealpha, 'XJitter', 'rand', 'XJitterWidth', 0.3); % draw data scatter
    bar(i, acc2show(i), 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1, 'BarWidth', 0.7); % draw bar
    errorbar(i, acc2show(i), accstd(i), 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1, 'CapSize', 6); % draw error bar
    % , 'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', colors(i, :), 'MarkerEdgeColor', 'k'
end
ylim([20 90]); set(gca, 'YTick', 0:10:100); ylabel('Accuracy (%)');
xlim([0.5 20.5]);
set(gca, 'XTick', [], 'Color', 'none');
set(gcf, 'Name', 'showresults_Delta_Acc_allwin_perwin_compar_scatter');
print(fig, fullfile(opt.groupresultdir, strcat(fig.Name, '.png')), '-dpng', '-r600', '-vector');
% print(fig, fullfile(opt.groupresultdir, strcat(fig.Name, '.eps')), '-depsc2', '-r600', '-vector');

%% Delta plot confusionMatrix
predictlabels = [];
truelabels = [];
for ss = 1:numel(meanACC_maxAI_allwin.mean)
    for k = 1:5
        predictlabels = cat(1, predictlabels, meanACC_maxAI_allwin.predictlabel{ss,k});
        truelabels = cat(1, truelabels, [ones(20,1); 2*ones(20,1); 3*ones(20,1); 4*ones(20,1)]);
    end
end

confusionMatrix = confusionmat(truelabels, predictlabels);
confusionMatrix = confusionMatrix/1600*100;
confusionMatrix_delta = confusionMatrix;
xname = {'Thumb','Index','Middle','Little'}; yname = {'Thumb','Index','Middle','Little'};

fig = figure('Units','centimeters','Position',[10 10 9 7], 'IntegerHandle', 'off');
tiledlayout(1, 1, 'TileSpacing', 'tight', 'Padding', 'tight');
nexttile;
f1 = heatmap(xname, yname, confusionMatrix);
set(f1, 'Title', 'Confusion matrix', 'XLabel', 'Predict class', 'YLabel', 'True class', 'FontName', 'Arial', 'FontSize', 10, 'GridVisible', 'off', 'CellLabelFormat', '%.2f'); % , 'CellLabelFormat', '%.4f'
colormap(ft_colormap('Purples'));
clim([5 75]);
set(gcf, 'Name', 'showresults_ConfusionMatrix_Delta_allwin');
print(fig, fullfile(opt.groupresultdir, strcat(fig.Name, '.png')), '-dpng', '-r600', '-vector');
% print(fig, fullfile(opt.groupresultdir, strcat(fig.Name, '.eps')), '-depsc2', '-r600', '-vector');

%% Theta
%% Load data
opt.freqname = 'Theta';
load(fullfile(opt.groupmatdir, strcat(opt.freqname, '_meanACC_allwin_group_16sub.mat')));
load(fullfile(opt.groupmatdir, strcat(opt.freqname, '_meanACC_perwin_group_16sub.mat')));

%% plot ACC
Acc_show = nan(16, 20);
Acc_show(:, [3:9 10:15 20]) = cat(2, meanACC_maxAI_perwindow.mean, meanACC_maxAI_allwin.mean'); % win_xx + allwin

predictlabels = [];
truelabels = [];
for ss = 1:numel(meanACC_maxAI_allwin.mean)
    for k = 1:5
        predictlabels = cat(1, predictlabels, meanACC_maxAI_allwin.predictlabel{ss,k});
        truelabels = cat(1, truelabels, [ones(20,1); 2*ones(20,1); 3*ones(20,1); 4*ones(20,1)]);
    end
end

% plot ACC
subname = arrayfun(@(x) strcat('win', num2str(x, '%02d')), 1:size(meanACC_maxAI_perwindow.mean,2), 'UniformOutput', false);
subname = [subname(:)', {'all-win'}];
acc2show = mean(Acc_show*100)';
accstd = std(Acc_show*100)';

colors = cat(1, cmap_grad, cmap_allvalue);
fig = figure('Units','centimeters','Position',[10 10 10 5], 'IntegerHandle', 'off');
tiledlayout(1, 1, 'TileSpacing','tight','Padding','tight');
nexttile; hold on;
plot([size(acc2show, 1)-0.5 size(acc2show, 1)-0.5], get(gca, 'YLim'), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.6);
plot([0 size(acc2show, 1)+0.5], [25 25], '--', 'Color', 'k', 'LineWidth', 0.6);
for i = 1:size(Acc_show, 2)
    facealpha = 0.6;
    edgealpha = 1;
    scatter(repmat(i, size(Acc_show, 1),1), Acc_show(:,i)*100, 15, 'filled', 'o', 'MarkerFaceColor', colors(i,:), 'MarkerFaceAlpha', facealpha, 'MarkerEdgeColor', colors(i,:), 'MarkerEdgeAlpha', edgealpha, 'XJitter', 'rand', 'XJitterWidth', 0.3); % draw data scatter
    bar(i, acc2show(i), 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1, 'BarWidth', 0.7); % draw bar
    errorbar(i, acc2show(i), accstd(i), 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1, 'CapSize', 6); % draw error bar
    % , 'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', colors(i, :), 'MarkerEdgeColor', 'k'
end
ylim([20 90]); set(gca, 'YTick', 0:10:100);
xlim([0.5 20.5]);
set(gca, 'XTick', [], 'Color', 'none');
set(gcf, 'Name', 'showresults_Theta_Acc_allwin_perwin_compar_scatter');
print(fig, fullfile(opt.groupresultdir, strcat(fig.Name, '.png')), '-dpng', '-r600', '-vector');
% print(fig, fullfile(opt.groupresultdir, strcat(fig.Name, '.eps')), '-depsc2', '-r600', '-vector');

%% Theta plot confusionMatrix
confusionMatrix = confusionmat(truelabels, predictlabels);
confusionMatrix = confusionMatrix/1600*100;
confusionMatrix_theta = confusionMatrix;
xname = {'Thumb','Index','Middle','Little'}; yname = {'Thumb','Index','Middle','Little'};

fig = figure('Units','centimeters','Position',[10 10 9 7], 'IntegerHandle', 'off');
tiledlayout(1, 1, 'TileSpacing', 'tight', 'Padding', 'tight');
nexttile;
f1 = heatmap(xname, yname, confusionMatrix);
set(f1, 'Title', 'Confusion matrix', 'XLabel', 'Predict class', 'YLabel', 'True class', 'FontName', 'Arial', 'FontSize', 10, 'GridVisible', 'off', 'CellLabelFormat', '%.2f'); % , 'CellLabelFormat', '%.4f'
colormap(wxyz_colormap(16)); % 16/17/flip(80)
clim([5 75]);
set(gcf, 'Name', 'showresults_ConfusionMatrix_Theta_allwin');
print(fig, fullfile(opt.groupresultdir, strcat(fig.Name, '.png')), '-dpng', '-r600', '-vector');
% print(fig, fullfile(opt.groupresultdir, strcat(fig.Name, '.eps')), '-depsc2', '-r600', '-vector');

%% SUBFUNCTIONS
function [cmap_grad] = genColormap()
    cmap_grad = [0.9425 0.7933 0.4179;
                0.9512 0.8481 0.5862;
                0.9384 0.8647 0.6765;
                0.9240 0.8787 0.7625;
                0.8715 0.8533 0.8073;
                0.8512 0.8410 0.8159;
                0.8242 0.8193 0.8087;
                0.7974 0.7980 0.8022;
                0.7706 0.7767 0.7957;
                0.7411 0.7517 0.7830;
                0.7082 0.7220 0.7622;
                0.6754 0.6923 0.7416;
                0.6426 0.6627 0.7210;
                0.6088 0.6312 0.6963;
                0.5741 0.5983 0.6685;
                0.5395 0.5654 0.6406;
                0.5048 0.5324 0.6128;
                0.4718 0.5005 0.5838;
                0.4394 0.4689 0.5543];
end