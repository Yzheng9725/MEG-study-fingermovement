%% This figure is shown in Supplementary Figure 4.
%%
clc; clear; close all;
ft_defaults;

%%
opt = [];
opt.filepath = fileparts(pwd);
opt.sub = 'S00';
opt.paradigm = 'Extension_4tasks';
opt.motion = 'Extension';
opt.mripath = fullfile('E:', 'MRI', opt.sub, 'forMEG');
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

load(fullfile(opt.groupmatdir, 'meanACC_pca_70var_TFR_ACTseg_-1_2s_group_16sub.mat'));
load(fullfile(opt.groupmatdir, 'meanACC_maxAI_TFR_ACTseg_-1_2s_group_16sub.mat'));

Acc_show(:,:,1) = meanACC_pca_group(:,2:7);
Acc_show(:,:,2) = meanACC_maxAI_group(:,2:7);

%% Plot using scatter points
subname = FN(2:7);
acc2show = Acc_show;
accmean = [squeeze(mean(Acc_show(:,:,1)*100))' squeeze(mean(Acc_show(:,:,2)*100))'];
accstd = [std(Acc_show(:,:,1)*100)' std(Acc_show(:,:,2)*100)'];
fig = figure('Units','centimeters','Position',[10 10 10 8], 'IntegerHandle', 'off');
tiledlayout(1, 1, 'TileSpacing', 'tight', 'Padding', 'tight');
nexttile; hold on;
plot([0 size(accmean, 1)+0.5], [25 25], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.6);
hbar = bar(accmean);
herrbar = gobjects(size(hbar));
marker = {'o', 's'};
hs = gobjects(1, numel(hbar));
for h = 1:numel(hbar)
    hbar(h).FaceColor = 'none';
    hbar(h).EdgeColor = 'k';
    hbar(h).BarWidth  = 0.9;
    % scatter individual points
    for k = 1:length(hbar(h).XEndPoints)
        hs(h) = scatter(repelem(hbar(h).XEndPoints(k), size(Acc_show, 1)), acc2show(:,k,h)*100, 15, k*ones(size(acc2show, 1), 1),...
                'Marker', marker{h}, 'MarkerEdgeColor', 'flat', 'MarkerEdgeAlpha', 1, 'MarkerFaceColor', 'flat', 'MarkerFaceAlpha', 0.3, 'XJitter', 'rand', 'XJitterWidth', 0.06);
    end
    herrbar(h) = errorbar(hbar(h).XEndPoints, hbar(h).YEndPoints, accstd(:,h), 'k', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 3);
end
colors = [65, 82, 84; 73, 159, 162; 162, 212, 210; 245, 215, 221; 245, 175, 191; 123, 105, 110]/255;
colormap(colors);
ylim([0 60]); set(gca, 'YTick', 0:10:60); ylabel('Accuracy (%)');
xlim([0.5 6.5]); set(gca, 'XTick', 1:6, 'XTickLabel', subname);
set(gca, 'Color', 'none');
legend(hs, {'PCA', 'A-epoch'}, 'Box', 'off', 'Location', 'northeast');
set(fig, 'Name', 'showresults_Acc_TFR_PCA_MaxAI_winappend_compar_barscatter_symbol');
print(fig, fullfile(opt.groupresultdir, strcat(fig.Name, '.png')), '-dpng', '-r600', '-image');
% print(fig, fullfile(opt.groupresultdir, strcat(fig.Name, '.eps')), '-depsc2', '-r600', '-vector');

%% Plot only Method comparison
accavg = squeeze(mean(Acc_show*100, 1));
accstd = squeeze(std(Acc_show*100, [], 1));
colors = [65, 82, 84; 73, 159, 162; 162, 212, 210; 245, 215, 221; 245, 175, 191; 123, 105, 110]/255;
fig = figure('Units','centimeters','Position',[10 10 6 8], 'IntegerHandle', 'off');
tiledlayout(1, 1, 'TileSpacing', 'tight', 'Padding', 'tight');
nexttile; hold on;
plot([0 size(accavg, 1)+0.5], [25 25], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.6);
herrbar = gobjects(1, size(accavg, 1));
margin = [-0.06 -0.04 0.02 0.02 0.04 0.06];
for i = 1:size(accavg, 1)
    herrbar(i) = plot(1:size(accavg,2), accavg(i,:), 'LineStyle', '-', 'LineWidth', 1, 'Color', colors(i,:), 'Marker', '^', 'MarkerSize', 6, 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'none');
    errorbar((1:size(accavg,2)), accavg(i,:), accstd(i,:), 'LineStyle', 'none', 'LineWidth', 1, 'Color', colors(i,:), 'CapSize', 3);
end
% set(herrbar(1), 'MarkerFaceAlpha', 0.3); % for legend display
ylim([0 60]); set(gca, 'YTick', 0:10:60); ylabel('Accuracy (%)');
xlim([0.5 2.5]); set(gca, 'XTick', 1:2, 'XTickLabel', {'PCA', 'A-epoch'});
set(gca, 'Color', 'none');
set(fig, 'Name', 'showresults_Acc_TFR_PCA_MaxAI_winappend_compar_method');
print(fig, fullfile(opt.groupresultdir, strcat(fig.Name, '.png')), '-dpng', '-r600', '-image');
% print(fig, fullfile(opt.groupresultdir, strcat(fig.Name, '.eps')), '-depsc2', '-r600', '-vector');