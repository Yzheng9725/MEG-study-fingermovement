%% This figure is shown in Supplementary Figure 5.
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
opt.ana = 'Visual'; % 'Visual', 'Source_Space'
opt.condition = 'Visual';
opt.subana = 'SourceLevel';
opt.groupmatdir = fullfile(opt.filepath, 'mat', opt.ana, opt.subana);
opt.groupresultdir = fullfile(opt.filepath, 'results', opt.ana, opt.subana);

if ~exist(opt.groupresultdir, 'dir')
    mkdir(opt.groupresultdir);
end

cm = [0 0.447 0.741; 0.85 0.325 0.098; 0.929 0.694 0.125; 0.466 0.674 0.188];

FN = {'Allband','Delta','Theta','Alpha','Beta','Gamma','High-gamma','LowFreqBand'};
FB = {[0.5 90],[0.5 4],[4 8],[8 13],[13 30],[30 60],[60 90],[0.5 8]};

%% Load sourcemodel and brainatlas
load(fullfile(opt.mripath, 'sourcemodel.mat'));
load(fullfile(opt.mripath, 'sourcemodel_inflated.mat'));
filename = fullfile(opt.filepath, 'mat', 'other', 'S00', '8k', 'template.L.inflated.8k_fs_LR.surf.gii');
brainsurf_8k = ft_read_headshape({filename strrep(filename, '.L.', '.R.')});
brainsurf_8k_lft = ft_read_headshape(filename);
brainsurf_8k_rgt = ft_read_headshape(strrep(filename, '.L.', '.R.'));
brainatlas_8k = importdata(fullfile(opt.filepath, 'mat', 'other', 'S00', 'forMEG', 'brainatlas.mat'));
cbasebrain = [229 229 229]/255;
croibrain  = [102 102 102]/255;

%% Load MVPA results
load(fullfile(opt.groupmatdir, strcat(opt.condition, '_MVPA_ROI_pca_70var_chantime_group_16sub.mat')));

ACC = ACC_ROI_Group;
ACC_all = cat(3, ACC{:});
meanACC_ROI_Group = ACC_all;
meanACC_ROI_Group = permute(meanACC_ROI_Group, [3 1 2]);
meanACC_ROI_Group = squeeze(mean(meanACC_ROI_Group,1));

load(fullfile(opt.mripath, 'brainatlas.mat'));
brainatlas = brainatlasL;
brainatlas.parcellation = [brainatlasL.parcellation; brainatlasR.parcellation+size(brainatlasL.parcellationlabel,1)];
brainatlas.parcellationlabel = [brainatlasL.parcellationlabel; brainatlasR.parcellationlabel];

stat2plot.time = stat_mvpa.time;
stat2plot.accuracy = meanACC_ROI_Group;
stat2plot.label = brainatlas.parcellationlabel;

label2plot = [10 10+34 7 7+34 28 28+34 21 21+34 23 23+34];
Label2legend = cell(1, numel(label2plot));
for i = 1:numel(label2plot)
    Label2legend{i} = strcat(stat2plot.label{label2plot(i)}(1), '-', stat2plot.label{label2plot(i)}(3:end));
end

%% Plot ACC timecourse with boundedline
Actseg = [-0.5 2];
cmap = [187, 79, 122; 29, 115, 58]/255;
fig = figure('Units','centimeters','Position', [10 10 6 12], 'IntegerHandle', 'off');
tiledlayout(5, 1, 'TileSpacing', 'tight', 'Padding', 'tight');  % ,'TileIndexing','columnmajor'
for i = numel(label2plot)-1:-2:1
    nexttile; hold on
    % plot(stat2plot.time, stat2plot.accuracy(label2plot(i),:),'LineWidth',1.0);
    plot(Actseg, [25 25], '--', 'Color', 'k', 'LineWidth', 1); % chance level
    plot([0 0], [0 100], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1); % time zero line
    acc_avg = mean(ACC_all(label2plot([i i+1]), :, :), 3)*100;
    acc_std = std(ACC_all(label2plot([i i+1]), :, :), 0, 3)*100;
    acc_sem = acc_std / sqrt(size(ACC_all, 3));
    hl = gobjects(2,1); hp = gobjects(2,1);
    for j = [2 1] % right first
        [hl(j), hp(j)] = boundedline(stat2plot.time, acc_avg(j,:), acc_sem(j,:), 'alpha', 'transparency', 0.2, 'orientation', 'vert', 'LineWidth', 1.2, 'cmap', cmap(j, :));
        % outlinebounds(hl, hp);
    end
    xlim(Actseg);
    ylim([20 40]); set(gca, 'YTick', [20 25 30 35 40], 'YTickLabel', {'20','25','30','35','40'});
    if i == 1
        xlabel('Time (s)');
        ylabel('Accuracy (%)');
        legend(hl, {'L','R'}, 'Box', 'off', 'Location', 'northeast');
    else
        set(gca, 'XTickLabel', []);
    end
    title(Label2legend{i}(3:end), 'FontSize', 10, 'FontName', 'Arial');
    set(gca, 'Color', 'none');
end
set(fig, 'Name', 'MVPA-ACC-timecourse-Visual-ROIs-boundline');
print(fig, fullfile(opt.groupresultdir, strcat(fig.Name, '.png')), '-dpng', '-r600', '-image');
% print(fig, fullfile(opt.groupresultdir, strcat(fig.Name, '.eps')), '-depsc2', '-r600', '-vector');

%% Plot ROI in brain surface
cmap = [195 82 127; 30 119 60]/255;
for i = 1:2:numel(label2plot)
    roilft = find(brainatlasL.parcellation == label2plot(i));
    roirgt = find(brainatlasR.parcellation == label2plot(i));
    croi = repmat([nan nan nan], [size(brainsurf_8k.pos, 1) 1]);
    croilft = repmat([nan nan nan], [size(brainsurf_8k_lft.pos, 1) 1]);
    croirgt = repmat([nan nan nan], [size(brainsurf_8k_rgt.pos, 1) 1]);
    croi(roilft, :) = repmat(cmap(1,:), [numel(roilft) 1]); croi(roirgt+numel(brainatlasL.parcellation), :) = repmat(cmap(2,:), [numel(roirgt) 1]);
    croilft(roilft, :) = repmat(cmap(1,:), [numel(roilft) 1]);
    croirgt(roirgt, :) = repmat(cmap(2,:), [numel(roirgt) 1]);

    fig = figure('Units', 'centimeters', 'Position', [10 10 4 4], 'IntegerHandle', 'off');
    tiledlayout(1, 1, 'TileSpacing','tight','Padding','compact');  % ,'TileIndexing','columnmajor'
    nexttile; hold on;
    ft_plot_mesh(brainsurf_8k_lft, 'vertexcolor', repmat(cbasebrain, [size(brainsurf_8k_lft.pos, 1) 1]));
    ft_plot_mesh(brainsurf_8k_lft, 'vertexcolor', croilft);
    view([-90 0]); material dull; lighting gouraud; camlight headlight; set(gca, 'PlotBoxAspectRatio', [1 1 1]);
    set(fig, 'Name', sprintf('BrainSurface_ROI_%s_Left', Label2legend{i}(3:end)));
    print(fig, fullfile(opt.groupresultdir, strcat(fig.Name, '.png')), '-dpng', '-r600', '-image');
    % print(fig, fullfile(opt.groupresultdir, strcat(fig.Name, '.eps')), '-depsc2', '-r600', '-image');

    fig = figure('Units', 'centimeters', 'Position', [10 10 4 4], 'IntegerHandle', 'off');
    tiledlayout(1, 1, 'TileSpacing','tight','Padding','compact');  % ,'TileIndexing','columnmajor'
    nexttile; hold on;
    ft_plot_mesh(brainsurf_8k_rgt, 'vertexcolor', repmat(cbasebrain, [size(brainsurf_8k_rgt.pos, 1) 1]));
    ft_plot_mesh(brainsurf_8k_rgt, 'vertexcolor', croirgt);
    view([90 0]); material dull; lighting gouraud; camlight headlight; set(gca, 'PlotBoxAspectRatio', [1 1 1]);
    set(fig, 'Name', sprintf('BrainSurface_ROI_%s_Right', Label2legend{i}(3:end)));
    print(fig, fullfile(opt.groupresultdir, strcat(fig.Name, '.png')), '-dpng', '-r600', '-image');
    % print(fig, fullfile(opt.groupresultdir, strcat(fig.Name, '.eps')), '-depsc2', '-r600', '-image');
end

%% Plot ACC topomap at specific timepoints. Using sourcemodel_inflated
T = stat_mvpa.time';
toi = 0:0.05:2.0-0.05;
cmap = ft_colormap('Purples');
for tt = 1:length(toi)
    tidx2plot = knnsearch(T,toi(tt)):knnsearch(T,toi(tt)+0.05);
    ACC2plot = nan(size(sourcemodel.thickness));
    acc2toi = mean(meanACC_ROI_Group(:, tidx2plot),2);
    for i = 1:length(brainatlas.parcellationlabel)
        ACC2plot(brainatlas.parcellation == i) = acc2toi(i);
    end
    if mod(tt, 10) == 1
        fig1 = figure('Units','centimeters','Position',[10 10 40 4], 'IntegerHandle', 'off');
        tile1 = tiledlayout(fig1, 1, 10, 'TileSpacing','tight','Padding','tight'); colormap(cmap);
        fig2 = figure('Units','centimeters','Position',[10 10 40 4], 'IntegerHandle', 'off');
        tile2 = tiledlayout(fig2, 1, 10, 'TileSpacing','tight','Padding','tight'); colormap(cmap);
        pause(0.1);
    end
    if mod(tt, 10) == 0
        tt_plot = 10;
    else
        tt_plot = mod(tt, 10);
    end
    nexttile(tile1, tt_plot); hold on;
    clim([0.25 0.35]);
    ft_plot_mesh(brainsurf_8k,'vertexcolor',ACC2plot);
    view(0, 90); material dull; camlight headlight; lighting gouraud;
    set(gca, 'PlotBoxAspectRatio', [1 1 1], 'Color', 'none');

    nexttile(tile2, tt_plot); hold on;
    clim([0.25 0.35]);
    ft_plot_mesh(brainsurf_8k, 'vertexcolor', ACC2plot);
    view(0, 0); material dull; camlight headlight; lighting gouraud;
    set(gca, 'PlotBoxAspectRatio', [1 1 1], 'Color', 'none');
    % set(fig, 'Color','none');          % figure 背景透明
    if mod(tt, 10) == 0
        set(fig1, 'Name', sprintf('MVPA_ACC_ROI_Group_TopView_%.2f_to_%.2f_s', toi(tt-9), toi(tt)));
        set(fig2, 'Name', sprintf('MVPA_ACC_ROI_Group_BackView_%.2f_to_%.2f_s', toi(tt-9), toi(tt)));
        print(fig1, fullfile(opt.groupresultdir, strcat(fig1.Name, '.png')), '-dpng', '-r600', '-image');
        print(fig2, fullfile(opt.groupresultdir, strcat(fig2.Name, '.png')), '-dpng', '-r600', '-image');
        close(fig1); close(fig2);
        clear fig1 fig2 tile1 tile2;
    end
end

fig = figure('Units','centimeters','Position',[10 10 4 2], 'IntegerHandle', 'off');
tiledlayout(1, 1, 'TileSpacing','compact','Padding','compact');
nexttile; hold on;
axis off; colormap(cmap);
clim([0.25 0.35]);
colorbar('Location', 'southoutside', 'ticks', [0.25 0.3 0.35], 'TickLabels', {'0.25', '0.3', '0.35'}, 'Box', 'on');
print(fig, fullfile(opt.groupresultdir, 'Colorbar_MVPA_ACC_ROI_Group.png'), '-dpng', '-r600', '-image');
% print(fig, fullfile(opt.groupresultdir, 'Colorbar_MVPA_ACC_ROI_Group.eps'), '-depsc2', '-r600', '-vector');
