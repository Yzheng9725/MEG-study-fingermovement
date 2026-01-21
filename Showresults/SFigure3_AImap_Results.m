%% This figure is shown in Supplementary Figure 3.
%%
clc; clear; close all;
ft_defaults;

%% Group
opt = [];
opt.filepath = fileparts(pwd);
opt.sub = 'S00';
opt.paradigm = 'Extension_4tasks';
opt.motion = 'Extension';
opt.mripath = fullfile(opt.filepath, 'mat', 'other', opt.sub, 'forMEG');
opt.task = {'Thumb','Index','Middle','Little'};
opt.ana = 'Source_Space';
opt.subana = 'Activation';
opt.groupmatdir = fullfile(opt.filepath, 'mat', opt.ana, opt.subana);
opt.groupresultdir = fullfile(opt.filepath, 'results', opt.ana, opt.subana);

if ~exist(opt.groupresultdir, 'dir')
    mkdir(opt.groupresultdir);
end

FN = {'Allband','Delta','Theta','Alpha','Beta','Gamma','High-gamma','LowFreqBand'};
FB = {[0.5 90],[0.5 4],[4 8],[8 13],[13 30],[30 60],[60 90],[0.5 8]};
timewin = 0.05;
timestep = 0.05;
Timeseg = -1:timestep:1.95;
Actseg = [-1 2.0];

%% plot AI - actseg (-1~2 s)
load(fullfile(opt.groupmatdir, 'showresults_AImap_Actseg.mat'));
load(fullfile(opt.mripath, 'sourcemodel_inflated.mat'));
brain_lft = ft_read_headshape(fullfile(opt.mripath, 'my.L.inflated_8k_fs_LR.surf.gii'));
brain_rgt = ft_read_headshape(fullfile(opt.mripath, 'my.R.inflated_8k_fs_LR.surf.gii'));
cbasebrain = [229 229 229]/255;
croibrain  = [102 102 102]/255;

climits = cell(size(FN));
climits{1} = [-0.15 0.15];
climits{2} = [-0.4 0.4];
climits{3} = [-0.2 0.2];
climits{4} = [-0.1 0.1];
climits{5} = [-0.15 0.15];
climits{6} = [-0.1 0.1];
climits{7} = [-0.1 0.1];
climits{8} = [-0.25 0.25];
pictype = {'png'};

fig = figure('Units', 'centimeters', 'Position', [10 10 4 24], 'IntegerHandle', 'off');
tiledlayout(6, 1, 'TileSpacing','tight','Padding','tight');
cmap = flip(ft_colormap('RdBu'));
for f = 2:7
    AI2plot = AI_map_show;
    v2plot = mean(cat(1, AI2plot{:, f}), 1)';
    v2plot_lft = v2plot(1:size(brain_lft.pos, 1)); % left hemi
    cdat = v2plot_lft;
    tile = nexttile; hold on;
    ft_plot_mesh(brain_lft, 'vertexcolor', repmat(cbasebrain, [size(brain_lft.pos, 1) 1]));
    ft_plot_mesh(brain_rgt, 'vertexcolor', repmat(cbasebrain, [size(brain_lft.pos, 1) 1]));
    if strcmpi(pictype{1}, 'png')
        ft_plot_mesh(brain_lft, 'vertexcolor', cdat);
    end
    view(-120, 40);
    material dull; lighting gouraud; camlight headlight;
    clim(climits{f});
    fprintf('Plotting %s. Climits = [%.4f %.4f]\n', FN{f}, min(cdat), max(cdat));
    cb = colorbar;
    colormap(cmap);
    title(FN{f}, 'Position', [0 -25 110]);
    set(gca, 'PlotBoxAspectRatio', [1 1 1], 'FontName', 'Arial', 'FontSize', 10);
    drawnow;
end
set(gcf, 'Name', 'showresults-AImap-AverageAcrossTime');
print(fig, fullfile(opt.groupresultdir, strcat(fig.Name, '.png')), '-dpng', '-r600', '-image');

%% show - slide window - AImap
load(fullfile(opt.mripath, 'sourcemodel_inflated.mat'));
v2plot = zeros(size(sourcemodel_inflated.thickness));

%% Plot necessary time windows
load(fullfile(opt.groupmatdir, 'showresults_AImap_50ms_window.mat'));
twin = {[], [0:0.05:0.45], [0:0.05:0.45], [0:0.05:0.2 1:0.05:1.2], [0:0.05:0.2 1:0.05:1.2], [0:0.05:0.2 1:0.05:1.2], [0:0.05:0.45], []};
twinidx = cell(size(twin));
for t = 1:numel(twin)
    if ~isempty(twin{t})
        twinidx{t} = knnsearch(Timeseg', twin{t}');
    end
end

% Get climtis
climits = cell(size(twin));
for f = 2:7
    AI2plot = AI_map_show;
    v2plot = mean(cat(3, AI2plot{:, f}), 3)'; % avg across tasks
    v2plot_lft = v2plot(1:size(brain_lft.pos, 1), :); % left hemi
    clim_min = nan(1, numel(twin{f}));
    clim_max = nan(1, numel(twin{f}));
    for c = 1:numel(twin{f})
        cdat = v2plot_lft(:, twinidx{f}(c));
        clim_min(c) = min(cdat);
        clim_max(c) = max(cdat);
    end
    climits{f} = [min(clim_min) max(clim_max)];
    if abs(climits{f}(1)) >= abs(climits{f}(2))
        climits{f}(2) = abs(climits{f}(1));
    else
        climits{f}(1) = -1 * abs(climits{f}(2));
    end
end

climits{2} = [-0.5 0.5];
climits{3} = [-0.3 0.3];
climits{4} = [-0.2 0.2];
climits{5} = [-0.3 0.3];
climits{6} = [-0.2 0.2];
climits{7} = [-0.3 0.3];

brain_lft = ft_read_headshape(fullfile(opt.mripath, 'my.L.inflated_8k_fs_LR.surf.gii'));
brain_rgt = ft_read_headshape(fullfile(opt.mripath, 'my.R.inflated_8k_fs_LR.surf.gii'));
cbasebrain = [229 229 229]/255;
fig = figure('Units', 'centimeters', 'Position', [10 10 40 24], 'IntegerHandle', 'off');
tiledlayout(6, 10, 'TileSpacing', 'compact', 'Padding', 'tight');  % ,'TileIndexing','columnmajor'
for f = 2:7 % Delta/Theta/Alpha/Beta/Gamma/High-gamma
    AI2plot = AI_map_show;
    v2plot = mean(cat(3, AI2plot{:, f}), 3)'; % avg across tasks
    v2plot_lft = v2plot(1:size(brain_lft.pos, 1), :); % left hemi
    fprintf('Plotting %s. Climits = [%.4f %.4f]', FN{f}, climits{f}(1), climits{f}(2));
    for c = 1:numel(twin{f})
        fprintf(', %.2fs', twin{f}(c));
        cdat = v2plot_lft(:, twinidx{f}(c));
        tile = nexttile; hold on;
        ft_plot_mesh(brain_lft, 'vertexcolor', repmat(cbasebrain, [size(brain_lft.pos, 1) 1]));
        ft_plot_mesh(brain_rgt, 'vertexcolor', repmat(cbasebrain, [size(brain_lft.pos, 1) 1]));
        ft_plot_mesh(brain_lft, 'vertexcolor', cdat);
        view(-120, 40);
        material dull; lighting gouraud; camlight headlight;
        clim(climits{f});
        if c == numel(twin{f})
            cb = colorbar;
            ylabel(cb, 'A');
        end
        colormap(tile, cmap);
        if c == 1
            title([FN{f} ' ' num2str(twin{f}(c)) ' s'],'Position',[0 -25 110]);
        else
            title([num2str(twin{f}(c)) ' s'],'Position',[0 -25 110]);
        end
        
        set(gca, 'PlotBoxAspectRatio', [1 1 1], 'FontName', 'Arial', 'FontSize', 10);
    end
    drawnow;
    fprintf('\n');
end
set(gcf, 'Name', 'showresults-AImap-50ms-selected-alltimewindow');
print(fig, fullfile(opt.groupresultdir, strcat(fig.Name, '.png')), '-dpng', '-r600', '-image');

