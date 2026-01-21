%% This figure is shown in Supplementary Figure 2.
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

%% Load flat surface and inflated surface
brainsurf = struct('surf_8k', struct(), 'surf_32k', struct());
s = {'32k'};
for h = {'L', 'R'}
    brainsurf.(strcat('surf_', s{1})).(h{1}).inflated = ft_read_headshape(fullfile(opt.filepath, 'mat', 'other', 'S00', s{1}, sprintf('template.%s.inflated.%s_fs_LR.surf.gii', h{1}, s{1})));
    brainsurf.(strcat('surf_', s{1})).(h{1}).flat = ft_read_headshape(fullfile(opt.filepath, 'mat', 'other', 'S00', s{1}, sprintf('template.%s.flat.%s_fs_LR.surf.gii', h{1}, s{1})));
    brainsurf.(strcat('surf_', s{1})).(h{1}).label = ft_read_atlas(fullfile(opt.filepath, 'mat', 'other', 'S00', s{1}, sprintf('template.%s.aparc.%s_fs_LR.label.gii', h{1}, s{1})));
    brainsurf.(strcat('surf_', s{1})).(h{1}).sulc = double(gifti(fullfile(opt.filepath, 'mat', 'other', 'S00', s{1}, sprintf('template.%s.sulc.%s_fs_LR.shape.gii', h{1}, s{1}))).cdata);
end
cbasebrain = [204 204 204]/255;
cprecentral = [212 185 218]/255;
cpostcentral = [223 101 176]/255;

%%
% inflated
s = {'32k'};
fig = figure('Units', 'centimeters', 'Position', [10 10 4 4], 'IntegerHandle', 'off');
tiledlayout(1, 1, 'TileSpacing', 'tight', 'Padding', 'tight');  % ,'TileIndexing','columnmajor'
nexttile; hold on;
ft_plot_mesh(brainsurf.(strcat('surf_', s{1})).L.inflated, 'vertexcolor', brainsurf.(strcat('surf_', s{1})).L.inflated.sulc);
ft_plot_mesh(brainsurf.(strcat('surf_', s{1})).R.inflated, 'vertexcolor', repmat(cbasebrain, size(brainsurf.(strcat('surf_', s{1})).R.inflated.pos, 1), 1), 'facealpha', 0.3);
view(-120, 40);
material dull; camlight headlight; lighting gouraud;
colormap(ft_colormap('Greys'));
clim([-15 15]);
set(fig, 'Name', sprintf('AI_pos_showresults_inflated_%s', s{1}));
print(fig, fullfile(opt.groupresultdir, strcat(fig.Name, '.png')), '-dpng', '-r600', '-image');
% print(fig, fullfile(opt.groupresultdir, strcat(fig.Name, '.eps')), '-depsc2', '-r600', '-vector');

% inflated highlight pre/post-central
fig = figure('Units', 'centimeters', 'Position', [10 10 4 4], 'IntegerHandle', 'off');
tiledlayout(1, 1, 'TileSpacing', 'tight', 'Padding', 'tight');  % ,'TileIndexing','columnmajor'
nexttile; hold on;
ft_plot_mesh(brainsurf.(strcat('surf_', s{1})).L.inflated, 'vertexcolor', brainsurf.(strcat('surf_', s{1})).L.inflated.sulc, 'facealpha', 0.3);
surf_label = brainsurf.(strcat('surf_', s{1})).L.label;
cvertex = nan(numel(surf_label.parcellation), 1);
cvertex(surf_label.parcellation == 21, :) = brainsurf.(strcat('surf_', s{1})).L.sulc(surf_label.parcellation == 21); % Postcentral
cvertex(surf_label.parcellation == 23, :) = brainsurf.(strcat('surf_', s{1})).L.sulc(surf_label.parcellation == 23); % Precentral
ft_plot_mesh(brainsurf.(strcat('surf_', s{1})).L.inflated, 'vertexcolor', cvertex, 'facealpha', 1);
ft_plot_mesh(brainsurf.(strcat('surf_', s{1})).R.inflated, 'vertexcolor', repmat(cbasebrain, size(brainsurf.(strcat('surf_', s{1})).R.inflated.pos, 1), 1), 'facealpha', 0.3);
view(-120, 40);
material dull; camlight headlight; lighting gouraud;
colormap(ft_colormap('Greys'));
clim([-15 15]);
set(fig, 'Name', sprintf('AI_pos_showresults_inflated_%s_roi', s{1}));
print(fig, fullfile(opt.groupresultdir, strcat(fig.Name, '.png')), '-dpng', '-r600', '-image');
% print(fig, fullfile(opt.groupresultdir, strcat(fig.Name, '.eps')), '-depsc2', '-r600', '-vector');

% flat surface
surf_tmp = brainsurf.(strcat('surf_', s{1})).L.flat;
fig = figure('Units', 'centimeters', 'Position', [10 10 14 10], 'IntegerHandle', 'off');
tiledlayout(1, 1, 'TileSpacing', 'tight', 'Padding', 'tight');
nexttile; hold on;
ft_plot_mesh(surf_tmp, 'vertexcolor', surf_tmp.sulc, 'facealpha', 0.3);
surf_label = brainsurf.(strcat('surf_', s{1})).L.label;
cvertex = nan(numel(surf_label.parcellation), 1);
cvertex(surf_label.parcellation == 21, :) = brainsurf.(strcat('surf_', s{1})).L.sulc(surf_label.parcellation == 21); % Postcentral
cvertex(surf_label.parcellation == 23, :) = brainsurf.(strcat('surf_', s{1})).L.sulc(surf_label.parcellation == 23); % Precentral
colormap(ft_colormap('Greys'));
clim([-15 15]);
ft_plot_mesh(surf_tmp, 'vertexcolor', cvertex, 'facealpha', 1);
% draw roi boundary
rectangle('Position', [-45 -40 120 190], 'EdgeColor', 'k', 'LineWidth', 1.5);
% axis normal;
axis off;
% xlim([-45 75]); ylim([-40 150]);
set(fig, 'Name', sprintf('AI_pos_showresults_flat_%s', s{1}));
print(fig, fullfile(opt.groupresultdir, strcat(fig.Name, '.png')), '-dpng', '-r600', '-image');
% print(fig, fullfile(opt.groupresultdir, strcat(fig.Name, '.eps')), '-depsc2', '-r600', '-vector');

% flat highlight pre/post-central
surf_tmp = brainsurf.(strcat('surf_', s{1})).L.flat;
fig = figure('Units', 'centimeters', 'Position', [10 10 8 10], 'IntegerHandle', 'off');
tiledlayout(1, 1, 'TileSpacing', 'tight', 'Padding', 'tight');
nexttile; hold on;
ft_plot_mesh(surf_tmp, 'vertexcolor', surf_tmp.sulc, 'facealpha', 0.3);
surf_label = brainsurf.(strcat('surf_', s{1})).L.label;
roiprecentral = find(surf_label.parcellation == 23);
roipostcentral = find(surf_label.parcellation == 21);
cvertex = nan(numel(surf_label.parcellation), 1);
cvertex(surf_label.parcellation == 21, :) = brainsurf.(strcat('surf_', s{1})).L.sulc(surf_label.parcellation == 21); % Postcentral
cvertex(surf_label.parcellation == 23, :) = brainsurf.(strcat('surf_', s{1})).L.sulc(surf_label.parcellation == 23); % Precentral
colormap(ft_colormap('Greys'));
clim([-15 15]);
ft_plot_mesh(surf_tmp, 'vertexcolor', cvertex, 'facealpha', 1);
axis normal;
axis off;
xlim([-45 75]); ylim([-40 150]);
set(fig, 'Name', sprintf('AI_pos_showresults_flat_%s_roi', s{1}));
print(fig, fullfile(opt.groupresultdir, strcat(fig.Name, '.png')), '-dpng', '-r600', '-image');
% print(fig, fullfile(opt.groupresultdir, strcat(fig.Name, '.eps')), '-depsc2', '-r600', '-vector');

%% Plot max/min AI pos (per time window)
% 50 ms time window
load(fullfile(opt.groupmatdir, 'showresults_maxAI_value_pos_50ms_window.mat'));
cm = [59 135 199; 242 120 115; 255 211 115; 0 138 69]/255; % thumb/index/middle/little

% map AI into size
psize = [10 150];
popc = [0 1];
maxAI_avg = cellfun(@(x) mean(x, 1), maxAI, 'UniformOutput', false); minAI_avg = cellfun(@(x) mean(x, 1), minAI, 'UniformOutput', false);
maxAI_avg_norm = nan(4, 8, 60); minAI_avg_norm = nan(4, 8, 60); maxAI_size = nan(4, 8, 60); minAI_size = nan(4, 8, 60);
maxAI_avg_opc = nan(4, 8, 60); minAI_avg_opc = nan(4, 8, 60); maxAI_opc = nan(4, 8, 60); minAI_opc = nan(4, 8, 60);
for t = 1:4
    for f = 2:7
        maxAI_avg_norm(t, f, :) = maxAI_avg{t, f};
        minAI_avg_norm(t, f, :) = minAI_avg{t, f};
        maxAI_avg_opc(t, f, :) = maxAI_avg{t, f};
        minAI_avg_opc(t, f, :) = minAI_avg{t, f};
        % normalize
        maxAI_size(t, f, :) = psize(1) + psize(2)*(maxAI_avg_norm(t, f, :) - min(squeeze(maxAI_avg_norm(t, f, :))))/(max(squeeze(maxAI_avg_norm(t, f, :))) - min(squeeze(maxAI_avg_norm(t, f, :))));
        minAI_size(t, f, :) = psize(2) - (psize(2)-psize(1))*(minAI_avg_norm(t, f, :) - min(squeeze(minAI_avg_norm(t, f, :))))/(max(squeeze(minAI_avg_norm(t, f, :))) - min(squeeze(minAI_avg_norm(t, f, :))));
        maxAI_opc(t, f, :) = popc(1) + popc(2)*(maxAI_avg_opc(t, f, :) - min(squeeze(maxAI_avg_opc(t, f, :))))/(max(squeeze(maxAI_avg_opc(t, f, :))) - min(squeeze(maxAI_avg_opc(t, f, :))));
        minAI_opc(t, f, :) = popc(2) - (popc(2)-popc(1))*(minAI_avg_opc(t, f, :) - min(squeeze(minAI_avg_opc(t, f, :))))/(max(squeeze(minAI_avg_opc(t, f, :))) - min(squeeze(minAI_avg_opc(t, f, :))));
    end
end

% map time into color
pcmap = GenColormap(60);
pcolor = nan(4, 60, 3);
pcolor(1, :, :) = pcmap.thumb/255;
pcolor(2, :, :) = pcmap.index/255;
pcolor(3, :, :) = pcmap.middle/255;
pcolor(4, :, :) = pcmap.little/255;

if ~exist(fullfile(opt.groupresultdir, 'AI_pos_50ms_window'), 'dir')
    mkdir(fullfile(opt.groupresultdir, 'AI_pos_50ms_window'));
end

for f = 2:7 % Delta/Theta/Alpha/Beta/Gamma/High-gamma
    surf_tmp = brainsurf.surf_32k.L.flat;
    surf_real = brainsurf.surf_8k.L.flat;
    for type = {'max', 'min'}
        for t = 1:4 % Thumb/Index/Middle/Little
            % draw flat surface
            fig = figure('Units','centimeters','Position',[10 10 8 10], 'IntegerHandle', 'off');
            tiledlayout(1, 1, 'TileSpacing', 'tight', 'Padding', 'compact'); nexttile; hold on;
            ft_plot_mesh(surf_tmp, 'vertexcolor', surf_tmp.sulc, 'facealpha', 0.3);
            surf_label = brainsurf.(strcat('surf_', s{1})).L.label;
            roiprecentral = find(surf_label.parcellation == 23); roipostcentral = find(surf_label.parcellation == 21);
            cvertex = nan(numel(surf_label.parcellation), 1);
            cvertex(roiprecentral, :) = brainsurf.(strcat('surf_', s{1})).L.sulc(roiprecentral); % Precentral
            cvertex(roipostcentral, :) = brainsurf.(strcat('surf_', s{1})).L.sulc(roipostcentral); % Postcentral
            ft_plot_mesh(surf_tmp, 'vertexcolor', cvertex, 'facealpha', 1);
            colormap(ft_colormap('Greys'));
            xlim([-45 75]); ylim([-40 150]); clim([-15 15]);
            axis normal; axis off;
            % draw pos max/min idx per time window
            for c = 1:1:60 % time windows
                if strcmpi(type{1}, 'max')
                    % draw pos maxidx
                    pos2plot = mean(surf_real.pos(squeeze(maxidx{t, f}(:, c)), :), 1);
                    scatter(pos2plot(1), pos2plot(2), maxAI_size(t, f, c), squeeze(pcolor(t, c, :))', 'o', 'filled', 'MarkerFaceAlpha', 1); %, 'MarkerEdgeColor', cm(t, :)
                    set(fig, 'Name', sprintf('%s_%s_maxAI_Pos_50ms_window_Lcortex_flat', opt.task{t}, FN{f}));
                elseif strcmpi(type{1}, 'min')
                    % draw pos minidx
                    pos2plot = mean(surf_real.pos(squeeze(minidx{t, f}(:, c)), :), 1);
                    scatter(pos2plot(1), pos2plot(2), minAI_size(t, f, c), squeeze(pcolor(t, c, :))', 'square', 'filled', 'MarkerFaceAlpha', 1); %, 'MarkerEdgeColor', cm(t, :)
                    set(fig, 'Name', sprintf('%s_%s_minAI_Pos_50ms_window_Lcortex_flat', opt.task{t}, FN{f}));
                end
            end
            drawnow;
            print(fig, fullfile(opt.groupresultdir, 'AI_pos_50ms_window', strcat(fig.Name, '.png')), '-dpng', '-r600', '-image');
            % print(fig, fullfile(opt.groupresultdir, 'AI_pos_50ms_window', strcat(fig.Name, '.eps')), '-depsc2', '-r600', '-vector');
        end
    end
    close all;
end

%% Plot colormap
fig = figure('Units','centimeters','Position',[10 10 16 2], 'IntegerHandle', 'off');
tiledlayout(1, 4, 'TileSpacing', 'tight', 'Padding', 'compact');
pcmap = GenColormap(60);
pcolor = nan(4, 60, 3);
pcolor(1, :, :) = pcmap.thumb/255;
pcolor(2, :, :) = pcmap.index/255;
pcolor(3, :, :) = pcmap.middle/255;
pcolor(4, :, :) = pcmap.little/255;

for t = 1:4
    nexttile; hold on;
    colormap(gca, squeeze(pcolor(t, :, :)));
    colorbar('Location', 'southoutside', 'Ticks', [], 'box', 'off');
    xlim([1 60]); ylim([0 1]);
    axis off;
    if t == 1
        title('Thumb');
    elseif t == 2
        title('Index');
    elseif t == 3
        title('Middle');
    elseif t == 4
        title('Little');
    end
end
set(fig, 'Name', 'AI_pos_50ms_window_Colormap');
wxyz_decorFigure(fig, 'FontName', 'Arial', 'FontSize', 8);
print(fig, fullfile(opt.groupresultdir, 'AI_pos_50ms_window', strcat(fig.Name, '.png')), '-dpng', '-r600', '-image');
print(fig, fullfile(opt.groupresultdir, 'AI_pos_50ms_window', strcat(fig.Name, '.eps')), '-depsc2', '-r600', '-vector');

%% Plot circle size
fig = figure('Units','centimeters','Position',[10 10 8 10], 'IntegerHandle', 'off');
tiledlayout(1, 1, 'TileSpacing', 'tight', 'Padding', 'compact');
psize = [10 120];
nexttile; hold on;
scatter(0, 0, psize(1), 'k', 'o', 'filled'); axis off; xlim([-45 75]); ylim([-40 150]); %title('min AI');
scatter(20, 0, psize(2), 'k', 'o', 'filled'); axis off; xlim([-45 75]); ylim([-40 150]); %title('max AI');
scatter(0, 20, psize(1), 'k', 'square', 'filled'); axis off; xlim([-45 75]); ylim([-40 150]); %title('min AI');
scatter(20, 20, psize(2), 'k', 'square', 'filled'); axis off; xlim([-45 75]); ylim([-40 150]); %title('max AI');
set(fig, 'Name', 'AI_pos_50ms_window_CircleSize');
wxyz_decorFigure(fig, 'FontName', 'Arial', 'FontSize', 8);
print(fig, fullfile(opt.groupresultdir, 'AI_pos_50ms_window', strcat(fig.Name, '.png')), '-dpng', '-r600', '-image');
print(fig, fullfile(opt.groupresultdir, 'AI_pos_50ms_window', strcat(fig.Name, '.eps')), '-depsc2', '-r600', '-vector');

%% SUBFUNCTION
function cmap = GenColormap(ncolor)
    % Thumb: Blue
    cmap.thumb = [
    169 203 232;  % #a9cbe8
    150 190 226;  % #96bee2
    130 177 219;  % #82b1db
    111 165 213;  % #6fa5d5
     90 152 207;  % #5a98cf
     77 143 202;  % #4d8fca
     63 133 196;  % #3f85c4
     47 124 191;  % #2f7cbf
     37 115 185;  % #2573b9
     27 107 179;  % #1b6bb3
     15  98 172;  % #0f62ac
      0  90 166   % #005aa6
    ];
    % Index: Red
    cmap.index = [
    246 214 213;  % #f6d6d5
    248 197 195;  % #f8c5c3
    248 179 177;  % #f8b3b1
    248 162 158;  % #f8a29e
    246 144 140;  % #f6908c
    244 132 127;  % #f4847f
    241 119 113;  % #f17771
    238 106 100;  % #ee6a64
    234  95  89;  % #ea5f59
    230  85  77;  % #e6554d
    226  73  66;  % #e24942
    221  61  54   % #dd3d36
    ];
    % Middle: Yellow
    cmap.middle = [
    255 246 230;  % #fff6e6
    255 239 208;  % #ffefd0
    255 233 187;  % #ffe9bb
    255 226 166;  % #ffe2a6
    255 220 144;  % #ffdc90
    255 214 128;  % #ffd680
    255 209 111;  % #ffd16f
    255 203  94;  % #ffcb5e
    255 197  78;  % #ffc54e
    255 190  61;  % #ffbe3d
    255 184  41;  % #ffb829
    255 177   8   % #ffb108
    ];
    % Little: Green
    cmap.little = [
    110 246 178;  % #6ef6b2
     93 225 156;  % #5de19c
     76 204 135;  % #4ccc87
     58 184 115;  % #3ab873
     39 164  95;  % #27a45f
     29 151  84;  % #1d9754
     16 138  73;  % #108a49
      0 125  62;  % #007d3e
      0 115  57;  % #007339
      0 105  52;  % #006934
      0  96  48;  % #006030
      0  86  43   % #00562b
    ];
    for t = {'thumb', 'index', 'middle', 'little'}
        x = 1:size(cmap.(t{1}), 1);
        xi = linspace(1, size(cmap.(t{1}), 1), ncolor);
        cmap.(t{1}) = interp1(x, double(cmap.(t{1})), xi, 'linear');
    end
end