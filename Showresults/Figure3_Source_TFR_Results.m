%% This figure is shown in Figure 3.
%%
clc; clear; close all;
ft_defaults;

%%
opt = [];
opt.filepath = fileparts(pwd);
opt.sub = 'S00';
opt.mripath = fullfile('E:', 'MRI', opt.sub, 'forMEG');
opt.paradigm = 'Extension_4tasks';
opt.motion = 'Extension';
opt.ana = 'Source_Space';
opt.task = {'Thumb','Index','Middle','Little'};
opt.subana = 'Process';
opt.groupmatdir = fullfile(opt.filepath, 'mat', opt.ana, opt.subana);
opt.groupresultdir = fullfile(opt.filepath, 'results', opt.ana, opt.subana);

if ~exist(opt.groupresultdir, 'dir')
    mkdir(opt.groupresultdir);
end

%% Load data - All task average TFR data and TFCE data
load(fullfile(opt.groupmatdir, 'TFR_Group_AllType_AllTaskAvg.mat'));

%% Show results - All task average TFR and stat
typeName = {'TFR_full', 'TFR_evoked', 'TFR_induced'};
climits = [-6 6]; % for t-value
for typeidx = 1:numel(typeName)
    tf_tmp = TFR_Group_AllType_AllTaskAvg_Stat_TFCE.(typeName{typeidx});
    % Show results
    fig = figure('Units','centimeters','Position', [10 10 8 6], 'IntegerHandle', 'off');
    tiledlayout(1, 1, 'TileSpacing','compact','Padding','tight');
    nexttile; hold on;
    imagesc(tf_tmp.time, tf_tmp.freq, squeeze(tf_tmp.stat));
    contour(tf_tmp.time, tf_tmp.freq, squeeze(tf_tmp.mask), 1, 'Color', 'k', 'LineWidth', 1);
    xlim([-1.5 2]); ylim([1 90]);
    xlabel('Time (s)'); ylabel('Frequency (Hz)');
    set(gca, 'YDir', 'normal', 'Color', 'none', 'xtick', -1:1:2, 'ytick', 10:10:90, 'Box', 'on');
    colormap(flip(ft_colormap('RdGy')));
    colorbar('ticks', [-6 0 6], 'TickLabels', {'-6', '0', '6'});
    clim(climits);
    set(fig, 'Name', sprintf('Stat_TFCE_%s', typeName{typeidx}));
    drawnow;
    print(fig, fullfile(opt.groupresultdir, strcat(fig.Name, '.png')), '-dpng', '-r600', '-image');
    % print(fig, fullfile(opt.groupresultdir, strcat(fig.Name, '.eps')), '-depsc2', '-r600', '-vector');
end

%% Show brain roi
filename = fullfile(opt.filepath, 'mat', 'other', 'S00', '8k', 'template.L.inflated.8k_fs_LR.surf.gii');
brainsurf_8k = ft_read_headshape({filename strrep(filename, '.L.', '.R.')});
brainatlas_8k = importdata(fullfile(opt.filepath, 'mat', 'other', 'S00', 'forMEG', 'brainatlas.mat'));
cbasebrain = [229 229 229]/255;
croibrain  = [102 102 102]/255;
atlasLabel = [21 23]; % Postcentral/Precentral
roi = find(brainatlas_8k.brainatlasL.parcellation == atlasLabel(1) | brainatlas_8k.brainatlasL.parcellation == atlasLabel(2));
croi = repmat([nan nan nan], [size(brainsurf_8k.pos, 1) 1]);
croi(roi, :) = repmat(croibrain, [numel(roi) 1]);

roiedge = [];
for i = 1:size(brainsurf_8k.tri, 1)
    tmp = sum(ismember(brainsurf_8k.tri(i, :), roi));
    if tmp == 1
        tmpidx = find(ismember(brainsurf_8k.tri(i, :), roi));
    end
end

fig = figure('Units', 'centimeters', 'Position', [10 10 4 4], 'IntegerHandle', 'off');
tiledlayout(1, 1, 'TileSpacing','compact','Padding','tight');
tile = nexttile; hold on;
ft_plot_mesh(brainsurf_8k, 'vertexcolor', repmat(cbasebrain, [size(brainsurf_8k.pos, 1) 1]));
ft_plot_mesh(brainsurf_8k, 'vertexcolor', croi);
view(-120, 40);
material dull; lighting gouraud; camlight('headlight');
set(gcf, 'Name', 'BrainSurface_SensoryMotorArea');
print(fig, fullfile(opt.groupresultdir, 'BrainSurface_SensoryMotorArea.png'), '-dpng', '-r600', '-image');
% print(fig, fullfile(opt.groupresultdir, 'BrainSurface_SensoryMotorArea.eps'), '-depsc2', '-r600', '-image');

