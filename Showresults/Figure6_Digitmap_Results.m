%% This figure is shown in Figure 6.
%%
clc; clear; close all;
ft_defaults;

%%
FN = {'Allband','Delta','Theta','Alpha','Beta','Gamma','High-gamma','LowFreqBand'};
FB = {[0.5 90],[0.5 4],[4 8],[8 13],[13 30],[30 60],[60 90],[0.5 8]};

for f = 2:3  % delta and theta
    %% cal Group - all sub intersection
    opt = [];
    opt.filepath = fileparts(pwd);
    opt.freqname = FN{f};
    opt.freqband = FB{f};
    opt.sub = 'S00';
    opt.paradigm = 'Extension_4tasks';
    opt.motion = 'Extension';
    opt.ana = 'Source_Space';
    opt.subana = 'DigitMap';
    opt.mripath = fullfile(opt.filepath, 'mat', 'other', opt.sub, 'forMEG');
    opt.task = {'Thumb','Index','Middle','Little'};
    opt.groupmatdir = fullfile(opt.filepath, 'mat', opt.ana, opt.subana);
    opt.groupresultdir = fullfile(opt.filepath, 'results', opt.ana, opt.subana);

    if ~exist(opt.groupresultdir, 'dir')
        mkdir(opt.groupresultdir);
    end

    timewin = 0.05;
    timestep = 0.05;
    Timeseg = -1:timestep:1.95;

    Subidx = [1:3 19 5:9 11 14:18 10];

    %%
    if f == 2
        load(fullfile(opt.groupmatdir, strcat(opt.freqname, '_DigitMAP_4tasks_Group_ACTseg_-1-2s_', num2str(numel(Subidx)), 'sub.mat')));
    elseif f == 3
        load(fullfile(opt.groupmatdir, strcat(opt.freqname, '_DigitMAP_4tasks_Group_ACTseg_-1-1_5s_', num2str(numel(Subidx)), 'sub.mat')));
    end

    %%
    brainsurf = struct('surf_8k', struct(), 'surf_32k', struct());
    s = {'8k'};
    for h = {'L', 'R'}
        brainsurf.(strcat('surf_', s{1})).(h{1}).inflated = ft_read_headshape(fullfile(opt.filepath, 'mat', 'other', 'S00', s{1}, sprintf('template.%s.inflated.%s_fs_LR.surf.gii', h{1}, s{1})));
        brainsurf.(strcat('surf_', s{1})).(h{1}).flat = ft_read_headshape(fullfile(opt.filepath, 'mat', 'other', 'S00', s{1}, sprintf('template.%s.flat.%s_fs_LR.surf.gii', h{1}, s{1})));
        brainsurf.(strcat('surf_', s{1})).(h{1}).label = ft_read_atlas(fullfile(opt.filepath, 'mat', 'other', 'S00', s{1}, sprintf('template.%s.aparc.%s_fs_LR.label.gii', h{1}, s{1})));
        brainsurf.(strcat('surf_', s{1})).(h{1}).sulc = double(gifti(fullfile(opt.filepath, 'mat', 'other', 'S00', s{1}, sprintf('template.%s.sulc.%s_fs_LR.shape.gii', h{1}, s{1}))).cdata);
    end
    filename = fullfile(opt.filepath, 'mat', 'other', 'S00', '8k', 'template.L.inflated.8k_fs_LR.surf.gii');
    brainsurf_8k = ft_read_headshape({filename strrep(filename, '.L.', '.R.')});
    cbasebrain = [229 229 229]/255;
    cprecentral = [212 185 218]/255;
    cpostcentral = [223 101 176]/255;

    %%
    color2PZ1 = [nan nan nan; 59 135 199; 242 120 115; 255 211 115; 54 151 88]/255;
    fig = figure('Units', 'centimeters', 'Position', [10 10 4 4], 'IntegerHandle', 'off');
    tiledlayout(1, 1, 'TileSpacing', 'tight', 'Padding', 'tight');
    nexttile; hold on;
    ft_plot_mesh(brainsurf_8k, 'facecolor', cbasebrain);
    digit_map = ft_plot_mesh(brainsurf_8k, 'vertexcolor', brainsurf_8k.sulc);
    view(-120, 40);
    material dull; lighting gouraud; camlight;
    set(gca,'FontSize', 10, 'FontName', 'Times New Roman')
    set(gca, 'PlotBoxAspectRatio', [1 1 1]); % 所有图窗大小一致
    twin = {[16:35], [19:32]};  % delta, theta
    for fband = {'delta', 'theta'}
        if ~exist(fullfile(opt.groupresultdir, fband{1}), 'dir')
            mkdir(fullfile(opt.groupresultdir, fband{1}));
        end
        if strcmp(fband{1}, 'delta')
            cidx = twin{1};
            load(fullfile(opt.groupmatdir, sprintf('%s_DigitMAP_4tasks_Group_ACTseg_-1-2s_%dsub.mat', fband{1}, numel(Subidx))));
        elseif strcmp(fband{1}, 'theta')
            cidx = twin{2};
            load(fullfile(opt.groupmatdir, sprintf('%s_DigitMAP_4tasks_Group_ACTseg_-1-1_5s_%dsub.mat', fband{1}, numel(Subidx))));
        end
        for c = cidx   % 16:35 %delta   % 19:32 % theta
            fprintf('Plotting %s time %.2fs\n', fband{1}, Timeseg(c));
            acc2plot = Map_4taskin1brain(c, :) + 1; % color2PZ1第一行为nan，实际颜色索引+1
            set(digit_map, 'FaceVertexCData', color2PZ1(acc2plot,:));
            set(gcf, 'Name', sprintf('%s_time_%.2fs', fband{1}, Timeseg(c)));
            drawnow;
            % pause(0.1);
            print(fig, fullfile(opt.groupresultdir, fband{1}, sprintf('%s_time_%.2fs_DigitMAP.png', fband{1}, Timeseg(c))), '-dpng', '-r600', '-image');
        end
    end
end

%% joy plot
toi = 11:41;
cm = [59 135 199; 242 120 115; 255 211 115; 54 151 88]/255;  % final color
for fband = {'delta', 'theta'}
    if ~exist(fullfile(opt.groupresultdir, fband{1}), 'dir')
        mkdir(fullfile(opt.groupresultdir, fband{1}));
    end
    if strcmp(fband{1}, 'delta')
        ylimits = [0 300];
        rectx   = [-0.25 0.75];
    elseif strcmp(fband{1}, 'theta')
        ylimits = [0 300];
        rectx   = [-0.1 0.55];
    end
    load(fullfile(opt.groupmatdir, sprintf('%s_DigitMAP_vertexnum_per_digit.mat', fband{1})));
    fprintf('Max digit num in %s: %d\n', fband{1}, max(digitnum(:, toi), [], 'all'));
    fig = figure('Units', 'centimeters', 'Position', [10 10 8 12], 'IntegerHandle', 'off');
    tiledlayout(4, 1, 'TileSpacing', 'tight', 'Padding', 'tight');
    for t = 1:4
        nexttile; hold on;
        rectangle('Position', [rectx(1) 0 rectx(2)-rectx(1) ylimits(2)], 'FaceColor', 0.9*[1 1 1], 'EdgeColor', 'none');
        plot([-0.5 1], numthres*[1 1], 'k--', 'LineWidth', 1);
        x = Timeseg;
        surf_x = [x; x];
        surf_y = [digitnum(t,:); zeros(1, length(x))];
        surf_z = zeros(size(surf_x));
        surf(surf_x, surf_y, surf_z, [digitnum(t,:); digitnum(t,:)], 'EdgeColor', 'none', 'FaceColor', cm(t,:), 'FaceAlpha', 0.3);
        plot(Timeseg, digitnum(t,:), 'Color', cm(t,:), 'LineWidth', 1.2);
        xlim([-0.5 1.0]);
        ylim(ylimits);
        if t ~= 4
            set(gca, 'XTickLabel', []);
        else
            set(gca, 'XTick', -0.5:0.5:1.0);
            xlabel('Time (s)', 'FontName', 'Arial', 'FontSize', 8);
        end
        set(gca, 'YTick', []);
        ylabel(opt.task{t});
        set(gca, 'TickDir', 'out', 'Color', 'none');
        set(fig, 'Name', sprintf('DigitNum_%s_task_%s', fband{1}, opt.task{t}));
    end
    print(fig, fullfile(opt.groupresultdir, fband{1}, strcat(fig.Name, '.png')), '-dpng', '-r600', '-vector');
    % print(fig, fullfile(opt.groupresultdir, fband{1}, strcat(fig.Name, '.eps')), '-depsc2', '-r600', '-vector');
end
