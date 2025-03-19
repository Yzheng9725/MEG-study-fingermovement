%%
% Finger extension movement
% Fig.5
% creator: WxyZ - Yu Zheng
% Date: 20250314

%%
ft_defaults;

baseline = [-1.5 -1.0];

% Colormap for 4 tasks  Blue-Red-Yellow-Green -> Thumb-Index-Middle-Little
cm = [59 135 199; 242 120 115; 255 211 115; 54 151 88]/255;

%% Source pos - select ROI & set colormap
% need MNI152 template
opt.mripath = '.\mat\S00\forMEG\';
opt.wbpath = '.\mat\S00\workbench\';

load([opt.mripath 'headmodel.mat'])
load([opt.mripath 'sourcemodel.mat'])
load([opt.mripath 'sourcemodel_inflated.mat'])
load([opt.mripath 'brainatlas.mat'])

atlas_wb_L = brainatlasL;
atlas_wb_R = brainatlasR;

atlasLabel = [21 23];
Atlasannot_Lshpere = [atlas_wb_L.parcellation];   % brain area label in atlas, only left

% select pre/post-all ROI pos index
ROIidx = [];
for i = 1:numel(atlasLabel)
    idx = find(Atlasannot_Lshpere == atlasLabel(i));
    ROIidx = [ROIidx;idx];
    clearvars idx
end

ROIidx = sort(ROIidx,'ascend');

%% fit a curved surface to find the projection of average source pos
% need s2009s
atlas_wb_L = ft_read_atlas([opt.wbpath 'S00.L.aparc.a2009s.8k_fs_LR.label.gii']);
atlas_wb_R = ft_read_atlas([opt.wbpath 'S00.R.aparc.a2009s.8k_fs_LR.label.gii']);
atlasLabel = [29 30 46];        % label:Postcentral/Precentral   aparc

% select pre/post-segment ROI pos index
[~, roiL_init, ~] = wxyz_selectbrainroi(atlas_wb_L, atlas_wb_R, 4, ...
    'sourcemodel',sourcemodel, 'isextend', true, 'roi4extend', atlasLabel, 'iteration', 1);

[~, roiL_end, ~] = wxyz_selectbrainroi(atlas_wb_L, atlas_wb_R, 4, ...
    'sourcemodel',sourcemodel, 'isextend', true, 'roi4extend', atlasLabel, 'iteration',20);

ROIidx_fit = setdiff(roiL_end, roiL_init);

figure
T = sourcemodel_inflated.pos(ROIidx_fit,:);
f = fit([T(:,1), T(:,2)],T(:,3),'linearinterp');
plot(f, [T(:,1), T(:,2)], T(:,3));

view(-100,45)

%% if need to plot on contral-sensorimotor cortex
% Only plot pre/post-central
sourcemodel_L = sourcemodel_inflated;
sourcemodel_L.pos = sourcemodel_inflated.pos(sourcemodel.brainstructure == 1,:);
sourcemodel_L.tri = sourcemodel_inflated.tri(1:end/2,:);
sourcemodel_L.sulc = sourcemodel_inflated.sulc(sourcemodel_inflated.brainstructure == 1,:);
sourcemodel_L.curv = sourcemodel_inflated.curv(sourcemodel_inflated.brainstructure == 1,:);
sourcemodel_L.thickness = sourcemodel_inflated.thickness(sourcemodel_inflated.brainstructure == 1,:);
sourcemodel_L.brainstructure = sourcemodel_inflated.brainstructure(sourcemodel_inflated.brainstructure == 1,:);
sourcemodel_L.inside = sourcemodel.inside(sourcemodel_inflated.brainstructure == 1,:);

sulc2plot = sourcemodel_L.sulc;

% create a circle
circle_center = [-42.8448 -23.4916 68.5311];   % avg source pos of LITTLE
radius = 18;  % plot in a bigger size 
vertices = sourcemodel_inflated.pos;
distances = sqrt(sum((vertices - circle_center).^2, 2)); % cal the distance betwwen very vertex and the center
circle_indices = distances <= radius; % obtain the index of the vertex in the circle
highlighted_vertices = vertices(circle_indices, :); % obtain the pos of the vertex in the circle
sulc2plot(~circle_indices(1:7842)) = NaN;  % only left
sulc2plot = -1*sulc2plot;  % to create the colormap


%%  fig.5 whole brain
h = figure('Units','centimeters','Position',[10 10 8 8], 'IntegerHandle', 'off');
tiledlayout(1, 1, 'TileSpacing','compact','Padding','tight');
nexttile
ft_plot_mesh(sourcemodel_inflated,'vertexcolor', -1*sourcemodel_inflated.sulc);
view(-100, 45);
colormap('gray');
clim([-20 10]);


%% Fig.5(a)
%% plot all sub all pos
opt.freqname = 'LowFreqBand';
load(['.\mat\SourceLevel\' opt.freqname '\Source_Group.mat'])

%%
circle_indices = distances <= 16;  % select the pos in the 16 mm-radius circle
vidx = find(circle_indices == 1);
TitleName = {'t1','t2','t3','t4'};

for l = 1:4
    %%
    h = figure('Units','centimeters','Position',[10 10 12 12], 'IntegerHandle', 'off');
    tiledlayout(1, 1,  'TileSpacing','compact','Padding','compact');  % ,'TileIndexing','columnmajor'

    nexttile
    ft_plot_mesh(sourcemodel_L,'vertexcolor', sulc2plot,'facealpha',0.8); % lighting gouraud; material dull;light
    hold on
    view(-100, 50)
    colormap('gray')
    clim([-20 10])

    %% plot all pos in circle
    for ss = 1:16   % [1 4 6 8 9 15 17 18]  % [4 8 9 14 15 16 17 18]
        for t = 1:4
            drawnow
            if ismember(SourcePosidx_Group(ss,t,l),vidx)
                Pos2plot = sourcemodel_inflated.pos(squeeze(SourcePosidx_Group(ss,t,l)),:);
                scatter3(Pos2plot(1), Pos2plot(2),Pos2plot(3),10,cm(t,:),'o','filled','MarkerFaceAlpha',0.9);  % ,'MarkerFaceAlpha',0.8
                hold on
            end
        end
    end

    % plot pos-avg per time
    o = 1;
    for t = 1:4
        %% cal source pos avg
        for ss = 1:16
            Pos2plot_all(ss,:) = sourcemodel_inflated.pos(squeeze(SourcePosidx_Group(ss,t,l)),:);
        end
        Pos2plot = mean(Pos2plot_all);
        Posavg_pertime(o,:) =  Pos2plot;

        %% project to ROI-curved surface
        objective_func = @(params) norm([params(1), params(2), f(params(1), params(2))] - Pos2plot);
        initial_guess = [Pos2plot(1,1),Pos2plot(1,2)];

        % mini using fminunc
        options = optimset('Display', 'off');
        [optimal_params, ~] = fminunc(objective_func, initial_guess, options);
        Posavg_pertime_proj = [optimal_params(1), optimal_params(2), f(optimal_params(1), optimal_params(2))];

        scatter3(Posavg_pertime_proj(1), Posavg_pertime_proj(2),Posavg_pertime_proj(3),90,cm(t,:),'x','LineWidth',3);

        %%
        o = o+1;
    end

    zoom(3.5)
    title(TitleName{l},'Position',[0 -31 50],'FontName', 'Calibri','FontSize', 12, 'FontWeight', 'bold')

    %% Show all surf
    zlim([-20 80]);

end



%% Show average pos low frequency
load(['.\mat\SourceLevel\' opt.freqname '\POSonMNI_pertask_perpeak.mat']);
%%  plot
for ss = 1:16
    for t = 1:4
        posavg_pertask_16sub(ss,t,:) = mean(squeeze(Pos_all(ss,t,:,:)));
    end
end
circle_center = [-42.8448 -23.4916 68.5311];   % center -> LITTLE average pos

%%
h = figure('Units','centimeters','Position',[10 10 12 12], 'IntegerHandle', 'off');
tiledlayout(1, 1,  'TileSpacing','compact','Padding','compact');  % ,'TileIndexing','columnmajor'

nexttile
ft_plot_mesh(sourcemodel_L,'vertexcolor', sulc2plot,'facealpha',0.8); % lighting gouraud; material dull;light
hold on
view(-100, 50)
colormap('gray')
clim([-20 10])

% plot all pos in circle
for ss = 1:16
    for t = 1:4
        Pos2plot = squeeze(posavg_pertask_16sub(ss,t,:))';
        distances = sqrt(sum((Pos2plot - circle_center).^2, 2));
        drawnow
        if distances <= 16
            Vnum2plot_pertask(ss,t) = 1;
            %% project to ROI
            objective_func = @(params) norm([params(1), params(2), f(params(1), params(2))] - Pos2plot);
            initial_guess = [Pos2plot(1,1),Pos2plot(1,2)];
            options = optimset('Display', 'off');
            [optimal_params, ~] = fminunc(objective_func, initial_guess, options);
            Posavg_all_pertask_proj = [optimal_params(1), optimal_params(2), f(optimal_params(1), optimal_params(2))];
            scatter3(Posavg_all_pertask_proj(1), Posavg_all_pertask_proj(2),Posavg_all_pertask_proj(3),10,cm(t,:),'o','filled','MarkerFaceAlpha',0.9);  % ,'filled','MarkerFaceAlpha',0.5
            hold on
        else
            Vnum2plot_pertask(ss,t) = 0;
        end
    end
end

% avg per task
Posavg_pertask= squeeze(mean(posavg_pertask_16sub));
for n = 1:size(Posavg_pertask,1)
    objective_func = @(params) norm([params(1), params(2), f(params(1), params(2))] - Posavg_pertask(n,:));
    initial_guess = [Posavg_pertask(n,1),Posavg_pertask(n,2)];
    options = optimset('Display', 'off');
    [optimal_params, ~] = fminunc(objective_func, initial_guess, options);
    Posavg_pertask_proj(n,:) = [optimal_params(1), optimal_params(2), f(optimal_params(1), optimal_params(2))];
end
for t = 4:-1:1
    scatter3(Posavg_pertask_proj(t,1), Posavg_pertask_proj(t,2),Posavg_pertask_proj(t,3),80,cm(t,:),'x' ,'LineWidth',3)  % ,'filled','LineWidth',3,
end

zoom(3.5)
title('Average','Position',[0 -31 50],'FontName', 'Calibri','FontSize', 12, 'FontWeight', 'bold')
zlim([-20 80]);




%% Fig.5(b)
%% avg per task
o = 1;
for t = 1:4
    for l = 1:4
        %% cal source pos avg
        for ss = 1:16
            Pos2plot_all(ss,:) = sourcemodel_inflated.pos(squeeze(SourcePosidx_Group(ss,t,l)),:);
        end
        Pos2plot = mean(Pos2plot_all);
        Posavg_pertime(o,:) =  Pos2plot;
        o = o+1;
    end
end

Posavg_pertask(1,:) = mean(Posavg_pertime(1:4,:));
Posavg_pertask(2,:) = mean(Posavg_pertime(5:8,:));
Posavg_pertask(3,:) = mean(Posavg_pertime(9:12,:));
Posavg_pertask(4,:) = mean(Posavg_pertime(13:16,:));

Dist_fingers(1) = norm(Posavg_pertask(2,:)-Posavg_pertask(1,:));
Dist_fingers(2) = norm(Posavg_pertask(3,:)-Posavg_pertask(2,:));
Dist_fingers(3) = norm(Posavg_pertask(4,:)-Posavg_pertask(3,:));
Dist_fingers(4) = norm(Posavg_pertask(4,:)-Posavg_pertask(1,:));

for n = 1:size(Posavg_pertask,1)
    objective_func = @(params) norm([params(1), params(2), f(params(1), params(2))] - Posavg_pertask(n,:));
    initial_guess = [Posavg_pertask(n,1),Posavg_pertask(n,2)];
    options = optimset('Display', 'off');
    [optimal_params, ~] = fminunc(objective_func, initial_guess, options);
    Posavg_pertask_proj(n,:) = [optimal_params(1), optimal_params(2), f(optimal_params(1), optimal_params(2))];
end

Dist_fingers_proj(1) = norm(Posavg_pertask_proj(2,:)-Posavg_pertask_proj(1,:));
Dist_fingers_proj(2) = norm(Posavg_pertask_proj(3,:)-Posavg_pertask_proj(2,:));
Dist_fingers_proj(3) = norm(Posavg_pertask_proj(4,:)-Posavg_pertask_proj(3,:));


h = figure('Units','centimeters','Position',[10 10 12 12], 'IntegerHandle', 'off');
tiledlayout(1, 1,  'TileSpacing','compact','Padding','compact');  % ,'TileIndexing','columnmajor'

nexttile
ft_plot_mesh(sourcemodel_L,'vertexcolor', sulc2plot,'facealpha',0.8); % lighting gouraud; material dull;light
hold on
view(-100, 50)
colormap('gray')
clim([-20 10])

for t = 4:-1:1
    scatter3(Posavg_pertask_proj(t,1), Posavg_pertask_proj(t,2),Posavg_pertask_proj(t,3),80,cm(t,:),'x' ,'LineWidth',3)  % ,'filled','LineWidth',3,
end

zoom(3.5)
zlim([-20 80]);
title('Average','Position',[0 -31 50],'FontName', 'Calibri','FontSize', 12, 'FontWeight', 'bold')



%% plot all sub all pos - alpha/beta
opt.freqname = 'alpha';    % 'beta'; 'gamma'
load(['.\mat\SourceLevel\' opt.freqname '\Source_Group.mat'])

mar = {'^','s'};  % alpha beta
% mar = {'s'};  % gamma

%% plot all sub pos
circle_indices = distances <= radius; % obtain the index of the vertex in the circle

circle_indices = distances <= 16;
vidx = find(circle_indices == 1);

h = figure('Units','centimeters','Position',[10 10 12 12], 'IntegerHandle', 'off');
tiledlayout(1, 1,  'TileSpacing','compact','Padding','compact');
nexttile
ft_plot_mesh(sourcemodel_L,'vertexcolor', sulc2plot,'facealpha',0.8);
hold on
view(-100, 50)
colormap('gray')
clim([-20 10])

for l = 1:size(SourcePosidx_Group,3)
    % plot all pos in circle
    for ss = 1:16
        for t = 1:4
            if ismember(SourcePosidx_Group(ss,t,l),vidx)
                 Pos2plot = sourcemodel_inflated.pos(squeeze(SourcePosidx_Group(ss,t,l)),:);
                %% project to ROI
                objective_func = @(params) norm([params(1), params(2), f(params(1), params(2))] - Pos2plot);
                initial_guess = [Pos2plot(1,1),Pos2plot(1,2)];
                options = optimset('Display', 'off');
                [optimal_params, ~] = fminunc(objective_func, initial_guess, options);
                Posavg_all_pertask_proj = [optimal_params(1), optimal_params(2), f(optimal_params(1), optimal_params(2))];
                scatter3(Posavg_all_pertask_proj(1), Posavg_all_pertask_proj(2),Posavg_all_pertask_proj(3),20,cm(t,:),mar{l},'MarkerFaceAlpha',0.75, 'LineWidth',1);  % ,'filled','MarkerFaceAlpha',0.5
                hold on
            end
        end
    end
end

% avg per task
o = 1;
for t = 1:4
    for l = 1:size(SourcePosidx_Group,3)
        %% cal source pos avg
        for ss = 1:16
            Pos2plot_all(ss,:) = sourcemodel_inflated.pos(squeeze(SourcePosidx_Group(ss,t,l)),:);
        end
        Pos2plot = mean(Pos2plot_all);
        Posavg_pertime(o,:) =  Pos2plot;
        o = o+1;
    end
end

for n = 1:size(Posavg_pertime,1)
    objective_func = @(params) norm([params(1), params(2), f(params(1), params(2))] - Posavg_pertime(n,:));
    initial_guess = [Posavg_pertime(n,1),Posavg_pertime(n,2)];
    options = optimset('Display', 'off');
    [optimal_params, ~] = fminunc(objective_func, initial_guess, options);
    Posavg_pertime_proj(n,:) = [optimal_params(1), optimal_params(2), f(optimal_params(1), optimal_params(2))];
end

o = 1;
for t = 1:4
    for l = 1: size(SourcePosidx_Group,3)
        scatter3(Posavg_pertime_proj(o,1), Posavg_pertime_proj(o,2),Posavg_pertime_proj(o,3),90,cm(t,:),mar{l},'filled')  % ,'filled','LineWidth',3,
        o = o+1;
    end
end

zoom(3.5)
zlim([-20 80]);
title(['\' opt.freqname],'Position',[0 -31 50],'FontName', 'Calibri','FontSize', 12, 'FontWeight', 'bold')




