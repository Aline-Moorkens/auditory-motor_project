%% Compare sit - stand - foam for PLV values
clc, clear
close all

%% Set parameters
subjectid = 'TH';
alpha = 0.05;
conditions = {'sit','firm','foam'};

% input
inputpath = fullfile('F:\Auditory-mov-data\control-no_plate\phase-analysis\plv_subsampled');


%% load plv data
PLV = load(fullfile(inputpath, 'PLV_table.mat'));
PLV = PLV.T_plv;

subT_sit = PLV(PLV.Subject == subjectid & PLV.Condition == 'sit' & PLV.Session == 'pre',:);
subT_firm = PLV(PLV.Subject == subjectid & PLV.Condition == 'firm' & PLV.Session == 'pre',:);
subT_foam = PLV(PLV.Subject == subjectid & PLV.Condition == 'foam' & PLV.Session == 'pre',:);

PLV_sit.plv = subT_sit{:,"PLV"}; PLV_sit.labels = subT_sit{:,"Channel"};
PLV_firm.plv = subT_firm{:,"PLV"}; PLV_firm.labels = subT_firm{:,"Channel"};
PLV_foam.plv = subT_foam{:,"PLV"}; PLV_foam.labels = subT_foam{:,"Channel"};

%% Only include channels with significant FFR
% non-siginificant chans are already left out in the previous step (PLV
% table computation)


%% Non-parametric test for 3 related conditions: Friedman test
plv_sit = PLV_sit.plv; plv_firm = PLV_firm.plv; plv_foam = PLV_foam.plv; 

data = [plv_foam, plv_firm, plv_sit];

[p_friedman, ~, ~] = friedman(data, 1, 'off');

disp('Friedman test p-value:')
disp(p_friedman)


%% Post-hoc pairwise comparisons (Wilcoxon signed-rank) + Bonferroni
[p_Foam_vs_NoFoam, h_Foam_vs_NoFoam, stats_Foam_vs_NoFoam] = ...
    signrank(plv_foam, plv_firm);

[p_Foam_vs_Sit, h_Foam_vs_Sit, stats_Foam_vs_Sit] = ...
    signrank(plv_foam, plv_sit);

[p_NoFoam_vs_Sit, h_NoFoam_vs_Sit, stats_NoFoam_vs_Sit] = ...
    signrank(plv_firm, plv_sit);

m         = numel(conditions);
alphaBonf = alpha/m;

p_raw  = [p_Foam_vs_NoFoam, p_Foam_vs_Sit, p_NoFoam_vs_Sit];
p_corr = min(p_raw * m, 1);

p_Foam_vs_NoFoam_corr = p_corr(1);
p_Foam_vs_Sit_corr    = p_corr(2);
p_NoFoam_vs_Sit_corr  = p_corr(3);

disp('Pairwise Wilcoxon p-values (raw and Bonferroni-corrected):')
disp(table( ...
    p_Foam_vs_NoFoam', p_Foam_vs_NoFoam_corr', ...
    p_Foam_vs_Sit',    p_Foam_vs_Sit_corr', ...
    p_NoFoam_vs_Sit',  p_NoFoam_vs_Sit_corr', ...
    'VariableNames', {'Foam_vs_NoFoam_raw','Foam_vs_NoFoam_corr', ...
                      'Foam_vs_Sit_raw','Foam_vs_Sit_corr', ...
                      'NoFoam_vs_Sit_raw','NoFoam_vs_Sit_corr'}))
fprintf('Bonferroni-corrected alpha: %.4f\n', alphaBonf);


%% Boxplot visualisation (3 conditions)
nr_conditions = numel(conditions);
nr_chans = size(plv_sit,1);
bp = zeros([nr_chans, nr_conditions]);
bp(:,1) = plv_foam;
bp(:,2) = plv_firm;
bp(:,3) = plv_sit;

sText = 14;
figure;
h = boxplot(bp, 'Notch','off', ...
              'Labels', {'STAND_FOAM','STAND_NO_FOAM','SIT'});

set(gca,'fontsize', sText);
set(h,{'linew'},{1});

c = findobj(gca,'Tag','Box');
patch(get(c(3),'XData'), get(c(3),'YData'), [0.4860 0.6940 0.1680], 'FaceAlpha', .26);
patch(get(c(2),'XData'), get(c(2),'YData'), [0.9290 0.6940 0.1250], 'FaceAlpha', .26);
patch(get(c(1),'XData'), get(c(1),'YData'), [0.8500 0.3250 0.0980], 'FaceAlpha', .26);

if p_friedman < 0.0001
    p_str = '<0.0001';
else
    p_str = sprintf('= %.4f', p_friedman);
end

title(['FFR PLV comparison [331–335 Hz] ',subjectid])
ylabel('PLV-value')

hold on
% parallelcoords(bp, 'Color', 0.7*[1 1 1], ...
%                    'Linestyle', '-', 'Marker', '.', 'Markersize', 7);


% Annotate vertex electrodes on boxplot
vertex_labels = {'Cz','C1','C2','CPz','FC1','FC2','CP1','CP2'};
labels = PLV_sit.labels;

% Logical/numeric mask: 1 where labels matches any entry in vertex_labels, 0 otherwise
mask = double(ismember(labels, vertex_labels));

% Same-size array of category names based on the mask
category = repmat({'non-vertex'}, size(mask));
category(mask == 1) = {'vertex'};


hold on
h = parallelcoords(bp, 'Color', 0.7*[1 1 1], ...
    'Linestyle', '-', 'Marker', '.', 'Markersize', 7);

% h is a vector of line handles, one per row/observation in bp
% Recolor/thicken the ones corresponding to vertex electrodes, and label them
for i = 1:numel(h)
    if mask(i) == 1
        set(h(i), 'Color', [1 0 0], 'LineWidth', 1);   % vertex: red, thick

        % Get this line's x/y data to find where to place the label
        xdata = get(h(i), 'XData');
        ydata = get(h(i), 'YData');

        % Place label just past the last point of the line
        text(xdata(end) + 0.1, ydata(end), labels(i), ...
            'Color', [1 0 0], 'FontSize', 8, 'FontWeight', 'bold');
    end
end

%% Significance annotation (Bonferroni-corrected)
yMax   = max(bp(:));
yMin   = min(bp(:));
yRange = yMax - yMin;
if yRange == 0
    yRange = 1;
end
yStep  = 0.02 * yRange;
startY = yMax + 0.01 * yRange;

ylim([0, yMax + 7*yStep]);

% 1) STAND_FOAM (1) vs STAND_NO_FOAM (2)
if p_Foam_vs_NoFoam_corr <= alphaBonf
    y1 = startY;
    plot([1 1 2 2], [y1 y1+yStep y1+yStep y1], 'k-', 'LineWidth', 1);
    if p_Foam_vs_NoFoam_corr <= 0.001
        stars = '***';
    elseif p_Foam_vs_NoFoam_corr <= 0.01
        stars = '**';
    elseif p_Foam_vs_NoFoam_corr <= 0.05
        stars = '*';
    else
        stars = '';
    end
    text(1.5, y1 + 1.1*yStep, stars, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment',   'bottom', ...
        'FontSize', sText);
end

% 2) STAND_FOAM (1) vs SIT (3)
if p_Foam_vs_Sit_corr <= alphaBonf
    y2 = startY + 1.6*yStep;
    plot([1 1 3 3], [y2 y2+yStep y2+yStep y2], 'k-', 'LineWidth', 1);
    if p_Foam_vs_Sit_corr <= 0.001
        stars = '***';
    elseif p_Foam_vs_Sit_corr <= 0.01
        stars = '**';
    elseif p_Foam_vs_Sit_corr <= 0.05
        stars = '*';
    else
        stars = '';
    end
    text(2, y2 + 1.1*yStep, stars, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment',   'bottom', ...
        'FontSize', sText);
end

% 3) STAND_NO_FOAM (2) vs SIT (3)
if p_NoFoam_vs_Sit_corr <= alphaBonf
    y3 = startY + 3.2*yStep;
    plot([2 2 3 3], [y3 y3+yStep y3+yStep y3], 'k-', 'LineWidth', 1);
    if p_NoFoam_vs_Sit_corr <= 0.001
        stars = '***';
    elseif p_NoFoam_vs_Sit_corr <= 0.01
        stars = '**';
    elseif p_NoFoam_vs_Sit_corr <= 0.05
        stars = '*';
    else
        stars = '';
    end
    text(2.5, y3 + 1.1*yStep, stars, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment',   'bottom', ...
        'FontSize', sText);
end



%% Only vertex electrodes included
% Keep only the vertex electrode rows
bp_vertex = bp(mask == 1, :);
plv_sit_vertex = bp(:,1); plv_firm_vertex = bp(:,2); plv_foam_vertex = bp(:,3);
% Non-parametric test for 3 related conditions: Friedman test

data = [plv_foam_vertex, plv_firm_vertex, plv_sit_vertex];

[p_friedman, tbl_friedman, stats_friedman] = friedman(data, 1, 'off');

disp('Friedman test p-value:')
disp(p_friedman)


% Post-hoc pairwise comparisons (Wilcoxon signed-rank) + Bonferroni
[p_1, h_1, stats_1] = ...
    signrank(plv_foam_vertex, plv_firm_vertex);

[p_2, h_2, stats_2] = ...
    signrank(plv_foam_vertex, plv_sit_vertex);

[p_3, h_3, stats_3] = ...
    signrank(plv_firm_vertex, plv_sit_vertex);

m         = numel(conditions);
alphaBonf = alpha/m;

p_raw_vertex  = [p_1, p_2, p_3];
p_corr_vertex = min(p_raw * m, 1);

p_1_corr = p_corr_vertex(1);
p_2_corr = p_corr_vertex(2);
p_3_corr = p_corr_vertex(3);

disp('Pairwise Wilcoxon p-values (raw and Bonferroni-corrected):')
disp(table( ...
    p_Foam_vs_NoFoam', p_1_corr', ...
    p_Foam_vs_Sit',    p_2_corr', ...
    p_NoFoam_vs_Sit',  p_3_corr', ...
    'VariableNames', {'Foam_vs_NoFoam_raw','Foam_vs_NoFoam_corr', ...
                      'Foam_vs_Sit_raw','Foam_vs_Sit_corr', ...
                      'NoFoam_vs_Sit_raw','NoFoam_vs_Sit_corr'}))
fprintf('Bonferroni-corrected alpha: %.4f\n', alphaBonf);

% Boxplot


sText = 14;

figure;
h = boxplot(bp_vertex, 'Notch','off', ...
    'Labels', {'STAND_FOAM','STAND_NO_FOAM','SIT'});

set(gca,'fontsize', sText);
set(h,{'linew'},{1});

c = findobj(gca,'Tag','Box');

patch(get(c(3),'XData'), get(c(3),'YData'), [0.4860 0.6940 0.1680], 'FaceAlpha', .26);
patch(get(c(2),'XData'), get(c(2),'YData'), [0.9290 0.6940 0.1250], 'FaceAlpha', .26);
patch(get(c(1),'XData'), get(c(1),'YData'), [0.8500 0.3250 0.0980], 'FaceAlpha', .26);

if p_friedman < 0.0001
    p_str = '<0.0001';
else
    p_str = sprintf('= %.4f', p_friedman);
end

title(['FFR PLV comparison [331–335 Hz] - Vertex electrodes ',subjectid])
ylabel('PLV-value')
hold on
parallelcoords(bp_vertex, 'Color', 0.7*[1 1 1], ...
                    'Linestyle', '-', 'Marker', '.', 'Markersize', 7);

% Significance annotation (Bonferroni-corrected)
yMax   = max(bp_vertex(:));
yMin   = min(bp_vertex(:));
yRange = yMax - yMin;
if yRange == 0
    yRange = 1;
end

yStep  = 0.02 * yRange;
startY = yMax + 0.01 * yRange;

ylim([yMin - 5*yStep, yMax + 7*yStep]);

% 1) STAND_FOAM (1) vs STAND_NO_FOAM (2)
if p_1_corr <= alphaBonf
    y1 = startY;
    plot([1 1 2 2], [y1 y1+yStep y1+yStep y1], 'k-', 'LineWidth', 1);
    if p_Foam_vs_NoFoam_corr <= 0.001
        stars = '***';
    elseif p_Foam_vs_NoFoam_corr <= 0.01
        stars = '**';
    elseif p_Foam_vs_NoFoam_corr <= 0.05
        stars = '*';
    else
        stars = '';
    end
    text(1.5, y1 + 1.1*yStep, stars, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment',   'bottom', ...
        'FontSize', sText);
end

% 2) STAND_FOAM (1) vs SIT (3)
if p_2_corr <= alphaBonf
    y2 = startY + 1.6*yStep;
    plot([1 1 3 3], [y2 y2+yStep y2+yStep y2], 'k-', 'LineWidth', 1);
    if p_Foam_vs_Sit_corr <= 0.001
        stars = '***';
    elseif p_Foam_vs_Sit_corr <= 0.01
        stars = '**';
    elseif p_Foam_vs_Sit_corr <= 0.05
        stars = '*';
    else
        stars = '';
    end
    text(2, y2 + 1.1*yStep, stars, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment',   'bottom', ...
        'FontSize', sText);
end

% 3) STAND_NO_FOAM (2) vs SIT (3)
if p_3_corr <= alphaBonf
    y3 = startY + 3.2*yStep;
    plot([2 2 3 3], [y3 y3+yStep y3+yStep y3], 'k-', 'LineWidth', 1);
    if p_NoFoam_vs_Sit_corr <= 0.001
        stars = '***';
    elseif p_NoFoam_vs_Sit_corr <= 0.01
        stars = '**';
    elseif p_NoFoam_vs_Sit_corr <= 0.05
        stars = '*';
    else
        stars = '';
    end
    text(2.5, y3 + 1.1*yStep, stars, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment',   'bottom', ...
        'FontSize', sText);
end
