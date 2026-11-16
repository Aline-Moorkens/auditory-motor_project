clc
% no "clear" / "close all" here - this script re-uses condOrder, plotSNR, nSubj, sigBandLo, sigBandHi
% that are already in the workspace from running the spaghettiplot script first.

% run first spaghettiplot script!

% LMM results!
p_Foam_vs_NoFoam  = 0.0001;
p_Foam_vs_Sit = 0.0001;
p_NoFoam_vs_Sit = 0.126;
alphaBonf = 0.05;


%%  Boxplot: one box per condition (Foam / Firm / Sit), subjects
%   overlaid as gray dots and connected subject-by-subject with gray
%   trajectory lines. Colors are PER CONDITION only
boxKeyOrder   = {'Foam','Firm','Sit'};        % display order, left to right
boxDispLabel  = {'Foam','Firm','Sit'};
nBox          = numel(boxKeyOrder);
 
boxColors = [224 159   0;           % Foam - amber/gold
             199  21 133;           % Firm - magenta/crimson
             106  27 154] / 255;    % Sit - purple/indigo
 
colDot  = [0.6 0.6 0.6];   % gray individual-subject dots
colLine = [0.6 0.6 0.6];   % gray subject trajectory lines


% map display keys to columns of plotSNR (condOrder = SIT,FIRM,FOAM,SITREPEAT)
colIdx = nan(1, nBox);
for b = 1:nBox
    colIdx(b) = find(strcmpi(condOrder, boxKeyOrder{b}));
end
boxVals = plotSNR(:, colIdx);   % nSubj x nBox, in display order


sText = 18;
figure('Color', 'w');
hold on; box on


% --- boxplots (median / IQR), one per condition ---
pBox = gobjects(1, nBox);
for b = 1:nBox
    vals = boxVals(:, b);
    vals = vals(~isnan(vals));
    bc = boxchart(b * ones(size(vals)), vals, ...
        'BoxFaceColor', boxColors(b,:), 'BoxFaceAlpha', 0.35, ...
        'BoxEdgeColor', boxColors(b,:), 'WhiskerLineColor', boxColors(b,:), ...
        'MarkerStyle', 'none', 'BoxWidth', 0.5);
    pBox(b) = bc;
end


% --- individual subject dots, gray ---
jitterWidth = 0.12;
pDot = gobjects(1, 1);
for b = 1:nBox
    vals = boxVals(:, b);
    keep = ~isnan(vals);
    %xj = b + (rand(sum(keep), 1) - 0.5) * 2 * jitterWidth;
    % hDot = scatter(xj, vals(keep), 30, colDot, 'filled', ...
    %     'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', 'none');
    hDot = scatter(b, vals(keep), 30, colDot, 'filled', ...
        'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', 'none');
    if b == 1, pDot = hDot; end
end
 
% --- subject trajectory lines across the 3 conditions (gray) ---
pLine = gobjects(1, 1);
firstLine = true;
for s = 1:nSubj
    vals = boxVals(s, :);
    if any(isnan(vals))
        continue   % only draw a trajectory where this subject has all 3 conditions
    end
    hLine = plot(1:nBox, vals, '-', 'Color', [colLine 0.5], 'LineWidth', 0.8, ...
        'Marker', 'none');
    if firstLine
        pLine = hLine;
        firstLine = false;
    end
end
 
% --- axes / cosmetics ---
xticks(1:nBox)
xticklabels(boxDispLabel)
xlim([0.4, nBox + 0.6])
 
xlabel('Condition')
ylabel('FFR SNR')
%ylabel('FFR SNR (signal power / noise power)')

 
set(gca, 'fontsize', sText)
% title(sprintf('FFR SNR [%d-%d Hz], subject-level (n = %d)', ...
%     sigBandLo, sigBandHi, sum(~all(isnan(boxVals), 2))))
title(sprintf('FFR SNR [%d-%d Hz], (N = %d)', ...
    sigBandLo, sigBandHi, sum(~all(isnan(boxVals), 2))),"FontSize",22)

% % Dedicated, neutral proxy handles for the legend only - not tied to any
% % one condition's color, and not to a single subject's actual line.
% pBoxProxy  = patch(nan, nan, 'w', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', [0.6 0.6 0.6], 'LineWidth', 1);
% pDotProxy  = scatter(nan, nan, 30, colDot, 'filled', 'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', 'none');
% pLineProxy = plot(nan, nan, '-', 'Color', colLine, 'LineWidth', 0.8);
% 
% legend([pBoxProxy, pDotProxy, pLineProxy], ...
%     {'Boxplot (median / IQR)', 'Individual subjects', 'Subject trajectory'}, ...
%     'AutoUpdate', 'off', ...
%     'Location', 'eastoutside', 'Orientation', 'vertical', 'Box', 'off')


%%% Significance annotation (Bonferroni-corrected)
yMax   = max(boxVals(:));
yMin   = min(boxVals(:));
yRange = yMax - yMin;
if yRange == 0
    yRange = 1;
end
yStep  = 0.02 * yRange;
startY = yMax + 0.01 * yRange;

ylim([yMin-5*yStep, yMax + 10*yStep]);

% 1) STAND_FOAM (1) vs STAND_NO_FOAM (2)
if p_Foam_vs_NoFoam <= alphaBonf
    y1 = startY;
    plot([1 1 2 2], [y1 y1+yStep y1+yStep y1], 'k-', 'LineWidth', 1);
    if p_Foam_vs_NoFoam <= 0.001
        stars = '***';
    elseif p_Foam_vs_NoFoam <= 0.01
        stars = '**';
    elseif p_Foam_vs_NoFoam <= 0.05
        stars = '*';
    else
        stars = '';
    end
    text(1.5, y1 + 0.1*yStep, stars, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment',   'bottom', ...
        'FontSize', sText);
end

% 2) STAND_FOAM (1) vs SIT (3)
if p_Foam_vs_Sit <= alphaBonf
    y2 = startY + 2.8*yStep;
    plot([1 1 3 3], [y2 y2+yStep y2+yStep y2], 'k-', 'LineWidth', 1);
    if p_Foam_vs_Sit <= 0.001
        stars = '***';
    elseif p_Foam_vs_Sit <= 0.01
        stars = '**';
    elseif p_Foam_vs_Sit <= 0.05
        stars = '*';
    else
        stars = '';
    end
    text(2, y2 + 0.1*yStep, stars, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment',   'bottom', ...
        'FontSize', sText);
end

% 3) STAND_NO_FOAM (2) vs SIT (3)
if p_NoFoam_vs_Sit <= alphaBonf
    y3 = startY + 3.5*yStep;
    plot([2 2 3 3], [y3 y3+yStep y3+yStep y3], 'k-', 'LineWidth', 1);
    if p_NoFoam_vs_Sit <= 0.001
        stars = '***';
    elseif p_NoFoam_vs_Sit <= 0.01
        stars = '**';
    elseif p_NoFoam_vs_Sit <= 0.05
        stars = '*';
    else
        stars = '';
    end
    text(2.5, y3 + 0.1*yStep, stars, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment',   'bottom', ...
        'FontSize', sText);
end


hold off

%% Save figure
dpi = 300;
exportTransparentPNG(fig1, 'SNR_333_boxplot.png', dpi);


%% Local function: export figure as PNG with transparent background + correct DPI metadata
function exportTransparentPNG(fig, filename, dpi)
tmpFile = [tempname, '.png'];
exportgraphics(fig, tmpFile, 'Resolution', dpi);

img = imread(tmpFile);
isWhite = all(img == 255, 3);
alphaMask = double(~isWhite);

pxPerMeter = round(dpi / 0.0254);

imwrite(img, filename, ...
    'Alpha', alphaMask, ...
    'ResolutionUnit', 'meter', ...
    'XResolution', pxPerMeter, ...
    'YResolution', pxPerMeter);

delete(tmpFile);
end