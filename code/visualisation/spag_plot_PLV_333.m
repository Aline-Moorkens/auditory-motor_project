%% PLV 333 Hz - SUBJECT-LEVEL pairwise Wilcoxon tests + spaghetti plot
%
% Mirrors the FFR333 SNR subject-level spaghetti script, but reads
% PRE-COMPUTED PLV values (per subject/condition/channel) from an Excel
% table instead of computing anything from raw EEG.
%
% Expected Excel columns: Subject, Condition, Frequency, Session,
% Channel, PLV, NumTrials, IsSigHere
%
% Condition mapping:
%   Sit       = Condition=='sit'  & Session=='pre'
%   SitRepeat = Condition=='sit'  & Session=='post'
%   Firm      = Condition=='firm'
%   Foam      = Condition=='foam'
%
% ASSUMPTIONS
%   1) Frequency filter fixed to 333 Hz (freqSelect below); change to 89
%      if you want the low-frequency PLV instead.
%   2) Firm/Foam are matched on Condition only (Session ignored for
%      those two) - only Sit needs Session to disambiguate pre vs post.
%   3) PLV is NOT converted to dB. Unlike SNR (a power ratio), PLV is
%      already a bounded [0,1] measure, so no log transform is applied.
%   4) No vertex ROI restriction - uses ALL electrodes with IsSigHere
%      true across whichever conditions are available for that subject,
%      matching your most recent SNR channel-handling change.
%   5) Same directional test logic as the SNR script: Sit-vs-Firm and
%      Firm-vs-Foam are one-tailed (testing Firm>Sit and Foam>Firm
%      respectively); Sit-vs-SitRepeat is two-sided.
%   6) Plain Bonferroni correction (m = 4, per subject) - matching the
%      current state of your SNR script. Swap in the Holm/FDR blocks
%      from that script if you want a less conservative correction here
%      too.

clc
clear
close all

%% 1. Parameters
excelPath  = 'F:\Auditory-mov-data\PLV_all\FFR_PLV_values.xlsx';
freqSelect = 333;

alpha  = 0.05;

condOrder = {'Sit', 'Firm', 'Foam', 'SitRepeat'};
nCond     = numel(condOrder);

contrastList = { {'Sit', 'SitRepeat'}; ...
                 {'Sit', 'Firm'}; ...
                 {'Sit', 'Foam'}; ...
                 {'Firm','Foam'}};
nContrasts    = numel(contrastList);
contrastNames = cellfun(@(x) sprintf('%s_vs_%s', x{1}, x{2}), contrastList, 'UniformOutput', false);

m_bonf = nContrasts;   % Bonferroni correction across the 4 contrasts, per subject

%% 2. Load table
T = readtable(excelPath);

% Normalise text columns for robust matching regardless of original case
T.Subject   = string(T.Subject);
T.Condition = lower(string(T.Condition));
T.Session   = lower(string(T.Session));
T.Channel   = string(T.Channel);

T = T(T.Frequency == freqSelect, :);

subjects = unique(T.Subject, 'stable');
nSubj    = numel(subjects);

fprintf('%d subjects found in the PLV table at %d Hz.\n', nSubj, freqSelect);

%% 3. Loop over subjects
p_raw_all     = nan(nSubj, nContrasts);
diffMedian_all = nan(nSubj, nContrasts);
nChan_all     = nan(nSubj, nContrasts);
dataMissing   = false(nSubj, nContrasts);
condAvailable = false(nSubj, nCond);

plotPLV = nan(nSubj, nCond);   % per-subject, per-condition mean PLV (for plotting)

for s = 1:nSubj
    subjectid = subjects(s);
    Tsubj = T(T.Subject == subjectid, :);

    labels_cond = struct('Sit', {{}}, 'Firm', {{}}, 'Foam', {{}}, 'SitRepeat', {{}});
    plv_subj    = struct('Sit', [], 'Firm', [], 'Foam', [], 'SitRepeat', []);
    sig_subj    = struct('Sit', {{}}, 'Firm', {{}}, 'Foam', {{}}, 'SitRepeat', {{}});

    for c = 1:nCond
        key = condOrder{c};
        switch key
            case 'Sit'
                rows = Tsubj.Condition == "sit"  & Tsubj.Session == "pre";
            case 'SitRepeat'
                rows = Tsubj.Condition == "sit"  & Tsubj.Session == "post";
            case 'Firm'
                rows = Tsubj.Condition == "firm";
            case 'Foam'
                rows = Tsubj.Condition == "foam";
        end

        Tcond = Tsubj(rows, :);
        if isempty(Tcond)
            continue   % condAvailable(s,c) stays false - e.g. no SitRepeat for this subject
        end

        labels_cond.(key) = cellstr(Tcond.Channel);
        plv_subj.(key)    = Tcond.PLV;
        sig_subj.(key)    = cellstr(Tcond.Channel(logical(Tcond.IsSigHere)));
        condAvailable(s, c) = true;
    end

    availKeys = condOrder(condAvailable(s, :));
    if isempty(availKeys)
        fprintf('Subject %-6s: no data found in table - skipped entirely.\n', subjectid);
        continue
    end

    % Channels significant across whichever conditions are actually
    % available for THIS subject (same logic as the SNR script)
    sig_labels = sig_subj.(availKeys{1});
    for kk = 2:numel(availKeys)
        sig_labels = intersect(sig_labels, sig_subj.(availKeys{kk}), 'stable');
    end
    final_labels = sig_labels;
    nr_chans     = numel(final_labels);


    % Per-condition mean PLV (for the plot), restricted to final_labels
    for c = 1:nCond
        key = condOrder{c};
        if ~condAvailable(s, c)
            continue   % leaves plotPLV(s,c) as NaN
        end

        [isMember, locInCond] = ismember(final_labels, labels_cond.(key));
        if ~all(isMember)
            warning('Subject %s, condition %s: %d of %d selected channels not found - excluded for this condition.', ...
                subjectid, key, sum(~isMember), numel(final_labels));
        end
        idx_this_cond = locInCond(isMember);

        plotPLV(s, c) = mean(plv_subj.(key)(idx_this_cond));
    end

    % Pairwise contrasts - ALL two-sided now
    for k = 1:nContrasts
        keyA = contrastList{k}{1};
        keyB = contrastList{k}{2};
        idxA = find(strcmp(condOrder, keyA));
        idxB = find(strcmp(condOrder, keyB));

        if ~condAvailable(s, idxA) || ~condAvailable(s, idxB)
            dataMissing(s, k) = true;
            continue
        end

        commonChans = intersect(labels_cond.(keyA), labels_cond.(keyB), 'stable');
        commonChans = intersect(commonChans, final_labels, 'stable');

        if numel(commonChans) < 2
            warning('Subject %s, contrast %s: fewer than 2 channels available - skipped.', ...
                subjectid, contrastNames{k});
            continue
        end

        [~, idxInA] = ismember(commonChans, labels_cond.(keyA));
        [~, idxInB] = ismember(commonChans, labels_cond.(keyB));

        plvA = plv_subj.(keyA)(idxInA);
        plvB = plv_subj.(keyB)(idxInB);

        p_raw_all(s, k)      = signrank(plvA, plvB);          % two-sided
        diffMedian_all(s, k) = median(plvA - plvB);            % sign -> direction (A vs B, per contrastList order)
        nChan_all(s, k)      = numel(commonChans);
    end

    fprintf('Subject %-6s processed (conditions available: %s).\n', ...
        subjectid, strjoin(condOrder(condAvailable(s,:)), ', '));
end

fprintf('\n%d / %d subjects have a SIT REPEAT recording.\n', sum(condAvailable(:,4)), nSubj);

p_corr_all = min(p_raw_all * m_bonf, 1);
sig_all    = p_corr_all <= alpha;

%%% Classify each subject/contrast into the SAME 4 categories as the bar plot
% Sign convention SWAPPED (< / <=, matching the SNR script's fix) so
% "Significant effect" means Firm>Sit / Foam>Firm / Foam>Sit, not the
% literal keyA>keyB ordering of contrastList.
category_all = strings(nSubj, nContrasts);
for s = 1:nSubj
    for k = 1:nContrasts
        if isnan(p_raw_all(s, k))
            continue
        end
        d = diffMedian_all(s, k);
        if sig_all(s, k)
            if d < 0
                category_all(s, k) = "Significant effect";
            else
                category_all(s, k) = "Significant opposite effect";
            end
        else
            if d <= 0
                category_all(s, k) = "Non-significant effect";
            else
                category_all(s, k) = "Non-significant opposite effect";
            end
        end
    end
end


%% 4. Save subject-level results
subjectLevelStatsPLV = struct();
subjectLevelStatsPLV.subjects      = cellstr(subjects);
subjectLevelStatsPLV.contrastNames = contrastNames(:)';
subjectLevelStatsPLV.condOrder      = condOrder(:)';
subjectLevelStatsPLV.condAvailable  = condAvailable;
subjectLevelStatsPLV.dataMissing    = dataMissing;
subjectLevelStatsPLV.p_raw          = p_raw_all;
subjectLevelStatsPLV.p_bonf         = p_corr_all;
subjectLevelStatsPLV.significant    = sig_all;
subjectLevelStatsPLV.diffMedian     = diffMedian_all;
subjectLevelStatsPLV.category       = category_all;
subjectLevelStatsPLV.nChannelsUsed  = nChan_all;
subjectLevelStatsPLV.plotPLV        = plotPLV;
subjectLevelStatsPLV.alpha          = alpha;
subjectLevelStatsPLV.bonferroni_m   = m_bonf;
subjectLevelStatsPLV.freqSelect     = freqSelect;

save('PLV333_subjectlevel_pairwise_stats.mat', 'subjectLevelStatsPLV');
fprintf('Saved subject-level PLV stats to PLV333_subjectlevel_pairwise_stats.mat\n');

sigTable = array2table(sig_all, 'VariableNames', contrastNames, 'RowNames', cellstr(subjects));
disp('Per-subject significance (Bonferroni-corrected, m = 4) for the 4 requested contrasts:')
disp(sigTable)
disp('(false may mean "not significant" OR "not tested" - check dataMissing / nChannelsUsed for the reason)')

%% Spaghetti plots - full, Sit-vs-SitRepeat only, and Sit-Firm-Foam only
colOrange      = [230 126  34]/255;   % Significant effect
colDarkPurple  = [ 80  10 120]/255;   % Significant opposite effect
colTan         = [214 197 165]/255;   % Non-significant effect
colLightPurple = [190 160 210]/255;   % Non-significant opposite effect
colDot         = [0.65 0.65 0.65];
colMean        = [0 0 0];
lineAlpha = 0.8;
lineWidth = 0.9;
xMargin   = 0.2;
sText     = 10;
dpi       = 300;

titleBase = sprintf('PLV [%d Hz]', freqSelect);   % FIX: was mislabeled "FFR SNR (dB)"

fig1 = plotSpaghettiSubset(plotPLV, condOrder, contrastList, category_all, ...
    {'SitRepeat','Sit','Firm','Foam'}, colOrange, colDarkPurple, colTan, colLightPurple, colDot, colMean, ...
    lineAlpha, lineWidth, xMargin, 'PLV', ...
    sprintf('%s, subject-level contrasts', titleBase), sText);
set(findobj(fig1, 'Type', 'axes'), 'Color', 'none');

fig2 = plotSpaghettiSubset(plotPLV, condOrder, contrastList, category_all, ...
    {'SitRepeat','Sit'}, colOrange, colDarkPurple, colTan, colLightPurple, colDot, colMean, ...
    lineAlpha, lineWidth, xMargin, 'PLV', ...
    sprintf('%s, Sit vs Sit-repeat only', titleBase), sText);
set(findobj(fig2, 'Type', 'axes'), 'Color', 'none');

fig3 = plotSpaghettiSubset(plotPLV, condOrder, contrastList, category_all, ...
    {'Sit','Firm','Foam'}, colOrange, colDarkPurple, colTan, colLightPurple, colDot, colMean, ...
    lineAlpha, lineWidth, xMargin, 'PLV', ...
    sprintf('%s, Sit-Firm-Foam only', titleBase), sText);
set(findobj(fig3, 'Type', 'axes'), 'Color', 'none');


%% Save the figures
%exportTransparentPNG(fig1, 'PLV_spaghetti_full.png', dpi);
%exportTransparentPNG(fig2, 'PLV_spaghetti_sitrepeat.png', dpi);
exportTransparentPNG(fig3, 'PLV_spaghetti_sitfirmfoam.png', dpi);

%% Local function: spaghetti plot for one subset of conditions
function fig = plotSpaghettiSubset(dataMat, condOrder, contrastList, category_all, ...
    displayOrder, colOrange, colDarkPurple, colTan, colLightPurple, colDot, colMean, ...
    lineAlpha, lineWidth, xMargin, yLabelStr, titlePrefix, sText)

nSubj      = size(dataMat, 1);
nDisp      = numel(displayOrder);
nContrasts = numel(contrastList);

xPos = struct();
for c = 1:nDisp
    xPos.(displayOrder{c}) = c;
end

colIdxDisp = zeros(1, nDisp);
for c = 1:nDisp
    colIdxDisp(c) = find(strcmp(condOrder, displayOrder{c}));
end
dataMat_disp = dataMat(:, colIdxDisp);

validSubj = ~all(isnan(dataMat_disp), 2);

fig = figure('Color', 'w');
hold on; grid on; box on

for c = 1:nDisp
    yVals = dataMat_disp(:, c);
    yVals = yVals(~isnan(yVals));
    scatter(repmat(c, numel(yVals), 1), yVals, 30, colDot, 'filled', ...
        'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', 'none');
end

usedSigEff  = false;
usedSigOpp  = false;
usedNSEff   = false;
usedNSOpp   = false;

for s = 1:nSubj
    if all(isnan(dataMat_disp(s,:))), continue, end

    for k = 1:nContrasts
        keyA = contrastList{k}{1};
        keyB = contrastList{k}{2};

        if ~isfield(xPos, keyA) || ~isfield(xPos, keyB)
            continue
        end

        xA = xPos.(keyA); xB = xPos.(keyB);
        if abs(xA - xB) > 1
            continue
        end

        yA = dataMat(s, strcmp(condOrder, keyA));
        yB = dataMat(s, strcmp(condOrder, keyB));
        if isnan(yA) || isnan(yB)
            continue
        end

        cat = category_all(s, k);
        if cat == ""
            continue
        end

        switch cat
            case "Significant effect"
                lineCol = colOrange;      usedSigEff = true;
            case "Significant opposite effect"
                lineCol = colDarkPurple;  usedSigOpp = true;
            case "Non-significant effect"
                lineCol = colTan;         usedNSEff  = true;
            case "Non-significant opposite effect"
                lineCol = colLightPurple; usedNSOpp  = true;
        end

        plot([xA xB], [yA yB], '-', 'Color', [lineCol lineAlpha], 'LineWidth', lineWidth);
    end
end

meanVals = nan(1, nDisp);
q1Vals   = nan(1, nDisp);
q3Vals   = nan(1, nDisp);
for c = 1:nDisp
    vals = dataMat_disp(:, c);
    vals = vals(~isnan(vals));
    if isempty(vals), continue, end
    meanVals(c) = mean(vals);
    q1Vals(c)   = prctile(vals, 25);
    q3Vals(c)   = prctile(vals, 75);
end

errLow  = meanVals - q1Vals;
errHigh = q3Vals - meanVals;

errorbar(1:nDisp, meanVals, errLow, errHigh, 'ko', 'LineWidth', 2.0, ...
    'MarkerSize', 10, 'MarkerFaceColor', 'k', 'CapSize', 12);

xticks(1:nDisp)
xticklabels(displayOrder)
xlim([1 - xMargin, nDisp + xMargin])

xlabel('Condition')
ylabel(yLabelStr)

set(gca, 'fontsize', sText)
title(sprintf('%s (n = %d)', titlePrefix, sum(validSubj)))

legHandles = [];
legLabels  = {};

if usedSigEff
    legHandles(end+1) = plot(nan, nan, '-', 'Color', colOrange, 'LineWidth', 2);
    legLabels{end+1}  = 'Significant effect';
end
if usedSigOpp
    legHandles(end+1) = plot(nan, nan, '-', 'Color', colDarkPurple, 'LineWidth', 2);
    legLabels{end+1}  = 'Significant opposite effect';
end
if usedNSEff
    legHandles(end+1) = plot(nan, nan, '-', 'Color', colTan, 'LineWidth', 2);
    legLabels{end+1}  = 'Non-significant effect';
end
if usedNSOpp
    legHandles(end+1) = plot(nan, nan, '-', 'Color', colLightPurple, 'LineWidth', 2);
    legLabels{end+1}  = 'Non-significant opposite effect';
end
if any(validSubj)
    legHandles(end+1) = plot(nan, nan, 'o', 'Color', colMean, 'MarkerFaceColor', colMean);
    legLabels{end+1}  = 'Group mean \pm IQR';
end

legend(legHandles, legLabels, 'Location', 'eastoutside', 'Orientation', 'vertical', 'Box', 'off')

hold off
end


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