%% FFR 333 Hz - SUBJECT-LEVEL pairwise Wilcoxon tests + spaghetti plot
%
% For each of the 20 subjects, and for each of 4 requested contrasts:
%   1) SIT       vs SIT REPEAT
%   2) SIT       vs FIRM
%   3) SIT       vs FOAM
%   4) FIRM      vs FOAM
%
% ...this script performs a Wilcoxon signed-rank test using that
% subject's own vertex electrodes (restricted to those significant in
% BOTH conditions of the given contrast) as the paired samples - i.e.
% the same logic as your original single-subject scripts, just looped
% over subjects and over these 4 specific contrasts.
%
% Results (raw p, Bonferroni-corrected p [m = 4], significance flag, and
% number of channels used) are saved per subject to:
%   FFR333_subjectlevel_pairwise_stats.mat
% so they can be reused later without re-running the analysis.
%
% The spaghetti plot at the end then uses each SUBJECT'S OWN saved
% result to color that subject's connecting lines - i.e. two different
% subjects can have differently colored lines for the same contrast,
% depending on whether the effect was significant for them individually.
%
% Noise sidebands: widened to 20 Hz on each side (6 Hz gap kept from the
% signal band edges, as before).
% ROI: vertex electrodes only {Cz,C1,C2,CPz,FC1,FC2,CP1,CP2}.

clc
clear
close all

%% 1. Parameters
subjects = {'NG','BVD','FD','EF','CX','HX','MCBG','LV','AM','FS', ...
            'LVB','RNR','VUR','LH','TF','MGZ','NSMT','MM','TH','TD'};
nSubj = numel(subjects);

alpha = 0.05;
basepath = 'F:\Auditory-mov-data\';

% condition keys -> file-name suffixes
condFileMap = containers.Map({'Sit', 'Firm', 'Foam', 'SitRepeat'}, ...
    {'_FFR_sit_333', '_FFR_firm_333', '_FFR_foam_333', '_FFR_sit_333_2'});

condOrder = {'Sit', 'Firm', 'Foam', 'SitRepeat'};   % x-axis / chronological order
nCond     = numel(condOrder);

% the 4 requested contrasts
contrastList = { {'Sit', 'SitRepeat'}; ...
                 {'Sit', 'Firm'}; ...
                 {'Sit', 'Foam'}; ...
                 {'Firm','Foam'}};
nContrasts = numel(contrastList);
contrastNames = cellfun(@(x) sprintf('%s_vs_%s', x{1}, x{2}), contrastList, 'UniformOutput', false);

m_bonf = nContrasts;   % Bonferroni correction across the 4 requested contrasts, per subject
%m_bonf = 3;

% vertex ROI
vertex_labels = {'Cz','C1','C2','CPz','FC1','FC2','CP1','CP2'};

% FFT parameters
fs   = 10000;
NFFT = ceil(fs/2);
df   = fs/NFFT;
freq_vec = (0:floor(NFFT/2))*df;
nFreqs    = numel(freq_vec);

% Signal band (333 Hz FFR)
sigBandLo = 331; sigBandHi = 335;

% Noise sidebands widened to 20 Hz each, same 6 Hz gap from signal edges
noiseLoLo = sigBandLo - 6 - 20; noiseLoHi = sigBandLo - 6;    % [305 325]
noiseHiLo = sigBandHi + 6; noiseHiHi = sigBandHi + 6 + 20; % [341 361]

sig_idx      = freq_vec >= sigBandLo & freq_vec <= sigBandHi;
noise_lo_idx = freq_vec >= noiseLoLo & freq_vec <= noiseLoHi;
noise_hi_idx = freq_vec >= noiseHiLo & freq_vec <= noiseHiHi;
noise_idx    = noise_lo_idx | noise_hi_idx;


%% 2. Loop over subjects
p_raw_all    = nan(nSubj, nContrasts);
diffMedian_all  = nan(nSubj, nContrasts);
nChan_all    = nan(nSubj, nContrasts);
dataMissing  = false(nSubj, nContrasts);   % true = a required condition wasn't recorded for this subject
condAvailable = false(nSubj, nCond);       % which of the 4 conditions loaded successfully, per subject

plotSNR = nan(nSubj, nCond);   % per-subject, per-condition mean SNR (for plotting)

for s = 1:nSubj
    subjectid = subjects{s};
    inputpath = fullfile(basepath, subjectid, 'preprocessed');

    % Initialise empty containers for this subject; conditions that fail
    % to load simply stay empty, they do not throw an error later on.
    psd_subj    = struct('Sit', [], 'Firm', [], 'Foam', [], 'SitRepeat', []);
    snr_subj    = struct('Sit', [], 'Firm', [], 'Foam', [], 'SitRepeat', []);
    labels_cond = struct('Sit', {{}}, 'Firm', {{}}, 'Foam', {{}}, 'SitRepeat', {{}});
    sig_subj    = struct('Sit', {{}}, 'Firm', {{}}, 'Foam', {{}}, 'SitRepeat', {{}});

    %%% load each condition independently
    for c = 1:nCond
        key    = condOrder{c};
        suffix = condFileMap(key);
        inputFile    = ['averagedFFR_', 'P[', subjectid, ']', suffix, '.mat'];
        inputSigFile = ['SigChansPadded_', 'P[', subjectid, ']', suffix, '.mat'];

        try
            tmp    = load(fullfile(inputpath, inputFile), 'final');
            tmpSig = load(fullfile(inputpath, inputSigFile), 'sigChansBoth');

            avg    = tmp.final.avg;
            labels = tmp.final.label;
            nChan  = size(avg, 1);
            N = size(avg,2);

            win       = hann(N)';
            win_power = sum(win.^2);
            psd_cond  = zeros(nChan, nFreqs);

            for ch = 1:nChan
                signal = avg(ch, 1:2500);

                S = fft(signal.*win, NFFT);
                S = S(1:floor(NFFT/2)+1);
                
                psd_full = (abs(S).^2)/(win_power*fs);
                psd_full(2:end-1) = 2 * psd_full(2:end-1);

                psd_cond(ch, 1:numel(psd_full)) = psd_full;
            end

            psd_subj.(key)    = psd_cond;
            labels_cond.(key) = labels;
            sig_subj.(key)    = tmpSig.sigChansBoth;
            condAvailable(s, c) = true;

        catch
            % File not found (e.g. no SITREPEAT recording for this
            % subject) - leave this condition empty and move on. This
            % does NOT abort processing of the subject's other
            % conditions/contrasts.
            fprintf('Subject %-6s: %-10s data not found - condition excluded for this subject.\n', ...
                subjectid, key);
        end
    end

    % Use ALL electrodes, significant across whichever conditions are
    % actually available for THIS subject - includes SitRepeat when present,
    % falls back to Sit/Firm/Foam only when it's missing (rather than
    % hardcoding just the three core conditions).
    availKeys = condOrder(condAvailable(s, :));

    labels_to_keep = labels_cond.(availKeys{1});
    for kk = 2:numel(availKeys)
        labels_to_keep = intersect(labels_to_keep, labels_cond.(availKeys{kk}), 'stable');
    end

    % Significant across all AVAILABLE conditions for this subject
    sig_labels = sig_subj.(availKeys{1});
    for kk = 2:numel(availKeys)
        sig_labels = intersect(sig_labels, sig_subj.(availKeys{kk}), 'stable');
    end

    final_labels = intersect(labels_to_keep, sig_labels, 'stable');
    nr_chans     = numel(final_labels);
    
    % fprintf('Subject %-6s: channels significant across %s: %d\n', ...
    %     subjectid, strjoin(availKeys, '/'), nr_chans);

    for c = 1:nCond
        key = condOrder{c};
        if ~condAvailable(s, c)
            continue   % leaves plotSNR(s,c) as NaN
        end

        % FIX: index into THIS condition's own label list (labels_cond.(key)),
        % which is what psd_subj.(key)'s rows actually correspond to - not
        % into labels_to_keep/final_labels, whose row order can differ.
        [isMember, locInCond] = ismember(final_labels, labels_cond.(key));
        if ~all(isMember)
            warning('Subject %s, condition %s: %d of %d selected channels not found in this montage - excluded for this condition.', ...
                subjectid, key, sum(~isMember), numel(final_labels));
        end
        
        idx_this_cond = locInCond(isMember);

        psd_subj.(key) = psd_subj.(key)(idx_this_cond, :);

        % SNR = mean signal power / mean noise power (per channel)
        snr_subj.(key) = mean(psd_subj.(key)(:, sig_idx), 2) ./ mean(psd_subj.(key)(:, noise_idx), 2);

        % Convert to decibels (per channel)
        snr_subj.(key) = 10*log10(snr_subj.(key));

        % Overall average SNR per condition
        plotSNR(s, c) = mean(snr_subj.(key));
    end
    

    % Wilcoxon signed-rank test - ALL 4 contrasts now two-sided
    for k = 1:nContrasts
        keyA = contrastList{k}{1};
        keyB = contrastList{k}{2};
        idxA = find(strcmp(condOrder, keyA));
        idxB = find(strcmp(condOrder, keyB));

        if ~condAvailable(s, idxA) || ~condAvailable(s, idxB)
            dataMissing(s, k) = true;
            continue
        end

        snrA = snr_subj.(keyA);
        snrB = snr_subj.(keyB);

        p_raw_all(s, k)      = signrank(snrA, snrB);        % two-sided
        diffMedian_all(s, k) = median(snrA - snrB);          % sign -> direction (A vs B, per contrastList order)
        nChan_all(s, k)      = numel(snrA);
    end
end

fprintf('\n%d / %d subjects have a SITREPEAT recording.\n', sum(condAvailable(:,4)), nSubj);

p_corr_all = min(p_raw_all*m_bonf, 1);
sig_all    = p_corr_all <= alpha;

%%% classify each subject/contrast into the SAME 4 categories as the bar plot
category_all = strings(nSubj, nContrasts);   % "" = not tested
for s = 1:nSubj
    for k = 1:nContrasts
        if isnan(p_raw_all(s, k))
            continue   % stays ""
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


% %% Holm-Bonferroni correction (per subject, across that subject's actually-tested contrasts)
% % More powerful than standard Bonferroni while still controlling the
% % family-wise error rate: sort each subject's raw p-values ascending,
% % then step down through them, comparing p_(i) to alpha/(m-i+1). NaN
% % (untested) contrasts are excluded from that subject's correction
% % family rather than counted as an extra test.
% p_corr_all = nan(nSubj, nContrasts);
% 
% for s = 1:nSubj
%     validK = find(~isnan(p_raw_all(s, :)));
%     nValid = numel(validK);
%     if nValid == 0
%         continue
%     end
% 
%     [sortedP, sortOrder] = sort(p_raw_all(s, validK));
% 
%     adjP = nan(1, nValid);
%     runningMax = 0;
%     for i = 1:nValid
%         stepP      = (nValid - i + 1) * sortedP(i);
%         runningMax = max(runningMax, stepP);   % enforces monotonicity
%         adjP(i)    = min(runningMax, 1);
%     end
% 
%     % map the Holm-adjusted p-values back to their original contrast columns
%     p_corr_all(s, validK(sortOrder)) = adjP;
% end
% 
% sig_all = p_corr_all <= alpha;   % NaN comparisons evaluate to false automatically
% 
% 
% %% Benjamini-Hochberg (FDR) correction, per subject, across that subject's
% % actually-tested contrasts (NaN/untested contrasts excluded from the
% %family rather than counted as an extra test).
% p_corr_all = nan(nSubj, nContrasts);
% 
% for s = 1:nSubj
%     validK = find(~isnan(p_raw_all(s, :)));
%     nValid = numel(validK);
%     if nValid == 0
%         continue
%     end
% 
%     [sortedP, sortOrder] = sort(p_raw_all(s, validK));
% 
%     % step-up BH-adjusted p-values (q-values): start from the largest
%     % p-value and work down, enforcing monotonicity
%     adjP = nan(1, nValid);
%     adjP(nValid) = sortedP(nValid);
%     for i = (nValid - 1):-1:1
%         adjP(i) = min(adjP(i + 1), (nValid / i) * sortedP(i));
%     end
%     adjP = min(adjP, 1);
% 
%     % map the BH-adjusted p-values back to their original contrast columns
%     p_corr_all(s, validK(sortOrder)) = adjP;
% end
% 
% sig_all = p_corr_all <= alpha;   % NaN comparisons evaluate to false automatically

%% 3. Assemble & save subject-level results
subjectLevelStats = struct();
subjectLevelStats.subjects       = subjects(:);
subjectLevelStats.contrastNames  = contrastNames(:)';
subjectLevelStats.condOrder      = condOrder(:)';
subjectLevelStats.condAvailable  = condAvailable;   % subjects x conditions, logical
subjectLevelStats.dataMissing    = dataMissing;     % subjects x contrasts, logical (required condition absent)
subjectLevelStats.p_raw          = p_raw_all;
subjectLevelStats.p_bonf         = p_corr_all;
subjectLevelStats.significant    = sig_all;
subjectLevelStats.diffMedian     = diffMedian_all;
subjectLevelStats.category       = category_all;
subjectLevelStats.nChannelsUsed  = nChan_all;
subjectLevelStats.plotSNR        = plotSNR;
subjectLevelStats.alpha          = alpha;
subjectLevelStats.bonferroni_m   = m_bonf;

save('FFR333_subjectlevel_pairwise_stats.mat', 'subjectLevelStats');
fprintf('Saved subject-level pairwise stats to FFR333_subjectlevel_pairwise_stats.mat\n');

% Human-readable summary table
sigTable = array2table(sig_all, 'VariableNames', contrastNames, 'RowNames', subjects);
disp('Per-subject significance (Bonferroni-corrected, m = 4) for the 4 requested contrasts:')
disp(sigTable)
disp('(false may mean "not significant" OR "not tested" - check dataMissing / nChannelsUsed for the reason)')

 
%% 4. Spaghetti plots - full, Sit-vs-SitRepeat only, and Sit-Firm-Foam only
% Colors now match the bar-plot script's 4-category scheme exactly.
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

titleBase = sprintf('FFR SNR (dB) [%d-%d Hz]', sigBandLo, sigBandHi); % adjust to however your title currently references frequency

% 1) Full plot
fig1 = plotSpaghettiSubset(plotSNR, condOrder, contrastList, category_all, ...
    {'SitRepeat','Sit','Firm','Foam'}, colOrange, colDarkPurple, colTan, colLightPurple, colDot, colMean, ...
    lineAlpha, lineWidth, xMargin, 'FFR SNR (dB)', ...
    sprintf('%s, subject-level contrasts', titleBase), sText);
set(findobj(fig1, 'Type', 'axes'), 'Color', 'none');

% 2) Sit vs SitRepeat only
fig2 = plotSpaghettiSubset(plotSNR, condOrder, contrastList, category_all, ...
    {'SitRepeat','Sit'}, colOrange, colDarkPurple, colTan, colLightPurple, colDot, colMean, ...
    lineAlpha, lineWidth, xMargin, 'FFR SNR (dB)', ...
    sprintf('%s, Sit vs Sit-repeat only', titleBase), sText);
set(findobj(fig2, 'Type', 'axes'), 'Color', 'none');

% 3) Sit-Firm-Foam only (no SitRepeat)
fig3 = plotSpaghettiSubset(plotSNR, condOrder, contrastList, category_all, ...
    {'Sit','Firm','Foam'}, colOrange, colDarkPurple, colTan, colLightPurple, colDot, colMean, ...
    lineAlpha, lineWidth, xMargin, 'FFR SNR (dB)', ...
    sprintf('%s, Sit-Firm-Foam only', titleBase), sText);
set(findobj(fig3, 'Type', 'axes'), 'Color', 'none');

%% Save figures
exportTransparentPNG(fig1, 'SNR_333_spaghetti_full.png', dpi);
exportTransparentPNG(fig2, 'SNR_333_spaghetti_sitrepeat.png', dpi);
exportTransparentPNG(fig3, 'SNR_333_spaghetti_sitfirmfoam.png', dpi);


%% Local function: spaghetti plot for one subset of conditions
function fig = plotSpaghettiSubset(dataMat, condOrder, contrastList, category_all, ...
    displayOrder, colOrange, colDarkPurple, colTan, colLightPurple, colDot, colMean, ...
    lineAlpha, lineWidth, xMargin, yLabelStr, titlePrefix, sText)
% Draws one spaghetti plot restricted to the conditions in displayOrder.
% Line color now encodes the subject's CATEGORY for that contrast
% (matching the bar-plot script's 4-category scheme), not just a
% significant/not-significant binary.

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
            continue   % not tested
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