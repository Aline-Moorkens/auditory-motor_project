%% PLV directionality summary bar plot (foam>firm, foam>sit, firm>sit)
%
% Mirrors SNR_directionality_barplot.m, but reads pre-computed PLV values
% from the Excel table instead of raw EEG, and does NOT apply a dB
% transform (PLV is already a bounded [0,1] measure).
%
% ASSUMPTIONS:
%   - alpha = 0.05, UNCORRECTED across the 3 contrasts (no Bonferroni/Holm/FDR)
%   --> reason: treating this as a descriptive per-subject summary rather than 
% a confirmatory family of tests, matching how this kind of figure is typically reported.
%
%   - Channels restricted to those with IsSigHere true AND shared across
%     Sit/Firm/Foam only - SitRepeat/Session is ignored entirely for this
%     chart (only the 'pre' Sit rows are used, to avoid double-counting
%     the repeat recording as if it were a separate "Sit").
%
%   - A subject is included only if Sit, Firm, AND Foam all have data.

clc
clear
close all

%% 1. Parameters
excelPath  = 'F:\Auditory-mov-data\PLV_all\FFR_PLV_values.xlsx';
freqSelect = 89;

alpha = 0.05;

comparisons = { {'Foam','Firm'}; {'Foam','Sit'}; {'Firm','Sit'} };
compLabels  = {'foam > firm', 'foam > sit', 'firm > sit'};
nComp       = numel(comparisons);

condOrder = {'Sit','Firm','Foam'};
nCond     = numel(condOrder);

%% 2. Load table
T = readtable(excelPath);
T.Subject   = string(T.Subject);
T.Condition = lower(string(T.Condition));
T.Channel   = string(T.Channel);
T = T(T.Frequency == freqSelect, :);

% Only the 'pre' Sit rows count as "Sit" here - excludes the repeat
% recording entirely, so it's never confused with the main Sit condition
if ismember('Session', T.Properties.VariableNames)
    T.Session   = lower(string(T.Session));
    isSitPre    = T.Condition == "sit" & T.Session == "pre";
    isOtherCond = T.Condition ~= "sit";
    T = T(isSitPre | isOtherCond, :);
end

subjects = unique(T.Subject, 'stable');
nSubj    = numel(subjects);

fprintf('%d subjects found in the PLV table at %d Hz.\n', nSubj, freqSelect);

%% 3. Loop over subjects
category  = strings(nSubj, nComp);
validSubj = false(nSubj, 1);

for s = 1:nSubj
    subjectid = subjects(s);
    Tsubj = T(T.Subject == subjectid, :);

    labels_cond = struct('Sit', {{}}, 'Firm', {{}}, 'Foam', {{}});
    plv_subj    = struct('Sit', [], 'Firm', [], 'Foam', []);
    sig_subj    = struct('Sit', {{}}, 'Firm', {{}}, 'Foam', {{}});
    condAvail   = false(1, nCond);

    for c = 1:nCond
        key  = condOrder{c};
        rows = Tsubj.Condition == lower(string(key));
        Tcond = Tsubj(rows, :);
        if isempty(Tcond)
            continue
        end

        labels_cond.(key) = cellstr(Tcond.Channel);
        plv_subj.(key)    = Tcond.PLV;
        sig_subj.(key)    = cellstr(Tcond.Channel(logical(Tcond.IsSigHere)));
        condAvail(c)      = true;
    end

    if ~all(condAvail)
        continue   % validSubj(s) stays false
    end
    validSubj(s) = true;

    labels_to_keep = intersect(labels_cond.Sit, labels_cond.Firm, 'stable');
    labels_to_keep = intersect(labels_to_keep, labels_cond.Foam, 'stable');

    sig_labels = intersect(sig_subj.Sit, sig_subj.Firm, 'stable');
    sig_labels = intersect(sig_labels, sig_subj.Foam, 'stable');

    final_labels = intersect(labels_to_keep, sig_labels, 'stable');
    if numel(final_labels) < 2
        warning('Subject %s: fewer than 2 qualifying channels across Sit/Firm/Foam - excluded.', subjectid);
        validSubj(s) = false;
        continue
    end

    plv_cond = struct();
    for c = 1:nCond
        key = condOrder{c};
        [~, locInCond] = ismember(final_labels, labels_cond.(key));
        plv_cond.(key) = plv_subj.(key)(locInCond);
    end

    for kk = 1:nComp
        keyA = comparisons{kk}{1};
        keyB = comparisons{kk}{2};

        valA = plv_cond.(keyA);
        valB = plv_cond.(keyB);

        p = signrank(valA, valB);   % TWO-SIDED
        d = median(valA - valB);

        if p < alpha
            if d > 0
                category(s, kk) = "Significant effect";
            else
                category(s, kk) = "Significant opposite effect";
            end
        else
            if d >= 0
                category(s, kk) = "Non-significant effect";
            else
                category(s, kk) = "Non-significant opposite effect";
            end
        end
    end

    fprintf('Subject %-6s processed.\n', subjectid);
end

nValid = sum(validSubj);
fprintf('\n%d / %d subjects included (require Sit, Firm, AND Foam).\n', nValid, nSubj);

%% 4. Tally counts per category, per comparison
catNames = ["Significant effect", "Significant opposite effect", ...
            "Non-significant effect", "Non-significant opposite effect"];
counts = zeros(nComp, numel(catNames));

for kk = 1:nComp
    for ci = 1:numel(catNames)
        counts(kk, ci) = sum(category(validSubj, kk) == catNames(ci));
    end
end

%% 5. Stacked bar plot, styled to match the reference figure
colOrange      = [230 126  34]/255;   % Significant effect
colDarkPurple  = [ 80  10 120]/255;   % Significant opposite effect
colTan         = [214 197 165]/255;   % Non-significant effect
colLightPurple = [190 160 210]/255;   % Non-significant opposite effect
barColors = [colOrange; colDarkPurple; colTan; colLightPurple];

figure('Color', 'w');
b = bar(counts, 'stacked','EdgeColor','none');
for ci = 1:numel(catNames)
    b(ci).FaceColor = barColors(ci, :);
end

hold on
yline(nValid, ':k', sprintf('N = %d', nValid), ...
    'LabelHorizontalAlignment', 'right', 'FontAngle', 'italic');

for kk = 1:nComp
    yBase = 0;
    for ci = 1:numel(catNames)
        v = counts(kk, ci);
        if v > 0
            yMid = yBase + v/2;
            txtColor = 'k';
            if ci == 1
                txtColor = 'w';
            end
            text(kk, yMid, sprintf('%d', v), 'HorizontalAlignment', 'center', ...
                'FontWeight', 'bold', 'Color', txtColor);
        end
        yBase = yBase + v;
    end
end

xticks(1:nComp)
xticklabels(compLabels)
ylabel('Number of participants')
ylim([0, nValid + 1])

legend(catNames, 'Location', 'eastoutside', 'Box', 'off','Color','none','Direction','normal')
box off
hold off


%% save the figures
% set(gca, 'Color', 'none');
% exportgraphics(gcf, 'SNR_directionality_barplot.png', ...
%     'BackgroundColor', 'none', 'Resolution', 300);

exportTransparentPNG(gcf, 'PLV_89_directionality_barplot.png', 300);

function exportTransparentPNG(fig, filename, dpi)
% Exports a figure as PNG with a transparent background AND correct DPI
% metadata (exportgraphics/imwrite's white-background workaround from
% before doesn't carry DPI info through the imread/imwrite round-trip,
% so viewers were defaulting to 96 DPI even though the pixel data itself
% was already rendered at the requested resolution).

    tmpFile = [tempname, '.png'];
    exportgraphics(fig, tmpFile, 'Resolution', dpi);   % white background - unavoidable here

    img = imread(tmpFile);
    isWhite = all(img == 255, 3);
    alphaMask = double(~isWhite);

    % DPI -> pixels-per-meter for the PNG's resolution metadata
    pxPerMeter = round(dpi / 0.0254);

    imwrite(img, filename, ...
        'Alpha', alphaMask, ...
        'ResolutionUnit', 'meter', ...
        'XResolution', pxPerMeter, ...
        'YResolution', pxPerMeter);

    delete(tmpFile);
end