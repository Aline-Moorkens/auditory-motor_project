%% SNR directionality summary bar plot (foam>firm, foam>sit, firm>sit)
%
% For each subject, and for each of the 3 pairwise comparisons among
% Sit/Firm/Foam, this script:
%   1) Computes per-channel SNR (dB), restricted to channels that are
%      BOTH shared across Sit/Firm/Foam AND significant across all three,
%      for that subject. SitRepeat is NOT required/used here at all -
%      this is a standalone script independent of the 4-contrast
%      subject-level pipeline.
%   2) Runs a TWO-SIDED Wilcoxon signed-rank test (channels as paired
%      samples) and determines the SIGN of the median per-channel
%      difference, classifying the subject's result for that contrast
%      into one of 4 categories:
%        - "Significant effect"            : p<alpha AND expected direction (A>B)
%        - "Significant opposite effect"   : p<alpha AND reversed direction (A<B)
%        - "Non-significant effect"        : p>=alpha AND expected direction
%        - "Non-significant opposite effect": p>=alpha AND reversed direction
%   3) A stacked bar chart (one bar per contrast) shows how many subjects
%      fall into each category.
%
% ASSUMPTIONS:
%   - alpha = 0.05, UNCORRECTED across the 3 contrasts - treated as a
%     descriptive per-subject summary, not a family of confirmatory
%     tests.
%
%   - A subject is included only if Sit, Firm, AND Foam all loaded
%     successfully (SitRepeat availability is irrelevant to this chart).
%
%   - SNR is converted to dB before testing.

clc
clear
close all

%% 1. Parameters
subjects = {'NG','BVD','FD','EF','CX','HX','MCBG','LV','AM','FS', ...
            'LVB','RNR','VUR','LH','TF','MGZ','NSMT','MM','TH','TD'};
nSubj = numel(subjects);

alpha    = 0.05;   % uncorrected - see note above
basepath = 'F:\Auditory-mov-data\';

condFileMap = containers.Map({'Sit','Firm','Foam'}, ...
    {'_FFR_sit_333','_FFR_firm_333','_FFR_foam_333'});
condOrder = {'Sit','Firm','Foam'};
nCond     = numel(condOrder);

% The 3 directional comparisons to summarize (A,B => testing A>B)
comparisons = { {'Foam','Firm'}; {'Foam','Sit'}; {'Firm','Sit'} };
compLabels  = {'foam > firm', 'foam > sit', 'firm > sit'};
nComp       = numel(comparisons);

% FFT / SNR parameters, matching the main SNR pipeline
fs   = 10000;
NFFT = ceil(fs/2);
df   = fs/NFFT;
freq_vec = (0:floor(NFFT/2))*df;

sigBandLo = 331; sigBandHi = 335;
noiseLoLo = sigBandLo - 6 - 20;  noiseLoHi = sigBandLo - 6;
noiseHiLo = sigBandHi + 6;       noiseHiHi = sigBandHi + 6 + 20;

sig_idx      = freq_vec >= sigBandLo & freq_vec <= sigBandHi;
noise_lo_idx = freq_vec >= noiseLoLo & freq_vec <= noiseLoHi;
noise_hi_idx = freq_vec >= noiseHiLo & freq_vec <= noiseHiHi;
noise_idx    = noise_lo_idx | noise_hi_idx;

%% 2. Loop over subjects
category  = strings(nSubj, nComp);   % category label per subject/comparison
validSubj = false(nSubj, 1);

for s = 1:nSubj
    subjectid = subjects{s};
    inputpath = fullfile(basepath, subjectid, 'preprocessed');

    psd_cond    = struct('Sit', [], 'Firm', [], 'Foam', []);
    labels_cond = struct('Sit', {{}}, 'Firm', {{}}, 'Foam', {{}});
    sig_cond    = struct('Sit', {{}}, 'Firm', {{}}, 'Foam', {{}});
    condAvail   = false(1, nCond);

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

            win       = hann(2500)';
            win_power = sum(win.^2);
            psd       = zeros(nChan, floor(NFFT/2)+1);

            for ch = 1:nChan
                signal = avg(ch, 1:2500);
                S = fft(signal .* win, NFFT);
                S = S(1:floor(NFFT/2)+1);
                psd_full = (abs(S).^2)/(win_power*fs);
                psd_full(2:end-1) = 2*psd_full(2:end-1);
                psd(ch,:) = psd_full;
            end

            psd_cond.(key)    = psd;
            labels_cond.(key) = labels;
            sig_cond.(key)    = tmpSig.sigChansBoth;
            condAvail(c)      = true;
        catch
            fprintf('Subject %-6s: %-6s data not found - subject excluded (needs Sit, Firm, AND Foam).\n', ...
                subjectid, key);
        end
    end

    if ~all(condAvail)
        continue   % validSubj(s) stays false
    end
    validSubj(s) = true;

    % Channels significant across all 3 (Sit, Firm, Foam) for this subject
    labels_to_keep = intersect(labels_cond.Sit, labels_cond.Firm, 'stable');
    labels_to_keep = intersect(labels_to_keep, labels_cond.Foam, 'stable');
    labels_to_keep = labels_to_keep(~ismember(labels_to_keep, {'TP9','TP10'}));

    sig_labels = intersect(sig_cond.Sit, sig_cond.Firm, 'stable');
    sig_labels = intersect(sig_labels, sig_cond.Foam, 'stable');

    final_labels = intersect(labels_to_keep, sig_labels, 'stable');
    if numel(final_labels) < 2
        warning('Subject %s: fewer than 2 qualifying channels across Sit/Firm/Foam - excluded.', subjectid);
        validSubj(s) = false;
        continue
    end

    % Per-condition, per-channel SNR (dB), restricted to final_labels
    snr_cond = struct();
    for c = 1:nCond
        key = condOrder{c};
        [~, locInCond] = ismember(final_labels, labels_cond.(key));
        psd_sel = psd_cond.(key)(locInCond, :);
        snr_lin = mean(psd_sel(:, sig_idx), 2) ./ mean(psd_sel(:, noise_idx), 2);
        snr_cond.(key) = 10 * log10(snr_lin);   % dB, matching main pipeline convention
    end

    for kk = 1:nComp
        keyA = comparisons{kk}{1};
        keyB = comparisons{kk}{2};

        valA = snr_cond.(keyA);
        valB = snr_cond.(keyB);

        p = signrank(valA, valB);        % TWO-SIDED TEST
        d = median(valA - valB);         % sign of the observed effect

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
fprintf('\n%d / %d subjects included (require Sit, Firm, AND Foam; SitRepeat not needed).\n', nValid, nSubj);

%% 3. Tally counts per category, per comparison
catNames = ["Significant effect", "Significant opposite effect", ...
            "Non-significant effect", "Non-significant opposite effect"];
counts = zeros(nComp, numel(catNames));

for kk = 1:nComp
    for ci = 1:numel(catNames)
        counts(kk, ci) = sum(category(validSubj, kk) == catNames(ci));
    end
end

%% 4. Stacked bar plot, styled to match the reference figure
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

% count labels centered within each non-empty segment
for kk = 1:nComp
    yBase = 0;
    for ci = 1:numel(catNames)
        v = counts(kk, ci);
        if v > 0
            yMid = yBase + v/2;
            txtColor = 'k';
            if ci == 1   % orange segment gets white text, matching reference
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

exportTransparentPNG(gcf, 'SNR_directionality_barplot.png', 300);

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