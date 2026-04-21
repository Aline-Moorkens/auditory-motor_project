%% Compare standing on foam vs standing no foam vs sitting
clc
clear
close all

%% Set parameters
subjectid = 'NG';

% input
inputpath = fullfile('F:\Auditory-mov-data',subjectid,'preprocessed');

nameSit  = '_FFR_sit_333';
nameFirm = '_FFR_firm_333';
nameFoam = '_FFR_foam_333';

inputSit  = ['averagedFFR_','P[',subjectid,']',nameSit, '.mat'];
inputFirm = ['averagedFFR_','P[',subjectid,']',nameFirm,'.mat'];
inputFoam = ['averagedFFR_','P[',subjectid,']',nameFoam,'.mat'];

%% Load data
load(fullfile(inputpath, inputSit),  'final'); ffr_sit  = final;
load(fullfile(inputpath, inputFirm), 'final'); ffr_firm = final;
load(fullfile(inputpath, inputFoam), 'final'); ffr_foam = final;

%% Frequency Analysis (Power Spectrum)
% FFT parameters
fs   = 10000;
N    = size(ffr_sit.avg, 2);
NFFT = ceil(fs / 2);

df        = fs / NFFT;
freq_vec  = (0 : floor(NFFT/2)) * df;
nFreqs    = numel(freq_vec);

conditions = {'sit', 'firm', 'foam'};
data_all   = {ffr_sit, ffr_firm, ffr_foam};

psd_all    = cell(1, 3);
labels_all = cell(1, 3);

for c = 1:3
    avg    = data_all{c}.avg;
    labels = data_all{c}.label;
    nChan  = size(avg, 1);
    N      = size(avg, 2);
    win       = hann(N)';
    win_power = sum(win .^ 2);
    psd_cond  = zeros(nChan, nFreqs);

    for ch = 1:nChan
        signal = avg(ch, :);

        S = fft(signal .* win, NFFT);
        S = S(1 : floor(NFFT/2) + 1);

        psd_full = (abs(S) .^ 2) / (win_power * fs);
        psd_full(2:end-1) = 2 * psd_full(2:end-1);

        psd_cond(ch, 1:numel(psd_full)) = psd_full;
    end

    psd_all{c}    = psd_cond;
    labels_all{c} = labels;

    fprintf('PSD computed: %s | %d channels | %d samples\n', ...
        conditions{c}, nChan, N);
end

psd_sit   = psd_all{1};   labels_sit  = labels_all{1};
psd_firm  = psd_all{2};   labels_firm = labels_all{2};
psd_foam  = psd_all{3};   labels_foam = labels_all{3};

fprintf('Frequency resolution: %.2f Hz\n', df);
fprintf('Frequency bins: %d\n', nFreqs);

%% Channel handling
chan_labels = intersect(labels_sit,  labels_firm, 'stable');
chan_labels = intersect(chan_labels, labels_foam,  'stable');

fprintf('Channels common to all conditions: %d\n', numel(chan_labels));

chans_to_remove = {'TP9', 'TP10'};
idx        = ismember(chan_labels, chans_to_remove);
chan_labels = chan_labels(~idx);
nr_chans   = numel(chan_labels);

fprintf('Channels after exclusion: %d\n', nr_chans);

idx_sit  = find(ismember(labels_sit,  chan_labels));
idx_firm = find(ismember(labels_firm, chan_labels));
idx_foam = find(ismember(labels_foam, chan_labels));

psd_sit  = psd_sit(idx_sit,   :);
psd_firm = psd_firm(idx_firm,  :);
psd_foam = psd_foam(idx_foam,  :);

%% Compute SNR per channel
% SNR = power in signal band / mean power in noise bands
% Signal band: 331-335 Hz (around 333 Hz FFR)
% Noise bands: equal-width flanking bands on either side, avoiding the signal
% Here: [321-325 Hz] and [341-345 Hz] (5 Hz bands, 6 Hz gap from signal)

sig_idx  = freq_vec >= 331 & freq_vec <= 335;

% Flanking noise bands — same width as signal band
noise_lo_idx = freq_vec >= 321 & freq_vec <= 325;
noise_hi_idx = freq_vec >= 341 & freq_vec <= 345;
noise_idx    = noise_lo_idx | noise_hi_idx;

% SNR = mean signal power / mean noise power (per channel)
snr_sit  = mean(psd_sit( :, sig_idx), 2) ./ mean(psd_sit( :, noise_idx), 2);
snr_firm = mean(psd_firm(:, sig_idx), 2) ./ mean(psd_firm(:, noise_idx), 2);
snr_foam = mean(psd_foam(:, sig_idx), 2) ./ mean(psd_foam(:, noise_idx), 2);

fprintf('SNR computed for %d channels across 3 conditions\n', nr_chans);

%% Non-parametric test for 3 related conditions: Friedman test
data = [snr_foam, snr_firm, snr_sit];

[p_friedman, tbl_friedman, stats_friedman] = friedman(data, 1, 'off');

disp('Friedman test p-value:')
disp(p_friedman)

%% Post-hoc pairwise comparisons (Wilcoxon signed-rank) + Bonferroni
[p_Foam_vs_NoFoam, h_Foam_vs_NoFoam, stats_Foam_vs_NoFoam] = ...
    signrank(snr_foam, snr_firm, 'tail', 'right');

[p_Foam_vs_Sit, h_Foam_vs_Sit, stats_Foam_vs_Sit] = ...
    signrank(snr_foam, snr_sit, 'tail', 'right');

[p_NoFoam_vs_Sit, h_NoFoam_vs_Sit, stats_NoFoam_vs_Sit] = ...
    signrank(snr_firm, snr_sit, 'tail', 'right');

m         = 3;
alpha     = 0.05;
alphaBonf = alpha / m;

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
nr_conditions = 3;
bp = zeros([nr_chans, nr_conditions]);
bp(:,1) = snr_foam;
bp(:,2) = snr_firm;
bp(:,3) = snr_sit;

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

title('FFR SNR comparison [331–335 Hz]')
ylabel('FFR SNR (signal power / noise power)')

hold on
parallelcoords(bp, 'Color', 0.7*[1 1 1], ...
                   'Linestyle', '-', 'Marker', '.', 'Markersize', 7);

%% Significance annotation (Bonferroni-corrected)
yMax   = max(bp(:));
yMin   = min(bp(:));
yRange = yMax - yMin;
if yRange == 0
    yRange = 1;
end
yStep  = 0.02 * yRange;
startY = yMax + 0.01 * yRange;

ylim([0, yMax + 4*yStep]);

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