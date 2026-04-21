%% Compare standing on foam vs standing no foam vs sitting
clc
clear
close all

%% Set parameters
subjectid = 'NG';

% input
inputpath = fullfile('F:\Auditory-mov-data',subjectid,'preprocessed');

nameSit = '_FFR_sit_333';
nameFirm = '_FFR_firm_333';
nameFoam = '_FFR_foam_333';

inputSit = ['averagedFFR_','P[',subjectid,']',nameSit,'.mat'];
inputFirm = ['averagedFFR_','P[',subjectid,']',nameFirm,'.mat'];
inputFoam = ['averagedFFR_','P[',subjectid,']',nameFoam,'.mat'];

%% Load data
load(fullfile(inputpath, inputSit),  'final'); ffr_sit  = final;
load(fullfile(inputpath, inputFirm), 'final'); ffr_firm = final;
load(fullfile(inputpath, inputFoam), 'final'); ffr_foam = final;

%% Frequency Analysis (Power Spectrum)
% FFT parameters
fs   = 10000;                       % sampling rate (Hz)
N    = size(ffr_sit.avg, 2);        % number of samples (from data)
NFFT = ceil(fs / 2);                % FFT length -> df ≈ 2 Hz

df        = fs / NFFT;                        % frequency resolution (Hz)
win       = hann(N)';                         % Hann window [1 x N]
win_power = sum(win .^ 2);                    % window energy correction

% One-sided frequency axis
freq_vec = (0 : floor(NFFT/2)) * df;          % [1 x (NFFT/2 + 1)] in Hz
nFreqs   = numel(freq_vec);

conditions = {'sit', 'firm', 'foam'};
data_all   = {ffr_sit, ffr_firm, ffr_foam};

% Preallocate outputs — cell arrays because channel count may differ
psd_all    = cell(1, 3);
labels_all = cell(1, 3);

for c = 1:3
     avg    = data_all{c}.avg;       % [nChan x nSamples]
     labels = data_all{c}.label;     % {nChan x 1} channel labels
     nChan  = size(avg, 1);
     psd_cond = zeros(nChan, nFreqs);

    for ch = 1:nChan
        signal = avg(ch, :);        % [1 x N]

        % Step 1: Apply Hann window and compute FFT
        S = fft(signal .* win, NFFT);

        % Step 2: Keep one-sided spectrum
        S = S(1 : floor(NFFT/2) + 1);

        % Step 3: Convert to PSD (µV²/Hz)
        psd_full = (abs(S) .^ 2) / (win_power * fs);

        % Step 4: Double all bins except DC and Nyquist
        psd_full(2:end-1) = 2 * psd_full(2:end-1);

        % Step 5: Store
        psd_cond(ch, 1:numel(psd_full)) = psd_full;
    end

    psd_all{c}    = psd_cond;
    labels_all{c} = labels;         % store channel labels for this condition

    fprintf('PSD computed: %s | %d channels | %d samples\n', ...
        conditions{c}, nChan, N);
end

% Unpack back to named variables
% Unpack to named variables
psd_sit   = psd_all{1};   labels_sit  = labels_all{1};
psd_firm  = psd_all{2};   labels_firm = labels_all{2};
psd_foam  = psd_all{3};   labels_foam = labels_all{3};

fprintf('Frequency resolution: %.2f Hz\n', df);
fprintf('Frequency bins: %d\n', nFreqs);


%% Channel handling
% Find channels present in all three conditions
chan_labels = intersect(labels_sit,  labels_firm, 'stable');
chan_labels = intersect(chan_labels, labels_foam,  'stable');

fprintf('Channels common to all conditions: %d\n', numel(chan_labels));

% Remove unwanted channels
chans_to_remove = {'TP9', 'TP10'};     % mastoids used for reference
idx        = ismember(chan_labels, chans_to_remove);
chan_labels = chan_labels(~idx);
nr_chans   = numel(chan_labels);

fprintf('Channels after exclusion: %d\n', nr_chans);

% Find indices of common channels in each condition's label list
idx_sit  = find(ismember(labels_sit,  chan_labels));
idx_firm = find(ismember(labels_firm, chan_labels));
idx_foam = find(ismember(labels_foam, chan_labels));

% Apply channel selection to PSD matrices
psd_sit  = psd_sit(idx_sit,   :);
psd_firm = psd_firm(idx_firm,  :);
psd_foam = psd_foam(idx_foam,  :);

%% Average power across frequency bins within 331-335 Hz
freq_idx = freq_vec >= 331 & freq_vec <= 335;

pow_sit  = mean(psd_sit( :, freq_idx), 2);   % [nChan x 1]
pow_firm = mean(psd_firm(:, freq_idx), 2);
pow_foam = mean(psd_foam(:, freq_idx), 2);


%% Non-parametric test for 3 related conditions: Friedman test

% Data matrix: rows = channels (subjects), columns = conditions
% Order: [STAND_FOAM, STAND_NO_FOAM, SIT] -> choose what makes sense for you
data = [pow_foam, pow_firm, pow_sit];

[p_friedman, tbl_friedman, stats_friedman] = friedman(data, 1, 'off');

disp('Friedman test p-value:')
disp(p_friedman)


%% Post-hoc pairwise comparisons (Wilcoxon signed-rank) + Bonferroni

% Pairwise Wilcoxon (right-tailed, as before)
[p_Foam_vs_NoFoam, h_Foam_vs_NoFoam, stats_Foam_vs_NoFoam] = ...
    signrank(pow_foam,   pow_firm, 'tail', 'right');

[p_Foam_vs_Sit, h_Foam_vs_Sit, stats_Foam_vs_Sit] = ...
    signrank(pow_foam,   pow_sit,         'tail', 'right');

[p_NoFoam_vs_Sit, h_NoFoam_vs_Sit, stats_NoFoam_vs_Sit] = ...
    signrank(pow_firm, pow_sit,         'tail', 'right');

% Bonferroni correction for 3 pairwise tests
m = 3;             % number of comparisons
alpha = 0.05;
alphaBonf = alpha / m;

p_raw  = [p_Foam_vs_NoFoam, p_Foam_vs_Sit, p_NoFoam_vs_Sit];
p_corr = min(p_raw * m, 1);   % Bonferroni-adjusted p-values

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
bp(:,1) = pow_foam;
bp(:,2) = pow_firm;
bp(:,3) = pow_sit;

sText = 14;
figure;
% Create a boxplot with notches and custom labels
h = boxplot(bp, 'Notch','off', ...
              'Labels', {'STAND_FOAM','STAND_NO_FOAM','SIT'});

set(gca,'fontsize',sText);
set(h,{'linew'},{1});

% Color boxes (findobj returns them in reverse order of plotting)
c = findobj(gca,'Tag','Box');  % c(1) = last box, etc.

% Make sure indexing matches 3 conditions:
% c(1) -> SIT
% c(2) -> STAND_NO_FOAM
% c(3) -> STAND_FOAM
patch(get(c(3),'XData'), get(c(3),'YData'), [0.4860 0.6940 0.1680], 'FaceAlpha', .26); % STAND_FOAM
patch(get(c(2),'XData'), get(c(2),'YData'), [0.9290 0.6940 0.1250], 'FaceAlpha', .26); % STAND_NO_FOAM
patch(get(c(1),'XData'), get(c(1),'YData'), [0.8500 0.3250 0.0980], 'FaceAlpha', .26); % SIT

% Choose what to show in title; e.g., global Friedman p-value
if p_friedman < 0.0001
    p_str = '<0.0001';
else
    p_str = sprintf('= %.4f', p_friedman);
end

%title(['FFR power comparison [331–335 Hz] (Friedman p ', p_str, ')'])
title('FFR power comparison [331–335 Hz]')
ylabel('FFR power')
ylim([0 0.05])

hold on
parallelcoords(bp, 'Color', 0.7*[1 1 1], ...
                   'Linestyle', '-', 'Marker', '.', 'Markersize', 7);

%% Significance annotation (Bonferroni-corrected, main axes)
% y-position settings
yMax   = max(bp(:));
yMin   = min(bp(:));
yRange = yMax - yMin;
if yRange == 0
    yRange = 1; % avoid zero-range issues
end
yStep  = 0.02 * yRange;      % was 0.08
startY = yMax + 0.01 * yRange;   % was 0.10

% yStep  = 0.08 * yRange;
% startY = yMax + 0.1 * yRange;

% startY = prctile(bp(:), 97);   % start above highest whisker
% yStep  = 0.5 * std(bp(:));     % small consistent spacing


% Helper to convert p-value to star string (based on corrected p)
getStars = @(p) ...
    (p <= 0.001) * "***" + ...
    (p > 0.001 & p <= 0.01) * "**"  + ...
    (p > 0.01  & p <= 0.05) * "*"   + ...
    (p > 0.05) * "";

% IMPORTANT: easier implementation without anonymous trick (MATLAB char handling):
% We'll just do it explicitly below.

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
        'VerticalAlignment', 'bottom', ...
        'FontSize', sText);
end

% 2) STAND_FOAM (1) vs SIT (3)
if p_Foam_vs_Sit_corr <= alphaBonf
    y2 = startY + 1.6*yStep; % a bit higher to avoid overlap
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
        'VerticalAlignment', 'bottom', ...
        'FontSize', sText);
end

% 3) STAND_NO_FOAM (2) vs SIT (3)
if p_NoFoam_vs_Sit_corr <= alphaBonf
    y3 = startY + 3.2*yStep; % highest line
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
        'VerticalAlignment', 'bottom', ...
        'FontSize', sText);
end