%% Compare standing on foam vs standing no foam vs sitting
clc
clear
close all

%% Set parameters
subjectid = 'AM';
alpha = 0.05;
conditions = {'sit','firm','foam'};

% input
inputpath = fullfile('F:\Auditory-mov-data\control-no_plate\',subjectid,'preprocessed');

nameSit = '_FFR_sit_2_333';
nameFirm = '_FFR_firm_333';
nameFoam = '_FFR_foam_333';

inputSit = ['preprocessed_','P[',subjectid,']',nameSit,'.mat'];
inputFirm = ['preprocessed_','P[',subjectid,']',nameFirm,'.mat'];
inputFoam = ['preprocessed_','P[',subjectid,']',nameFoam,'.mat'];

input_sigchans_sit = ['SigChans_','P[',subjectid,']',nameSit,'.mat'];
input_sigchans_firm = ['Sigchans_','P[',subjectid,']',nameFirm,'.mat'];
input_sigchans_foam = ['Sigchans_','P[',subjectid,']',nameFoam,'.mat'];

%% load trial level data

load(fullfile(inputpath, inputSit),  'data_final'); ffr_sit  = data_final;
load(fullfile(inputpath, inputFirm), 'data_final'); ffr_firm = data_final;
load(fullfile(inputpath, inputFoam), 'data_final'); ffr_foam = data_final;

load(fullfile(inputpath, input_sigchans_sit),  'sigChansBoth'); sigchans_sit  = sigChansBoth;
load(fullfile(inputpath, input_sigchans_firm), 'sigChansBoth'); sigchans_firm  = sigChansBoth;
load(fullfile(inputpath, input_sigchans_firm), 'sigChansBoth'); sigchans_foam  = sigChansBoth;


%% PLV calculation

data_all   = {ffr_sit, ffr_firm, ffr_foam};

plv_all    = cell(1, 3);
labels_all = cell(1, 3);

for c = 1:numel(conditions)
    data   = data_all{c};
    labels = data_all{c}.label;

    nchan  = size(data.trial, 1);
    n_trials = size(data.trial, 2);
    NFFT = 10000;
    fs = data.fsample;
    f0 = 333;
    toi = [0.050 0.220];   % 20-180 ms post stimulus onset

    polarity_labels = {'FFR', '-FFR'};   % the two polarity labels to use
    exclude_labels  = {'notrial'};       % any labels to exclude entirely

    % Identify polarities --> see trialinfo
    % Convert trialinfo to a string array for uniform handling,
    % whether it came in as cell, char, or string.
    pol_labels = string(data.trialinfo);
    
    % --- Compute PLV per channel, per polarity ---
    plv_per_pol = zeros(nchan, numel(polarity_labels));
    n_trials_per_pol = zeros(1, numel(polarity_labels));% Frequency bin index for F0
    df = fs/NFFT;
    k_f0 = round(f0/df)+1;

    for pidx = 1:numel(polarity_labels)
        pol_str = string(polarity_labels{pidx});
        trial_mask = (pol_labels == pol_str);
        n_trials_per_pol(pidx) = sum(trial_mask);

         trial_idx = find(trial_mask);

         % Per-trial complex FFT coefficient at F0
         complex_coeffs = zeros(nchan, numel(trial_idx));

         for t = 1:numel(trial_idx)
             ti = trial_idx(t);
             tr = data.trial{ti}(ch_idx, :);
             tt = data.time{ti};

             if ~isempty(toi)
                 t_mask = tt >= toi(1) & tt <= toi(2);
                 tr = tr(:, t_mask);
             end

                F = fft(tr, NFFT, 2);
                complex_coeffs(:, t) = F(:, k_f0);
          end

         % PLV per channel = | mean over trials of (z / |z|) |
         unit_vectors = complex_coeffs ./ abs(complex_coeffs);
         unit_vectors(isnan(unit_vectors)) = 0;
         plv_per_pol(:, pidx) = abs(mean(unit_vectors, 2));
         
     end

        % Average PLV across polarities (ignoring missing polarities)
        plv_cond = mean(plv_per_pol, 2, 'omitnan');
        total_trials = sum(n_trials_per_pol);

        plv_all{c}    = plv_cond;
        labels_all{c} = labels;

        fprintf('PSD computed: %s | %d channels | %d trials\n', ...
        conditions{c}, nchan, n_trials);
end


plv_sit   = plv_all{1};   labels_sit  = labels_all{1};
plv_firm  = plv_all{2};   labels_firm = labels_all{2};
plv_foam  = plv_all{3};   labels_foam = labels_all{3};

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

plv_sit  = plv_sit(idx_sit,   :);
plv_firm = plv_firm(idx_firm,  :);
plv_foam = plv_foam(idx_foam,  :);

%% Only include channels with significant FFR
sig_labels = intersect(sigchans_sit, sigchans_firm, 'stable');
sig_labels = intersect(sig_labels, sigchans_foam,  'stable');

fprintf('Significant channels common to all conditions: %d\n', numel(sig_labels));

% Find indices of common channels in each condition's label list
idx_sig_sit  = find(ismember(chan_labels, sig_labels));
idx_sig_firm = find(ismember(chan_labels, sig_labels));
idx_sig_foam = find(ismember(chan_labels, sig_labels));

% Apply channel selection to PSD matrices
psd_sig_sit  = plv_sit(idx_sig_sit,:);
psd_sig_firm = plv_firm(idx_sig_firm,:);
psd_sig_foam = plv_foam(idx_sig_foam,:);

%% Non-parametric test for 3 related conditions: Friedman test
data = [plv_foam, plv_firm, plv_sit];

[p_friedman, tbl_friedman, stats_friedman] = friedman(data, 1, 'off');

disp('Friedman test p-value:')
disp(p_friedman)

%% Post-hoc pairwise comparisons (Wilcoxon signed-rank) + Bonferroni
[p_Foam_vs_NoFoam, h_Foam_vs_NoFoam, stats_Foam_vs_NoFoam] = ...
    signrank(plv_foam, plv_firm, 'tail', 'right');

[p_Foam_vs_Sit, h_Foam_vs_Sit, stats_Foam_vs_Sit] = ...
    signrank(plv_foam, plv_sit, 'tail', 'right');

[p_NoFoam_vs_Sit, h_NoFoam_vs_Sit, stats_NoFoam_vs_Sit] = ...
    signrank(plv_firm, plv_sit, 'tail', 'right');

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