%% FFR Subsampling - Significance as a function of trial count
% Sweeps trial counts, computes FFR SNR per channel, thresholds on 95th
% percentile of the noise distribution across iterations.

%% Set parameters
subjectid = 'AM';

% Input paths
inputpath = fullfile('F:\Auditory-mov-data', subjectid, 'preprocessed');
nameSit   = '_FFR_sit_333';
nameFirm  = '_FFR_firm_333';
nameFoam  = '_FFR_foam_333';

inputSit  = ['preprocessed_', 'P[', subjectid, ']', nameSit,  '.mat'];
inputFirm = ['preprocessed_', 'P[', subjectid, ']', nameFirm, '.mat'];
inputFoam = ['preprocessed_', 'P[', subjectid, ']', nameFoam, '.mat'];

%% Load data
load(fullfile(inputpath, inputSit),  'final'); ffr_sit  = final;
load(fullfile(inputpath, inputFirm), 'final'); ffr_firm = final;
load(fullfile(inputpath, inputFoam), 'final'); ffr_foam = final;

%% Analysis parameters
N_iterations    = 500;
trial_counts    = [500, 800, 1100, 1400, 1700, 2000, 2400, 2800, 3000];
freq_of_interest = 333;   % Hz
fs              = 10000;  % sampling rate (Hz)
freq_tol        = 3;      % ±3 Hz around FFR frequency = signal bins
noise_halfwidth = 20;     % ±20 Hz noise floor window (excluding signal bins)
noise_pct       = 95;     % percentile threshold for significance

% Choose which condition(s) to analyse — loop over all three or pick one
conditions      = {'sit', 'firm', 'foam'};
datasets        = {ffr_sit, ffr_firm, ffr_foam};


%% Main loop over conditions
for cond_idx = 1:numel(conditions)

    cond_name = conditions{cond_idx};
    data      = datasets{cond_idx};

    fprintf('\n=== Condition: %s ===\n', cond_name);

    % Polarity indices (only FFR and -FFR)
    % trialinfo is a cell array of strings: 'FFR', '-FFR', 'notrial'
    pos_idx = find(strcmp(data, 'FFR'));
    neg_idx = get_polarity_idx(data, '-FFR');


    % time axis & FFT setup (use first trial's time vector)
    t           = data.time{1};
    toi  = [0.03 0.230]
    N_samples   = numel(data.trial{1,1});
    N_chans     = numel(data.label);

    freqs       = (0:N_samples-1) * (fs / N_samples);

    signal_bins = (freqs >= freq_of_interest - freq_tol) & ...
                  (freqs <= freq_of_interest + freq_tol);
    noise_bins  = (freqs >= freq_of_interest - noise_halfwidth) & ...
                  (freqs <= freq_of_interest + noise_halfwidth) & ...
                  ~signal_bins;

    % --- Pre-extract trial matrices for speed ---
    % Shape: [N_chans x N_samples_stim x N_trials]
    all_pos_trials = cat(3, data.trial{pos_idx});
    all_neg_trials = cat(3, data.trial{neg_idx});
    all_pos_trials = all_pos_trials(:, toi:end, :);
    all_neg_trials = all_neg_trials(:, toi:end, :);

    % --- Preallocate results ---
    % SNR_matrix : [N_chans x N_trial_counts x N_iterations]
    SNR_matrix    = zeros(N_chans, numel(trial_counts), N_iterations);
    % noise_dist  : [N_chans x N_trial_counts x N_iterations]
    % (mean noise power per iteration, used for percentile threshold)
    noise_dist    = zeros(N_chans, numel(trial_counts), N_iterations);

    %% Sweep over trial counts
    for t_idx = 1:numel(trial_counts)
        N_FFR  = trial_counts(t_idx);
        N_draw = min(N_FFR, min(numel(pos_idx), numel(neg_idx)));

        if N_draw < N_FFR
            warning('Condition %s: only %d trials per polarity available (requested %d).', ...
                cond_name, N_draw, N_FFR);
        end

        fprintf('  Trial count: %d/%d ...\n', N_FFR, trial_counts(end));

        for iter = 1:N_iterations

            % Random subsets (indices into the already-extracted matrices)
            rp = randperm(numel(pos_idx), N_draw);
            rn = randperm(numel(neg_idx), N_draw);

            % Average each polarity: [N_chans x N_samples_stim]
            avg_pos = mean(all_pos_trials(:, :, rp), 3);
            avg_neg = mean(all_neg_trials(:, :, rn), 3);

            % FFR via polarity subtraction: cancels stimulus artefact
            ffr_signal = (avg_pos - avg_neg) / 2;   % [N_chans x N_samples_stim]

            % Power spectrum for every channel simultaneously
            fft_power = abs(fft(ffr_signal, [], 2)).^2 / N_samples;
            % fft_power: [N_chans x N_samples_stim]

            % Signal and noise power (mean across respective bins)
            sig_power   = mean(fft_power(:, signal_bins), 2);  % [N_chans x 1]
            noise_power = mean(fft_power(:, noise_bins),  2);  % [N_chans x 1]

            SNR_matrix(:, t_idx, iter) = 10 * log10(sig_power ./ noise_power);
            noise_dist(:, t_idx, iter) = noise_power;
        end
    end

    %% Significance threshold: 95th percentile of noise-based SNR distribution
    % Under H0 (no signal), SNR would reflect noise variability only.
    % We use the empirical 95th percentile across iterations as threshold.
    % threshold: [N_chans x N_trial_counts]
    threshold = prctile(SNR_matrix, noise_pct, 3);

    % Mean and CI across iterations
    mean_SNR = mean(SNR_matrix, 3);         % [N_chans x N_trial_counts]
    std_SNR  = std(SNR_matrix, 0, 3);
    ci95_SNR = 1.96 * std_SNR / sqrt(N_iterations);

    % Significant = mean SNR exceeds 95th percentile threshold
    is_significant = mean_SNR > threshold;  % [N_chans x N_trial_counts]

    
    %% Summary: first trial count at which each channel becomes significant
    min_trials_per_chan = NaN(N_chans, 1);
    for ch = 1:N_chans
        sig_tc = trial_counts(is_significant(ch, :));
        if ~isempty(sig_tc)
            min_trials_per_chan(ch) = sig_tc(1);
        end
    end

    fprintf('\n  Channels significant at each trial count:\n');
    for t_idx = 1:numel(trial_counts)
        n_sig = sum(is_significant(:, t_idx));
        fprintf('    %4d trials/polarity: %d / %d channels significant\n', ...
            trial_counts(t_idx), n_sig, N_chans);
    end

    %% Save results
    results.(cond_name).SNR_matrix        = SNR_matrix;
    results.(cond_name).mean_SNR          = mean_SNR;
    results.(cond_name).ci95_SNR          = ci95_SNR;
    results.(cond_name).threshold         = threshold;
    results.(cond_name).is_significant    = is_significant;
    results.(cond_name).min_trials_per_chan = min_trials_per_chan;
    results.(cond_name).chan_labels        = data.label;
    results.(cond_name).trial_counts       = trial_counts;


    %% Plot 1: SNR curve averaged across channels
    mean_SNR_across_chans     = mean(mean_SNR,  1);   % [1 x N_trial_counts]
    mean_thresh_across_chans  = mean(threshold, 1);
    ci95_across_chans         = mean(ci95_SNR,  1);

    figure('Name', sprintf('FFR SNR - %s - %s', subjectid, cond_name), ...
           'Units', 'normalized', 'Position', [0.05 0.1 0.55 0.45]);

    x = trial_counts;
    fill([x, fliplr(x)], ...
         [mean_SNR_across_chans + ci95_across_chans, ...
          fliplr(mean_SNR_across_chans - ci95_across_chans)], ...
         [0.7 0.85 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    hold on;
    plot(x, mean_SNR_across_chans,    'b-o', 'LineWidth', 2, 'MarkerSize', 6);
    plot(x, mean_thresh_across_chans, 'r--', 'LineWidth', 1.5);
    xlabel('Number of trials per polarity');
    ylabel('FFR SNR (dB)');
    title(sprintf('FFR SNR vs trial count | %s | %s | %d Hz', ...
        subjectid, cond_name, freq_of_interest));
    legend({'95% CI', 'Mean SNR (all channels)', ...
            sprintf('%dth pct threshold', noise_pct)}, 'Location', 'southeast');
    grid on; hold off;

end
