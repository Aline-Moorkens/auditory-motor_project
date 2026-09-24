%% NOISE FLOOR + SNR CHECK — averaged FFR waveforms (multi-subject)
% Loops over all subjects, for each subject and condition:
%   - Loads the polarity-subtracted averaged FFR waveform
%   - Computes the power spectrum
%   - Extracts the broadband noise floor in the user-defined bands
%   - Computes the SNR at the stimulus frequency, using narrow flanking
%     noise bands (per-condition local noise estimate — the standard FFR
%     SNR definition)
% Saves per-subject figures + asymmetry summary as .txt (one per subject).
% Saves a group-level summary figure and .txt with Wilcoxon signed-rank
% tests on both the noise floor and the SNR.
%
% The SNR computation is the key analysis: if SNR(foam) > SNR(firm) at the
% group level after normalizing each condition by its own local noise
% floor, the effect is robust to the contamination identified in the
% noise-floor analysis.

clc; clear; close all

% restoredefaultpath
% addpath 'C:\Users\Aline Moorkens\OneDrive - KU Leuven\Documents\fieldtrip-20250106'
% ft_defaults


%% Data input settings

subjects = {'AM' 
'LVB'
'TH'
'TF' 
'NSMT'
'MGZ'
'LH'
'MM'
'FS'
'TD'
'RNR'
'HX'
'VUR'
'FD'
'EF'
'MCBG'
'LV'
'CX'
'NG'
'BVD'
};

conditions = {'sit', 'firm', 'foam'};
stim_freq  = 333;                                   % 333, 433, or 89

% Input: <basepath>/<subjectid>/preprocessed/averagedFFR_P[<subjectid>]_FFR_<cond>_<stim>.mat
basepath = 'F:\Auditory-mov-data';

% Name of the FieldTrip variable inside each .mat file.
fname_var = 'final';

% Output directory.
outdir = 'F:\Auditory-mov-data\noisefloor_check\SNR';

% Channels to use.
channels = {'Cz', 'FCz', 'Fz'};

% Time window of the averaged waveform to use for the spectrum.
time_window = [0 0.2];


%% Noise bands (broadband check) — automatically chosen based on stim_freq

switch stim_freq
    case 333
        noise_bands = [260 290; 360 390];
    case 433
        noise_bands = [360 390; 460 490];
    case 89
        noise_bands = [60 80; 105 125];
    otherwise
        error('Unrecognized stim_freq %g — set noise_bands manually.', stim_freq);
end


%% SNR definition — narrow signal and flanking noise bands
% Matches your existing SNR definition:
%   Signal band: 331-335 Hz (around 333 Hz FFR)
%   Noise bands: 321-325 Hz and 341-345 Hz (5 Hz bands, 6 Hz gap from signal)
% Edit if you want different definitions for 89 or 433 Hz stim.

switch stim_freq
    case 333
        snr_signal_band = [331 335];
        snr_noise_bands = [321 325; 341 345];
    case 433
        snr_signal_band = [431 435];
        snr_noise_bands = [421 425; 441 445];
    case 89
        snr_signal_band = [87 91];
        snr_noise_bands = [77 81; 97 101];
    otherwise
        error('No SNR bands defined for stim_freq %g — set manually.', stim_freq);
end


%% Setup

if ~exist(outdir, 'dir'), mkdir(outdir); end

n_sub  = numel(subjects);
n_cond = numel(conditions);
n_band = size(noise_bands, 1);

fprintf('\n=== Noise floor + SNR check (averaged FFR, %d subjects) ===\n', n_sub);
fprintf('Stimulus frequency: %g Hz\n', stim_freq);
fprintf('Broadband noise bands: ');
for b = 1:n_band
    fprintf('[%g-%g] ', noise_bands(b,1), noise_bands(b,2));
end
fprintf('Hz\n');
fprintf('SNR signal band: [%g-%g] Hz\n', snr_signal_band(1), snr_signal_band(2));
fprintf('SNR noise bands: [%g-%g] and [%g-%g] Hz\n', ...
    snr_noise_bands(1,1), snr_noise_bands(1,2), ...
    snr_noise_bands(2,1), snr_noise_bands(2,2));
if ~isempty(time_window)
    fprintf('Time window: [%.3f %.3f] s\n', time_window(1), time_window(2));
end
fprintf('\n');


%% Pre-allocate result containers (group-level)

noise_power      = nan(n_sub, n_cond, n_band);   % (sub, cond, band)
sig_power        = nan(n_sub, n_cond);            % (sub, cond) signal power at stim freq
snr_noise_power  = nan(n_sub, n_cond);            % (sub, cond) mean power across SNR noise bands
snr_values       = nan(n_sub, n_cond);            % (sub, cond) signal / noise = SNR
chanlabels       = {};
freq_axis        = [];
spectra_all      = [];                             % (sub, cond, freq)
ratios_all       = nan(n_sub, n_band);             % (sub, band) — (sit+firm)/foam noise ratio


%% Main loop

for s = 1:n_sub
    subjectid = subjects{s};
    inputpath = fullfile(basepath, subjectid, 'preprocessed');

    fprintf('--- Processing subject: %s ---\n', subjectid);

    sub_noise_power     = nan(n_cond, n_band);
    sub_sig_power       = nan(n_cond, 1);
    sub_snr_noise_power = nan(n_cond, 1);
    sub_snr_values      = nan(n_cond, 1);

    for c = 1:n_cond

        % --- Build filename ---
        fname = fullfile(inputpath, ...
            sprintf('averagedFFR_P[%s]_FFR_%s_%d.mat', subjectid, conditions{c}, stim_freq));

        if ~exist(fname, 'file')
            warning('Missing file: %s — skipping.', fname);
            continue;
        end

        % --- Load ---
        tmp = load(fname);
        if ~isfield(tmp, fname_var)
            warning('Variable "%s" not found in %s — skipping.', fname_var, fname);
            continue;
        end
        data = tmp.(fname_var);

        % --- Channel selection ---
        if ~strcmp(channels, 'all')
            sel = []; sel.channel = channels;
            data = ft_selectdata(sel, data);
        end

        % --- Extract waveform (timelock or raw) ---
        if isfield(data, 'avg')
            ffr_wave  = data.avg;
            ffr_time  = data.time;
            ffr_label = data.label;
        elseif isfield(data, 'trial') && iscell(data.trial)
            if numel(data.trial) ~= 1
                warning('Expected 1 trial for %s/%s, got %d. Using the first.', ...
                    subjectid, conditions{c}, numel(data.trial));
            end
            ffr_wave  = data.trial{1};
            ffr_time  = data.time{1};
            ffr_label = data.label;
        else
            warning('Unrecognized data format for %s/%s — skipping.', subjectid, conditions{c});
            continue;
        end

        % --- Time window restriction ---
        if ~isempty(time_window)
            t_idx = ffr_time >= time_window(1) & ffr_time <= time_window(2);
            ffr_wave = ffr_wave(:, t_idx);
            ffr_time = ffr_time(t_idx);
        end

        % --- Wrap as one-trial raw struct ---
        raw = [];
        raw.label   = ffr_label;
        raw.fsample = round(1 / mean(diff(ffr_time)));
        raw.trial   = {ffr_wave};
        raw.time    = {ffr_time};

        % --- Spectrum ---
        % Extend foilim a bit beyond max to ensure we capture all band edges.
        max_band = max([noise_bands(:); snr_noise_bands(:); snr_signal_band(:)]);
        cfg_ft = [];
        cfg_ft.method     = 'mtmfft';
        cfg_ft.taper      = 'hanning';
        cfg_ft.output     = 'pow';
        cfg_ft.foilim     = [10 max_band + 50];
        cfg_ft.keeptrials = 'no';
        cfg_ft.pad        = 'nextpow2';
        freq = ft_freqanalysis(cfg_ft, raw);

        % --- Lazy-init group spectrum array ---
        if isempty(freq_axis)
            chanlabels  = freq.label;
            freq_axis   = freq.freq;
            spectra_all = nan(n_sub, n_cond, numel(freq_axis));
        end

        % --- Channel alignment safeguard ---
        if ~isequal(freq.label, chanlabels)
            warning('Channel mismatch for %s / %s — aligning by intersection.', ...
                subjectid, conditions{c});
            [~, ia, ib] = intersect(chanlabels, freq.label, 'stable');
            tmp_pow = nan(numel(chanlabels), numel(freq.freq));
            tmp_pow(ia, :) = freq.powspctrm(ib, :);
            freq.powspctrm = tmp_pow;
        end

        % --- Save spectrum (averaged across channels) ---
        spectra_all(s, c, :) = mean(freq.powspctrm, 1, 'omitnan');

        % --- Broadband noise floor in each noise band ---
        for b = 1:n_band
            f_idx = freq.freq >= noise_bands(b,1) & freq.freq <= noise_bands(b,2);
            if ~any(f_idx), continue; end
            pow_per_chan = mean(freq.powspctrm(:, f_idx), 2, 'omitnan');
            sub_noise_power(c, b) = mean(pow_per_chan, 'omitnan');
        end

        % --- SNR computation ---
        % Signal: mean power across signal band, then across channels.
        % Noise: mean power across both flanking noise bands, then across channels.
        % SNR = signal / noise. Computed PER condition with PER condition local noise.
        sig_idx = freq.freq >= snr_signal_band(1) & freq.freq <= snr_signal_band(2);
        n1_idx  = freq.freq >= snr_noise_bands(1,1) & freq.freq <= snr_noise_bands(1,2);
        n2_idx  = freq.freq >= snr_noise_bands(2,1) & freq.freq <= snr_noise_bands(2,2);

        if any(sig_idx) && any(n1_idx) && any(n2_idx)
            sig_pow_per_chan = mean(freq.powspctrm(:, sig_idx), 2, 'omitnan');
            n1_pow_per_chan  = mean(freq.powspctrm(:, n1_idx),  2, 'omitnan');
            n2_pow_per_chan  = mean(freq.powspctrm(:, n2_idx),  2, 'omitnan');

            sig_pow = mean(sig_pow_per_chan, 'omitnan');
            % Average power across both flanking noise bands
            noise_pow = mean([mean(n1_pow_per_chan, 'omitnan'), mean(n2_pow_per_chan, 'omitnan')]);

            sub_sig_power(c)       = sig_pow;
            sub_snr_noise_power(c) = noise_pow;
            sub_snr_values(c)      = sig_pow / noise_pow;
        end

        fprintf('  %s / %-6s done.\n', subjectid, conditions{c});
    end

    % Push to group arrays
    noise_power(s, :, :)   = sub_noise_power;
    sig_power(s, :)        = sub_sig_power;
    snr_noise_power(s, :)  = sub_snr_noise_power;
    snr_values(s, :)       = sub_snr_values;

    % Compute (sit+firm)/foam ratio per band (noise floor)
    for b = 1:n_band
        plate_mean = mean(sub_noise_power(1:2, b), 'omitnan');
        ratios_all(s, b) = plate_mean / sub_noise_power(3, b);
    end


    % =================================================================
    % Per-subject outputs (figures + asymmetry text summary)
    % =================================================================

    % --- Figure 1: noise power by condition ---
    fig1 = figure('Position', [100 100 350*n_band 400], 'Color', 'w', 'Visible', 'off');
    for b = 1:n_band
        subplot(1, n_band, b);
        bar(sub_noise_power(:, b));
        set(gca, 'XTickLabel', conditions);
        title(sprintf('Noise %g-%g Hz', noise_bands(b,1), noise_bands(b,2)));
        ylabel('Power (\muV^2/Hz)');
        set(gca, 'YScale', 'log');
        grid on;
    end
    sgtitle(sprintf('%s — broadband noise floor (stim %g Hz)', subjectid, stim_freq));
    saveas(fig1, fullfile(outdir, sprintf('%s_noise_floor_by_condition.png', subjectid)));
    close(fig1);

    % --- Figure 2: mean spectrum per condition ---
    fig2 = figure('Position', [100 550 1000 400], 'Color', 'w', 'Visible', 'off');
    hold on;
    colors = lines(n_cond);
    for c = 1:n_cond
        mean_spec = squeeze(spectra_all(s, c, :));
        plot(freq_axis, mean_spec, 'Color', colors(c, :), 'LineWidth', 1.5, ...
            'DisplayName', conditions{c});
    end
    yl = ylim;
    plot([stim_freq stim_freq], yl, 'k--', 'HandleVisibility', 'off');
    text(stim_freq, yl(2)*0.9, sprintf(' stim %g Hz', stim_freq), 'HorizontalAlignment', 'left');
    for h = 50:50:max(freq_axis)
        plot([h h], yl, 'Color', [0.8 0.8 0.8], 'LineStyle', ':', 'HandleVisibility', 'off');
    end
    % Shade broadband noise bands
    for b = 1:n_band
        patch([noise_bands(b,1) noise_bands(b,2) noise_bands(b,2) noise_bands(b,1)], ...
              [yl(1) yl(1) yl(2) yl(2)], [1 0.9 0.7], ...
              'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
    set(gca, 'YScale', 'log');
    xlabel('Frequency (Hz)');
    ylabel('Power (\muV^2/Hz)');
    title(sprintf('%s — mean spectrum per condition', subjectid));
    legend('Location', 'best');
    grid on;
    saveas(fig2, fullfile(outdir, sprintf('%s_mean_spectrum_per_condition.png', subjectid)));
    close(fig2);

    % --- Figure 3: SNR by condition ---
    % Single bar chart — directly visible whether SNR differs by condition.
    fig3 = figure('Position', [100 100 500 400], 'Color', 'w', 'Visible', 'off');
    bar(sub_snr_values);
    set(gca, 'XTickLabel', conditions);
    ylabel('SNR (signal / local noise)');
    title(sprintf('%s — SNR at %g Hz by condition', subjectid, stim_freq));
    grid on;
    saveas(fig3, fullfile(outdir, sprintf('%s_SNR_by_condition.png', subjectid)));
    close(fig3);

    % --- Asymmetry + SNR summary: print to console AND save to .txt ---
    txt_path = fullfile(outdir, sprintf('%s_asymmetry_summary.txt', subjectid));
    fid = fopen(txt_path, 'w');

    fprintf('\n=== Summary for %s ===\n', subjectid);
    fprintf(fid, '=== Summary for %s ===\n', subjectid);

    % Broadband noise floor asymmetry
    fprintf('\n-- Broadband noise floor --\n');
    fprintf(fid, '\n-- Broadband noise floor --\n');
    for b = 1:n_band
        sit_pow  = sub_noise_power(1, b);
        firm_pow = sub_noise_power(2, b);
        foam_pow = sub_noise_power(3, b);
        plate_mean = mean([sit_pow, firm_pow]);
        ratio = plate_mean / foam_pow;
        excess_pct = 100 * (plate_mean - foam_pow) / foam_pow;

        if ratio < 1.2
            interp_str = 'asymmetry essentially gone. Polarity subtraction worked.';
        elseif ratio < 1.5
            interp_str = 'small residual asymmetry, probably acceptable.';
        elseif ratio < 3
            interp_str = 'moderate residual asymmetry, may bias raw-power SNR contrasts.';
        else
            interp_str = 'substantial residual asymmetry.';
        end

        for fid_target = [1, fid]   % 1 = stdout
            fprintf(fid_target, '\nBand %g-%g Hz:\n', noise_bands(b,1), noise_bands(b,2));
            fprintf(fid_target, '  sit:  %.4e\n', sit_pow);
            fprintf(fid_target, '  firm: %.4e\n', firm_pow);
            fprintf(fid_target, '  foam: %.4e\n', foam_pow);
            fprintf(fid_target, '  mean(sit, firm) / foam = %.2fx\n', ratio);
            fprintf(fid_target, '  Plate-contact conditions have %.1f%% higher noise than foam.\n', excess_pct);
            fprintf(fid_target, '  Interpretation: %s\n', interp_str);
        end
    end

    % SNR section
    fprintf('\n-- SNR at %g Hz (signal: %g-%g, noise: %g-%g + %g-%g) --\n', ...
        stim_freq, snr_signal_band(1), snr_signal_band(2), ...
        snr_noise_bands(1,1), snr_noise_bands(1,2), ...
        snr_noise_bands(2,1), snr_noise_bands(2,2));
    fprintf(fid, '\n-- SNR at %g Hz (signal: %g-%g, noise: %g-%g + %g-%g) --\n', ...
        stim_freq, snr_signal_band(1), snr_signal_band(2), ...
        snr_noise_bands(1,1), snr_noise_bands(1,2), ...
        snr_noise_bands(2,1), snr_noise_bands(2,2));

    for fid_target = [1, fid]
        for c = 1:n_cond
            fprintf(fid_target, '  %-6s: signal = %.4e, noise = %.4e, SNR = %.2f\n', ...
                conditions{c}, sub_sig_power(c), sub_snr_noise_power(c), sub_snr_values(c));
        end
        if ~any(isnan(sub_snr_values))
            fprintf(fid_target, '  SNR(foam) / SNR(firm) = %.2fx\n', ...
                sub_snr_values(3) / sub_snr_values(2));
            fprintf(fid_target, '  SNR(foam) / SNR(sit)  = %.2fx\n', ...
                sub_snr_values(3) / sub_snr_values(1));
            fprintf(fid_target, '  SNR(firm) / SNR(sit)  = %.2fx\n', ...
                sub_snr_values(2) / sub_snr_values(1));
        end
    end

    fclose(fid);
    fprintf('\nSaved per-subject outputs for %s.\n', subjectid);
end


%% Save the combined group-level results

results = struct();
results.subjects        = subjects;
results.conditions      = conditions;
results.noise_bands     = noise_bands;
results.snr_signal_band = snr_signal_band;
results.snr_noise_bands = snr_noise_bands;
results.stim_freq       = stim_freq;
results.channels        = channels;
results.time_window     = time_window;
results.noise_power     = noise_power;
results.sig_power       = sig_power;
results.snr_noise_power = snr_noise_power;
results.snr_values      = snr_values;
results.spectra_all     = spectra_all;
results.freq_axis       = freq_axis;
results.ratios_all      = ratios_all;
results.note            = 'Polarity-subtracted averaged FFR. Group-level. Includes SNR.';

save(fullfile(outdir, 'group_noise_floor_results.mat'), 'results', '-v7.3');
fprintf('\nGroup results saved to %s\n', fullfile(outdir, 'group_noise_floor_results.mat'));


%% Group summary figure

fig_group = figure('Position', [50 50 1500 1000], 'Color', 'w');

% --- Row 1: (sit+firm)/foam noise-floor ratio per band ---
for b = 1:n_band
    subplot(4, n_band, b);
    hold on;

    boxplot(ratios_all(:, b), 'Labels', {sprintf('%g-%g Hz', noise_bands(b,1), noise_bands(b,2))}, ...
        'Colors', 'k', 'Symbol', '');
    jitter = 0.05 * randn(n_sub, 1);
    scatter(1 + jitter, ratios_all(:, b), 30, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.6);
    plot(xlim, [1 1], 'k:', 'HandleVisibility', 'off');

    ylabel('(sit+firm) / foam ratio');
    title(sprintf('Noise floor ratio: %g-%g Hz', noise_bands(b,1), noise_bands(b,2)));
    grid on;
end

% --- Row 2: per-subject noise power across conditions, one panel per band ---
for b = 1:n_band
    subplot(4, n_band, n_band + b);
    hold on;
    pow_band = squeeze(noise_power(:, :, b));
    for s = 1:n_sub
        if all(~isnan(pow_band(s, :)))
            plot(1:n_cond, pow_band(s, :), '-', 'Color', [0.5 0.5 0.5 0.3]);
            plot(1:n_cond, pow_band(s, :), 'o', 'MarkerSize', 4, ...
                'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', 'none');
        end
    end
    boxplot(pow_band, 'Labels', conditions, 'Colors', 'k', 'Symbol', '', 'Whisker', 1.5);
    set(gca, 'YScale', 'log');
    ylabel('Power (\muV^2/Hz)');
    title(sprintf('Noise floor: %g-%g Hz', noise_bands(b,1), noise_bands(b,2)));
    grid on;
end

% --- Row 3: SNR per subject across conditions ---
subplot(4, n_band, (2*n_band + 1):(3*n_band));
hold on;
for s = 1:n_sub
    if all(~isnan(snr_values(s, :)))
        plot(1:n_cond, snr_values(s, :), '-', 'Color', [0.5 0.5 0.5 0.3]);
        plot(1:n_cond, snr_values(s, :), 'o', 'MarkerSize', 4, ...
            'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', 'none');
    end
end
boxplot(snr_values, 'Labels', conditions, 'Colors', 'k', 'Symbol', '', 'Whisker', 1.5);
plot(xlim, [1 1], 'k:', 'HandleVisibility', 'off');
ylabel('SNR (signal / local noise)');
title(sprintf('SNR at %g Hz per subject per condition', stim_freq));
grid on;

% --- Row 4: mean spectrum across subjects, per condition ---
subplot(4, n_band, (3*n_band + 1):(4*n_band));
hold on;
colors = lines(n_cond);
for c = 1:n_cond
    mean_spec = squeeze(mean(spectra_all(:, c, :), 1, 'omitnan'));
    plot(freq_axis, mean_spec, 'Color', colors(c, :), 'LineWidth', 1.5, ...
        'DisplayName', conditions{c});
end
yl = ylim;
plot([stim_freq stim_freq], yl, 'k--', 'HandleVisibility', 'off');
text(stim_freq, yl(2)*0.9, sprintf(' stim %g Hz', stim_freq), 'HorizontalAlignment', 'left');
for h = 50:50:max(freq_axis)
    plot([h h], yl, 'Color', [0.8 0.8 0.8], 'LineStyle', ':', 'HandleVisibility', 'off');
end
for b = 1:n_band
    patch([noise_bands(b,1) noise_bands(b,2) noise_bands(b,2) noise_bands(b,1)], ...
          [yl(1) yl(1) yl(2) yl(2)], [1 0.9 0.7], ...
          'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
% Also shade SNR signal + noise bands in a different color
patch([snr_signal_band(1) snr_signal_band(2) snr_signal_band(2) snr_signal_band(1)], ...
      [yl(1) yl(1) yl(2) yl(2)], [0.7 0.9 1], ...
      'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
for b = 1:size(snr_noise_bands, 1)
    patch([snr_noise_bands(b,1) snr_noise_bands(b,2) snr_noise_bands(b,2) snr_noise_bands(b,1)], ...
          [yl(1) yl(1) yl(2) yl(2)], [0.9 0.7 1], ...
          'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
set(gca, 'YScale', 'log');
xlabel('Frequency (Hz)');
ylabel('Power (\muV^2/Hz)');
title(sprintf('Mean spectrum per condition (n=%d). Orange = broadband noise bands; blue = SNR signal band; pink = SNR noise bands.', n_sub));
legend('Location', 'best');
grid on;

sgtitle(sprintf('GROUP SUMMARY — noise floor & SNR (n=%d, stim %g Hz)', n_sub, stim_freq));

saveas(fig_group, fullfile(outdir, 'GROUP_summary.png'));
fprintf('Group summary figure saved to %s\n', fullfile(outdir, 'GROUP_summary.png'));


%% Group-level statistical summary

fprintf('\n=== GROUP-LEVEL SUMMARY (n=%d) ===\n', n_sub);

group_txt = fullfile(outdir, 'GROUP_summary.txt');
fid = fopen(group_txt, 'w');
fprintf(fid, '=== GROUP-LEVEL SUMMARY (n=%d) ===\n', n_sub);

% --- Noise floor stats ---
fprintf('\n-- Broadband noise floor --\n');
fprintf(fid, '\n-- Broadband noise floor --\n');

for b = 1:n_band
    fprintf('\nBand %g-%g Hz:\n', noise_bands(b,1), noise_bands(b,2));
    fprintf(fid, '\nBand %g-%g Hz:\n', noise_bands(b,1), noise_bands(b,2));

    pow_band = squeeze(noise_power(:, :, b));

    for c = 1:n_cond
        vals = pow_band(:, c);
        vals = vals(~isnan(vals));
        if isempty(vals), continue; end
        line_str = sprintf('  %-6s: median = %.4e, IQR = [%.4e, %.4e], n = %d', ...
            conditions{c}, median(vals), quantile(vals, 0.25), quantile(vals, 0.75), numel(vals));
        fprintf('%s\n', line_str);
        fprintf(fid, '%s\n', line_str);
    end

    r = ratios_all(:, b);
    r = r(~isnan(r));
    if ~isempty(r)
        line_str = sprintf('  Ratio (sit+firm)/foam: median = %.2f, IQR = [%.2f, %.2f]', ...
            median(r), quantile(r, 0.25), quantile(r, 0.75));
        fprintf('%s\n', line_str);
        fprintf(fid, '%s\n', line_str);
    end

    fprintf('  Pairwise Wilcoxon signed-rank:\n');
    fprintf(fid, '  Pairwise Wilcoxon signed-rank:\n');
    for ci = 1:n_cond
        for cj = ci+1:n_cond
            v1 = pow_band(:, ci);
            v2 = pow_band(:, cj);
            valid = ~isnan(v1) & ~isnan(v2);
            if sum(valid) >= 5
                try
                    p = signrank(v1(valid), v2(valid));
                catch
                    p = NaN;
                end
                md = median(v1(valid) - v2(valid));
                line_str = sprintf('    %-6s vs %-6s: median diff = %+.4e, p = %.4g, n = %d', ...
                    conditions{ci}, conditions{cj}, md, p, sum(valid));
                fprintf('%s\n', line_str);
                fprintf(fid, '%s\n', line_str);
            end
        end
    end
end

% --- SNR stats ---
% This is the analysis that matters: does SNR differ across conditions
% after each condition is normalized by its own local noise floor?

fprintf('\n-- SNR at %g Hz --\n', stim_freq);
fprintf(fid, '\n-- SNR at %g Hz --\n', stim_freq);

for c = 1:n_cond
    vals = snr_values(:, c);
    vals = vals(~isnan(vals));
    if isempty(vals), continue; end
    line_str = sprintf('  %-6s: median SNR = %.2f, IQR = [%.2f, %.2f], n = %d', ...
        conditions{c}, median(vals), quantile(vals, 0.25), quantile(vals, 0.75), numel(vals));
    fprintf('%s\n', line_str);
    fprintf(fid, '%s\n', line_str);
end

fprintf('\n  Pairwise Wilcoxon signed-rank on SNR (within-subject):\n');
fprintf(fid, '\n  Pairwise Wilcoxon signed-rank on SNR (within-subject):\n');
for ci = 1:n_cond
    for cj = ci+1:n_cond
        v1 = snr_values(:, ci);
        v2 = snr_values(:, cj);
        valid = ~isnan(v1) & ~isnan(v2);
        if sum(valid) >= 5
            try
                p = signrank(v1(valid), v2(valid));
            catch
                p = NaN;
            end
            md = median(v1(valid) - v2(valid));
            line_str = sprintf('    %-6s vs %-6s: median SNR diff = %+.3f, p = %.4g, n = %d', ...
                conditions{ci}, conditions{cj}, md, p, sum(valid));
            fprintf('%s\n', line_str);
            fprintf(fid, '%s\n', line_str);
        end
    end
end

% --- Interpretation guide ---
fprintf('\n  Interpretation guide:\n');
fprintf(fid, '\n  Interpretation guide:\n');
guide_lines = {
    '    If SNR(foam) > SNR(firm) significantly at the group level, the effect'
    '    is robust to the plate noise contamination. The SNR metric normalizes'
    '    out the broadband asymmetry visible in the noise floor analysis.'
    '    '
    '    If SNR(foam) is NOT significantly different from SNR(firm) after this'
    '    correction, the apparent foam>firm effect in raw power may have been'
    '    driven by the noise asymmetry rather than a real neurophysiological'
    '    effect, and would not survive the correction.'
};
for k = 1:numel(guide_lines)
    fprintf('%s\n', guide_lines{k});
    fprintf(fid, '%s\n', guide_lines{k});
end

fclose(fid);
fprintf('\nGroup summary text saved to %s\n', group_txt);
fprintf('\nDone.\n');


%%
% Median signal power per condition
fprintf('Median signal power at 333 Hz:\n');
for c = 1:n_cond
    vals = sig_power(:, c);
    vals = vals(~isnan(vals));
    fprintf('  %-6s: median = %.4e, IQR = [%.4e, %.4e]\n', ...
        conditions{c}, median(vals), quantile(vals, 0.25), quantile(vals, 0.75));
end

% Wilcoxon signed-rank on signal power
fprintf('\nWilcoxon signed-rank on signal power (within-subject):\n');
for ci = 1:n_cond
    for cj = ci+1:n_cond
        v1 = sig_power(:, ci);
        v2 = sig_power(:, cj);
        valid = ~isnan(v1) & ~isnan(v2);
        if sum(valid) >= 5
            p = signrank(v1(valid), v2(valid));
            md = median(v1(valid) - v2(valid));
            fprintf('  %-6s vs %-6s: median diff = %+.4e, p = %.4g\n', ...
                conditions{ci}, conditions{cj}, md, p);
        end
    end
end
