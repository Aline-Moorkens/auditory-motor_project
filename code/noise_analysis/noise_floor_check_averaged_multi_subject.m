%% NOISE FLOOR CHECK — averaged FFR waveforms (multi-subject)
% Loops over all subjects, computes the noise floor by condition from the
% polarity-subtracted averaged FFR waveforms, and saves per-subject results
% (figures + asymmetry summary as .txt) plus a group-level summary figure.
%
% Output naming convention: every per-subject output file is prefixed with
% the subject ID so the output directory stays flat and scrollable.
%
% Run section by section (Ctrl+Enter) to inspect intermediate variables.

clc; clear; close all

% restoredefaultpath
% addpath 'C:\Users\Aline Moorkens\OneDrive - KU Leuven\Documents\fieldtrip-20250106'
% ft_defaults


%% Data input settings

% List of subject IDs to process. Edit this to match your 20 participants.
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

% Output directory — all subjects' results go here, distinguished by prefix.
outdir = 'F:\Auditory-mov-data\noisefloor_check';

% Channels to use.
channels = {'Cz', 'FCz', 'Fz'};

% Time window of the averaged waveform to use for the spectrum.
% [0 0.2] = stimulus-on portion only (excludes 30 ms pre-stim baseline).
time_window = [0 0.2];


%% Noise bands — automatically chosen based on stim_freq, edit if needed

switch stim_freq
    case 333
       % noise_bands = [260 290; 360 390];
        noise_bands = [321 325; 341 345]; % more narrow bands that I used
    case 433
        noise_bands = [360 390; 460 490];
    case 89
        noise_bands = [60 80; 105 125];
    otherwise
        error('Unrecognized stim_freq %g — set noise_bands manually.', stim_freq);
end

if ~exist(outdir, 'dir'), mkdir(outdir); end

n_sub  = numel(subjects);
n_cond = numel(conditions);
n_band = size(noise_bands, 1);

fprintf('\n=== Broadband noise floor check (averaged FFR, %d subjects) ===\n', n_sub);
fprintf('Stimulus frequency: %g Hz\n', stim_freq);
fprintf('Noise bands: ');
for b = 1:n_band
    fprintf('[%g-%g] ', noise_bands(b,1), noise_bands(b,2));
end
fprintf('Hz\n');
if ~isempty(time_window)
    fprintf('Time window: [%.3f %.3f] s\n', time_window(1), time_window(2));
end
fprintf('\n');


%% Pre-allocate result containers (group-level)

noise_power = nan(n_sub, n_cond, n_band);   % (sub, cond, band) averaged across channels
chanlabels  = {};
freq_axis   = [];
spectra_all = [];                            % (sub, cond, freq)
ratios_all  = nan(n_sub, n_band);            % (sub, band) — (sit+firm)/foam per subject


%% Main loop — process each subject, save per-subject outputs

for s = 1:n_sub
    subjectid = subjects{s};
    inputpath = fullfile(basepath, subjectid, 'preprocessed');

    fprintf('--- Processing subject: %s ---\n', subjectid);

    % Per-subject storage (will be written back into the group arrays after)
    sub_noise_power = nan(n_cond, n_band);
    sub_spectra     = [];

    for c = 1:n_cond

        % --- Build filename ---
        fname = fullfile(inputpath, ...
            sprintf('averagedFFR_P[%s]_FFR_%s_%d.mat', subjectid, conditions{c}, stim_freq));

        if ~exist(fname, 'file')
            warning('Missing file: %s — skipping.', fname);
            continue;
        end

        % --- Load FieldTrip structure ---
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

        % --- Extract the FFR waveform (timelock or raw format) ---
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

        % --- Wrap as one-trial raw struct for ft_freqanalysis ---
        raw = [];
        raw.label   = ffr_label;
        raw.fsample = round(1 / mean(diff(ffr_time)));
        raw.trial   = {ffr_wave};
        raw.time    = {ffr_time};

        % --- Spectrum ---
        cfg_ft = [];
        cfg_ft.method     = 'mtmfft';
        cfg_ft.taper      = 'hanning';
        cfg_ft.output     = 'pow';
        cfg_ft.foilim     = [10 max(noise_bands(:))+50];
        cfg_ft.keeptrials = 'no';
        cfg_ft.pad        = 'nextpow2';
        freq = ft_freqanalysis(cfg_ft, raw);

        % --- Lazy-init group spectrum array on first successful load ---
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

        % --- Extract mean power in each noise band ---
        for b = 1:n_band
            f_idx = freq.freq >= noise_bands(b,1) & freq.freq <= noise_bands(b,2);
            if ~any(f_idx), continue; end
            pow_per_chan = mean(freq.powspctrm(:, f_idx), 2, 'omitnan');
            sub_noise_power(c, b) = mean(pow_per_chan, 'omitnan');
        end

        fprintf('  %s / %-6s done.\n', subjectid, conditions{c});
    end

    % Push subject results into the group array
    noise_power(s, :, :) = sub_noise_power;

    % Compute (sit+firm)/foam ratio per band
    for b = 1:n_band
        plate_mean = mean(sub_noise_power(1:2, b), 'omitnan');
        ratios_all(s, b) = plate_mean / sub_noise_power(3, b);
    end


    % =================================================================
    % Per-subject outputs (figures + asymmetry text summary)
    % =================================================================

    % --- Figure 1: noise power by condition (bar chart, log scale) ---
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
    sgtitle(sprintf('%s — averaged-data noise floor (stim %g Hz)', subjectid, stim_freq));
    saveas(fig1, fullfile(outdir, sprintf('%s_noise_floor_narrow_per_condition.png', subjectid)));
    % saveas(fig1, fullfile(outdir, sprintf('%s_noise_floor_per_condition.png', subjectid)));
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
    % saveas(fig2, fullfile(outdir, sprintf('%s_mean_spectrum_per_condition.png', subjectid)));
    saveas(fig2, fullfile(outdir, sprintf('%s_mean_spectrum_narrow_per_condition.png', subjectid)));
    close(fig2);

    % --- Asymmetry summary: print to console AND save to .txt ---
    txt_path = fullfile(outdir, sprintf('%s_asymmetry_summary.txt', subjectid));
    fid = fopen(txt_path, 'w');

    % Helper-style: write to both stdout and the file. We avoid a nested
    % function because the user prefers flat scripts; just call fprintf twice.

    fprintf('\n=== Asymmetry summary (averaged data, %s) ===\n', subjectid);
    fprintf(fid, '=== Asymmetry summary (averaged data, %s) ===\n', subjectid);

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
            interp_str = 'moderate residual asymmetry, may bias SNR contrasts.';
        else
            interp_str = 'substantial residual asymmetry, likely problematic.';
        end

        % Print to console
        fprintf('\nBand %g-%g Hz:\n', noise_bands(b,1), noise_bands(b,2));
        fprintf('  sit:  %.4e\n', sit_pow);
        fprintf('  firm: %.4e\n', firm_pow);
        fprintf('  foam: %.4e\n', foam_pow);
        fprintf('  mean(sit, firm) / foam = %.2fx\n', ratio);
        fprintf('  Plate-contact conditions have %.1f%% higher noise than foam.\n', excess_pct);
        fprintf('  Interpretation: %s\n', interp_str);

        % Print to file
        fprintf(fid, '\nBand %g-%g Hz:\n', noise_bands(b,1), noise_bands(b,2));
        fprintf(fid, '  sit:  %.4e\n', sit_pow);
        fprintf(fid, '  firm: %.4e\n', firm_pow);
        fprintf(fid, '  foam: %.4e\n', foam_pow);
        fprintf(fid, '  mean(sit, firm) / foam = %.2fx\n', ratio);
        fprintf(fid, '  Plate-contact conditions have %.1f%% higher noise than foam.\n', excess_pct);
        fprintf(fid, '  Interpretation: %s\n', interp_str);
    end

    fclose(fid);
    fprintf('\nSaved per-subject outputs for %s.\n', subjectid);
end


%% Save the combined group-level results

results = struct();
results.subjects     = subjects;
results.conditions   = conditions;
results.noise_bands  = noise_bands;
results.stim_freq    = stim_freq;
results.channels     = channels;
results.time_window  = time_window;
results.noise_power  = noise_power;     % (sub, cond, band)
results.spectra_all  = spectra_all;     % (sub, cond, freq)
results.freq_axis    = freq_axis;
results.ratios_all   = ratios_all;      % (sub, band) — (sit+firm)/foam
results.note         = 'Polarity-subtracted averaged FFR. Group-level.';

%save(fullfile(outdir, 'group_noise_floor_results.mat'), 'results', '-v7.3');
save(fullfile(outdir, 'group_noise_floor_narrow_band_results.mat'), 'results', '-v7.3');
fprintf('\nGroup results saved to %s\n', fullfile(outdir, 'group_noise_floor_results.mat'));


%% Group summary figure
% Three panels per band:
%   - Boxplot of the (sit+firm)/foam ratio across subjects (shows
%     consistency of the asymmetry across the group)
%   - Per-subject noise power across the three conditions with within-
%     subject lines (shows the raw pattern)
%   - Mean spectrum across subjects, per condition (shows the contamination
%     in frequency-resolved view)

fig_group = figure('Position', [50 50 1400 900], 'Color', 'w');

% --- Row 1: ratio boxplot per band ---
for b = 1:n_band
    subplot(3, n_band, b);
    hold on;

    boxplot(ratios_all(:, b), 'Labels', {sprintf('%g-%g Hz', noise_bands(b,1), noise_bands(b,2))}, ...
        'Colors', 'k', 'Symbol', '');
    % Overlay individual subjects as dots
    jitter = 0.05 * randn(n_sub, 1);
    scatter(1 + jitter, ratios_all(:, b), 30, [0.5 0.5 0.5], 'filled', ...
        'MarkerFaceAlpha', 0.6);

    % Reference line at ratio = 1 (no asymmetry)
    plot(xlim, [1 1], 'k:', 'HandleVisibility', 'off');

    ylabel('(sit+firm) / foam ratio');
    title(sprintf('Band %g-%g Hz', noise_bands(b,1), noise_bands(b,2)));
    grid on;
end

% --- Row 2: per-subject noise power across conditions, one panel per band ---
for b = 1:n_band
    subplot(3, n_band, n_band + b);
    hold on;
    pow_band = squeeze(noise_power(:, :, b));   % (sub, cond)
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
    title(sprintf('Noise %g-%g Hz', noise_bands(b,1), noise_bands(b,2)));
    grid on;
end

% --- Row 3: mean spectrum across subjects, per condition (spans bottom row) ---
subplot(3, n_band, (2*n_band + 1):(3*n_band));
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
set(gca, 'YScale', 'log');
xlabel('Frequency (Hz)');
ylabel('Power (\muV^2/Hz)');
title(sprintf('Mean spectrum per condition (across %d subjects)', n_sub));
legend('Location', 'best');
grid on;

sgtitle(sprintf('GROUP SUMMARY — averaged-data noise floor by condition (n=%d, stim %g Hz)', ...
    n_sub, stim_freq));

%saveas(fig_group, fullfile(outdir, 'GROUP_summary_narrow_band.png'));
saveas(fig_group, fullfile(outdir, 'GROUP_summary_narrow_band.png'));
fprintf('Group summary figure saved to %s\n', fullfile(outdir, 'GROUP_summary.png'));


%% Group-level statistical summary
% Wilcoxon signed-rank tests per band, across subjects. Tests whether
% the noise asymmetry is consistent at the group level.

fprintf('\n=== GROUP-LEVEL SUMMARY (n=%d) ===\n', n_sub);

% Also save group summary to a .txt
group_txt = fullfile(outdir, 'GROUP_summary_narrow-band_analysis.txt');
fid = fopen(group_txt, 'w');
fprintf(fid, '=== GROUP-LEVEL SUMMARY (n=%d) ===\n', n_sub);

for b = 1:n_band
    fprintf('\nBand %g-%g Hz:\n', noise_bands(b,1), noise_bands(b,2));
    fprintf(fid, '\nBand %g-%g Hz:\n', noise_bands(b,1), noise_bands(b,2));

    pow_band = squeeze(noise_power(:, :, b));   % (sub, cond)

    % Descriptives per condition
    for c = 1:n_cond
        vals = pow_band(:, c);
        vals = vals(~isnan(vals));
        if isempty(vals), continue; end
        line_str = sprintf('  %-6s: median = %.4e, IQR = [%.4e, %.4e], n = %d', ...
            conditions{c}, median(vals), quantile(vals, 0.25), quantile(vals, 0.75), numel(vals));
        fprintf('%s\n', line_str);
        fprintf(fid, '%s\n', line_str);
    end

    % Ratio descriptives across subjects
    r = ratios_all(:, b);
    r = r(~isnan(r));
    if ~isempty(r)
        line_str = sprintf('  Ratio (sit+firm)/foam across subjects: median = %.2f, IQR = [%.2f, %.2f]', ...
            median(r), quantile(r, 0.25), quantile(r, 0.75));
        fprintf('%s\n', line_str);
        fprintf(fid, '%s\n', line_str);
    end

    % Pairwise Wilcoxon signed-rank
    fprintf('  Pairwise Wilcoxon signed-rank (within-subject):\n');
    fprintf(fid, '  Pairwise Wilcoxon signed-rank (within-subject):\n');
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

fclose(fid);
fprintf('\nGroup summary text saved to %s\n', group_txt);
fprintf('\nDone.\n');
