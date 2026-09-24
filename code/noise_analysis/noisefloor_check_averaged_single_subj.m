%% NOISE FLOOR CHECK — averaged FFR waveforms (single-participant pilot)
% Companion to noise_floor_check_single_subject.m. Where that script
% operates on trial-level data (preserving uncorrelated noise), this one
% operates on the polarity-subtracted averaged waveforms — the data that
% feeds your downstream FFR analysis. Comparing the two tells you whether
% any condition-level noise asymmetry survives polarity subtraction.
%
% INTERPRETATION:
%   trial-level shows asymmetry, averaged does NOT  -->  polarity
%       subtraction killed the contamination; final FFRs are clean
%   both show asymmetry  -->  contamination survived your pipeline and
%       may be biasing condition contrasts
%
% Run section by section (Ctrl+Enter) to inspect intermediate variables.

clc; clear; close all

% restoredefaultpath
% addpath 'C:\Users\Aline Moorkens\OneDrive - KU Leuven\Documents\fieldtrip-20250106'
% ft_defaults

%% Data input settings
subjectid  = 'AM';
conditions = {'sit', 'firm', 'foam'};
stim_freq  = 333;                                   % 333, 433, or 89

% Input: one file per condition for this subject, containing the
% pre-computed polarity-subtracted FFR waveform.
inputpath = fullfile('F:\Auditory-mov-data', subjectid, 'preprocessed');

% Name of the FieldTrip variable inside each .mat file.
% Change this if your variable is called something else.
fname_var = 'final';

% Output directory.
outdir = 'F:\Auditory-mov-data\noisefloor_check';

% Channels to use (same as trial-level script for consistency).
channels = {'Cz', 'FCz', 'Fz'};

% Time window of the averaged waveform to use for the spectrum.
% Cutting out the pre-stimulus baseline gives a cleaner FFR spectrum.
% Your protocol uses 30 ms prestim + 200 ms stim, so [0 0.2] uses just the
% stimulus-on portion. Set to [] to use the full waveform.
time_window = [0 0.2];


%% Noise bands — automatically chosen based on stim_freq, edit if needed
% Bands avoid the stimulus and fall BETWEEN mains harmonics (50, 100, ...).

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

if ~exist(outdir, 'dir'), mkdir(outdir); end

n_sub  = 1;
n_cond = numel(conditions);
n_band = size(noise_bands, 1);

fprintf('\n=== Broadband noise floor check (averaged FFR, single participant: %s) ===\n', subjectid);
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


%% Pre-allocate result containers

noise_power         = nan(n_sub, n_cond, n_band);
noise_power_perchan = [];
chanlabels          = {};
freq_axis           = [];
spectra_all         = [];


%% Main loop — load each condition, compute spectrum of the FFR waveform

for s = 1:n_sub
    for c = 1:n_cond

        % --- Build filename ---
        % Same naming convention as the trial-level data, just in a different
        % folder. Adjust the sprintf pattern if your averaged-file names
        % differ (e.g., 'averaged_P[...]' instead of 'preprocessed_P[...]').
        fname = fullfile(inputpath, ...
            sprintf('averagedFFR_P[%s]_FFR_%s_%d.mat', subjectid, conditions{c}, stim_freq));

        if ~exist(fname, 'file')
            warning('Missing file: %s — skipping.', fname);
            continue;
        end

        % --- Load FieldTrip timelock/raw structure ---
        tmp = load(fname);
        if ~isfield(tmp, fname_var)
            % If you don't know the variable name:
            % disp(fieldnames(tmp));
            warning('Variable "%s" not found in %s — skipping.', fname_var, fname);
            continue;
        end
        data = tmp.(fname_var);

        % --- Channel selection ---
        if ~strcmp(channels, 'all')
            sel = []; sel.channel = channels;
            data = ft_selectdata(sel, data);
        end

        % --- Extract the FFR waveform ---
        % FieldTrip can store averaged data as either a 'timelock' (with .avg
        % field, [chan x time]) or as 'raw' (with .trial as a cell array).
        % Handle both possibilities.
        if isfield(data, 'avg')
            % Timelock format
            ffr_wave = data.avg;      % [chan x time]
            ffr_time = data.time;     % [1 x time]
            ffr_label = data.label;
        elseif isfield(data, 'trial') && iscell(data.trial)
            % Raw format with (typically) one "trial" = the FFR waveform
            if numel(data.trial) ~= 1
                warning('Expected 1 trial in averaged data for %s/%s, got %d. Using the first.', ...
                    subjectid, conditions{c}, numel(data.trial));
            end
            ffr_wave = data.trial{1};
            ffr_time = data.time{1};
            ffr_label = data.label;
        else
            warning('Unrecognized data format for %s/%s — skipping.', subjectid, conditions{c});
            continue;
        end

        % --- Time window restriction ---
        % Cutting away the pre-stimulus baseline gives a cleaner FFR
        % spectrum since the baseline contributes only noise.
        if ~isempty(time_window)
            t_idx = ffr_time >= time_window(1) & ffr_time <= time_window(2);
            ffr_wave = ffr_wave(:, t_idx);
            ffr_time = ffr_time(t_idx);
        end

        % --- Wrap as FieldTrip raw data with one trial for ft_freqanalysis ---
        % This is the standard idiom for taking a spectrum of an already-
        % averaged waveform. mtmfft with keeptrials='no' on a one-trial
        % raw struct just returns that trial's spectrum.
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

        % --- Lazy-init result arrays ---
        if isempty(noise_power_perchan)
            chanlabels          = freq.label;
            noise_power_perchan = nan(n_sub, n_cond, n_band, numel(chanlabels));
            freq_axis           = freq.freq;
            spectra_all         = nan(n_sub, n_cond, numel(freq_axis));
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

        % --- Save full spectrum (averaged across channels) ---
        spectra_all(s, c, :) = mean(freq.powspctrm, 1, 'omitnan');

        % --- Extract mean power in each noise band ---
        for b = 1:n_band
            f_idx = freq.freq >= noise_bands(b,1) & freq.freq <= noise_bands(b,2);
            if ~any(f_idx), continue; end
            pow_per_chan = mean(freq.powspctrm(:, f_idx), 2, 'omitnan');
            noise_power_perchan(s, c, b, :) = pow_per_chan;
            noise_power(s, c, b)            = mean(pow_per_chan, 'omitnan');
        end

        fprintf('  %s / %-6s done.\n', subjectid, conditions{c});
    end
end


%% Save results

results = struct();
results.subjectid           = subjectid;
results.conditions          = conditions;
results.noise_bands         = noise_bands;
results.stim_freq           = stim_freq;
results.channels            = channels;
results.time_window         = time_window;
results.noise_power         = noise_power;
results.noise_power_perchan = noise_power_perchan;
results.chanlabels          = chanlabels;
results.freq_axis           = freq_axis;
results.spectra_all         = spectra_all;
results.note                = 'Spectra computed from polarity-subtracted averaged FFR waveforms.';

save(fullfile(outdir, sprintf('noise_floor_results_averaged_%s.mat', subjectid)), 'results', '-v7.3');
fprintf('\nResults saved to %s\n', fullfile(outdir, sprintf('noise_floor_results_averaged_%s.mat', subjectid)));


%% Figure 1 — noise power by condition (bar chart for one subject)

fig1 = figure('Position', [100 100 350*n_band 400], 'Color', 'w');
for b = 1:n_band
    subplot(1, n_band, b);

    pow_band = squeeze(noise_power(:, :, b));   % (1, cond) for one subject
    bar(pow_band);
    set(gca, 'XTickLabel', conditions);

    title(sprintf('Noise %g-%g Hz', noise_bands(b,1), noise_bands(b,2)));
    ylabel('Power (\muV^2/Hz)');
    set(gca, 'YScale', 'log');
    grid on;
end
sgtitle(sprintf('Averaged-data noise floor by condition — %s (stim %g Hz)', subjectid, stim_freq));
saveas(fig1, fullfile(outdir, sprintf('noise_floor_averaged_by_condition_%s.png', subjectid)));


%% Figure 2 — mean spectrum per condition (averaged waveform)

fig2 = figure('Position', [100 550 1000 400], 'Color', 'w');
hold on;
colors = lines(n_cond);

for c = 1:n_cond
    mean_spec = squeeze(spectra_all(1, c, :));
    plot(freq_axis, mean_spec, 'Color', colors(c, :), 'LineWidth', 1.5, ...
        'DisplayName', conditions{c});
end

yl = ylim;
plot([stim_freq stim_freq], yl, 'k--', 'HandleVisibility', 'off');
text(stim_freq, yl(2)*0.9, sprintf(' stim %g Hz', stim_freq), ...
    'HorizontalAlignment', 'left');

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
title(sprintf('Mean spectrum per condition (polarity-subtracted FFR) — %s', subjectid));
legend('Location', 'best');
grid on;
saveas(fig2, fullfile(outdir, sprintf('mean_spectrum_averaged_per_condition_%s.png', subjectid)));


%% Asymmetry summary — quantify what survives polarity subtraction
% The key comparison: if trial-level showed (sit ≈ firm) > foam, does the
% averaged data still show that, or did polarity subtraction kill it?

fprintf('\n=== Asymmetry summary (averaged data, %s) ===\n', subjectid);

for b = 1:n_band
    fprintf('\nBand %g-%g Hz:\n', noise_bands(b,1), noise_bands(b,2));

    sit_pow  = noise_power(1, 1, b);
    firm_pow = noise_power(1, 2, b);
    foam_pow = noise_power(1, 3, b);

    plate_mean = mean([sit_pow, firm_pow]);
    ratio = plate_mean / foam_pow;

    fprintf('  sit:  %.4e\n', sit_pow);
    fprintf('  firm: %.4e\n', firm_pow);
    fprintf('  foam: %.4e\n', foam_pow);
    fprintf('  mean(sit, firm) / foam = %.2fx\n', ratio);

    excess_pct = 100 * (plate_mean - foam_pow) / foam_pow;
    fprintf('  Plate-contact conditions have %.1f%% higher noise than foam.\n', excess_pct);

    if ratio < 1.2
        fprintf('  Interpretation: asymmetry essentially gone. Polarity subtraction worked.\n');
    elseif ratio < 1.5
        fprintf('  Interpretation: small residual asymmetry, probably acceptable.\n');
    elseif ratio < 3
        fprintf('  Interpretation: moderate residual asymmetry, may bias SNR contrasts.\n');
    else
        fprintf('  Interpretation: substantial residual asymmetry, likely problematic.\n');
    end
end

fprintf('\nDone.\n');
