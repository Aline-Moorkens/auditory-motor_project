%% NOISE FLOOR CHECK — trial-level data (single-participant pilot)
% Computes per-channel, per-condition broadband EEG power in clean
% frequency bands away from the FFR stimulus and away from mains harmonics.
% Tests whether sit / firm / foam have systematically different noise floors
% that could bias condition-level SNR comparisons.
%
% Operates on the UNAVERAGED trial-level data because polarity subtraction
% in the averaged waveform deliberately suppresses uncorrelated noise — the
% contamination you're trying to detect.

clc; clear; close all

restoredefaultpath
addpath 'C:\Users\Aline Moorkens\OneDrive - KU Leuven\Documents\fieldtrip-20250106'
ft_defaults

%% Data input settings
subjectid  = 'AM';
conditions = {'sit', 'firm', 'foam'};
stim_freq  = 333;                                   % 333, 433, or 89

% input
inputpath = fullfile('F:\Auditory-mov-data', subjectid, 'preprocessed');
fname_var = 'data_final';                                        % Name of the FieldTrip variable inside each .mat file.
outdir = 'F:\Auditory-mov-data\noisefloor_check';           % Output directory for results and figures.

% Channels to use. For FFR, restrict to your FFR-relevant channels.
% Set to 'all' to include every channel (will dilute the signal).
channels = {'Cz', 'FCz', 'Fz'};


%% Noise bands — automatically chosen based on stim_freq, edit if needed
% Bands avoid the stimulus and fall BETWEEN mains harmonics (50, 100, ...).
% Edit these if your line frequency is 60 Hz instead of 50 Hz.

switch stim_freq
    case 333
        % Stimulus 333. Mains harmonics nearby: 300, 350.
        % Use 260-290 (between 250 and 300) and 360-390 (between 350 and 400).
        noise_bands = [260 290; 360 390];
    case 433
        % Stimulus 433. Mains harmonics nearby: 400, 450.
        noise_bands = [360 390; 460 490];
    case 89
        % Stimulus 89. Mains harmonics nearby: 50, 100.
        noise_bands = [60 80; 105 125];
    otherwise
        error('Unrecognized stim_freq %g — set noise_bands manually.', stim_freq);
end

if ~exist(outdir, 'dir'), mkdir(outdir); end

n_sub  = 1;
n_cond = numel(conditions);
n_band = size(noise_bands, 1);

fprintf('\n=== Broadband noise floor check (trial-level, single participant: %s) ===\n', subjectid);
fprintf('Stimulus frequency: %g Hz\n', stim_freq);
fprintf('Noise bands: ');
for b = 1:n_band
    fprintf('[%g-%g] ', noise_bands(b,1), noise_bands(b,2));
end
fprintf('Hz\n\n');


%% Pre-allocate result containers
% Per-channel array is allocated lazily once we know the channel count
% (after the first successful file load).

noise_power         = nan(n_sub, n_cond, n_band);   % (sub, cond, band) averaged over channels
noise_power_perchan = [];                            % (sub, cond, band, chan) — lazy alloc
chanlabels          = {};
freq_axis           = [];
spectra_all         = [];                            % (sub, cond, freq) for diagnostic plot


%% Main loop — load each condition, compute spectrum, extract band powers
% For one subject this is really just a loop over conditions. The outer
% subject loop is kept so the structure scales to many subjects later.

for s = 1:n_sub
    for c = 1:n_cond

        % --- Build filename for this condition ---
        % Matches: preprocessed_P[<subjectid>]_FFR_<condition>_<stim>.mat
        fname = fullfile(inputpath, ...
            sprintf('preprocessed_P[%s]_FFR_%s_%d.mat', subjectid, conditions{c}, stim_freq));

        if ~exist(fname, 'file')
            warning('Missing file: %s — skipping.', fname);
            continue;
        end

        % --- Load FieldTrip trial structure ---
        tmp = load(fname);
        if ~isfield(tmp, fname_var)
            % If you don't know what the variable is called, uncomment the
            % next line on first run to see what's inside:
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

        % --- Spectrum: mtmfft with Hanning taper, averaged across trials ---
        % keeptrials='no' returns mean power across trials, one spectrum per channel.
        % pad='nextpow2' zero-pads for stable frequency bins.
        cfg_ft = [];
        cfg_ft.method     = 'mtmfft';
        cfg_ft.taper      = 'hanning';
        cfg_ft.output     = 'pow';
        cfg_ft.foilim     = [10 max(noise_bands(:))+50];
        cfg_ft.keeptrials = 'no';
        cfg_ft.pad        = 'nextpow2';
        freq = ft_freqanalysis(cfg_ft, data);
        % freq.powspctrm is [chan x freq]

        % --- Lazy-init the per-channel array now that we know the channel count ---
        if isempty(noise_power_perchan)
            chanlabels          = freq.label;
            noise_power_perchan = nan(n_sub, n_cond, n_band, numel(chanlabels));
            freq_axis           = freq.freq;
            spectra_all         = nan(n_sub, n_cond, numel(freq_axis));
        end

        % --- Handle channel mismatch (in case some files have different montages) ---
        if ~isequal(freq.label, chanlabels)
            warning('Channel mismatch for %s / %s — aligning by intersection.', ...
                subjectid, conditions{c});
            [~, ia, ib] = intersect(chanlabels, freq.label, 'stable');
            tmp_pow = nan(numel(chanlabels), numel(freq.freq));
            tmp_pow(ia, :) = freq.powspctrm(ib, :);
            freq.powspctrm = tmp_pow;
        end

        % --- Save full spectrum (averaged across channels) for the diagnostic plot ---
        spectra_all(s, c, :) = mean(freq.powspctrm, 1, 'omitnan');

        % --- Extract mean power in each noise band ---
        for b = 1:n_band
            f_idx = freq.freq >= noise_bands(b,1) & freq.freq <= noise_bands(b,2);
            if ~any(f_idx)
                warning('No frequency bins in band [%g-%g] for %s/%s', ...
                    noise_bands(b,1), noise_bands(b,2), subjectid, conditions{c});
                continue;
            end
            pow_per_chan = mean(freq.powspctrm(:, f_idx), 2, 'omitnan');
            noise_power_perchan(s, c, b, :) = pow_per_chan;
            noise_power(s, c, b)            = mean(pow_per_chan, 'omitnan');
        end

        fprintf('  %s / %-6s done.\n', subjectid, conditions{c});
    end
end


%% Save results
% Stored variables can be loaded back later without re-running the loop.

results = struct();
results.subjectid           = subjectid;
results.conditions          = conditions;
results.noise_bands         = noise_bands;
results.stim_freq           = stim_freq;
results.channels            = channels;
results.noise_power         = noise_power;
results.noise_power_perchan = noise_power_perchan;
results.chanlabels          = chanlabels;
results.freq_axis           = freq_axis;
results.spectra_all         = spectra_all;

save(fullfile(outdir, sprintf('noise_floor_results_%s.mat', subjectid)), 'results', '-v7.3');
fprintf('\nResults saved to %s\n', fullfile(outdir, sprintf('noise_floor_results_%s.mat', subjectid)));


%% Figure 1 — noise power by condition (bar chart for one subject)
% With one subject there's no distribution to box, so this is just a bar
% chart per band. When you scale up to multiple subjects, switch back to
% the boxplot version.

fig1 = figure('Position', [100 100 350*n_band 400], 'Color', 'w');
for b = 1:n_band
    subplot(1, n_band, b);

    pow_band = squeeze(noise_power(:, :, b));   % (1, cond) for one subject
    bar(pow_band);
    set(gca, 'XTickLabel', conditions);

    title(sprintf('Noise %g-%g Hz', noise_bands(b,1), noise_bands(b,2)));
    ylabel('Power (\muV^2/Hz)');
    set(gca, 'YScale', 'log');   % power spans orders of magnitude
    grid on;
end
sgtitle(sprintf('Noise floor by condition — %s (stim %g Hz)', subjectid, stim_freq));
saveas(fig1, fullfile(outdir, sprintf('noise_floor_by_condition_%s.png', subjectid)));


%% Figure 2 — mean spectrum per condition
% Shows where the energy actually lives. Mains harmonics and the stimulus
% frequency are marked; the noise bands are shaded for visual verification.

fig2 = figure('Position', [100 550 1000 400], 'Color', 'w');
hold on;
colors = lines(n_cond);

for c = 1:n_cond
    mean_spec = squeeze(spectra_all(1, c, :));   % only one subject
    plot(freq_axis, mean_spec, 'Color', colors(c, :), 'LineWidth', 1.5, ...
        'DisplayName', conditions{c});
end

% Mark stimulus
yl = ylim;
plot([stim_freq stim_freq], yl, 'k--', 'HandleVisibility', 'off');
text(stim_freq, yl(2)*0.9, sprintf(' stim %g Hz', stim_freq), ...
    'HorizontalAlignment', 'left');

% Mark mains harmonics
for h = 50:50:max(freq_axis)
    plot([h h], yl, 'Color', [0.8 0.8 0.8], 'LineStyle', ':', 'HandleVisibility', 'off');
end

% Shade noise bands
for b = 1:n_band
    patch([noise_bands(b,1) noise_bands(b,2) noise_bands(b,2) noise_bands(b,1)], ...
          [yl(1) yl(1) yl(2) yl(2)], [1 0.9 0.7], ...
          'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

set(gca, 'YScale', 'log');
xlabel('Frequency (Hz)');
ylabel('Power (\muV^2/Hz)');
title(sprintf('Mean spectrum per condition — %s (averaged across channels)', subjectid));
legend('Location', 'best');
grid on;
saveas(fig2, fullfile(outdir, sprintf('mean_spectrum_per_condition_%s.png', subjectid)));


%% Summary — just print the noise band powers
% For one subject there's no statistical test to run (signrank needs paired
% observations across subjects). The single-subject result is just a
% sanity check that the contamination shows up at all.

fprintf('\n=== Noise band powers for %s ===\n', subjectid);
for b = 1:n_band
    fprintf('\nBand %g-%g Hz:\n', noise_bands(b,1), noise_bands(b,2));
    for c = 1:n_cond
        fprintf('  %-6s: %.4e\n', conditions{c}, noise_power(1, c, b));
    end

    % Quick ratio summary — useful for spotting asymmetry at a glance
    if n_cond == 3
        fprintf('  Ratios: firm/sit = %.2f, foam/sit = %.2f, foam/firm = %.2f\n', ...
            noise_power(1,2,b)/noise_power(1,1,b), ...
            noise_power(1,3,b)/noise_power(1,1,b), ...
            noise_power(1,3,b)/noise_power(1,2,b));
    end
end


%% Asymmetry summary — quantify the foam-vs-plate-contact difference
% The contamination hypothesis predicts foam to be cleaner than sit/firm
% because the foam pad insulates the participant's feet from the plate
% surface. This block prints the ratios explicitly in both directions.

fprintf('\n=== Asymmetry summary (%s) ===\n', subjectid);

for b = 1:n_band
    fprintf('\nBand %g-%g Hz:\n', noise_bands(b,1), noise_bands(b,2));

    sit_pow  = noise_power(1, 1, b);
    firm_pow = noise_power(1, 2, b);
    foam_pow = noise_power(1, 3, b);

    % Plate-contact conditions (sit, firm) vs insulated condition (foam)
    plate_mean = mean([sit_pow, firm_pow]);  % average noise across plate-contact conditions
    ratio = plate_mean / foam_pow;

    fprintf('  sit:  %.4e\n', sit_pow);
    fprintf('  firm: %.4e\n', firm_pow);
    fprintf('  foam: %.4e\n', foam_pow);
    fprintf('  mean(sit, firm) / foam = %.2fx\n', ratio);

    % How much of the plate-contact noise is "extra" relative to foam?
    excess_pct = 100 * (plate_mean - foam_pow) / foam_pow;
    fprintf('  Plate-contact conditions have %.1f%% higher noise than foam.\n', excess_pct);

    % Interpretation guide
    if ratio < 1.5
        fprintf('  Interpretation: small asymmetry, unlikely to dominate SNR.\n');
    elseif ratio < 3
        fprintf('  Interpretation: moderate asymmetry, worth checking averaged data.\n');
    elseif ratio < 10
        fprintf('  Interpretation: large asymmetry, definitely check averaged data.\n');
    else
        fprintf('  Interpretation: very large asymmetry, likely to bias SNR contrasts.\n');
    end
end

fprintf('\nDone.\n');
