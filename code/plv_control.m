%% PLV computation (control band version) — averages PLV across a frequency range
% Adapted from the original single-frequency PLV script. Computes PLV at every
% FFT bin in a specified range and averages, giving the mean phase-locking
% strength in that band. Useful as a control analysis to rule out generic
% phase-locking inflation from contamination at non-stimulus frequencies.
%
% For the 270-290 Hz control band with fs=10000 and NFFT=10000, this covers
% bins 271-291 (21 bins at 1 Hz spacing). PLV is computed per bin per polarity,
% averaged across polarities first, then averaged across bins.

clc; clear; close all
 
addpath 'C:\Users\Aline Moorkens\OneDrive - KU Leuven\Documents\fieldtrip-20250106'
ft_defaults

%% Configuration
cfg = struct();

% Paths
cfg.basepath  = 'F:\Auditory-mov-data';
cfg.outpath   = 'F:\Auditory-mov-data\LMM_PLV\control_band_270_290';
cfg.subfolder = 'preprocessed';

% Subjects
cfg.subjects   = {'NG', 'BVD', 'FD', 'EF', 'CX', 'HX', 'MCBG', 'LV', 'AM', ...
                  'FS', 'LVB', 'RNR', 'VUR', 'LH', 'TF', 'MGZ', 'NSMT', ...
                  'MM', 'TH', 'TD'};

% Data file naming
cfg.trial_fname_pattern       = 'preprocessed_P[%s]_FFR_%s_%d.mat';
cfg.trial_fname_pattern_sit2  = 'preprocessed_P[%s]_FFR_sit_%d_2.mat';
cfg.trial_var = 'data_final';

% Polarity
cfg.polarity_labels = {'FFR', '-FFR'};
cfg.exclude_labels  = {'notrial'};

% Channels
cfg.channels      = 'vertex';
cfg.vertex_chans  = {'Cz','C1','C2','CPz','FC1','FC2','CP1','CP2'};

cfg.fs = 10000;

% --- PLV computation settings ---
cfg.NFFT = 10000;       % gives df=1 Hz
cfg.toi  = [0.030 0.230];

% --- Control band specification ---
% Frequency range over which to compute and average PLV.
% This is the key change from the original script.
cfg.f_range = [270 290];   % Hz, inclusive

% Aggregation across bins: 'mean' or 'median'. Mean is more standard;
% median is more robust if a single bin within the range has an artifact.
cfg.bin_aggregation = 'mean';

% --- Note on stimulus frequency in the combo table ---
% We still need a "Frequency" column in the output table to match MLM.m's
% schema, but the value represents the *centre* of the control band rather
% than a true stimulus frequency. Set to the band midpoint.
cfg.control_band_label = round(mean(cfg.f_range));   % 280 Hz

% Sig-channel restriction
cfg.restrict_to_sigchans = true;

% Which combinations to compute
cfg.run.A_within_333  = true;   % sit/firm/foam at the 270-290 control band
cfg.run.B_prepost_333 = false;
cfg.run.C_within_89   = false;
cfg.run.D_between_freq= false;

% File-name templates
cfg.condition     = {'sit','firm','foam'};
cfg.prepost_names = {'sit','sit_2'};

% Output
if ~exist(cfg.outpath,'dir'), mkdir(cfg.outpath); end
cfg.plv_matfile   = fullfile(cfg.outpath, 'PLV_table_control_band.mat');
cfg.plv_xlsxsheet = 'plv_control';
cfg.write_xlsx    = true;
cfg.xlsx_file     = fullfile(cfg.outpath, 'FFR_PLV_control_band.xlsx');

cfg.verbose_missing = false;


%% Build the combo list
% IMPORTANT: when loading the trial-level data, we still need to load the
% files associated with the actual stimulus frequency (333 Hz), because
% the trial-level data exists there. The control-band PLV is computed from
% those same trials, just at different FFT bins.
%
% So the combo table tracks two things:
%   - load_freq:   the stimulus frequency used to find the data file (333 Hz)
%   - report_freq: the label written into the output table (280 Hz, the band centre)

need_333_condition = cfg.run.A_within_333 || cfg.run.B_prepost_333 || cfg.run.D_between_freq;
need_sit2          = cfg.run.B_prepost_333;
need_89_condition  = cfg.run.C_within_89  || cfg.run.D_between_freq;

% Columns: {condition, load_freq, report_freq, session}
combos_template = {};

if need_333_condition
    for k = 1:numel(cfg.condition)
        combos_template(end+1,:) = {cfg.condition{k}, 333, cfg.control_band_label, 'pre'};
    end
end
if need_sit2
    combos_template(end+1,:) = {'sit_2', 333, cfg.control_band_label, 'post'};
end
if need_89_condition
    for k = 1:numel(cfg.condition)
        combos_template(end+1,:) = {cfg.condition{k}, 89, cfg.control_band_label, 'pre'};
    end
end

n_combos = size(combos_template, 1);
n_sub    = numel(cfg.subjects);

fprintf('Computing PLV in control band %d-%d Hz (n_subjects=%d, n_combos=%d)\n', ...
    cfg.f_range(1), cfg.f_range(2), n_sub, n_combos);

missing_count = zeros(n_combos, 1);


%% Precompute the FFT bin range
df = cfg.fs / cfg.NFFT;
k_lo = round(cfg.f_range(1)/df) + 1;
k_hi = round(cfg.f_range(2)/df) + 1;
k_band = k_lo:k_hi;
fprintf('Control band covers FFT bins %d-%d (%d bins, df=%.2f Hz).\n', ...
    k_lo, k_hi, numel(k_band), df);


%% Main loop

rows = {};
rowi = 0;

for s = 1:n_sub
    subj   = cfg.subjects{s};
    inpath = fullfile(cfg.basepath, subj, cfg.subfolder);

    for r = 1:n_combos
        cond        = combos_template{r,1};
        load_freq   = combos_template{r,2};
        report_freq = combos_template{r,3};
        sess        = combos_template{r,4};

        % Build filename
        if strcmp(cond, 'sit_2')
            fname = fullfile(inpath, ...
                sprintf(cfg.trial_fname_pattern_sit2, subj, load_freq));
        else
            fname = fullfile(inpath, ...
                sprintf(cfg.trial_fname_pattern, subj, cond, load_freq));
        end

        cond_factor = cond;
        if strcmp(cond, 'sit_2'), cond_factor = 'sit'; end

        if ~exist(fname, 'file')
            if cfg.verbose_missing
                warning('Missing: %s', fname);
            end
            missing_count(r) = missing_count(r) + 1;
            continue
        end

        S = load(fname);
        if ~isfield(S, cfg.trial_var)
            if cfg.verbose_missing
                warning('Variable "%s" not found in %s', cfg.trial_var, fname);
            end
            missing_count(r) = missing_count(r) + 1;
            continue
        end
        data = S.(cfg.trial_var);

        % Channel selection
        if ischar(cfg.channels) && strcmpi(cfg.channels, 'all')
            ch_idx = 1:numel(data.label);
        elseif ischar(cfg.channels) && strcmpi(cfg.channels, 'vertex')
            [~, ch_idx] = ismember(cfg.vertex_chans, data.label);
            ch_idx = ch_idx(ch_idx > 0);
        else
            [~, ch_idx] = ismember(cfg.channels, data.label);
            ch_idx = ch_idx(ch_idx > 0);
        end

        ref_chans = {'TP9','TP10','Trigger'};
        is_ref = ismember(data.label, ref_chans);
        ch_idx = ch_idx(~is_ref(ch_idx));
        chan_labels = data.label(ch_idx);

        if isempty(ch_idx)
            warning('No channels left after selection for %s/%s', subj, cond);
            continue
        end

        % Polarity handling
        if ~isfield(data, 'trialinfo') || isempty(data.trialinfo)
            warning('No trialinfo for %s/%s — skipping.', subj, cond);
            continue
        end
        pol_labels = string(data.trialinfo);

        % --- Compute PLV per channel, per polarity, averaged across band ---
        n_chan = numel(ch_idx);
        n_pol  = numel(cfg.polarity_labels);
        n_band = numel(k_band);

        % plv_per_pol_per_bin: n_chan x n_pol x n_band
        plv_per_pol_per_bin = nan(n_chan, n_pol, n_band);
        n_trials_per_pol    = zeros(1, n_pol);

        for pidx = 1:n_pol
            pol_str = string(cfg.polarity_labels{pidx});
            trial_mask = (pol_labels == pol_str);
            n_trials_per_pol(pidx) = sum(trial_mask);

            if n_trials_per_pol(pidx) < 5
                continue
            end

            trial_idx = find(trial_mask);

            % complex_coeffs: n_chan x n_trials x n_band
            complex_coeffs = nan(n_chan, numel(trial_idx), n_band);

            for t = 1:numel(trial_idx)
                ti = trial_idx(t);
                tr = data.trial{ti}(ch_idx, :);
                tt = data.time{ti};

                if ~isempty(cfg.toi)
                    t_mask = tt >= cfg.toi(1) & tt <= cfg.toi(2);
                    tr = tr(:, t_mask);
                end

                F = fft(tr, cfg.NFFT, 2);
                complex_coeffs(:, t, :) = F(:, k_band);   % grab whole band at once
            end

            % Unit-normalize and average across trials, per bin
            unit_vectors = complex_coeffs ./ abs(complex_coeffs);
            unit_vectors(isnan(unit_vectors)) = 0;

            % PLV: |mean over trials| → result is n_chan x 1 x n_band
            plv_chan_bin = abs(mean(unit_vectors, 2));
            plv_per_pol_per_bin(:, pidx, :) = plv_chan_bin;
        end

        % Average PLV across polarities (n_chan x n_band)
        plv_per_chan_per_bin = squeeze(mean(plv_per_pol_per_bin, 2, 'omitnan'));

        % Aggregate across bins in the control band
        switch lower(cfg.bin_aggregation)
            case 'mean'
                plv_per_chan = mean(plv_per_chan_per_bin, 2, 'omitnan');
            case 'median'
                plv_per_chan = median(plv_per_chan_per_bin, 2, 'omitnan');
            otherwise
                error('Unknown bin_aggregation: %s', cfg.bin_aggregation);
        end

        total_trials = sum(n_trials_per_pol);

        % Sig-channel restriction
        % NB: sig channels were defined for the *stimulus* frequency, not the
        % control band. They're still useful as a channel selector (same vertex
        % set the LMM uses), but conceptually they no longer mean "significant
        % response at this frequency". You may want to set
        % cfg.restrict_to_sigchans = false for the control analysis, or keep
        % it on for direct comparability with the stimulus-frequency LMM.
        is_sig_here = true(numel(chan_labels), 1);
        if cfg.restrict_to_sigchans
            if strcmp(cond, 'sit_2')
                name_part = sprintf('_FFR_sit_%d_2', load_freq);
            else
                name_part = sprintf('_FFR_%s_%d', cond, load_freq);
            end
            sigfile = fullfile(inpath, ['SigChansPadded_P[', subj, ']', name_part, '.mat']);
            if exist(sigfile, 'file')
                Ssig = load(sigfile, 'sigChansBoth');
                if isfield(Ssig, 'sigChansBoth') && ~isempty(Ssig.sigChansBoth)
                    is_sig_here = ismember(chan_labels, Ssig.sigChansBoth);
                end
            end
        end

        % Append rows
        for c = 1:numel(chan_labels)
            rowi = rowi + 1;
            rows(rowi, :) = { string(subj), string(cond_factor), report_freq, ...
                              string(sess), string(chan_labels{c}), ...
                              plv_per_chan(c), total_trials, ...
                              is_sig_here(c) };
        end

        fprintf('  %s / %-6s / %d Hz band (load %d Hz) / %s : %d trials, mean PLV = %.3f\n', ...
            subj, cond, report_freq, load_freq, sess, total_trials, ...
            mean(plv_per_chan, 'omitnan'));
    end
end


%% Assemble and save the table

T_plv = cell2table(rows, 'VariableNames', ...
    {'Subject','Condition','Frequency','Session','Channel','PLV','NumTrials','IsSigHere'});

T_plv.Subject   = categorical(T_plv.Subject);
T_plv.Condition = categorical(T_plv.Condition);
T_plv.Session   = categorical(T_plv.Session);
T_plv.Channel   = categorical(T_plv.Channel);
T_plv.IsSigHere = logical(T_plv.IsSigHere);

save(cfg.plv_matfile, 'T_plv');
fprintf('\nControl-band PLV table saved to %s\n', cfg.plv_matfile);
fprintf('  %d rows, %d unique subjects, %d unique channels.\n', ...
    height(T_plv), numel(unique(T_plv.Subject)), numel(unique(T_plv.Channel)));

if cfg.write_xlsx
    writetable(T_plv, cfg.xlsx_file, 'Sheet', cfg.plv_xlsxsheet, ...
        'WriteMode', 'overwritesheet');
    fprintf('Also written to %s (sheet: %s)\n', cfg.xlsx_file, cfg.plv_xlsxsheet);
end


%% Missing-file summary

fprintf('\n--- Missing-file summary (out of %d subjects) ---\n', n_sub);
for r = 1:n_combos
    fprintf('  %-6s load=%d Hz report=%d Hz %-4s : %d/%d subjects missing\n', ...
        combos_template{r,1}, combos_template{r,2}, combos_template{r,3}, ...
        combos_template{r,4}, missing_count(r), n_sub);
end


%% Sanity-check summary

fprintf('\n--- Median PLV per Condition / Frequency / Session ---\n');
G = findgroups(T_plv.Condition, T_plv.Frequency, T_plv.Session);
[grp_cond, grp_freq, grp_sess] = findgroups_to_labels(T_plv, G);
plv_medians = splitapply(@(x) median(x, 'omitnan'), T_plv.PLV, G);
plv_iqr_lo  = splitapply(@(x) quantile(x, 0.25), T_plv.PLV, G);
plv_iqr_hi  = splitapply(@(x) quantile(x, 0.75), T_plv.PLV, G);
n_obs       = splitapply(@(x) numel(x), T_plv.PLV, G);

for k = 1:numel(plv_medians)
    fprintf('  %-6s %3d Hz %-4s : median PLV = %.3f, IQR = [%.3f, %.3f], n = %d\n', ...
        grp_cond(k), grp_freq(k), grp_sess(k), ...
        plv_medians(k), plv_iqr_lo(k), plv_iqr_hi(k), n_obs(k));
end

fprintf('\nDone.\n');


%% Helper

function [condL, freqL, sessL] = findgroups_to_labels(T, G)
[~, idx] = unique(G, 'stable');
condL = T.Condition(idx);
freqL = T.Frequency(idx);
sessL = T.Session(idx);
end