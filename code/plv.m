%% PLV computation — produces a long-format table for direct use in MLM.m
% Computes phase-locking value (PLV) at the stimulus frequency per
% Subject x Condition x Frequency x Session x Channel
%
% PLV is computed on TRIAL-LEVEL data (not averaged waveform), because
% it measures phase consistency across trials. For each trial, the complex
% FFT coefficient at F0 is normalized to unit magnitude, and PLV is the
% magnitude of the mean of these unit vectors across trials. PLV ranges
% from 0 (no phase locking) to 1 (perfect phase locking).
%
% POLARITY HANDLING:
% FFR uses alternating polarity. Polarity 2 trials have inverted signal,
% which inverts the phase at F0. Naive pooling gives a bimodal phase
% distribution and a near-zero PLV. This script computes PLV separately
% per polarity and averages them, which is mathematically equivalent to
% sign-flipping and pooling and is the standard approach.
%
% OUTPUT:
% A long-format table with columns matching MLM.m's SNR table:
%   Subject, Condition, Frequency, Session, Channel, PLV, NumTrials
% Saved as a .mat file and (optionally) as an Excel sheet.
%
% HOW TO USE THE OUTPUT IN MLM.m:
% 1. Run this script first to generate PLV_table.mat.
% 2. In MLM.m, after the line `T = buildTable(cfg);`, add:
%
%       PLV = load(fullfile(cfg.outpath, 'PLV_table.mat'));
%       T = outerjoin(T, PLV.T_plv, ...
%           'Keys', {'Subject','Condition','Frequency','Session','Channel'}, ...
%           'MergeKeys', true, 'Type','left');
%
% 3. Then set cfg.outcome = 'PLV' to use PLV as the outcome variable.
%    The LMM formula becomes:  PLV ~ Condition + (1|Subject)
%
% Run section by section (Ctrl+Enter) to inspect intermediate variables.
close all
clear
clc

% addpath 'C:\Users\Aline Moorkens\OneDrive - KU Leuven\Documents\fieldtrip-20250106'
% ft_defaults

%% Configuration
cfg = struct();

% Paths
cfg.basepath  = 'F:\Auditory-mov-data';
%cfg.basepath  = 'F:\Auditory-mov-data\control-no_plate';
%cfg.basepath  = 'F:\Auditory-mov-data\control-no_plate\phase-analysis\plv_subsampled';

cfg.outpath   = 'F:\Auditory-mov-data\PLV_all';
%cfg.outpath   = 'F:\Auditory-mov-data\control-no_plate';
%cfg.outpath   = 'F:\Auditory-mov-data\control-no_plate\phase-analysis\plv_subsampled';

cfg.subfolder = 'preprocessed';

% Subjects
cfg.subjects  = {'NG', 'BVD', 'FD', 'EF', 'CX', 'HX', 'MCBG', 'LV',  ...
                 'AM', 'FS', 'LVB', 'RNR', 'VUR', 'LH', 'TF', 'MGZ', ...
                 'NSMT', 'MM', 'TH', 'TD'};

% data file naming
% Use %s for subject, %c for condition, and %d for frequency.
% If the file has a "_2" suffix for post-session, that's handled below.
cfg.trial_fname_pattern       = 'preprocessed_P[%s]_FFR_%s_%d.mat';     % subj, cond, freq
cfg.trial_fname_pattern_sit2  = 'preprocessed_P[%s]_FFR_sit_%d_2.mat';   % subj, freq
% cfg.trial_fname_pattern       = 'subsampled_P[%s]_FFR_%s_%d.mat';     % subj, cond, freq
% cfg.trial_fname_pattern_sit2  = 'subsampled_P[%s]_FFR_sit_%d_2.mat';   % subj, freq

% Name of the FieldTrip variable inside each .mat file.
cfg.trial_var = 'data_final';
% cfg.trial_var = 'data_sub';

% polarity
cfg.polarity_labels = {'FFR', '-FFR'};   % the two polarity labels to use

% channels
%cfg.channels      = 'vertex';
cfg.channels      = 'all';
%cfg.vertex_chans  = {'Cz','C1','C2','CPz','FC1','FC2','CP1','CP2'};

cfg.fs = 10000;     % fs

% --- PLV computation settings ---
% NFFT controls frequency-bin resolution. With fs=10000 and NFFT=10000,
% df=1 Hz, which puts a bin exactly at 333 Hz. For 89 Hz, the closest bin
% with NFFT=10000 is bin 90 (89 Hz exactly). Both stimulus frequencies are
% therefore covered cleanly with NFFT=10000.
cfg.NFFT = 10000;

% Time window of the trial to use for the FFT, in TRIAL time (t=0 = trigger,
% not stimulus onset - the stimulus itself starts 30 ms into the trial).
% This selects 20-180 ms POST STIMULUS ONSET, i.e. trigger-relative
% [0.030+0.020, 0.030+0.180] = [0.050, 0.210].
cfg.toi = [0.050 0.210];   % 20-180 ms post stimulus onset (trigger-relative)

% Sig-channel restriction
cfg.restrict_to_sigchans = true;

% Which combinations to compute
cfg.run.A_within_333  = true;
cfg.run.B_prepost_333 = true;
cfg.run.C_within_89   = true;

% File-name templates
cfg.condition     = {'sit','firm','foam'};
cfg.prepost_names = {'sit','sit_2'};

% Output handling
if ~exist(cfg.outpath,'dir'), mkdir(cfg.outpath); end
cfg.plv_matfile   = fullfile(cfg.outpath, 'PLV_table.mat');
cfg.plv_xlsxsheet = 'plv';   % sheet name in the same Excel file MLM uses
cfg.write_xlsx    = true;
cfg.xlsx_file     = fullfile(cfg.outpath, 'FFR_PLV_values.xlsx');

% Logging
cfg.verbose_missing = false;

%% Build the combo list
% This mirrors buildTable() in MLM.m so the resulting PLV table has rows
% that match the SNR table 1-to-1.

%need_333_condition = cfg.run.A_within_333 && cfg.run.B_prepost_333; %% In this case 333 Hz analysis is only ran when both sit-firm-foam AND pre-post are necessary.
need_333_condition = cfg.run.A_within_333 || cfg.run.B_prepost_333;
need_sit2          = cfg.run.B_prepost_333;
need_89_condition  = cfg.run.C_within_89;

combos_template = {};

if need_333_condition
    for k = 1:numel(cfg.condition)
        combos_template(end+1,:) = {cfg.condition{k}, 333, 'pre'}; %#ok<SAGROW>
    end
end
if need_sit2
    combos_template(end+1,:) = {'sit_2', 333, 'post'};
end
if need_89_condition
    for k = 1:numel(cfg.condition)
        combos_template(end+1,:) = {cfg.condition{k}, 89, 'pre'}; %#ok<SAGROW>
    end
end

n_combos = size(combos_template, 1);
n_sub    = numel(cfg.subjects);

fprintf('Computing PLV (n_subjects=%d, n_combos=%d)', n_sub, n_combos);

missing_count = zeros(n_combos, 1);


%% Main loop 
% For each subject x condition x frequency x session, load trial-level
% data, compute PLV per channel per polarity, average across polarities,
% and append rows to the long-format table.

rows = {};
rowi = 0;

for s = 1:n_sub
    subj   = cfg.subjects{s};
    inpath = fullfile(cfg.basepath, subj, cfg.subfolder);
    %inpath = fullfile(cfg.basepath);   % for subsampled plv

    for r = 1:n_combos
        cond = combos_template{r,1};
        f0   = combos_template{r,2};
        sess = combos_template{r,3};

        % Build filename
        % '_FFR_sit_<freq>_2' for sit_2 (post session)
        if strcmp(cond, 'sit_2')
            fname = fullfile(inpath, ...
                sprintf(cfg.trial_fname_pattern_sit2, subj, f0));
        else
            fname = fullfile(inpath, ...
                sprintf(cfg.trial_fname_pattern, subj, cond, f0));
        end


        % Map sit_2 -> sit for the condition factor (session carries pre/post)
        cond_factor = cond;
        if strcmp(cond, 'sit_2'), cond_factor = 'sit'; end

        if ~exist(fname, 'file')
            if cfg.verbose_missing
                warning('Missing: %s', fname);
            end
            missing_count(r) = missing_count(r) + 1;
            continue
        end

        % Load trial-level data
        S = load(fname);

        if ~isfield(S, cfg.trial_var)
            if cfg.verbose_missing
                warning('Variable "%s" not found in %s', cfg.trial_var, fname);
            end
            missing_count(r) = missing_count(r) + 1;
            continue
        end
        
        data = S.(cfg.trial_var); % we now have the trial level data

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

        % Identify polarities --> see trialinfo
        n_trials = numel(data.trial);

        if ~isfield(data, 'trialinfo') || isempty(data.trialinfo)
            warning('No trialinfo for %s/%s — skipping.', subj, cond);
            continue
        end

        % Convert trialinfo to a string array for uniform handling,
        % whether it came in as cell, char, or string.
        pol_labels = string(data.trialinfo);

        % % If trialinfo has multiple columns, take the first one (adjust
        % % here if your polarity label is in a different column).
        % if size(pol_labels, 2) > 1
        %     pol_labels = pol_labels(:, 1);
        % end
        % 
        % % Sanity check: report what labels we see (helpful first run)
        % unique_labels = unique(pol_labels);
        % unrecognized = setdiff(unique_labels, [string(cfg.polarity_labels),string(cfg.exclude_labels)]);
        % if ~isempty(unrecognized)
        %     warning('Unrecognized trialinfo labels in %s/%s: %s', ...
        %         subj, cond, strjoin(unrecognized, ', '));
        % end

        %%%% Compute PLV per channel, per polarity %%%%
        n_chan = numel(ch_idx);
        n_trials_per_pol = zeros(1, numel(cfg.polarity_labels));
        all_unit_vectors = [];   %must be reset here, once per subject/condition/frequency combo

        % Frequency bin index for F0
        df = cfg.fs/cfg.NFFT;
        k_f0 = round(f0/df)+1;

        for pidx = 1:numel(cfg.polarity_labels)
            pol_str = string(cfg.polarity_labels{pidx});
            trial_mask = (pol_labels == pol_str);
            n_trials_per_pol(pidx) = sum(trial_mask);

            if n_trials_per_pol(pidx) < 5
                % Not enough trials in this polarity to estimate PLV reliably
                continue
            end

            trial_idx = find(trial_mask);

            % Per-trial complex FFT coefficient at F0
            complex_coeffs = nan(n_chan, numel(trial_idx));

            for t = 1:numel(trial_idx)
                ti = trial_idx(t);
                tr = data.trial{ti}(ch_idx, :);
                tt = data.time{ti};

                if ~isempty(cfg.toi)
                    t_mask = tt >= cfg.toi(1) & tt <= cfg.toi(2);
                    tr = tr(:, t_mask);
                end

                % Apply a taper before the FFT, for consistency with the SNR script and
                % to reduce spectral leakage into the F0 bin (this reduces but does not eliminate leakage, 
                % since the window length is not an integer number of cycles at 333 Hz).
                win = hann(size(tr,2))';
                tr  = tr.*win;

                F = fft(tr, cfg.NFFT, 2);
                complex_coeffs(:,t) = F(:,k_f0);
            end

            % PLV per channel = | mean over trials of (z / |z|) |
            unit_vectors = complex_coeffs ./ abs(complex_coeffs);
            unit_vectors(isnan(unit_vectors)) = 0;

            % Sign-flip the inverted polarity BEFORE pooling, rather than
            % computing a separate PLV magnitude for this polarity alone. Matched
            % by the actual label string (not loop index), so this stays correct
            % even if cfg.polarity_labels is ever reordered.
            if strcmp(cfg.polarity_labels{pidx}, '-FFR')
                unit_vectors = -unit_vectors;
            end

            all_unit_vectors = [all_unit_vectors, unit_vectors]; %#ok<AGROW>
        end

        % ONE PLV computed from all trials, pooled after sign-correction - not
        % an average of two separately-computed magnitudes.
        plv_per_chan = abs(mean(all_unit_vectors, 2));
        total_trials = sum(n_trials_per_pol);

        % Sig-channel restriction 
        is_sig_here = true(numel(chan_labels), 1); % default when restriction is OFF entirely
        if cfg.restrict_to_sigchans
            % Once restriction is actually requested, default set to FALSE, not
            % TRUE --> a missing or empty sig-channel file should conservatively
            % exclude channels, not silently count them all as significant.
             is_sig_here = false(numel(chan_labels), 1);

            if strcmp(cond, 'sit_2')
                name_part = sprintf('_FFR_sit_%d_2', f0);
            else
                name_part = sprintf('_FFR_%s_%d', cond, f0);
            end
            
            sigfile = fullfile(inpath, ['SigChansPadded_P[', subj, ']', name_part, '.mat']);
            
            if exist(sigfile, 'file')
                Ssig = load(sigfile, 'sigChansBoth');
                if isfield(Ssig, 'sigChansBoth') && ~isempty(Ssig.sigChansBoth)
                    is_sig_here = ismember(chan_labels, Ssig.sigChansBoth);
                else
                    warning('Sig-channel file %s exists but sigChansBoth is missing/empty - marking ALL channels as NOT significant for %s/%s/%dHz/%s.', ...
                        sigfile, subj, cond, f0, sess);
                end
            else
                warning('Sig-channel file not found: %s - marking ALL channels as NOT significant for %s/%s/%dHz/%s.', ...
                    sigfile, subj, cond, f0, sess);
            end
        end

        % --- Append rows to the long-format table ---
        for c = 1:numel(chan_labels)
            rowi = rowi + 1;
            rows(rowi, :) = { string(subj), string(cond_factor), f0, ...
                              string(sess), string(chan_labels{c}), ...
                              plv_per_chan(c), total_trials, ...
                              is_sig_here(c) }; %#ok<SAGROW>
        end

        fprintf('  %s / %-6s / %d Hz / %s : %d trials [%s], mean PLV = %.3f\n', ...
            subj, cond, f0, sess, total_trials, ...
            strjoin(arrayfun(@(p,n) sprintf('%s:%d', cfg.polarity_labels{p}, n), ...
                             1:numel(cfg.polarity_labels), n_trials_per_pol, ...
                             'UniformOutput', false), ', '), ...
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
fprintf('\nPLV table saved to %s\n', cfg.plv_matfile);
fprintf('  %d rows, %d unique subjects, %d unique channels.\n', ...
    height(T_plv), numel(unique(T_plv.Subject)), numel(unique(T_plv.Channel)));

if cfg.write_xlsx
    writetable(T_plv, cfg.xlsx_file, 'Sheet', cfg.plv_xlsxsheet, ...
        'WriteMode', 'overwritesheet');
    fprintf('PLV table also written to %s (sheet: %s)\n', cfg.xlsx_file, cfg.plv_xlsxsheet);
end


%% Missing file summary

fprintf('\n--- Missing-file summary (out of %d subjects) ---\n', n_sub);
for r = 1:n_combos
    fprintf('  %-6s %3d Hz %-4s : %d/%d subjects missing\n', ...
        combos_template{r,1}, combos_template{r,2}, combos_template{r,3}, ...
        missing_count(r), n_sub);
end


%% Sanity check summary
% Print median PLV per condition/frequency so you can eyeball whether the
% computation produced sensible values before going to the LMM.

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
% Inverse of findgroups: get the unique label triples in group order.
[~, idx] = unique(G, 'stable');
condL = T.Condition(idx);
freqL = T.Frequency(idx);
sessL = T.Session(idx);
end