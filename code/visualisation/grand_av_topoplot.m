%% Grand-averaged topoplot
%
% Requires: preprocessed and averaged FFR data (the 'final' structure
% saved at the end of the preprocessing script)
% FieldTrip must be on the path

clc; clear; close all;

%% define parameters and shit
subjects   = {'NG', 'BVD', 'FD', 'EF', 'CX', 'HX', 'MCBG', 'LV', ...
              'AM', 'FS', 'LVB', 'RNR', 'VUR', 'LH', 'TF', 'MGZ', ...
              'NSMT', 'MM', 'TH', 'TD'};
conditions  = {'sit', 'firm', 'foam'};
frequencies = [89, 333];

basepath           = 'F:\Auditory-mov-data';
subfolder          = 'preprocessed';
trial_fname_pattern = 'averagedFFR_P[%s]_FFR_%s_%d.mat';  % subj, cond, freq
templ_file          = 'easycapM11.mat';
fz_label            = 'Fz';

n_sub  = numel(subjects);
n_cond = numel(conditions);
n_freq = numel(frequencies);

%% 2a. First pass: find union of all channel labels
all_labels = {};

for si = 1:n_sub
    subj   = subjects{si};
    inpath = fullfile(basepath, subj, subfolder);
    fname  = fullfile(inpath, sprintf(trial_fname_pattern, subj, conditions{1}, frequencies(1)));
    S = load(fname);
    all_labels = union(all_labels, S.final.label, 'stable');
end

n_chan = numel(all_labels);
fprintf('Total unique channels across subjects: %d\n', n_chan);

%% 2. Load data and compute per-subject peak/trough maps
% pos_maps{fi, ci} and neg_maps{fi, ci} are [nChan x n_sub] matrices

pos_maps = cell(n_freq, n_cond);
neg_maps = cell(n_freq, n_cond);

for fi = 1:n_freq
    f0 = frequencies(fi);
    for ci = 1:n_cond
        cond = conditions{ci};

        pos_stack = [];   % will grow to [nChan x n_valid_subjects]
        neg_stack = [];

        for si = 1:n_sub
            subj   = subjects{si};
            inpath = fullfile(basepath, subj, subfolder);
            fname  = fullfile(inpath, sprintf(trial_fname_pattern, subj, cond, f0));

            if ~isfile(fname)
                warning('Missing: %s', fname);
                continue
            end

            S     = load(fname);
            final = S.final;

            fz_idx = find(strcmp(final.label, fz_label));
            if isempty(fz_idx)
                warning('Fz not found for %s %s %dHz', subj, cond, f0);
                continue
            end

            waveform = final.avg(fz_idx, :);

            min_prom      = 0.7 * max(abs(waveform));
            [~, pos_locs] = findpeaks( waveform, 'MinPeakProminence', min_prom);
            [~, neg_locs] = findpeaks(-waveform, 'MinPeakProminence', min_prom);

            if isempty(pos_locs) || isempty(neg_locs)
                warning('No peaks found for %s %s %dHz', subj, cond, f0);
                continue
            end

            % Average at peak/trough samples (in subject's own channel space)
            pos_map_subj = mean(final.avg(:, pos_locs), 2);   % [subj_nChan x 1]
            neg_map_subj = mean(final.avg(:, neg_locs), 2);

           % Map into full channel space, NaN for missing channels
            pos_map = NaN(n_chan, 1);
            neg_map = NaN(n_chan, 1);
            [~, idx] = ismember(final.label, all_labels);
            valid   = idx > 0;                        % guard against unmatched labels
            pos_map(idx(valid)) = pos_map_subj(valid);
            neg_map(idx(valid)) = neg_map_subj(valid);

            pos_stack = [pos_stack, pos_map]; %#ok<AGROW>
            neg_stack = [neg_stack, neg_map]; %#ok<AGROW>
        end

        pos_maps{fi, ci} = pos_stack;
        neg_maps{fi, ci} = neg_stack;
    end
end

%% 3. Grand average: mean across subjects -> build ft timelock structs
% Load one file just to get the label list and a template struct
S_templ     = load(fullfile(basepath, subjects{1}, subfolder, ...
    sprintf(trial_fname_pattern, subjects{1}, conditions{1}, frequencies(1))));
templ       = S_templ.final;
templ.label = all_labels;

grandAvg_pos = cell(n_freq, n_cond);
grandAvg_neg = cell(n_freq, n_cond);

for fi = 1:n_freq
    for ci = 1:n_cond
        if isempty(pos_maps{fi, ci})
            continue
        end

        % Mean across subjects (columns)
        ga_pos_map = nanmean(pos_maps{fi, ci}, 2);  % [nChan x 1]
        ga_neg_map = nanmean(neg_maps{fi, ci}, 2);

        % Build minimal ft timelock structs (same approach as single-subject script)
        topo_pos          = templ;
        topo_pos.avg      = ga_pos_map;
        topo_pos.time     = 0;
        topo_pos.dimord   = 'chan_time';

        topo_neg          = templ;
        topo_neg.avg      = ga_neg_map;
        topo_neg.time     = 0;
        topo_neg.dimord   = 'chan_time';

        grandAvg_pos{fi, ci} = topo_pos;
        grandAvg_neg{fi, ci} = topo_neg;
    end
end


%% 4. Plot: one figure per frequency, rows = pos/neg, columns = conditions
for fi = 1:n_freq
    f0 = frequencies(fi);

    % Shared colour scale across all conditions for this frequency
    all_vals = [];
    for ci = 1:n_cond
        if ~isempty(grandAvg_pos{fi, ci})
            all_vals = [all_vals; grandAvg_pos{fi, ci}.avg; ...
                                  grandAvg_neg{fi, ci}.avg]; %#ok<AGROW>
        end
    end
    clim_val = max(abs(all_vals(~isnan(all_vals))));   % ignore NaNs in colour scale
    clim     = [-clim_val, clim_val];

    figure('Color', 'w', 'Position', [100 100 400*n_cond 750], ...
           'Name', sprintf('Grand avg topoplot — %d Hz', f0));

    for ci = 1:n_cond
        if isempty(grandAvg_pos{fi, ci})
            continue
        end

        cfg          = [];
        cfg.layout   = templ_file;
        cfg.marker   = 'on';
        cfg.colorbar = 'yes';
        cfg.zlim     = clim;
        cfg.colormap = 'RdBu';
        %cfg.colormap = 'plasma';
        cfg.style   = 'both';          % contour lines + colour fill
        %cfg.style    = 'straight';      % colormap only
        cfg.comment  = 'no';

        % Positive phase
        subplot(2, n_cond, ci);
        cfg.figure = gca;
        ft_topoplotER(cfg, grandAvg_pos{fi, ci});
        title(sprintf('%s — positive', conditions{ci}), 'FontSize', 11);

        % Negative phase
        subplot(2, n_cond, n_cond + ci);
        cfg.figure = gca;
        ft_topoplotER(cfg, grandAvg_neg{fi, ci});
        title(sprintf('%s — negative', conditions{ci}), 'FontSize', 11);
    end

    sgtitle(sprintf('Grand Average FFR Topography — %d Hz  (N=%d)', f0, n_sub), ...
        'FontSize', 13, 'FontWeight', 'bold');
end