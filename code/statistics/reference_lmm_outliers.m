%% LMM: Multi-subject 3 conditions x 2 references on FFR SNR & Power
clc
clear
close all

%% Set parameters
subjectids = {'NG', 'BVD', 'FD', 'EF', 'CX', 'HX', 'MCBG', 'LV', 'NG', 'BVD', 'AM', 'FS', 'LVB', 'RNR', 'VUR', 'LH', 'TF', 'MGZ', 'NSMT', 'MM', 'TH', 'TD'}
basepath_root = 'D:\Auditory-mov-data';

conditions = {'sit', 'firm', 'foam', 'sit_2'};
condFiles  = {'_FFR_sit_333', '_FFR_firm_333', '_FFR_foam_333', '_FFR_sit_333_2'};
refs       = {'mastoid', 'CAR'};
refFolders = {'preprocessed', 'preprocessed_CAR'};

nC = numel(conditions);
nR = numel(refs);
nS = numel(subjectids);

% --- Outlier handling (Cook's distance, two-pass refit) ---
remove_outliers = true;
outlier_thresh  = 'auto';   % 'auto' = 4/n, or a numeric value (e.g. 0.01)

% --- Outcome variable for the SNR LMM (matches main LMM convention) ---
%   'SNR_dB'  : log-transformed SNR (recommended; more normal residuals)
%   'SNR_lin' : raw linear SNR ratio
snr_outcome = 'SNR_lin';

%% FFT / SNR parameters (from your existing script)
fs   = 10000;
NFFT = ceil(fs / 2);
df   = fs / NFFT;
freq_vec = (0 : floor(NFFT/2)) * df;
nFreqs   = numel(freq_vec);

% Signal & noise bands (333 Hz F0)
sig_idx      = freq_vec >= 331 & freq_vec <= 335;
noise_lo_idx = freq_vec >= 321 & freq_vec <= 325;
noise_hi_idx = freq_vec >= 341 & freq_vec <= 345;
noise_idx    = noise_lo_idx | noise_hi_idx;

fprintf('Frequency resolution: %.2f Hz\n', df);
fprintf('Subjects: %s\n\n', strjoin(subjectids, ', '));

%% Loop: load + compute PSD/SNR/Power per subject, condition, reference
rows = {};

for s = 1:nS
    sid = subjectids{s};
    fprintf('=== Subject %s ===\n', sid);
    basepath = fullfile(basepath_root, sid);

    psd_subj    = cell(nC, nR);
    labels_subj = cell(nC, nR);
    sig_subj    = cell(nC, nR);

    for c = 1:nC
        for r = 1:nR
            folder  = fullfile(basepath, refFolders{r});
            avgFile = ['averagedFFR_P[', sid, ']', condFiles{c}, '.mat'];
            sigFile1 = ['SigChansPadded_P[', sid, ']', condFiles{c}, '.mat'];
            sigFile2 = ['SigchansPadded_P[', sid, ']', condFiles{c}, '.mat'];

            avgPath = fullfile(folder, avgFile);
            if      isfile(fullfile(folder, sigFile1)),  sigPath = fullfile(folder, sigFile1);
            elseif  isfile(fullfile(folder, sigFile2)),  sigPath = fullfile(folder, sigFile2);
            else,                                        sigPath = '';
            end

            if ~isfile(avgPath) || isempty(sigPath)
                fprintf('  [SKIP] %s / %s -- missing file\n', conditions{c}, refs{r});
                continue
            end

            S    = load(avgPath, 'final');         ffr  = S.final;
            S    = load(sigPath, 'sigChansBoth');  sigC = S.sigChansBoth;

            avg    = ffr.avg;
            labels = ffr.label;
            nChan  = size(avg, 1);
            N      = size(avg, 2);
            win    = hann(N)';
            win_pw = sum(win .^ 2);

            psd_cond = zeros(nChan, nFreqs);
            for ch = 1:nChan
                Sft = fft(avg(ch,:) .* win, NFFT);
                Sft = Sft(1 : floor(NFFT/2) + 1);
                psd_full = (abs(Sft).^2) / (win_pw * fs);
                psd_full(2:end-1) = 2 * psd_full(2:end-1);
                psd_cond(ch, 1:numel(psd_full)) = psd_full;
            end

            psd_subj{c,r}    = psd_cond;
            labels_subj{c,r} = labels;
            sig_subj{c,r}    = sigC;

            fprintf('  Loaded %s / %-7s: %d chans, %d sig\n', ...
                conditions{c}, refs{r}, nChan, numel(sigC));
        end
    end

    for c = 1:nC
        if isempty(labels_subj{c,1}) || isempty(labels_subj{c,2})
            continue
        end

        common_chans = intersect(labels_subj{c,1}, labels_subj{c,2}, 'stable');
        % Always exclude reference + trigger channels (same as main LMM script)
        common_chans = setdiff(common_chans, {'TP9','TP10','Trigger'}, 'stable');

        sig_keep = intersect(sig_subj{c,1}, sig_subj{c,2}, 'stable');
        sig_keep = intersect(common_chans, sig_keep, 'stable');

        if isempty(sig_keep)
            fprintf('  [WARN] %s / %s: no channels significant in both refs -- using all common\n', sid, conditions{c});
            sig_keep = common_chans;
        end

        for r = 1:nR
            [~, idx] = ismember(sig_keep, labels_subj{c,r});
            psd_cell = psd_subj{c,r}(idx, :);

            % SNR computation identical to main LMM 'custom' method:
            %   signal = mean PSD over 331-335 Hz (= F0 +/- 2)
            %   noise  = mean PSD over 321-325 and 341-345 Hz (F0 +/- [8,12])
            %   snr_lin = signal / max(noise, eps)  -- guard against zero
            %   snr_db  = 10*log10(snr_lin)
            sig_power   = mean(psd_cell(:, sig_idx),   2);
            noise_power = mean(psd_cell(:, noise_idx), 2);
            snr_lin     = sig_power ./ max(noise_power, eps);
            snr_db      = 10 * log10(snr_lin);

            for k = 1:numel(sig_keep)
                rows(end+1, :) = { sid, conditions{c}, refs{r}, sig_keep{k}, ...
                                   snr_lin(k), snr_db(k), sig_power(k) }; %#ok<SAGROW>
            end
        end

        fprintf('  %s: %d paired channels (mast<->CAR) kept\n', conditions{c}, numel(sig_keep));
    end
    fprintf('\n');
end

if isempty(rows)
    error('No data loaded -- check paths and filenames.');
end

%% Build long-format table
T = cell2table(rows, 'VariableNames', ...
    {'Subject','Condition','Reference','Channel','SNR_lin','SNR_dB','Power'});
% Outcome columns: SNR_dB (recommended) and SNR_lin are computed directly.
% Power gets log-transformed for the second LMM as before.
T.logPower = log(T.Power);
T.Subject   = categorical(T.Subject);
T.Channel   = categorical(T.Channel);
T.Condition = categorical(T.Condition, conditions);
T.Reference = categorical(T.Reference, refs);
T.Condition = reordercats(T.Condition, conditions);
T.Reference = reordercats(T.Reference, {'mastoid','CAR'});

fprintf('Long-format table: %d observations from %d subjects\n', height(T), numel(unique(T.Subject)));
disp(grpstats(T, {'Subject','Condition','Reference'}, 'numel', 'DataVars', snr_outcome))

%% Fit LMMs -- two-pass with optional Cook's-distance outlier removal
re_struct = '+ (1|Subject) + (1|Subject:Channel)';

fprintf('\n========== LMM 1: %s ==========\n', snr_outcome);
[lmm_snr, T_snr_trim, snr_outliers] = fit_with_outliers( ...
    T, [snr_outcome ' ~ Condition * Reference ' re_struct], ...
    remove_outliers, outlier_thresh, snr_outcome);

fprintf('\n========== LMM 2: log(Power at 333 Hz) ==========\n');
[lmm_pow, T_pow_trim, pow_outliers] = fit_with_outliers( ...
    T, ['logPower ~ Condition * Reference ' re_struct], ...
    remove_outliers, outlier_thresh, 'Power');

%% Likelihood-ratio tests for the interaction (run on the trimmed tables)
lmm_snr_main = fitlme(T_snr_trim, [snr_outcome ' ~ Condition + Reference ' re_struct]);
lmm_pow_main = fitlme(T_pow_trim, ['logPower ~ Condition + Reference '   re_struct]);

fprintf('\n=== LRT: Condition x Reference interaction (SNR) ===\n');
disp(compare(lmm_snr_main, lmm_snr))
fprintf('\n=== LRT: Condition x Reference interaction (Power) ===\n');
disp(compare(lmm_pow_main, lmm_pow))

%% Residual diagnostics for both models
figure('Position',[100 100 900 600],'Name','LMM residual diagnostics');
subplot(2,2,1); histogram(residuals(lmm_snr), 30); title('SNR residuals');   xlabel('residual')
subplot(2,2,2); qqplot(residuals(lmm_snr));        title('SNR Q-Q')
subplot(2,2,3); histogram(residuals(lmm_pow), 30); title('Power residuals'); xlabel('residual')
subplot(2,2,4); qqplot(residuals(lmm_pow));        title('Power Q-Q')

%% Post-hoc contrasts on the trimmed models
print_contrasts(lmm_snr, T_snr_trim, conditions, refs, snr_outcome);
print_contrasts(lmm_pow, T_pow_trim, conditions, refs, 'Power');

%% Visualisations
plot_cond_x_ref(T_snr_trim, snr_outcome, conditions, refs, lmm_snr, ['FFR ' snr_outcome ' [331-335 Hz]']);
plot_cond_x_ref(T_pow_trim, 'Power',     conditions, refs, lmm_pow, 'FFR signal power at 331-335 Hz');

%% =================== HELPER FUNCTIONS ===================

function [lmm, T_used, info] = fit_with_outliers(T, formula, do_remove, thr, label)
% Fit LMM. If do_remove is true, perform Cook's-distance based outlier
% detection on the first fit and refit on the trimmed table. Returns the
% (possibly refit) model, the table actually used, and an info struct.

% --- First (full-data) fit ---
lmm_full = fitlme(T, formula);
fprintf('--- %s, full data: n=%d obs ---\n', label, height(T));
disp(lmm_full);
disp(anova(lmm_full, 'DFMethod','satterthwaite'))

info = struct('method','none','threshold',NaN,'flagged_rows',[], ...
              'n_flagged',0,'flagged_table',table());

if ~do_remove
    lmm = lmm_full; T_used = T; return
end

% --- Cook's distance ---
res_raw = residuals(lmm_full, 'Conditional', true);
res_std = res_raw ./ std(res_raw, 'omitnan');
Xfx = designMatrix(lmm_full, 'Fixed');
try
    H = Xfx * ((Xfx' * Xfx) \ Xfx');
    h = diag(H);
catch
    h = sum(Xfx.^2, 2) ./ max(sum(Xfx(:).^2), eps);
end
p     = size(Xfx, 2);
cooks = (h ./ max((1 - h).^2, eps)) .* (res_std.^2 ./ max(p, 1));

n = height(T);
if ischar(thr) && strcmpi(thr,'auto'), thr = 4/n; end
keep = cooks < thr;

info.method        = 'cooks';
info.threshold     = thr;
info.scores        = cooks;          % full per-observation Cook's distances
info.flagged_rows  = find(~keep);
info.n_flagged     = sum(~keep);
if info.n_flagged > 0
    Tf = T(~keep, {'Subject','Condition','Reference','Channel'});
    Tf.CooksDistance = cooks(~keep);
    info.flagged_table = Tf;
end

fprintf('  Cook''s distance (threshold = %.4g): flagged %d / %d obs (%.1f%%)\n', ...
        thr, info.n_flagged, n, 100*info.n_flagged/n);

% --- Diagnostic plot ---
plot_cooks(cooks, thr, info.flagged_rows, label);

if info.n_flagged == 0
    lmm = lmm_full; T_used = T; return
end

% --- Refit on trimmed table ---
T_used = T(keep, :);
% Drop unused categorical levels so the design matrix stays full rank
for v = T_used.Properties.VariableNames
    if iscategorical(T_used.(v{1}))
        T_used.(v{1}) = removecats(T_used.(v{1}));
    end
end

try
    lmm = fitlme(T_used, formula);
    fprintf('--- %s, outliers removed: n=%d obs ---\n', label, height(T_used));
    disp(lmm);
    disp(anova(lmm, 'DFMethod','satterthwaite'))
catch ME
    warning('Refit after outlier removal failed (%s); using full-data fit.', ME.message);
    lmm = lmm_full; T_used = T;
end
end


function print_contrasts(lmm, T, conditions, refs, label)
    fprintf('\n=== Post-hoc contrasts (%s, log-scale) ===\n', label);

    fprintf('\nReference effect (mastoid - CAR) within each condition:\n');
    for c = 1:numel(conditions)
        H = build_contrast(lmm, T, 'Condition', conditions{c}, 'Reference', 'mastoid', 'CAR');
        [p, F, df1, df2] = coefTest(lmm, H);
        est = H * fixedEffects(lmm);
        fprintf('  %-6s:  est = %+.3f log-units (%.2fx),  F(%d, %.1f) = %.2f,  p = %.4f\n', ...
            conditions{c}, est, exp(est), df1, df2, F, p);
    end

    % All unique pairwise condition contrasts (e.g. for 4 conditions: 6 pairs)
    condPairs = {};
    for a = 1:numel(conditions)
        for b = a+1:numel(conditions)
            condPairs(end+1,:) = {conditions{a}, conditions{b}};  %#ok<AGROW>
        end
    end
    nPairs = size(condPairs, 1);

    fprintf('\nPairwise condition contrasts within each reference (Bonferroni-corrected across %d pairs):\n', nPairs);
    for r = 1:numel(refs)
        fprintf('  Reference = %s:\n', refs{r});
        p_raw = zeros(nPairs,1); est = zeros(nPairs,1);
        for k = 1:nPairs
            H = build_contrast(lmm, T, 'Reference', refs{r}, 'Condition', condPairs{k,1}, condPairs{k,2});
            [p_raw(k), ~, ~, ~] = coefTest(lmm, H);
            est(k) = H * fixedEffects(lmm);
        end
        p_corr = min(p_raw * nPairs, 1);
        for k = 1:nPairs
            fprintf('    %-6s vs %-6s:  est = %+.3f,  p_raw = %.4f,  p_Bonf = %.4f\n', ...
                condPairs{k,1}, condPairs{k,2}, est(k), p_raw(k), p_corr(k));
        end
    end
end


function plot_cond_x_ref(T, varname, conditions, refs, lmm, ttl)
    nC = numel(conditions);
    nR = numel(refs);

    sm = grpstats(T, {'Subject','Condition','Reference'}, 'mean', 'DataVars', varname);

    subjects = categories(removecats(sm.Subject));
    nS = numel(subjects);
    mat = nan(nS, nC, nR);
    for s = 1:nS
        for c = 1:nC
            for r = 1:nR
                mask = sm.Subject == subjects{s} & ...
                       sm.Condition == conditions{c} & ...
                       sm.Reference == refs{r};
                if any(mask), mat(s,c,r) = sm.(['mean_' varname])(mask); end
            end
        end
    end

    nBoxes = nC * nR;
    use_channel_level = (nS < 2);

    if use_channel_level
        fprintf('  [%s plot] Only %d subject(s); using channel-level values for boxplot.\n', varname, nS);
        cell_vals = cell(1, nBoxes);
        col = 1;
        for c = 1:nC
            for r = 1:nR
                mask = T.Condition == conditions{c} & T.Reference == refs{r};
                cell_vals{col} = T.(varname)(mask);
                col = col + 1;
            end
        end
        maxLen = max(cellfun(@numel, cell_vals));
        bp_plot = nan(maxLen, nBoxes);
        for col = 1:nBoxes
            v = cell_vals{col};
            bp_plot(1:numel(v), col) = v;
        end
    else
        bp_plot = nan(nS, nBoxes);
        col = 1;
        for c = 1:nC
            for r = 1:nR
                bp_plot(:, col) = mat(:, c, r);
                col = col + 1;
            end
        end
    end

    % X positions: nR adjacent boxes per condition, with a gap between conditions
    within_step = 0.7;                       % gap between refs inside a condition
    between_gap = 1.5;                       % gap between conditions
    xpos = zeros(1, nBoxes);
    box_labels = cell(1, nBoxes);
    col = 1;
    for c = 1:nC
        x0 = (c-1)*(nR*within_step + between_gap) + 1;
        for r = 1:nR
            xpos(col)       = x0 + (r-1)*within_step;
            box_labels{col} = abbrevRef(refs{r});
            col = col + 1;
        end
    end

    figure('Position', [100 100 max(900, 250*nC) 550], 'Name', ttl);
    h = boxplot(bp_plot, 'Notch','off', 'Positions', xpos, ...
        'Labels', box_labels, 'Widths', 0.5);
    set(gca,'fontsize', 13); set(h, {'linew'}, {1});

    % Color per reference
    if nR == 2
        colors = [0.4860 0.6940 0.1680;
                  0.8500 0.3250 0.0980];
    else
        colors = lines(nR);
    end
    bh = findobj(gca,'Tag','Box');
    % boxplot returns boxes in reverse order
    for k = 1:numel(bh)
        ri = mod(numel(bh)-k, nR) + 1;
        patch(get(bh(k),'XData'), get(bh(k),'YData'), colors(ri,:), 'FaceAlpha', .26);
    end

    hold on
    % Connect mast<->CAR per condition (only meaningful when nR == 2)
    if nR == 2
        if use_channel_level
            for c = 1:nC
                idxL = (c-1)*nR + 1;
                xL = xpos(idxL); xR = xpos(idxL+1);
                mask_m = T.Condition == conditions{c} & T.Reference == refs{1};
                mask_c = T.Condition == conditions{c} & T.Reference == refs{2};
                chans_m = T.Channel(mask_m);  vals_m = T.(varname)(mask_m);
                chans_c = T.Channel(mask_c);  vals_c = T.(varname)(mask_c);
                [~, ia, ib] = intersect(chans_m, chans_c, 'stable');
                for k = 1:numel(ia)
                    plot([xL xR], [vals_m(ia(k)) vals_c(ib(k))], ...
                        'Color',[0.7 0.7 0.7], 'LineWidth', 0.5);
                end
            end
        else
            for c = 1:nC
                idxL = (c-1)*nR + 1;
                for s = 1:nS
                    plot([xpos(idxL) xpos(idxL+1)], [mat(s,c,1) mat(s,c,2)], ...
                        'Color',[0.6 0.6 0.6], 'Marker','.', 'MarkerSize', 12);
                end
            end
        end
    end

    yL = ylim;
    for c = 1:nC
        idxs = (c-1)*nR + (1:nR);
        text(mean(xpos(idxs)), yL(1) - 0.08*range(yL), ...
            upper(conditions{c}), 'HorizontalAlignment','center', ...
            'FontWeight','bold', 'FontSize', 13);
    end

    if use_channel_level
        title(sprintf('%s  (n = 1 subject, channel-level)', ttl))
    else
        title(sprintf('%s  (n = %d subjects)', ttl, nS))
    end
    ylabel(varname)
    xlim([min(xpos)-0.7 max(xpos)+0.7])

    yMax = max(bp_plot(:),[],'omitnan');
    yR   = range(bp_plot(:),'omitnan');
    if yR == 0 || isnan(yR), yR = 1; end
    yStep = 0.025 * yR;
    ylim([yL(1) - 0.12*range(yL), yMax + 4*yStep]);

    % Significance brackets for the within-condition reference effect
    % (only meaningful when nR == 2; mast vs CAR)
    if nR == 2
        for c = 1:nC
            H = build_contrast(lmm, T, 'Condition', conditions{c}, 'Reference', 'mastoid', 'CAR');
            p_c = coefTest(lmm, H);
            if p_c <= 0.05
                idxL = (c-1)*nR + 1;
                xL = xpos(idxL); xR = xpos(idxL+1); y = yMax + yStep;
                plot([xL xL xR xR], [y y+0.5*yStep y+0.5*yStep y], 'k-', 'LineWidth', 1);
                if      p_c <= 0.001, st = '***';
                elseif  p_c <= 0.01,  st = '**';
                else,                 st = '*';
                end
                text(mean([xL xR]), y + 0.7*yStep, st, 'HorizontalAlignment','center', 'FontSize', 13);
            end
        end
    end

    % Legend
    hLeg = gobjects(1, nR);
    for r = 1:nR
        hLeg(r) = patch(NaN, NaN, colors(r,:), 'FaceAlpha', .26);
    end
    legend(hLeg, refs, 'Location', 'northeast')
end


function s = abbrevRef(refname)
% Short label for reference name in boxplot ticks
switch lower(refname)
    case 'mastoid', s = 'Mast';
    case 'car',     s = 'CAR';
    otherwise,      s = refname;
end
end


function H = build_contrast(lmm, T, fixVar, fixLevel, contrastVar, levelA, levelB)
    rowA = mean_design_row(lmm, T, fixVar, fixLevel, contrastVar, levelA);
    rowB = mean_design_row(lmm, T, fixVar, fixLevel, contrastVar, levelB);
    H = reshape(rowA - rowB, 1, []);
end


function row = mean_design_row(lmm, T, varA, levA, varB, levB)
    mask = (T.(varA) == levA) & (T.(varB) == levB);
    X    = designMatrix(lmm, 'Fixed');
    row  = mean(X(mask,:), 1);
end


function plot_cooks(cooks, thr, flagged, label)
% Two-panel diagnostic figure for Cook's distance:
%   left  = stem plot per observation with threshold line + flagged in red
%   right = log-scale histogram of distances
figure('Position', [120 120 1000 450], 'Color','w', ...
       'Name', sprintf('Cook''s distance (%s)', label));
tiledlayout(1, 2, 'Padding','compact', 'TileSpacing','compact');

% Left panel
nexttile;
hold on; grid on; box on
plot(1:numel(cooks), cooks, '.', 'MarkerSize', 7, 'Color', [0.4 0.4 0.4]);
yline(thr, 'r--', sprintf('threshold = %.4g', thr), ...
      'LineWidth', 1.2, 'LabelHorizontalAlignment','left');
if ~isempty(flagged)
    plot(flagged, cooks(flagged), 'ro', 'MarkerSize', 7, 'LineWidth', 1.2);
end
xlabel('Observation index')
ylabel('Cook''s distance')
title(sprintf('%s  |  %d / %d flagged', label, numel(flagged), numel(cooks)))

% Right panel: log-scale histogram
nexttile;
histogram(cooks, 60, 'FaceColor',[0.4 0.6 0.85], 'EdgeColor','w');
hold on
xline(thr, 'r--', 'LineWidth', 1.2);
xlabel('Cook''s distance')
ylabel('Count')
title('Distribution (log-y to show tails)')
set(gca, 'YScale', 'log')
grid on; box on
end