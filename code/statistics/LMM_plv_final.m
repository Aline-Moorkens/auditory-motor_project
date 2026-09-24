%% Linear Mixed Models for FFR phase-locking value (PLV)
% Outcome:  PLV at F0, precomputed by compute_PLV.m and saved as
%           <outpath>/PLV_table.mat (a long-format table with
%           Subject, Condition, Frequency, Session, Channel, PLV,
%           NumTrials, IsSigHere).
%
% PLV is bounded in [0, 1]. Two outcome choices:
%   'PLV'       : raw PLV (simplest interpretation, but residuals can be
%                  skewed near the 0 and 1 boundaries)
%   'PLV_logit' : log(PLV / (1 - PLV)), unbounded and usually closer to
%                  normal-distributed residuals; preferred for LMMs.
%
% Comparisons (toggle individually in the CONFIG block, or run any %% section
% on its own with Ctrl+Enter):
%   (A) 333 Hz, surface conditions:  sit vs firm vs foam
%   (B) 333 Hz, pre vs post:         sit vs sit_2
%   (C) 89  Hz, surface conditions:  sit vs firm vs foam
%   (D) 89 vs 333 Hz across surfaces
%
% LMM (per comparison):
%   <outcome> ~ <fixed effects> + (1 | Subject) + (1 | Channel)
% Effect sizes reported:
%   - Fixed-effect betas with 95% CIs
%   - Marginal & conditional R^2
%   - Cohen's d (paired) for each pairwise contrast
%
% Aline / Auditory-mov-data pipeline -- PLV LMM module

clc; clear; close all

%% Configuration
cfg = struct();

% --- Paths ---
cfg.basepath   = 'F:\Auditory-mov-data';
cfg.outpath    = 'F:\Auditory-mov-data\LMM_PLV_robust';

% Where PLV.m saved the long-format PLV table.
cfg.plv_table_file = fullfile(cfg.basepath,'PLV_all','PLV_table.mat');

% Subjects (not used in this script)
cfg.subjects   = {'NG', 'BVD', 'FD', 'EF', 'CX', 'HX', 'MCBG', 'LV', 'AM', ...
                  'FS', 'LVB', 'RNR', 'VUR', 'LH', 'TF', 'MGZ', 'NSMT', ...
                  'MM', 'TH', 'TD'};

% Channel restriction
% If true, drop rows where IsSigHere==false using a combine rule:
%   'union'     : channel must be significant in AT LEAST ONE condition
%                 (default; keeps channels that drop below threshold in
%                 some conditions, so the LMM can quantify the actual drop)
%   'intersect' : channel must be significant in ALL conditions
%                 (more conservative; can hide the effect you want to see)
cfg.restrict_to_sigchans = true;
cfg.sig_combine          = 'intersect';

% Model structure
% Whether to include (1|Channel) as a random effect.
% Set false for small ROIs (e.g. vertex with ~8 electrodes) where channel
% variance is undefined or parked at the boundary.
cfg.channel_random_effect = true;

% Outcome variable
% 'PLV'       : raw PLV (in [0, 1])
% 'PLV_logit' : log(PLV / (1 - PLV)), unbounded — usually better for LMMs
%
% If you pick PLV_logit, the script derives that column from PLV
% automatically (with a small epsilon to handle exactly-0 / exactly-1).
cfg.outcome = 'PLV_logit';

cfg.cohensd = 'PLV';
% NOTE: Cohen's d is maybe best computed on raw PLV (not cfg.outcome), for
% interpretability -- regardless of which scale the LMM itself uses.

% Epsilon used for the logit transform to avoid -Inf/+Inf at PLV=0 or 1.
cfg.logit_eps = 1e-4;

% Outlier handling
% Two-pass fit: fit once, flag influential observations, refit without them.
%   method  : 'cooks'      -> Cook's distance (recommended)
%             'residual_z' -> |standardised residuals| > threshold
%   thresh  : numeric value, or 'auto' (cooks: 4/n, residual_z: 3.0)
cfg.remove_outliers = false;
cfg.outlier_method  = 'cooks';
cfg.outlier_thresh  = 'auto';

% Robust estimation (IRLS, Tukey's bisquare)
cfg.robust_estimation = true;
cfg.robust_tuning      = 4.685;
cfg.robust_maxiter     = 20;
cfg.robust_tol         = 1e-3;

% Which of the three fits (full / cooks / robust) feeds the headline
% contrasts, LOSO, assumption checks, and the main plotResults figure.
cfg.headline_fit = 'robust';      % 'full' | 'cooks' | 'robust'

%  Plot labelling
cfg.label_subjects      = true;
cfg.label_subject_mode  = 'all';   % 'all' or 'outliers'

% Which analyses to run
cfg.run.A_within_333  = true;   % sit vs firm vs foam at 333 Hz
cfg.run.B_prepost_333 = true;   % sit vs sit_2 at 333 Hz
cfg.run.C_within_89   = true;   % sit vs firm vs foam at 89 Hz
cfg.run.D_between_freq= true;   % 89 vs 333 across conditions

% Plotting
cfg.make_plots = true;
cfg.figpath    = fullfile(cfg.outpath, 'figures_PLV');
if cfg.make_plots && ~exist(cfg.figpath,'dir'), mkdir(cfg.figpath); end

% Assumption checks
cfg.checkAssumptions = true;

% Leave-one-subject-out sensitivity
cfg.run_loso = true;

% --- Output ---
if ~exist(cfg.outpath,'dir'), mkdir(cfg.outpath); end
cfg.logfile = fullfile(cfg.outpath, 'FFR_LMM_results_PLV.xlsx');


%% LOAD THE PLV TABLE 
% compute_PLV.m already did the per-trial PLV computation and saved a
% long-format table. Here we just load it, optionally derive the logit
% transform, and prepare it for the LMMs.

fprintf('\n=== Loading PLV table from %s ===\n', cfg.plv_table_file);
if ~exist(cfg.plv_table_file, 'file')
    error(['PLV table not found at %s. ', ...
           'Run compute_PLV.m first.'], cfg.plv_table_file);
end
S = load(cfg.plv_table_file, 'T_plv');
T = S.T_plv;

% Make sure key columns are the expected types
if ~iscategorical(T.Subject),   T.Subject   = categorical(T.Subject);   end
if ~iscategorical(T.Condition), T.Condition = categorical(T.Condition); end
if ~iscategorical(T.Session),   T.Session   = categorical(T.Session);   end
if ~iscategorical(T.Channel),   T.Channel   = categorical(T.Channel);   end
if ~islogical(T.IsSigHere),     T.IsSigHere = logical(T.IsSigHere);     end

% Derive the logit-transformed PLV if requested. We keep both columns in
% the table so plots can show the raw PLV scale even when the LMM uses
% the logit transform.
plv_clipped = max(min(T.PLV, 1 - cfg.logit_eps), cfg.logit_eps);
T.PLV_logit = log(plv_clipped ./ (1 - plv_clipped));

fprintf('Loaded %d rows, %d unique subjects, %d unique channels.\n', ...
        height(T), numel(unique(T.Subject)), numel(unique(T.Channel)));

% --- Per-subject/condition summary (sanity check) ---
fprintf('\n--- Median PLV per Subject x Condition x Frequency x Session ---\n');
G = findgroups(T.Subject, T.Condition, T.Frequency, T.Session);
sub_meds = splitapply(@(x) median(x,'omitnan'), T.PLV, G);
[subL, condL, freqL, sessL] = groupLabels(T, G);
% (Don't dump 60+ rows to console; just print the per-condition aggregates)

fprintf('\n--- Median PLV per Condition x Frequency x Session ---\n');
G2 = findgroups(T.Condition, T.Frequency, T.Session);
plv_meds = splitapply(@(x) median(x,'omitnan'), T.PLV, G2);
plv_iqr_lo = splitapply(@(x) quantile(x,0.25), T.PLV, G2);
plv_iqr_hi = splitapply(@(x) quantile(x,0.75), T.PLV, G2);
n_obs      = splitapply(@(x) numel(x), T.PLV, G2);
[condL2, freqL2, sessL2] = condFreqSessLabels(T, G2);
for k = 1:numel(plv_meds)
    fprintf('  %-6s %3d Hz %-4s : median PLV = %.3f, IQR = [%.3f, %.3f], n = %d\n', ...
        condL2(k), freqL2(k), sessL2(k), ...
        plv_meds(k), plv_iqr_lo(k), plv_iqr_hi(k), n_obs(k));
end


%% (A) 333 Hz: sit vs firm vs foam
if cfg.run.A_within_333
    fprintf('\n=== (A) 333 Hz: sit vs firm vs foam ===\n');
    Ta = T(T.Frequency==333 & T.Session=="pre" & ...
           ismember(T.Condition, {'sit','firm','foam'}), :);
    Ta = dropUnused(Ta, {'Condition','Subject','Channel'});
    if cfg.restrict_to_sigchans
        Ta.Condition = categorical(Ta.Condition);
        Ta = selectSigChans(Ta, 'Condition', cfg.sig_combine);
    end

    tag = 'A_333';
    Ta.Condition = reorderCategoricalRef(Ta.Condition, 'sit', {'sit','firm','foam'});

    if checkFactor(Ta, 'Condition', 2, tag)
        chanRE = '';
        if cfg.channel_random_effect, chanRE = ' + (1|Channel)'; end
        formula = [cfg.outcome ' ~ Condition + (1|Subject)' chanRE];
        res_A   = fitLMM(Ta, formula, tag, cfg);

        lme = res_A.lme;
        contrasts = struct( ...
               'Label', {'firm_vs_sit','foam_vs_sit','foam_vs_firm'}, ...
               'L', { makeContrast(lme, {'Condition_firm'}), ...
                      makeContrast(lme, {'Condition_foam'}), ...
                      makeContrast(lme, {'Condition_foam'}, {'Condition_firm'}) });

        pvals = nan(numel(contrasts),1);
        for k = 1:numel(contrasts)
            r = contrastTest(lme, contrasts(k).L);
            contrasts(k).Estimate = r.Estimate;
            contrasts(k).SE       = r.SE;
            contrasts(k).tStat    = r.tStat;
            contrasts(k).DF       = r.DF;
            contrasts(k).pValue   = r.pValue;
            contrasts(k).CI_lower = r.CI_lower;
            contrasts(k).CI_upper = r.CI_upper;
            pvals(k) = r.pValue;
        end
        padj = bonferroniAdjust(pvals);
        res_A.contrasts = buildContrastsTable(contrasts, 'pValue_bonf', padj);
        disp(res_A.contrasts);

        res_A.cohens_d = cohensD(Ta, 'Condition', cfg.cohensd);
        writeResults(res_A, tag, cfg);
        if cfg.run_loso
            leaveOneSubjectOut(res_A, tag, cfg);
        end
        if cfg.checkAssumptions
            checkAssumptions(res_A, Ta, 'Condition', tag, cfg);
        end
        if cfg.make_plots
            plotOutlierComparison(res_A, tag, cfg);
            plotRobustWeights(res_A, tag, cfg);
            plotResults(Ta, res_A, 'Condition', tag, '333 Hz: sit vs firm vs foam', cfg);
        end
    end
end



%% (B) 333 Hz: pre vs post (sit vs sit_2)
if cfg.run.B_prepost_333
    fprintf('\n=== (B) 333 Hz: sit (pre) vs sit_2 (post) ===\n');
    Tb = T(T.Frequency==333 & T.Condition=="sit", :);
    Tb.Session = categorical(Tb.Session, {'pre','post'});
    Tb = dropUnused(Tb, {'Session','Subject','Channel'});
    Tb = keepPaired(Tb, 'Session');
    if cfg.restrict_to_sigchans
        Tb = selectSigChans(Tb, 'Session', cfg.sig_combine);
    end

    tag = 'B_333_prepost';

    if checkFactor(Tb, 'Session', 2, tag)
        chanRE = '';
        if cfg.channel_random_effect, chanRE = ' + (1|Channel)'; end
        formula = [cfg.outcome ' ~ Session + (1|Subject)' chanRE];
        res_B   = fitLMM(Tb, formula, tag, cfg);

        lme = res_B.lme;
        L   = makeContrast(lme, {'Session_post'});
        r   = contrastTest(lme, L);
        contrasts = struct('Label', {'post_vs_pre'}, 'L', {L}, ...
            'Estimate', r.Estimate, 'SE', r.SE, 'tStat', r.tStat, 'DF', r.DF, ...
            'pValue', r.pValue, 'CI_lower', r.CI_lower, 'CI_upper', r.CI_upper);
        res_B.contrasts = buildContrastsTable(contrasts, 'pValue_bonf', contrasts.pValue);
        disp(res_B.contrasts);

        res_B.cohens_d = cohensD(Tb, 'Session', cfg.cohensd);
        writeResults(res_B, tag, cfg);
        if cfg.run_loso
            leaveOneSubjectOut(res_B, tag, cfg);
        end
        if cfg.checkAssumptions
            checkAssumptions(res_B, Tb, 'Session', tag, cfg);
        end
        if cfg.make_plots
            plotOutlierComparison(res_B, tag, cfg);
            plotRobustWeights(res_B, tag, cfg);
            plotResults(Tb, res_B, 'Session', tag, '333 Hz: sit-pre vs sit-post', cfg);
        end
    end
end



%% (C) 89 Hz: sit vs firm vs foam
if cfg.run.C_within_89
    fprintf('\n=== (C) 89 Hz: sit vs firm vs foam ===\n');
    Tc = T(T.Frequency==89 & T.Session=="pre" & ...
           ismember(T.Condition, {'sit','firm','foam'}), :);
    Tc = dropUnused(Tc, {'Condition','Subject','Channel'});
    if cfg.restrict_to_sigchans
        Tc.Condition = categorical(Tc.Condition);
        Tc = selectSigChans(Tc, 'Condition', cfg.sig_combine);
    end

    tag = 'C_89';
    Tc.Condition = reorderCategoricalRef(Tc.Condition, 'sit', {'sit','firm','foam'});

    if checkFactor(Tc, 'Condition', 2, tag)
        chanRE = '';
        if cfg.channel_random_effect, chanRE = ' + (1|Channel)'; end
        formula = [cfg.outcome ' ~ Condition + (1|Subject)' chanRE];
        res_C   = fitLMM(Tc, formula, tag, cfg);

        lme = res_C.lme;
        contrasts = struct( ...
               'Label', {'firm_vs_sit','foam_vs_sit','foam_vs_firm'}, ...
               'L', { makeContrast(lme, {'Condition_firm'}), ...
                      makeContrast(lme, {'Condition_foam'}), ...
                      makeContrast(lme, {'Condition_foam'}, {'Condition_firm'}) });

        pvals = nan(numel(contrasts),1);
        for k = 1:numel(contrasts)
            r = contrastTest(lme, contrasts(k).L);
            contrasts(k).Estimate = r.Estimate;
            contrasts(k).SE       = r.SE;
            contrasts(k).tStat    = r.tStat;
            contrasts(k).DF       = r.DF;
            contrasts(k).pValue   = r.pValue;
            contrasts(k).CI_lower = r.CI_lower;
            contrasts(k).CI_upper = r.CI_upper;
            pvals(k) = r.pValue;
        end
        padj = bonferroniAdjust(pvals);
        res_C.contrasts = buildContrastsTable(contrasts, 'pValue_bonf', padj);
        disp(res_C.contrasts);

        res_C.cohens_d = cohensD(Tc, 'Condition', cfg.cohensd);
        writeResults(res_C, tag, cfg);
        if cfg.run_loso
            leaveOneSubjectOut(res_C, tag, cfg);
        end
        if cfg.checkAssumptions
            checkAssumptions(res_C, Tc, 'Condition', tag, cfg);
        end
        if cfg.make_plots
            plotOutlierComparison(res_C, tag, cfg);
            plotRobustWeights(res_C, tag, cfg);
            plotResults(Tc, res_C, 'Condition', tag, '89 Hz: sit vs firm vs foam', cfg);
        end
    end
end



%% (D) 89 vs 333 across conditions
if cfg.run.D_between_freq
    fprintf('\n=== (D) Frequency: 89 vs 333 (across conditions) ===\n');
    Td = T(T.Session=="pre" & ismember(T.Condition, {'sit','firm','foam'}), :);
    Td.Frequency = categorical(Td.Frequency);
    Td = dropUnused(Td, {'Condition','Frequency','Subject','Channel'});
    if cfg.restrict_to_sigchans
        Td.Condition = categorical(Td.Condition);
        Td.CondFreq = categorical(strcat(string(Td.Condition), '_', string(Td.Frequency)));
        Td = selectSigChans(Td, 'CondFreq', cfg.sig_combine);
        Td.CondFreq = [];
    end

    tag = 'D_FxC';
    Td.Condition = reorderCategoricalRef(Td.Condition, 'sit', {'sit','firm','foam'});
    Td.Frequency = reorderCategoricalRef(Td.Frequency, '89', {'89','333'});

    if checkFactor(Td, 'Frequency', 2, tag) && checkFactor(Td, 'Condition', 2, tag)
        chanRE = '';
        if cfg.channel_random_effect, chanRE = ' + (1|Channel)'; end
        formula = [cfg.outcome ' ~ Frequency * Condition + (1|Subject)' chanRE];
        res_D   = fitLMM(Td, formula, tag, cfg);

        lme = res_D.lme;
        % IMPORTANT: confirm interaction coefficient name order via
        % disp(lme.CoefficientNames) once -- may come out as
        % 'Condition_firm:Frequency_333' rather than the reverse (this
        % was the exact issue we hit in the SNR script's Section D).
        contrasts = struct( ...
                    'Label', {'freq333v89_at_sit','freq333v89_at_firm','freq333v89_at_foam'}, ...
                    'L',{ makeContrast(lme, {'Frequency_333'}), ...
                          makeContrast(lme, {'Frequency_333','Condition_firm:Frequency_333'}), ...
                          makeContrast(lme, {'Frequency_333','Condition_foam:Frequency_333'}) });

        pvals = nan(numel(contrasts),1);
        for k = 1:numel(contrasts)
            r = contrastTest(lme, contrasts(k).L);
            contrasts(k).Estimate = r.Estimate;
            contrasts(k).SE       = r.SE;
            contrasts(k).tStat    = r.tStat;
            contrasts(k).DF       = r.DF;
            contrasts(k).pValue   = r.pValue;
            contrasts(k).CI_lower = r.CI_lower;
            contrasts(k).CI_upper = r.CI_upper;
            pvals(k) = r.pValue;
        end
        padj = fdrBH(pvals);
        res_D.contrasts = buildContrastsTable(contrasts, 'pValue_fdr', padj);
        disp(res_D.contrasts);

        res_D.cohens_d = cohensD(Td, 'Frequency', cfg.cohensd);
        writeResults(res_D, tag, cfg);
        if cfg.run_loso
            leaveOneSubjectOut(res_D, tag, cfg);
        end
        if cfg.checkAssumptions
            checkAssumptions(res_D, Td, 'Frequency', tag, cfg);
        end
        if cfg.make_plots
            plotOutlierComparison(res_D, tag, cfg);
            plotRobustWeights(res_D, tag, cfg);
            plotResults(Td, res_D, 'Frequency', tag, '89 vs 333 Hz across conditions', cfg);
        end
    end
end

fprintf('\nAll selected PLV analyses complete. Results in:\n  %s\n', cfg.logfile);


%% HELPER FUNCTIONS
function L = makeContrast(lme, posNames, negNames)
if nargin < 3, negNames = {}; end
cNames = lme.CoefficientNames;
L = zeros(1, numel(cNames));
for i = 1:numel(posNames)
    idx = strcmp(cNames, posNames{i});
    if ~any(idx)
        error('makeContrast: coefficient "%s" not found. Available: %s', ...
              posNames{i}, strjoin(cNames, ', '));
    end
    L(idx) = L(idx) + 1;
end
for i = 1:numel(negNames)
    idx = strcmp(cNames, negNames{i});
    if ~any(idx)
        error('makeContrast: coefficient "%s" not found. Available: %s', ...
              negNames{i}, strjoin(cNames, ', '));
    end
    L(idx) = L(idx) - 1;
end
end


function out = contrastTest(lme, L, alpha)
if nargin < 3, alpha = 0.05; end
beta = fixedEffects(lme);
covB = lme.CoefficientCovariance;
est = L * beta;
se  = sqrt(L * covB * L');
[pval, Fstat, df] = coefTest(lme, L, 0, 'DFMethod', 'satterthwaite');
tstat = sign(est) * sqrt(max(Fstat, 0));
tcrit = tinv(1 - alpha/2, df);
ci_lo = est - tcrit*se;
ci_hi = est + tcrit*se;
out = struct('Estimate', est, 'SE', se, 'tStat', tstat, 'DF', df, ...
             'pValue', pval, 'CI_lower', ci_lo, 'CI_upper', ci_hi);
end


function padj = bonferroniAdjust(pvals)
m = numel(pvals);
padj = min(pvals * m, 1);
end


function padj = fdrBH(pvals)
m = numel(pvals);
[psort, ord] = sort(pvals(:));
padj_sorted = psort .* m ./ (1:m)';
padj_sorted = flipud(cummin(flipud(padj_sorted)));
padj = nan(size(pvals));
padj(ord) = min(padj_sorted, 1);
end


function T = buildContrastsTable(contrasts, padjName, padjVals)
% Forces the L column to stay a cell array (see the struct2table bug we
% found in the SNR script -- same-length numeric fields get silently
% merged into a matrix column otherwise).
Label    = {contrasts.Label}';
L        = {contrasts.L}';
Estimate = [contrasts.Estimate]';
SE       = [contrasts.SE]';
tStat    = [contrasts.tStat]';
DF       = [contrasts.DF]';
pValue   = [contrasts.pValue]';
CI_lower = [contrasts.CI_lower]';
CI_upper = [contrasts.CI_upper]';
T = table(Label, L, Estimate, SE, tStat, DF, pValue, CI_lower, CI_upper, ...
    'VariableNames', {'Label','L','Estimate','SE','tStat','DF','pValue','CI_lower','CI_upper'});
T.(padjName) = padjVals(:);
end


function ok = checkFactor(T, varName, minLevels, tag)
ok = true;
if ~ismember(varName, T.Properties.VariableNames)
    fprintf('  [SKIP %s] factor "%s" not present in table.\n', tag, varName);
    ok = false; return
end
v = T.(varName);
if iscategorical(v), lv = categories(removecats(v));
else,                lv = unique(string(v));
end
nObs = height(T);
fprintf('  Factor "%s": %d level(s) (%s), %d observations\n', ...
        varName, numel(lv), strjoin(string(lv), ', '), nObs);
if numel(lv) < minLevels
    fprintf('  [SKIP %s] need >= %d levels of "%s", found %d.\n', ...
            tag, minLevels, varName, numel(lv));
    ok = false;
end
if nObs < 10
    fprintf('  [SKIP %s] too few observations (%d).\n', tag, nObs);
    ok = false;
end
end


function T = dropUnused(T, vars)
for i = 1:numel(vars)
    if ismember(vars{i}, T.Properties.VariableNames) && ...
       iscategorical(T.(vars{i}))
        T.(vars{i}) = removecats(T.(vars{i}));
    end
end
end


function T = keepPaired(T, factorVar)
if ~ismember(factorVar, T.Properties.VariableNames), return; end
lv = categories(removecats(T.(factorVar)));
subs = unique(T.Subject);
keep = false(numel(subs),1);
for i = 1:numel(subs)
    have = unique(string(T.(factorVar)(T.Subject==subs(i))));
    keep(i) = all(ismember(string(lv), have));
end
sub_keep = subs(keep);
sub_drop = subs(~keep);
if ~isempty(sub_drop)
    fprintf('  Excluding subjects without data in all "%s" levels: %s\n', ...
            factorVar, strjoin(string(sub_drop), ', '));
end
T = T(ismember(T.Subject, sub_keep), :);
end


function T = selectSigChans(T, factorVar, mode)
if ~ismember(factorVar, T.Properties.VariableNames) || isempty(T), return; end
if nargin < 3 || isempty(mode), mode = 'union'; end

if ~ismember('IsSigHere', T.Properties.VariableNames)
    warning('PLV:NoSigCol', ...
            'selectSigChans: IsSigHere column missing; skipping.');
    return
end

lv   = categories(removecats(T.(factorVar)));
subs = unique(T.Subject);
keep_rows = false(height(T), 1);

for i = 1:numel(subs)
    sub_rows = T.Subject == subs(i);
    if ~any(sub_rows), continue; end

    sub_chans = unique(string(T.Channel(sub_rows)));
    sig_count = zeros(numel(sub_chans), 1);
    for cidx = 1:numel(sub_chans)
        ch_rows = sub_rows & T.Channel == sub_chans(cidx);
        sig_lv  = unique(string(T.(factorVar)(ch_rows & T.IsSigHere)));
        sig_count(cidx) = numel(sig_lv);
    end

    switch lower(mode)
        case 'intersect'
            chan_ok = sub_chans(sig_count >= numel(lv));
        case 'union'
            chan_ok = sub_chans(sig_count >= 1);
        otherwise
            error('Unknown sig_combine mode: %s', mode);
    end

    keep_rows = keep_rows | (sub_rows & ismember(string(T.Channel), chan_ok));
end

fprintf('  Combine sig-channels across "%s" (mode=%s): %d -> %d rows\n', ...
        factorVar, mode, height(T), sum(keep_rows));
T = T(keep_rows, :);
T = dropUnused(T, {'Subject','Channel', factorVar});
end


function res = fitLMM(Tin, formula, tag, cfg)
res = struct();
res.tag     = tag;
res.formula = formula;

res.full = fitOneLMM(Tin, formula);
printLMM(res.full, 'full data');

if isfield(cfg,'remove_outliers') && cfg.remove_outliers
    [keepMask, flagInfo] = flagOutliers(res.full.lme, Tin, cfg);
    res.outlier_info = flagInfo;
    n_flag = sum(~keepMask);
    fprintf('  Cook''s-trim (%s): flagged %d / %d obs (%.1f%%)\n', ...
            flagInfo.method, n_flag, numel(keepMask), 100*n_flag/numel(keepMask));
    if n_flag == 0
        res.trimmed = res.full;
    else
        Ttrim = dropUnused(Tin(keepMask,:), Tin.Properties.VariableNames);
        try
            res.trimmed = fitOneLMM(Ttrim, formula);
            printLMM(res.trimmed, 'cooks-trim');
        catch ME
            warning('Cook''s-trim refit failed: %s', ME.message);
            res.trimmed = res.full;
        end
    end
else
    res.trimmed = res.full;
end

if isfield(cfg,'robust_estimation') && cfg.robust_estimation
    try
        res.robust = fitRobustLMM(Tin, formula, cfg);
        printLMM(res.robust, sprintf('robust (bisquare, %d iter)', res.robust.n_iter));
    catch ME
        warning('Robust fit failed: %s', ME.message);
        res.robust = res.full;
    end
else
    res.robust = res.full;
end

headline = 'robust';
if isfield(cfg,'headline_fit') && ~isempty(cfg.headline_fit)
    headline = lower(cfg.headline_fit);
end
switch headline
    case 'full',   main = res.full;
    case 'cooks',  main = res.trimmed;
    case 'robust', main = res.robust;
    otherwise
        warning('Unknown cfg.headline_fit "%s"; defaulting to robust.', headline);
        main = res.robust; headline = 'robust';
end
res.headline_name = headline;

res.lme            = main.lme;
res.fixedEffects   = main.fixedEffects;
res.R2_marginal    = main.R2_marginal;
res.R2_conditional = main.R2_conditional;
res.n_obs          = main.n_obs;
res.n_subj         = main.n_subj;
end


function out = fitOneLMM(T, formula)
out = struct();
out.n_obs   = height(T);
out.n_subj  = numel(unique(T.Subject));

lme = fitlme(T, formula, 'FitMethod','REML');
out.lme = lme;

[~,~,FE] = fixedEffects(lme, 'DFMethod','satterthwaite');
out.fixedEffects = FE;

[psi, mse] = covarianceParameters(lme);
out.randomEffects_cov = psi;
out.residual_var      = mse;

yhat_fixed = predict(lme, T, 'Conditional', false);
var_fixed  = var(yhat_fixed, 0, 'omitnan');
var_rand = 0;
for k = 1:numel(psi)
    M = psi{k};
    if ~isempty(M), var_rand = var_rand + M(1,1); end
end
var_res = mse;

out.R2_marginal    = var_fixed / (var_fixed + var_rand + var_res);
out.R2_conditional = (var_fixed + var_rand) / (var_fixed + var_rand + var_res);
end


function out = fitRobustLMM(T, formula, cfg)
maxIter = cfg.robust_maxiter;
tol     = cfg.robust_tol;
tuning  = cfg.robust_tuning;

n = height(T);
w = ones(n,1);
lme = [];

for iter = 1:maxIter
    lme_new = fitlme(T, formula, 'FitMethod','REML', 'Weights', w);
    res_raw = residuals(lme_new, 'Conditional', true);
    s = median(abs(res_raw - median(res_raw))) / 0.6745;
    s = max(s, eps);
    u = res_raw ./ (tuning * s);
    w_new = (1 - u.^2).^2;
    w_new(abs(u) >= 1) = 0;
    w_new = max(w_new, 1e-6);
    delta = max(abs(w_new - w));
    w   = w_new;
    lme = lme_new;
    if delta < tol, break; end
end

out = struct();
out.lme     = lme;
out.weights = w;
out.n_iter  = iter;
out.tuning  = tuning;
out.n_obs   = n;
out.n_subj  = numel(unique(T.Subject));

[~,~,FE] = fixedEffects(lme, 'DFMethod','satterthwaite');
out.fixedEffects = FE;

[psi, mse] = covarianceParameters(lme);
var_rand = 0;
for k = 1:numel(psi)
    M = psi{k};
    if ~isempty(M), var_rand = var_rand + M(1,1); end
end
yhat_fixed = predict(lme, T, 'Conditional', false);
var_fixed  = var(yhat_fixed, 0, 'omitnan');

out.R2_marginal    = var_fixed / (var_fixed + var_rand + mse);
out.R2_conditional = (var_fixed + var_rand) / (var_fixed + var_rand + mse);
end


function printLMM(out, label)
fprintf('--- LMM (%s): n=%d obs, %d subjects ---\n', label, out.n_obs, out.n_subj);
disp(out.lme);
fprintf('R^2 marginal    = %.3f\n', out.R2_marginal);
fprintf('R^2 conditional = %.3f\n', out.R2_conditional);
end


function [keep, info] = flagOutliers(lme, T, cfg)
info = struct();
n = height(T);

method = lower(cfg.outlier_method);
thr    = cfg.outlier_thresh;

res_raw = residuals(lme, 'Conditional', true);
res_std = res_raw ./ std(res_raw, 'omitnan');

switch method
    case 'cooks'
        Xfx = designMatrix(lme, 'Fixed');
        try
            H = Xfx * ((Xfx' * Xfx) \ Xfx');
            h = diag(H);
        catch
            h = sum(Xfx.^2, 2) ./ max(sum(Xfx(:).^2), eps);
        end
        p = size(Xfx, 2);
        cooks = (h ./ max((1 - h).^2, eps)) .* (res_std.^2 ./ max(p, 1));
        score = cooks;

        if ischar(thr) && strcmpi(thr,'auto'), thr = 4/n; end
        keep = score < thr;

        info.method     = 'cooks';
        info.threshold  = thr;
        info.score      = cooks;

    case 'residual_z'
        score = abs(res_std);
        if ischar(thr) && strcmpi(thr,'auto'), thr = 3.0; end
        keep = score < thr;

        info.method     = 'residual_z';
        info.threshold  = thr;
        info.score      = res_std;

    otherwise
        error('Unknown outlier method: %s', method);
end

info.flagged_rows = find(~keep);
if ~isempty(info.flagged_rows)
    Tf = T(~keep, {'Subject','Condition','Frequency','Session','Channel'});
    Tf.OutlierScore = info.score(~keep);
    info.flagged_table = Tf;
else
    info.flagged_table = table();
end
end


function d_tbl = cohensD(Tin, groupVar, depVar)
% Paired Cohen's d on subject-level means (averaged across channels).
G  = Tin.(groupVar);
lv = categories(G);
sub = unique(Tin.Subject);
M = nan(numel(sub), numel(lv));
for i = 1:numel(sub)
    for j = 1:numel(lv)
        sel = Tin.Subject==sub(i) & G==lv{j};
        if any(sel), M(i,j) = mean(Tin.(depVar)(sel), 'omitnan'); end
    end
end

have_mes = exist('meanEffectSize','file') == 2;

rows = {};
for a = 1:numel(lv)
    for b = a+1:numel(lv)
        valid = ~isnan(M(:,a)) & ~isnan(M(:,b));
        xa = M(valid, a);
        xb = M(valid, b);
        if numel(xa) < 2
            d = NaN; ci_lo = NaN; ci_hi = NaN;
        elseif have_mes
            try
                S = meanEffectSize(xa, xb, 'Paired', true, ...
                    'Effect','cohen', 'ConfidenceIntervalType','bootstrap');
                d     = S.Effect;
                ci_lo = S.ConfidenceIntervals(1);
                ci_hi = S.ConfidenceIntervals(2);
            catch
                diffs = xa - xb;
                d  = mean(diffs) / std(diffs);
                ci_lo = NaN; ci_hi = NaN;
            end
        else
            diffs = xa - xb;
            d  = mean(diffs) / std(diffs);
            ci_lo = NaN; ci_hi = NaN;
        end
        rows(end+1,:) = { string(lv{a}), string(lv{b}), ...
                          mean(xa,'omitnan'), mean(xb,'omitnan'), ...
                          d, ci_lo, ci_hi, numel(xa) }; %#ok<AGROW>
    end
end
d_tbl = cell2table(rows, 'VariableNames', ...
    {'Level_A','Level_B','Mean_A','Mean_B','Cohens_d_paired', ...
     'CI_lower','CI_upper','N_pairs'});
disp(d_tbl);
end


function writeResults(res, tag, cfg)
file = cfg.logfile;

headlineName = 'full';
if isfield(res,'headline_name') && ~isempty(res.headline_name)
    headlineName = res.headline_name;
end

FE = dataset2table(res.fixedEffects);
writetable(FE, file, 'Sheet', sheetName(tag, sprintf('FE_headline_%s', headlineName)), ...
           'WriteMode','overwritesheet');

hasTrim   = isfield(res,'trimmed') && isfield(res,'full') && ...
            ~isequal(res.full.n_obs, res.trimmed.n_obs);
hasRobust = isfield(res,'robust') && isfield(res.robust,'weights');

fitNames = {"full"};
n_obs    = res.full.n_obs;
n_subj   = res.full.n_subj;
R2m      = res.full.R2_marginal;
R2c      = res.full.R2_conditional;

if hasTrim
    fitNames{end+1} = "cooks-trim";
    n_obs(end+1)    = res.trimmed.n_obs;
    n_subj(end+1)   = res.trimmed.n_subj;
    R2m(end+1)      = res.trimmed.R2_marginal;
    R2c(end+1)      = res.trimmed.R2_conditional;
end
if hasRobust
    fitNames{end+1} = "robust";
    n_obs(end+1)    = res.robust.n_obs;
    n_subj(end+1)   = res.robust.n_subj;
    R2m(end+1)      = res.robust.R2_marginal;
    R2c(end+1)      = res.robust.R2_conditional;
end

isHeadline = strcmpi(fitNames, headlineName);
S = table(string(fitNames)', n_obs', n_subj', R2m', R2c', isHeadline', ...
    'VariableNames', {'Fit','N_obs','N_subjects','R2_marginal','R2_conditional','IsHeadline'});
writetable(S, file, 'Sheet', sheetName(tag,'summary'), 'WriteMode','overwritesheet');

if hasTrim
    FE_full = dataset2table(res.full.fixedEffects);
    writetable(FE_full, file, 'Sheet', sheetName(tag,'FE_full'), 'WriteMode','overwritesheet');
    if isfield(res, 'outlier_info') && ~isempty(res.outlier_info.flagged_table)
        writetable(res.outlier_info.flagged_table, file, ...
                   'Sheet', sheetName(tag,'outliers'), 'WriteMode','overwritesheet');
    end
end

if hasRobust
    FE_robust = dataset2table(res.robust.fixedEffects);
    writetable(FE_robust, file, 'Sheet', sheetName(tag,'FE_robust'), 'WriteMode','overwritesheet');
    Wtab = table((1:numel(res.robust.weights))', res.robust.weights, ...
        'VariableNames', {'ObsIndex','RobustWeight'});
    writetable(Wtab, file, 'Sheet', sheetName(tag,'robust_w'), 'WriteMode','overwritesheet');
end

if isfield(res,'cohens_d') && ~isempty(res.cohens_d)
    writetable(res.cohens_d, file, 'Sheet', sheetName(tag,'d'), 'WriteMode','overwritesheet');
end

if isfield(res,'contrasts') && ~isempty(res.contrasts)
    writetable(res.contrasts, file, 'Sheet', sheetName(tag,'contrasts'), 'WriteMode','overwritesheet');
end

fprintf('Wrote results for %s to %s (headline fit: %s)\n', tag, file, headlineName);
end


function name = sheetName(tag, suffix)
forbidden = {':','\','/','?','*','[',']'};
for k = 1:numel(forbidden)
    tag    = strrep(tag,    forbidden{k}, '-');
    suffix = strrep(suffix, forbidden{k}, '-');
end
maxLen = 31;
fullName = [tag '_' suffix];
if numel(fullName) > maxLen
    keepTag = maxLen - numel(suffix) - 1;
    if keepTag < 1, keepTag = 1; end
    tag = tag(1:min(numel(tag), keepTag));
    fullName = [tag '_' suffix];
end
name = fullName;
end


function plotResults(Tin, res, groupVar, tag, ttl, cfg)
G  = Tin.(groupVar);
lv = categories(removecats(G));
nL = numel(lv);

% Outcome variable name + axis label
outcomeVar = cfg.outcome;
switch outcomeVar
    case 'PLV',       outcomeLabel = 'PLV';
    case 'PLV_logit', outcomeLabel = 'logit(PLV)';
    otherwise,        outcomeLabel = outcomeVar;
end

sub = unique(Tin.Subject);
M = nan(numel(sub), nL);
for i = 1:numel(sub)
    for j = 1:nL
        sel = Tin.Subject==sub(i) & G==lv{j};
        if any(sel), M(i,j) = mean(Tin.(outcomeVar)(sel),'omitnan'); end
    end
end

fig = figure('Name', ['PLV LMM: ' tag], ...
             'Position', [80 80 1500 900], 'Color','w');
tiledlayout(fig, 2, 2, 'Padding','compact', 'TileSpacing','compact');

% =========================================================
% Panel 1: paired-subject plot (spaghetti) + group means
% =========================================================
nexttile;
hold on; grid on; box on
xpos = 1:nL;

do_label = isfield(cfg,'label_subjects') && cfg.label_subjects;
label_mode = 'all';
if isfield(cfg,'label_subject_mode') && ~isempty(cfg.label_subject_mode)
    label_mode = lower(cfg.label_subject_mode);
end
label_this = false(numel(sub),1);
if do_label
    switch label_mode
        case 'outliers'
            for j = 1:nL
                col = M(:,j);
                isout = isoutlier(col, 'quartiles');
                label_this = label_this | isout;
            end
        otherwise
            label_this = ~all(isnan(M),2);
    end
end

is_outlier = false(numel(sub),1);
for j = 1:nL
    col = M(:,j);
    is_outlier = is_outlier | isoutlier(col, 'quartiles');
end

for i = 1:numel(sub)
    if is_outlier(i)
        lc = [0.85 0.33 0.10];
        lw = 1.4;
    else
        lc = [0.6 0.6 0.6 0.6];
        lw = 0.7;
    end
    hL = plot(xpos, M(i,:), '-o', 'Color', lc, ...
        'MarkerSize', 5, 'MarkerFaceColor', [0.85 0.85 0.85], ...
        'LineWidth', lw);
    try
        hL.DataTipTemplate.DataTipRows(1) = ...
            dataTipTextRow('Subject', repmat({char(string(sub(i)))}, 1, nL));
        hL.DataTipTemplate.DataTipRows(2) = ...
            dataTipTextRow(groupVar, cellstr(string(lv(:)'))) ;
        hL.DataTipTemplate.DataTipRows(3) = ...
            dataTipTextRow(outcomeLabel, M(i,:));
    catch
    end
end

grp_mean = mean(M, 1, 'omitnan');
grp_sd   = std(M, 0, 1, 'omitnan');
nValid   = sum(~isnan(M), 1);
grp_se   = grp_sd ./ sqrt(max(nValid,1));
grp_ci   = 1.96 * grp_se;

errorbar(xpos, grp_mean, grp_ci, 'ko', 'LineWidth', 2.0, ...
    'MarkerSize', 10, 'MarkerFaceColor','k', 'CapSize', 12);

xticks(xpos); xticklabels(lv);
xlim([0.5 nL + 0.9])
xlabel(groupVar); ylabel(outcomeLabel)

if do_label && any(label_this)
    pts = struct('y_data',{},'y_lab',{},'id',{},'is_out',{});
    for i = 1:numel(sub)
        if ~label_this(i), continue; end
        lastCol = find(~isnan(M(i,:)), 1, 'last');
        if isempty(lastCol), continue; end
        pts(end+1).y_data = M(i, lastCol); %#ok<AGROW>
        pts(end).y_lab    = M(i, lastCol);
        pts(end).id       = char(string(sub(i)));
        pts(end).is_out   = is_outlier(i);
    end

    if ~isempty(pts)
        [~, ord] = sort([pts.y_data]);
        pts = pts(ord);
        yl  = ylim;
        min_spacing = (yl(2) - yl(1)) * 0.032;
        for k = 2:numel(pts)
            if pts(k).y_lab < pts(k-1).y_lab + min_spacing
                pts(k).y_lab = pts(k-1).y_lab + min_spacing;
            end
        end
        if pts(end).y_lab > yl(2)
            pts(end).y_lab = yl(2);
            for k = numel(pts)-1:-1:1
                if pts(k).y_lab > pts(k+1).y_lab - min_spacing
                    pts(k).y_lab = pts(k+1).y_lab - min_spacing;
                end
            end
        end

        x_data_edge = nL;
        x_leader    = nL + 0.18;
        x_text      = nL + 0.25;
        for k = 1:numel(pts)
            if abs(pts(k).y_lab - pts(k).y_data) > min_spacing * 0.2
                plot([x_data_edge, x_leader, x_text - 0.02], ...
                     [pts(k).y_data, pts(k).y_data, pts(k).y_lab], ...
                     '-', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.3, ...
                     'HandleVisibility', 'off');
            end
            if pts(k).is_out
                tcol = [0.65 0.20 0.05]; tw = 'bold';
            else
                tcol = [0.20 0.20 0.45]; tw = 'normal';
            end
            text(x_text, pts(k).y_lab, pts(k).id, ...
                'Color', tcol, 'FontSize', 6.5, 'FontWeight', tw, ...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
        end
    end
end

if do_label && strcmp(label_mode,'outliers')
    title('Subject-level means (orange = outliers, black = mean \pm 95% CI)')
else
    title('Subject means (grey/orange = subjects [orange = outlier], black = mean \pm 95% CI)')
end

% =========================================================
% Panel 2: per-channel mean PLV by condition
% =========================================================
nexttile;
chans = categories(removecats(Tin.Channel));
H = nan(numel(chans), nL);
for ic = 1:numel(chans)
    for j = 1:nL
        sel = Tin.Channel==chans{ic} & G==lv{j};
        if any(sel), H(ic,j) = mean(Tin.(outcomeVar)(sel),'omitnan'); end
    end
end
imagesc(H); colorbar
xticks(1:nL); xticklabels(lv); xlabel(groupVar)
set(gca,'YDir','normal','FontSize',7)
title(sprintf('Mean %s per channel \\times condition', outcomeLabel))

if numel(chans) <= 32
    yt = 1:numel(chans);
    yticks(yt); yticklabels(chans);
else
    yt = unique(round(linspace(1, numel(chans), 16)));
    yticks(yt); yticklabels(chans(yt));
    set(gca,'FontSize',8)
end

% =========================================================
% Panel 3: forest plot of fixed-effect estimates
% =========================================================
nexttile;
hold on; grid on; box on
FE = res.fixedEffects;
isInt = strcmp(string(FE.Name), '(Intercept)');
FEp = FE(~isInt, :);
nFE = height(FEp);

ypos = 1:nFE;
FEp  = flipud(FEp);

errorbar(FEp.Estimate, ypos, FEp.Estimate - FEp.Lower, FEp.Upper - FEp.Estimate, ...
    'horizontal', 'ko', 'LineWidth', 1.6, 'MarkerSize', 8, 'MarkerFaceColor','k', 'CapSize', 10);
xline(0, 'r--', 'LineWidth', 1.2);
if nFE >= 1
    yticks(ypos); yticklabels(cellstr(FEp.Name));
    ylim([0.5 nFE+0.5])
end
xlabel(sprintf('Estimate (%s), 95%% CI', outcomeLabel))
title(sprintf('Fixed effects (relative to ref. level)\nR^2_{marg}=%.3f, R^2_{cond}=%.3f', ...
              res.R2_marginal, res.R2_conditional))

for i = 1:nFE
    pv = FEp.pValue(i);
    if pv < 0.001, sig = '***';
    elseif pv < 0.01, sig = '**';
    elseif pv < 0.05, sig = '*';
    else, sig = 'n.s.';
    end
    text(FEp.Upper(i)+0.05, ypos(i), sprintf(' %s (p=%.3g)', sig, pv), ...
         'FontSize', 9, 'VerticalAlignment','middle');
end
xlim('padded')

% =========================================================
% Panel 4: variance decomposition
% =========================================================
nexttile;
[psi, mse] = covarianceParameters(res.lme);
var_sub = 0; var_chn = 0;
for k = 1:numel(psi)
    M2 = psi{k};
    if ~isempty(M2)
        if k==1, var_sub = M2(1,1);
        elseif k==2, var_chn = M2(1,1);
        end
    end
end
var_res = mse;
yhat_fixed = predict(res.lme, Tin, 'Conditional', false);
var_fix    = var(yhat_fixed, 0, 'omitnan');

components = [var_fix, var_sub, var_chn, var_res];
total = sum(components);
pct   = 100 * components ./ total;
labels = {'Fixed effects','Subject (random)','Channel (random)','Residual'};
b = bar(pct, 'FaceColor','flat');
b.CData = [0.20 0.45 0.85;
           0.85 0.55 0.20;
           0.45 0.75 0.45;
           0.55 0.55 0.55];
xticks(1:4); xticklabels(labels);
xtickangle(20)
ylabel('% of total variance')
title('Where does the variance live?')
grid on; box on
for i = 1:4
    text(i, pct(i)+1, sprintf('%.1f%%', pct(i)), ...
        'HorizontalAlignment','center', 'FontSize',10, 'FontWeight','bold')
end
ylim([0 max(pct)*1.2 + 5])

sgtitle(sprintf('%s   |   n=%d obs, %d subjects, %d channels', ...
        ttl, res.n_obs, res.n_subj, numel(chans)), ...
        'FontWeight','bold', 'FontSize',13)

outpng = fullfile(cfg.figpath, [tag '.png']);
exportgraphics(fig, outpng, 'Resolution', 200);
savefig(fig, fullfile(cfg.figpath, [tag '.fig']));
fprintf('Saved figure: %s\n', outpng);
end


function checkAssumptions(res, Tin, groupVar, tag, cfg)
lme   = res.lme;
res_raw = residuals(lme, 'Conditional', true);
res_std = res_raw ./ std(res_raw, 'omitnan');
yfit    = fitted(lme, 'Conditional', true);

outcomeVar = cfg.outcome;
switch outcomeVar
    case 'PLV',       unitLabel = 'PLV';
    case 'PLV_logit', unitLabel = 'logit';
    otherwise,        unitLabel = outcomeVar;
end

Tmodel = lme.Variables;

try
    [h_lil, p_lil] = lillietest(res_raw, 'Alpha', 0.05);
catch
    h_lil = NaN; p_lil = NaN;
end

if ismember(groupVar, Tmodel.Properties.VariableNames)
    G = Tmodel.(groupVar);
else
    G = Tin.(groupVar);
end
try
    p_lev = vartestn(res_raw, G, 'TestType','LeveneAbsolute', 'Display','off');
    F_lev = NaN;
catch
    p_lev = NaN; F_lev = NaN;
end

try
    Xfx = designMatrix(lme, 'Fixed');
    if all(Xfx(:,1) == 1), Xfx = Xfx(:, 2:end); end
    nc = size(Xfx, 2);
    if nc >= 2
        vif_vals = nan(nc, 1);
        for j = 1:nc
            mdl = fitlm(Xfx(:, setdiff(1:nc, j)), Xfx(:, j));
            vif_vals(j) = 1 / max(1 - mdl.Rsquared.Ordinary, eps);
        end
    else
        vif_vals = NaN;
    end
catch
    vif_vals = NaN;
end

fprintf('--- Assumption checks for %s ---\n', tag);
if isnan(p_lil),     msg_lil = '(n/a)';
elseif p_lil < 0.05, msg_lil = '(NORMALITY VIOLATED)';
else,                msg_lil = '(OK)';
end
fprintf('  Residual normality (Lilliefors):  p = %.4g  %s\n', p_lil, msg_lil);

if isnan(p_lev),     msg_lev = '(n/a)';
elseif p_lev < 0.05, msg_lev = '(EQUAL VARIANCE VIOLATED)';
else,                msg_lev = '(OK)';
end
fprintf('  Homoscedasticity (Levene by %s): p = %.4g  %s\n', groupVar, p_lev, msg_lev);

if all(isnan(vif_vals))
    fprintf('  VIF: not applicable (only one fixed predictor)\n');
else
    if any(vif_vals > 5), msg_vif = '(MULTICOLLINEARITY: VIF > 5)';
    else,                 msg_vif = '(OK)';
    end
    fprintf('  VIF values: %s  %s\n', sprintf('%.2f ', vif_vals), msg_vif);
end

fig = figure('Name', ['PLV LMM diagnostics: ' tag], ...
             'Position', [60 60 1500 900], 'Color','w');
tiledlayout(fig, 2, 3, 'Padding','compact', 'TileSpacing','compact');

nexttile;
qqplot(res_raw);
title(sprintf('Q-Q residuals\nLilliefors p=%.3g', p_lil))
grid on; box on
xlabel('Theoretical quantiles'); ylabel('Sample quantiles')

nexttile;
histogram(res_raw, 40, 'Normalization','pdf', ...
    'FaceColor',[0.4 0.6 0.85], 'EdgeColor','w'); hold on
xs = linspace(min(res_raw), max(res_raw), 200);
ys = normpdf(xs, mean(res_raw,'omitnan'), std(res_raw,'omitnan'));
plot(xs, ys, 'r-', 'LineWidth', 1.5);
title('Residual distribution (red = normal)')
xlabel(sprintf('Residual (%s)', unitLabel)); ylabel('Density')
grid on; box on

nexttile;
plot(yfit, res_raw, '.', 'MarkerSize', 6, 'Color',[0.4 0.4 0.4]); hold on
yline(0, 'r--', 'LineWidth', 1);
[fs_sort, ord] = sort(yfit);
rs_sort        = res_raw(ord);
w              = max(20, round(length(fs_sort)/30));
smooth         = movmean(rs_sort, w, 'omitnan');
plot(fs_sort, smooth, 'b-', 'LineWidth', 1.5);
xlabel('Fitted values'); ylabel('Residuals')
title(sprintf('Residuals vs fitted\nLevene p=%.3g (by %s)', p_lev, groupVar))
grid on; box on

nexttile;
plot(yfit, sqrt(abs(res_std)), '.', 'MarkerSize', 6, 'Color',[0.4 0.4 0.4]); hold on
smooth2 = movmean(sqrt(abs(res_std(ord))), w, 'omitnan');
plot(fs_sort, smooth2, 'b-', 'LineWidth', 1.5);
xlabel('Fitted values'); ylabel('sqrt(|standardised residuals|)')
title('Scale-location (homoscedasticity)')
grid on; box on

nexttile;
try
    Bsub = randomEffects(lme);
    if ismember('Subject', Tmodel.Properties.VariableNames)
        nSub = numel(unique(Tmodel.Subject));
    else
        nSub = numel(unique(Tin.Subject));
    end
    if numel(Bsub) >= nSub
        sub_blup = Bsub(1:nSub);
        qqplot(sub_blup);
        title('Q-Q of subject random intercepts')
    else
        text(0.1,0.5,'Random-effect BLUPs unavailable','Units','normalized');
        axis off
    end
catch
    text(0.1,0.5,'Random-effect BLUPs unavailable','Units','normalized');
    axis off
end
grid on; box on

nexttile;
if numel(G) == numel(res_raw)
    boxplot(res_raw, G, 'Symbol','.');
    yline(0, 'r--', 'LineWidth', 1);
    ylabel(sprintf('Residual (%s)', unitLabel)); xlabel(groupVar)
    title('Residuals by group (look for unequal spread)')
else
    text(0.1,0.5, sprintf('Length mismatch: %d residuals vs %d group labels', ...
        numel(res_raw), numel(G)), 'Units','normalized');
    axis off
end
grid on; box on

verdict = '';
if ~isnan(p_lil) && p_lil < 0.05, verdict = [verdict ' [normality !]']; end
if ~isnan(p_lev) && p_lev < 0.05, verdict = [verdict ' [variance !]']; end
if all(isnan(vif_vals))
    vifstr = 'VIF: n/a';
elseif any(vif_vals > 5)
    vifstr = sprintf('VIF max=%.1f', max(vif_vals));
    verdict = [verdict ' [VIF !]'];
else
    vifstr = sprintf('VIF max=%.1f', max(vif_vals));
end
sgtitle(sprintf('PLV LMM assumption checks: %s   |   %s%s', tag, vifstr, verdict), ...
        'FontWeight','bold', 'FontSize',13)

outpng = fullfile(cfg.figpath, [tag '_diagnostics.png']);
exportgraphics(fig, outpng, 'Resolution', 200);
fprintf('Saved diagnostics figure: %s\n', outpng);

if isfield(res, 'outlier_info') && ~isempty(res.outlier_info)
    plotOutlierScores(res, tag, cfg);
end

if ~isnan(vif_vals(1))
    vif_str = strjoin(arrayfun(@(x) sprintf('%.2f',x), vif_vals, 'uni', 0), ', ');
else
    vif_str = 'n/a';
end
A = table( ...
    string(tag), p_lil, h_lil, p_lev, F_lev, string(vif_str), ...
    'VariableNames', {'Tag','Lilliefors_p','Lilliefors_reject', ...
                      'Levene_p','Levene_F','VIF_values'});
writetable(A, cfg.logfile, 'Sheet', sheetName(tag,'assumptions'), 'WriteMode','overwritesheet');
end


function [subL, condL, freqL, sessL] = groupLabels(T, G)
nG = max(G);
subL  = strings(nG,1);
condL = strings(nG,1);
freqL = nan(nG,1);
sessL = strings(nG,1);
for k = 1:nG
    sel = find(G == k, 1, 'first');
    subL(k)  = string(T.Subject(sel));
    condL(k) = string(T.Condition(sel));
    freqL(k) = double(string(T.Frequency(sel)));
    sessL(k) = string(T.Session(sel));
end
end


function [condL, freqL, sessL] = condFreqSessLabels(T, G)
% Reduced version of groupLabels for the (Condition x Frequency x Session) grouping.
nG = max(G);
condL = strings(nG,1);
freqL = nan(nG,1);
sessL = strings(nG,1);
for k = 1:nG
    sel = find(G == k, 1, 'first');
    condL(k) = string(T.Condition(sel));
    freqL(k) = double(string(T.Frequency(sel)));
    sessL(k) = string(T.Session(sel));
end
end


function plotOutlierScores(res, tag, cfg)
info = res.outlier_info;
if isempty(info), return; end

score = info.score;
thr   = info.threshold;
flag  = info.flagged_rows;

fig = figure('Name', ['Outliers: ' tag], ...
             'Position', [120 120 1100 500], 'Color','w');
tiledlayout(fig, 1, 2, 'Padding','compact', 'TileSpacing','compact');

nexttile;
hold on; grid on; box on
plot(1:numel(score), score, '.', 'MarkerSize', 6, 'Color', [0.4 0.4 0.4]);
yline(thr, 'r--', sprintf('threshold = %.4g', thr), ...
      'LineWidth', 1.2, 'LabelHorizontalAlignment','left');
if ~isempty(flag)
    plot(flag, score(flag), 'ro', 'MarkerSize', 6, 'LineWidth', 1.2);
end
xlabel('Observation index'); ylabel(['Score (' info.method ')'])
title(sprintf('Outlier scores  |  %d / %d flagged', numel(flag), numel(score)))

nexttile;
histogram(score, 60, 'FaceColor',[0.4 0.6 0.85], 'EdgeColor','w');
hold on
xline(thr, 'r--', 'LineWidth', 1.2);
xlabel(['Score (' info.method ')']); ylabel('Count')
title('Distribution of outlier scores')
grid on; box on
set(gca, 'YScale', 'log')

sgtitle(sprintf('Outlier diagnostics: %s', tag), 'FontWeight','bold', 'FontSize',12)

outpng = fullfile(cfg.figpath, [tag '_outliers.png']);
exportgraphics(fig, outpng, 'Resolution', 200);
fprintf('Saved outlier figure: %s\n', outpng);
end


function plotOutlierComparison(res, tag, cfg)
% Comparison of fitting approaches (full / [cooks-trim] / robust):
%   1) Forest plot: each stored contrast, re-evaluated under every fit
%   2) Marginal & conditional R^2 by fit
%   3) Robust (bisquare) weight per observation, flagging w < 0.5
%      (only drawn if a robust fit is present)
%
% Which fits are shown adapts to cfg.remove_outliers: when false,
% the redundant "cooks-trim == full" series is omitted entirely.

if ~isfield(res,'contrasts') || isempty(res.contrasts), return; end
if ~(isfield(cfg,'make_plots') && cfg.make_plots), return; end

% ---------- Decide which fits to show ----------
fits   = {};
fitRes = {};
allColors = struct( ...
    'full',       [0.00 0.4470 0.7410], ...
    'cookstrim',  [0.8500 0.3250 0.0980], ...
    'robust',     [0.4660 0.6740 0.1880]);

fits{end+1} = 'full';        fitRes{end+1} = res.full;

if isfield(cfg,'remove_outliers') && cfg.remove_outliers && isfield(res,'trimmed')
    fits{end+1} = 'cooks-trim';  fitRes{end+1} = res.trimmed;
end

hasRobust = isfield(cfg,'robust_estimation') && cfg.robust_estimation && ...
            isfield(res,'robust') && isfield(res.robust,'weights');
if hasRobust
    fits{end+1} = 'robust';      fitRes{end+1} = res.robust;
end

nFits  = numel(fits);
colors = zeros(nFits,3);
for f = 1:nFits
    key = matlab.lang.makeValidName(strrep(fits{f},'-',''));
    if isfield(allColors, key)
        colors(f,:) = allColors.(key);
    else
        colors(f,:) = lines(1);   % fallback, shouldn't normally trigger
    end
end

% ---------- Evaluate every contrast under every fit ----------
C  = res.contrasts;
nC = height(C);
est   = nan(nC, nFits);
ci_lo = nan(nC, nFits);
ci_hi = nan(nC, nFits);
for f = 1:nFits
    for c = 1:nC
        try
            r = contrastTest(fitRes{f}.lme, C.L{c});
            est(c,f)   = r.Estimate;
            ci_lo(c,f) = r.CI_lower;
            ci_hi(c,f) = r.CI_upper;
        catch
        end
    end
end

% ---------- Figure layout: 3 panels if robust present, else 2 ----------
nPanels = 2 + double(hasRobust);
fig = figure('Name', ['Outlier comparison: ' tag], ...
             'Position', [60 60 min(1700, 600*nPanels) 550], 'Color','w');
tiledlayout(fig, 1, nPanels, 'Padding','compact', 'TileSpacing','compact');

% =========================================================
% Panel 1: forest plot, one row of markers per fit, offset vertically
% =========================================================
nexttile; hold on; grid on; box on

% Symmetric offsets around each contrast's y-position, spread across nFits
if nFits == 1
    offsets = 0;
else
    offsets = linspace(-0.18, 0.18, nFits);
end
ypos_base = nC:-1:1;   % top-to-bottom matches table order

for f = 1:nFits
    errorbar(est(:,f), ypos_base + offsets(f), ...
        est(:,f) - ci_lo(:,f), ci_hi(:,f) - est(:,f), ...
        'horizontal', 'o', 'Color', colors(f,:), 'MarkerFaceColor', colors(f,:), ...
        'MarkerSize', 7, 'LineWidth', 1.6, 'CapSize', 6, 'DisplayName', fits{f});
end
xline(0, 'k--', 'LineWidth', 1);
yticks(1:nC); yticklabels(flipud(cellstr(string(C.Label))));
ylim([0.5 nC+0.5]); xlabel('Estimate (95% CI)');
legend('Location','best');
title('Fixed effects by fit');

% =========================================================
% Panel 2: variance explained, grouped bars per fit
% =========================================================
nexttile;
R2m = nan(1,nFits); R2c = nan(1,nFits);
for f = 1:nFits
    R2m(f) = fitRes{f}.R2_marginal;
    R2c(f) = fitRes{f}.R2_conditional;
end
b = bar([R2m; R2c]');
b(1).FaceColor = [0.20 0.45 0.85];
b(2).FaceColor = [0.45 0.75 0.45];
xticks(1:nFits); xticklabels(fits);
ylabel('R^2'); ylim([0 1]); grid on; box on
legend({'marginal','conditional'}, 'Location','northoutside', 'Orientation','horizontal');
title('Variance explained');
for f = 1:nFits
    text(f-0.15, R2m(f)+0.02, sprintf('%.2f',R2m(f)), 'HorizontalAlignment','center','FontSize',9);
    text(f+0.15, R2c(f)+0.02, sprintf('%.2f',R2c(f)), 'HorizontalAlignment','center','FontSize',9);
end

% =========================================================
% Panel 3: robust weights -- only if a robust fit exists
% =========================================================
if hasRobust
    nexttile; hold on; grid on; box on
    w = res.robust.weights;
    n = numel(w);
    robustColor = allColors.robust;
    stem(1:n, w, 'Marker','none', 'Color', robustColor, 'LineWidth', 0.5);
    yline(0.5, 'r--', 'LineWidth', 1.2);
    lowIdx = find(w < 0.5);
    plot(lowIdx, w(lowIdx), 'ro', 'MarkerSize', 5, 'LineWidth', 1.2);
    xlabel('Observation index'); ylabel('Robust weight');
    ylim([0 1.05]); xlim([0 n+1]);
    title(sprintf('Robust weights (bisquare): %d / %d obs w<0.5', numel(lowIdx), n));
end

headlineName = 'full';
if isfield(res,'headline_name'), headlineName = res.headline_name; end
sgtitle(sprintf('Outlier-handling comparison: %s   |   headline fit: %s', ...
        tag, headlineName), 'FontWeight','bold', 'FontSize',13);

outpng = fullfile(cfg.figpath, [tag '_fitcompare.png']);
exportgraphics(fig, outpng, 'Resolution', 200);
fprintf('Saved outlier-comparison figure: %s\n', outpng);
end


function plotRobustWeights(res, tag, cfg)
% Stand-alone robust-weight diagnostics figure:
%   1) Per-observation bisquare weight (stem plot), threshold line at 0.5,
%      flagged points (w < 0.5) circled in red
%   2) Histogram of the weight distribution across all observations
%
% Mirrors plotOutlierScores (used for Cook's-distance diagnostics) but
% for the robust IRLS fit's final bisquare weights.

if ~isfield(res,'robust') || ~isfield(res.robust,'weights')
    fprintf('  [robust weights %s] no robust fit/weights found; skipping.\n', tag);
    return
end
if ~(isfield(cfg,'make_plots') && cfg.make_plots), return; end

w = res.robust.weights;
n = numel(w);
lowIdx = find(w < 0.5);

fig = figure('Name', ['Robust weights: ' tag], ...
             'Position', [100 100 1500 500], 'Color','w');
tiledlayout(fig, 1, 2, 'Padding','compact', 'TileSpacing','compact');

% =========================================================
% Panel 1: per-observation weight stem plot
% =========================================================
nexttile;
hold on; grid on; box on
stem(1:n, w, 'Marker','none', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 0.5);
yline(0.5, 'r--', 'LineWidth', 1.2);
if ~isempty(lowIdx)
    plot(lowIdx, w(lowIdx), 'ro', 'MarkerSize', 5, 'LineWidth', 1.2);
end
xlabel('Observation index'); ylabel('Robust weight');
ylim([0 1.05]); xlim([0 n+1]);
title(sprintf('bisquare weights: %d / %d obs downweighted (w<0.5)', numel(lowIdx), n));

% =========================================================
% Panel 2: weight distribution (histogram)
% =========================================================
nexttile;
histogram(w, 40, 'FaceColor', [0.6 0.85 0.6], 'EdgeColor', [0.4 0.6 0.4]);
xlabel('Robust weight'); ylabel('Count');
xlim([0 1]);
grid on; box on
title('Weight distribution');

sgtitle(sprintf('Robust estimation diagnostics: %s', tag), ...
        'FontWeight','bold', 'FontSize',13);

outpng = fullfile(cfg.figpath, [tag '_robustw.png']);
exportgraphics(fig, outpng, 'Resolution', 200);
fprintf('Saved robust-weight diagnostics figure: %s\n', outpng);
end


function loso_tbl = leaveOneSubjectOut(res, tag, cfg)
% LOSO sensitivity check, now consistent with whichever fit is the
% headline (full / cooks-trim / robust). Each per-subject refit uses
% the SAME fitting procedure as res.headline_name, so LOSO is actually
% stress-testing the estimator being reported, not a different one.

loso_tbl = table();

try
    Tin     = res.lme.Variables;
    formula = res.formula;
catch
    fprintf('  [LOSO %s] could not extract data from res.lme; skipping.\n', tag);
    return
end

if ~ismember('Subject', Tin.Properties.VariableNames)
    fprintf('  [LOSO %s] no Subject column in model data; skipping.\n', tag);
    return
end
if ~isfield(res, 'contrasts') || isempty(res.contrasts)
    fprintf('  [LOSO %s] no contrasts table found; skipping.\n', tag);
    return
end

headline = 'full';
if isfield(res,'headline_name') && ~isempty(res.headline_name)
    headline = res.headline_name;
end
fprintf('  [LOSO %s] refitting with headline method: %s\n', tag, headline);

subs  = unique(Tin.Subject);
n_sub = numel(subs);
if n_sub < 3
    fprintf('  [LOSO %s] only %d subject(s); LOSO is not meaningful.\n', tag, n_sub);
    return
end

C            = res.contrasts;
coefNames    = string(C.Label);
nCoef        = height(C);
refCoefNames = res.lme.CoefficientNames;

ref_est = C.Estimate;
ref_p   = C.pValue;

fprintf('--- Leave-one-subject-out for %s (n=%d subjects, %d contrasts, method=%s) ---\n', ...
        tag, n_sub, nCoef, headline);

est_mat   = nan(n_sub, nCoef);
p_mat     = nan(n_sub, nCoef);
ci_lo_mat = nan(n_sub, nCoef);
ci_hi_mat = nan(n_sub, nCoef);
fit_ok    = false(n_sub, 1);

for i = 1:n_sub
    Tsub = Tin(Tin.Subject ~= subs(i), :);
    Tsub = dropUnused(Tsub, Tsub.Properties.VariableNames);
    try
        switch headline
            case 'robust'
                out_i = fitRobustLMM(Tsub, formula, cfg);
                lme_i = out_i.lme;
            case 'cooks'
                % Re-run the same exclude-then-refit logic used headline-side.
                lme0 = fitlme(Tsub, formula, 'FitMethod','REML');
                [keepMask_i, ~] = flagOutliers(lme0, Tsub, cfg);
                if all(keepMask_i)
                    lme_i = lme0;
                else
                    Ttrim_i = dropUnused(Tsub(keepMask_i,:), Tsub.Properties.VariableNames);
                    lme_i = fitlme(Ttrim_i, formula, 'FitMethod','REML');
                end
            otherwise  % 'full'
                lme_i = fitlme(Tsub, formula, 'FitMethod','REML');
        end

        % Coefficient-space safety check -- L vectors are positional.
        if ~isequal(lme_i.CoefficientNames, refCoefNames)
            warning('LOSO %s: coefficient names changed when excluding %s -- skipping this subject.', ...
                    tag, char(string(subs(i))));
            continue
        end

        for c = 1:nCoef
            r = contrastTest(lme_i, C.L{c});
            est_mat(i,c)   = r.Estimate;
            p_mat(i,c)     = r.pValue;
            ci_lo_mat(i,c) = r.CI_lower;
            ci_hi_mat(i,c) = r.CI_upper;
        end
        fit_ok(i) = true;
    catch ME
        fprintf('  LOSO refit failed excluding %s (%s): %s\n', ...
                char(string(subs(i))), headline, ME.message);
    end
end

% --- Console summary: range of estimates/p-values, flag "flippers" ---
for c = 1:nCoef
    e  = est_mat(fit_ok, c);
    pv = p_mat(fit_ok, c);
    if isempty(e), continue; end
    fprintf('  %-22s : beta [%+.4g, %+.4g] (ref %+.4g),  p [%.3g, %.3g] (ref %.3g)\n', ...
        coefNames(c), min(e), max(e), ref_est(c), min(pv), max(pv), ref_p(c));

    if ref_p(c) < 0.05 && any(pv >= 0.05)
        flippers = subs(fit_ok & p_mat(:,c) >= 0.05);
        fprintf('    *** Significance LOST when excluding: %s\n', strjoin(string(flippers), ', '));
    elseif ref_p(c) >= 0.05 && any(pv < 0.05)
        flippers = subs(fit_ok & p_mat(:,c) < 0.05);
        fprintf('    Significance gained when excluding: %s\n', strjoin(string(flippers), ', '));
    end
end

% --- Wide-format table for Excel ---
T_wide = table(subs, fit_ok, 'VariableNames', {'Excluded_Subject','Fit_OK'});
for c = 1:nCoef
    safe = matlab.lang.makeValidName(char(coefNames(c)));
    T_wide.(['Est_' safe]) = est_mat(:,c);
    T_wide.(['p_'   safe]) = p_mat(:,c);
end
ref_row = T_wide(1,:);
ref_row.Excluded_Subject = categorical("<none>");
ref_row.Fit_OK = true;
for c = 1:nCoef
    safe = matlab.lang.makeValidName(char(coefNames(c)));
    ref_row.(['Est_' safe]) = ref_est(c);
    ref_row.(['p_'   safe]) = ref_p(c);
end
try
    T_wide = [ref_row; T_wide];
catch
end

writetable(T_wide, cfg.logfile, 'Sheet', sheetName(tag,'LOSO'), 'WriteMode','overwritesheet');
loso_tbl = T_wide;

% --- Forest plot per contrast ---
if ~(isfield(cfg,'make_plots') && cfg.make_plots), return; end

panel_w = max(380, 220 * nCoef);
fig = figure('Name', ['LOSO: ' tag], ...
             'Position', [80 80 panel_w 520], 'Color','w');
tiledlayout(fig, 1, nCoef, 'Padding','compact', 'TileSpacing','compact');

for c = 1:nCoef
    nexttile;
    hold on; grid on; box on
    ypos = 1:n_sub;
    for i = 1:n_sub
        if ~fit_ok(i), continue; end
        if p_mat(i,c) < 0.05
            col = [0.20 0.20 0.20];
        else
            col = [0.85 0.33 0.10];
        end
        errorbar(est_mat(i,c), ypos(i), ...
                 est_mat(i,c) - ci_lo_mat(i,c), ci_hi_mat(i,c) - est_mat(i,c), ...
                 'horizontal', 'o', 'Color', col, 'MarkerFaceColor', col, ...
                 'MarkerSize', 5, 'LineWidth', 0.9, 'CapSize', 4);
    end
    xline(ref_est(c), 'b-', 'LineWidth', 1.5);
    xline(0, 'k:', 'LineWidth', 0.8);
    yticks(ypos); yticklabels(cellstr(string(subs)));
    ylim([0.5 n_sub+0.5]); set(gca, 'FontSize', 7);
    xlabel('Estimate (95% CI)')
    n_lost = sum(fit_ok & ref_p(c) < 0.05 & p_mat(:,c) >= 0.05);
    title(sprintf('%s\nref beta=%.3g, p=%.3g\n%d / %d exclusions -> p>=0.05', ...
          coefNames(c), ref_est(c), ref_p(c), n_lost, sum(fit_ok)), ...
          'Interpreter','none', 'FontSize',8)
end
sgtitle(sprintf('LOSO sensitivity: %s  (method=%s)  |  blue=full-fit estimate, orange=p>=0.05 when excluded', ...
        tag, headline), 'FontWeight','bold', 'FontSize',11)

outpng = fullfile(cfg.figpath, [tag '_LOSO.png']);
exportgraphics(fig, outpng, 'Resolution', 200);
fprintf('Saved LOSO figure: %s\n', outpng);
end


function c = reorderCategoricalRef(c, refLvl, allLvls)
if ~iscategorical(c), c = categorical(c); end
present = categories(removecats(c));
if ~ismember(refLvl, present)
    warning('PLV:RefLevelMissing', ...
            'reorderCategoricalRef: reference level "%s" not present in data.', refLvl);
    return
end
others   = setdiff(allLvls, refLvl, 'stable');
others   = intersect(others, present, 'stable');
newOrder = [{refLvl}, others(:)'];
c = categorical(c, newOrder);
end
