%% Linear Mixed Models for FFR strength
% Outcome:  spectral SNR (dB) at F0, computed from averagedFFR_*.mat files
%
% Comparisons (toggle individually in the CONFIG block, or run any %% section
% on its own with Ctrl+Enter):
%   (A) 333 Hz, surface conditions:  sit vs firm vs foam
%   (B) 333 Hz, pre vs post:         sit vs sit_2
%   (C) 89  Hz, surface conditions:  sit vs firm vs foam
%   (D) 89 vs 333 Hz across surfaces
%
% LMM (per comparison):
%   SNR_dB ~ <fixed effects> + (1 | Subject) + (1 | Channel)
% 
% Effect sizes reported:
%   - Fixed-effect betas with 95% CIs
%   - Marginal & conditional R^2
%   - Cohen's d (paired) for each pairwise contrast
%
% Aline / Auditory-mov-data pipeline -- LMM module
% Requires: Statistics & Machine Learning Toolbox, FieldTrip on path (for ft_selectdata)

clc; clearvars -except ft_defaults_done; close all

%% CONFIG
CFG = struct();

% Paths
CFG.basepath   = 'F:\Auditory-mov-data';
CFG.outpath    = 'F:\Auditory-mov-data\LMM_SNR_robust\final';
CFG.subfolder  = 'preprocessed';

% Subjects
CFG.subjects   = {'NG', 'BVD', 'FD', 'EF', 'CX', 'HX', 'MCBG', 'LV', 'AM', 'FS', 'LVB', 'RNR', 'VUR', 'LH', 'TF', 'MGZ', 'NSMT', 'MM', 'TH', 'TD'};

%  Channels
% Either 'all' to use every EEG channel, a cell-array of labels to restrict
% to a specific ROI, or the special string 'vertex' (uses CFG.vertex_chans).
% Reference channels (TP9/TP10) and the trigger channel are always excluded.
CFG.channels      = 'all';
CFG.vertex_chans  = {'Cz','C1','C2','CPz','FC1','FC2','CP1','CP2'};   % vertex ROI definition

% SNR computation
CFG.fs              = 10000;          % sample rate (Hz)
CFG.snr_method      = 'custom';       % 'custom' matches reference script; 'simple' is the original LMM default
CFG.toi             = [0.050 0.210];  % steady-state window (used by 'simple' method only)
CFG.detrend_segment = true;           % linear detrend (used by 'simple' method only)

% Settings for 'custom' method:
CFG.custom_sig_band       = [-2 2];    % Hz around F0 for signal band (e.g. 331-335 at F0=333)

%CFG.custom_noise_lo_band  = [-12 -8];  % Hz offsets for lower noise band
CFG.custom_noise_lo_band  = [-30 -8];

%CFG.custom_noise_hi_band  = [ 8  12];  % Hz offsets for upper noise band
CFG.custom_noise_hi_band  = [8  30];

% Settings for 'simple' method:
CFG.noise_bw_hz     = 10;             % +/- Hz around F0 used to estimate noise floor
CFG.noise_gap_hz    = 2;              % +/- Hz around F0 to exclude from noise band
 
% Channel restriction
% If true, only include channels that are in SigChansPadded_*.mat (sigChansBoth)
% for the conditions contributing to a comparison.
% sig_combine = how to combine sig-channel lists across conditions per subject:
%   'union'     : channel must be significant in AT LEAST ONE condition (default)
%                 -> Keeps channels that drop below threshold in some conditions,
%                    so the LMM can quantify the actual drop-off.
%   'intersect' : channel must be significant in ALL conditions
%                 --> More conservative; can hide the effect you want to measure
%                    if a channel's response drops out in one condition.
CFG.restrict_to_sigchans = true;
CFG.sig_combine          = 'intersect';
 
% Model structure
% Whether to include (1|Channel) as a random effect.
% Set false for small ROIs (e.g. vertex with ~8 electrodes) where channel
% variance is undefined or parked at the boundary.
CFG.channel_random_effect = true;
 
% Outcome variable for the LMM
% 'SNR_dB'  : log-transformed SNR (dB), more normal residuals (recommended)
% 'SNR_lin' : raw linear SNR ratio (signal power / noise power)
CFG.outcome = 'SNR_dB';

% Outlier handling
% Two-pass fit: fit once, flag influential observations, refit without them.
% Both fits are reported (with/without outliers) for transparency.
%   method  : 'cooks'      -> Cook's distance (recommended)
%             'residual_z' -> |standardised residuals| > threshold
%   thresh  : numeric value, or 'auto' (cooks: 4/n,  residual_z: 3.0)
CFG.remove_outliers = true;    % Set to true to also get Cook's distance trimmed results
CFG.outlier_method  = 'cooks';
CFG.outlier_thresh  = 'auto';

%  Robust estimation (IRLS, Tukey's bisquare)
CFG.robust_estimation = true;     % fit a bisquare-weighted robust LMM alongside full/cooks-trim
CFG.robust_tuning      = 4.685;   % Tukey's bisquare tuning constant (95% efficiency under normality)
CFG.robust_maxiter     = 20;
CFG.robust_tol         = 1e-3;    % convergence: max change in any weight between iterations

% Which of the three fits (full / cooks / robust) feeds the headline
% contrasts, LOSO, assumption checks, and the main plotResults figure.
CFG.headline_fit = 'robust';      % 'full' | 'cooks' | 'robust'

% Plot labelling
% If true, the spaghetti panel labels each subject's line with its ID at the
% leftmost condition (mirrors the channel labels in the boxplot script).
% 'all'      : label every subject
% 'outliers' : label only subjects whose line is a statistical outlier in any
%              condition (1.5*IQR rule) -- keeps busy plots readable
CFG.label_subjects      = true;
CFG.label_subject_mode  = 'all';   % 'all' or 'outliers'
 
% Which analyses to run (set to false to skip)
CFG.run.A_within_333  = true;   % sit vs firm vs foam at 333 Hz
CFG.run.B_prepost_333 = true;   % sit vs sit_2 at 333 Hz
CFG.run.C_within_89   = true;   % sit vs firm vs foam at 89 Hz
CFG.run.D_between_freq= true;   % 89 vs 333 across conditions
 
% Logging
CFG.verbose_missing = false;    % set true to see each missing file (default: summary only)
 
% Plotting
CFG.make_plots = true;          % create explanatory figures per analysis
CFG.figpath    = fullfile('F:\Auditory-mov-data', 'LMM_SNR_robust', 'figures-final');
if CFG.make_plots && ~exist(CFG.figpath,'dir'), mkdir(CFG.figpath); end
 
% Assumption checks
CFG.checkAssumptions = true;   % LMM diagnostics figure + numerical tests per analysis

% Leave-one-subject-out sensitivity
% Refits the model once per subject (each time with that subject excluded)
% to test whether the fixed-effect conclusions depend on any single subject.
% Adds ~n_subj extra LMM fits per analysis, so it roughly doubles runtime.
CFG.run_loso = true;
 
% File-name templates
% Conditions A & C use these surface names; pre/post uses sit and sit_2
CFG.surfaces      = {'sit','firm','foam'};
CFG.prepost_names = {'sit','sit_2'};   % first = pre, second = post
 
% Output
if ~exist(CFG.outpath,'dir'), mkdir(CFG.outpath); end
CFG.logfile = fullfile(CFG.outpath, 'FFR_LMM_results.xlsx');

%% BUILD THE DATA TABLE 
% This crawls every needed file once and computes SNR + PLV per
% Subject x Condition x Frequency x Session x Channel.
 
fprintf('\n=== Building SNR table ===\n');
T = buildTable(CFG);
 
% Quick save of the raw long-format table
writetable(T, CFG.logfile, 'Sheet','data', 'WriteMode','overwritesheet');
fprintf('Saved long-format data to %s (sheet: data)\n', CFG.logfile);


%% (A) 333 Hz: sit vs firm vs foam
if CFG.run.A_within_333
    fprintf('\n=== (A) 333 Hz: sit vs firm vs foam ===\n');
    Ta = T(T.Frequency == 333 & T.Session == "pre" & ...
            ismember(T.Condition, {'sit','firm','foam'}), :);
    Ta = dropUnused(Ta, {'Condition','Subject','Channel'});
    if CFG.restrict_to_sigchans
        Ta.Condition = categorical(Ta.Condition);
        Ta = selectSigChans(Ta, 'Condition', CFG.sig_combine);
    end

    tag = 'A_333';
    % Reference level is arbitrary now -- sit is just a sensible default
    Ta.Condition = reorderCategoricalRef(Ta.Condition, 'sit', {'sit','firm','foam'});

    if checkFactor(Ta, 'Condition', 2, tag)
        chanRE = '';
        if CFG.channel_random_effect, chanRE = ' + (1|Channel)'; end
        formula = [CFG.outcome ' ~ Condition + (1|Subject)' chanRE];
        res_A   = fitLMM(Ta, formula, tag, CFG);

        % ---- Three pairwise contrasts from the ONE fitted model ----
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
        padj = bonferroniAdjust(pvals);   % family of 3, within-frequency surface contrasts
        for k = 1:numel(contrasts)
            contrasts(k).pValue_bonf = padj(k);
        end
        res_A.contrasts = buildContrastsTable(contrasts, 'pValue_bonf', padj);
        disp(res_A.contrasts);

        res_A.cohens_d = cohensD(Ta, 'Condition', CFG.outcome);
        writeResults(res_A, tag, CFG);
        if CFG.run_loso
            leaveOneSubjectOut(res_A, tag, CFG);
        end
        if CFG.checkAssumptions
            checkAssumptions(res_A, Ta, 'Condition', tag, CFG);
        end
        if CFG.make_plots
                plotOutlierComparison(res_A, tag, CFG);
                plotRobustWeights(res_A, tag, CFG);
                plotResults(Ta, res_A, 'Condition', tag, '333 Hz: sit vs firm vs foam', CFG);
        end
    end
end


%%  (B) 333 Hz: pre vs post (sit vs sit_2)
if CFG.run.B_prepost_333
    fprintf('\n=== (B) 333 Hz: sit (pre) vs sit_2 (post) ===\n');
    Tb = T(T.Frequency==333 & T.Condition=="sit", :);
    Tb.Session = categorical(Tb.Session, {'pre','post'});
    Tb = dropUnused(Tb, {'Session','Subject','Channel'});

    Tb = keepPaired(Tb, 'Session');
    if CFG.restrict_to_sigchans
        Tb = selectSigChans(Tb, 'Session', CFG.sig_combine);
    end

    tag = 'B_333_prepost';   % <-- NEW: own tag, used consistently below

    if checkFactor(Tb, 'Session', 2, tag)
        chanRE = '';
        if CFG.channel_random_effect, chanRE = ' + (1|Channel)'; end
        formula = [CFG.outcome ' ~ Session + (1|Subject)' chanRE];
        res_B   = fitLMM(Tb, formula, tag, CFG);

        lme = res_B.lme;
        L   = makeContrast(lme, {'Session_post'});
        r   = contrastTest(lme, L);

        contrasts = struct( ...
            'Label', {'post_vs_pre'}, ...
            'L',     {L}, ...
            'Estimate', r.Estimate, 'SE', r.SE, 'tStat', r.tStat, 'DF', r.DF, ...
            'pValue', r.pValue, 'CI_lower', r.CI_lower, 'CI_upper', r.CI_upper);

        % Single comparison -- no real correction needed, but pValue_bonf column
        % kept for schema consistency across all four analyses' _contrasts sheets.
        res_B.contrasts = buildContrastsTable(contrasts, 'pValue_bonf', contrasts.pValue);
        disp(res_B.contrasts);

        res_B.cohens_d = cohensD(Tb, 'Session', CFG.outcome);
        writeResults(res_B, tag, CFG);
        if CFG.run_loso
            leaveOneSubjectOut(res_B, tag, CFG);
        end
        if CFG.checkAssumptions
            checkAssumptions(res_B, Tb, 'Session', tag, CFG);
        end
        if CFG.make_plots
            plotOutlierComparison(res_B, tag, CFG);
            plotRobustWeights(res_B, tag, CFG);
            plotResults(Tb, res_B, 'Session', tag, '333 Hz: sit-pre vs sit-post', CFG);   % <-- fixed: 'Session', not 'Condition'
        end
    end
end


%% (C) 89 Hz: sit vs firm vs foam
if CFG.run.C_within_89
    fprintf('\n=== (C) 89 Hz: sit vs firm vs foam ===\n');
    Tc = T(T.Frequency == 89 & T.Session == "pre" & ...
            ismember(T.Condition, {'sit','firm','foam'}), :);
    Tc = dropUnused(Tc, {'Condition','Subject','Channel'});
    if CFG.restrict_to_sigchans
        Tc.Condition = categorical(Tc.Condition);
        Tc = selectSigChans(Tc, 'Condition', CFG.sig_combine);
    end

    tag = 'C_89';
    % Reference level is arbitrary now -- sit is just a sensible default
    Tc.Condition = reorderCategoricalRef(Tc.Condition, 'sit', {'sit','firm','foam'});

    if checkFactor(Tc, 'Condition', 2, tag)
        chanRE = '';
        if CFG.channel_random_effect, chanRE = ' + (1|Channel)'; end
        formula = [CFG.outcome ' ~ Condition + (1|Subject)' chanRE];
        res_C   = fitLMM(Tc, formula, tag, CFG);

        % ---- Three pairwise contrasts from the ONE fitted model ----
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
        padj = bonferroniAdjust(pvals);   % family of 3, within-frequency surface contrasts
        for k = 1:numel(contrasts)
            contrasts(k).pValue_bonf = padj(k);
        end
        res_C.contrasts = buildContrastsTable(contrasts, 'pValue_bonf', padj);
        disp(res_C.contrasts);

        res_C.cohens_d = cohensD(Tc, 'Condition', CFG.outcome);
        writeResults(res_C, tag, CFG);
        if CFG.run_loso
            leaveOneSubjectOut(res_C, tag, CFG);
        end
        if CFG.checkAssumptions
            checkAssumptions(res_C, Tc, 'Condition', tag, CFG);
        end
        if CFG.make_plots
            plotOutlierComparison(res_C, tag, CFG);
            plotRobustWeights(res_C, tag, CFG);
            plotResults(Tc, res_C, 'Condition', tag, '89 Hz: sit vs firm vs foam', CFG);
        end
    end
end



%% (D) 89 vs 333 across conditions
if CFG.run.D_between_freq
    fprintf('\n=== (D) Frequency: 89 vs 333 (across conditions) ===\n');
    Td = T(T.Session=="pre" & ismember(T.Condition, {'sit','firm','foam'}), :);
    Td.Frequency = categorical(Td.Frequency);
    Td = dropUnused(Td, {'Condition','Frequency','Subject','Channel'});
    if CFG.restrict_to_sigchans
        Td.Condition = categorical(Td.Condition);
        Td.CondFreq = categorical(strcat(string(Td.Condition), '_', string(Td.Frequency)));
        Td = selectSigChans(Td, 'CondFreq', CFG.sig_combine);
        Td.CondFreq = [];
    end

    tag = 'D_FxC';
    Td.Condition = reorderCategoricalRef(Td.Condition, 'sit', {'sit','firm','foam'});
    Td.Frequency = reorderCategoricalRef(Td.Frequency, '89', {'89','333'});

    if checkFactor(Td, 'Frequency', 2, tag) && checkFactor(Td, 'Condition', 2, tag)
        chanRE = '';
        if CFG.channel_random_effect, chanRE = ' + (1|Channel)'; end
        formula = [CFG.outcome ' ~ Frequency * Condition + (1|Subject)' chanRE];
        res_D   = fitLMM(Td, formula, tag, CFG);

        lme = res_D.lme;
        % Frequency effect (333-89) evaluated AT EACH surface condition:
        %   at sit  : Frequency_333 alone (sit is the Condition reference)
        %   at firm : Frequency_333 + Frequency_333:Condition_firm
        %   at foam : Frequency_333 + Frequency_333:Condition_foam
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
        padj = fdrBH(pvals);   % FDR family: 333-vs-89, tested once per surface
        for k = 1:numel(contrasts)
            contrasts(k).pValue_fdr = padj(k);
        end
        res_D.contrasts = buildContrastsTable(contrasts, 'pValue_fdr', padj);        disp(res_D.contrasts);

        res_D.cohens_d = cohensD(Td, 'Frequency', CFG.outcome);
        writeResults(res_D, tag, CFG);
        if CFG.run_loso
            leaveOneSubjectOut(res_D, tag, CFG);
        end
        if CFG.checkAssumptions
            checkAssumptions(res_D, Td, 'Frequency', tag, CFG);
        end
        if CFG.make_plots
            plotOutlierComparison(res_D, tag, CFG);
            plotRobustWeights(res_D, tag, CFG);
            plotResults(Td, res_D, 'Frequency', tag, '89 vs 333 Hz across conditions', CFG);   % <-- fixed: 'Frequency', not 'Condition'; title corrected
        end
    end
end



%% Helper funcitons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function L = makeContrast(lme, posNames, negNames)
% Build a contrast row vector over the model's fixed-effect coefficient space.
% L*beta = sum(coef for posNames) - sum(coef for negNames)
% Coefficient names must match lme.CoefficientNames exactly (case-sensitive).
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


function T = buildContrastsTable(contrasts, padjName, padjVals)
% Convert a contrasts struct array into a table, forcing the L column to
% remain a cell array regardless of whether all L vectors happen to be
% the same length. (struct2table silently merges same-length numeric
% fields into a matrix column instead, which breaks any downstream code
% that indexes L with C.L{c}.)
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


function out = contrastTest(lme, L, alpha)
% Test one linear combination of fixed effects: H0: L*beta = 0.
% Returns estimate, SE, t, df (Satterthwaite), p, and 95% CI --
% all from the SAME fitted model, regardless of which level was the reference.
if nargin < 3, alpha = 0.05; end

beta = fixedEffects(lme);
covB = lme.CoefficientCovariance;

est = L * beta;
se  = sqrt(L * covB * L');

[pval, Fstat, df] = coefTest(lme, L, 0, 'DFMethod', 'satterthwaite');
tstat = sign(est) * sqrt(max(Fstat, 0));   % single-row contrast: F = t^2

tcrit = tinv(1 - alpha/2, df);
ci_lo = est - tcrit*se;
ci_hi = est + tcrit*se;

out = struct('Estimate', est, 'SE', se, 'tStat', tstat, 'DF', df, ...
             'pValue', pval, 'CI_lower', ci_lo, 'CI_upper', ci_hi);
end


function padj = bonferroniAdjust(pvals)
% Classic Bonferroni: multiply by family size, cap at 1.
m = numel(pvals);
padj = min(pvals * m, 1);
end


function padj = fdrBH(pvals)
% Benjamini-Hochberg FDR correction (manual implementation, no
% Bioinformatics Toolbox dependency).
m = numel(pvals);
[psort, ord] = sort(pvals(:));
padj_sorted = psort .* m ./ (1:m)';
% Enforce monotonicity (step-up procedure)
padj_sorted = flipud(cummin(flipud(padj_sorted)));
padj = nan(size(pvals));
padj(ord) = min(padj_sorted, 1);
end


function ok = checkFactor(T, varName, minLevels, tag)
% Verify a factor column has enough levels and observations to fit a model.
% Returns true if OK, false (with a clear message) if not.
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
% Remove unused categorical levels so the design matrix is full rank.
for i = 1:numel(vars)
    if ismember(vars{i}, T.Properties.VariableNames) && ...
       iscategorical(T.(vars{i}))
        T.(vars{i}) = removecats(T.(vars{i}));
    end
end
end
 
 
function T = keepPaired(T, factorVar)
% For paired designs: keep only subjects that have data in ALL levels of factorVar.
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
% For each subject, keep only channels that pass a significance criterion
% across the levels of factorVar.
%   mode = 'intersect' : channel must be in sigChans for EVERY level (conservative)
%   mode = 'union'     : channel must be in sigChans for AT LEAST ONE level (default)
%
% NOTE: the rows in T already come from sig-channel-filtered loading
% (CFG.restrict_to_sigchans). This second pass decides how to combine the
% per-condition lists *within* a subject. With 'union', a channel can appear
% in the table for a condition where its own Rayleigh test was non-significant,
% as long as it was significant in at least one other condition for the same
% subject -- and we then re-attach that channel's data from the underlying
% averaged files.
%
% IMPORTANT: rows in T have already been pre-filtered, so for 'union' to
% actually recover dropped channel-condition cells, buildTable must
% have included those rows. We handle that by leaving buildTable to
% emit rows for all channels first (sig + non-sig), and using a marker
% column 'IsSigHere' (boolean) that tells us whether the original Rayleigh
% test flagged this channel in this condition. Below, we filter by:
%   - intersect: keep rows where IsSigHere==true AND channel passed in all levels
%   - union:     keep rows where channel passed in at least one level
% In both cases, channels never significant for this subject are dropped.
 
if ~ismember(factorVar, T.Properties.VariableNames) || isempty(T), return; end
if nargin < 3 || isempty(mode), mode = 'union'; end
 
% Need the IsSigHere column to do this properly
if ~ismember('IsSigHere', T.Properties.VariableNames)
    warning('selectSigChans: IsSigHere column missing; skipping.');
    return
end
 
lv   = categories(removecats(T.(factorVar)));
subs = unique(T.Subject);
keep_rows = false(height(T), 1);
 
for i = 1:numel(subs)
    sub_rows = T.Subject == subs(i);
    if ~any(sub_rows), continue; end
 
    % Tally significance per channel across the levels of factorVar
    sub_chans = unique(string(T.Channel(sub_rows)));
    sig_count = zeros(numel(sub_chans), 1);
    for cidx = 1:numel(sub_chans)
        ch_rows = sub_rows & T.Channel == sub_chans(cidx);
        % count distinct levels where this channel was significant
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
 
 
function T = buildTable(CFG)
% Crawl files and compute SNR (dB) + PLV per subject/condition/channel/freq
 
rows = {};
rowi = 0;
 
% Track missing files per combo (initialized once we know the combo list)
missing_count = [];
combos_template = {};
 
for s = 1:numel(CFG.subjects)
    subj = CFG.subjects{s};
    inpath = fullfile(CFG.basepath, subj, CFG.subfolder);
 
    % Build the surface x frequency x session list, but ONLY for analyses
    % that are actually enabled in CFG.run. This avoids loading (and warning
    % about) files for conditions that wouldn't be used downstream anyway --
    % e.g. no need to crawl 89 Hz files when only A_within_333 is true.
    need_333_surf = CFG.run.A_within_333 || CFG.run.B_prepost_333 || CFG.run.D_between_freq;
    need_sit2     = CFG.run.B_prepost_333;
    need_89_surf  = CFG.run.C_within_89  || CFG.run.D_between_freq;

    combos = {};
    if need_333_surf
        for k = 1:numel(CFG.surfaces)
            combos(end+1,:) = {CFG.surfaces{k}, 333, 'pre'};   %#ok<AGROW>
        end
    end
    if need_sit2
        combos(end+1,:) = {'sit_2', 333, 'post'};
    end
    if need_89_surf
        for k = 1:numel(CFG.surfaces)
            combos(end+1,:) = {CFG.surfaces{k}, 89, 'pre'};   %#ok<AGROW>
        end
    end
 
    if isempty(missing_count)
        missing_count   = zeros(size(combos,1),1);
        combos_template = combos;
    end
 
    for r = 1:size(combos,1)
        surf = combos{r,1};
        f0   = combos{r,2};
        sess = combos{r,3};
 
        % File naming follows the preprocessing script.
        % Convention: averagedFFR_P[<subj>]_FFR_<surf>_<freq>[_2].mat
        % The "_2" suffix (post measurement) comes AFTER the frequency.
        if strcmp(surf, 'sit_2')
            name = ['_FFR_sit_', num2str(f0), '_2'];
        else
            name = ['_FFR_', surf, '_', num2str(f0)];
        end
        inputname = ['averagedFFR_P[', subj, ']', name, '.mat'];
        fpath     = fullfile(inpath, inputname);
 
        % Map sit_2 back to 'sit' for the Condition factor; Session
        % carries the pre/post info
        cond_factor = surf;
        if strcmp(surf, 'sit_2'), cond_factor = 'sit'; end
 
        if ~exist(fpath,'file')
            if CFG.verbose_missing
                warning('Missing: %s', fpath);
            end
            missing_count(r) = missing_count(r) + 1; %#ok<AGROW>
            continue
        end
 
        S = load(fpath);
        if ~isfield(S,'final')
            if CFG.verbose_missing
                warning('No "final" variable in %s', fpath);
            end
            missing_count(r) = missing_count(r) + 1; %#ok<AGROW>
            continue
        end
        D = S.final;   % FieldTrip timelock struct
 
        % Channel selection (initial)
        if ischar(CFG.channels) && strcmpi(CFG.channels,'all')
            ch_idx = 1:numel(D.label);
        elseif ischar(CFG.channels) && strcmpi(CFG.channels,'vertex')
            [~, ch_idx] = ismember(CFG.vertex_chans, D.label);
            ch_idx = ch_idx(ch_idx>0);
            missing_vx = setdiff(CFG.vertex_chans, D.label);
            if ~isempty(missing_vx)
                fprintf('  [%s/%s] Vertex channels not present in data: %s\n', ...
                        subj, surf, strjoin(missing_vx, ', '));
            end
        else
            [~, ch_idx] = ismember(CFG.channels, D.label);
            ch_idx = ch_idx(ch_idx>0);
        end
 
        % Always exclude reference electrodes from analysis channels
        % (TP9/TP10 may be flat after mastoid rereference)
        ref_chans = {'TP9','TP10','Trigger'};
        is_ref = ismember(D.label, ref_chans);
        ch_idx = ch_idx(~is_ref(ch_idx));
 
        % -------- Compute SNR --------
        fs = CFG.fs;
        switch lower(CFG.snr_method)
            case 'custom'
                % Replicates the reference script exactly
                avg     = D.avg(ch_idx, :);
                Nfull   = size(avg, 2);
                NFFT    = ceil(fs / 2);                % 2 Hz bin
                win     = hann(Nfull)';
                winPow  = sum(win .^ 2);
 
                S = fft(avg .* win, NFFT, 2);
                S = S(:, 1:floor(NFFT/2)+1);
                P = (abs(S).^2) / (winPow * fs);
                P(:,2:end-1) = 2 * P(:,2:end-1);
 
                df    = fs / NFFT;
                faxis = (0:floor(NFFT/2)) * df;
 
                sig_mask     = faxis >= (f0 + CFG.custom_sig_band(1)) & ...
                               faxis <= (f0 + CFG.custom_sig_band(2));
                noise_lo     = faxis >= (f0 + CFG.custom_noise_lo_band(1)) & ...
                               faxis <= (f0 + CFG.custom_noise_lo_band(2));
                noise_hi     = faxis >= (f0 + CFG.custom_noise_hi_band(1)) & ...
                               faxis <= (f0 + CFG.custom_noise_hi_band(2));
                noise_mask   = noise_lo | noise_hi;
 
                snr_lin = mean(P(:, sig_mask), 2) ./ ...
                          max(mean(P(:, noise_mask), 2), eps);
                snr_db  = 10 * log10(snr_lin);
 
            case 'simple'
                % Original LMM-script method
                t   = D.time;
                idx = t >= CFG.toi(1) & t <= CFG.toi(2);
                X   = D.avg(ch_idx, idx);
                if CFG.detrend_segment
                    X = detrend(X','linear')';
                end
                N    = size(X,2);
                nfft = 2^nextpow2(N*4);
                Y    = fft(X, nfft, 2);
                P    = abs(Y(:,1:nfft/2+1)).^2;
                faxis = (0:nfft/2) * fs/nfft;
 
                [~, kF0] = min(abs(faxis - f0));
                sig_pow  = P(:, kF0);
                in_band  = abs(faxis - f0) <= CFG.noise_bw_hz;
                in_guard = abs(faxis - f0) <= CFG.noise_gap_hz;
                noise_mask = in_band & ~in_guard;
                noise_pow  = mean(P(:, noise_mask), 2);
 
                snr_lin = sig_pow ./ max(noise_pow, eps);
                snr_db  = 10*log10(snr_lin);
 
            otherwise
                error('Unknown CFG.snr_method: %s (use ''custom'' or ''simple'')', ...
                      CFG.snr_method);
        end
 
        % -------- Optional: restrict to channels flagged as significant --------
        chan_labels = D.label(ch_idx);
        is_sig      = true(numel(chan_labels), 1);   % default: keep all
        if CFG.restrict_to_sigchans
            sigfile = fullfile(inpath, ['SigChansPadded_P[', subj, ']', name, '.mat']);
            if exist(sigfile, 'file')
                Ssig = load(sigfile, 'sigChansBoth');
                if isfield(Ssig,'sigChansBoth') && ~isempty(Ssig.sigChansBoth)
                    is_sig = ismember(chan_labels, Ssig.sigChansBoth);
                else
                    if CFG.verbose_missing
                        warning('SigChansBoth empty for %s -- keeping no channels.', sigfile);
                    end
                    is_sig = false(size(is_sig));
                end
            else
                if CFG.verbose_missing
                    warning('Missing sig-channel file %s -- keeping all channels for this condition.', sigfile);
                end
            end
        end
 
        % Append rows for ALL channels with an IsSigHere flag.
        % This lets downstream code choose 'union' (keep channels significant
        % in at least one condition) vs 'intersect' (significant in all).
        for c = 1:numel(chan_labels)
            rowi = rowi + 1;
            rows(rowi,:) = { string(subj), string(cond_factor), f0, ...
                             string(sess), string(chan_labels{c}), ...
                             snr_db(c), snr_lin(c), is_sig(c) }; %#ok<AGROW>
        end
    end
end
 
T = cell2table(rows, 'VariableNames', ...
    {'Subject','Condition','Frequency','Session','Channel','SNR_dB','SNR_lin','IsSigHere'});
T.Subject   = categorical(T.Subject);
T.Condition   = categorical(T.Condition);
T.Session   = categorical(T.Session);
T.Channel   = categorical(T.Channel);
T.IsSigHere = logical(T.IsSigHere);
 
% --- Missing file summary ---
nSubj = numel(CFG.subjects);
fprintf('\n--- Missing-file summary (out of %d subjects) ---\n', nSubj);
for r = 1:size(combos_template,1)
    fprintf('  %-6s %3d Hz %-4s : %d/%d subjects missing\n', ...
        combos_template{r,1}, combos_template{r,2}, combos_template{r,3}, ...
        missing_count(r), nSubj);
end
fprintf('Table built: %d rows, %d unique subjects.\n', ...
        height(T), numel(unique(T.Subject)));
 
% --- Per-subject/condition channel summary (useful for vertex ROI) ---
fprintf('\n--- Channels per subject x condition (after channel + sig-channel filter) ---\n');
G = findgroups(T.Subject, T.Condition, T.Frequency, T.Session);
S = splitapply(@(c) {sort(unique(string(c)))}, T.Channel, G);
[subL, condL, freqL, sessL] = groupLabels(T, G);
for k = 1:numel(S)
    chans_str = strjoin(S{k}, ',');
    fprintf('  %-6s %-6s %3d %-4s : %2d ch  [%s]\n', ...
        subL(k), condL(k), freqL(k), sessL(k), numel(S{k}), chans_str);
end
end
 
 
function res = fitLMM(Tin, formula, tag, CFG)
res = struct();
res.tag     = tag;
res.formula = formula;

% ---------- Fit 1: full data, unweighted ----------
res.full = fitOneLMM(Tin, formula);
printLMM(res.full, 'full data');

% ---------- Fit 2: Cook's-distance exclude & refit ----------
if isfield(CFG,'remove_outliers') && CFG.remove_outliers
    [keepMask, flagInfo] = flagOutliers(res.full.lme, Tin, CFG);
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

% ---------- Fit 3: robust (IRLS bisquare) ----------
if isfield(CFG,'robust_estimation') && CFG.robust_estimation
    try
        res.robust = fitRobustLMM(Tin, formula, CFG);
        printLMM(res.robust, sprintf('robust (bisquare, %d iter)', res.robust.n_iter));
    catch ME
        warning('Robust fit failed: %s', ME.message);
        res.robust = res.full;
    end
else
    res.robust = res.full;
end

% ---------- Select headline fit ----------
headline = 'robust';
if isfield(CFG,'headline_fit') && ~isempty(CFG.headline_fit)
    headline = lower(CFG.headline_fit);
end
switch headline
    case 'full',   main = res.full;
    case 'cooks',  main = res.trimmed;
    case 'robust', main = res.robust;
    otherwise
        warning('Unknown CFG.headline_fit "%s"; defaulting to robust.', headline);
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
% Helper: fit a single LMM and compute all summary stats
out = struct();
out.n_obs   = height(T);
out.n_subj  = numel(unique(T.Subject));
 
lme = fitlme(T, formula, 'FitMethod','REML');
out.lme = lme;
 
% Fixed-effects table with 95% CI
[~,~,FE] = fixedEffects(lme, 'DFMethod','satterthwaite');
out.fixedEffects = FE;
 
% Random-effects variances
[psi, mse] = covarianceParameters(lme);
out.randomEffects_cov = psi;
out.residual_var      = mse;
 
% R^2 (Nakagawa & Schielzeth approximation)
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
 

function out = fitRobustLMM(T, formula, CFG)
% Iteratively reweighted LMM using Tukey's bisquare weight function.
% Approximates the Dornheim & Brazauskas (2011) robust-estimation idea:
% observations with large residuals are progressively downweighted
% across iterations, rather than excluded outright as in Cook's-trim.
%
% NOTE: this is an IRLS approximation built on top of fitlme's WLS-style
% 'Weights' option, not a dedicated joint robust-MLE mixed-model solver
% (MATLAB has no built-in equivalent to R's robustlmm). Worth stating
% this explicitly in the methods text so the estimator is described
% accurately.

maxIter = CFG.robust_maxiter;
tol     = CFG.robust_tol;
tuning  = CFG.robust_tuning;

n = height(T);
w = ones(n,1);   % start: equal weights = ordinary REML fit
lme = [];

for iter = 1:maxIter
    lme_new = fitlme(T, formula, 'FitMethod','REML', 'Weights', w);

    res_raw = residuals(lme_new, 'Conditional', true);
    % Robust scale via MAD (median absolute deviation), normalised to be
    % consistent with the SD under a normal distribution (/0.6745)
    s = median(abs(res_raw - median(res_raw))) / 0.6745;
    s = max(s, eps);

    u = res_raw ./ (tuning * s);
    w_new = (1 - u.^2).^2;
    w_new(abs(u) >= 1) = 0;
    w_new = max(w_new, 1e-6);   % tiny floor for numerical stability, not zero

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
 
 
function [keep, info] = flagOutliers(lme, T, CFG)
% Return boolean keep-mask + diagnostics struct.
% Method: 'cooks' (Cook's-distance analogue for LMM) or 'residual_z'.
info = struct();
n = height(T);
 
method = lower(CFG.outlier_method);
thr    = CFG.outlier_thresh;
 
% Get residuals and standardised residuals
res_raw = residuals(lme, 'Conditional', true);
res_std = res_raw ./ std(res_raw, 'omitnan');
 
switch method
    case 'cooks'
        % Cook's-distance analogue for LMM:
        % approx D_i = (h_i / (1-h_i)^2) * (r_std_i^2 / p_fix)
        % where h_i are diagonal entries of the hat-style projection on the
        % fixed-effects design and p_fix is the number of fixed coefficients.
        Xfx = designMatrix(lme, 'Fixed');
        try
            H = Xfx * ((Xfx' * Xfx) \ Xfx');
            h = diag(H);
        catch
            h = sum(Xfx.^2, 2) ./ max(sum(Xfx(:).^2), eps);  % crude fallback
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
 
% Record which observations got flagged (subject/condition/channel)
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
% Uses built-in meanEffectSize (R2022a+, Statistics Toolbox) when available
% so we also get a bootstrap 95% CI. Falls back to manual computation otherwise.
 
G  = Tin.(groupVar);
lv = categories(G);
% subject-level means
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
            % Built-in: paired Cohen's d with bootstrap 95% CI
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
 

function writeResults(res, tag, CFG)
% Write LMM results to Excel: full, cooks-trim (if run), and robust (if run)
% fits are all reported side-by-side for transparency. The HEADLINE fit
% (res.headline_name, set in fitLMM per CFG.headline_fit) gets its own
% clearly-named sheet -- this is the one to cite in the manuscript.
file = CFG.logfile;

headlineName = 'full';
if isfield(res,'headline_name') && ~isempty(res.headline_name)
    headlineName = res.headline_name;
end

% --- Fixed-effects: HEADLINE fit (the one to report) ---
FE = dataset2table(res.fixedEffects);
writetable(FE, file, 'Sheet', sheetName(tag, sprintf('FE_headline_%s', headlineName)), ...
           'WriteMode','overwritesheet');

% --- Determine which fits actually ran (beyond the always-present 'full') ---
hasTrim   = isfield(res,'trimmed') && isfield(res,'full') && ...
            ~isequal(res.full.n_obs, res.trimmed.n_obs);
hasRobust = isfield(res,'robust') && isfield(res.robust,'weights');  % weights = genuine robust fit, not a full-fit fallback

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

isHeadline = strcmpi(fitNames, headlineName);   % marks which row is the headline

S = table(string(fitNames)', n_obs', n_subj', R2m', R2c', isHeadline', ...
    'VariableNames', {'Fit','N_obs','N_subjects','R2_marginal','R2_conditional','IsHeadline'});
writetable(S, file, 'Sheet', sheetName(tag,'summary'), 'WriteMode','overwritesheet');

% --- Cook's-trim details (only if it actually ran and differed) ---
if hasTrim
    FE_full = dataset2table(res.full.fixedEffects);
    writetable(FE_full, file, 'Sheet', sheetName(tag,'FE_full'), 'WriteMode','overwritesheet');

    if isfield(res, 'outlier_info') && ~isempty(res.outlier_info.flagged_table)
        writetable(res.outlier_info.flagged_table, file, ...
                   'Sheet', sheetName(tag,'outliers'), 'WriteMode','overwritesheet');
    end
end

% --- Robust-fit details ---
if hasRobust
    FE_robust = dataset2table(res.robust.fixedEffects);
    writetable(FE_robust, file, 'Sheet', sheetName(tag,'FE_robust'), 'WriteMode','overwritesheet');

    Wtab = table((1:numel(res.robust.weights))', res.robust.weights, ...
        'VariableNames', {'ObsIndex','RobustWeight'});
    writetable(Wtab, file, 'Sheet', sheetName(tag,'robust_w'), 'WriteMode','overwritesheet');
end

% --- Cohen's d ---
if isfield(res,'cohens_d') && ~isempty(res.cohens_d)
    writetable(res.cohens_d, file, 'Sheet', sheetName(tag,'d'), 'WriteMode','overwritesheet');
end

% --- Contrasts (from the single headline fit) ---
if isfield(res,'contrasts') && ~isempty(res.contrasts)
    writetable(res.contrasts, file, 'Sheet', sheetName(tag,'contrasts'), 'WriteMode','overwritesheet');
end

fprintf('Wrote results for %s to %s (headline fit: %s)\n', tag, file, headlineName);
end


function name = sheetName(tag, suffix)
% Build an Excel-safe sheet name: <31 chars, no forbidden chars.
% Replaces forbidden chars (: \ / ? * [ ]) with '-' and truncates the
% tag if combined with the suffix it would exceed 31 chars.
forbidden = {':','\','/','?','*','[',']'};
for k = 1:numel(forbidden)
    tag    = strrep(tag,    forbidden{k}, '-');
    suffix = strrep(suffix, forbidden{k}, '-');
end
maxLen = 31;
fullName = [tag '_' suffix];
if numel(fullName) > maxLen
    keepTag = maxLen - numel(suffix) - 1;   % -1 for the underscore
    if keepTag < 1, keepTag = 1; end
    tag = tag(1:min(numel(tag), keepTag));
    fullName = [tag '_' suffix];
end
name = fullName;
end
 

function plotResults(Tin, res, groupVar, tag, ttl, CFG)
% Create four explanatory figures for one LMM analysis:
%   1) Subject-level means with paired lines (spaghetti) + group means
%   2) Per-channel mean SNR by condition (heatmap-like)
%   3) Fixed-effect estimates with 95% CIs (forest plot)
%   4) Random-effect spread (subject vs channel) + R^2 bar
%
% Saves a single combined figure as PNG in CFG.figpath.
 
G  = Tin.(groupVar);
lv = categories(removecats(G));
nL = numel(lv);
 
% Outcome variable name + axis label
outcomeVar = 'SNR_dB';
if isfield(CFG,'outcome') && ~isempty(CFG.outcome)
    outcomeVar = CFG.outcome;
end
switch outcomeVar
    case 'SNR_dB',  outcomeLabel = 'SNR (dB)';
    case 'SNR_lin', outcomeLabel = 'SNR (linear ratio)';
    otherwise,      outcomeLabel = outcomeVar;
end
 
% Subject-level means matrix: rows = subjects, cols = condition levels
sub = unique(Tin.Subject);
M = nan(numel(sub), nL);
for i = 1:numel(sub)
    for j = 1:nL
        sel = Tin.Subject==sub(i) & G==lv{j};
        if any(sel), M(i,j) = mean(Tin.(outcomeVar)(sel),'omitnan'); end
    end
end
 
fig = figure('Name', ['LMM results: ' tag], ...
             'Position', [80 80 1500 900], 'Color','w');
tiledlayout(fig, 2, 2, 'Padding','compact', 'TileSpacing','compact');
 
% =========================================================
% Panel 1: paired-subject plot (spaghetti) + group means
% =========================================================
nexttile;
hold on; grid on; box on
xpos = 1:nL;
 
% Decide which subjects to label (mirrors the channel-label logic from the
% boxplot script). 'all' labels every subject; 'outliers' labels only those
% whose subject-mean is a 1.5*IQR outlier in at least one condition.
do_label = isfield(CFG,'label_subjects') && CFG.label_subjects;
label_mode = 'all';
if isfield(CFG,'label_subject_mode') && ~isempty(CFG.label_subject_mode)
    label_mode = lower(CFG.label_subject_mode);
end
label_this = false(numel(sub),1);
if do_label
    switch label_mode
        case 'outliers'
            for j = 1:nL
                col = M(:,j);
                isout = isoutlier(col, 'quartiles');   % 1.5*IQR, NaNs ignored
                label_this = label_this | isout;
            end
        otherwise   % 'all'
            label_this = ~all(isnan(M),2);   % label any subject with data
    end
end
 
% Pre-compute which subjects are statistical outliers (1.5*IQR in any
% condition). Outliers get a distinct line colour AND are always labelled,
% so the most informative IDs stand out from the grey mass.
is_outlier = false(numel(sub),1);
for j = 1:nL
    col = M(:,j);
    is_outlier = is_outlier | isoutlier(col, 'quartiles');
end

% Spaghetti per subject. Each line is its own object carrying a hover data
% tip with the subject ID, so the saved .fig is interactive too.
for i = 1:numel(sub)
    if is_outlier(i)
        lc = [0.85 0.33 0.10];   % orange-red highlight
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
        % older MATLAB without DataTipTemplate: skip silently
    end
end
 
% Group mean +/- 95% CI
grp_mean = mean(M, 1, 'omitnan');
grp_sd   = std(M, 0, 1, 'omitnan');
nValid   = sum(~isnan(M), 1);
grp_se   = grp_sd ./ sqrt(max(nValid,1));
grp_ci   = 1.96 * grp_se;
 
errorbar(xpos, grp_mean, grp_ci, 'ko', 'LineWidth', 2.0, ...
    'MarkerSize', 10, 'MarkerFaceColor','k', 'CapSize', 12);
 
xticks(xpos); xticklabels(lv);
% Extend x-axis on the right to make room for subject labels
xlim([0.5 nL + 0.9])
xlabel(groupVar); ylabel(outcomeLabel)

% -------- Right-side subject labels with vertical de-overlap --------
% Each label is anchored at the subject's LAST non-NaN value, then a greedy
% pass enforces a minimum vertical spacing so IDs never overlap. When a
% label gets pushed away from its true y, a thin leader line connects them.
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
        % Sort by data y, then push upward to enforce min spacing
        [~, ord] = sort([pts.y_data]);
        pts = pts(ord);
        yl  = ylim;
        min_spacing = (yl(2) - yl(1)) * 0.032;   % ~3% of axis range per label
        for k = 2:numel(pts)
            if pts(k).y_lab < pts(k-1).y_lab + min_spacing
                pts(k).y_lab = pts(k-1).y_lab + min_spacing;
            end
        end
        % If the top one got pushed above the axis, walk back down
        if pts(end).y_lab > yl(2)
            pts(end).y_lab = yl(2);
            for k = numel(pts)-1:-1:1
                if pts(k).y_lab > pts(k+1).y_lab - min_spacing
                    pts(k).y_lab = pts(k+1).y_lab - min_spacing;
                end
            end
        end

        x_data_edge = nL;            % rightmost data column
        x_leader    = nL + 0.18;     % elbow of the leader line
        x_text      = nL + 0.25;     % where text starts
        for k = 1:numel(pts)
            % Leader line only if the label was actually displaced
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
% Panel 2: per-channel mean SNR by condition (heatmap)
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
 
% Y ticks: show all if few channels, otherwise a subset of unique positions
if numel(chans) <= 32
    yt = 1:numel(chans);
    yticks(yt); yticklabels(chans);
else
    yt = unique(round(linspace(1, numel(chans), 16)));   % unique => strictly increasing
    yticks(yt); yticklabels(chans(yt));
    set(gca,'FontSize',8)
end
 
% =========================================================
% Panel 3: forest plot of fixed-effect estimates
% =========================================================
nexttile;
hold on; grid on; box on
FE = res.fixedEffects;
% Skip the intercept for visibility
isInt = strcmp(string(FE.Name), '(Intercept)');
FEp = FE(~isInt, :);
nFE = height(FEp);
 
% Use ascending positions (yticks requires strictly increasing values),
% then flip the table so the first effect appears at the top.
ypos = 1:nFE;
FEp  = flipud(FEp);
 
errorbar(FEp.Estimate, ypos, FEp.Estimate - FEp.Lower, FEp.Upper - FEp.Estimate, ...
    'horizontal', 'ko', 'LineWidth', 1.6, 'MarkerSize', 8, 'MarkerFaceColor','k', 'CapSize', 10);
xline(0, 'r--', 'LineWidth', 1.2);
if nFE >= 1
    yticks(ypos); yticklabels(cellstr(FEp.Name));
    ylim([0.5 nFE+0.5])
end
xlabel(sprintf('Estimate (%s), 95%% CI', strrep(outcomeLabel,'SNR ','')))
title(sprintf('Fixed effects (relative to ref. level)\nR^2_{marg}=%.3f, R^2_{cond}=%.3f', ...
              res.R2_marginal, res.R2_conditional))
 
% Annotate significance
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
% Panel 4: variance decomposition + Cohen's d bar
% =========================================================
nexttile;
 
% Compute the three variance components
[psi, mse] = covarianceParameters(res.lme);
var_sub = 0; var_chn = 0;
for k = 1:numel(psi)
    M2 = psi{k};
    if ~isempty(M2)
        % Heuristic: first random group is Subject, second is Channel
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
b.CData = [0.20 0.45 0.85;   % fixed
           0.85 0.55 0.20;   % subject
           0.45 0.75 0.45;   % channel
           0.55 0.55 0.55];  % residual
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
 
% Overall title
sgtitle(sprintf('%s   |   n=%d obs, %d subjects, %d channels', ...
        ttl, res.n_obs, res.n_subj, numel(chans)), ...
        'FontWeight','bold', 'FontSize',13)
 
% Save
outpng = fullfile(CFG.figpath, [tag '.png']);
exportgraphics(fig, outpng, 'Resolution', 200);
% Interactive .fig so the subject-line hover data tips remain usable
savefig(fig, fullfile(CFG.figpath, [tag '.fig']));
fprintf('Saved figure: %s\n', outpng);
end
 

function checkAssumptions(res, Tin, groupVar, tag, CFG)
% Run a panel of LMM assumption checks and write results to Excel + a PNG.
%   1) Residual normality      : Q-Q + histogram + Lilliefors test
%   2) Homoscedasticity        : residuals vs fitted, scale-location, Levene
%   3) Linearity / influence   : Cook's-distance-style index, residuals vs group
%   4) Random-effect normality : Q-Q of subject and channel BLUPs
%   5) Multicollinearity (VIF) : on the fixed-effects design matrix
%
% Saves figure: <figpath>/<tag>_diagnostics.png
% Writes numerical summary to sheet: <tag>_assumptions in CFG.logfile
 
lme   = res.lme;
res_raw = residuals(lme, 'Conditional', true);    % residuals after random effects
res_std = res_raw ./ std(res_raw, 'omitnan');     % crude standardisation
yfit    = fitted(lme, 'Conditional', true);
 
% Outcome label (used for axis labels)
outcomeVar = 'SNR_dB';
if isfield(CFG,'outcome') && ~isempty(CFG.outcome), outcomeVar = CFG.outcome; end
switch outcomeVar
    case 'SNR_dB',  unitLabel = 'dB';
    case 'SNR_lin', unitLabel = 'linear';
    otherwise,      unitLabel = '';
end
 
% Use the model's own variables table so G aligns with residual length.
% After outlier removal, lme.Variables has fewer rows than the caller's Tin.
Tmodel = lme.Variables;
 
% ---------- Numerical tests ----------
% Lilliefors test on residuals (null = normal)
try
    [h_lil, p_lil] = lillietest(res_raw, 'Alpha', 0.05);
catch
    h_lil = NaN; p_lil = NaN;
end
 
% Levene's test for equal variance across the main grouping factor
% (built-in: vartestn, Brown-Forsythe absolute-deviation variant)
if ismember(groupVar, Tmodel.Properties.VariableNames)
    G = Tmodel.(groupVar);
else
    G = Tin.(groupVar);   % fallback (shouldn't normally happen)
end
try
    p_lev = vartestn(res_raw, G, 'TestType','LeveneAbsolute', 'Display','off');
    F_lev = NaN;   % vartestn returns p only; F-stat not exposed
catch
    p_lev = NaN; F_lev = NaN;
end
 
% VIF on the fixed-effects design matrix (excluding the intercept).
% Uses fitlm to regress each predictor on the others; VIF = 1 / (1 - R^2).
try
    Xfx = designMatrix(lme, 'Fixed');
    if all(Xfx(:,1) == 1), Xfx = Xfx(:, 2:end); end   % drop intercept column
    nc = size(Xfx, 2);
    if nc >= 2
        vif_vals = nan(nc, 1);
        for j = 1:nc
            mdl = fitlm(Xfx(:, setdiff(1:nc, j)), Xfx(:, j));
            vif_vals(j) = 1 / max(1 - mdl.Rsquared.Ordinary, eps);
        end
    else
        vif_vals = NaN;   % VIF undefined with a single predictor
    end
catch
    vif_vals = NaN;
end
 
% Print to console
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
 
% ---------- Figure ----------
fig = figure('Name', ['LMM diagnostics: ' tag], ...
             'Position', [60 60 1500 900], 'Color','w');
tiledlayout(fig, 2, 3, 'Padding','compact', 'TileSpacing','compact');
 
% Panel 1: Q-Q of residuals
nexttile;
qqplot(res_raw);
title(sprintf('Q-Q residuals\nLilliefors p=%.3g', p_lil))
grid on; box on
xlabel('Theoretical quantiles'); ylabel('Sample quantiles')
 
% Panel 2: histogram of residuals + normal overlay
nexttile;
histogram(res_raw, 40, 'Normalization','pdf', ...
    'FaceColor',[0.4 0.6 0.85], 'EdgeColor','w'); hold on
xs = linspace(min(res_raw), max(res_raw), 200);
ys = normpdf(xs, mean(res_raw,'omitnan'), std(res_raw,'omitnan'));
plot(xs, ys, 'r-', 'LineWidth', 1.5);
title('Residual distribution (red = normal)')
xlabel(sprintf('Residual (%s)', unitLabel)); ylabel('Density')
grid on; box on
 
% Panel 3: residuals vs fitted
nexttile;
plot(yfit, res_raw, '.', 'MarkerSize', 6, 'Color',[0.4 0.4 0.4]); hold on
yline(0, 'r--', 'LineWidth', 1);
% loess-style trend (moving average for visual aid)
[fs_sort, ord] = sort(yfit);
rs_sort        = res_raw(ord);
w              = max(20, round(length(fs_sort)/30));
smooth         = movmean(rs_sort, w, 'omitnan');
plot(fs_sort, smooth, 'b-', 'LineWidth', 1.5);
xlabel('Fitted values'); ylabel('Residuals')
title(sprintf('Residuals vs fitted\nLevene p=%.3g (by %s)', p_lev, groupVar))
grid on; box on
 
% Panel 4: scale-location plot
nexttile;
plot(yfit, sqrt(abs(res_std)), '.', 'MarkerSize', 6, 'Color',[0.4 0.4 0.4]); hold on
smooth2 = movmean(sqrt(abs(res_std(ord))), w, 'omitnan');
plot(fs_sort, smooth2, 'b-', 'LineWidth', 1.5);
xlabel('Fitted values'); ylabel('sqrt(|standardised residuals|)')
title('Scale-location (homoscedasticity)')
grid on; box on
 
% Panel 5: Q-Q of subject random effects (BLUPs)
nexttile;
try
    [Bsub] = randomEffects(lme);          % concatenated random effects
    % Best-effort: take the first len(unique subjects) elements as subject BLUPs
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
 
% Panel 6: residuals by group (boxplot)
nexttile;
% Defensive: ensure G and res_raw have the same length
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
 
% Title
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
sgtitle(sprintf('LMM assumption checks: %s   |   %s%s', tag, vifstr, verdict), ...
        'FontWeight','bold', 'FontSize',13)
 
outpng = fullfile(CFG.figpath, [tag '_diagnostics.png']);
exportgraphics(fig, outpng, 'Resolution', 200);
fprintf('Saved diagnostics figure: %s\n', outpng);
 
% ---------- Separate outlier figure (only if outlier removal was run) ----------
if isfield(res, 'outlier_info') && ~isempty(res.outlier_info)
    plotOutlierScores(res, tag, CFG);
end
 
% ---------- Write numerical summary to Excel ----------
if ~isnan(vif_vals(1))
    vif_str = strjoin(arrayfun(@(x) sprintf('%.2f',x), vif_vals, 'uni', 0), ', ');
else
    vif_str = 'n/a';
end
A = table( ...
    string(tag), p_lil, h_lil, p_lev, F_lev, string(vif_str), ...
    'VariableNames', {'Tag','Lilliefors_p','Lilliefors_reject', ...
                      'Levene_p','Levene_F','VIF_values'});
writetable(A, CFG.logfile, 'Sheet', sheetName(tag,'assumptions'), 'WriteMode','overwritesheet');
end


function [subL, condL, freqL, sessL] = groupLabels(T, G)
% Helper: for each group id in G, return the representative
% Subject/Condition/Frequency/Session label as parallel arrays.
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
 

function plotOutlierScores(res, tag, CFG)
% Stand-alone figure: outlier scores per observation + threshold line.
info = res.outlier_info;
if isempty(info), return; end
 
score = info.score;
thr   = info.threshold;
flag  = info.flagged_rows;
 
fig = figure('Name', ['Outliers: ' tag], ...
             'Position', [120 120 1100 500], 'Color','w');
tiledlayout(fig, 1, 2, 'Padding','compact', 'TileSpacing','compact');
 
% Panel 1: outlier score by observation index
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
 
% Panel 2: distribution
nexttile;
histogram(score, 60, 'FaceColor',[0.4 0.6 0.85], 'EdgeColor','w');
hold on
xline(thr, 'r--', 'LineWidth', 1.2);
xlabel(['Score (' info.method ')']); ylabel('Count')
title('Distribution of outlier scores')
grid on; box on
set(gca, 'YScale', 'log')   % outlier tails are easier to see on log
 
sgtitle(sprintf('Outlier diagnostics: %s', tag), 'FontWeight','bold', 'FontSize',12)
 
outpng = fullfile(CFG.figpath, [tag '_outliers.png']);
exportgraphics(fig, outpng, 'Resolution', 200);
fprintf('Saved outlier figure: %s\n', outpng);
end
 

function plotOutlierComparison(res, tag, CFG)
% Comparison of fitting approaches (full / [cooks-trim] / robust):
%   1) Forest plot: each stored contrast, re-evaluated under every fit
%   2) Marginal & conditional R^2 by fit
%   3) Robust (bisquare) weight per observation, flagging w < 0.5
%      (only drawn if a robust fit is present)
%
% Which fits are shown adapts to CFG.remove_outliers: when false,
% the redundant "cooks-trim == full" series is omitted entirely.

if ~isfield(res,'contrasts') || isempty(res.contrasts), return; end
if ~(isfield(CFG,'make_plots') && CFG.make_plots), return; end

% ---------- Decide which fits to show ----------
fits   = {};
fitRes = {};
allColors = struct( ...
    'full',       [0.00 0.4470 0.7410], ...
    'cookstrim',  [0.8500 0.3250 0.0980], ...
    'robust',     [0.4660 0.6740 0.1880]);

fits{end+1} = 'full';        fitRes{end+1} = res.full;

if isfield(CFG,'remove_outliers') && CFG.remove_outliers && isfield(res,'trimmed')
    fits{end+1} = 'cooks-trim';  fitRes{end+1} = res.trimmed;
end

hasRobust = isfield(CFG,'robust_estimation') && CFG.robust_estimation && ...
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

outpng = fullfile(CFG.figpath, [tag '_fitcompare.png']);
exportgraphics(fig, outpng, 'Resolution', 200);
fprintf('Saved outlier-comparison figure: %s\n', outpng);
end


function plotRobustWeights(res, tag, CFG)
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
if ~(isfield(CFG,'make_plots') && CFG.make_plots), return; end

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

outpng = fullfile(CFG.figpath, [tag '_robustw.png']);
exportgraphics(fig, outpng, 'Resolution', 200);
fprintf('Saved robust-weight diagnostics figure: %s\n', outpng);
end


function loso_tbl = leaveOneSubjectOut(res, tag, CFG)
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
                out_i = fitRobustLMM(Tsub, formula, CFG);
                lme_i = out_i.lme;
            case 'cooks'
                % Re-run the same exclude-then-refit logic used headline-side.
                lme0 = fitlme(Tsub, formula, 'FitMethod','REML');
                [keepMask_i, ~] = flagOutliers(lme0, Tsub, CFG);
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

writetable(T_wide, CFG.logfile, 'Sheet', sheetName(tag,'LOSO'), 'WriteMode','overwritesheet');
loso_tbl = T_wide;

% --- Forest plot per contrast ---
if ~(isfield(CFG,'make_plots') && CFG.make_plots), return; end

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

outpng = fullfile(CFG.figpath, [tag '_LOSO.png']);
exportgraphics(fig, outpng, 'Resolution', 200);
fprintf('Saved LOSO figure: %s\n', outpng);
end


function c = reorderCategoricalRef(c, refLvl, allLvls)
% Put refLvl first in the category ordering so it becomes the reference
% level for fitlme. Any other levels follow in the order given by allLvls
% (with levels not present in the data quietly skipped).
if ~iscategorical(c), c = categorical(c); end
present = categories(removecats(c));
if ~ismember(refLvl, present)
    warning('reorderCategoricalRef: reference level "%s" not present in data.', refLvl);
    return
end
% Build the new order: ref first, then the remaining allLvls that are present
others   = setdiff(allLvls, refLvl, 'stable');
others   = intersect(others, present, 'stable');
newOrder = [{refLvl}, others(:)'];
c = categorical(c, newOrder);
end