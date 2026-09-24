%% Linear Mixed Models for FFR PLV — CONTROL BAND (270-290 Hz)
%
% This is the control-band counterpart to MLM.m. It runs the same LMM
% analyses on PLV values computed in a non-stimulus frequency band
% (270-290 Hz), loaded from PLV_table_control_band.mat.
%
% PURPOSE:
% If the foam > firm/sit PLV effect at the stimulus frequency (333 Hz) is
% driven by contamination-induced phase-locking, the same effect should
% appear at non-stimulus frequencies. If the effect is specific to the
% stimulus frequency, this control analysis should show no condition
% effect (or a much smaller one) in the 270-290 Hz band.
%
% Interpretation of results:
%   - Flat / near-zero foam effect here  -> stimulus-specific effect, main finding defended
%   - Foam effect of similar magnitude   -> contamination is likely driving the main result
%   - Reduced but non-zero foam effect   -> mixed; partial contamination contribution
%
% NOTE: The control-band PLV table has Frequency=280 (band centre) instead
% of 333. This script remaps 280 -> 333 right after loading so all the
% existing filters (T.Frequency==333) work without further changes.
%
% Aline / Auditory-mov-data pipeline -- PLV LMM control-band module

clc; clear; close all

%% Configuration
cfg = struct();

% --- Paths ---
cfg.basepath   = 'F:\Auditory-mov-data';
cfg.outpath    = 'F:\Auditory-mov-data\LMM_PLV\control_band_270_290';

% Where compute_PLV.m saved the control-band long-format PLV table.
cfg.plv_table_file = fullfile(cfg.outpath, 'PLV_table_control_band.mat');

% --- Control band specification (used only for labelling) ---
cfg.control_band_label_str = '270-290 Hz';

% --- Subjects ---
cfg.subjects   = {'NG', 'BVD', 'FD', 'EF', 'CX', 'HX', 'MCBG', 'LV', 'AM', ...
                  'FS', 'LVB', 'RNR', 'VUR', 'LH', 'TF', 'MGZ', 'NSMT', ...
                  'MM', 'TH', 'TD'};

% --- Channel restriction ---
cfg.restrict_to_sigchans = true;
cfg.sig_combine          = 'union';

% --- Model structure ---
cfg.channel_random_effect = false;

% --- Outcome variable ---
cfg.outcome = 'PLV_logit';
cfg.logit_eps = 1e-4;

% --- Reference levels for the Condition factor ---
cfg.condition_refs = {'sit', 'firm'};

% --- Outlier handling ---
cfg.remove_outliers = true;
cfg.outlier_method  = 'cooks';
cfg.outlier_thresh  = 'auto';

% --- Plot labelling ---
cfg.label_subjects      = true;
cfg.label_subject_mode  = 'all';

% --- Which analyses to run ---
cfg.run.A_within_333  = true;   % sit vs firm vs foam, control band on 333 Hz data
cfg.run.B_prepost_333 = false;  % sit vs sit_2, control band on 333 Hz data
cfg.run.C_within_89   = false;  % sit vs firm vs foam, control band on 89 Hz data
cfg.run.D_between_freq= false;

% --- Plotting ---
cfg.make_plots = true;
cfg.figpath    = fullfile(cfg.outpath, 'figures_PLV_controlband');
if cfg.make_plots && ~exist(cfg.figpath,'dir'), mkdir(cfg.figpath); end

% --- Assumption checks ---
cfg.checkAssumptions = true;

% --- Leave-one-subject-out sensitivity ---
cfg.run_loso = true;

% --- Output ---
if ~exist(cfg.outpath,'dir'), mkdir(cfg.outpath); end
cfg.logfile = fullfile(cfg.outpath, 'FFR_LMM_results_PLV_controlband.xlsx');


%% Load PLV table

fprintf('\n=== Loading CONTROL-BAND PLV table from %s ===\n', cfg.plv_table_file);
if ~exist(cfg.plv_table_file, 'file')
    error(['Control-band PLV table not found at %s. ', ...
           'Run compute_PLV_controlband.m first.'], cfg.plv_table_file);
end
S = load(cfg.plv_table_file, 'T_plv');
T = S.T_plv;

% --- Remap control-band frequency label to match LMM filters ---
% The control-band PLV table uses Frequency=280 (band centre) but the
% LMM analyses filter on Frequency==333 (or 89). Remap so the existing
% filters work without further changes. This is purely a relabelling --
% the actual PLV values were computed from the 270-290 Hz band.
n_remapped_333 = sum(T.Frequency == 280);
T.Frequency(T.Frequency == 280) = 333;
fprintf('Remapped %d rows from Frequency=280 to Frequency=333 (for filter compatibility).\n', ...
        n_remapped_333);

% Make sure key columns are the expected types
if ~iscategorical(T.Subject),   T.Subject   = categorical(T.Subject);   end
if ~iscategorical(T.Condition), T.Condition = categorical(T.Condition); end
if ~iscategorical(T.Session),   T.Session   = categorical(T.Session);   end
if ~iscategorical(T.Channel),   T.Channel   = categorical(T.Channel);   end
if ~islogical(T.IsSigHere),     T.IsSigHere = logical(T.IsSigHere);     end

% Derive logit-transformed PLV
plv_clipped = max(min(T.PLV, 1 - cfg.logit_eps), cfg.logit_eps);
T.PLV_logit = log(plv_clipped ./ (1 - plv_clipped));

fprintf('Loaded %d rows, %d unique subjects, %d unique channels.\n', ...
        height(T), numel(unique(T.Subject)), numel(unique(T.Channel)));

% --- Per-condition summary ---
fprintf('\n--- Median PLV (control band) per Condition x Frequency x Session ---\n');
G2 = findgroups(T.Condition, T.Frequency, T.Session);
plv_meds = splitapply(@(x) median(x,'omitnan'), T.PLV, G2);
plv_iqr_lo = splitapply(@(x) quantile(x,0.25), T.PLV, G2);
plv_iqr_hi = splitapply(@(x) quantile(x,0.75), T.PLV, G2);
n_obs      = splitapply(@(x) numel(x), T.PLV, G2);
[condL2, freqL2, sessL2] = condFreqSessLabels(T, G2);
for k = 1:numel(plv_meds)
    fprintf('  %-6s (recording label %3d Hz) %-4s : median PLV = %.3f, IQR = [%.3f, %.3f], n = %d\n', ...
        condL2(k), freqL2(k), sessL2(k), ...
        plv_meds(k), plv_iqr_lo(k), plv_iqr_hi(k), n_obs(k));
end


%% (A) Control band on 333 Hz data: sit vs firm vs foam
if cfg.run.A_within_333
    fprintf('\n=== (A) Control band %s on 333 Hz recordings: sit vs firm vs foam ===\n', ...
            cfg.control_band_label_str);
    Ta = T(T.Frequency==333 & T.Session=="pre" & ...
           ismember(T.Condition, {'sit','firm','foam'}), :);
    Ta = dropUnused(Ta, {'Condition','Subject','Channel'});
    if cfg.restrict_to_sigchans
        Ta.Condition = categorical(Ta.Condition);
        Ta = selectSigChans(Ta, 'Condition', cfg.sig_combine);
    end

    for refIdx = 1:numel(cfg.condition_refs)
        refLvl = cfg.condition_refs{refIdx};
        tag    = sprintf('A_ctrlband_333_ref-%s', refLvl);
        fprintf('  -- reference = %s --\n', refLvl);

        Ta_r = Ta;
        Ta_r.Condition = reorderCategoricalRef(Ta_r.Condition, refLvl, ...
                                               {'sit','firm','foam'});

        if checkFactor(Ta_r, 'Condition', 2, tag)
            chanRE = '';
            if cfg.channel_random_effect, chanRE = ' + (1|Channel)'; end
            formula = [cfg.outcome ' ~ Condition + (1|Subject)' chanRE];
            res_A   = fitLMM(Ta_r, formula, tag, cfg);
            res_A.cohens_d = cohensD(Ta_r, 'Condition', cfg.outcome);
            writeResults(res_A, tag, cfg);
            if cfg.run_loso
                leaveOneSubjectOut(res_A, tag, cfg);
            end
            if cfg.checkAssumptions
                checkAssumptions(res_A, Ta_r, 'Condition', tag, cfg);
            end
            if cfg.make_plots
                plotResults(Ta_r, res_A, 'Condition', tag, ...
                    sprintf('CONTROL BAND %s (on 333 Hz recordings): sit vs firm vs foam (ref = %s)', ...
                            cfg.control_band_label_str, refLvl), cfg);
            end
        end
    end
end


%% =================== (B) Control band on 333 Hz data: pre vs post (sit vs sit_2) ===================
if cfg.run.B_prepost_333
    fprintf('\n=== (B) Control band %s on 333 Hz recordings: sit (pre) vs sit_2 (post) ===\n', ...
            cfg.control_band_label_str);
    Tb = T(T.Frequency==333 & T.Condition=="sit", :);
    Tb.Session = categorical(Tb.Session, {'pre','post'});
    Tb = dropUnused(Tb, {'Session','Subject','Channel'});

    Tb = keepPaired(Tb, 'Session');
    if cfg.restrict_to_sigchans
        Tb = selectSigChans(Tb, 'Session', cfg.sig_combine);
    end

    tag = 'B_ctrlband_333_prepost';
    if checkFactor(Tb, 'Session', 2, tag)
        chanRE = '';
        if cfg.channel_random_effect, chanRE = ' + (1|Channel)'; end
        formula = [cfg.outcome ' ~ Session + (1|Subject)' chanRE];
        res_B   = fitLMM(Tb, formula, tag, cfg);
        res_B.cohens_d = cohensD(Tb, 'Session', cfg.outcome);
        writeResults(res_B, tag, cfg);
        if cfg.run_loso
            leaveOneSubjectOut(res_B, tag, cfg);
        end
        if cfg.checkAssumptions
            checkAssumptions(res_B, Tb, 'Session', tag, cfg);
        end
        if cfg.make_plots
            plotResults(Tb, res_B, 'Session', tag, ...
                sprintf('CONTROL BAND %s (on 333 Hz recordings): pre (sit) vs post (sit\\_2)', ...
                        cfg.control_band_label_str), cfg);
        end
    end
end


%% =================== (C) Control band on 89 Hz data: sit vs firm vs foam ===================
if cfg.run.C_within_89
    fprintf('\n=== (C) Control band %s on 89 Hz recordings: sit vs firm vs foam ===\n', ...
            cfg.control_band_label_str);
    Tc = T(T.Frequency==89 & T.Session=="pre" & ...
           ismember(T.Condition, {'sit','firm','foam'}), :);
    Tc = dropUnused(Tc, {'Condition','Subject','Channel'});
    if cfg.restrict_to_sigchans
        Tc.Condition = categorical(Tc.Condition);
        Tc = selectSigChans(Tc, 'Condition', cfg.sig_combine);
    end

    for refIdx = 1:numel(cfg.condition_refs)
        refLvl = cfg.condition_refs{refIdx};
        tag    = sprintf('C_ctrlband_89_ref-%s', refLvl);
        fprintf('  -- reference = %s --\n', refLvl);

        Tc_r = Tc;
        Tc_r.Condition = reorderCategoricalRef(Tc_r.Condition, refLvl, ...
                                               {'sit','firm','foam'});

        if checkFactor(Tc_r, 'Condition', 2, tag)
            chanRE = '';
            if cfg.channel_random_effect, chanRE = ' + (1|Channel)'; end
            formula = [cfg.outcome ' ~ Condition + (1|Subject)' chanRE];
            res_C   = fitLMM(Tc_r, formula, tag, cfg);
            res_C.cohens_d = cohensD(Tc_r, 'Condition', cfg.outcome);
            writeResults(res_C, tag, cfg);
            if cfg.run_loso
                leaveOneSubjectOut(res_C, tag, cfg);
            end
            if cfg.checkAssumptions
                checkAssumptions(res_C, Tc_r, 'Condition', tag, cfg);
            end
            if cfg.make_plots
                plotResults(Tc_r, res_C, 'Condition', tag, ...
                    sprintf('CONTROL BAND %s (on 89 Hz recordings): sit vs firm vs foam (ref = %s)', ...
                            cfg.control_band_label_str, refLvl), cfg);
            end
        end
    end
end


fprintf('\nAll selected CONTROL-BAND PLV analyses complete. Results in:\n  %s\n', cfg.logfile);


%% =================== HELPER FUNCTIONS ===================
% These are identical to MLM.m -- copied here so the script is standalone.

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

if isfield(cfg, 'remove_outliers') && cfg.remove_outliers
    [keepMask, flagInfo] = flagOutliers(res.full.lme, Tin, cfg);
    res.outlier_info = flagInfo;

    n_flag = sum(~keepMask);
    fprintf('  Outlier removal (%s): flagged %d / %d obs (%.1f%%)\n', ...
            flagInfo.method, n_flag, numel(keepMask), 100*n_flag/numel(keepMask));

    if n_flag == 0
        res.trimmed = res.full;
    else
        Ttrim = Tin(keepMask, :);
        Ttrim = dropUnused(Ttrim, Ttrim.Properties.VariableNames);
        try
            res.trimmed = fitOneLMM(Ttrim, formula);
            printLMM(res.trimmed, 'outliers removed');
        catch ME
            warning('PLV:RefitFailed', ...
                    'Refit after outlier removal failed: %s', ME.message);
            res.trimmed = res.full;
        end
    end
else
    res.trimmed = res.full;
end

main = res.trimmed;
res.lme           = main.lme;
res.fixedEffects  = main.fixedEffects;
res.R2_marginal   = main.R2_marginal;
res.R2_conditional= main.R2_conditional;
res.n_obs         = main.n_obs;
res.n_subj        = main.n_subj;
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
                          d, ci_lo, ci_hi, numel(xa) };
    end
end
d_tbl = cell2table(rows, 'VariableNames', ...
    {'Level_A','Level_B','Mean_A','Mean_B','Cohens_d_paired', ...
     'CI_lower','CI_upper','N_pairs'});
disp(d_tbl);
end


function writeResults(res, tag, cfg)
file = cfg.logfile;

FE = dataset2table(res.fixedEffects);
writetable(FE, file, 'Sheet', sheetName(tag,'FE'), 'WriteMode','overwritesheet');

hasTrim = isfield(res,'full') && isfield(res,'trimmed') && ...
          ~isequal(res.full.n_obs, res.trimmed.n_obs);
if hasTrim
    S = table( ...
        ["full"; "trimmed"], ...
        [res.full.n_obs;     res.trimmed.n_obs], ...
        [res.full.n_subj;    res.trimmed.n_subj], ...
        [res.full.R2_marginal;    res.trimmed.R2_marginal], ...
        [res.full.R2_conditional; res.trimmed.R2_conditional], ...
        'VariableNames', {'Fit','N_obs','N_subjects','R2_marginal','R2_conditional'});
else
    S = table( string(res.tag), string(res.formula), res.n_obs, res.n_subj, ...
               res.R2_marginal, res.R2_conditional, ...
               'VariableNames', {'Tag','Formula','N_obs','N_subjects', ...
                                 'R2_marginal','R2_conditional'});
end
writetable(S, file, 'Sheet', sheetName(tag,'summary'), 'WriteMode','overwritesheet');

if hasTrim
    FE_full = dataset2table(res.full.fixedEffects);
    writetable(FE_full, file, 'Sheet', sheetName(tag,'FE_full'), ...
               'WriteMode','overwritesheet');

    if isfield(res, 'outlier_info') && ~isempty(res.outlier_info.flagged_table)
        writetable(res.outlier_info.flagged_table, file, ...
                   'Sheet', sheetName(tag,'outliers'), 'WriteMode','overwritesheet');
    end
end

if isfield(res,'cohens_d') && ~isempty(res.cohens_d)
    writetable(res.cohens_d, file, 'Sheet', sheetName(tag,'d'), 'WriteMode','overwritesheet');
end

fprintf('Wrote results for %s to %s\n', tag, file);
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

outcomeVar = cfg.outcome;
switch outcomeVar
    case 'PLV',       outcomeLabel = 'PLV (control band)';
    case 'PLV_logit', outcomeLabel = 'logit(PLV) (control band)';
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

fig = figure('Name', ['PLV LMM control band: ' tag], ...
             'Position', [80 80 1500 900], 'Color','w');
tiledlayout(fig, 2, 2, 'Padding','compact', 'TileSpacing','compact');

% Panel 1: paired-subject plot
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
        pts(end+1).y_data = M(i, lastCol);
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

% Panel 2: per-channel mean PLV (control band) by condition
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

% Panel 3: forest plot
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

% Panel 4: variance decomposition
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
    case 'PLV',       unitLabel = 'PLV (control band)';
    case 'PLV_logit', unitLabel = 'logit (control band)';
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

fig = figure('Name', ['PLV LMM control-band diagnostics: ' tag], ...
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
sgtitle(sprintf('PLV LMM CONTROL-BAND assumption checks: %s   |   %s%s', tag, vifstr, verdict), ...
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


function [condL, freqL, sessL] = condFreqSessLabels(T, G)
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

fig = figure('Name', ['Outliers control band: ' tag], ...
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

sgtitle(sprintf('Outlier diagnostics (control band): %s', tag), 'FontWeight','bold', 'FontSize',12)

outpng = fullfile(cfg.figpath, [tag '_outliers.png']);
exportgraphics(fig, outpng, 'Resolution', 200);
fprintf('Saved outlier figure: %s\n', outpng);
end


function loso_tbl = leaveOneSubjectOut(res, tag, cfg)
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

subs  = unique(Tin.Subject);
n_sub = numel(subs);
if n_sub < 3
    fprintf('  [LOSO %s] only %d subject(s); LOSO is not meaningful.\n', tag, n_sub);
    return
end

FE_ref    = res.fixedEffects;
isInt     = strcmp(string(FE_ref.Name), '(Intercept)');
coefNames = string(FE_ref.Name(~isInt));
nCoef     = numel(coefNames);
if nCoef == 0
    fprintf('  [LOSO %s] no non-intercept coefficients to track.\n', tag);
    return
end

ref_est = FE_ref.Estimate(~isInt);
ref_p   = FE_ref.pValue(~isInt);

fprintf('--- Leave-one-subject-out for %s (n=%d subjects, %d coef) ---\n', ...
        tag, n_sub, nCoef);

est_mat   = nan(n_sub, nCoef);
p_mat     = nan(n_sub, nCoef);
ci_lo_mat = nan(n_sub, nCoef);
ci_hi_mat = nan(n_sub, nCoef);
fit_ok    = false(n_sub, 1);

for i = 1:n_sub
    Tsub = Tin(Tin.Subject ~= subs(i), :);
    Tsub = dropUnused(Tsub, Tsub.Properties.VariableNames);
    try
        lme_i = fitlme(Tsub, formula, 'FitMethod','REML');
        [~,~,FE_i] = fixedEffects(lme_i, 'DFMethod','satterthwaite');
        for c = 1:nCoef
            idx = strcmp(string(FE_i.Name), coefNames(c));
            if any(idx)
                est_mat(i,c)   = FE_i.Estimate(idx);
                p_mat(i,c)     = FE_i.pValue(idx);
                ci_lo_mat(i,c) = FE_i.Lower(idx);
                ci_hi_mat(i,c) = FE_i.Upper(idx);
            end
        end
        fit_ok(i) = true;
    catch ME
        fprintf('  LOSO refit failed excluding %s: %s\n', ...
                char(string(subs(i))), ME.message);
    end
end

for c = 1:nCoef
    e  = est_mat(fit_ok, c);
    pv = p_mat(fit_ok, c);
    if isempty(e), continue; end
    fprintf('  %-22s : beta [%+.4g, %+.4g] (ref %+.4g),  p [%.3g, %.3g] (ref %.3g)\n', ...
        coefNames(c), min(e), max(e), ref_est(c), min(pv), max(pv), ref_p(c));

    if ref_p(c) < 0.05 && any(pv >= 0.05)
        flippers = subs(fit_ok & p_mat(:,c) >= 0.05);
        fprintf('    *** Significance LOST when excluding: %s\n', ...
                strjoin(string(flippers), ', '));
    elseif ref_p(c) >= 0.05 && any(pv < 0.05)
        flippers = subs(fit_ok & p_mat(:,c) < 0.05);
        fprintf('    Significance gained when excluding: %s\n', ...
                strjoin(string(flippers), ', '));
    end
end

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

if ~(isfield(cfg,'make_plots') && cfg.make_plots), return; end

panel_w = max(380, 220 * nCoef);
fig = figure('Name', ['LOSO control band: ' tag], ...
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
                 est_mat(i,c) - ci_lo_mat(i,c), ...
                 ci_hi_mat(i,c) - est_mat(i,c), ...
                 'horizontal', 'o', 'Color', col, ...
                 'MarkerFaceColor', col, 'MarkerSize', 5, ...
                 'LineWidth', 0.9, 'CapSize', 4);
    end
    xline(ref_est(c), 'b-', 'LineWidth', 1.5);
    xline(0, 'k:', 'LineWidth', 0.8);

    yticks(ypos);
    yticklabels(cellstr(string(subs)));
    ylim([0.5 n_sub+0.5])
    set(gca, 'FontSize', 7)
    xlabel('Estimate (95% CI)')

    n_lost = sum(fit_ok & ref_p(c) < 0.05 & p_mat(:,c) >= 0.05);
    title(sprintf('%s\nref beta=%.3g, p=%.3g\n%d / %d exclusions -> p>=0.05', ...
          coefNames(c), ref_est(c), ref_p(c), n_lost, sum(fit_ok)), ...
          'Interpreter','none', 'FontSize',8)
end
sgtitle(sprintf('LOSO sensitivity (control band): %s   |   blue = full-fit estimate, orange = p>=0.05 when subject excluded', tag), ...
        'FontWeight','bold', 'FontSize',11)

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