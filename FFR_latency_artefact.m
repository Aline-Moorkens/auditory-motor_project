%% FFR latency, script 1 of 2: neural response vs artefact (onset/offset)
%
%  TIME AXIS. The epochs run from -10 to 240 ms re the TRIGGER. The tone
%  starts 30 ms after the trigger (all frequencies), so the script shifts the
%  time axis to t = 0 at tone onset (toneOnset = 0.030):
%       trigger = -30 ms | tone onset = 0 | tone offset = 200 ms | epoch end = 210 ms
%  Onset latencies in the output are re tone onset (add 30 ms for re trigger);
%  offset latencies are re tone offset (230 ms re trigger).
%
%  LOGIC. The difference waveform (pos - neg)/2, i.e. the FFR as defined in
%  the paper, also contains the cochlear microphonic (CM) and any electro-
%  magnetic artefact. These start at the electrical/acoustic arrival of the
%  tone and stop at its offset. A neural response starts several ms later and
%  outlasts the tone by about the same lag. Per recording (subject x frequency
%  x condition), on the polarity-combined average:
%   - OnsetNoise : the envelope of the difference waveform rises above the
%                  baseline noise (mean + kSD x SD) and stays above for >= minAbove
%   - OnsetSum   : the same on the sum waveform (CM and artefact cancel), as a
%                  neural onset reference
%   - Offset     : time after tone offset at which the envelope drops back
%                  below threshold (NaN with OffsetCensored = true means it is
%                  still above at the end of the epoch, 10 ms after offset)
%
%  FILTER (AND RAMP) CORRECTION. The zero-phase Butterworth band-pass used in
%  preprocessing has no delay, but it smears the response both ways in time:
%  onsets look earlier and offsets later. The 5 ms raised-cosine ramps also
%  delay the point where any threshold is crossed. Both effects are removed
%  with a zero-lag MODEL RESPONSE: the half-wave rectified stimulus train
%  (alternating polarity, same SOA and ramps), filtered with the same filter,
%  epoched and averaged like the data. Its difference waveform is the
%  filtered tone, and its sum is the filtered rectified envelope.
%  The model is measured at the SAME threshold level as the data: threshold
%  as a fraction of each recording's plateau envelope. Then:
%       corrected latency = data latency - model latency
%  This is the lag of the response relative to an identically filtered copy
%  of the stimulus.
%  It assumes that the response envelope has the shape of the stimulus
%  envelope. Neural onset responses that overshoot or adapt break that
%  assumption a little. The model assumes the filter ran on the continuous
%  recording (before epoching), as in the methods notes.
%  Hardware delays come on top: fill in earDelay_ms / elecDelay_ms from the
%  control recordings. Then run the PHANTOM recordings through this script
%  (phantomFiles): they give artefact-only values under identical processing.
%
%  File naming follows the preprocessing script:
%     preprocessed_P<SUBJ>_FFR_<cond>_<freq>(_2).mat   (variable 'data_final')
%  Requires the Signal Processing Toolbox (butter, filtfilt, hilbert).
clc; clear; close all;
rng(1);

%% Settings
basepath  = 'F:\Auditory-mov-project\master-thesis';
subfolder = 'preprocessed';
subjects  = {'NG', 'BVD', 'FD', 'EF', 'CX', 'HX', 'MCBG', 'LV',  ...
             'AM', 'FS', 'LVB', 'RNR', 'VUR', 'LH', 'TF', 'MGZ', ...
             'NSMT', 'MM', 'TH', 'TD'};

freqs           = [89 333];               % add 433 for the 333-433 protocol version
conditions      = {'sit','firm','foam'};
includeRestarts = true;                   % *_2.mat epochs are pooled with the first file

% Phantom control recordings, preprocessed like the data (full paths). The
% names must contain _FFR_<cond>_<freq>.mat or _FFR_<cond>_<freq>_2.mat.
phantomFiles = {};
% phantomFiles = {'F:\Auditory-mov-data\phantom\preprocessed\preprocessed_Pphantom_FFR_sit_333.mat'};

chanNames = {'Cz'};                       % or the a priori ROI of the paper:
% chanNames = {'Fz','F1','F2','FC1','FCz','FC2','C1','Cz','C2'};

% Timing (s)
toneOnset    = 0.030;                     % trigger -> tone onset. Epoch time is re TRIGGER.
toneDur      = 0.200;
rampDur      = 0.005;                     % raised-cosine on/off ramps (methods notes)
% Windows below: re TONE ONSET (offsetSearch: re tone offset)
baseWin      = [-0.030 -0.005];           % noise window (= 0-25 ms re trigger)
onsetSearch  = [-0.005 0.040];
steadyWin    = [0.020 0.180];             % plateau window (= 50-210 ms re trigger)
offsetSearch = [0 0.010];                 % the epoch ends 10 ms after tone offset

% Preprocessing filter and stimulus period (SOA) per frequency. These MUST
% match the preprocessing (FieldTrip bpfilttype 'but', bpfreq, bpfiltord).
%            freq   lo    hi   order  SOA (s)
filtTab = [   89    60   190    4    0.650
             333   200   900    6    0.240
             433   200   900    6    0.240 ];
filtDir = 'twopass';                      % 'twopass' = zero-phase (filtfilt), as in preprocessing;
                                          % 'onepass' if the data were filtered causally

% Measured delays from the control recordings, in ms re the NOMINAL tone
% onset (measured trigger-to-transducer or trigger-to-ear delay minus 30 ms).
% Leave NaN until measured.
elecDelay_ms = NaN;                       % electrical arrival (electromagnetic artefact)
earDelay_ms  = NaN;                       % acoustic arrival at the ER-3A tip (CM)

% Epochs per polarity: [] = all (the two polarities are always equalised
% within a recording); or a number to cap every recording at the same count
nPerPolarity = [];

% Detection
kSD      = 3;                             % threshold = baseline envelope mean + kSD x SD
minAbove = 0.002;                         % a crossing must last >= 2 ms

qcPlot  = true;
outpath = 'F:\Auditory-mov-data\FFR_latency_artefact.xlsx';

P = struct('baseWin', baseWin, 'onsetSearch', onsetSearch, 'steadyWin', steadyWin, ...
           'offsetSearch', offsetSearch, 'toneDur', toneDur, 'rampDur', rampDur, ...
           'kSD', kSD, 'minAbove', minAbove, 'filtDir', filtDir);

%% Files -> recordings
F = listFiles(basepath, subfolder, subjects, freqs, conditions, includeRestarts);
F.IsPhantom = false(height(F), 1);
for i = 1:numel(phantomFiles)
    [fo, nm, ex] = fileparts(phantomFiles{i});
    tok = regexp([nm ex], '_FFR_(sit|firm|foam)_(\d+)(?:_2)?\.mat$', 'tokens', 'once');
    if isempty(tok), fprintf('Phantom file name not recognised: %s\n', phantomFiles{i}); continue; end
    F = [F; table("PHANTOM_" + string(nm), string(fo), string([nm ex]), str2double(tok{2}), ...
                  string(tok{1}), false, true, 'VariableNames', F.Properties.VariableNames)]; %#ok<AGROW>
end
if isempty(F), error('No files found. Check paths and file patterns.'); end
[G, gSubj, gFreq, gCond] = findgroups(F.Subject, F.Frequency, F.Condition);
gPhantom = splitapply(@any, F.IsPhantom, G);

%% Measures per recording
modelCache = cell(size(filtTab, 1), 1);   % zero-lag model per frequency (same time axis)
rows    = cell(max(G), 1);
curSubj = "";
for g = 1:max(G)
    [Y, pol, t] = loadGroup(F(G == g, :), chanNames, toneOnset);
    if isempty(Y), continue; end
    f  = gFreq(g);
    fi = find(filtTab(:, 1) == f, 1);
    if isempty(fi), fprintf('%d Hz not in filtTab, skipped\n', f); continue; end

    % polarity averages with equal epochs per polarity
    iP = find(pol == 1); iN = find(pol == -1);
    n  = min(numel(iP), numel(iN));
    if ~isempty(nPerPolarity), n = min(n, nPerPolarity); end
    iP = iP(randperm(numel(iP), n)); iN = iN(randperm(numel(iN), n));
    mP = mean(Y(iP, :), 1); mN = mean(Y(iN, :), 1);
    yD = (mP - mN) / 2;                   % FFR + CM + artefact
    yS = (mP + mN) / 2;                   % CM and artefact cancel

    % zero-lag model, filtered like the data (computed once per frequency)
    if isempty(modelCache{fi})
        [mD, mS] = modelResponse(t, f, filtTab(fi, :), P);
        modelCache{fi} = struct('mD', mD, 'mS', mS);
        reportModel(t, f, mD, P);
    end
    mdl = modelCache{fi};

    [M, Q] = measureArtefact(t, yD, yS, mdl.mD, mdl.mS, P);
    rows{g} = [table(gSubj(g), gPhantom(g), f, gCond(g), n, ...
                     'VariableNames', {'Subject','IsPhantom','Frequency','Condition','nPerPolarity'}), ...
               struct2table(M)];

    fprintf('%-14s %4d Hz %-5s | onset raw %6.2f corr %6.2f | sum corr %6.2f | offset raw %5.2f corr %5.2f ms\n', ...
            gSubj(g), f, gCond(g), M.OnsetNoise_ms, M.OnsetCorr_ms, M.OnsetSumCorr_ms, ...
            M.Offset_ms, M.OffsetCorr_ms);

    % --- QC plot ---
    if qcPlot
        if gSubj(g) ~= curSubj
            curSubj = gSubj(g);
            fig = figure('Name', curSubj, 'Color', 'w');
            tl  = tiledlayout(fig, 'flow', 'TileSpacing', 'compact');
            title(tl, {sprintf('%s | %s', curSubj, strjoin(chanNames, '+')), ...
                'grey: difference waveform | black: envelope | green: zero-lag filtered model (scaled) | red: onset/offset'}, ...
                'Interpreter', 'none', 'FontSize', 9);
        end
        ax = nexttile(tl); hold(ax, 'on');
        tm = t * 1000;
        plot(ax, tm, Q.yD, 'Color', [0.75 0.75 0.75]);
        plot(ax, tm, Q.envModel, 'Color', [0 0.6 0], 'LineWidth', 1);
        plot(ax, tm, Q.envD, 'k');
        yline(ax, Q.thr, '--r');
        xline(ax, [0 toneDur*1000], ':k');
        if ~isnan(earDelay_ms),     xline(ax, earDelay_ms, '--b'); end
        if ~isnan(M.OnsetNoise_ms), xline(ax, M.OnsetNoise_ms, 'r'); end
        if ~isnan(M.Offset_ms),     xline(ax, toneDur*1000 + M.Offset_ms, 'r'); end
        title(ax, sprintf('%s %d Hz | onset %.1f (corr %.1f) | off +%.1f (corr %.1f) ms', ...
              gCond(g), f, M.OnsetNoise_ms, M.OnsetCorr_ms, M.Offset_ms, M.OffsetCorr_ms), ...
              'Interpreter', 'none', 'FontSize', 8);
        xlim(ax, [tm(1) tm(end)]); xlabel(ax, 'ms re tone onset'); box(ax, 'off');
    end
end
R = vertcat(rows{:});
if isempty(R), error('No recordings could be processed.'); end

% relative to the measured hardware delays (NaN until filled in)
R.OnsetCorr_reEar_ms  = R.OnsetCorr_ms - earDelay_ms;
R.OnsetCorr_reElec_ms = R.OnsetCorr_ms - elecDelay_ms;

%% Summary per group (participants per frequency x condition; phantom separately)
sumRows = {};
grp = unique(R(:, {'IsPhantom','Frequency','Condition'}));
for k = 1:height(grp)
    r = R(R.IsPhantom == grp.IsPhantom(k) & R.Frequency == grp.Frequency(k) & ...
          R.Condition == grp.Condition(k), :);
    sumRows(end+1, :) = {grp.IsPhantom(k), grp.Frequency(k), grp.Condition(k), height(r), ...
        nnz(r.ResponsePresent), ...
        median(r.OnsetNoise_ms, 'omitnan'), median(r.ModelOnset_ms, 'omitnan'), ...
        median(r.OnsetCorr_ms, 'omitnan'), prctile(r.OnsetCorr_ms, 25), prctile(r.OnsetCorr_ms, 75), ...
        nnz(r.OnsetCensored), median(r.OnsetSumCorr_ms, 'omitnan'), ...
        median(r.OffsetCorr_ms, 'omitnan'), nnz(r.OffsetCensored), ...
        median(r.OnsetCorr_reEar_ms, 'omitnan')}; %#ok<SAGROW>
end
A = cell2table(sumRows, 'VariableNames', {'IsPhantom','Frequency','Condition','nRecordings', ...
    'nResponsePresent','OnsetRaw_median_ms','ModelOnset_median_ms', ...
    'OnsetCorr_median_ms','OnsetCorr_Q1_ms','OnsetCorr_Q3_ms','nOnsetCensored', ...
    'OnsetSumCorr_median_ms','OffsetCorr_median_ms','nOffsetCensored','OnsetCorr_reEar_median_ms'});

writetable(R, outpath, 'Sheet', 'Recordings', 'WriteMode', 'overwritesheet');
writetable(A, outpath, 'Sheet', 'Summary',    'WriteMode', 'overwritesheet');

fprintf('\n===== Neural vs artefact: corrected latencies (onset re tone onset, offset re tone offset) =====\n');
disp(A);
fprintf('Written to %s (sheets: Recordings, Summary)\n', outpath);




%% ======================================================================
function [M, Q] = measureArtefact(t, yD, yS, mD, mS, P)
% Data and zero-lag model measured at the same threshold level (fraction of plateau).
fs   = 1 / median(diff(t));
nMin = max(1, round(P.minAbove * fs));
b    = t >= P.baseWin(1)     & t <= P.baseWin(2);
st   = t >= P.steadyWin(1)   & t <= P.steadyWin(2);
on   = t >= P.onsetSearch(1) & t <= P.onsetSearch(2);
off  = t >= P.toneDur + P.offsetSearch(1) & t <= P.toneDur + P.offsetSearch(2);
tOn  = t(on);
tOff = t(off) - P.toneDur;

% --- difference waveform ---
yD    = yD - mean(yD(b));
envD  = envelopeOf(yD);
thr   = mean(envD(b)) + P.kSD * std(envD(b));
plat  = median(envD(st));
frac  = thr / plat;                        % threshold as a fraction of the plateau
envMD = envelopeOf(mD);
platM = median(envMD(st));

M.BaselineSD      = std(yD(b));
M.PlateauEnv      = plat;
M.ThrFrac         = frac;
M.ResponsePresent = frac < 1;

i0 = firstSustained(envD(on) > thr, nMin);
M.OnsetNoise_ms  = msOrNaN(tOn, i0);
M.OnsetCensored  = isequal(i0, 1);         % already above threshold at the window start
M.ModelOnset_ms  = msOrNaN(tOn, firstSustained(envMD(on) > frac * platM, nMin));
M.OnsetCorr_ms   = M.OnsetNoise_ms - M.ModelOnset_ms;

i0 = firstSustained(envD(off) < thr, nMin);
M.Offset_ms      = msOrNaN(tOff, i0);
M.OffsetCensored = isempty(i0);
M.ModelOffset_ms = msOrNaN(tOff, firstSustained(envMD(off) < frac * platM, nMin));
M.OffsetCorr_ms  = M.Offset_ms - M.ModelOffset_ms;

if ~M.ResponsePresent
    [M.ModelOnset_ms, M.OnsetCorr_ms, M.ModelOffset_ms, M.OffsetCorr_ms] = deal(NaN);
end

% --- sum waveform (CM and artefact cancel) ---
yS     = yS - mean(yS(b));
envS   = envelopeOf(yS);
thrS   = mean(envS(b)) + P.kSD * std(envS(b));
platS  = median(envS(st));
fracS  = thrS / platS;
envMS  = envelopeOf(mS);
platMS = median(envMS(st));

M.SumResponsePresent = fracS < 1;
M.OnsetSum_ms        = msOrNaN(tOn, firstSustained(envS(on) > thrS, nMin));
M.ModelOnsetSum_ms   = msOrNaN(tOn, firstSustained(envMS(on) > fracS * platMS, nMin));
M.OnsetSumCorr_ms    = M.OnsetSum_ms - M.ModelOnsetSum_ms;
if ~M.SumResponsePresent, [M.ModelOnsetSum_ms, M.OnsetSumCorr_ms] = deal(NaN); end

Q = struct('yD', yD, 'envD', envD, 'thr', thr, 'envModel', envMD * plat / platM);
end

function reportModel(t, f, mD, P)
% Print how far filter + ramps alone move the onset of a zero-lag response
fs   = 1 / median(diff(t));
nMin = max(1, round(P.minAbove * fs));
on   = t >= P.onsetSearch(1) & t <= P.onsetSearch(2);
st   = t >= P.steadyWin(1) & t <= P.steadyWin(2);
env  = envelopeOf(mD); pl = median(env(st)); tOn = t(on);
lv   = [0.1 0.25 0.5];
v    = arrayfun(@(x) msOrNaN(tOn, firstSustained(env(on) > x * pl, nMin)), lv);
fprintf('Zero-lag model %d Hz (%s): onset at 10/25/50%% of plateau = %.2f / %.2f / %.2f ms\n', ...
        f, P.filtDir, v);
end

function [mD, mS] = modelResponse(t, f, spec, P)
% Zero-lag model response, filtered like the data. Half-wave rectified
% alternating-polarity tone train (5 tones, same SOA, raised-cosine ramps),
% band-pass filtered on the "continuous" train, epoched at tones 3 (+) and 4 (-).
% Difference -> filtered tone; sum -> filtered rectified envelope.
fs  = 1 / median(diff(t));
pad = 0.5; nT = 5; soa = spec(5);
tt  = (0 : round((2*pad + nT*soa) * fs)) / fs - pad;
s   = zeros(size(tt));
for k = 0:nT-1
    tk = tt - k*soa;
    in = tk >= 0 & tk < P.toneDur;
    s(in) = s(in) + (-1)^k * sin(2*pi*f*tk(in)) .* rampEnv(tk(in), P.toneDur, P.rampDur);
end
r = max(s, 0);                                        % half-wave rectification, zero lag
[B, A] = butter(spec(4), spec(2:3) / (fs/2));         % as FieldTrip ft_preproc_bandpassfilter
if any(abs(roots(A)) >= 1)
    warning('Model filter %g-%g Hz (order %d) is unstable in [b,a] form.', spec(2), spec(3), spec(4));
end
if strcmpi(P.filtDir, 'twopass'), rf = filtfilt(B, A, r); else, rf = filter(B, A, r); end
e2 = interp1(tt, rf, 2*soa + t, 'linear', 0);
e3 = interp1(tt, rf, 3*soa + t, 'linear', 0);
mD = (e2 - e3) / 2;
mS = (e2 + e3) / 2;
end

function e = rampEnv(tk, dur, rampDur)
% Raised-cosine (cos^2) onset and offset ramps
e = ones(size(tk));
if rampDur <= 0, return; end
up = tk < rampDur;        e(up) = sin(pi/2 * tk(up) / rampDur).^2;
dn = tk > dur - rampDur;  e(dn) = sin(pi/2 * (dur - tk(dn)) / rampDur).^2;
end

function F = listFiles(basepath, subfolder, subjects, freqs, conditions, includeRestarts)
% One row per file: Subject, Folder, Name, Frequency, Condition, Restart
Subject = strings(0,1); Folder = strings(0,1); Name = strings(0,1);
Frequency = zeros(0,1); Condition = strings(0,1); Restart = false(0,1);
for s = 1:numel(subjects)
    folder = fullfile(basepath, subjects{s}, subfolder);
    nFound = 0;
    for ff = freqs
        for cc = 1:numel(conditions)
            d = dir(fullfile(folder, sprintf('preprocessed_P*_FFR_%s_%d*.mat', conditions{cc}, ff)));
            for k = 1:numel(d)
                nm  = d(k).name;
                tok = regexp(nm, '_FFR_(sit|firm|foam)_(\d+)(?:_2)?\.mat$', 'tokens', 'once');
                if isempty(tok) || str2double(tok{2}) ~= ff, continue; end
                isR = ~isempty(regexp(nm, '_2\.mat$', 'once'));
                if isR && ~includeRestarts, continue; end
                Subject(end+1,1) = string(subjects{s}); Folder(end+1,1) = string(folder); %#ok<AGROW>
                Name(end+1,1) = string(nm); Frequency(end+1,1) = ff;                      %#ok<AGROW>
                Condition(end+1,1) = string(conditions{cc}); Restart(end+1,1) = isR;      %#ok<AGROW>
                nFound = nFound + 1;
            end
        end
    end
    if nFound == 0, fprintf('No matching files for %s in %s\n', subjects{s}, folder); end
end
F = table(Subject, Folder, Name, Frequency, Condition, Restart);
end

function [Y, pol, t] = loadGroup(Fg, chanNames, toneOnset)
% Load all files of one recording; epochs averaged over chanNames.
% t is re TONE ONSET (epoch time minus toneOnset).
Y = []; pol = []; t = [];
for i = 1:height(Fg)
    nm = Fg.Name(i);
    S  = load(fullfile(Fg.Folder(i), nm), 'data_final');
    D  = S.data_final;
    ci = find(ismember(D.label, chanNames));
    if numel(ci) ~= numel(chanNames)
        fprintf('%s: channel(s) %s not found, skipped\n', nm, strjoin(setdiff(chanNames, D.label), ', '));
        continue;
    end
    nsamp = cellfun(@(x) size(x, 2), D.trial);
    if any(nsamp ~= nsamp(1)), fprintf('%s: epochs differ in length, skipped\n', nm); continue; end
    Yi = cell2mat(cellfun(@(x) mean(x(ci, :), 1), D.trial(:), 'UniformOutput', false));
    p  = getPolarity(D, size(Yi, 1));
    if any(isnan(p)), fprintf('%s: no 2-level polarity label in trialinfo, skipped\n', nm); continue; end
    ti = D.time{1}(:)' - toneOnset;
    if isempty(t)
        t = ti;
    elseif numel(ti) ~= numel(t) || max(abs(ti - t)) > 1e-9
        fprintf('%s: time axis differs from the first file, skipped\n', nm); continue;
    end
    Y = [Y; Yi]; pol = [pol; p]; %#ok<AGROW>
end
end

function pol = getPolarity(D, ntr)
% +1 / -1 per epoch from trialinfo. The lower label after sorting becomes -1
% (correct for -1/+1 coding; otherwise the sign of the difference waveform is
% arbitrary, which does not affect latencies).
pol = nan(ntr, 1);
if ~isfield(D, 'trialinfo') || isempty(D.trialinfo), return; end
if iscell(D.trialinfo), lab = string(D.trialinfo(:));
else,                   lab = string(D.trialinfo(:, 1));
end
u = unique(lab);
if numel(u) ~= 2 || numel(lab) ~= ntr, return; end
pol(lab == u(1)) = -1;
pol(lab == u(2)) = +1;
end

function env = envelopeOf(y)
% Hilbert envelope with mirror padding, to limit edge effects at the epoch ends
y   = y(:)';
n   = numel(y);
e   = abs(hilbert([fliplr(y) y fliplr(y)]));
env = e(n+1 : 2*n);
end

function i0 = firstSustained(m, nMin)
% First index where m is true for nMin consecutive samples ([] if none)
c  = conv(double(m(:)'), ones(1, nMin), 'valid');
i0 = find(c == nMin, 1);
end

function v = msOrNaN(tt, i0)
if isempty(i0), v = NaN; else, v = tt(i0) * 1000; end
end
