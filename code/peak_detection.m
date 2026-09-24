%% FFR per-epoch first-peak latency (Cz, mastoid, single trials)
%  Loads preprocessed single-trial files, detects the first peak after
%  stimulus onset in each epoch (regardless of polarity), records every
%  epoch's latency to Excel, and reports mean/median across all epochs.
%
%  File naming follows the preprocessing script:
%     preprocessed_P[SUBJ]_FFR_<cond>_<freq>(_2).mat   (variable 'data_final')
clc; clear; close all;

%% Settings
basepath  = 'F:\Auditory-mov-data';
subfolder = 'preprocessed';
subjects  = {'NG', 'BVD', 'FD', 'EF', 'CX', 'HX', 'MCBG', 'LV',  ...
                 'AM', 'FS', 'LVB', 'RNR', 'VUR', 'LH', 'TF', 'MGZ', ...
                 'NSMT', 'MM', 'TH', 'TD'};

% Which data to include
freqs      = [89 333];                 % 89, 333, or both [89 333]
conditions = {'sit','firm','foam'};
chanName  = 'Cz';                      % channel to read peaks from

% Peak-search options
searchMin = 0.030;                       % start search at stimulus onset (in sec)
searchMax = 0.100;                        % [] = whole post-onset epoch, or e.g. 0.050
minProm   = true;                        % [] = no constraint; else prominence (data units)

outpath = 'F:\Auditory-mov-data\FFR_epoch_peak_latencies.xlsx';

% build file names
patterns = {};
for ff = freqs
    for cc = 1:numel(conditions)
        patterns{end+1} = sprintf('preprocessed_P[*]_FFR_%s_%d.mat',conditions{cc}, ff); 
    end
end

%% Collect per-epoch results
rows = {};   % Subject, File, Frequency, Condition, Restart, TrialIdx, CondLabel, Latency_ms, Amp, Sign

for s = 1:numel(subjects)
    subj   = subjects{s};
    folder = fullfile(basepath, subj, subfolder);

    % gather files for all requested freq/condition patterns
    files = [];

    for p = 1:numel(patterns)
        files = [files; dir(fullfile(folder, patterns{p}))]; 
    end

    if isempty(files)
        fprintf('No matching files for %s in %s\n', subj, folder); continue;
    end

    [~, ui] = unique({files.name}); files = files(ui);   % dedupe

    for f = 1:numel(files)
        nm = files(f).name;

        % --- parse condition / frequency from filename ---
        tokC = regexp(nm, '_FFR_(sit|firm|foam)_', 'tokens', 'once');
        condStr = ""; if ~isempty(tokC), condStr = string(tokC{1}); end
        tokF = regexp(nm, '_(89|333)\.mat$', 'tokens', 'once');
        freqVal = NaN; if ~isempty(tokF), freqVal = str2double(tokF{1}); end

        % --- load single-trial data ---
        S = load(fullfile(folder, nm));
        D = S.data_final;

        ci = find(strcmp(D.label, chanName), 1);
        if isempty(ci)
            fprintf('%s: channel %s not found\n', nm, chanName); continue;
        end

        % --- per-epoch first-peak detection ---
        ntr = numel(D.trial);
        for k = 1:ntr
            t = D.time{k}(:)';
            y = D.trial{k}(ci,:);

            if isempty(searchMax), mask = t >= searchMin;
            else,                  mask = t >= searchMin & t <= searchMax; end
            tt = t(mask); yy = y(mask);
            if isempty(yy), continue; end

            % first peak regardless of sign: detect on +yy and -yy, take earliest
            if ~ (minProm)
                [~, lPos] = findpeaks(yy);
                [~, lNeg] = findpeaks(-yy);
            else
                % MinPeakProminence threshold: adjust based on your data amplitude
                min_prom = 0.8*max(abs(yy));   % 80% of max amplitude

                [~, lPos] = findpeaks(yy,  'MinPeakProminence', min_prom);
                [~, lNeg] = findpeaks(-yy, 'MinPeakProminence', min_prom);
            end

            cand = [];
            if ~isempty(lPos), cand = [cand; lPos(1) +1]; end 
            if ~isempty(lNeg), cand = [cand; lNeg(1) -1]; end
            if isempty(cand), continue; end

            [idx, ~] = min(cand(:,1));            % earliest sample index wins
            sgn = cand(cand(:,1)==idx, 2); sgn = sgn(1);

            lat = tt(idx) * 1000;                 % ms
            amp = yy(idx);

            % condition label (cell-string or numeric trialinfo)
            if iscell(D.trialinfo)
                clab = string(D.trialinfo{k});
            else
                clab = string(D.trialinfo(k));
            end

            rows = [rows; {string(subj), string(nm), freqVal, condStr, ...
                           k, clab, lat, amp, sgn}]; 
        end

        fprintf('%-6s %-45s  freq=%g cond=%s  %d epochs\n', ...
            subj, nm, freqVal, condStr, ntr);
    end
end

if isempty(rows)
    error('No epochs collected. Check paths, file patterns, and channel name.');
end


%% Build table & write to Excel
T = cell2table(rows, 'VariableNames', ...
    {'Subject','File','Frequency','Condition', ...
     'TrialIdx','CondLabel','Latency_ms','Amp','Sign'});

writetable(T, outpath, 'Sheet', 'Epochs', 'WriteMode','overwritesheet');


%% Summary: separate per frequency (89 vs 333)
fprintf('\n================ Summary (%s, %s) ================\n', chanName, subfolder);

uFreq = unique(T.Frequency);
freqRows = {};
for q = 1:numel(uFreq)
    fv  = uFreq(q);
    sub = T(T.Frequency == fv, :);
    lat = sub.Latency_ms;
    nP  = sum(sub.Sign ==  1);
    nN  = sum(sub.Sign == -1);

    fprintf('\n-- %g Hz --\n', fv);
    fprintf('  Epochs with a peak: %d\n', height(sub));
    fprintf('  Positive: %d   Negative: %d\n', nP, nN);
    fprintf('  Mean   onset latency: %.3f ms\n', mean(lat,'omitnan'));
    fprintf('  Median onset latency: %.3f ms\n', median(lat,'omitnan'));

    freqRows = [freqRows; {fv, mean(lat,'omitnan'), median(lat,'omitnan'), ...
                           std(lat,'omitnan'), height(sub), nP, nN}]; 
end

summary = cell2table(freqRows, 'VariableNames', ...
    {'Frequency','MeanLatency_ms','MedianLatency_ms','SD_ms', ...
     'nEpochs','nPositive','nNegative'});
writetable(summary, outpath, 'Sheet', 'Summary', 'WriteMode','overwritesheet');

% Per frequency x condition breakdown
byGroup = groupsummary(T, {'Frequency','Condition'}, {'mean','median','std'}, 'Latency_ms');
writetable(byGroup, outpath, 'Sheet', 'ByGroup', 'WriteMode','overwritesheet');

fprintf('\nWritten to %s (sheets: Epochs, Summary, ByGroup)\n', outpath);
