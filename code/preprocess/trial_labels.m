function conditionlabels = trial_labels(trig, fs)

% TRIAL_DEF_MISSING_TRIG  Assign condition labels to the triggers of a
% 333 Hz FFR recording, correcting for missing and extra triggers.
%
%   conditionlabels = trial_def_missing_trig(trig, fs)
%
%   trig  vector with the sample of each trigger, e.g. trig = [evt.sample]
%   fs    sampling rate in Hz (10000)
%
%   conditionlabels  1 x n cell, one label per trigger:
%                    'FFR', '-FFR' or 'notrials' (START/MID/END markers
%                    and extra triggers)
%
%   Expected structure: START + 1500 trials (block 1) + MID + 1500 trials
%   (block 2) + END. Polarity alternates within each block, starting with FFR.
%
%   Handled automatically:
%   - missing triggers: the number of trials in an interval is counted, so
%     the polarity of the next trigger stays correct
%   - extra triggers before START or after END: labelled 'notrials'
%   - extra triggers within a block (interval that does not fit the 240 ms
%     trial grid): labelled 'notrials' and skipped, so the next interval is
%     measured from the last valid trigger
%   Not handled: EC sit_333_2, which is split over two files.


%% Timing (in samples)
trial_dur = 0.24 * fs;        % ~240 ms between consecutive triggers
tol_trial = 0.004 * fs;       % tolerance for trial intervals (4 ms; must stay
                              % below the 6.4 ms offset of known extra triggers)
tol_gap   = 0.15 * fs;        % tolerance for the gaps around the markers

% Boundary gap signatures
gap_start_mid  = 1.97 * fs;   % Gap after START or MID markers
gap_block1_end = 3.27 * fs;   % Gap after last trial of block 1
gap_block2_end = 2.27 * fs;   % Gap after last trial of block 2

n_per_block = 1500;

%% Initialise
n_trigs = length(trig);
conditionlabels = repmat({'notrials'}, 1, n_trigs);

fprintf('\n=== Assigning Condition Labels (1500-trial blocks) ===\n');
fprintf('Total triggers: %d (expected 3003)\n', n_trigs);

%% Find the START marker: the first trigger followed by the start gap
start_idx = find(abs(diff(trig) - gap_start_mid) < tol_gap, 1);
if isempty(start_idx)
    error('No START marker found: no gap of %.2f s between triggers.', gap_start_mid/fs);
end
if start_idx > 1
    fprintf('Trigger 1-%d: extra trigger(s) before START, labelled notrials\n', start_idx-1);
end
fprintf('Trigger %4d: START marker\n', start_idx);

%% State tracking
current_block = 1;
trial_count_in_block = 0;
block1_count = 0;
block2_count = 0;
last_trig = trig(start_idx);  % sample of the last valid trigger
max_dev = 0;                  % largest deviation of an accepted trial interval

%% Process all triggers after START
for i = start_idx+1 : n_trigs
    % Interval to the last VALID trigger, so extra triggers are skipped
    interval = trig(i) - last_trig;

    % Regular trial interval (possibly with missing triggers in between)?
    n_trials_in_interval = round(interval / trial_dur);
    is_trial_interval = n_trials_in_interval > 0 && ...
        abs(interval - n_trials_in_interval*trial_dur) < tol_trial;

    % Gap to the end marker of this block. Normally the block is complete;
    % if its last trigger(s) are missing, those trials are added to the gap.
    n_missing_at_end = max(0, n_per_block - trial_count_in_block);
    if current_block == 1
        end_gap = gap_block1_end + n_missing_at_end*trial_dur;
    else
        end_gap = gap_block2_end + n_missing_at_end*trial_dur;
    end
    is_end_gap = trial_count_in_block > 0 && abs(interval - end_gap) < tol_gap;

    if trial_count_in_block > 0 && is_trial_interval
        % Trial within a block: count all trials in this interval
        % (including missing ones) and label the recorded trigger
        trial_count_in_block = trial_count_in_block + n_trials_in_interval;
        max_dev = max(max_dev, abs(interval - n_trials_in_interval*trial_dur));

        % Polarity: odd trials in block -> FFR, even -> -FFR
        if mod(trial_count_in_block, 2) == 1
            conditionlabels{i} = 'FFR';
        else
            conditionlabels{i} = '-FFR';
        end

        if n_trials_in_interval > 1
            fprintf('Trigger %4d: %s (block %d, trial %4d, %d skipped)\n', ...
                    i, conditionlabels{i}, current_block, ...
                    trial_count_in_block, n_trials_in_interval-1);
        end

    elseif trial_count_in_block == 0 && abs(interval - gap_start_mid) < tol_gap
        % Gap after START or MID marker: first trial of a block (always FFR)
        trial_count_in_block = 1;
        conditionlabels{i} = 'FFR';
        fprintf('Trigger %4d: FFR (block %d, trial 1, start gap)\n', i, current_block);

    elseif is_end_gap && current_block == 1
        % Gap after last trial of block 1: this trigger is the MID marker
        block1_count = trial_count_in_block + n_missing_at_end;
        fprintf('Trigger %4d: MID marker (block 1 complete: %d trials, %d missing at end)\n', ...
                i, block1_count, n_missing_at_end);
        current_block = 2;
        trial_count_in_block = 0;

    elseif is_end_gap && current_block == 2
        % Gap after last trial of block 2: this trigger is the END marker
        block2_count = trial_count_in_block + n_missing_at_end;
        fprintf('Trigger %4d: END marker (block 2 complete: %d trials, %d missing at end)\n', ...
                i, block2_count, n_missing_at_end);
        if i < n_trigs
            fprintf('Trigger %d-%d: extra trigger(s) after END, labelled notrials\n', i+1, n_trigs);
        end
        break   % everything after END stays 'notrials'

    else
        % Extra trigger: keep 'notrials' and do not use it as reference
        fprintf('Trigger %4d: EXTRA trigger (%.3f s after last valid trigger), labelled notrials\n', ...
                i, interval/fs);
        continue
    end

    last_trig = trig(i);
end

% If no END marker was found, block 2 ends at the last counted trial
if current_block == 2 && block2_count == 0
    block2_count = trial_count_in_block;
    warning('No END marker found.');
end

%% Final validation
fprintf('\n=== Summary ===\n');

n_ffr = sum(strcmp(conditionlabels, 'FFR'));
n_mffr = sum(strcmp(conditionlabels, '-FFR'));
n_notrials = sum(strcmp(conditionlabels, 'notrials'));
n_valid = n_ffr + n_mffr;

fprintf('Total triggers: %d\n', n_trigs);
fprintf('FFR:           %d (expected 1500)\n', n_ffr);
fprintf('-FFR:          %d (expected 1500)\n', n_mffr);
fprintf('notrials:      %d (expected 3: START, MID, END)\n', n_notrials);
fprintf('Valid trials:  %d (expected 3000)\n', n_valid);
fprintf('Largest deviation of a trial interval: %.2f ms (tolerance %.1f ms)\n', ...
        max_dev/fs*1000, tol_trial/fs*1000);
fprintf('\nBlock 1: %d trials\n', block1_count);
fprintf('Block 2: %d trials\n', block2_count);

% Validation checks
if block1_count ~= n_per_block
    warning('Block 1 size unexpected: %d (expected %d)', block1_count, n_per_block);
end
if block2_count ~= n_per_block
    warning('Block 2 size unexpected: %d (expected %d)', block2_count, n_per_block);
end

fprintf('\nCondition labeling complete!\n');
end
