% Trial def of files with missing triggers
subjectid = 'HX';

% input
inputpath = fullfile('F:\Auditory-mov-data\',subjectid,'eeg'); % in case of normal data with pressure plate
name = '_FFR_sit_333_2';
inputname = ['P[',subjectid,']',name,'.dat'];

% read events once (same for both)
evt = ft_read_event(fullfile(inputpath,inputname));
codes = {evt.type};
trig_code = codes{3};

% inspect the triggers
trig = [evt.sample];
xx = diff(trig);
plot(xx)

% At some point a few triggers are skipped --> adjust the data labels accordingly

fs = 10000;

%% Assign condition labels
% Handles missing triggers correctly for 1500-trial blocks
% Total structure: 3003 triggers (1 START + 1500 block1 + 1 MID + 1500 block2 + 1 END)

n_trigs = length(trig);
conditionlabels = repmat({'notrials'}, 1, n_trigs);

% Define timing thresholds (in samples)
trial_dur = 0.24 * fs;      % ~240 ms between consecutive triggers
tol = 0.15 * fs;            % tolerance for interval matching

% Boundary gap signatures
gap_start_mid = 1.97*fs;  % Gap after START or MID markers
gap_block1_end = 3.27*fs; % Gap after last trial of block 1
gap_block2_end = 2.27*fs; % Gap after last trial of block 2

% State tracking
current_block = 1;
trial_count_in_block = 0;
block1_count = 0;
block2_count = 0;

fprintf('\n=== Assigning Condition Labels (1500-trial blocks) ===\n');
fprintf('Total triggers: %d (expected 3003)\n', n_trigs);

%% Label trigger 1 as START marker
conditionlabels{1} = 'notrials';


%% Process all intervals
for i = 1:length(xx)
    interval = xx(i);
    next_trig_idx = i + 1;  % Index of trigger after this interval
    
    %%% Classify the interval
    
    % Check for boundary gaps
    is_start_mid_gap = abs(interval - gap_start_mid) < 0.15*fs;
    is_block1_end_gap = abs(interval - gap_block1_end) < 0.15*fs;
    is_block2_end_gap = abs(interval - gap_block2_end) < 0.15*fs;
    
    if is_block1_end_gap
        % This is the gap after last trial of block 1
        % Next trigger is MID marker
        conditionlabels{next_trig_idx} = 'notrials';
        fprintf('Trigger %4d: MID marker (block 1 complete: %d trials)\n', ...
                next_trig_idx, trial_count_in_block);
        
        block1_count = trial_count_in_block;
        current_block = 2;
        trial_count_in_block = 0;
        
    elseif is_start_mid_gap
        % This is the gap after START or MID marker
        % Next trigger is first trial of a block
        current_block_temp = current_block;  % Which block are we entering?
        trial_count_in_block = 1;  % RESET counter
        
        % First trial of each block is always FFR
        conditionlabels{next_trig_idx} = 'FFR';
        fprintf('Trigger %4d: FFR (block %d, trial 1, start gap)\n', ...
                next_trig_idx, current_block_temp);
        
    elseif is_block2_end_gap
        % This is the gap after last trial of block 2
        % Next trigger is END marker
        conditionlabels{next_trig_idx} = 'notrials';
        fprintf('Trigger %4d: END marker (block 2 complete: %d trials)\n', ...
                next_trig_idx, trial_count_in_block);
        
        block2_count = trial_count_in_block;
        
    else
        % Regular trial interval (possibly with missing triggers)
        n_trials_in_interval = round(interval/trial_dur);
        expected_interval = n_trials_in_interval*trial_dur;
        deviation = abs(interval - expected_interval);
        
        if deviation < tol && n_trials_in_interval > 0
            % Valid trial interval
            % Count all trials in this interval (including missing ones)
            for k = 1:n_trials_in_interval
                trial_count_in_block = trial_count_in_block + 1;
                
                % Only label the actual recorded trigger (last in interval)
                if k == n_trials_in_interval
                    % Polarity: odd trials in block → FFR, even → -FFR
                    if mod(trial_count_in_block, 2) == 1
                        conditionlabels{next_trig_idx} = 'FFR';
                    else
                        conditionlabels{next_trig_idx} = '-FFR';
                    end
                    
                    if n_trials_in_interval > 1
                        fprintf('Trigger %4d: %s (block %d, trial %4d, %d skipped)\n', ...
                                next_trig_idx, conditionlabels{next_trig_idx}, ...
                                current_block, trial_count_in_block, n_trials_in_interval-1);
                    end
                end
            end
        else
            % Unrecognized interval
            conditionlabels{next_trig_idx} = 'notrials';
            fprintf('Trigger %4d: UNRECOGNIZED (gap: %.3fs)\n', ...
                    next_trig_idx, interval/fs);
        end
    end
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
fprintf('\nBlock 1: %d trials\n', block1_count);
fprintf('Block 2: %d trials\n', block2_count);

% Validation checks
if abs(block1_count - 1500) > 10
    warning('Block 1 size unexpected: %d (expected 1500)', block1_count);
end

if abs(block2_count - 1500) > 10
    warning('Block 2 size unexpected: %d (expected 1500)', block2_count);
end

if n_ffr ~= 1500 || n_mffr ~= 1500
    warning('Polarity distribution unexpected: FFR=%d, -FFR=%d', n_ffr, n_mffr);
end

if n_notrials ~= 3
    warning('Marker count unexpected: %d (expected 3)', n_notrials);
end

fprintf('\nCondition labeling complete!\n');

% Assign to data structure
data_segm.trialinfo = conditionlabels';



%% DIAGNOSTIC: Check for event type mismatch
% 
% evt = ft_read_event(fullfile(inputpath,inputname));
% codes = {evt.type};
% trig_code = codes{3};
% 
% % Check what event types exist
% unique_types = unique(codes);
% fprintf('=== Event Types Found ===\n');
% for i = 1:length(unique_types)
%     n_this_type = sum(strcmp(codes, unique_types{i}));
%     fprintf('Type: "%s"  -->  Count: %d\n', unique_types{i}, n_this_type);
% end
% 
% fprintf('\nSelected trig_code: "%s"\n', trig_code);
% fprintf('Total events (ALL types): %d\n', length(evt));
% fprintf('Events matching trig_code: %d\n', sum(strcmp(codes, trig_code)));
% 
% % Compare to the actual trial definition
% cfg2 = [];
% cfg2.dataset = fullfile(inputpath,inputname);
% cfg2.trialdef.eventtype = trig_code;
% cfg2.trialdef.prestim  = 0.01;
% cfg2.trialdef.poststim = 0.240;
% cfg_trl_check = ft_definetrial(cfg2);
% 
% fprintf('\nNumber of trials in cfg_trl.trl: %d\n', size(cfg_trl_check.trl,1));
% fprintf('Number of elements in "trig" (all events): %d\n', length(evt));
% 
% if length(evt) ~= size(cfg_trl_check.trl,1)
%     fprintf(2, '\n*** MISMATCH DETECTED! ***\n');
%     fprintf(2, 'trig (%d) does NOT match number of actual trials (%d)\n', ...
%             length(evt), size(cfg_trl_check.trl,1));
%     fprintf(2, 'This means your conditionlabels are likely MISALIGNED with your epoched trials!\n');
% else
%     fprintf('\nCounts match - no obvious misalignment from event type mixing.\n');
% end