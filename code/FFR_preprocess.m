%% Preprocessing of the FFR data
clc
clear
close all

restoredefaultpath 
addpath 'C:\Users\nicol\AppData\Roaming\MathWorks\MATLAB Add-Ons\Collections\FieldTrip' 
ft_defaults

%% Set parameters
filepath = 'C:\Users\nicol\OneDrive\Desktop\Leuven\TESI\Reading_FFR\PD\PD12';
filename = 'PwPD12_FFR.eeg';
subjectID = 'PD12_FFR';

templ_file = 'easycapM11.mat';


% % Load data
cfg            = [];
cfg.dataset    = fullfile(filepath,filename);
cfg.continuous = 'yes';
cfg.channel    = 'all';
data           = ft_preprocessing(cfg);

OrigEEGchan = ft_channelselection('eeg', data);

fs = data.fsample;
%% Filtering
% Using hp=300/lp=400 (Butterworth, onepass) --> very narrow, causal band-pass that can distort phase and dynamics.
% For FFRs, a safer default is a broader, zero-phase band-pass (e.g.,[250–1000] Hz for 333 Hz F0)
% Use twopass (acausal) unless you need causal filtering.
% --- Band-pass for FFR (333 Hz target, keep harmonics) ---
cfg = [];
cfg.bpfilter  = 'yes';
cfg.channel   = OrigEEGchan;    % don't filter trigger channel
cfg.bpfreq    = [250 500];     % or [80 2000] if you want very broad  
cfg.bpfiltord = 6;
cfg.bpfilttype= 'but';
cfg.bpfiltdir = 'twopass';     % zero-phase!
data_bp = ft_preprocessing(cfg, data);


% Optional: DFT at 50 Hz or harmonics only if needed
% Prefer:
%   - No explicit notch (the band-pass may already suppress 50 Hz), or
%   - A DFT filter at 50 Hz only if line noise is visible (avoid harmonics in this paradigm).
% cfg = [];
% cfg.dftfilter = 'yes';
% cfg.dftfreq   = [300];             % avoid harmonics like 350
% data_bp = ft_preprocessing(cfg, data_bp);


%% Detection & remove bad channels
% Quick automatic screening (robust RMS / flatline) (continuous, filtered data)

% Taking only the usefull part of signal
% cfg = [];
% cfg.dataset = fullfile(filepath,filename);
% evt = ft_read_event(cfg.dataset);
% sample = [evt.sample];
% cfg = [];
% cfg.latency = [sample(3)-2*fs sample(end)+2*fs] / fs;
% data_bp = ft_selectdata(cfg, data_bp);

X = data_bp.trial{1};                  % [nChan x nSamples]  
rms_ch = sqrt(mean(X.^2, 2));          % per-channel RMS

robust_c = median(rms_ch);             % robust center
robust_s = mad(rms_ch, 1);             % robust spread         

z_rms = (rms_ch - robust_c)./max(robust_s, eps);

% Flatline detection (very low variance)
ptp_ch = max(X,[],2) - min(X,[],2);     % peak-to-peak
ptp_med = median(ptp_ch);

% Mark channels with extreme RMS or flatline
bad_rms  = abs(z_rms) > 3.1;                       % robust outlier
bad_flat_abs = ptp_ch < 1;                       % ONLY IF units are µV
bad_flat_rel = ptp_ch < 0.02 * ptp_med;          % unit-free fallback

bad_idx = find(bad_rms | bad_flat_abs | bad_flat_rel);
badlabels = data_bp.label(bad_idx);

fprintf('Bad channels flagged: %d/%d\n', numel(badlabels), numel(data_bp.label));

%%
% Decide: interpolate vs exclude bad channels (based on fraction bad)
% Set the treshold for discarding vs interpolation ---
interp_max_frac  = 0.10;   % <= 10% bad -> interpolate (keeps montage consistent)
exclude_max_frac = 0.20;   % 10-15% bad -> exclude (still salvageable)
% >15% bad -> stop / mark dataset as poor quality

badlabels = intersect(badlabels, OrigEEGchan);   % only EEG bad channels
nEEG   = numel(OrigEEGchan);
nBad   = numel(badlabels);
badFrac = nBad / max(nEEG,1);

fprintf('EEG bad channels: %d/%d (%.1f%%)\n', nBad, nEEG, 100*badFrac);

% Protect reference electrodes
% refChans = {'TP9','TP10'};
% badRef = intersect(badlabels, refChans);
% if ~isempty(badRef)
%     warning('Reference channel(s) flagged as bad: %s. Consider changing reference strategy or repairing first.', ...
%         strjoin(badRef, ', '));
% end

% Visual inspection of bad channels quality using ft_databrowser
cfg = [];
cfg.viewmode = 'vertical'; % 'vertical' or 'butterfly' view
cfg.continuous = 'yes'; % 'yes' for continuous data
cfg.blocksize = 10; % Set the blocksize (number of seconds per window)
cfg.channel = badlabels; % 'all' or specify channels of interest
ft_databrowser(cfg, data_bp);


%%
% This next part of the code interpolates or rejects the bad channels. In
% this case the interpolation is based on the personalized EEG electrode
% locations (registered during the EEG session) but our dataset does not
% have this so interpolation is not the best option then. We could
% interpolate based on the EEG electrode template but this is not always
% equally reliable.
 
if nBad == 0
    data_nobad = data_bp;
    excluded = [];
    
elseif badFrac <= interp_max_frac
    fprintf('Action: INTERPOLATE bad channels.\n');
    
    % Loading layout from template easycapM11.mat
    cfg = [];
    cfg.layout = templ_file;          
    layout = ft_prepare_layout(cfg);
    
    cfg = [];
    cfg.method = 'triangulation';
    cfg.layout = layout;
    cfg.channel = OrigEEGchan;
    neighbours = ft_prepare_neighbours(cfg, data_bp);
    
    % Interpolation of bad channels
    cfg = [];
    cfg.badchannel = badlabels;
    cfg.neighbours = neighbours;
    cfg.method = 'average';            
    data_nobad = ft_channelrepair(cfg, data_bp);
    excluded = [];
    
    fprintf('Interpolated channels: %s\n', strjoin(badlabels, ', '));
    
elseif badFrac <= exclude_max_frac
    fprintf('Action: EXCLUDE bad channels.\n');
    
    cfg = [];
    cfg.channel = setdiff(data_bp.label, badlabels, 'stable');
    data_nobad = ft_selectdata(cfg, data_bp);
    excluded = badlabels;
    
    fprintf('Excluded channels: %s\n', strjoin(badlabels, ', '));
    
else
    error('Too many bad channels (%.1f%%). Consider excluding this file/block.', 100*badFrac);
end
    

% Don’t over-interpolate. 
% If many channels are bad (>10–15%), results become unreliable—consider excluding that run/subject.

%% Re-referencing
% Do NOT interpolate reference electrodes you'll use for rereferencing (e.g., TP9/TP10). 
% If they are bad, pick a different reference strategy (e.g., average reference after repairing the rest).

% mastoid rereference: offline linked-ears reference
% haveTP9  = any(strcmp(data_bp.label,'TP9'));
% haveTP10 = any(strcmp(data_bp.label,'TP10'));
% 
% % Visual inspection of mastoid quality using ft_databrowser
% cfg = [];
% cfg.viewmode = 'vertical'; % 'vertical' or 'butterfly' view
% cfg.continuous = 'yes'; % 'yes' for continuous data
% cfg.blocksize = 10; % Set the blocksize (number of seconds per window)
% cfg.channel = {'TP9','TP10'}; % 'all' or specify channels of interest
% ft_databrowser(cfg, data_bp);
% 
% 
% if ~isempty(badRef)
%     cfg = [];
%     cfg.reref     = 'yes';
%     cfg.refmethod = 'avg';
%     cfg.channel   = 'all';
%     cfg.refchannel = setdiff(OrigEEGchan, badlabels);
%     data_ref = ft_preprocessing(cfg, data_nobad);
%     fprintf('Reference channel(s) flagged as bad: %s. Average reference is used\n',  strjoin(badRef, ', '));
% 
% elseif haveTP9 && haveTP10 && isempty(badRef)
%     cfg = [];
%     cfg.reref      = 'yes';
%     cfg.refmethod  = 'avg';
%     cfg.refchannel = {'TP9','TP10'};
%     cfg.channel    = 'all';
%     %cfg.implicitref = 'FCz';               % Bob
%     %cfg.implicitref = 'AFz';               % others
%     data_ref = ft_preprocessing(cfg, data_nobad);
%     fprintf('Using TP9/TP10 as reference\n');
% else
%     error('TP9/TP10 not found. Check mastoid labels (e.g., A1/A2, M1/M2).');
% end

cfg = [];
cfg.reref     = 'yes';
cfg.refmethod = 'avg';
cfg.channel   = 'all';
cfg.implicitref = 'FCz';
cfg.refchannel = setdiff(OrigEEGchan, excluded);
data_ref = ft_preprocessing(cfg, data_nobad);
fprintf('Average reference is used\n');

%% Define trials
% Epoch the data according to the trial definition

% % check for the trigger codes
cfg.dataset = fullfile(filepath,filename);
evt = ft_read_event(cfg.dataset);
codes = {evt.type};
trig_code = codes{3};

% Defining trials based on event type
cfg = [];
cfg.dataset = fullfile(filepath,filename);
cfg.trialdef.eventtype = trig_code; % Trigger event name
cfg.trialdef.prestim = 0.01;
cfg.trialdef.poststim = 0.240;
cfg_trl = ft_definetrial(cfg);

% Segment the continuous data based on trials
data_segmented = ft_redefinetrial(cfg_trl, data_ref); 

%% MANUAL TRIALS SEGMENTATION
% cfg.dataset = fullfile(filepath,filename);
% evt = ft_read_event(cfg.dataset);
% trig = [evt.sample];
% trigger_sample = trig(3);
% 
% nTrials = 2003;
% 
% pretrig  = round(0.01  * fs);
% posttrig = round(0.240 * fs);
% 
% step_std = round(0.24 * fs);
% step_long = round(1.97 * fs);
% gap_1001_1002 = round(3.27 * fs);
% gap_1002_1003 = round(1.97 * fs);
% 
% steps = step_std * ones(nTrials-1,1);
% 
% steps(1) = step_long;         % 1 -> 2
% steps(1001) = gap_1001_1002;  % 1001 -> 1002
% steps(1002) = gap_1002_1003;  % 1002 -> 1003
% steps(2002) = step_long;      % 2002 -> 2003
% 
% trig = zeros(nTrials,1);
% trig(1) = trigger_sample;
% 
% for i = 2:nTrials
%     trig(i) = trig(i-1) + steps(i-1);
% end
% 
% trl = zeros(nTrials,3);
% 
% for i = 1:nTrials
%     
%     begsample = trig(i) - pretrig;
%     endsample = trig(i) + posttrig;
%     offset    = -pretrig;
%     
%     trl(i,:) = [begsample endsample offset];
% end
% 
% cfg = [];
% cfg.trl = trl;
% 
% data_segmented = ft_redefinetrial(cfg, data_ref);

%% Assign the conditionlabels
% trig = [evt.sample];
% trig = trig(3:end);     % First 2 triggers are system initialisation

xx = diff(trig);

n_trigs = length(trig);
conditionlabels = repmat({'notrials'}, 1, n_trigs);

toggle_ind = 1; % 1 = FFR, -1 = -FFR


% Define timing thresholds (in samples)
trial_dur  = 0.24 * fs;  % ~230 ms: expected gap between consecutive triggers
tol        = 0.15 * fs;  % tolerance window around each expected multiple

% Boundary gaps: mid-block marker and start marker are ~1.97s
% These are distinct from end-of-block-1 (~3.27s) and end-of-block-2 (~2.27s)
% which are just long inter-trial gaps with many missed triggers.
% So we define a narrow boundary window around 1.97s only.
boundary_low  = 1.85 * fs;  % lower bound for true boundary gap
boundary_high = 2.1 * fs;  % upper bound for true boundary gap

% Structural event gaps in samples
first_gap = 1.97 * fs;     % start/mid marker to first real trial

% Variables for counting blocks size
tot_blocksize = 0;
block1_count = 0;
block2_count = 0;
current_block = 1;



for i = 1:length(xx)
    interval = xx(i);
    
    % Compute how far the current interval deviates from the expected boundary gap,
    % estimate how many trial durations (k) fit into this deviation,
    % and classify the interval as a boundary either if it falls within the
    % predefined boundary window or if it matches a multiple of the trial duration.
    delta = interval - first_gap;
    k = round(delta / trial_dur);
    is_boundary = (interval > boundary_low && interval < boundary_high) || ...
    (k >= 0 && abs(delta - k * trial_dur) < 0.04*fs);
    
    % Handle very large gaps separately: if the mid-measurement trigger is missing,
    % the time interval becomes unusually long (approximately 3.27 s + 1.97 s,
    % plus possible multiples of the trial duration ~0.24 s). In this case,
    % recompute delta and k assuming a longer expected gap and force boundary classification.
    if interval >= 5.24 * fs
        delta = interval - 5.24 * fs;
        k = round(delta / trial_dur);
        is_boundary = 1;
    end


    % Number of triggers skipped in this interval
    n_skipped = round(interval / trial_dur) - 1;
    is_valid_trial = n_skipped >= 0 && ...
                     abs(interval - (n_skipped + 1) * trial_dur) < tol;


    % Enter boundary handling only if the interval is classified as a boundary
    % and we are not within the central portion of the block (to avoid false boundaries).
    if is_boundary && ~(tot_blocksize>10 && tot_blocksize<800)
        % Mid-block or end marker — notrials, reset toggle
        increment =  1;
        conditionlabels{i} = 'notrials';
        toggle_ind = 1;
        
        % If multiple trials are spanned, update increment accordingly
        % and flip toggle if the number of skipped trials is odd.
        if k > 0 
            increment = k + 1;
            if  mod(k, 2) == 1
                toggle_ind = -toggle_ind;
            end
        end

    elseif is_valid_trial
        % Label current trigger
        if toggle_ind == 1
            conditionlabels{i} = 'FFR';
        else
            conditionlabels{i} = '-FFR';
        end


        % If there are skipped triggers in this interval
        if n_skipped > 0
            % Ignore skipped count if interval matches start/mid boundary (~1.97s)
            if (interval > 1.95 * fs && interval < 2.05 * fs)
                increment =  1;

            % Ignore skipped count if interval matches end of block 1 (~3.27s)    
            elseif abs(interval - 3.27*fs) < 0.04*fs
                increment =  1;

            % A gap of ~2.27 s is expected only at the final trigger (end of the experiment);
            % if the same interval occurs earlier, it should be interpreted as skipped trials.  
            elseif abs(interval - 2.27*fs) < 0.04*fs && tot_blocksize > 1990
                increment =  1;

            % Otherwise, count skipped trials normally   
            else
                increment = n_skipped + 1;
            end

        % No skipped triggers: count single trial    
        else
            increment =  1;
        end


        % Advance toggle by (n_skipped + 1) steps
        if mod(n_skipped + 1, 2) == 1
            toggle_ind = -toggle_ind;
        end

    else
        % Unrecognised interval — notrials, reset toggle
        conditionlabels{i} = 'notrials'; 
        increment =  0;
    end

    % Detect transition to second block using a large gap (~3.27s)
     if interval > 3.26 * fs  
        current_block = 2;
        block1_count = block1_count + increment; % Assign this increment to block 1 (last event before transition)
        tot_blocksize = block1_count + block2_count;
        %disp([block2_count, block1_count, tot_blocksize, increment])
        continue
     end

     % Update trial count based on current block
     if current_block == 1
         block1_count = block1_count + increment;
     else 
         block2_count = block2_count + increment;
     end
     
     tot_blocksize = block1_count + block2_count;
        
     %disp([block2_count, block1_count, tot_blocksize,increment])
end


% Last trigger is end marker — already labelled 'notrials' by initialisation

data_segmented.trialinfo = conditionlabels';

% tot_blocksize should be 2002 (1001 for first block and 1001 for second
% one), the central trigger is not taken into consideration.
if tot_blocksize ~= 2002
    warning('Block size mismatch: B1=%d, B2=%d', block1_count, block2_count);
end


disp(sum(strcmp(conditionlabels, 'FFR')))
disp(sum(strcmp(conditionlabels, '-FFR')))


%% Baseline + artifact rejection

% Baseline first (on epoched data)
% cfg = [];
% cfg.demean         = 'yes';
% cfg.keeptrialinfo  = 'yes';
% cfg.baselinewindow = [-0.010 0.020];  %0.020
% data_bl = ft_preprocessing(cfg, data_segmented);

data_bl = data_segmented;

% Fast per-trial absolute-amplitude rejection
cfg = [];
cfg.latency = [0.030 0.230];             % check only the response window (in seconds!)   
data_toi = ft_selectdata(cfg, data_bl);    


ntr = numel(data_toi.trial);
trial_rms = zeros(ntr,1);

for j = 1:ntr
    Z = data_toi.trial{j};                 % [nChan x nSamples]
    rms_ch = sqrt(mean(Z.^2, 2));          % RMS for channel
    trial_rms(j) = median(rms_ch);         
end

% robust statistics tra trial
robust_c = median(trial_rms);
robust_s = mad(trial_rms, 1);
z_trial = (trial_rms - robust_c) ./ max(robust_s, eps);

% threshold
rej = abs(z_trial) > 3;

% keep the trials that are NOT rejected
keep = find(~rej);

% apply the trial selection to the FULL baselined dataset (not only TOI)
cfg = [];
cfg.trials = keep;
data_clean = ft_selectdata(cfg, data_bl);


% reporting: counts + percentage
n_before = numel(data_bl.trial);
n_after  = numel(data_clean.trial);
n_rej    = n_before - n_after;
pct_rej  = 100 * n_rej / max(n_before,1);

%fprintf('Fast amplitude rejection (±%g):\n', thr);
fprintf('Trials before: %d\n', n_before);
fprintf('Trials after : %d\n', n_after);
fprintf('Rejected     : %d (%.1f%%)\n', n_rej, pct_rej);


%% Possibly additional artifact rejection for the EEG channels
cfg        = [];
cfg.channel = 'all';
cfg.metric = 'var';  % use by default zvalue method
cfg.method = 'summary'; % use by default summary method
data_final = ft_rejectvisual(cfg,data_clean);


% Save the preprocessed data
name = ['preprocessed_', subjectID];
outputdir = 'C:\Users\nicol\OneDrive\Desktop\Leuven\TESI\Reading_FFR\Results\PD\PD12';
path = fullfile(outputdir,name);
save(path,"data_final","-v7.3");


%% Averaging of the different trials
% Select the data per data label
cfg = [];
cfg.trials = strcmp(data_final.trialinfo, 'FFR');
dataFFR = ft_selectdata(cfg, data_final);

cfg = [];
cfg.trials = strcmp(data_final.trialinfo, '-FFR');
dataMinFFR = ft_selectdata(cfg, data_final);

% Average the data per condition
cfg = [];
cfg.trials = 'all';
avgFFR = ft_timelockanalysis(cfg, dataFFR);

cfg = [];
cfg.trials = 'all';
avgMinFFR = ft_timelockanalysis(cfg, dataMinFFR);

% Compute the difference between the polarities, FFR and -FFR and divide by 2
final_avg = (avgFFR.avg - avgMinFFR.avg)/2;

final = avgFFR;
final.avg = final_avg;

% Save the result
name = ['averagedFFR_', subjectID];
outputdir = 'C:\Users\nicol\OneDrive\Desktop\Leuven\TESI\Reading_FFR\Results\PD\PD12';
path = fullfile(outputdir,name);
save(path,"final","-v7.3");

