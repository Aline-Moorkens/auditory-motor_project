%% Preprocessing of the FFR data
clc
clear
close all

% restoredefaultpath 
% addpath 'C:\Users\Aline Moorkens\OneDrive - KU Leuven\Documents\fieldtrip-20250106'
% ft_defaults

%% Set path
subjectid = 'AM';

% input
inputpath = fullfile('F:\Auditory-mov-data',subjectid,'eeg');
%name = '_FFR_sit_333';
%name = '_FFR_firm_333';
name = '_FFR_foam_333';
inputname = ['P[',subjectid,']',name,'.dat'];

% output
outputpath = fullfile('F:\Auditory-mov-data',subjectid,'preprocessed');
outputname = ['P[',subjectid,']',name];

% layout
templ_file = 'easycapM11.mat';


% % Load data
cfg            = [];
cfg.dataset    = fullfile(inputpath,inputname);
cfg.continuous = 'yes';
cfg.channel    = 'all';
data           = ft_preprocessing(cfg);

OrigEEGchan = ft_channelselection('eeg', data);

fs = data.fsample;

% % Filtering
cfg = [];
cfg.bpfilter  = 'yes';
cfg.channel   = OrigEEGchan;    % don't filter trigger channel
cfg.bpfreq    = [200 900];     % or [80 2000] if you want very broad  
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


%% Detection of bad channels (FASTER / PREP approach)
% Quick automatic screening (robust RMS / flatline) (continuous, filtered data)
%%% 1. Variance (FASTER)
X = data_bp.trial{1};               % [nChan x nSamples]  
var_ch = var(X, 0, 2);              % variance per channel, [nChan x 1]

robust_c = median(var_ch);
robust_s = mad(var_ch, 1);
z_var = (var_ch - robust_c) ./ max(robust_s, eps);

bad_var = abs(z_var) > 4;

%%% 2. Flatline detection (very low variance)
ptp_ch = max(X,[],2) - min(X,[],2);     % peak-to-peak
ptp_med = median(ptp_ch);

bad_flat_abs = ptp_ch < 1;                       % ONLY IF units are µV
bad_flat_rel = ptp_ch < 0.02 * ptp_med;          % unit-free fallback

%%% 3. correlation (FASTER)
% Compute pairwise correlations on bandpass-filtered data
C = corr(X');                         % [nChan x nChan]
C(logical(eye(size(C)))) = NaN;       % remove self-correlations
mean_corr = mean(C, 2, 'omitnan');    % mean correlation per channel

robust_mc = median(mean_corr);
robust_ms = mad(mean_corr, 1);
z_corr = (mean_corr - robust_mc) ./ max(robust_ms, eps);

bad_corr = z_corr < -4;              % low correlation = bad channel

%%% 4. HF noisiness criterion (PREP)
% % Ratio of upper-band power (500-900 Hz) to lower-band power (200-500 Hz)
% % A noisy channel will have disproportionately high power in the upper band
% b_low  = fir1(200, [200 500]/(fs/2), 'bandpass');
% b_high = fir1(200, [500 900]/(fs/2), 'bandpass');
% 
% X_low  = filtfilt(b_low,  1, X')';
% X_high = filtfilt(b_high, 1, X')';
% 
% pow_low  = mean(X_low.^2,  2);
% pow_high = mean(X_high.^2, 2);
% 
% noise_ratio = pow_high ./ max(pow_low, eps);
% 
% robust_cn = median(noise_ratio);
% robust_sn = mad(noise_ratio, 1);
% z_noise   = (noise_ratio - robust_cn) ./ max(robust_sn, eps);
% 
% bad_hf_noise = z_noise > 3;

%%% Find the bad channel labels
bad_idx = find(bad_var | bad_flat_abs | bad_flat_rel|bad_corr);
%bad_idx = find(bad_var | bad_flat_abs | bad_flat_rel|bad_corr|bad_hf_noise);
badlabels = data_bp.label(bad_idx);

fprintf('Bad channels flagged: %d/%d\n', numel(badlabels), numel(data_bp.label));

fprintf('Bad by variance:  %d channels\n', sum(bad_var));
fprintf('Bad by flatline:  %d channels\n', sum(bad_flat_abs | bad_flat_rel));
fprintf('Bad by low corr:  %d channels\n', sum(bad_corr));
%fprintf('Bad by HF noise:  %d channels\n', sum(bad_hf_noise));

%% Visual inspection of bad channels

badlabels = intersect(badlabels, OrigEEGchan);   % only EEG bad channels
nEEG   = numel(OrigEEGchan);
nBad   = numel(badlabels);
badFrac = nBad / max(nEEG,1);

fprintf('EEG bad channels: %d/%d (%.1f%%)\n', nBad, nEEG, 100*badFrac);

% Protect reference electrodes
refChans = {'TP9','TP10'};
badRef = intersect(badlabels, refChans);
if ~isempty(badRef)
    warning('Reference channel(s) flagged as bad: %s. Consider changing reference strategy or repairing first.', ...
        strjoin(badRef, ', '));
end

% Visual inspection of bad channels quality using ft_databrowser
cfg = [];
cfg.viewmode = 'vertical'; % 'vertical' or 'butterfly' view
cfg.continuous = 'yes'; % 'yes' for continuous data
cfg.blocksize = 10; % Set the blocksize (number of seconds per window)
cfg.channel = badlabels; % 'all' or specify channels of interest
ft_databrowser(cfg, data_bp);

%% Define the badchannels
badlabels = {'P8','C5','FT7','F2', 'FC5', 'T7'};
nEEG   = numel(OrigEEGchan);
nBad   = numel(badlabels);
badFrac = nBad / max(nEEG,1);
fprintf('EEG bad channels: %d/%d (%.1f%%)\n', nBad, nEEG, 100*badFrac);

%% Bad channel interpolation/ removal
% Decide: interpolate vs exclude bad channels (based on fraction bad)
% Set the treshold for discarding vs interpolation ---
interp_max_frac  = 0.20;   % <= 10% bad -> interpolate (keeps montage consistent)
exclude_max_frac = 0.25;   % 10-15% bad -> exclude (still salvageable)
% >20% bad -> stop / mark dataset as poor quality

% This next part of the code interpolates or rejects the bad channels.
if nBad == 0
    data_nobad = data_bp;
    excluded = [];
    
elseif badFrac <= interp_max_frac
    fprintf('Action: INTERPOLATE bad channels.\n');
    % Use 2D layout-based neighbours + simple interpolation
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
haveTP9  = any(strcmp(data_nobad.label,'TP9'));
haveTP10 = any(strcmp(data_nobad.label,'TP10'));

% % Visual inspection of mastoid quality using ft_databrowser
cfg = [];
cfg.viewmode = 'vertical'; % 'vertical' or 'butterfly' view
cfg.continuous = 'yes'; % 'yes' for continuous data
cfg.blocksize = 10; % Set the blocksize (number of seconds per window)
cfg.channel = {'TP9','TP10'}; % 'all' or specify channels of interest
ft_databrowser(cfg, data_nobad);


cfg = [];
cfg.reref      = 'yes';
cfg.refmethod  = 'avg';
cfg.refchannel = {'TP9','TP10'};
cfg.channel    = 'all';
data_ref = ft_preprocessing(cfg, data_nobad);
fprintf('Using TP9/TP10 as reference\n');

% OR CAR
% cfg = [];
% cfg.reref     = 'yes';
% cfg.refmethod = 'avg';
% cfg.channel   = 'all';
% cfg.implicitref = 'FCz';
% cfg.refchannel = setdiff(OrigEEGchan, excluded);
% data_ref = ft_preprocessing(cfg, data_nobad);
% fprintf('Average reference is used\n');

%% Define trials
% Epoch the data according to the trial definition

% % check for the trigger codes
cfg.dataset = fullfile(inputpath,inputname);
evt = ft_read_event(cfg.dataset);
codes = {evt.type};
trig_code = codes{3};

% Defining trials based on event type
cfg = [];
cfg.dataset = fullfile(inputpath,inputname);
cfg.trialdef.eventtype = trig_code; % Trigger event name
cfg.trialdef.prestim = 0.01;
cfg.trialdef.poststim = 0.240;
cfg_trl = ft_definetrial(cfg);

% Segment the continuous data based on trials
data_segm = ft_redefinetrial(cfg_trl, data_ref); 


%% Assign the conditionlabels
n_trigs = size(cfg_trl.trl, 1);  % Should be 3003

% Initialize all as 'notrials'
conditionlabels = repmat({'notrials'}, n_trigs, 1);
 
fprintf('Total triggers: %d (expected 3003)\n', n_trigs);

% Block 1: Triggers 2-1501 (1500 trials)
for i = 2:1501
    trial_num_in_block = i - 1;  % Trial 1, 2, 3, ...
    if mod(trial_num_in_block, 2) == 1
        conditionlabels{i} = 'FFR';   % Odd trials: FFR
    else
        conditionlabels{i} = '-FFR';  % Even trials: -FFR
    end
end

% Trigger 1502: MID marker (already 'notrials')
 
% Block 2: Triggers 1503-3002 (1500 trials)
for i = 1503:3002
    trial_num_in_block = i - 1502;  % Trial 1, 2, 3, ... (resets for block 2)
    if mod(trial_num_in_block, 2) == 1
        conditionlabels{i} = 'FFR';   % Odd trials: FFR
    else
        conditionlabels{i} = '-FFR';  % Even trials: -FFR
    end
end
 
% Trigger 3003: END marker (already 'notrials')


% Assign the trial labels
data_segm.trialinfo = conditionlabels;

%% Baseline + artifact rejection

% Baseline first (on epoched data)
% cfg = [];
% cfg.demean         = 'yes';
% cfg.keeptrialinfo  = 'yes';
% cfg.baselinewindow = [-0.010 0.020];  %0.020
% data_bl = ft_preprocessing(cfg, data_segm);

data_bl = data_segm;

% Fast per-trial absolute-amplitude rejection
cfg = [];
cfg.latency = [0.030 0.230];             % check only the response window (in seconds!)   
data_toi = ft_selectdata(cfg, data_bl);    

ntr = numel(data_toi.trial);

% trial_rms = zeros(ntr,1);
% 
% for j = 1:ntr
%     Z = data_toi.trial{j};                 % [nChan x nSamples]
%     rms_ch = sqrt(mean(Z.^2, 2));          % RMS for channel
%     trial_rms(j) = median(rms_ch);         
% end
% 
% % robust statistics tra trial
% robust_c = median(trial_rms);
% robust_s = mad(trial_rms, 1);
% z_trial = (trial_rms - robust_c) ./ max(robust_s, eps);

trial_var = zeros(ntr,1);

for j = 1:ntr
    Z = data_toi.trial{j};
    var_ch = var(Z, 0, 2);
    trial_var(j) = median(var_ch);
end

robust_c = median(trial_var);
robust_s = mad(trial_var, 1);
z_trial  = (trial_var - robust_c) ./ max(robust_s, eps);

% threshold
rej = abs(z_trial) > 4;

% keep the trials that are NOT rejected
keep = find(~rej);

% apply the trial selection to the FULL baselined dataset (not only TOI)
cfg = [];
cfg.trials = keep;
cfg.keeptrialinfo = 'yes';
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
name = ['preprocessed_', outputname];
output = fullfile(outputpath,name);
save(output,"data_final","-v7.3");


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
name = ['averagedFFR_', outputname];
output = fullfile(outputpath,name);
save(output,"final","-v7.3");