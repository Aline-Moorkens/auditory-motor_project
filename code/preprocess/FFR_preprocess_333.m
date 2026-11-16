%% Preprocessing of the FFR data
clc
clear
close all

addpath 'C:\Users\Aline Moorkens\OneDrive - KU Leuven\Documents\fieldtrip-20250106'
ft_defaults

%% Settings

% SubjectID

%%% Batch 1: 333 - 89
% subjectid = 'NSMT';
% subjectid = 'AM';
% subjectid = 'FS';
% subjectid = 'HX';     % triggers missing in sit_2 file
% subjectid = 'LVB';
% subjectid = 'RNR';
% subjectid = 'VUR';
% subjectid = 'LH';
% subjectid = 'TF';
% subjectid = 'MGZ';
% subjectid = 'MM';
% subjectid = 'TH';
% subjectid = 'TD';
% subjectid = 'FD'; 
% subjectid = 'EF';     % triggers missing in sit_2 file
% subjectid = 'MCBG';
% subjectid = 'LV';
% subjectid = 'CX';     % triggers missing in sit_2 file
% subjectid = 'NG';
%subjectid = 'BVD';

%%%  Batch 2: 333 - 433
% subjectid = 'RR';
% subjectid = 'AF';
% subjectid = 'LM';
% subjectid = 'ZB';
% subjectid = 'CF';
% subjectid = 'CB';
% subjectid = 'ECM';
% subjectid = 'ML';
% subjectid = 'NS';
% subjectid = 'TM';
% subjectid = 'LPB';
% subjectid = 'SK';
% subjectid = 'SG';
% subjectid = 'EC';
% subjectid = 'SP';

% Postural condition
%name = '_FFR_sit_333';
%name = '_FFR_firm_333';
%name = '_FFR_foam_333';
name = '_FFR_sit_333_2';


% Paths
inputpath  = fullfile('F:\Auditory-mov-data', subjectid, 'eeg');
basepath   = fullfile('F:\Auditory-mov-data', subjectid);
inputname  = ['P[', subjectid, ']', name, '.dat'];

outputpath = fullfile('F:\Auditory-mov-project\main',subjectid, 'preprocessed');
outputname = ['P[', subjectid, ']', name];
if ~exist(outputpath, 'dir'), mkdir(outputpath); end  % make sure folders exist


logpath = 'C:\Users\Aline Moorkens\OneDrive - KU Leuven\PhD\auditory-motor_project\subjects';
logfile = fullfile(logpath, sprintf('FFR_preprocessing_333.xlsx'));
badFile = fullfile(logpath, 'FFR_bad_channels.csv');        % one row per file

% electrode layout
templ_file = 'easycapM11.mat';


%% Main analysis
%%% Load data
cfg            = [];
cfg.dataset    = fullfile(inputpath, inputname);
cfg.continuous = 'yes';
cfg.channel    = 'all';
data           = ft_preprocessing(cfg);

OrigEEGchan = ft_channelselection('eeg', data);
fs = data.fsample;

% The 1 uV flatline criterion assumes microvolts: print the units once
if isfield(data, 'hdr') && isfield(data.hdr, 'chanunit')
    fprintf('Channel units: %s\n', strjoin(unique(data.hdr.chanunit(:))', ', '));
end


%%% Filtering: zero-phase Butterworth band-pass, 200-900 Hz, order 6
cfg = [];
cfg.bpfilter  = 'yes';
cfg.channel   = OrigEEGchan;    % don't filter trigger channel
cfg.bpfreq    = [200 900];     % or [80 2000] if you want very broad  
cfg.bpfiltord = 6;
cfg.bpfilttype= 'but';
cfg.bpfiltdir = 'twopass';     % zero-phase!
data_bp = ft_preprocessing(cfg, data);


%% Automatic bad-channel detection (FASTER-style, continuous filtered data)
X = data_bp.trial{1};               % [nChan x nSamples]

% 1. Variance, robust z (two-sided)
var_ch = var(X, 0, 2);              % variance per channel, [nChan x 1]
robust_c = median(var_ch);
robust_s = mad(var_ch, 1);
z_var = (var_ch - robust_c) ./ max(robust_s, eps);
bad_var = abs(z_var) > 4;

% 2. Flatline: peak-to-peak < 1 uV or < 2% of the median peak-to-peak
ptp_ch       = max(X, [], 2) - min(X, [], 2);   % peak-to-peak
bad_flat_abs = ptp_ch < 1;
bad_flat_rel = ptp_ch < 0.02 * median(ptp_ch);


% 3. Mean correlation with the other channels, robust z < -4
% Compute pairwise correlations on bandpass-filtered data
C = corr(X');                         % [nChan x nChan]
C(logical(eye(size(C)))) = NaN;       % remove self-correlations
mean_corr = mean(C, 2, 'omitnan');    % mean correlation per channel
z_corr    = (mean_corr - median(mean_corr)) ./ max(mad(mean_corr, 1), eps);
bad_corr = z_corr < -4;              % low correlation = bad channel

% 4. HF noisiness criterion (PREP)
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
bad_idx       = find(bad_var | bad_flat_abs | bad_flat_rel | bad_corr);
suggested_bad = intersect(data_bp.label(bad_idx), OrigEEGchan);

fprintf('Suggested bad channels: %s\n', strjoin(suggested_bad, ', '));
fprintf('  by variance: %d | flatline: %d | low correlation: %d\n', ...
        sum(bad_var), sum(bad_flat_abs | bad_flat_rel), sum(bad_corr));


%%% Visual inspection of bad channels

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


%% Define the bad channels
%badlabels = {'POz', 'P3', 'C1', 'FC5', 'FC6', 'F5'};  %'O1','Pz', 'FT10', 'P7', 'P1'};%, 'TP7', 'C2'};%'F7'}; %{'Oz', 'PO7', 'PO8', 'Pz', 'P2', 'FC1', 'F1', 'Fz', ' AF4'};%;%'F2', 'FC5', 'T7'};
%badlabels = {'Fz','FC3', 'C4','FT7', 'FT8','F1','F3','P8','CPz','F4'};
badlabels = {'PO7', 'PO3', 'FC3', 'FC1', 'AF3', 'F3'};

nEEG   = numel(OrigEEGchan);
nBad   = numel(badlabels);
badFrac = nBad / max(nEEG,1);
fprintf('EEG bad channels: %d/%d (%.1f%%)\n', nBad, nEEG, 100*badFrac);

% Bad channel interpolation/ removal
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

% Re-referencing
% Do NOT interpolate reference electrodes you'll use for rereferencing (e.g., TP9/TP10). 
% If they are bad, pick a different reference strategy (e.g., average reference after repairing the rest).

% mastoid rereference: offline linked-ears reference
% haveTP9  = any(strcmp(data_nobad.label,'TP9'));
% haveTP10 = any(strcmp(data_nobad.label,'TP10'));
% 
% % % Visual inspection of mastoid quality using ft_databrowser
% cfg = [];
% cfg.viewmode = 'vertical'; % 'vertical' or 'butterfly' view
% cfg.continuous = 'yes'; % 'yes' for continuous data
% cfg.blocksize = 10; % Set the blocksize (number of seconds per window)
% cfg.channel = {'TP9','TP10'}; % 'all' or specify channels of interest
% ft_databrowser(cfg, data_nobad);


% --- Version 1: Mastoid reference (TP9/TP10) ---
cfg = [];
cfg.reref      = 'yes';
cfg.refmethod  = 'avg';
cfg.implicitref = 'FCz';
cfg.refchannel = {'TP9','TP10'};
cfg.channel    = 'all';
data_ref = ft_preprocessing(cfg, data_nobad);
fprintf('Mastoid (TP9/TP10) reference computed\n');


% read events once (same for both)
evt = ft_read_event(fullfile(inputpath,inputname));
codes = {evt.type};
trig_code = codes{3};

% %% different trigger definition for P[AM]_FFR_sit_333 for control_no-plate
%inspect the triggers
% samples = [evt.sample];
% xx = diff(samples);
% plot(xx)
% 
% trig_samp = samples(3:end);
% fsample = data.fsample;
% 
% timewin = [-70, 180];    % Prestimulus and poststimulus (amount of ms to include before and after trigger)
% trlshift_files = 0;      % ms between trigger and stimulus
% 
% % Manual trial definition
% trl = [trig_samp'+(trlshift_files+timewin(1))*(fsample/1000), ...
%        trig_samp'+(trlshift_files+timewin(2))*(fsample/1000), ...
%        repmat((-10*fsample)/1000, fliplr(size(trig_samp)))];
% 
% conditionlabels = repmat({'notrials'}, size(trig_samp));
% 
% for i = size(xx,1)-1:-1:1
%     if xx(i)<1.0000e-03*fsample
%         trl(i+1,:) = [];
%         conditionlabels(i+1) = [];
%     end
% end
% 
% cfg = [];
% cfg.trl = trl;
% data_segm = ft_redefinetrial(cfg, data_ref_mastoid);
% 
% data_segm.trialinfo = conditionlabels';


% Define trials
cfg = [];
cfg.dataset = fullfile(inputpath,inputname);
cfg.trialdef.eventtype = trig_code;
cfg.trialdef.prestim   = 0.01;
cfg.trialdef.poststim  = 0.240;
cfg_trl = ft_definetrial(cfg);
    
data_segm = ft_redefinetrial(cfg_trl, data_ref);
    
% Assign condition labels
n_trigs = size(cfg_trl.trl, 1);
%n_trigs = size(trl, 1); % only for AM sit333-control
conditionlabels = repmat({'notrials'}, n_trigs, 1);
    
for i = 2:1501
    trial_num_in_block = i - 1;
     if mod(trial_num_in_block, 2) == 1
         conditionlabels{i} = 'FFR';
     else
         conditionlabels{i} = '-FFR';
     end
end

for i = 1503:3002
    trial_num_in_block = i - 1502;
    if mod(trial_num_in_block, 2) == 1
        conditionlabels{i} = 'FFR';
    else
        conditionlabels{i} = '-FFR';
    end
end
    
data_segm.trialinfo = conditionlabels;

% Baseline + automatic artifact rejection
data_bl = data_segm;
    
cfg = [];
cfg.latency = [0.030 0.230];
data_toi = ft_selectdata(cfg, data_bl);
    
ntr = numel(data_toi.trial);
trial_var = zeros(ntr,1);

for j = 1:ntr
    Z = data_toi.trial{j};
    var_ch = var(Z, 0, 2);
    trial_var(j) = median(var_ch);
end
    
robust_c = median(trial_var);
robust_s = mad(trial_var, 1);
z_trial  = (trial_var - robust_c) ./ max(robust_s, eps);
rej  = abs(z_trial) > 4;
keep = find(~rej);
    
cfg = [];
cfg.trials = keep;
cfg.keeptrialinfo = 'yes';
data_clean = ft_selectdata(cfg, data_bl);
    
n_before = numel(data_bl.trial);
n_after  = numel(data_clean.trial);

fprintf('[%s] Trials: %d -> %d (%.1f%% rejected)\n', 'mastoid', n_before, n_after, ...
        100*(n_before-n_after)/max(n_before,1));

    
% Manual visual artifact rejection
cfg = [];
cfg.channel = 'all';
cfg.metric  = 'var';
cfg.method  = 'summary';
data_final = ft_rejectvisual(cfg, data_clean);
    
% Log to Excel (one logfile per reference)
logpath = 'C:\Users\Aline Moorkens\OneDrive - KU Leuven\PhD\auditory-motor_project\subjects';
logfile = fullfile(logpath, 'FFR_preprocessing_333_batch2.xlsx');

suggested_bad = data_bp.label(bad_idx);
suggested_bad = intersect(suggested_bad, OrigEEGchan);
 
if exist('excluded','var') && ~isempty(excluded)
    affected_chans = excluded;
elseif nBad > 0
    affected_chans = badlabels;
else
    affected_chans = {};
end
 
n_final     = numel(data_final.trial);
n_rej_total = n_before - n_final;
pct_rej_tot = 100 * n_rej_total / max(n_before,1);

newRow = table( ...
     string(subjectid), string(name), ...
     string(strjoin(suggested_bad, ', ')), numel(suggested_bad), ...
     string(strjoin(affected_chans, ', ')), numel(affected_chans), ...
     n_before, n_rej_total, n_final, pct_rej_tot, ...
     'VariableNames', {'Subject','Condition', ...
                       'SuggestedBadChans','NumSuggestedBad', ...
                       'AffectedChans','NumAffected', ...
                       'TrialsBefore','TrialsRejected','TrialsAfter','PctRejected'});

sheet = 'Log';

if isfile(logfile)
     existing = readtable(logfile, 'Sheet', sheet, 'TextType','string');
     existing = existing(~ismember(existing.Subject, ["AVERAGE","TOTAL"]), :);
     combined = [existing; newRow];
else
    combined = newRow;
end

writetable(combined, logfile, 'Sheet', sheet, 'WriteMode','overwritesheet');

nRows     = height(combined);
firstData = 2;
lastData  = firstData + nRows - 1;
avgRow    = lastData + 2;
totRow    = lastData + 3;

writecell({'AVERAGE'}, logfile, 'Sheet', sheet, 'Range', sprintf('A%d', avgRow));
    
for col = ["D","F","G","H","I","J"]
     writecell({sprintf('=AVERAGE(%s%d:%s%d)', col, firstData, col, lastData)}, ...
         logfile, 'Sheet', sheet, 'Range', sprintf('%s%d', col, avgRow));
end

writecell({'TOTAL'}, logfile, 'Sheet', sheet, 'Range', sprintf('A%d', totRow));

for col = ["F","G","H","I"]
     writecell({sprintf('=SUM(%s%d:%s%d)', col, firstData, col, lastData)}, ...
         logfile, 'Sheet', sheet, 'Range', sprintf('%s%d', col, totRow));
end
    
% Save preprocessed data
save(fullfile(outputpath, ['preprocessed_', outputname]), 'data_final', '-v7.3');
    
    
% Average per condition and save
cfg = []; cfg.trials = strcmp(data_final.trialinfo, 'FFR');
dataFFR = ft_selectdata(cfg, data_final);
cfg = []; cfg.trials = strcmp(data_final.trialinfo, '-FFR');
dataMinFFR = ft_selectdata(cfg, data_final);
    
cfg = []; cfg.trials = 'all';
avgFFR    = ft_timelockanalysis(cfg, dataFFR);
avgMinFFR = ft_timelockanalysis(cfg, dataMinFFR);
    
final_avg = (avgFFR.avg - avgMinFFR.avg)/2;

final = avgFFR;
final.avg = final_avg;

save(fullfile(outputpath, ['averagedFFR_', outputname]), 'final', '-v7.3');
    