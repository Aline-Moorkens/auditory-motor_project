%% Preprocessing of the FFR data
clc
clear
close all

addpath 'C:\Users\Aline Moorkens\OneDrive - KU Leuven\Documents\fieldtrip-20250106'
ft_defaults

%% Settings

% SubjectID

%%% Batch 1: 333 - 89
% subjectid = 'LVB';    % no sit_2 file
% subjectid = 'AM';     
% subjectid = 'TH';     % 1 trigger too many in foam_333 file, skip the first trigger
% subjectid = 'TF';
% subjectid = 'NSMT';   % no sit_2 file
%
% subjectid = 'MGZ';
% subjectid = 'LH';     % no sit_2 file & 1 trigger missing sit_333 file
% subjectid = 'MM';      
% subjectid = 'FS';     
% subjectid = 'TD';
% subjectid = 'RNR';    % 2 triggers too many at end of sit_333 file
% subjectid = 'HX';     % multiple triggers missing in sit_2 file & 1 trigger missing in foam_333 file
% subjectid = 'VUR';
% subjectid = 'FD'; 
% subjectid = 'EF';     % 4 consecutive triggers missing in sit_2 file (block 1)
% subjectid = 'MCBG';
% subjectid = 'LV';     % no sit_2 file
% subjectid = 'CX';     % 4 triggers missing in sit_2 file
% subjectid = 'NG';
%subjectid = 'BVD';


%%%  Batch 2: 333 - 433
% subjectid = 'RR';     % no sit_2 file & 4 triggers missing in foam file
% subjectid = 'AF';     % multiple triggers missing in sit_2 file
% subjectid = 'LM';     % 2 triggers too many at end of sit_333 file; 1 trigger missing in foam_333
% subjectid = 'ZB';
% subjectid = 'CF';     % 1 trigger missing in foam_333
% subjectid = 'CB';     % 1 trigger missing in sit_333, 4 trig miss block 1
% & 4 trigg miss block 2 of foam_333 & multiple triggers missing in sit_333_2 plus there is also an extra trigger that needs to be removed
% subjectid = 'ECM';    % triggers missing in foam_333
% subjectid = 'ML';     % trigger missing in firm_333
% subjectid = 'NS';
% subjectid = 'TM';
% subjectid = 'LPB';
% subjectid = 'SK';
% subjectid = 'SG';
% subjectid = 'SP';     % triggers missing in foam_333
% subjectid = 'EC';     % sit_2a and sit_2b recording here instead of 1 sit_2 recording (both are approx 1500 trials long)

% Postural condition
name = '_FFR_sit_333';
%name = '_FFR_firm_333';
%name = '_FFR_foam_333';
%name = '_FFR_sit_333_2';

doVisualRejection = true;          % set true to also run ft_rejectvisual (saved separately)

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


% Missing markers?
nMarkers = 3003;
nMissing = 0;

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
% confirmed.

%%% Filtering: zero-phase Butterworth band-pass, 200-900 Hz, order 6
cfg = [];
cfg.bpfilter  = 'yes';
cfg.channel   = OrigEEGchan;    % don't filter trigger channel
cfg.bpfreq    = [200 900];     % or [80 2000] if you want very broad  
cfg.bpfiltord = 6;
cfg.bpfilttype= 'but';
cfg.bpfiltdir = 'twopass';     % zero-phase!
data_bp = ft_preprocessing(cfg, data);


%%% Automatic bad-channel detection (FASTER-style, continuous filtered data)
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
% Final bad-channel list: read from the per-file CSV, or decide once after visual inspection
% The list is stored, so a rerun gives exactly the same result and no list from a
% previous file can be applied by mistake.
if isfile(badFile)
    T = readtable(badFile, 'TextType', 'string', 'Delimiter', ',');
else
    T = table(strings(0,1), strings(0,1), 'VariableNames', {'File', 'BadChannels'});
end
row = T.File == string(outputname);

if any(row)
    bc = T.BadChannels(find(row, 1));
    if ismissing(bc) || strtrim(bc) == ""
        badlabels = {};
    else
        badlabels = cellstr(split(strtrim(bc), ' '))';
    end
    fprintf('Bad channels from %s: %s\n', badFile, strjoin(badlabels, ', '));
else
    cfg = [];
    cfg.viewmode   = 'vertical';
    cfg.continuous = 'yes';
    cfg.blocksize  = 10;
    cfg.channel    = suggested_bad;
    if isempty(cfg.channel), cfg.channel = 'all'; end
    ft_databrowser(cfg, data_bp);
    answer    = input('Final bad channels for this file (space-separated, Enter for none): ', 's');
    badlabels = strsplit(strtrim(answer), ' ');
    badlabels = badlabels(~cellfun(@isempty, badlabels));
    T = [T; {string(outputname), string(strjoin(badlabels, ' '))}];
    writetable(T, badFile);
end

% Mastoids are the primary reference and may never be interpolated
if any(ismember(badlabels, {'TP9', 'TP10'}))
    error(['TP9/TP10 marked bad in %s. The analysis does not allow interpolating the mastoids: ', ...
           'exclude this run (and log the reason).'], outputname);
end
assert(all(ismember(badlabels, OrigEEGchan)), 'Unknown channel label in the bad-channel list.');

nEEG    = numel(OrigEEGchan);
nBad    = numel(badlabels);
badFrac = nBad / max(nEEG, 1);
fprintf('EEG bad channels: %d/%d (%.1f%%)\n', nBad, nEEG, 100*badFrac);


%%% Bad-channel repair: <= 20% interpolate, 20-25% remove, > 25% exclude the run
interp_max_frac  = 0.20;
exclude_max_frac = 0.25;
excluded = {};

if nBad == 0
    data_nobad = data_bp;
elseif badFrac <= interp_max_frac
    cfg = [];
    cfg.layout = templ_file;
    layout = ft_prepare_layout(cfg);

    cfg = [];
    cfg.method  = 'triangulation';
    cfg.layout  = layout;
    cfg.channel = OrigEEGchan;
    neighbours  = ft_prepare_neighbours(cfg, data_bp);

    cfg = [];
    cfg.badchannel = badlabels;
    cfg.neighbours = neighbours;
    cfg.method     = 'average';                 % average of neighbours (report as such in the methods)
    data_nobad = ft_channelrepair(cfg, data_bp);
    fprintf('Interpolated channels: %s\n', strjoin(badlabels, ', '));
elseif badFrac <= exclude_max_frac
    cfg = [];
    cfg.channel = setdiff(data_bp.label, badlabels, 'stable');
    data_nobad  = ft_selectdata(cfg, data_bp);
    excluded    = badlabels;
    warning('Channels removed (not interpolated): %s. Check whether ROI channels are affected.', strjoin(badlabels, ', '));
else
    error('Too many bad channels (%.1f%%): exclude this run.', 100*badFrac);
end


%%% Re-referencing: linked mastoids (primary) and common average (robustness), FCz recovered in both
cfg = [];
cfg.reref       = 'yes';
cfg.refmethod   = 'avg';
cfg.implicitref = 'FCz';
cfg.refchannel  = {'TP9', 'TP10'};
cfg.channel     = 'all';
data_ref = ft_preprocessing(cfg, data_nobad);

cfg = [];
cfg.reref       = 'yes';
cfg.refmethod   = 'avg';
cfg.implicitref = 'FCz';
cfg.refchannel  = 'all';                         % all EEG channels, FCz included
cfg.channel     = 'all';
data_ref_car = ft_preprocessing(cfg, data_nobad);

% Sanity check: FCz must exist and not be flat after re-referencing
iFCz = strcmp(data_ref.label, 'FCz');
assert(any(iFCz) && std(data_ref.trial{1}(iFCz, :)) > 0, 'FCz was not recovered correctly.');


%%% Define trials: -10 to 240 ms around the trigger (trigger = 30 ms before tone onset)
% read events once (same for both)
evt   = ft_read_event(fullfile(inputpath, inputname));
codes = {evt.type};
trig_code = codes{3};                            % check once per batch that this is the stimulus trigger type
fprintf('Trigger event type used: %s (%d events of this type)\n', ...
        string(trig_code), sum(strcmp(codes, trig_code)));

cfg = [];
cfg.dataset            = fullfile(inputpath, inputname);
cfg.trialdef.eventtype = trig_code;
cfg.trialdef.prestim   = 0.01;
cfg.trialdef.poststim  = 0.240;
cfg_trl = ft_definetrial(cfg);

data_segm     = ft_redefinetrial(cfg_trl, data_ref);
data_segm_car = ft_redefinetrial(cfg_trl, data_ref_car);


%%% Polarity labels (NOT robust to missing triggers!! See separate code for this)
n_trigs = size(cfg_trl.trl, 1);
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
data_segm_car.trialinfo = conditionlabels;


%%% Automatic trial rejection (decided once, on the primary reference)data_bl = data_segm;
cfg = [];
cfg.latency = [0.030 0.230];                     % tone period relative to the trigger
data_toi = ft_selectdata(cfg, data_segm);

ntr = numel(data_toi.trial);
trial_var = zeros(ntr, 1);

for j = 1:ntr
    Z = data_toi.trial{j};
    trial_var(j) = median(var(Z, 0, 2));
end
    
robust_c = median(trial_var);
robust_s = mad(trial_var, 1);
z_trial  = (trial_var - robust_c) ./ max(robust_s, eps);
keep    = find(abs(z_trial) <= 4);

    
cfg = [];
cfg.trials        = keep;
cfg.keeptrialinfo = 'yes';
data_clean     = ft_selectdata(cfg, data_segm);       % primary: linked mastoids
data_clean_car = ft_selectdata(cfg, data_segm_car);   % robustness: same trials, average reference

n_before = numel(data_segm.trial);                    % tone epochs only
n_after  = numel(data_clean.trial);
pctKept  = 100 * n_after / max(n_before, 1);
nPerPol  = [sum(strcmp(data_clean.trialinfo, 'FFR')), sum(strcmp(data_clean.trialinfo, '-FFR'))];
fprintf('Trials: %d -> %d (%.1f%% kept); per polarity: %d / %d\n', n_before, n_after, pctKept, nPerPol);

% Runs keeping < 70% of their trials are excluded from the analyses
excludeRun = pctKept < 70;
if excludeRun
    warning('Only %.1f%% of trials kept: this run is EXCLUDED under the 70%% rule.', pctKept);
end


%%% Optional visual rejection (robustness or blinded primary: decision pending)
data_visual = [];
if doVisualRejection
    cfg = [];
    cfg.channel      = 'all';
    cfg.metric       = 'var';
    cfg.method       = 'summary';
    cfg.keepchannel  = 'yes';                        % never remove channels here: the ROI must stay identical
    data_visual = ft_rejectvisual(cfg, data_clean);
    fprintf('Visual rejection: %d -> %d trials\n', n_after, numel(data_visual.trial));
end

%% Log to Excel file
affected_chans = badlabels;
 
newRow = table( ...
    string(subjectid), string(name), ...
    string(strjoin(suggested_bad, ', ')), numel(suggested_bad), ...
    string(strjoin(affected_chans, ', ')), numel(affected_chans), string(strjoin(excluded, ', ')), ...
    nMarkers, nMissing, n_before, n_after, pctKept, nPerPol(1), nPerPol(2), excludeRun, ...
    numel(data_visual), numel(data_visual.trial), sum(strcmp(data_clean.trialinfo, 'FFR')), ...
    sum(strcmp(data_clean.trialinfo, '-FFR')),...
    'VariableNames', {'Subject', 'Condition', 'SuggestedBadChans', 'NumSuggestedBad', ...
                      'BadChans', 'NumBad', 'RemovedChans', 'Markers', 'MissingTriggers', ...
                      'TrialsBefore', 'TrialsAfterAuto', 'PctKept', 'KeptFFR', 'KeptMinFFR', ...
                      'ExcludeRun70', 'VisualRun','TrialsAfterVisual','KeptFFRvisual', 'KeptMinFFRvisual'});

sheet = 'Log';
if isfile(logfile)
    existing = readtable(logfile, 'Sheet', sheet, 'TextType', 'string');
    existing = existing(~ismember(existing.Subject, ["AVERAGE", "TOTAL"]), :);
    existing = existing(~(existing.Subject == string(subjectid) & existing.Condition == string(name)), :);  % rerun replaces the old row
    combined = [existing; newRow];
else
    combined = newRow;
end
writetable(combined, logfile, 'Sheet', sheet, 'WriteMode', 'overwritesheet');

    
%% Save (single trials; equalisation, ROI and spectral measures follow in the analysis script)
save(fullfile(outputpath, ['preprocessed_', outputname]),     'data_clean',     '-v7.3');   % primary
save(fullfile(outputpath, ['preprocessed_CAR_', outputname]), 'data_clean_car', '-v7.3');   % robustness

if doVisualRejection
    save(fullfile(outputpath, ['preprocessed_visual_', outputname]), 'data_visual', '-v7.3');
end


%% Average: FFR = (mean(+) - mean(-)) / 2
cfg = []; cfg.trials = strcmp(data_final.trialinfo, 'FFR');
avgFFR = ft_timelockanalysis(cfg, data_clean);
cfg = []; cfg.trials = strcmp(data_final.trialinfo, '-FFR');
avgMinFFR = ft_timelockanalysis(cfg, data_clean);

final = avgFFR;
final.avg = (avgFFR.avg - avgMinFFR.avg) / 2;

save(fullfile(outputpath, ['averagedFFR_', outputname]), 'final', '-v7.3');
fprintf('[%s] done.\n', outputname);
