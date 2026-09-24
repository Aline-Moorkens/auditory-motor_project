%% Test the significance of the FFR data with Rayleigh test
clc
clear
close all

%% Set path
subjectid = 'LVB';

% input
inputpath = fullfile('F:\Auditory-mov-data',subjectid,'preprocessed');
%name = '_FFR_sit_333.mat';
%name = '_FFR_firm_333.mat';
name = '_FFR_foam_333.mat';
inputname = ['P[',subjectid,']',name];

input = fullfile(inputpath, ['preprocessed_', inputname]);

load(input)

%% Rayleigh test
% (Rayleigh test and Inter-trial Phase Coherence look at the same thing, nl:
% uniformity of the phase of the data)

% % Phase calculation
% For this I need the non-averaged data of all channels, the EEG channels

data = data_final;
fs = 10000;

% % Throw away the trigger channel (if still present)
% Specify the channels you want to keep, excluding the one you want to discard
cfg = [];
cfg.channel = {'all', '-Trigger'};
final_phase_data = ft_selectdata(cfg, data);


% Definition of some parameters;
% Nconditions = 2;
% conditions = {'FFR', '-FFR'};
stim_freq = 333;
Ntrials1 = numel(find(strcmp(final_phase_data.trialinfo,'FFR'))); % Number of trials with postive polarity
Ntrials2 = numel(find(strcmp(final_phase_data.trialinfo,'-FFR'))); % Number of trials with negative polarity (180 degrees out of phase with FFR)

% Interval around the frequency of interest
lower_freq = stim_freq-2; 
upper_freq = stim_freq+2;

% Filter definition
% A 4th-order Butterworth bandpass filter is designed to isolate the frequency 
% range centered around the stimulation frequency (stim_freq), with a ±2 Hz window.
order = 4; 
[b, a] = butter(order, [lower_freq/(fs/2), upper_freq/(fs/2)], 'bandpass'); % filter coeff of the bandpass filter

% Zero-phase filtering of the data for both conditions in order to extract
% the activity in the frequency band of interest without distorting the
% phase.
for n_tr = 1:size(final_phase_data.trial,2)
    for n_ch = 1:size(final_phase_data.label,1)
        data_target(n_ch, :, n_tr) = filtfilt(b, a, final_phase_data.trial{n_tr}(n_ch,:));
    end
end

FFR = data_target(:,:,strcmp(final_phase_data.trialinfo,'FFR')); % channels x time x trials for the FFR trials
minFFR = data_target(:,:,strcmp(final_phase_data.trialinfo,'-FFR'));  % channels x time x trials for the minus FFR trials

% The instantaneous phase of the filtered signals is extracted for 
% each condition (FFR and -FFR) using the Hilbert transform.

for i = 1:size(FFR,3)
    %phase1 = mod(angle(hilbert(FFR(:,:,i)')),2*pi);
    phase1 = angle(hilbert(FFR(:,:,i)')); 
    phase_FFR(:,:,i) = phase1';
end

for i = 1:size(minFFR,3)
    %phase2 = mod(angle(hilbert(minFFR(:,:,i)')),2*pi);
    phase2 = angle(hilbert(minFFR(:,:,i)'));
    phase_minFFR(:,:,i) = phase2'; % channels x time x trials
end

e1 = exp(1i*(phase_FFR));
e2 = exp(1i*(phase_minFFR));

% part for combining the phase angles, check this maybe becaus not entirely
% sure about this
e2_aligned = e2 * exp(1i*pi);  % rotate them 180 degrees
e_combined = cat(3, e1, e2_aligned);
combined_angles = angle(e_combined);


%% Rayleigh test 
% (we need the phase/angle in RADIANS for each channel)
addpath('C:\Users\Aline Moorkens\OneDrive - KU Leuven\Documents\CircStat2012a');

% Apply the Rayleigh test per timepoint, over all channels 
% Test is applied to the phase data per timepoint for all the timepoints in the specified interval of the trials.
%
%   [pval, z] = circ_rtest(alpha,w)
%   
%   Computes Rayleigh test for non-uniformity of circular data.
%       H0: the population is uniformly distributed around the circle
%       HA: the populatoin is not distributed uniformly around the circle
%   
%   Assumption: the distribution has maximally one mode and the data is 
%   sampled from a von Mises distribution!
%
%   Input:
%     alpha	sample of angles in radians
%     [w = number of incidences in case of binned angle data]

% Do Rayleigh test for a specific time window (for the steady state part of the FFR response)

%test_values = phase_FFR(:,dsearchn(data.time{1,1}', 0.03):dsearchn(data.time{1,1}', 0.18),:);
%test_values = phase_minFFR(:,dsearchn(data.time{1,1}', 0.03):dsearchn(data.time{1,1}', 0.18),:);

test_values = combined_angles(:, dsearchn(data.time{1}', 0.05):dsearchn(data.time{1}', 0.21), :);

p_rayleigh = zeros(size(test_values,1),size(test_values,2));  % contains the p_values per channel and per timepoint (channels x timepoints)

for j = 1:size(test_values,1)
    for k = 1:size(test_values,2)
    p_rayleigh(j,k) = circ_rtest(squeeze(test_values(j,k,:)));
    end
end


%% Visualise the p-values
% Significance level
alpha = 0.05;
% Create a new figure
figure;
tiledlayout(8,8);

for b = 1:size(test_values,1)

    % Sample p-values from the Rayleigh test
    p_values = p_rayleigh(b,:);

    % Plot the p-values
    %subplot(8,8,b)
    nexttile(b)
    plot(p_values, 'bo-', 'LineWidth', 1, 'MarkerSize', 4); % All p-values in blue
    hold on;
    % Highlight the p-values below the significance level
    significant_indices = find(p_values < alpha);
    plot(significant_indices, p_values(significant_indices), 'ro', 'LineWidth', 1, 'MarkerSize', 6);

    % Add a horizontal line representing the significance level
    %yline(alpha, 'r--', 'LineWidth', 1, 'Label', '\alpha = 0.05', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'bottom');

    % Labels and title
    xlabel('Sample number');
    ylabel('p-value');
    channel = final_phase_data.label{b};
    str = 'p-value  ';
    subplot_title = strcat(str,channel);
    title(subplot_title);

    % Adjust axes limits
    %ylim([0 max(max(p_values), alpha) * 1.1]);
    xlim([0 length(p_values) + 1]);

    if b == 1
    % Add a legend
    legend('p-values', 'Significant p-values', 'Location', 'northwest');
    end
    
    % Add general title
    sgtitle('Rayleigh Test p-values Visualization')

end


%% Collect significant channels for both conditions
alpha = 0.05;  % significance level

% --- Run Rayleigh test for both conditions ---
conds = {'FFR','minFFR'};
sigChans = struct();

for ci = 1:numel(conds)
    which_cond = conds{ci};

    switch which_cond
        case 'FFR'
            phase_cond = phase_FFR;    % channels x time x trials
        case 'minFFR'
            phase_cond = phase_minFFR; % channels x time x trials
    end

    % Time window indices
    ix1 = dsearchn(data.time{1,1}', 0.05);
    ix2 = dsearchn(data.time{1,1}', 0.21);

    % Compute Rayleigh p-values per channel × time
    nCh = size(phase_cond,1);
    nT  = ix2 - ix1 + 1;
    p_rayleigh = nan(nCh,nT);

    for ch = 1:nCh
        for ti = 1:nT
            p_rayleigh(ch,ti) = circ_rtest(squeeze(phase_cond(ch, ix1+ti-1, :)));
        end
    end

    % Mark channels as significant if ANY time point passes alpha
    sig_mask = any(p_rayleigh < alpha, 2);
    sigChans.(which_cond) = final_phase_data.label(sig_mask);
end

% --- Find overlap between the two conditions ---
sigChansBoth = intersect(sigChans.FFR, sigChans.('minFFR'));


%% Save significant channel names
output = fullfile(inputpath, ['SigChans_', inputname]);
save(output, 'sigChans', 'sigChansBoth');
fprintf('Saved significant channels to:\n%s\n', output);
disp('Channels significant in BOTH conditions:');
disp(sigChansBoth);
