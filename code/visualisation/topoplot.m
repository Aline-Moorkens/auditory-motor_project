%% FFR Scalp Topography Plot
% Replicates the structure of Bidelman (2014) Figure 1:
% Panel A: FFR waveform at Fz with marked peaks and troughs
% Panel B: 2D scalp topographies for averaged positive and negative phases
%
% Requires: preprocessed and averaged FFR data (the 'final' structure
% saved at the end of the preprocessing script)
% FieldTrip must be on the path

clc
clear
close all

%% 1. LOAD DATA
subjectid = 'LVB';

% input
inputpath = fullfile('F:\Auditory-mov-data',subjectid,'preprocessed');
%name = '_FFR_sit_333.mat';
%name = '_FFR_firm_333.mat';
name = '_FFR_foam_333.mat';
inputname = ['P[',subjectid,']',name];

input = fullfile(inputpath, ['averagedFFR_', inputname]);

load(input)

templ_file = 'easycapM11.mat';

%% 2. EXTRACT Fz WAVEFORM AND TIME VECTOR
% Find Fz channel index
fz_idx = find(strcmp(final.label, 'Fz'));
if isempty(fz_idx)
    error('Fz channel not found. Check channel labels.');
end

time     = final.time;                  % [1 x nSamples] in seconds
waveform = final.avg(fz_idx, :);        % [1 x nSamples] FFR at Fz

% Restrict to the response window (after stimulus onset, exclude onset artifact)
t_start = 0.030;    % 30 ms post-stimulus: skip onset artifact
t_end   = 0.250;    % end of stimulus

t_mask   = time >= t_start & time <= t_end;
% t_plot   = time(t_mask) * 1000;         % convert to ms for plotting
t_plot = time*1000;
% w_plot   = waveform(t_mask);
w_plot = waveform;

%% 3. DETECT PEAKS AND TROUGHS IN THE Fz WAVEFORM
% These correspond to the starred time points in Bidelman Fig 1A

% MinPeakProminence threshold: adjust based on your data amplitude
min_prom = 0.8 * max(abs(w_plot));   % 80% of max amplitude

[pos_vals, pos_locs] = findpeaks(w_plot, 'MinPeakProminence', min_prom);
[neg_vals, neg_locs] = findpeaks(-w_plot, 'MinPeakProminence', min_prom);
neg_vals = -neg_vals;   % restore sign

% Time points of peaks and troughs in ms
t_pos = t_plot(pos_locs);
t_neg = t_plot(neg_locs);

fprintf('Detected %d positive peaks and %d negative troughs\n', ...
    numel(pos_locs), numel(neg_locs));

%% 4. COMPUTE AVERAGE SCALP MAPS AT PEAK AND TROUGH TIME POINTS
%     Average across all detected positive peaks -> positive phase map
%     Average across all detected negative troughs -> negative phase map

% Convert peak/trough times back to sample indices in the full time vector
fs = 10000;   % sampling frequency in Hz

pos_samples = dsearchn(time', (t_pos/1000)');   % nearest sample indices
neg_samples = dsearchn(time', (t_neg/1000)');

% Average across all channels at each peak/trough sample
% final.avg is [nChan x nSamples]
% pos_map = mean(final.avg(:, pos_samples), 2);   % [nChan x 1]
% neg_map = mean(final.avg(:, neg_samples), 2);   % [nChan x 1]
pos_map = mean(final.avg(:, pos_locs), 2);   % [nChan x 1]
neg_map = mean(final.avg(:, neg_locs), 2);   % [nChan x 1]

%% 5. BUILD FIELDTRIP TIMELOCK STRUCTURES FOR TOPOPLOTS
topo_pos        = final;
topo_pos.avg    = pos_map;
topo_pos.time   = 0;        % single time point
topo_pos.dimord = 'chan_time';

topo_neg        = final;
topo_neg.avg    = neg_map;
topo_neg.time   = 0;
topo_neg.dimord = 'chan_time';

%% 6. DETERMINE COLOUR SCALE (symmetric, shared across both maps)
clim_val = max(abs([pos_map; neg_map]));
clim     = [-clim_val clim_val];

%% 7. PLOT
% ---- Colour scheme ----
col_pos  = [0.85 0.15 0.15];   % red for positive peaks
col_neg  = [0.15 0.35 0.75];   % blue for negative troughs
col_wave = [0.15 0.15 0.15];   % near-black for waveform

%  Figure 1: Fz waveform with peaks and lows
figure('Color', 'w', 'Position', [100 100 1100 500]);
ax_wave = gca;   % capture the current axes handle

hold on;

% Stimulus onset dotted line
xline(t_start*1000, 'k:', 'LineWidth', 1.2);

% Zero amplitude line
yline(0, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8);

% Waveform
plot(t_plot, w_plot, 'Color', col_wave, 'LineWidth', 1.4);

% Mark positive peaks with filled circles
scatter(t_pos, pos_vals, 60, col_pos, 'filled', 'Marker', 'o');

% Mark negative troughs with filled circles
scatter(t_neg, neg_vals, 60, col_neg, 'filled', 'Marker', 'o');

xlabel('Time (ms)', 'FontSize', 11);
ylabel('Amplitude (\muV)', 'FontSize', 11);
title('(A)  FFR at Fz', 'FontSize', 12, 'FontWeight', 'bold');

% Extend x-axis slightly to include pre-stimulus
xlim([time(1)*1000 t_end*1000]);
set(ax_wave, 'FontSize', 10, 'Box', 'off', 'TickDir', 'out');

% Legend
legend({'Stimulus onset', '', 'FFR waveform', 'Positive peaks', 'Negative troughs'}, ...
    'Location', 'northeast', 'FontSize', 8, 'Box', 'off');


%  Figure 2 left: Positive phase topography
figure('Color', 'w', 'Position', [100 100 1100 500]);
subplot(1, 2, 1);
ax_pos = gca;

cfg                  = [];
cfg.layout           = templ_file;
%cfg.channel         = {'all', '-F4', '-F5', '-FC2','-O2','-C3','-PO3','-P6','-C2','-FC3','-PO7'};
cfg.marker           = 'on';
cfg.colorbar         = 'yes';
cfg.zlim             = clim;
%cfg.colormap         = 'RdBu';          % red-blue diverging
%cfg.colormap         = 'magma';
cfg.colormap         = 'plasma';          % red-blue diverging
%cfg.style            = 'both';          % contour lines + colour fill
cfg.style            = 'straight';      % colormap only
cfg.comment          = 'no';
cfg.title            = '(B)  Positive phase';
cfg.figure           = ax_pos;
ft_topoplotER(cfg, topo_pos);

title('(B)  Positive phase', 'FontSize', 12, 'FontWeight', 'bold');

%  PANEL B right: Negative phase topography
subplot(1, 2, 2);
ax_neg = gca;

cfg.title  = 'Negative phase';
cfg.figure = ax_neg;
ft_topoplotER(cfg, topo_neg);

title('Negative phase', 'FontSize', 12, 'FontWeight', 'bold');

%  8. GLOBAL TITLE AND SAVE

sgtitle(sprintf('FFR Scalp Topography — Subject %s — 333 Hz', subjectid), ...
    'FontSize', 13, 'FontWeight', 'bold');

% % Save figure
% saveas(gcf, fullfile(outputpath, sprintf('FFR_topography_%s_333Hz.png', subjectid)));
% fprintf('Figure saved to: %s\n', outputpath);