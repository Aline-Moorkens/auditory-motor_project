%% Compare standing versus sitting using SNR (331–335 Hz) instead of band power
clc
clear
close all


%% Define input data
filepath  = 'C:\Users\Aline Moorkens\OneDrive - KU Leuven\PhD\Pilot_audiory-motor';
subject   = 'Aline';
subfolder = 'preprocessed/mastoid_reref';

file_1 = 'averagedFFR_SIT333';     % sit
file_2 = 'averagedFFR_Stand333';   % stand - foam

% Channels to exclude
chans_to_remove = {'FT9', 'FT10','F7'};

% Frequency band settings (Hz)
sigBand    = [331 335];   % signal band around 333 Hz
noiseBand1 = [321 329];   % lower sideband noise
noiseBand2 = [337 345];   % upper sideband noise

% Plot text size
sText = 14;

%% Load data
% SIT
load(fullfile(filepath, subject, subfolder, file_1));
ffr_sit = final;

% STAND
load(fullfile(filepath, subject, subfolder, file_2));
ffr_stand = final;

%% Frequency analysis (FFT power spectrum)
cfg = [];
cfg.method  = 'mtmfft';
cfg.output  = 'pow';
cfg.taper   = 'hanning';
cfg.pad     = 'nextpow2';     % zero-padding to increase frequency resolution
cfg.channel = 'all';

freq_sit   = ft_freqanalysis(cfg, ffr_sit);
freq_stand = ft_freqanalysis(cfg, ffr_stand);


%% Match channels + remove unwanted channels
% Use channel intersection to guarantee consistent ordering
chan_labels = intersect(ffr_stand.label, ffr_sit.label, 'stable');

% Remove unwanted channels
idx_rm    = ismember(chan_labels, chans_to_remove);
chans_adj = chan_labels(~idx_rm);

nr_chans      = numel(chans_adj);
nr_conditions = 2;

% Select only wanted channels in both frequency structs
cfg = [];
cfg.channel  = chans_adj;
freq_sit     = ft_selectdata(cfg, freq_sit);
freq_stand   = ft_selectdata(cfg, freq_stand);


%% Compute SNR in the 331–335 Hz band
% Indices for signal band
sig_idx_sit   = freq_sit.freq   >= sigBand(1)   & freq_sit.freq   <= sigBand(2);
sig_idx_stand = freq_stand.freq >= sigBand(1)   & freq_stand.freq <= sigBand(2);

% Indices for noise sidebands
n1_idx_sit   = freq_sit.freq   >= noiseBand1(1) & freq_sit.freq   <= noiseBand1(2);
n2_idx_sit   = freq_sit.freq   >= noiseBand2(1) & freq_sit.freq   <= noiseBand2(2);

n1_idx_stand = freq_stand.freq >= noiseBand1(1) & freq_stand.freq <= noiseBand1(2);
n2_idx_stand = freq_stand.freq >= noiseBand2(1) & freq_stand.freq <= noiseBand2(2);

% Signal power (mean over bins in signal band)
sig_sit   = mean(freq_sit.powspctrm(:, sig_idx_sit), 2);      % [chan x 1]
sig_stand = mean(freq_stand.powspctrm(:, sig_idx_stand), 2);  % [chan x 1]

% Noise power (mean of mean-power in each sideband)
noise_sit = mean([ ...
    mean(freq_sit.powspctrm(:, n1_idx_sit), 2), ...
    mean(freq_sit.powspctrm(:, n2_idx_sit), 2)], 2);

noise_stand = mean([ ...
    mean(freq_stand.powspctrm(:, n1_idx_stand), 2), ...
    mean(freq_stand.powspctrm(:, n2_idx_stand), 2)], 2);

% Guard against division by zero (unlikely, but safe)
eps_val = 1e-20;
noise_sit(noise_sit < eps_val)     = eps_val;
noise_stand(noise_stand < eps_val) = eps_val;

% SNR (linear + dB)
snr_sit_lin   = sig_sit   ./ noise_sit;
snr_stand_lin = sig_stand ./ noise_stand;

snr_sit_dB   = 10*log10(snr_sit_lin);
snr_stand_dB = 10*log10(snr_stand_lin);


%% Normality checks
figure()
subplot(1,2,1)
qqplot(snr_sit_dB)
title('Q-Q plot SNR (dB) during SITTING')

subplot(1,2,2)
qqplot(snr_stand_dB)
title('Q-Q plot SNR (dB) during STANDING')

% Shapiro-Wilk tests (requires swtest in path)
[H0_sit,   psw_sit,   W_sit]   = swtest(snr_sit_dB,   0.05);
[H0_stand, psw_stand, W_stand] = swtest(snr_stand_dB, 0.05);


%% Paired hypothesis test (stand > sit)
% If you want the same nonparametric approach as before:
[p_val_nonparam, h, stats] = signrank(snr_stand_dB, snr_sit_dB, 'tail', 'right');

% (Optional) Parametric paired t-test, if distributions look normal:
% [~, p_value_t, ~, stats_t] = ttest(snr_stand_dB, snr_sit_dB, "Tail","right");


%% Boxplot visualisation
bp = zeros(nr_chans, nr_conditions);
bp(:,1) = snr_stand_dB;
bp(:,2) = snr_sit_dB;

figure(1);
hbox = boxplot(bp, 'Notch','off', 'Labels', {'FFR_STAND', 'FFR_SIT'});

set(gca,'fontsize', sText);
set(hbox, {'linew'}, {1})

% Color boxes
c = findobj(gca,'Tag','Box');
if numel(c) >= 2
    patch(get(c(2), 'XData'), get(c(2), 'YData'), [0.4860 0.6940 0.1680], 'FaceAlpha', .26) % stand
    patch(get(c(1), 'XData'), get(c(1), 'YData'), [0.8500 0.3250 0.0980], 'FaceAlpha', .26) % sit
end

% p-value string
if p_val_nonparam < 0.0001
    p_str = '<0.0001';
else
    p_str = sprintf('= %.4f', p_val_nonparam);
end

title(['FFR SNR comparison [', num2str(sigBand(1)),'–', num2str(sigBand(2)), ' Hz] (p ', p_str, ')'])
ylabel('SNR (dB)')

% Add paired lines (parallel coordinates)
hold on
parallelcoords(bp, 'Color', 0.7*[1 1 1], 'Linestyle', '-', 'Marker', '.', 'Markersize', 7);

