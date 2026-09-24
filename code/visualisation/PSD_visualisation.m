%% PSD_visualisation_vertx.m
%
% DESCRIPTION:
%   Visualises the Frequency-Following Response (FFR) per channel in both
%   the time domain and the frequency domain (Power Spectral Density).
%
% INPUT (must be present in the workspace before running):
%
%   final   - FieldTrip timelock structure (output of ft_timelockanalysis)
%               Required fields:
%                 .avg    [nChannels × nSamples]  trial-averaged EEG (µV)
%                 .label  {nChannels × 1}         channel labels
%                 .time   [1 × nSamples]          time axis (seconds)
%
%   fs      - Sampling rate in Hz (scalar)
%               e.g. fs = 5000;
%
% OUTPUT:
%
%   One figure per channel, each containing:
%     - Top panel:    Time-domain waveform                 (µV vs ms)
%     - Bottom panel: Power Spectral Density (PSD)         (µV²/Hz vs Hz)
%
%   No variables are written to the workspace.
%
% SETTINGS — adjust these before running:

channels_of_interest = {'Fz','Cz','FCz','CPz','C1','C2','FC1','FC2','AFz','CP1','CP2'};
freq_max = 1000;   % Hz — upper limit of the frequency axis in the plot


%% Extract data
chan_idx  = find(ismember(final.label, channels_of_interest));
time_vec  = final.time;   % seconds
N         = length(time_vec);

% Warn if any requested channels are missing
missing = setdiff(channels_of_interest, final.label);
if ~isempty(missing)
    warning('The following channels were not found and will be skipped: %s', ...
        strjoin(missing, ', '));
end

%% FFT parameters
% NFFT is chosen so that the frequency resolution df = fs/NFFT ≈ 2 Hz.
% Zero-padding (NFFT >> N) interpolates the spectrum for a smoother plot
% without changing the underlying frequency resolution of the data.
NFFT      = ceil(fs / 2);          % gives df ≈ 2 Hz
df        = fs / NFFT;             % actual frequency resolution (Hz)

% Hanning window matched to the actual signal length N (not the padded length).
% Applied to the signal before the FFT to suppress spectral leakage caused
% by the sharp edges of the epoch.
win = hann(N)';

% Window power normalisation factor (equivalent noise bandwidth correction).
% Dividing by win_power in the PSD formula corrects for the energy removed
% from the signal by the window, ensuring correct power scaling.
win_power = sum(win.^2);

% Frequency axis and plot mask
freq_vec  = (0 : floor(NFFT/2)) * df;    % one-sided frequency axis (Hz)
f_idx     = freq_vec <= freq_max;        % logical mask to limit plot to freq_max


%% Compute PSD and plot — one figure per channel
for ch = 1:length(chan_idx)

    signal    = final.avg(chan_idx(ch), :);   % [1 × N], µV
    chan_name = final.label{chan_idx(ch)};

    % PSD computation 
    % Step 1: Apply window and compute zero-padded FFT.
    %         Zero-padding to NFFT points gives finer frequency sampling.
    S = fft(signal.*win, NFFT);

    % Step 2: Keep only the one-sided (positive-frequency) spectrum.
    S = S(1 : floor(NFFT/2) + 1);

    % Step 3: Convert to PSD in µV²/Hz.
    %         abs(S).^2  →  raw power (µV²)
    %         / win_power →  correct for window energy loss
    %         / fs        →  normalise for sampling rate, giving µV²/Hz
    psd = (abs(S).^2)/(win_power * fs);

    % Step 4: Double all bins except DC (index 1) and Nyquist (last index).
    %         Because we discarded the negative-frequency half of the FFT,
    %         we must multiply by 2 to conserve total signal power.
    psd(2:end-1) = 2*psd(2:end-1);

    % Plot
    figure('Color','w', 'Position',[100 100 1200 500], ...
           'Name', sprintf('FFR – %s', chan_name));

    % Time domain
    subplot(2,1,1);
    plot(time_vec * 1e3, signal, 'Color','#05802C', 'LineWidth',1.4);
    xlabel('Time (ms)');
    ylabel('Amplitude (µV)');
    title(sprintf('Time Domain – %s', chan_name));
    grid on; box off;

    % Frequency domain (PSD)
    subplot(2,1,2);
    plot(freq_vec(f_idx), psd(f_idx), 'Color','#410580', 'LineWidth',1.4);
    xlabel('Frequency (Hz)');
    ylabel('Power Spectral Density (µV²/Hz)');
    title(sprintf('Frequency Domain – %s', chan_name));
    grid on; box off;

    sgtitle(sprintf('Frequency-Following Response (FFR) – %s', chan_name));

end
