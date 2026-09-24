%% Compare sit pre- and post-session for power, SNR and PLV
clc, clear
close all

%% Set parameters
subjectid = 'TH';
alpha = 0.05;
conditions = {'sit pre-session', 'sit post-session'};


%% ---------- Power ----------
%%% 1. Input
inputpath = fullfile('F:\Auditory-mov-data\control-no_plate\',subjectid,'preprocessed');

nameSit1 = '_FFR_sit_333';
nameSit2 = '_FFR_sit_333_2';

inputSit1 = ['averagedFFR_','P[',subjectid,']',nameSit1,'.mat'];
inputSit2 = ['averagedFFR_','P[',subjectid,']',nameSit2,'.mat'];

input_sigchans_sit1 = ['SigChans_','P[',subjectid,']',nameSit1,'.mat'];
input_sigchans_sit2 = ['SigChans_','P[',subjectid,']',nameSit2,'.mat'];


%%% 2. Load data
load(fullfile(inputpath, inputSit1),  'final'); ffr_sit1  = final;
load(fullfile(inputpath, inputSit2),  'final'); ffr_sit2  = final;

load(fullfile(inputpath, input_sigchans_sit1),  'sigChansBoth'); sigchans_sit1  = sigChansBoth;
load(fullfile(inputpath, input_sigchans_sit2),  'sigChansBoth'); sigchans_sit2  = sigChansBoth;

%%% 3. Frequency Analysis (Power Spectrum)
% FFT parameters
fs   = 10000;                     % sampling rate (Hz)
N    = size(ffr_sit2.avg, 2);     % number of samples (from data)
NFFT = ceil(fs/2);                % FFT length -> df ≈ 2 Hz

df        = fs/NFFT;                % frequency resolution (Hz)
win       = hann(N)';               % Hann window [1 x N]
win_power = sum(win .^ 2);          % window energy correction

% One-sided frequency axis
freq_vec = (0 : floor(NFFT/2)) * df;          % [1 x (NFFT/2 + 1)] in Hz
nFreqs   = numel(freq_vec);

data_all   = {ffr_sit1, ffr_sit2};

% Preallocate outputs — cell arrays because channel count may differ
psd_all    = cell(1, 3);
labels_all = cell(1, 3);

for c = 1:numel(conditions)
     avg    = data_all{c}.avg(:,1:2500);       % [nChan x nSamples]
     labels = data_all{c}.label;               % {nChan x 1} channel labels
     nChan  = size(avg, 1);
     psd_cond = zeros(nChan, nFreqs);

    for ch = 1:nChan
        signal = avg(ch, :);        % [1 x N]

        % Step 1: Apply Hann window and compute FFT
        S = fft(signal .* win, NFFT);

        % Step 2: Keep one-sided spectrum
        S = S(1 : floor(NFFT/2) + 1);

        % Step 3: Convert to PSD (µV²/Hz)
        psd_full = (abs(S) .^ 2) / (win_power * fs);

        % Step 4: Double all bins except DC and Nyquist
        psd_full(2:end-1) = 2 * psd_full(2:end-1);

        % Step 5: Store
        psd_cond(ch, 1:numel(psd_full)) = psd_full;
    end

    psd_all{c}    = psd_cond;
    labels_all{c} = labels;         % store channel labels for this condition

    fprintf('PSD computed: %s | %d channels | %d samples\n', ...
        conditions{c}, nChan, N);
end

% Unpack back to named variables
% Unpack to named variables
psd_sit1   = psd_all{1};  labels_sit1  = labels_all{1};
psd_sit2  = psd_all{2};   labels_sit2 = labels_all{2};

fprintf('Frequency resolution: %.2f Hz\n', df);
fprintf('Frequency bins: %d\n', nFreqs);


%%% 4. Channel handling
% Find channels present in both conditions
chan_labels = intersect(labels_sit1,  labels_sit2, 'stable');

fprintf('Channels common to all conditions: %d\n', numel(chan_labels));

% Remove unwanted channels
chans_to_remove = {'TP9', 'TP10'};     % mastoids used for reference
idx             = ismember(chan_labels, chans_to_remove);
chan_labels     = chan_labels(~idx);
nr_chans        = numel(chan_labels);

fprintf('Channels after exclusion: %d\n', nr_chans);

% Find indices of common channels in each condition's label list
idx_sit1  = find(ismember(labels_sit1,  chan_labels));
idx_sit2  = find(ismember(labels_sit2,  chan_labels));

% Apply channel selection to PSD matrices
psd_sit1  = psd_sit1(idx_sit1,:);
psd_sit2  = psd_sit2(idx_sit2,:);


% --> Only include channels with significant FFR
sig_labels = intersect(sigchans_sit1, sigchans_sit2, 'stable');

fprintf('Significant channels common to all conditions: %d\n', numel(sig_labels));

% Find indices of common channels in each condition's label list
idx_sig_sit1  = find(ismember(chan_labels, sig_labels));
idx_sig_sit2  = find(ismember(chan_labels, sig_labels));

% Apply channel selection to PSD matrices
psd_sig_sit1  = psd_sit1(idx_sig_sit1,:);
psd_sig_sit2  = psd_sit2(idx_sig_sit2,:);


%%% 5. Average power across frequency bins within 331-335 Hz
freq_idx = freq_vec >= 331 & freq_vec <= 335;

pow_sit1  = mean(psd_sig_sit1( :, freq_idx), 2);   % [nChan x 1]
pow_sit2  = mean(psd_sig_sit2( :, freq_idx), 2);   % [nChan x 1]


%%% 6. Wilcoxon signed-rank
% Pairwise Wilcoxon
[p_power, ~, ~] = signrank(pow_sit1, pow_sit2);


%%% 7. Boxplot visualisation power comparison
nr_conditions = 2;
nr_chans = size(pow_sit1,1);
bp_power = zeros([nr_chans, nr_conditions]);
bp_power(:,1) = pow_sit1;
bp_power(:,2) = pow_sit2;

sText = 14;
figure;
h = boxplot(bp_power, 'Notch','off', ...
              'Labels', {conditions{1,1},conditions{1,2}});

set(gca,'fontsize',sText);
set(h,{'linew'},{1});

c = findobj(gca,'Tag','Box');  % c(1) = last box, etc.

% Make sure indexing matches 2 conditions:
patch(get(c(2),'XData'), get(c(2),'YData'), [0.9290 0.6940 0.1250], 'FaceAlpha', .26); % Sit pre-session
patch(get(c(1),'XData'), get(c(1),'YData'), [0.8500 0.3250 0.0980], 'FaceAlpha', .26); % Sit post-session


title('FFR power comparison [331–335 Hz]')
ylabel('FFR power')

hold on
parallelcoords(bp_power, 'Color', 0.7*[1 1 1], ...
                   'Linestyle', '-', 'Marker', '.', 'Markersize', 7);

% Significance annotation
yMax   = max(bp_power(:));
yMin   = min(bp_power(:));
yRange = yMax - yMin;
if yRange == 0
    yRange = 1;
end
yStep  = 0.02 * yRange;
startY = yMax + 0.01 * yRange;

ylim([0, yMax + 7*yStep]);

% Sit pre-session (1) vs Sit post-session (2)
if p_power <= alpha
    y1 = startY;
    plot([1 1 2 2], [y1 y1+yStep y1+yStep y1], 'k-', 'LineWidth', 1);
    if p_power <= 0.001
        stars = '***';
    elseif p_power <= 0.01
        stars = '**';
    elseif p_power <= 0.05
        stars = '*';
    else
        stars = '';
    end
    text(1.5, y1 + 1.1*yStep, stars, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment',   'bottom', ...
        'FontSize', sText);
end




%% ---------- SNR -----------
%%% 1. Compute SNR per channel
% SNR = power in signal band / mean power in noise bands
% Signal band: 331-335 Hz (around 333 Hz FFR)
% Noise bands: equal-width flanking bands on either side, avoiding the signal
% Here: [321-325 Hz] and [341-345 Hz] (5 Hz bands, 6 Hz gap from signal)

sig_idx  = freq_vec >= 331 & freq_vec <= 335;

% Flanking noise bands — same width as signal band
noise_lo_idx = freq_vec >= 321 & freq_vec <= 325;
noise_hi_idx = freq_vec >= 341 & freq_vec <= 345;
noise_idx    = noise_lo_idx | noise_hi_idx;

% SNR = mean signal power / mean noise power (per channel)
snr_sit1  = mean(psd_sit1( :, sig_idx), 2) ./ mean(psd_sit1( :, noise_idx), 2);
snr_sit2  = mean(psd_sit2( :, sig_idx), 2) ./ mean(psd_sit2( :, noise_idx), 2);

%%% 2. Wilcoxon signed-rank
% Pairwise Wilcoxon
[p_snr, ~, ~] = signrank(snr_sit1, snr_sit2);


%%% 3. Boxplot visualisation SNR comparison
nr_conditions = 2;
nr_chans = size(snr_sit1,1);
bp_snr = zeros([nr_chans, nr_conditions]);
bp_snr(:,1) = snr_sit1;
bp_snr(:,2) = snr_sit2;

sText = 14;
figure;

h = boxplot(bp_snr, 'Notch','off', ...
              'Labels', {conditions{1,1},conditions{1,2}});

set(gca,'fontsize',sText);
set(h,{'linew'},{1});

% Color boxes (findobj returns them in reverse order of plotting)
c = findobj(gca,'Tag','Box');  % c(1) = last box, etc.

% Make sure indexing matches 3 conditions:
% c(1) -> SIT
% c(2) -> SIT2
patch(get(c(2),'XData'), get(c(2),'YData'), [0.9290 0.6940 0.1250], 'FaceAlpha', .26); % Sit pre-session
patch(get(c(1),'XData'), get(c(1),'YData'), [0.8500 0.3250 0.0980], 'FaceAlpha', .26); % Sit post-session


title('FFR SNR comparison [331–335 Hz]')
ylabel('FFR SNR')
%ylim([0 0.05])

hold on
parallelcoords(bp_snr, 'Color', 0.7*[1 1 1], ...
                   'Linestyle', '-', 'Marker', '.', 'Markersize', 7);

% Significance annotation
yMax   = max(bp_snr(:));
yMin   = min(bp_snr(:));
yRange = yMax - yMin;
if yRange == 0
    yRange = 1;
end
yStep  = 0.02 * yRange;
startY = yMax + 0.01 * yRange;

ylim([0, yMax + 7*yStep]);

% Sit pre-session (1) vs Sit post-session (2)
if p_snr <= alpha
    y1 = startY;
    plot([1 1 2 2], [y1 y1+yStep y1+yStep y1], 'k-', 'LineWidth', 1);
    if p_snr <= 0.001
        stars = '***';
    elseif p_snr <= 0.01
        stars = '**';
    elseif p_snr <= 0.05
        stars = '*';
    else
        stars = '';
    end
    text(1.5, y1 + 1.1*yStep, stars, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment',   'bottom', ...
        'FontSize', sText);
end


%%% 4. SNR in dB
% --> convert SNR to dB
snr_sit1_dB = 10*log10(snr_sit1);
snr_sit2_dB = 10*log10(snr_sit2);

% --> Wilcoxon signed-rank
% Pairwise Wilcoxon
[p_snr_dB, ~, ~] = signrank(snr_sit1_dB, snr_sit2_dB);


% Boxplot visualisation SNR comparison
nr_conditions = 2;
nr_chans = size(snr_sit1_dB,1);
bp_snr_dB = zeros([nr_chans, nr_conditions]);
bp_snr_dB(:,1) = snr_sit1_dB;
bp_snr_dB(:,2) = snr_sit2_dB;

sText = 14;
figure;
% Create a boxplot with notches and custom labels
h = boxplot(bp_snr_dB, 'Notch','off', ...
              'Labels', {conditions{1,1},conditions{1,2}});

set(gca,'fontsize',sText);
set(h,{'linew'},{1});

% Color boxes (findobj returns them in reverse order of plotting)
c = findobj(gca,'Tag','Box');  % c(1) = last box, etc.

% Make sure indexing matches conditions
patch(get(c(2),'XData'), get(c(2),'YData'), [0.9290 0.6940 0.1250], 'FaceAlpha', .26);
patch(get(c(1),'XData'), get(c(1),'YData'), [0.8500 0.3250 0.0980], 'FaceAlpha', .26);

title('FFR SNR comparison [331–335 Hz]')
ylabel('FFR SNR')
%ylim([0 0.05])

hold on
parallelcoords(bp_snr_dB, 'Color', 0.7*[1 1 1], ...
                   'Linestyle', '-', 'Marker', '.', 'Markersize', 7);

% Significance annotation
yMax   = max(bp_snr_dB(:));
yMin   = min(bp_snr_dB(:));
yRange = yMax - yMin;
if yRange == 0
    yRange = 1;
end
yStep  = 0.02 * yRange;
startY = yMax + 0.01 * yRange;

ylim([0, yMax + 7*yStep]);

% Sit pre-session (1) vs Sit post-session (2)
if p_snr_dB <= alpha
    y1 = startY;
    plot([1 1 2 2], [y1 y1+yStep y1+yStep y1], 'k-', 'LineWidth', 1);
    if p_snr_dB <= 0.001
        stars = '***';
    elseif p_snr_dB <= 0.01
        stars = '**';
    elseif p_snr_dB <= 0.05
        stars = '*';
    else
        stars = '';
    end
    text(1.5, y1 + 1.1*yStep, stars, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment',   'bottom', ...
        'FontSize', sText);
end


%% ---------- PLV analysis ----------------
%%% 1. Load PLV table
PLV = load(fullfile('F:\Auditory-mov-data\control-no_plate\plv', 'PLV_table.mat'));
PLV = PLV.T_plv;

subT_sit1 = PLV(PLV.Subject == subjectid & PLV.Condition == 'sit' & PLV.Session == 'pre',:);
subT_sit2 = PLV(PLV.Subject == subjectid & PLV.Condition == 'sit' & PLV.Session == 'post',:);

PLV_sit1.plv = subT_sit1{:,"PLV"}; PLV_sit1.labels = subT_sit1{:,"Channel"};
PLV_sit2.plv = subT_sit2{:,"PLV"}; PLV_sit2.labels = subT_sit2{:,"Channel"};

%%% 2. Only include channels with significant FFR
% --> non-siginificant chans are already left out in the previous step (PLV
% table computation)

%%% 3. Wilcoxon signed-rank
plv_sit1 = PLV_sit1.plv; plv_sit2 = PLV_sit2.plv;

[p_plv, h_plv, stats_plv] = signrank(plv_sit1, plv_sit2);

%%% 4.Boxplot visualisation
nr_conditions = numel(conditions);
nr_chans = size(plv_sit1,1);
bp = zeros([nr_chans, nr_conditions]);
bp(:,1) = plv_sit1;
bp(:,2) = plv_sit2;

sText = 14;
figure;
h = boxplot(bp, 'Notch','off', ...
              'Labels', {conditions{1,1},conditions{1,2}});

set(gca,'fontsize', sText);
set(h,{'linew'},{1});

c = findobj(gca,'Tag','Box');
patch(get(c(2),'XData'), get(c(2),'YData'), [0.4860 0.6940 0.1680], 'FaceAlpha', .26);
patch(get(c(1),'XData'), get(c(1),'YData'), [0.9290 0.6940 0.1250], 'FaceAlpha', .26);

if p_plv < 0.0001
    p_str = '<0.0001';
else
    p_str = sprintf('= %.4f', p_plv);
end

title(['FFR PLV comparison [331–335 Hz] ',subjectid])
ylabel('PLV-value')

hold on
% parallelcoords(bp, 'Color', 0.7*[1 1 1], ...
%                    'Linestyle', '-', 'Marker', '.', 'Markersize', 7);


% Annotate vertex electrodes on boxplot
vertex_labels = {'Cz','C1','C2','CPz','FC1','FC2','CP1','CP2'};
labels = PLV_sit1.labels;

% Logical/numeric mask: 1 where labels matches any entry in vertex_labels, 0 otherwise
mask = double(ismember(labels, vertex_labels));

% Same-size array of category names based on the mask
category = repmat({'non-vertex'}, size(mask));
category(mask == 1) = {'vertex'};


hold on
h = parallelcoords(bp, 'Color', 0.7*[1 1 1], ...
    'Linestyle', '-', 'Marker', '.', 'Markersize', 7);

% h is a vector of line handles, one per row/observation in bp
% Recolor/thicken the ones corresponding to vertex electrodes, and label them
for i = 1:numel(h)
    if mask(i) == 1
        set(h(i), 'Color', [1 0 0], 'LineWidth', 1);   % vertex: red, thick

        % Get this line's x/y data to find where to place the label
        xdata = get(h(i), 'XData');
        ydata = get(h(i), 'YData');

        % Place label just past the last point of the line
        text(xdata(end) + 0.1, ydata(end), labels(i), ...
            'Color', [1 0 0], 'FontSize', 8, 'FontWeight', 'bold');
    end
end

% Significance annotation
yMax   = max(bp(:));
yMin   = min(bp(:));
yRange = yMax - yMin;
if yRange == 0
    yRange = 1;
end
yStep  = 0.02 * yRange;
startY = yMax + 0.01 * yRange;

ylim([0, yMax + 7*yStep]);

% Sit pre-session (1) vs Sit post-session (2)
if p_plv <= alpha
    y1 = startY;
    plot([1 1 2 2], [y1 y1+yStep y1+yStep y1], 'k-', 'LineWidth', 1);
    if p_plv <= 0.001
        stars = '***';
    elseif p_plv <= 0.01
        stars = '**';
    elseif p_plv <= 0.05
        stars = '*';
    else
        stars = '';
    end
    text(1.5, y1 + 1.1*yStep, stars, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment',   'bottom', ...
        'FontSize', sText);
end


% Only vertex electrodes included
% Keep only the vertex electrode rows
bp_vertex = bp(mask == 1, :);
plv_sit1_vertex = bp(:,1); plv_sit2_vertex = bp(:,2);

% Wilcoxon signed-rank for vertex
[p_plv_vertex, h_plv_vertex, stats_plv_vertex] = signrank(plv_sit1_vertex, plv_sit2_vertex);

% Boxplot
sText = 14;

figure;
h = boxplot(bp_vertex, 'Notch','off', ...
    'Labels', {conditions{1,1},conditions{1,2}});

set(gca,'fontsize', sText);
set(h,{'linew'},{1});

c = findobj(gca,'Tag','Box');
patch(get(c(2),'XData'), get(c(2),'YData'), [0.4860 0.6940 0.1680], 'FaceAlpha', .26);
patch(get(c(1),'XData'), get(c(1),'YData'), [0.9290 0.6940 0.1250], 'FaceAlpha', .26);

if p_plv_vertex < 0.0001
    p_str = '<0.0001';
else
    p_str = sprintf('= %.4f', p_plv_vertex);
end

title(['FFR PLV comparison [331–335 Hz] - Vertex electrodes ',subjectid])
ylabel('PLV-value')
hold on
parallelcoords(bp_vertex, 'Color', 0.7*[1 1 1], ...
                    'Linestyle', '-', 'Marker', '.', 'Markersize', 7);

% Significance annotation
yMax   = max(bp_vertex(:));
yMin   = min(bp_vertex(:));
yRange = yMax - yMin;
if yRange == 0
    yRange = 1;
end
yStep  = 0.02 * yRange;
startY = yMax + 0.01 * yRange;

ylim([0, yMax + 7*yStep]);

% Sit pre-session (1) vs Sit post-session (2)
if p_plv_vertex <= alpha
    y1 = startY;
    plot([1 1 2 2], [y1 y1+yStep y1+yStep y1], 'k-', 'LineWidth', 1);
    if p_plv_vertex <= 0.001
        stars = '***';
    elseif p_plv_vertex <= 0.01
        stars = '**';
    elseif p_plv_vertex <= 0.05
        stars = '*';
    else
        stars = '';
    end
    text(1.5, y1 + 1.1*yStep, stars, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment',   'bottom', ...
        'FontSize', sText);
end
