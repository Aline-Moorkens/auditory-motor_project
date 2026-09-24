%% Compare standing versus sitting
clc
clear
close all

% % Set parameters
filepath = 'C:\Users\Aline Moorkens\OneDrive - KU Leuven\PhD\Pilot_audiory-motor';
subject = 'Aline';
subfolder = 'preprocessed/mastoid_reref'

file_1 = 'averagedFFR_SIT333'; %sit
file_2 = 'averagedFFR_Stand333'; %stand - foam

% SIT
load(fullfile(filepath,subject,subfolder,file_1));
ffr_sit = final;

% STAND
load(fullfile(filepath,subject,subfolder,file_2));
ffr_stand = final;

% % Frequency Analysis (Power Spectrum): to test if the ffr_sit data has an overall 
% higher power at freq 333 Hz than the ffr_stand data
cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'pow';
cfg.taper = 'hanning';
%cfg.foilim = [331 335];       % band of interest
cfg.pad = 'nextpow2';          % zero-padding to increase frequency resolution
cfg.channel = 'all';

freq_sit = ft_freqanalysis(cfg, ffr_sit);
freq_stand = ft_freqanalysis(cfg, ffr_stand);

% Extract the channel labels
%chan_labels = freq_sit.label;
chan_labels = intersect(ffr_stand.label, ffr_sit.label, 'stable');

% channels to include/exclude
%chans_to_remove = {'TP9', 'TP10'}; % No mastoids
%chans_to_remove = {'TP9', 'TP10', 'FT9', 'FT10', 'T8', 'T7'}; % No pam channels
%chans_to_remove = {'T8', 'T7'}; % bad chans for Acq 1-4 Bob
%chans_to_remove = {'P7', 'FC2', 'CP2', 'CP1','F4'}; % bad chan for acquisition 5&7 Bob
% chans_to_remove = {'FC1', 'FT10', 'PO9', 'TP10', 'TP9'}; % Mansoureh
chans_to_remove = {'Fp1', 'Fp2'}; % Neurocom

% Find the indices of the channels to be removed
idx = ismember(chan_labels, chans_to_remove);
% Remove the channels using logical indexing
chans_adj = chan_labels(~idx);

nr_chans = size(chans_adj,1);
nr_conditions = 2;

% Remove the channels that are not plotted
cfg = [];
cfg.channel = chans_adj;
freq_sit = ft_selectdata(cfg, freq_sit);
freq_stand = ft_selectdata(cfg,freq_stand);


% Average across frequency bins within 331–335 Hz
freq_idx = freq_sit.freq >= 331 & freq_sit.freq <= 335;
pow_sit = mean(freq_sit.powspctrm(:, freq_idx), 2);  % [chan x 1]
pow_stand = mean(freq_stand.powspctrm(:, freq_idx), 2);


%% Visualisation

% %%  Visual comparison for a representative channel
% cfg         = [];
% cfg.output  = 'pow';
% cfg.channel = EEG_chans(29);
% cfg.method  = 'mtmfft';
% cfg.taper   = 'hanning';
% cfg.foi     = 200:400;
% 
% spectr_stand  = ft_freqanalysis(cfg, FFR_stand);
% 
% spectr_sit = ft_freqanalysis(cfg, FFR_sit);
% 
% 
% figure;
% hold on;
% plot(spectr_stand.freq, (spectr_stand.powspctrm), 'linewidth', 2)
% plot(spectr_sit.freq, (spectr_sit.powspctrm), 'linewidth', 2)
% legend('stand', 'sit')
% xlabel('Frequency (Hz)')
% ylabel('Power (\mu V^2)')
% title('Power spectrum channel Fz in frequency range 200-400 Hz')

% % Paired t-test or Wilcoxin rank test
% Check for normality of the data
figure()
subplot(1,2,1)
qqplot(pow_sit)
title('Q-Q plot power during SITTING')
subplot(1,2,2)
qqplot(pow_stand)
title('Q-Q plot power during STANDING')

% Shapiro_wilk test
[H0_sit, psw_sit, W_sit] = swtest(pow_sit, 0.05);
[H0_stand, psw_stand, W_stand] = swtest(pow_stand, 0.05);

% Get p-value with a paired t-test
% p_val = ttest(stand_power,sit_power);

% Perform a paired t-test --> Only if the power-distribution is
% approximately normal! Not the case here
[~, p_value, ~, stats] = ttest(stand_power, sit_power,"Tail","right");

% Plot the corresponding Student's t-distribution, the returned t-statistic, and the critical t-value. 
% Calculate the critical t-value for the default confidence level of 95% by using tinv.

nu = stats.df;
k = linspace(-15,15,300);
tdistpdf = tpdf(k,nu);
tval = stats.tstat

tvalpdf = tpdf(tval,nu);
tcrit = tinv(0.95,nu)

plot(k,tdistpdf)
hold on
scatter(tval,tvalpdf,"filled")
xline(tcrit,"--")
legend(["Student's t pdf", "t-Statistic", ...
    "Critical Cutoff"])


% h = ttest(x,y) returns a test decision for the null hypothesis that the 
% data in x – y comes from a normal distribution with mean equal to zero and 
% unknown variance, using the paired-sample t-test.

% So the null hypothesis actually states there is no significant difference
% between the mean of these two distributions of values.

% In case of a very skewed ditribution: Use Wilcoxen paired test
% returns the p-value of a paired, two-sided test for the null hypothesis 
% that x – y comes from a distribution with zero median.
%
% right tailed test: the alternate hypothesis states the data in x – y 
% come from a distribution with median greater than 0.
%[p_val_nonparam, h, stats] = signrank(pow_stand,pow_sit);
[p_val_nonparam, h, stats]= signrank(pow_stand,pow_sit,'tail','right');

%  h = 1 indicates a rejection of the null hypothesis, and h = 0 indicates 
% a failure to reject the null hypothesis at the 5% significance level. 

%% Boxplot visualisation
bp = zeros([nr_chans, nr_conditions]);
bp(:,1) = pow_stand;
bp(:,2) = pow_sit;

sText = 14;
figure(1);
% Create a boxplot with notches and custom labels
h = boxplot(bp, 'Notch','off','Labels', {'FFR_STAND_FOAM', 'FFR_SIT'});

% Set the font size for the axes labels
set(gca,'fontsize',sText);

% Set the line width for the boxplot
set(h,{'linew'},{1})

% Find and color the boxplot elements
c = findobj(gca,'Tag','Box');
patch(get(c(2), 'XData'), get(c(2),'YData'),[0.4860 0.6940 0.1680], 'FaceAlpha', .26) % Coloring
patch(get(c(1), 'XData'), get(c(1),'YData'),[0.8500 0.3250 0.0980], 'FaceAlpha', .26) % Coloring

% Add title and label to the y-axis
if p_val_nonparam < 0.0001
    p_str = '<0.0001';
else
    p_str = sprintf('= %.4f', p_val_nonparam);
end

title(['FFR power comparison [331–335 Hz] (p ', p_str, ')'])
ylabel('FFR power')

% Add parallel coordinates for the data points
hold on
parallelcoords(bp, 'Color', 0.7*[1 1 1], 'Linestyle', '-', 'Marker', '.', 'Markersize', 7); % The grey lines connecting the data

% Adjust the y-axis limits if necessary
% ylim([-0.02 1]); % Adjust this range based on your data

% Apply the breaky function to create a broken y-axis
% Specify the ranges to be displayed
% breakyaxis([0.07 0.8],0.005);

%%%%% Create an inset for zoomed view
axes('Position',[.65 .65 .25 .25])
h2 = boxplot(bp, 'Notch', 'off', 'Labels', {'FFR_STAND', 'FFR_SIT'});
ylim([-0.00001 0.00002]); % Adjust this range based on your data
% Set the line width for the boxplot
set(h2,{'linew'},{1})

% Find and color the boxplot elements
c2 = findobj(gca,'Tag','Box');
patch(get(c2(2), 'XData'), get(c2(2),'YData'),[0.4860 0.6940 0.1680], 'FaceAlpha', .26) % Coloring
patch(get(c2(1), 'XData'), get(c2(1),'YData'),[0.8500 0.3250 0.0980], 'FaceAlpha', .26) % Coloring
% Add title and label to the y-axis
%title('FFR comparison')
ylabel('FFR power')
hold on
parallelcoords(bp, 'Color', 0.7*[1 1 1], 'Linestyle', '-', 'Marker', '.', 'Markersize', 7); % The grey lines connecting the data
