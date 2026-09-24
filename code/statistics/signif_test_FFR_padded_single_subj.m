%% Test the significance of the FFR data with Rayleigh test
clc
clear
close all

%% User input - participants and conditions
subjectid = 'HX';
conditions = {'sit', 'firm', 'foam'};
%conditions = {'sit'};
subfolder = 'preprocessed';

basepath   = 'F:\Auditory-mov-data';
stim_freq  = 333;
fs         = 10000;
alpha      = 0.05;

% Toolbox path (Circular Statistics)
addpath('C:\Users\Aline Moorkens\OneDrive - KU Leuven\Documents\CircStat2012a');


%% Filter design (same for all participants/conditions)
lower_freq = stim_freq - 5;
upper_freq = stim_freq + 5;
order = 4;   % now safe to crank up
[z, p, k] = butter(order, [lower_freq/(fs/2), upper_freq/(fs/2)], 'bandpass');
[sos, g]  = zp2sos(z, p, k);


%% Filter response sanity check (run once)
filt_fig = figure('Name', 'Filter response', 'Position', [200 200 900 600]);

% Magnitude response
[H, F] = freqz(sos, 8192, fs);
H = H * g;
subplot(2,1,1)
plot(F, 20*log10(abs(H)), 'b', 'LineWidth', 1.5); hold on
xline(stim_freq, 'r--', sprintf('%d Hz', stim_freq))
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)')
title(sprintf('Butterworth bandpass  (order %d, %.1f-%.1f Hz)', ...
      order, lower_freq, upper_freq))
xlim([stim_freq-20, stim_freq+20])
ylim([-60 5])
grid on

% Zoomed phase response (for filtfilt this should be ~zero across passband, but worth showing)
subplot(2,1,2)
plot(F, unwrap(angle(H)), 'b', 'LineWidth', 1.5); hold on
xline(stim_freq, 'r--')
xlabel('Frequency (Hz)'); ylabel('Phase (rad)')
title('Phase response (one-pass)')
xlim([stim_freq-20, stim_freq+20])
grid on

saveas(filt_fig, fullfile(basepath, subjectid,'FilterResponse.png'));

%% Loop over the conditions

    for c = 1:numel(conditions)
        
        %%% input
        cond = conditions{c};
        name = ['_FFR_', cond, '_', num2str(stim_freq),'.mat'];
        %name = ['_FFR_', cond, '_', num2str(stim_freq),'_2', '.mat'];
        inputname = ['P[', subjectid, ']', name];
        input = fullfile(basepath, subjectid, subfolder, ['preprocessed_', inputname]);

        load(input)  % expects variable: data_final

        %%% Select EEG channels (drop trigger if present)
        data = data_final;
        cfg = [];
        cfg.channel = {'all', '-Trigger'};
        final_phase_data = ft_selectdata(cfg, data);

        Ntrials1 = numel(find(strcmp(final_phase_data.trialinfo, 'FFR')));
        Ntrials2 = numel(find(strcmp(final_phase_data.trialinfo, '-FFR')));

        if Ntrials1 == 0 || Ntrials2 == 0
            warning('Subject %s, condition %s: missing FFR or -FFR trials. Skipping.', ...
                    subjectid, cond);
            continue
        end

        %%% Padding length
        nTrials       = numel(final_phase_data.trial);
        trial_lengths = cellfun(@(x) size(x,2), final_phase_data.trial);
        min_trial_len = min(trial_lengths);

        min_cycles       = 3;
        pad_len_required = ceil(min_cycles * fs / stim_freq);
        max_pad_allowed  = floor(min_trial_len * 0.4);

        if pad_len_required > max_pad_allowed
            warning(['Required padding (%d samples, %.1f ms) exceeds 40%% of signal. ' ...
                     'Using maximum: %d samples (%.1f ms). Edge effects may be present.'], ...
                    pad_len_required, 1000*pad_len_required/fs, ...
                    max_pad_allowed,  1000*max_pad_allowed/fs);
            pad_len = max_pad_allowed;
        else
            pad_len = pad_len_required;
        end

        %%% Trial-by-trial filtering + Hilbert with reflection padding
        pad_idx_pre  = 1:pad_len;
        pad_idx_post = (min_trial_len - pad_len + 1):min_trial_len;

        nChans       = size(final_phase_data.label, 1);
        data_target  = zeros(nChans, min_trial_len, nTrials);
        phase_target = zeros(nChans, min_trial_len, nTrials);

        for n_tr = 1:nTrials
            trial_data = final_phase_data.trial{n_tr}(:, 1:min_trial_len);
            for n_ch = 1:nChans
                sig_raw    = trial_data(n_ch, :);
                sig_padded = [fliplr(sig_raw(pad_idx_pre)), sig_raw, fliplr(sig_raw(pad_idx_post))];

                % Filter on padded signal
                sig_filt_padded = filtfilt(sos, g, sig_padded);
                % Hilbert on padded filtered signal — KEY: do this BEFORE trimming
                phase_padded = angle(hilbert(sig_filt_padded));

                % Now trim padding off both
                data_target(n_ch, :, n_tr)  = sig_filt_padded(pad_len+1 : pad_len+min_trial_len);
                phase_target(n_ch, :, n_tr) = phase_padded(pad_len+1   : pad_len+min_trial_len);
            end
        end

        %%% Diagnostic plot: filter & padding sanity check 
        diag_chan  = 33;
        diag_trial = 50;

        sig_raw_diag    = final_phase_data.trial{diag_trial}(diag_chan, 1:min_trial_len);
        sig_padded_diag = [fliplr(sig_raw_diag(pad_idx_pre)), sig_raw_diag, fliplr(sig_raw_diag(pad_idx_post))];
        sig_filt_diag = filtfilt(sos, g, sig_padded_diag);
        phase_diag      = angle(hilbert(sig_filt_diag));

        % Time vectors
        t_trial  = final_phase_data.time{diag_trial}(1:min_trial_len) * 1000;
        dt       = 1/fs;
        t_padded = ((1:length(sig_padded_diag)) - pad_len - 1) * dt * 1000;

        diag_fig = figure('Name', sprintf('Filter diagnostic - %s - %s', subjectid, cond), ...
                          'Position', [150 150 1200 800]);

        % Panel 1: raw signal
        subplot(4,1,1)
        plot(t_trial, sig_raw_diag, 'k')
        xlabel('Time (ms)'); ylabel('Amplitude (\muV)')
        title(sprintf('Raw signal - channel %s, trial %d', ...
              final_phase_data.label{diag_chan}, diag_trial))
        xlim([t_trial(1) t_trial(end)])

        % Panel 2: padded + filtered, with pad regions shaded
        subplot(4,1,2)
        plot(t_padded, sig_filt_diag, 'b'); hold on
        yl = ylim;
        patch([t_padded(1) t_padded(pad_len) t_padded(pad_len) t_padded(1)], ...
              [yl(1) yl(1) yl(2) yl(2)], [1 0.8 0.8], 'EdgeColor','none','FaceAlpha',0.4)
        patch([t_padded(end-pad_len+1) t_padded(end) t_padded(end) t_padded(end-pad_len+1)], ...
              [yl(1) yl(1) yl(2) yl(2)], [1 0.8 0.8], 'EdgeColor','none','FaceAlpha',0.4)
        xline(t_padded(pad_len),       'r--', 'LineWidth', 1)
        xline(t_padded(end-pad_len+1), 'r--', 'LineWidth', 1)
        xlabel('Time (ms)'); ylabel('Amplitude')
        title(sprintf('Filtered signal w/ reflection padding (red = padded regions, %d samples each)', pad_len))
        xlim([t_padded(1) t_padded(end)])

        % Panel 3: trimmed filtered signal
        subplot(4,1,3)
        sig_filt_trimmed = sig_filt_diag(pad_len+1 : pad_len+min_trial_len);
        plot(t_trial, sig_filt_trimmed, 'b')
        xlabel('Time (ms)'); ylabel('Amplitude')
        title('Filtered signal after trimming (this is what enters the analysis)')
        xlim([t_trial(1) t_trial(end)])

        % Panel 4: instantaneous phase (trimmed)
        subplot(4,1,4)
        phase_trimmed = phase_diag(pad_len+1 : pad_len+min_trial_len);
        plot(t_trial, phase_trimmed, 'm')
        xlabel('Time (ms)'); ylabel('Phase (rad)')
        title('Instantaneous phase (Hilbert, after trim)')
        xlim([t_trial(1) t_trial(end)])
        ylim([-pi pi])

        sgtitle(sprintf('Filter diagnostic - %s - %s - channel %s', ...
                subjectid, cond, final_phase_data.label{diag_chan}))

        diag_name = fullfile(basepath,subjectid, ['FilterDiag_', inputname(1:end-4), '.png']);
        saveas(diag_fig, diag_name);
        fprintf('Saved filter diagnostic: %s\n', diag_name);

        %%% Split phase by condition
        phase_FFR    = phase_target(:, :, strcmp(final_phase_data.trialinfo, 'FFR'));
        phase_minFFR = phase_target(:, :, strcmp(final_phase_data.trialinfo, '-FFR'));


        %%% Combine polarities (rotate -FFR by 180°) for the visualization Rayleigh test
        e1 = exp(1i * phase_FFR);
        e2 = exp(1i * phase_minFFR);
        e2_aligned     = e2 * exp(1i * pi);
        e_combined     = cat(3, e1, e2_aligned);
        combined_angles = angle(e_combined);

        %%% Rayleigh test on combined data (visualization)
        ix1 = dsearchn(final_phase_data.time{1}', 0.05);
        ix2 = dsearchn(final_phase_data.time{1}', 0.21);
        test_values = combined_angles(:, ix1:ix2, :);

        p_rayleigh = zeros(size(test_values, 1), size(test_values, 2));
        for j = 1:size(test_values, 1)
            for k = 1:size(test_values, 2)
                p_rayleigh(j, k) = circ_rtest(squeeze(test_values(j, k, :)));
            end
        end

        %%% Visualise the p-values
        fig = figure('Name', sprintf('%s - %s', subjectid, cond), ...
                     'Visible', 'on', 'Position', [100 100 1400 900]);
        tiledlayout(8, 8);

        for ch = 1:size(test_values, 1)
            p_values = p_rayleigh(ch, :);
            nexttile(ch)
            plot(p_values, 'bo-', 'LineWidth', 1, 'MarkerSize', 4);
            hold on;
            sig_idx = find(p_values < alpha);
            plot(sig_idx, p_values(sig_idx), 'ro', 'LineWidth', 1, 'MarkerSize', 6);
            xlabel('Sample number');
            ylabel('p-value');
            channel = final_phase_data.label{ch};
            title(['p-value ', channel]);
            xlim([0 length(p_values) + 1]);
            if ch == 1
                legend('p-values', 'Significant p-values', 'Location', 'northwest');
            end
        end
        sgtitle(sprintf('Rayleigh Test p-values - %s - %s', subjectid, cond))


        %%% Save figure
        fig_name = fullfile(basepath, subjectid, ['Rayleigh_', inputname(1:end-4), '.png']);
        saveas(fig, fig_name);
        fprintf('Saved figure: %s\n', fig_name);


        %%% Per-condition Rayleigh tests for sigChans (FFR and -FFR separately)
        condsLocal = {'FFR', 'minFFR'};
        sigChans   = struct();

        for ci = 1:numel(condsLocal)
            switch condsLocal{ci}
                case 'FFR',    phase_cond = phase_FFR;
                case 'minFFR', phase_cond = phase_minFFR;
            end

            nCh = size(phase_cond, 1);
            nT  = ix2 - ix1 + 1;
            p_ray = nan(nCh, nT);

            for ch = 1:nCh
                for ti = 1:nT
                    p_ray(ch, ti) = circ_rtest(squeeze(phase_cond(ch, ix1+ti-1, :)));
                end
            end

            sig_mask = any(p_ray < alpha, 2);
            sigChans.(condsLocal{ci}) = final_phase_data.label(sig_mask);
        end

        sigChansBoth = intersect(sigChans.FFR, sigChans.minFFR);

        %%% Save significant channel names
        output = fullfile(basepath,subjectid,subfolder, ['SigChansPadded_', inputname]);
        save(output, 'sigChans', 'sigChansBoth');
        fprintf('Saved significant channels to:\n%s\n', output);
        disp('Channels significant in BOTH conditions:');
        disp(sigChansBoth);

    end  % end conditions loop
