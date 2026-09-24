%% FFR Single subsample - all conditions, 2000 trials (1000 per polarity)
clc, clear, close all

subjectid  = 'TH';
%conditions = {'sit','firm','foam'};
conditions = {'sit'}; % for 2nd sit file

inputpath  = fullfile('F:\Auditory-mov-data\control-no_plate', subjectid, 'preprocessed');

outputpath = fullfile('F:\Auditory-mov-data\control-no_plate\phase-analysis\plv_subsampled');

N_draw = 1000;
% rng(1);   % uncomment for a reproducible draw

for c = 1:numel(conditions)
    cond  = conditions{c};
    input = ['preprocessed_', 'P[', subjectid, ']', '_FFR_', cond, '_333.mat'];
    %input = ['preprocessed_', 'P[', subjectid, ']', '_FFR_', cond, '_333_2.mat']; % for 2nd sit file
    load(fullfile(inputpath, input), 'data_final');

    % Separate trials by polarity
    pos_idx = find(strcmp(data_final.trialinfo, 'FFR'));
    neg_idx = find(strcmp(data_final.trialinfo, '-FFR'));

    fprintf('%s | Positive polarity: %d trials | Negative polarity: %d trials\n', ...
        cond, numel(pos_idx), numel(neg_idx));

    if numel(pos_idx) < N_draw || numel(neg_idx) < N_draw
        warning('Condition %s: not enough trials to draw %d per polarity — skipping.', cond, N_draw);
        continue
    end

    % Draw at random 1000 trials from each polarity, without replacement
    rand_pos = pos_idx(randperm(numel(pos_idx), N_draw));
    rand_neg = neg_idx(randperm(numel(neg_idx), N_draw));

    % With replacement (alternative)
    % rand_pos = pos_idx(randi(numel(pos_idx), N_draw, 1));
    % rand_neg = neg_idx(randi(numel(neg_idx), N_draw, 1));

    % Select the subsampled trials via FieldTrip, keeping struct consistent
    cfg_pos        = [];
    cfg_pos.trials = rand_pos;
    data_pos        = ft_selectdata(cfg_pos, data_final);

    cfg_neg        = [];
    cfg_neg.trials = rand_neg;
    data_neg        = ft_selectdata(cfg_neg, data_final);

    % Recombine into a single preprocessed structure (2000 trials) ---
    data_sub = ft_appenddata([], data_pos, data_neg);

    % ft_appenddata drops/rebuilds some fields (e.g. trialinfo, cfg) —
    % restore trialinfo explicitly so polarity labels are preserved
    data_sub.trialinfo = [data_pos.trialinfo; data_neg.trialinfo];

    % Save
    outputname = ['subsampled_', 'P[', subjectid, ']', '_FFR_', cond, '_333.mat'];
    %outputname = ['subsampled_', 'P[', subjectid, ']', '_FFR_', cond, '_333_2.mat'];
    save(fullfile(outputpath, outputname), 'data_sub', '-v7.3');
    fprintf('Saved: %s (%d trials)\n\n', outputname, numel(data_sub.trial));

    clear data_final data_pos data_neg data_sub
end