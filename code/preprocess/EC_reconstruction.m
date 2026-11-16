%% Restore sit_2 file for subject EC
clc
clear
close all

% addpath 'C:\Users\Aline Moorkens\OneDrive - KU Leuven\Documents\fieldtrip-20250106'
% ft_defaults

%% Load data
subjectid = 'EC';

% input
inputpath = fullfile('F:\Auditory-mov-data\',subjectid,'eeg'); % in case of normal data with pressure plate
name = '_FFR_sit_333_2';
subfiles = ['a','b'];

for s = 1:numel(subfiles)
    inputname = ['P[',subjectid,']',name,subfiles(s),'.dat'];
    filenames{s,1} = inputname; 
end


% Load files
cfg            = [];
cfg.dataset    = fullfile(inputpath, filenames{1,1});
cfg.continuous = 'yes';
cfg.channel    = 'all';
dummy1      = ft_preprocessing(cfg);


cfg            = [];
cfg.dataset    = fullfile(inputpath, filenames{2,1});
cfg.continuous = 'yes';
cfg.channel    = 'all';
dummy2      = ft_preprocessing(cfg);


%% Inspect trigger
% read events once (same for both)
evt1 = ft_read_event(fullfile(inputpath,filenames{1,1}));
evt2 = ft_read_event(fullfile(inputpath,filenames{2,1}));

