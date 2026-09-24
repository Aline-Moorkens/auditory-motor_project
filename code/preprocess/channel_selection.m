%% save subset of the HX files
close all
clear
clc

%% Set path & load data
subjectid = 'HX';
basepath = fullfile('F:\Auditory-mov-data',subjectid); 
name = '_FFR_foam_333';
inputname = ['preprocessed_P[',subjectid,']',name,'.mat'];
outputname = ['P[',subjectid,']',name];

load(fullfile('F:\Auditory-mov-data\HX\preprocessed_old',inputname));

%% Select channels in previously preprocessed data
chans_remove = {'F1','FC1','TP7','P5','PO8','Oz','CP1','PO3'};
chans_keep = data_final.label(~ismember(data_final.label,chans_remove));

cfg = []; 
cfg.channel = chans_keep;
data_final = ft_selectdata(cfg, data_final);

% Save preprocessed data
save(fullfile('F:\Auditory-mov-data\HX\preprocessed', ['preprocessed_', outputname]), 'data_final', '-v7.3');


%% Average per condition and save
cfg = []; cfg.trials = strcmp(data_final.trialinfo, 'FFR');
dataFFR = ft_selectdata(cfg, data_final);
cfg = []; cfg.trials = strcmp(data_final.trialinfo, '-FFR');
dataMinFFR = ft_selectdata(cfg, data_final);
    
cfg = []; cfg.trials = 'all';
avgFFR    = ft_timelockanalysis(cfg, dataFFR);
avgMinFFR = ft_timelockanalysis(cfg, dataMinFFR);
    
final_avg = (avgFFR.avg - avgMinFFR.avg)/2;

final = avgFFR;
final.avg = final_avg;

save(fullfile('F:\Auditory-mov-data\HX\preprocessed', ['averagedFFR_', outputname]), 'final', '-v7.3');