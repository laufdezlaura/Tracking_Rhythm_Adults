%% Preprocessing script for the RhythAdults_1 data
% Author: Laura Fernández-Merino 
% Last change date: 3/04/23
% Matlab version: R2021a on Mac
% Fieldtrip version: 20211001
% This script contains information about preprocessing data from the project RhythAdults.  
% NOTE: All the changes made to the script can be found in my github.* 

clear all;clc;

%% Add paths
restoredefaultpath
pathhome='/Users/laura/Desktop/Midterm_Paper/Analisis/Experiment1/DataRhythAdults_corregida';
%pathfldtrip='/Users/laura/Documents/fieldtrip-20211001';
%addpath(genpath(pathhome));
addpath '/Users/laura/Documents/fieldtrip-20211001'
ft_defaults

%% Define variables
subj={'1','2','3','4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '35'};
cond={'S 11','S 12','S 13','S 14','S 15','S 16', 'S 17', 'S 18',...
    'S 21','S 22','S 23','S 24','S 25','S 26', 'S 27', 'S 28',...
    'S 31','S 32','S 33','S 34','S 35','S 36', 'S 37', 'S 38',...
    'S 41','S 42','S 43','S 44','S 45','S 46', 'S 47', 'S 48',...
    'S 51','S 52','S 53','S 54','S 55','S 56', 'S 57', 'S 58',...
    'S 61','S 62','S 63','S 64','S 65','S 66', 'S 67', 'S 68'};

%% Read data
s=33; %select subject
% trial info:
%Spanish Regular: 11:18 (odd melody, even sentence).
%Spanish Irregular: 21:28 (odd melody, even sentence).
%Spanish Mismatch: 31:38 (odd melody, even sentence).
%Basque Regular: 41:48 (odd melody, even sentence).
%Basque Irregular: 51:58 (odd melody, even sentence).
%Basque Mismatch: 61:68 (odd melody, even sentence).
cfg = [];
cfg.dataset = [pathhome,'/subj',subj{1,s},'/0',subj{1,s},'.eeg']; %file directory
cfg.trialdef.eventtype='Stimulus'; %name of the trigger channel
cfg.trialdef.eventvalue=cond;
cfg.trialdef.prestim=-0.084;
cfg.trialdef.poststim=2.584;
cfg=ft_definetrial(cfg);
cfg.channel    = 'EEG';
cfg.reref='yes';
cfg.refchannel={'A2'};
%cfg.continuous = 'yes';
cfg.lpfilter = 'yes';
cfg.lpfreq   = 30;
dataEEG= ft_preprocessing(cfg);

% Change the label names
Nsens=size(dataEEG.label,1);
for ch=1:Nsens
a=strfind(dataEEG.label{ch},'FP1');
if a == 1
   dataEEG.label{ch}='Fp1';
end
end

for ch=1:Nsens
a=strfind(dataEEG.label{ch},'FP2');
if a == 1
   dataEEG.label{ch}='Fp2';
end
end

for ch=1:Nsens
a=strfind(dataEEG.label{ch},'FZ');
if a == 1
   dataEEG.label{ch}='Fz';
end
end

for ch=1:Nsens
a=strfind(dataEEG.label{ch},'CZ');
if a == 1
   dataEEG.label{ch}='Cz';
end
end

for ch=1:Nsens
a=strfind(dataEEG.label{ch},'PZ');
if a == 1
   dataEEG.label{ch}='Pz';
end
end

%% Resample data
cfg = [];
cfg.resamplefs = 200;
dataEEG_rsmpl= ft_resampledata(cfg,dataEEG);

%% Compute ICA using FastICA
cfg = [];
cfg.channel= 'EEG';
cfg.method = 'fastica';
cfg.numcomponent=20;
comp = ft_componentanalysis(cfg,dataEEG_rsmpl);

%Plot ICA component dynamics
cfg = [];
cfg.viewmode = 'component';
cfg.continuous = 'yes';
cfg.layout   = 'easycapM1.lay';
ft_databrowser(cfg,comp);

%% Remove ica components
cfg = [];
cfg.component = [10 11 16];
cfg.demean = 'no';
dataEEG_ica= ft_rejectcomponent(cfg,comp,dataEEG_rsmpl);

%% Reject EEG sensors when the mean z-score exceeded a value of 3
cfg = [];
cfg.method = 'summary';
cfg.metric = 'zvalue';
cfg.keepchannel = 'yes';
cfg.channel = 'EEG';
data_clean = ft_rejectvisual(cfg, dataEEG_ica);

%% Save the data
save([pathhome,'/subj',subj{1,s},'/',subj{1,s},'_data_clean'],'data_clean');