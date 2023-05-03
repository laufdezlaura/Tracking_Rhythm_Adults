%% Preprocessing pipeline for RhythAdults_2
% Author: Laura Fernández-Merino
% Matlab version R2021a 
% Fieldtrip version 20211001
% Last update: 9/Feb/2023. 

% trial info:

%Basque Regular: 71:78 (odd melody, even sentence).
%Basque Irregular: 81:88 (odd melody, even sentence).
%Basque Mismatch: 91:98 (odd melody, even sentence).

clear all;clc;

%% Add paths
restoredefaultpath
pathhome='/Users/laura/Desktop/Midterm_Paper/Analisis/Experiment2/Data_Experiment2_corregida';
%pathfldtrip='/Users/laura/Documents/fieldtrip-20211001';
%addpath(genpath(pathhome));
addpath '/Users/laura/Documents/fieldtrip-20211001'
ft_defaults
%% Define variables
subj={'1','2','3','4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20'};
cond  = {'S 71','S 72','S 73','S 74','S 75','S 76', 'S 77', 'S 78',...
         'S 81','S 82','S 83','S 84','S 85','S 86', 'S 87', 'S 88',...
         'S 91','S 92','S 93','S 94','S 95','S 96', 'S 97', 'S 98'};

  
%% Read data
s=1; %select subject
cfg = [];
cfg.dataset = [pathhome,'/subj',subj{1,s},'/00',subj{1,s},'.eeg']; %file directory
cfg.trialdef.eventtype='Stimulus'; %name of the trigger channel
cfg.trialdef.eventvalue=cond;
cfg.trialdef.prestim=-0.084;
cfg.trialdef.poststim = 3.834; % estos trials son 3.75 + el delay 
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
cfg.component = [1];
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