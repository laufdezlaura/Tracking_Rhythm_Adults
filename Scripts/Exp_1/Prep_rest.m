%% Preprocessing script resting state data for the RhythAdults_1 data
% Author: Laura Fernández-Merino and Mikel Lizarazu
% Last change date: 30/JAN/23
% Matlab version: R2021a on Mac
% Fieldtrip version: 20211001
% This script contains information about preprocessing resting state data from the project RhythAdults.  
% NOTE: All the changes made to the script can be found in my github.* 

%% Add paths
restoredefaultpath
pathhome='/Users/laura/Desktop/Midterm_Paper/Analisis/Experiment1/DataRhythAdults_corregida';
%pathfldtrip='/Users/laura/Documents/fieldtrip-20211001';
%addpath(genpath(pathhome));
addpath '/Users/laura/Documents/fieldtrip-20211001'
ft_defaults

%% Define variables
subj={'1','2','3','4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '35'};
cond={'S  2'}; 
%%
clearvars -except pathhome subj s cond

%% Read data
s=33; %select subject
% trial info:
% S 2 - 5 minutes of resting state
cfg = [];
cfg.dataset = [pathhome,'/subj',subj{1,s},'/0',subj{1,s},'.eeg']; %file directory
cfg.trialdef.eventtype = 'Stimulus'; %name of the trigger channel
cfg.trialdef.eventvalue = cond;
cfg.trialdef.prestim = 0;
cfg.trialdef.poststim = 300; % five minutes of resting state
cfg = ft_definetrial(cfg);
cfg.channel    = 'EEG';
cfg.reref = 'yes';
cfg.refchannel = {'A2'};
%cfg.continuous = 'yes';
cfg.lpfilter = 'yes';
cfg.lpfreq   = 30;
dataEEG = ft_preprocessing(cfg);


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

% Plot ICA component dynamics
cfg = [];
cfg.viewmode = 'component';
cfg.continuous = 'yes';
cfg.layout   = 'easycapM1.lay';
ft_databrowser(cfg,comp);


%% Remove ica components
cfg = [];
cfg.component = [5];
cfg.demean = 'no';
dataEEG_ica= ft_rejectcomponent(cfg,comp,dataEEG_rsmpl);

%% Segment the data
% We need to transform the time window into points. 
data_clean_rest = [];
fs        = 200;
window    = 2.5*fs;
L=size(dataEEG_ica.trial{1},2);
cnt=1;
for i=1:window:L
    data_clean_rest.trial{1,cnt}=dataEEG_ica.trial{1,1}(:,i:(i+window-1));
    data_clean_rest.time{1,cnt}=dataEEG_ica.time{1,1}(:,i:(i+window-1));
    cnt=cnt+1;
end
data_clean_rest.fsample=fs;
data_clean_rest.label=dataEEG_ica.label;

%%% Reject EEG sensors when the mean z-score exceeded a value of 3
fs=200;
   window=2.5*fs; %time window
   noverlap=1.25*fs; %Overlap


   Nsensor=size(data_clean_rest.label,1); %number of sensors
   count=1;
   segm_data=[]; m=[];
   for i=1:size(data_clean_rest.trial,2)
       for m=1:window-noverlap:size(data_clean_rest.trial{1,i},2)-window+1
           for j=1:Nsensor
               n=round(m);
               segm_data.trial{count}(j,:)=data_clean_rest.trial{1,i}(j,n:n+round(window)-1);
               
           end
           segm_data.time{count}(1,:)=1/fs:1/fs:window/fs;
           count=count+1;
       end
   end
   
   data_clean_rest.trial=segm_data.trial;
   data_clean_rest.time=segm_data.time;

%% Save the data
save([pathhome,'/subj',subj{1,s},'/',subj{1,s},'_data_clean_rest'],'data_clean_rest');