%% Coherence pipeline for RhythAdults_2
% Author: Laura Fernández-Merino and Mikel Lizarazu
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
subj= [2:20];

s=1;    %:length(subj)
clearvars -except pathhome subj s
%% Load data
filename=[pathhome,'/subj',num2str(subj(s)),'/',num2str(subj(s)),'_data_clean'];
load(filename,'data_clean');

%% 

% Read all EEG data
cfg=[];
cfg.trials=[find(data_clean.trialinfo(:)==71);find(data_clean.trialinfo(:)==72);find(data_clean.trialinfo(:)==73);...
    find(data_clean.trialinfo(:)==74);find(data_clean.trialinfo(:)==75);find(data_clean.trialinfo(:)==76);...
    find(data_clean.trialinfo(:)==77);find(data_clean.trialinfo(:)==78);...
    find(data_clean.trialinfo(:)==81);find(data_clean.trialinfo(:)==82);find(data_clean.trialinfo(:)==83);...
    find(data_clean.trialinfo(:)==84);find(data_clean.trialinfo(:)==85);find(data_clean.trialinfo(:)==86);...
    find(data_clean.trialinfo(:)==87);find(data_clean.trialinfo(:)==88);...
    find(data_clean.trialinfo(:)==91);find(data_clean.trialinfo(:)==92);find(data_clean.trialinfo(:)==93);...
    find(data_clean.trialinfo(:)==94);find(data_clean.trialinfo(:)==95);find(data_clean.trialinfo(:)==96);...
    find(data_clean.trialinfo(:)==97);find(data_clean.trialinfo(:)==98);];
Sp_Reg=ft_selectdata(cfg,data_clean);

filename = ['/Users/laura/Desktop/Midterm_Paper/Analisis/Experiment2/SentencesOrder_RH2/subj',num2str(subj(s)),'.xlsx'];
[num,txt,raw] =  xlsread(filename);

V=[];
k=1;
for n=1:length(num)
    if k<=size(Sp_Reg.trialinfo,1)
    if num(n)==Sp_Reg.trialinfo(k)
        V(n)=num(n);
        k=k+1;
    else
        V(n)=1;
    end
    end
end
V=V';

Fs=200;Nsens=size(Sp_Reg.label,1);
music=[];
music.fsample=Fs;
music.label=Sp_Reg.label;
music.label{Nsens+1,1}='Audio';
music.trialinfo=V;
Ntime=size(Sp_Reg.time{1,1},2);

k=1;
for v=1:length(V)
    if V(v) == 1
        music.trial{1,v}=zeros(Nsens, Ntime);
        music.time{1,v}=Sp_Reg.time{1,1};
    else 
        music.trial{1,v}=Sp_Reg.trial{1,k};
        music.time{1,v}=Sp_Reg.time{1,1};
        k=k+1;
    end
   
end

WAVS={};cnt=1;
txt_size=size(txt,1);
for tx=2:length(txt)
    k = strfind(txt{tx,1},'\');
    wav=txt{tx,1}(k(end)+1:end);
    WAVS{1,cnt}=wav;
    cnt=cnt+1;
end

for k=1:size(music.trial,2)
    Y=[];yMono_reg =[]; amp=[];
    [Y, Fs]=audioread(WAVS{1,k});
    yMono_reg = sum(Y, 2) / size(Y, 2);
    M=hilbert(yMono_reg);
    amp=abs(M);
    t = 1/Fs:1/Fs:size(Y,1)/Fs;
    audio_file = [];
    S = ['Audio'];
    audio_file.label = cellstr(S);
    audio_file.fsample = Fs;
    audio_file.trial{1,1} = amp';
    audio_file.time{1,1} = t;
    
    % Resample audio data
    audio_file_resample=[];
    cfg            = [];
    cfg.resamplefs = 200;
    cfg.detrend    = 'no';
    audio_file_resample = ft_resampledata(cfg,audio_file);   
    
    music.trial{1,k}(Nsens+1,:)=audio_file_resample.trial{1,1};
    
end

%% Remove bad trials
% Remove the audio label from the original data

selchan=ft_channelselection({'all','-Audio'},music.label);
cfg=[];
cfg.channel=selchan;
music_nonaudio=ft_selectdata(cfg,music);

%% Reject EEG sensors when the mean z-score exceeded a value of 3
cfg = [];
cfg.method = 'summary';
cfg.metric = 'zvalue';
cfg.keepchannel = 'yes';
cfg.channel = 'EEG';
data_clean_music = ft_rejectvisual(cfg, music_nonaudio);

%% copy the index of the clean trials
v_post=data_clean_music.cfg.trials;

% select only the clean trials from the original data
cfg=[];
cfg.trials=v_post;
music_clean=ft_selectdata(cfg,music);

%% Create music conditions separated by lang and rhythm condition

    % Music BQ regular condition

cfg=[];
cfg.trials=[find(music_clean.trialinfo(:)==71);find(music_clean.trialinfo(:)==73);find(music_clean.trialinfo(:)==75);find(music_clean.trialinfo(:)==77)];
   music_reg_basque=ft_selectdata(cfg,music_clean);
   
% Create new trials with overlapping windows

   fs=200;
   window=2.5*fs; %time window
   noverlap=1.25*fs; %Overlap

    Nsensor=size(music_reg_basque.label,1); %number of sensors
   count=1;
   segm_data=[]; m=[];
   for i=1:size(music_reg_basque.trial,2)
       for m=1:window-noverlap:size(music_reg_basque.trial{1,i},2)-window+1
           for j=1:Nsensor
               n=round(m);
               segm_data.trial{count}(j,:)=music_reg_basque.trial{1,i}(j,n:n+round(window)-1);
               
           end
           segm_data.time{count}(1,:)=1/fs:1/fs:window/fs;
           count=count+1;
       end
   end
   
   music_reg_basque.trial=segm_data.trial;
   music_reg_basque.time=segm_data.time;

   % Music BQ irregular condition

cfg=[];
cfg.trials=[find(music_clean.trialinfo(:)==81);find(music_clean.trialinfo(:)==83);find(music_clean.trialinfo(:)==85);find(music_clean.trialinfo(:)==87)];
   music_irreg_basque=ft_selectdata(cfg,music_clean);
 
   % Create new trials with overlapping windows

   fs=200;
   window=2.5*fs; %time window
   noverlap=1.25*fs; %Overlap

   
     Nsensor=size(music_irreg_basque.label,1); %number of sensors
   count=1;
   segm_data=[]; m=[];
   for i=1:size(music_irreg_basque.trial,2)
       for m=1:window-noverlap:size(music_irreg_basque.trial{1,i},2)-window+1
           for j=1:Nsensor
               n=round(m);
               segm_data.trial{count}(j,:)=music_irreg_basque.trial{1,i}(j,n:n+round(window)-1);
               
           end
           segm_data.time{count}(1,:)=1/fs:1/fs:window/fs;
           count=count+1;
       end
   end
   
   music_irreg_basque.trial=segm_data.trial;
   music_irreg_basque.time=segm_data.time;

      % Music BQ mismatching condition

cfg=[];
cfg.trials=[find(music_clean.trialinfo(:)==91);find(music_clean.trialinfo(:)==93);find(music_clean.trialinfo(:)==95);find(music_clean.trialinfo(:)==97)];
   music_mm_basque=ft_selectdata(cfg,music_clean);

 
   % Create new trials with overlapping windows

   fs=200;
   window=2.5*fs; %time window
   noverlap=1.25*fs; %Overlap


   Nsensor=size(music_mm_basque.label,1); %number of sensors
   count=1;
   segm_data=[]; m=[];
   for i=1:size(music_mm_basque.trial,2)
       for m=1:window-noverlap:size(music_mm_basque.trial{1,i},2)-window+1
           for j=1:Nsensor
               n=round(m);
               segm_data.trial{count}(j,:)=music_mm_basque.trial{1,i}(j,n:n+round(window)-1);
               
           end
           segm_data.time{count}(1,:)=1/fs:1/fs:window/fs;
           count=count+1;
       end
   end
   
   music_mm_basque.trial=segm_data.trial;
   music_mm_basque.time=segm_data.time;
   
   
%% Obtain the Fourier Transform of the EEG music data and the audio files

% Basque
Fdata_music_reg_basque=[];Faudio_music_reg_basque=[];count=0;
for i=1:size(music_reg_basque.trial,2)
    for j=1:Nsens
        Fdata_music_reg_basque(j,:,count+1)=fft(music_reg_basque.trial{1,i}(j,:),[],2);
    end
    Faudio_music_reg_basque(:,count+1)=fft(music_reg_basque.trial{1,i}(Nsens+1,:),[],2);
    count=count+1;
end

Fdata_music_irreg_basque=[];Faudio_music_irreg_basque=[];count=0;
for i=1:size(music_irreg_basque.trial,2)
    for j=1:Nsens
        Fdata_music_irreg_basque(j,:,count+1)=fft(music_irreg_basque.trial{1,i}(j,:),[],2);
    end
    Faudio_music_irreg_basque(:,count+1)=fft(music_irreg_basque.trial{1,i}(Nsens+1,:),[],2);
    count=count+1;
end

Fdata_music_mm_basque=[];Faudio_music_mm_basque=[];count=0;
for i=1:size(music_mm_basque.trial,2)
    for j=1:Nsens
        Fdata_music_mm_basque(j,:,count+1)=fft(music_mm_basque.trial{1,i}(j,:),[],2);
    end
    Faudio_music_mm_basque(:,count+1)=fft(music_mm_basque.trial{1,i}(Nsens+1,:),[],2);
    count=count+1;
end


%Compute coherence for the audio 

Faudio_music_reg_basque = permute(Faudio_music_reg_basque,[3 1 2]);
Fxx_reg_basque=[];Fyy_reg_basque=[];Fxy_reg_basque=[];coh_audio_music_reg_basque=[];
Fxx_reg_basque=mean(Fdata_music_reg_basque.*conj(Fdata_music_reg_basque),3);
Fyy_reg_basque=mean(Faudio_music_reg_basque.*conj(Faudio_music_reg_basque),3);
Fxy_reg_basque=mean(Fdata_music_reg_basque.*repmat(conj(Faudio_music_reg_basque),[Nsens 1]),3); 
coh_audio_music_reg_basque(:,:)=Fxy_reg_basque.*conj(Fxy_reg_basque)./(repmat(Fyy_reg_basque,[Nsens 1]).*Fxx_reg_basque);

Faudio_music_irreg_basque = permute(Faudio_music_irreg_basque,[3 1 2]);
Fxx_irreg_basque=[];Fyy_irreg_basque=[];Fxy_irreg_basque=[];coh_audio_music_irreg_basque=[];
Fxx_irreg_basque=mean(Fdata_music_irreg_basque.*conj(Fdata_music_irreg_basque),3);
Fyy_irreg_basque=mean(Faudio_music_irreg_basque.*conj(Faudio_music_irreg_basque),3);
Fxy_irreg_basque=mean(Fdata_music_irreg_basque.*repmat(conj(Faudio_music_irreg_basque),[Nsens 1]),3); 
coh_audio_music_irreg_basque(:,:)=Fxy_irreg_basque.*conj(Fxy_irreg_basque)./(repmat(Fyy_irreg_basque,[Nsens 1]).*Fxx_irreg_basque);

Faudio_music_mm_basque = permute(Faudio_music_mm_basque,[3 1 2]);
Fxx_mm_basque=[];Fyy_mm_basque=[];Fxy_mm_basque=[];coh_audio_music_mm_basque=[];
Fxx_mm_basque=mean(Fdata_music_mm_basque.*conj(Fdata_music_mm_basque),3);
Fyy_mm_basque=mean(Faudio_music_mm_basque.*conj(Faudio_music_mm_basque),3);
Fxy_mm_basque=mean(Fdata_music_mm_basque.*repmat(conj(Faudio_music_mm_basque),[Nsens 1]),3); 
coh_audio_music_mm_basque(:,:)=Fxy_mm_basque.*conj(Fxy_mm_basque)./(repmat(Fyy_mm_basque,[Nsens 1]).*Fxx_mm_basque);


filename=[pathhome,'/Coherence/subj',num2str(subj(s)),'_coherence_music'];
save(filename,'coh_audio_music_reg_basque','coh_audio_music_irreg_basque','coh_audio_music_mm_basque')

%% Create the speech conditions separated by lang and rhythm condition for ALL TRIALS 

% Speech Basque regular condition 

cfg=[];
cfg.trials=[find(music_clean.trialinfo(:)==72);find(music_clean.trialinfo(:)==74);find(music_clean.trialinfo(:)==76);find(music_clean.trialinfo(:)==78)];
   speech_reg_basque=ft_selectdata(cfg,music_clean);
   
   % Create new trials with overlapping windows

   fs=200;
   window=2.5*fs; %time window
   noverlap=1.25*fs; %Overlap


   Nsensor=size(speech_reg_basque.label,1); %number of sensors
   count=1;
   segm_data=[]; m=[];
   for i=1:size(speech_reg_basque.trial,2)
       for m=1:window-noverlap:size(speech_reg_basque.trial{1,i},2)-window+1
           for j=1:Nsensor
               n=round(m);
               segm_data.trial{count}(j,:)=speech_reg_basque.trial{1,i}(j,n:n+round(window)-1);
               
           end
           segm_data.time{count}(1,:)=1/fs:1/fs:window/fs;
           count=count+1;
       end
   end
   
   speech_reg_basque.trial=segm_data.trial;
   speech_reg_basque.time=segm_data.time;

   % Speech Basque irregular condition 
   
cfg=[];
cfg.trials=[find(music_clean.trialinfo(:)==82);find(music_clean.trialinfo(:)==84);find(music_clean.trialinfo(:)==86);find(music_clean.trialinfo(:)==88)];
   speech_irreg_basque=ft_selectdata(cfg,music_clean);
   
    % Create new trials with overlapping windows

   fs=200;
   window=2.5*fs; %time window
   noverlap=1.25*fs; %Overlap


   Nsensor=size(speech_irreg_basque.label,1); %number of sensors
   count=1;
   segm_data=[]; m=[];
   for i=1:size(speech_irreg_basque.trial,2)
       for m=1:window-noverlap:size(speech_irreg_basque.trial{1,i},2)-window+1
           for j=1:Nsensor
               n=round(m);
               segm_data.trial{count}(j,:)=speech_irreg_basque.trial{1,i}(j,n:n+round(window)-1);
               
           end
           segm_data.time{count}(1,:)=1/fs:1/fs:window/fs;
           count=count+1;
       end
   end
   
   speech_irreg_basque.trial=segm_data.trial;
   speech_irreg_basque.time=segm_data.time;

      % Speech Basque mismatching condition 

cfg=[];
cfg.trials=[find(music_clean.trialinfo(:)==92);find(music_clean.trialinfo(:)==94);find(music_clean.trialinfo(:)==96);find(music_clean.trialinfo(:)==98)];
   speech_mm_basque=ft_selectdata(cfg,music_clean);
   
    % Create new trials with overlapping windows

   fs=200;
   window=2.5*fs; %time window
   noverlap=1.25*fs; %Overlap

   Nsensor=size(speech_mm_basque.label,1); %number of sensors
   count=1;
   segm_data=[]; m=[];
   for i=1:size(speech_mm_basque.trial,2)
       for m=1:window-noverlap:size(speech_mm_basque.trial{1,i},2)-window+1
           for j=1:Nsensor
               n=round(m);
               segm_data.trial{count}(j,:)=speech_mm_basque.trial{1,i}(j,n:n+round(window)-1);
               
           end
           segm_data.time{count}(1,:)=1/fs:1/fs:window/fs;
           count=count+1;
       end
   end
   
   speech_mm_basque.trial=segm_data.trial;
   speech_mm_basque.time=segm_data.time;

%% Obtain the Fourier Transform of the EEG speech data and the audio files 

% Basque
Fdata_reg_basque=[];Faudio_reg_basque=[];count=0;
for i=1:size(speech_reg_basque.trial,2)
    for j=1:Nsens
        Fdata_reg_basque(j,:,count+1)=fft(speech_reg_basque.trial{1,i}(j,:),[],2);
    end
    Faudio_reg_basque(:,count+1)=fft(speech_reg_basque.trial{1,i}(Nsens+1,:),[],2);
    count=count+1;
end

Fdata_irreg_basque=[];Faudio_irreg_basque=[];count=0;
for i=1:size(speech_irreg_basque.trial,2)
    for j=1:Nsens
        Fdata_irreg_basque(j,:,count+1)=fft(speech_irreg_basque.trial{1,i}(j,:),[],2);
    end
    Faudio_irreg_basque(:,count+1)=fft(speech_irreg_basque.trial{1,i}(Nsens+1,:),[],2);
    count=count+1;
end

Fdata_mm_basque=[];Faudio_mm_basque=[];count=0;
for i=1:size(speech_mm_basque.trial,2)
    for j=1:Nsens
        Fdata_mm_basque(j,:,count+1)=fft(speech_mm_basque.trial{1,i}(j,:),[],2);
    end
    Faudio_mm_basque(:,count+1)=fft(speech_mm_basque.trial{1,i}(Nsens+1,:),[],2);
    count=count+1;
end


% Compute coherence for the audio 

Faudio_reg_basque = permute(Faudio_reg_basque,[3 1 2]);
Fxx_reg_basque=[];Fyy_reg_basque=[];Fxy_reg_basque=[];coh_audio_reg_basque=[];
Fxx_reg_basque=mean(Fdata_reg_basque.*conj(Fdata_reg_basque),3);
Fyy_reg_basque=mean(Faudio_reg_basque.*conj(Faudio_reg_basque),3);
Fxy_reg_basque=mean(Fdata_reg_basque.*repmat(conj(Faudio_reg_basque),[Nsens 1]),3); 
coh_audio_reg_basque(:,:)=Fxy_reg_basque.*conj(Fxy_reg_basque)./(repmat(Fyy_reg_basque,[Nsens 1]).*Fxx_reg_basque);

Faudio_irreg_basque = permute(Faudio_irreg_basque,[3 1 2]);
Fxx_irreg_basque=[];Fyy_irreg_basque=[];Fxy_irreg_basque=[];coh_audio_irreg_basque=[];
Fxx_irreg_basque=mean(Fdata_irreg_basque.*conj(Fdata_irreg_basque),3);
Fyy_irreg_basque=mean(Faudio_irreg_basque.*conj(Faudio_irreg_basque),3);
Fxy_irreg_basque=mean(Fdata_irreg_basque.*repmat(conj(Faudio_irreg_basque),[Nsens 1]),3); 
coh_audio_irreg_basque(:,:)=Fxy_irreg_basque.*conj(Fxy_irreg_basque)./(repmat(Fyy_irreg_basque,[Nsens 1]).*Fxx_irreg_basque);

Faudio_mm_basque = permute(Faudio_mm_basque,[3 1 2]);
Fxx_mm_basque=[];Fyy_mm_basque=[];Fxy_mm_basque=[];coh_audio_mm_basque=[];
Fxx_mm_basque=mean(Fdata_mm_basque.*conj(Fdata_mm_basque),3);
Fyy_mm_basque=mean(Faudio_mm_basque.*conj(Faudio_mm_basque),3);
Fxy_mm_basque=mean(Fdata_mm_basque.*repmat(conj(Faudio_mm_basque),[Nsens 1]),3); 
coh_audio_mm_basque(:,:)=Fxy_mm_basque.*conj(Fxy_mm_basque)./(repmat(Fyy_mm_basque,[Nsens 1]).*Fxx_mm_basque);


filename=[pathhome,'/Coherence/subj',num2str(subj(s)),'_coherence_speech'];
save(filename,'coh_audio_reg_basque','coh_audio_irreg_basque','coh_audio_mm_basque')