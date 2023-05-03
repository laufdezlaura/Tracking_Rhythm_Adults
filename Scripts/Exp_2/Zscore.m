%% Zscore pipeline for RhythAdults_2
% Author: Laura Fernández-Merino and Mikel Lizarazu
% Matlab version R2021a 
% Fieldtrip version 20211001
% Date of last change: 3/april/23

clear all;clc;

%% Paths
%% Add paths
restoredefaultpath
pathhome='/Users/laura/Desktop/Midterm_Paper/Analisis/Experiment2/Data_Experiment2_corregida';
%pathfldtrip='/Users/laura/Documents/fieldtrip-20211001';
%addpath(genpath(pathhome));
addpath '/Users/laura/Documents/fieldtrip-20211001'
ft_defaults
%% Define variables
subj= [2:20];

s=19;    %:length(subj)
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
speech=ft_selectdata(cfg,music);

%% ZSCORE
coh_real=load([pathhome,'/Coherence/subj',num2str(subj(s)),'_coherence_speech']);

%% SPEECH CONDITION - BASQUE REGULAR

    % speech BQ regular condition

cfg=[];w=[];p=[];data_clean_p=[];
cfg.trials=[find(speech.trialinfo(:)==72);find(speech.trialinfo(:)==74);find(speech.trialinfo(:)==76);find(speech.trialinfo(:)==78)];
speech_reg_basque_all=ft_selectdata(cfg,speech);

% create new trials 

   fs=200;
   window=2.5*fs; %time window
   noverlap=1.25*fs; %Overlap


   Nsensor=size(speech_reg_basque_all.label,1); %number of sensors
   count=1;
   segm_data=[]; m=[];
   for i=1:size(speech_reg_basque_all.trial,2)
       for m=1:window-noverlap:size(speech_reg_basque_all.trial{1,i},2)-window+1
           for j=1:Nsensor
               n=round(m);
               segm_data.trial{count}(j,:)=speech_reg_basque_all.trial{1,i}(j,n:n+round(window)-1);
               
           end
           segm_data.time{count}(1,:)=1/fs:1/fs:window/fs;
           count=count+1;
       end
   end
   
   speech_reg_basque_all.trial=segm_data.trial;
   speech_reg_basque_all.time=segm_data.time;
   
N_reg_basque=size(speech_reg_basque_all.trial,2);  
load([pathhome,'/subj',num2str(subj(s)),'/',num2str(subj(s)),'_data_clean_rest'])
w=1:size(data_clean_rest.trial,2);

for perm=1:500
    rng('shuffle');
    p=w(randperm(length(w)));

    cfg=[];
    cfg.trials=p(1:N_reg_basque);
    data_clean_p=ft_selectdata(cfg,data_clean_rest);
    for trl=1:length(w)
        data_clean_p.trial{1,trl} = data_clean_rest.trial{1,p(trl)};
    end
    Fdata_reg_basque_perm=[];Faudio_reg_basque_perm=[];count=0;
    for i=1:size(speech_reg_basque_all.trial,2)
        for j=1:Nsens
            Fdata_reg_basque_perm(j,:,count+1)=fft(data_clean_p.trial{1,i}(j,:),[],2);
        end
        Faudio_reg_basque_perm(:,count+1)=fft(speech_reg_basque_all.trial{1,i}(Nsens+1,:),[],2);
        count=count+1;
    end
    
    %Compute coherence for the audio
    
    Faudio_reg_basque_perm = permute(Faudio_reg_basque_perm,[3 1 2]);
    Fxx_reg_basque_perm=[];Fyy_reg_basque_perm=[];Fxy_reg_basque_perm=[];
    Fxx_reg_basque_perm=mean(Fdata_reg_basque_perm.*conj(Fdata_reg_basque_perm),3);
    Fyy_reg_basque_perm=mean(Faudio_reg_basque_perm.*conj(Faudio_reg_basque_perm),3);
    Fxy_reg_basque_perm=mean(Fdata_reg_basque_perm.*repmat(conj(Faudio_reg_basque_perm),[Nsens 1]),3);
    coh_speech_reg_basque_all_perm(:,:,perm)=Fxy_reg_basque_perm.*conj(Fxy_reg_basque_perm)./(repmat(Fyy_reg_basque_perm,[Nsens 1]).*Fxx_reg_basque_perm);
    
    
end


Nfreq=500;

for ch=1:Nsens-1
    for f=1:Nfreq
        kk=[];
        kk=zscore([coh_real.coh_audio_reg_basque(ch,f);squeeze(coh_speech_reg_basque_all_perm(ch,f,:))]);
        z_score_coh_speech_reg_basque_all(ch,f)=kk(1);
    end
end

%% SPEECH CONDITION - BASQUE IRREGULAR


cfg=[];w=[];p=[];data_clean_p=[];
cfg.trials=[find(speech.trialinfo(:)==82);find(speech.trialinfo(:)==84);find(speech.trialinfo(:)==86);find(speech.trialinfo(:)==88)];
   speech_irreg_basque_all=ft_selectdata(cfg,speech);
   
% create new trials 

   fs=200;
   window=2.5*fs; %time window
   noverlap=1.25*fs; %Overlap


   Nsensor=size(speech_irreg_basque_all.label,1); %number of sensors
   count=1;
   segm_data=[]; m=[];
   for i=1:size(speech_irreg_basque_all.trial,2)
       for m=1:window-noverlap:size(speech_irreg_basque_all.trial{1,i},2)-window+1
           for j=1:Nsensor
               n=round(m);
               segm_data.trial{count}(j,:)=speech_irreg_basque_all.trial{1,i}(j,n:n+round(window)-1);
               
           end
           segm_data.time{count}(1,:)=1/fs:1/fs:window/fs;
           count=count+1;
       end
   end
   
   speech_irreg_basque_all.trial=segm_data.trial;
   speech_irreg_basque_all.time=segm_data.time;

N_irreg_basque=size(speech_irreg_basque_all.trial,2);  
load([pathhome,'/subj',num2str(subj(s)),'/',num2str(subj(s)),'_data_clean_rest']) 
w=1:size(data_clean_rest.trial,2);

for perm=1:500
    rng('shuffle');
    p=w(randperm(length(w)));
    
    cfg=[];
    cfg.trials=p(1:N_irreg_basque);
    data_clean_p=ft_selectdata(cfg,data_clean_rest);
    for trl=1:length(w)
        data_clean_p.trial{1,trl} = data_clean_rest.trial{1,p(trl)};
    end
    
    Fdata_irreg_basque_perm=[];Faudio_irreg_basque_perm=[];count=0;
    for i=1:size(speech_irreg_basque_all.trial,2)
        for j=1:Nsens
            Fdata_irreg_basque_perm(j,:,count+1)=fft(data_clean_p.trial{1,i}(j,:),[],2);
        end
        Faudio_irreg_basque_perm(:,count+1)=fft(speech_irreg_basque_all.trial{1,i}(Nsens+1,:),[],2);
        count=count+1;
    end
    
    %Compute coherence for the audio
    
    Faudio_irreg_basque_perm = permute(Faudio_irreg_basque_perm,[3 1 2]);
    Fxx_irreg_basque_perm=[];Fyy_irreg_basque_perm=[];Fxy_irreg_basque_perm=[];
    Fxx_irreg_basque_perm=mean(Fdata_irreg_basque_perm.*conj(Fdata_irreg_basque_perm),3);
    Fyy_irreg_basque_perm=mean(Faudio_irreg_basque_perm.*conj(Faudio_irreg_basque_perm),3);
    Fxy_irreg_basque_perm=mean(Fdata_irreg_basque_perm.*repmat(conj(Faudio_irreg_basque_perm),[Nsens 1]),3);
    coh_speech_irreg_basque_all_perm(:,:,perm)=Fxy_irreg_basque_perm.*conj(Fxy_irreg_basque_perm)./(repmat(Fyy_irreg_basque_perm,[Nsens 1]).*Fxx_irreg_basque_perm);
    
    
end


Nfreq=500;

for ch=1:Nsens-1
    for f=1:Nfreq
        kk=[];
        kk=zscore([coh_real.coh_audio_irreg_basque(ch,f);squeeze(coh_speech_irreg_basque_all_perm(ch,f,:))]);
        z_score_coh_speech_irreg_basque_all(ch,f)=kk(1);
    end
end

%% SPEECH CONDITION - BASQUE MISMATCHING


cfg=[];w=[];p=[];data_clean_p=[];
cfg.trials=[find(speech.trialinfo(:)==92);find(speech.trialinfo(:)==94);find(speech.trialinfo(:)==96);find(speech.trialinfo(:)==98)];
   speech_mm_basque_all=ft_selectdata(cfg,speech);
   
    % Create new trials with overlapping windows

   fs=200;
   window=2.5*fs; %time window
   noverlap=1.25*fs; %Overlap

   Nsensor=size(speech_mm_basque_all.label,1); %number of sensors
   count=1;
   segm_data=[]; m=[];
   for i=1:size(speech_mm_basque_all.trial,2)
       for m=1:window-noverlap:size(speech_mm_basque_all.trial{1,i},2)-window+1
           for j=1:Nsensor
               n=round(m);
               segm_data.trial{count}(j,:)=speech_mm_basque_all.trial{1,i}(j,n:n+round(window)-1);
               
           end
           segm_data.time{count}(1,:)=1/fs:1/fs:window/fs;
           count=count+1;
       end
   end
   
   speech_mm_basque_all.trial=segm_data.trial;
   speech_mm_basque_all.time=segm_data.time;

N_mm_basque=size(speech_mm_basque_all.trial,2);  
load([pathhome,'/subj',num2str(subj(s)),'/',num2str(subj(s)),'_data_clean_rest']) 
w=1:size(data_clean_rest.trial,2);

for perm=1:500
    rng('shuffle');
    p=w(randperm(length(w)));
    
    cfg=[];
    cfg.trials=p(1:N_mm_basque);
    data_clean_p=ft_selectdata(cfg,data_clean_rest);
    
     for trl=1:length(w)
        data_clean_p.trial{1,trl} = data_clean_rest.trial{1,p(trl)};
      end 
    
    Fdata_mm_basque_perm=[];Faudio_mm_basque_perm=[];count=0;
    for i=1:size(speech_mm_basque_all.trial,2)
        for j=1:Nsens
            Fdata_mm_basque_perm(j,:,count+1)=fft(data_clean_p.trial{1,i}(j,:),[],2);
        end
        Faudio_mm_basque_perm(:,count+1)=fft(speech_mm_basque_all.trial{1,i}(Nsens+1,:),[],2);
        count=count+1;
    end
    
    %Compute coherence for the audio
    
    Faudio_mm_basque_perm = permute(Faudio_mm_basque_perm,[3 1 2]);
    Fxx_mm_basque_perm=[];Fyy_mm_basque_perm=[];Fxy_mm_basque_perm=[];
    Fxx_mm_basque_perm=mean(Fdata_mm_basque_perm.*conj(Fdata_mm_basque_perm),3);
    Fyy_mm_basque_perm=mean(Faudio_mm_basque_perm.*conj(Faudio_mm_basque_perm),3);
    Fxy_mm_basque_perm=mean(Fdata_mm_basque_perm.*repmat(conj(Faudio_mm_basque_perm),[Nsens 1]),3);
    coh_speech_mm_basque_all_perm(:,:,perm)=Fxy_mm_basque_perm.*conj(Fxy_mm_basque_perm)./(repmat(Fyy_mm_basque_perm,[Nsens 1]).*Fxx_mm_basque_perm);
    
    
end

Nfreq=500;

for ch=1:Nsens-1
    for f=1:Nfreq
        kk=[];
        kk=zscore([coh_real.coh_audio_mm_basque(ch,f);squeeze(coh_speech_mm_basque_all_perm(ch,f,:))]);
        z_score_coh_speech_mm_basque_all(ch,f)=kk(1);
    end
end

filename=[pathhome, '/Zscore/', num2str(subj(s)),'_coherence_speech_zscore'];
save(filename,'z_score_coh_speech_reg_basque_all','z_score_coh_speech_irreg_basque_all','z_score_coh_speech_mm_basque_all')

%% ZSCORE
coh_real_music=load([pathhome,'/Coherence/subj',num2str(subj(s)),'_coherence_music']);

%% MUSIC CONDITION - BASQUE REGULAR

cfg=[];w=[];p=[];data_clean_p=[];
cfg.trials=[find(speech.trialinfo(:)==71);find(speech.trialinfo(:)==73);find(speech.trialinfo(:)==75);find(speech.trialinfo(:)==77)];
   music_reg_basque_all=ft_selectdata(cfg,speech);
   
% Create new trials with overlapping windows

   fs=200;
   window=2.5*fs; %time window
   noverlap=1.25*fs; %Overlap

    Nsensor=size(music_reg_basque_all.label,1); %number of sensors
   count=1;
   segm_data=[]; m=[];
   for i=1:size(music_reg_basque_all.trial,2)
       for m=1:window-noverlap:size(music_reg_basque_all.trial{1,i},2)-window+1
           for j=1:Nsensor
               n=round(m);
               segm_data.trial{count}(j,:)=music_reg_basque_all.trial{1,i}(j,n:n+round(window)-1);
               
           end
           segm_data.time{count}(1,:)=1/fs:1/fs:window/fs;
           count=count+1;
       end
   end
   
   music_reg_basque_all.trial=segm_data.trial;
   music_reg_basque_all.time=segm_data.time;
   
N_music_reg_basque=size(music_reg_basque_all.trial,2);  
load([pathhome,'/subj',num2str(subj(s)),'/',num2str(subj(s)),'_data_clean_rest']) 
w=1:size(data_clean_rest.trial,2);

for perm=1:500
    rng('shuffle');
    p=w(randperm(length(w)));
    
    cfg=[];
    cfg.trials=p(1:N_music_reg_basque);
    data_clean_p=ft_selectdata(cfg,data_clean_rest);
    
    for trl=1:length(w)
        data_clean_p.trial{1,trl} = data_clean_rest.trial{1,p(trl)};
    end
    
    Fdata_music_reg_basque_perm=[];Faudio_music_reg_basque_perm=[];count=0;
    for i=1:size(music_reg_basque_all.trial,2)
        for j=1:Nsens
            Fdata_music_reg_basque_perm(j,:,count+1)=fft(data_clean_p.trial{1,i}(j,:),[],2);
        end
        Faudio_music_reg_basque_perm(:,count+1)=fft(music_reg_basque_all.trial{1,i}(Nsens+1,:),[],2);
        count=count+1;
    end
    
    %Compute coherence for the audio
    
    Faudio_music_reg_basque_perm = permute(Faudio_music_reg_basque_perm,[3 1 2]);
    Fxx_music_reg_basque_perm=[];Fyy_music_reg_basque_perm=[];Fxy_music_reg_basque_perm=[];
    Fxx_music_reg_basque_perm=mean(Fdata_music_reg_basque_perm.*conj(Fdata_music_reg_basque_perm),3);
    Fyy_music_reg_basque_perm=mean(Faudio_music_reg_basque_perm.*conj(Faudio_music_reg_basque_perm),3);
    Fxy_music_reg_basque_perm=mean(Fdata_music_reg_basque_perm.*repmat(conj(Faudio_music_reg_basque_perm),[Nsens 1]),3);
    coh_music_reg_basque_all_perm(:,:,perm)=Fxy_music_reg_basque_perm.*conj(Fxy_music_reg_basque_perm)./(repmat(Fyy_music_reg_basque_perm,[Nsens 1]).*Fxx_music_reg_basque_perm);    
    
end

Nfreq=500;

for ch=1:Nsens-1
    for f=1:Nfreq
        kk=[];
        kk=zscore([coh_real_music.coh_audio_music_reg_basque(ch,f);squeeze(coh_music_reg_basque_all_perm(ch,f,:))]);
        z_score_coh_music_reg_basque_all(ch,f)=kk(1);
    end
end

%% MUSIC CONDITION - BASQUE IRREGULAR

cfg=[];w=[];p=[];data_clean_p=[];
cfg.trials=[find(speech.trialinfo(:)==81);find(speech.trialinfo(:)==83);find(speech.trialinfo(:)==85);find(speech.trialinfo(:)==87)];
   music_irreg_basque_all=ft_selectdata(cfg,speech);
 
   % Create new trials with overlapping windows

   fs=200;
   window=2.5*fs; %time window
   noverlap=1.25*fs; %Overlap

   
     Nsensor=size(music_irreg_basque_all.label,1); %number of sensors
   count=1;
   segm_data=[]; m=[];
   for i=1:size(music_irreg_basque_all.trial,2)
       for m=1:window-noverlap:size(music_irreg_basque_all.trial{1,i},2)-window+1
           for j=1:Nsensor
               n=round(m);
               segm_data.trial{count}(j,:)=music_irreg_basque_all.trial{1,i}(j,n:n+round(window)-1);
               
           end
           segm_data.time{count}(1,:)=1/fs:1/fs:window/fs;
           count=count+1;
       end
   end
   
   music_irreg_basque_all.trial=segm_data.trial;
   music_irreg_basque_all.time=segm_data.time;
   
N_music_irreg_basque=size(music_irreg_basque_all.trial,2);  
load([pathhome,'/subj',num2str(subj(s)),'/',num2str(subj(s)),'_data_clean_rest']) 
w=1:size(data_clean_rest.trial,2);

for perm=1:500
    rng('shuffle');
    p=w(randperm(length(w)));
    
    cfg=[];
    cfg.trials=p(1:N_music_irreg_basque);
    data_clean_p=ft_selectdata(cfg,data_clean_rest);
    
    for trl=1:length(w)
        data_clean_p.trial{1,trl} = data_clean_rest.trial{1,p(trl)};
    end
    
    Fdata_music_irreg_basque_perm=[];Faudio_music_irreg_basque_perm=[];count=0;
    for i=1:size(music_irreg_basque_all.trial,2)
        for j=1:Nsens
            Fdata_music_irreg_basque_perm(j,:,count+1)=fft(data_clean_p.trial{1,i}(j,:),[],2);
        end
        Faudio_music_irreg_basque_perm(:,count+1)=fft(music_irreg_basque_all.trial{1,i}(Nsens+1,:),[],2);
        count=count+1;
    end
    
    %Compute coherence for the audio
    
    Faudio_music_irreg_basque_perm = permute(Faudio_music_irreg_basque_perm,[3 1 2]);
    Fxx_music_irreg_basque_perm=[];Fyy_music_irreg_basque_perm=[];Fxy_music_irreg_basque_perm=[];
    Fxx_music_irreg_basque_perm=mean(Fdata_music_irreg_basque_perm.*conj(Fdata_music_irreg_basque_perm),3);
    Fyy_music_irreg_basque_perm=mean(Faudio_music_irreg_basque_perm.*conj(Faudio_music_irreg_basque_perm),3);
    Fxy_music_irreg_basque_perm=mean(Fdata_music_irreg_basque_perm.*repmat(conj(Faudio_music_irreg_basque_perm),[Nsens 1]),3);
    coh_music_irreg_basque_all_perm(:,:,perm)=Fxy_music_irreg_basque_perm.*conj(Fxy_music_irreg_basque_perm)./(repmat(Fyy_music_irreg_basque_perm,[Nsens 1]).*Fxx_music_irreg_basque_perm);   
    
end

Nfreq=500;

for ch=1:Nsens-1
    for f=1:Nfreq
        kk=[];
        kk=zscore([coh_real_music.coh_audio_music_irreg_basque(ch,f);squeeze(coh_music_irreg_basque_all_perm(ch,f,:))]);
        z_score_coh_music_irreg_basque_all(ch,f)=kk(1);
    end
end

%% MUSIC CONDITION - BASQUE MISMATCHING

cfg=[];w=[];p=[];data_clean_p=[];
cfg.trials=[find(speech.trialinfo(:)==91);find(speech.trialinfo(:)==93);find(speech.trialinfo(:)==95);find(speech.trialinfo(:)==97)];
   music_mm_basque_all=ft_selectdata(cfg,speech);

 
   % Create new trials with overlapping windows

   fs=200;
   window=2.5*fs; %time window
   noverlap=1.25*fs; %Overlap


   Nsensor=size(music_mm_basque_all.label,1); %number of sensors
   count=1;
   segm_data=[]; m=[];
   for i=1:size(music_mm_basque_all.trial,2)
       for m=1:window-noverlap:size(music_mm_basque_all.trial{1,i},2)-window+1
           for j=1:Nsensor
               n=round(m);
               segm_data.trial{count}(j,:)=music_mm_basque_all.trial{1,i}(j,n:n+round(window)-1);
               
           end
           segm_data.time{count}(1,:)=1/fs:1/fs:window/fs;
           count=count+1;
       end
   end
   
   music_mm_basque_all.trial=segm_data.trial;
   music_mm_basque_all.time=segm_data.time;
   
   
N_music_mm_basque=size(music_mm_basque_all.trial,2);  
load([pathhome,'/subj',num2str(subj(s)),'/',num2str(subj(s)),'_data_clean_rest']) 
w=1:size(data_clean_rest.trial,2);

for perm=1:500
    rng('shuffle');
    p=w(randperm(length(w)));
    
    cfg=[];
    cfg.trials=p(1:N_music_mm_basque);
    data_clean_p=ft_selectdata(cfg,data_clean_rest);
    
    for trl=1:length(w)
        data_clean_p.trial{1,trl} = data_clean_rest.trial{1,p(trl)};
    end 
    
    Fdata_music_mm_basque_perm=[];Faudio_music_mm_basque_perm=[];count=0;
    for i=1:size(music_mm_basque_all.trial,2)
        for j=1:Nsens
            Fdata_music_mm_basque_perm(j,:,count+1)=fft(data_clean_p.trial{1,i}(j,:),[],2);
        end
        Faudio_music_mm_basque_perm(:,count+1)=fft(music_mm_basque_all.trial{1,i}(Nsens+1,:),[],2);
        count=count+1;
    end
    
    %Compute coherence for the audio
    
    Faudio_music_mm_basque_perm = permute(Faudio_music_mm_basque_perm,[3 1 2]);
    Fxx_music_mm_basque_perm=[];Fyy_music_mm_basque_perm=[];Fxy_music_mm_basque_perm=[];
    Fxx_music_mm_basque_perm=mean(Fdata_music_mm_basque_perm.*conj(Fdata_music_mm_basque_perm),3);
    Fyy_music_mm_basque_perm=mean(Faudio_music_mm_basque_perm.*conj(Faudio_music_mm_basque_perm),3);
    Fxy_music_mm_basque_perm=mean(Fdata_music_mm_basque_perm.*repmat(conj(Faudio_music_mm_basque_perm),[Nsens 1]),3);
    coh_music_mm_basque_all_perm(:,:,perm)=Fxy_music_mm_basque_perm.*conj(Fxy_music_mm_basque_perm)./(repmat(Fyy_music_mm_basque_perm,[Nsens 1]).*Fxx_music_mm_basque_perm);
    
end

Nfreq=500;

for ch=1:Nsens-1
    for f=1:Nfreq
        kk=[];
        kk=zscore([coh_real_music.coh_audio_music_mm_basque(ch,f);squeeze(coh_music_mm_basque_all_perm(ch,f,:))]);
        z_score_coh_music_mm_basque_all(ch,f)=kk(1);
    end
end

filename=[pathhome, '/Zscore/', num2str(subj(s)),'_coherence_music_zscore'];
save(filename,'z_score_coh_music_reg_basque_all','z_score_coh_music_irreg_basque_all','z_score_coh_music_mm_basque_all')
