%% Coherence script for the RhythAdults_1 data
% Author: Laura Fernández-Merino and Mikel Lizarazu
% Last change date: 3/04/23
% Matlab version: R2021a on Mac
% Fieldtrip version: 20211001
% This script contains information about computing coherence from the project RhythAdults.  
% NOTE: All the changes made to the script can be found in my github.* 

%% Coherence script for the RhythAdults_1 data
clear all;clc;

%% Add paths
restoredefaultpath
pathhome='/Users/laura/Desktop/Midterm_Paper/Analisis/Experiment1/DataRhythAdults_corregida';
%pathfldtrip='/Users/laura/Documents/fieldtrip-20211001';
%addpath(genpath(pathhome));
addpath '/Users/laura/Documents/fieldtrip-20211001'
ft_defaults
%% Define variables
subj = [1:35];
s=35;    %:length(subj)
clearvars -except pathhome subj s
%% Load data
filename=[pathhome,'/subj',num2str(subj(s)),'/',num2str(subj(s)),'_data_clean'];
load(filename,'data_clean');

%%
% Read all EEG data
cfg=[];
cfg.trials=[find(data_clean.trialinfo(:)==11);find(data_clean.trialinfo(:)==12);find(data_clean.trialinfo(:)==13);...
    find(data_clean.trialinfo(:)==14);find(data_clean.trialinfo(:)==15);find(data_clean.trialinfo(:)==16);...
    find(data_clean.trialinfo(:)==17);find(data_clean.trialinfo(:)==18);...
    find(data_clean.trialinfo(:)==21);find(data_clean.trialinfo(:)==22);find(data_clean.trialinfo(:)==23);...
    find(data_clean.trialinfo(:)==24);find(data_clean.trialinfo(:)==25);find(data_clean.trialinfo(:)==26);...
    find(data_clean.trialinfo(:)==27);find(data_clean.trialinfo(:)==28);...
    find(data_clean.trialinfo(:)==31);find(data_clean.trialinfo(:)==32);find(data_clean.trialinfo(:)==33);...
    find(data_clean.trialinfo(:)==34);find(data_clean.trialinfo(:)==35);find(data_clean.trialinfo(:)==36);...
    find(data_clean.trialinfo(:)==37);find(data_clean.trialinfo(:)==38);...
    find(data_clean.trialinfo(:)==41);find(data_clean.trialinfo(:)==42);find(data_clean.trialinfo(:)==43);...
    find(data_clean.trialinfo(:)==44);find(data_clean.trialinfo(:)==45);find(data_clean.trialinfo(:)==46);...
    find(data_clean.trialinfo(:)==47);find(data_clean.trialinfo(:)==48);...
    find(data_clean.trialinfo(:)==51);find(data_clean.trialinfo(:)==52);find(data_clean.trialinfo(:)==53);...
    find(data_clean.trialinfo(:)==54);find(data_clean.trialinfo(:)==55);find(data_clean.trialinfo(:)==56);...
    find(data_clean.trialinfo(:)==57);find(data_clean.trialinfo(:)==58);...
    find(data_clean.trialinfo(:)==61);find(data_clean.trialinfo(:)==62);find(data_clean.trialinfo(:)==63);...
    find(data_clean.trialinfo(:)==64);find(data_clean.trialinfo(:)==65);find(data_clean.trialinfo(:)==66);...
    find(data_clean.trialinfo(:)==67);find(data_clean.trialinfo(:)==68);];
Sp_Reg=ft_selectdata(cfg,data_clean);

filename = ['/Users/laura/Desktop/Midterm_Paper/Analisis/Experiment1/SentencesOrder_DataRhythAdults_corregida/subj',num2str(subj(s)),'.xlsx'];
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

cfg=[];
cfg.trials=[find(music_clean.trialinfo(:)==11);find(music_clean.trialinfo(:)==13);find(music_clean.trialinfo(:)==15);find(music_clean.trialinfo(:)==17)];
    music_reg_spanish=ft_selectdata(cfg,music_clean);
    
cfg=[];
cfg.trials=[find(music_clean.trialinfo(:)==21);find(music_clean.trialinfo(:)==23);find(music_clean.trialinfo(:)==25);find(music_clean.trialinfo(:)==27)];
   music_irreg_spanish=ft_selectdata(cfg,music_clean);

cfg=[];
cfg.trials=[find(music_clean.trialinfo(:)==31);find(music_clean.trialinfo(:)==33);find(music_clean.trialinfo(:)==35);find(music_clean.trialinfo(:)==37)];
   music_mm_spanish=ft_selectdata(cfg,music_clean);

cfg=[];
cfg.trials=[find(music_clean.trialinfo(:)==41);find(music_clean.trialinfo(:)==43);find(music_clean.trialinfo(:)==45);find(music_clean.trialinfo(:)==47)];
   music_reg_basque=ft_selectdata(cfg,music_clean);

cfg=[];
cfg.trials=[find(music_clean.trialinfo(:)==51);find(music_clean.trialinfo(:)==53);find(music_clean.trialinfo(:)==55);find(music_clean.trialinfo(:)==57)];
   music_irreg_basque=ft_selectdata(cfg,music_clean);

cfg=[];
cfg.trials=[find(music_clean.trialinfo(:)==61);find(music_clean.trialinfo(:)==63);find(music_clean.trialinfo(:)==65);find(music_clean.trialinfo(:)==67)];
   music_mm_basque=ft_selectdata(cfg,music_clean);

%% Obtain the Fourier Transform of the EEG music data and the audio files
%Spanish 
Fdata_music_reg_spanish=[];Faudio_music_reg_spanish=[];count=0;
for i=1:size(music_reg_spanish.trial,2)
    for j=1:Nsens
        Fdata_music_reg_spanish(j,:,count+1)=fft(music_reg_spanish.trial{1,i}(j,:),[],2);
    end
    Faudio_music_reg_spanish(:,count+1)=fft(music_reg_spanish.trial{1,i}(Nsens+1,:),[],2);
    count=count+1;
end

Fdata_music_irreg_spanish=[];Faudio_music_irreg_spanish=[];count=0;
for i=1:size(music_irreg_spanish.trial,2)
    for j=1:Nsens
        Fdata_music_irreg_spanish(j,:,count+1)=fft(music_irreg_spanish.trial{1,i}(j,:),[],2);
    end
    Faudio_music_irreg_spanish(:,count+1)=fft(music_irreg_spanish.trial{1,i}(Nsens+1,:),[],2);
    count=count+1;
end

Fdata_music_mm_spanish=[];Faudio_music_mm_spanish=[];count=0;
for i=1:size(music_mm_spanish.trial,2)
    for j=1:Nsens
        Fdata_music_mm_spanish(j,:,count+1)=fft(music_mm_spanish.trial{1,i}(j,:),[],2);
    end
    Faudio_music_mm_spanish(:,count+1)=fft(music_mm_spanish.trial{1,i}(Nsens+1,:),[],2);
    count=count+1;
end
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

Faudio_music_reg_spanish = permute(Faudio_music_reg_spanish,[3 1 2]);
Fxx_reg_spanish=[];Fyy_reg_spanish=[];Fxy_reg_spanish=[];coh_audio_music_reg_spanish=[];
Fxx_reg_spanish=mean(Fdata_music_reg_spanish.*conj(Fdata_music_reg_spanish),3);
Fyy_reg_spanish=mean(Faudio_music_reg_spanish.*conj(Faudio_music_reg_spanish),3);
Fxy_reg_spanish=mean(Fdata_music_reg_spanish.*repmat(conj(Faudio_music_reg_spanish),[Nsens 1]),3); 
coh_audio_music_reg_spanish(:,:)=Fxy_reg_spanish.*conj(Fxy_reg_spanish)./(repmat(Fyy_reg_spanish,[Nsens 1]).*Fxx_reg_spanish);

Faudio_music_irreg_spanish = permute(Faudio_music_irreg_spanish,[3 1 2]);
Fxx_irreg_spanish=[];Fyy_irreg_spanish=[];Fxy_irreg_spanish=[];coh_audio_music_irreg_spanish=[];
Fxx_irreg_spanish=mean(Fdata_music_irreg_spanish.*conj(Fdata_music_irreg_spanish),3);
Fyy_irreg_spanish=mean(Faudio_music_irreg_spanish.*conj(Faudio_music_irreg_spanish),3);
Fxy_irreg_spanish=mean(Fdata_music_irreg_spanish.*repmat(conj(Faudio_music_irreg_spanish),[Nsens 1]),3); 
coh_audio_music_irreg_spanish(:,:)=Fxy_irreg_spanish.*conj(Fxy_irreg_spanish)./(repmat(Fyy_irreg_spanish,[Nsens 1]).*Fxx_irreg_spanish);

Faudio_music_mm_spanish = permute(Faudio_music_mm_spanish,[3 1 2]);
Fxx_mm_spanish=[];Fyy_mm_spanish=[];Fxy_mm_spanish=[];coh_audio_music_mm_spanish=[];
Fxx_mm_spanish=mean(Fdata_music_mm_spanish.*conj(Fdata_music_mm_spanish),3);
Fyy_mm_spanish=mean(Faudio_music_mm_spanish.*conj(Faudio_music_mm_spanish),3);
Fxy_mm_spanish=mean(Fdata_music_mm_spanish.*repmat(conj(Faudio_music_mm_spanish),[Nsens 1]),3); 
coh_audio_music_mm_spanish(:,:)=Fxy_mm_spanish.*conj(Fxy_mm_spanish)./(repmat(Fyy_mm_spanish,[Nsens 1]).*Fxx_mm_spanish);

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
save(filename,'coh_audio_music_reg_spanish','coh_audio_music_irreg_spanish','coh_audio_music_mm_spanish','coh_audio_music_reg_basque','coh_audio_music_irreg_basque','coh_audio_music_mm_basque')

%% Create the speech conditions separated by lang and rhythm condition for ALL TRIALS 


cfg=[];
cfg.trials=[find(music_clean.trialinfo(:)==12);find(music_clean.trialinfo(:)==14);find(music_clean.trialinfo(:)==16);find(music_clean.trialinfo(:)==18)];
    speech_reg_spanish=ft_selectdata(cfg,music_clean);
    
cfg=[];
cfg.trials=[find(music_clean.trialinfo(:)==22);find(music_clean.trialinfo(:)==24);find(music_clean.trialinfo(:)==26);find(music_clean.trialinfo(:)==28)];
   speech_irreg_spanish=ft_selectdata(cfg,music_clean);

cfg=[];
cfg.trials=[find(music_clean.trialinfo(:)==32);find(music_clean.trialinfo(:)==34);find(music_clean.trialinfo(:)==36);find(music_clean.trialinfo(:)==38)];
   speech_mm_spanish=ft_selectdata(cfg,music_clean);

cfg=[];
cfg.trials=[find(music_clean.trialinfo(:)==42);find(music_clean.trialinfo(:)==44);find(music_clean.trialinfo(:)==46);find(music_clean.trialinfo(:)==48)];
   speech_reg_basque=ft_selectdata(cfg,music_clean);

cfg=[];
cfg.trials=[find(music_clean.trialinfo(:)==52);find(music_clean.trialinfo(:)==54);find(music_clean.trialinfo(:)==56);find(music_clean.trialinfo(:)==58)];
   speech_irreg_basque=ft_selectdata(cfg,music_clean);

cfg=[];
cfg.trials=[find(music_clean.trialinfo(:)==62);find(music_clean.trialinfo(:)==64);find(music_clean.trialinfo(:)==66);find(music_clean.trialinfo(:)==68)];
   speech_mm_basque=ft_selectdata(cfg,music_clean);

%% Obtain the Fourier Transform of the EEG speech data and the audio files
%Spanish 
Fdata_reg_spanish=[];Faudio_reg_spanish=[];count=0;
for i=1:size(speech_reg_spanish.trial,2)
    for j=1:Nsens
        Fdata_reg_spanish(j,:,count+1)=fft(speech_reg_spanish.trial{1,i}(j,:),[],2);
    end
    Faudio_reg_spanish(:,count+1)=fft(speech_reg_spanish.trial{1,i}(Nsens+1,:),[],2);
    count=count+1;
end

Fdata_irreg_spanish=[];Faudio_irreg_spanish=[];count=0;
for i=1:size(speech_irreg_spanish.trial,2)
    for j=1:Nsens
        Fdata_irreg_spanish(j,:,count+1)=fft(speech_irreg_spanish.trial{1,i}(j,:),[],2);
    end
    Faudio_irreg_spanish(:,count+1)=fft(speech_irreg_spanish.trial{1,i}(Nsens+1,:),[],2);
    count=count+1;
end

Fdata_mm_spanish=[];Faudio_mm_spanish=[];count=0;
for i=1:size(speech_mm_spanish.trial,2)
    for j=1:Nsens
        Fdata_mm_spanish(j,:,count+1)=fft(speech_mm_spanish.trial{1,i}(j,:),[],2);
    end
    Faudio_mm_spanish(:,count+1)=fft(speech_mm_spanish.trial{1,i}(Nsens+1,:),[],2);
    count=count+1;
end
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

Faudio_reg_spanish = permute(Faudio_reg_spanish,[3 1 2]);
Fxx_reg_spanish=[];Fyy_reg_spanish=[];Fxy_reg_spanish=[];coh_audio_reg_spanish=[];
Fxx_reg_spanish=mean(Fdata_reg_spanish.*conj(Fdata_reg_spanish),3);
Fyy_reg_spanish=mean(Faudio_reg_spanish.*conj(Faudio_reg_spanish),3);
Fxy_reg_spanish=mean(Fdata_reg_spanish.*repmat(conj(Faudio_reg_spanish),[Nsens 1]),3); 
coh_audio_reg_spanish(:,:)=Fxy_reg_spanish.*conj(Fxy_reg_spanish)./(repmat(Fyy_reg_spanish,[Nsens 1]).*Fxx_reg_spanish);

Faudio_irreg_spanish = permute(Faudio_irreg_spanish,[3 1 2]);
Fxx_irreg_spanish=[];Fyy_irreg_spanish=[];Fxy_irreg_spanish=[];coh_audio_irreg_spanish=[];
Fxx_irreg_spanish=mean(Fdata_irreg_spanish.*conj(Fdata_irreg_spanish),3);
Fyy_irreg_spanish=mean(Faudio_irreg_spanish.*conj(Faudio_irreg_spanish),3);
Fxy_irreg_spanish=mean(Fdata_irreg_spanish.*repmat(conj(Faudio_irreg_spanish),[Nsens 1]),3); 
coh_audio_irreg_spanish(:,:)=Fxy_irreg_spanish.*conj(Fxy_irreg_spanish)./(repmat(Fyy_irreg_spanish,[Nsens 1]).*Fxx_irreg_spanish);

Faudio_mm_spanish = permute(Faudio_mm_spanish,[3 1 2]);
Fxx_mm_spanish=[];Fyy_mm_spanish=[];Fxy_mm_spanish=[];coh_audio_mm_spanish=[];
Fxx_mm_spanish=mean(Fdata_mm_spanish.*conj(Fdata_mm_spanish),3);
Fyy_mm_spanish=mean(Faudio_mm_spanish.*conj(Faudio_mm_spanish),3);
Fxy_mm_spanish=mean(Fdata_mm_spanish.*repmat(conj(Faudio_mm_spanish),[Nsens 1]),3); 
coh_audio_mm_spanish(:,:)=Fxy_mm_spanish.*conj(Fxy_mm_spanish)./(repmat(Fyy_mm_spanish,[Nsens 1]).*Fxx_mm_spanish);

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
save(filename,'coh_audio_reg_spanish','coh_audio_irreg_spanish','coh_audio_mm_spanish','coh_audio_reg_basque','coh_audio_irreg_basque','coh_audio_mm_basque')

