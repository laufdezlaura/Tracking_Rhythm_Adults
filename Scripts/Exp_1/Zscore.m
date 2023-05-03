%% Z-scored coherence for the RhythAdults_1 data
% Author: Laura Fernández-Merino and Mikel Lizarazu
% Date of last change: 3/abril/23
% Matlab version: R2021a on Mac
% Fieldtrip version: 20211001
% This script contains information about computing z-scored coherence from the project RhythAdults.  
% NOTE: All the changes made to the script can be found in my github.* 

clear all;clc;

%% Paths
%% Add paths
restoredefaultpath
pathhome='/Users/laura/Desktop/Midterm_Paper/Analisis/Experiment1/DataRhythAdults_corregida';
%pathfldtrip='/Users/laura/Documents/fieldtrip-20211001';
%addpath(genpath(pathhome));
addpath '/Users/laura/Documents/fieldtrip-20211001'
ft_defaults
%% Define variables
 subj=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 35];

s=1;    %:length(subj)
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
speech=ft_selectdata(cfg,music);

%% ZSCORE
coh_real=load([pathhome,'/Coherence/subj',num2str(subj(s)),'_coherence_speech']);

%% SPEECH CONDITION - SPANISH REGULAR

cfg=[];
cfg.trials=[find(speech.trialinfo(:)==12);find(speech.trialinfo(:)==14);find(speech.trialinfo(:)==16);find(speech.trialinfo(:)==18)];
speech_reg_spanish_all=ft_selectdata(cfg,speech);

N_reg_spanish=size(speech_reg_spanish_all.trial,2);  
load([pathhome,'/subj',num2str(subj(s)),'/',num2str(subj(s)),'_data_clean_rest']) 
w=1:size(data_clean_rest.trial,2);

for perm=1:500
    
    p=w(randperm(length(w)));
    cfg=[];
    cfg.trials=p(1:N_reg_spanish);
    data_clean_p=ft_selectdata(cfg,data_clean_rest);
    
    
    Fdata_reg_spanish_perm=[];Faudio_reg_spanish_perm=[];count=0;
    for i=1:size(speech_reg_spanish_all.trial,2)
        for j=1:Nsens
            Fdata_reg_spanish_perm(j,:,count+1)=fft(data_clean_p.trial{1,i}(j,:),[],2);
        end
        Faudio_reg_spanish_perm(:,count+1)=fft(speech_reg_spanish_all.trial{1,i}(Nsens+1,:),[],2);
        count=count+1;
    end
    
    %Compute coherence for the audio
    
    Faudio_reg_spanish_perm = permute(Faudio_reg_spanish_perm,[3 1 2]);
    Fxx_reg_spanish_perm=[];Fyy_reg_spanish_perm=[];Fxy_reg_spanish_perm=[];
    Fxx_reg_spanish_perm=mean(Fdata_reg_spanish_perm.*conj(Fdata_reg_spanish_perm),3);
    Fyy_reg_spanish_perm=mean(Faudio_reg_spanish_perm.*conj(Faudio_reg_spanish_perm),3);
    Fxy_reg_spanish_perm=mean(Fdata_reg_spanish_perm.*repmat(conj(Faudio_reg_spanish_perm),[Nsens 1]),3);
    coh_speech_reg_spanish_all_perm(:,:,perm)=Fxy_reg_spanish_perm.*conj(Fxy_reg_spanish_perm)./(repmat(Fyy_reg_spanish_perm,[Nsens 1]).*Fxx_reg_spanish_perm);
    

end

Nfreq=500;

for ch=1:Nsens-1
    for f=1:Nfreq
        kk=[];
        kk=zscore([coh_real.coh_audio_reg_spanish(ch,f);squeeze(coh_speech_reg_spanish_all_perm(ch,f,:))]);
        z_score_coh_speech_reg_spanish_all(ch,f)=kk(1);
    end
end

%% SPEECH CONDITION - SPANISH IRREGULAR


cfg=[];w=[];p=[];data_clean_p=[];
cfg.trials=[find(speech.trialinfo(:)==22);find(speech.trialinfo(:)==24);find(speech.trialinfo(:)==26);find(speech.trialinfo(:)==28)];
speech_irreg_spanish_all=ft_selectdata(cfg,speech);

N_irreg_spanish=size(speech_irreg_spanish_all.trial,2);  
load([pathhome,'/subj',num2str(subj(s)),'/',num2str(subj(s)),'_data_clean_rest']) 
w=1:size(data_clean_rest.trial,2);

for perm=1:500
    
    p=w(randperm(length(w)));
    cfg=[];
    cfg.trials=p(1:N_irreg_spanish);
    data_clean_p=ft_selectdata(cfg,data_clean_rest);
    
    
    Fdata_irreg_spanish_perm=[];Faudio_irreg_spanish_perm=[];count=0;
    for i=1:size(speech_irreg_spanish_all.trial,2)
        for j=1:Nsens
            Fdata_irreg_spanish_perm(j,:,count+1)=fft(data_clean_p.trial{1,i}(j,:),[],2);
        end
        Faudio_irreg_spanish_perm(:,count+1)=fft(speech_irreg_spanish_all.trial{1,i}(Nsens+1,:),[],2);
        count=count+1;
    end
    
    %Compute coherence for the audio
    
    Faudio_irreg_spanish_perm = permute(Faudio_irreg_spanish_perm,[3 1 2]);
    Fxx_irreg_spanish_perm=[];Fyy_irreg_spanish_perm=[];Fxy_irreg_spanish_perm=[];
    Fxx_irreg_spanish_perm=mean(Fdata_irreg_spanish_perm.*conj(Fdata_irreg_spanish_perm),3);
    Fyy_irreg_spanish_perm=mean(Faudio_irreg_spanish_perm.*conj(Faudio_irreg_spanish_perm),3);
    Fxy_irreg_spanish_perm=mean(Fdata_irreg_spanish_perm.*repmat(conj(Faudio_irreg_spanish_perm),[Nsens 1]),3);
    coh_speech_irreg_spanish_all_perm(:,:,perm)=Fxy_irreg_spanish_perm.*conj(Fxy_irreg_spanish_perm)./(repmat(Fyy_irreg_spanish_perm,[Nsens 1]).*Fxx_irreg_spanish_perm);
    
    
end

Nfreq=500;

for ch=1:Nsens-1
    for f=1:Nfreq
        kk=[];
        kk=zscore([coh_real.coh_audio_irreg_spanish(ch,f);squeeze(coh_speech_irreg_spanish_all_perm(ch,f,:))]);
        z_score_coh_speech_irreg_spanish_all(ch,f)=kk(1);
    end
end

%% SPEECH CONDITION - SPANISH MISMATCHING


cfg=[];w=[];p=[];data_clean_p=[];
cfg.trials=[find(speech.trialinfo(:)==32);find(speech.trialinfo(:)==34);find(speech.trialinfo(:)==36);find(speech.trialinfo(:)==38)];
speech_mm_spanish_all=ft_selectdata(cfg,speech);

N_mm_spanish=size(speech_mm_spanish_all.trial,2);  
load([pathhome,'/subj',num2str(subj(s)),'/',num2str(subj(s)),'_data_clean_rest'])
w=1:size(data_clean_rest.trial,2);

for perm=1:500
    
    p=w(randperm(length(w)));
    cfg=[];
    cfg.trials=p(1:N_mm_spanish);
    data_clean_p=ft_selectdata(cfg,data_clean_rest);
    
    
    Fdata_mm_spanish_perm=[];Faudio_mm_spanish_perm=[];count=0;
    for i=1:size(speech_mm_spanish_all.trial,2)
        for j=1:Nsens
            Fdata_mm_spanish_perm(j,:,count+1)=fft(data_clean_p.trial{1,i}(j,:),[],2);
        end
        Faudio_mm_spanish_perm(:,count+1)=fft(speech_mm_spanish_all.trial{1,i}(Nsens+1,:),[],2);
        count=count+1;
    end
    
    %Compute coherence for the audio
    
    Faudio_mm_spanish_perm = permute(Faudio_mm_spanish_perm,[3 1 2]);
    Fxx_mm_spanish_perm=[];Fyy_mm_spanish_perm=[];Fxy_mm_spanish_perm=[];
    Fxx_mm_spanish_perm=mean(Fdata_mm_spanish_perm.*conj(Fdata_mm_spanish_perm),3);
    Fyy_mm_spanish_perm=mean(Faudio_mm_spanish_perm.*conj(Faudio_mm_spanish_perm),3);
    Fxy_mm_spanish_perm=mean(Fdata_mm_spanish_perm.*repmat(conj(Faudio_mm_spanish_perm),[Nsens 1]),3);
    coh_speech_mm_spanish_all_perm(:,:,perm)=Fxy_mm_spanish_perm.*conj(Fxy_mm_spanish_perm)./(repmat(Fyy_mm_spanish_perm,[Nsens 1]).*Fxx_mm_spanish_perm);
    
    
end


Nfreq=500;

for ch=1:Nsens-1
    for f=1:Nfreq
        kk=[];
        kk=zscore([coh_real.coh_audio_mm_spanish(ch,f);squeeze(coh_speech_mm_spanish_all_perm(ch,f,:))]);
        z_score_coh_speech_mm_spanish_all(ch,f)=kk(1);
    end
end

%% SPEECH CONDITION - BASQUE REGULAR


cfg=[];w=[];p=[];data_clean_p=[];
cfg.trials=[find(speech.trialinfo(:)==42);find(speech.trialinfo(:)==44);find(speech.trialinfo(:)==46);find(speech.trialinfo(:)==48)];
speech_reg_basque_all=ft_selectdata(cfg,speech);


N_reg_basque=size(speech_reg_basque_all.trial,2); 
load([pathhome,'/subj',num2str(subj(s)),'/',num2str(subj(s)),'_data_clean_rest']) 
w=1:size(data_clean_rest.trial,2);

for perm=1:500
    
    p=w(randperm(length(w)));
    cfg=[];
    cfg.trials=p(1:N_reg_basque);
    data_clean_p=ft_selectdata(cfg,data_clean_rest);
    
    
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
cfg.trials=[find(speech.trialinfo(:)==52);find(speech.trialinfo(:)==54);find(speech.trialinfo(:)==56);find(speech.trialinfo(:)==58)];
speech_irreg_basque_all=ft_selectdata(cfg,speech);


N_irreg_basque=size(speech_irreg_basque_all.trial,2);  
load([pathhome,'/subj',num2str(subj(s)),'/',num2str(subj(s)),'_data_clean_rest']) 
w=1:size(data_clean_rest.trial,2);

for perm=1:500
    
    p=w(randperm(length(w)));
    cfg=[];
    cfg.trials=p(1:N_irreg_basque);
    data_clean_p=ft_selectdata(cfg,data_clean_rest);
    
    
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
cfg.trials=[find(speech.trialinfo(:)==62);find(speech.trialinfo(:)==64);find(speech.trialinfo(:)==66);find(speech.trialinfo(:)==68)];
speech_mm_basque_all=ft_selectdata(cfg,speech);


N_mm_basque=size(speech_mm_basque_all.trial,2); 
load([pathhome,'/subj',num2str(subj(s)),'/',num2str(subj(s)),'_data_clean_rest']) 
w=1:size(data_clean_rest.trial,2);

for perm=1:500
    
    p=w(randperm(length(w)));
    cfg=[];
    cfg.trials=p(1:N_mm_basque);
    data_clean_p=ft_selectdata(cfg,data_clean_rest);
    
    
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
save(filename,'z_score_coh_speech_reg_spanish_all','z_score_coh_speech_irreg_spanish_all','z_score_coh_speech_mm_spanish_all','z_score_coh_speech_reg_basque_all','z_score_coh_speech_irreg_basque_all','z_score_coh_speech_mm_basque_all')

%% ZSCORE
coh_real_music=load([pathhome,'/Coherence/subj',num2str(subj(s)),'_coherence_music']);

%% MUSIC CONDITION - SPANISH REGULAR

cfg=[];w=[];p=[];data_clean_p=[];
cfg.trials=[find(speech.trialinfo(:)==11);find(speech.trialinfo(:)==13);find(speech.trialinfo(:)==15);find(speech.trialinfo(:)==17)];
music_reg_spanish_all=ft_selectdata(cfg,speech);

N_music_reg_spanish=size(music_reg_spanish_all.trial,2);  
load([pathhome,'/subj',num2str(subj(s)),'/',num2str(subj(s)),'_data_clean_rest'])
w=1:size(data_clean_rest.trial,2);

for perm=1:500
    
    p=w(randperm(length(w)));
    cfg=[];
    cfg.trials=p(1:N_music_reg_spanish);
    data_clean_p=ft_selectdata(cfg,data_clean_rest);
    
    Fdata_music_reg_spanish_perm=[];Faudio_music_reg_spanish_perm=[];count=0;
    for i=1:size(music_reg_spanish_all.trial,2)
        for j=1:Nsens
            Fdata_music_reg_spanish_perm(j,:,count+1)=fft(data_clean_p.trial{1,i}(j,:),[],2);
        end
        Faudio_music_reg_spanish_perm(:,count+1)=fft(music_reg_spanish_all.trial{1,i}(Nsens+1,:),[],2);
        count=count+1;
    end
    
    %Compute coherence for the audio
    
    Faudio_music_reg_spanish_perm = permute(Faudio_music_reg_spanish_perm,[3 1 2]);
    Fxx_music_reg_spanish_perm=[];Fyy_music_reg_spanish_perm=[];Fxy_music_reg_spanish_perm=[];
    Fxx_music_reg_spanish_perm=mean(Fdata_music_reg_spanish_perm.*conj(Fdata_music_reg_spanish_perm),3);
    Fyy_music_reg_spanish_perm=mean(Faudio_music_reg_spanish_perm.*conj(Faudio_music_reg_spanish_perm),3);
    Fxy_music_reg_spanish_perm=mean(Fdata_music_reg_spanish_perm.*repmat(conj(Faudio_music_reg_spanish_perm),[Nsens 1]),3);
    coh_music_reg_spanish_all_perm(:,:,perm)=Fxy_music_reg_spanish_perm.*conj(Fxy_music_reg_spanish_perm)./(repmat(Fyy_music_reg_spanish_perm,[Nsens 1]).*Fxx_music_reg_spanish_perm);
    
end

Nfreq=500;

for ch=1:Nsens-1
    for f=1:Nfreq
        kk=[];
        kk=zscore([coh_real_music.coh_audio_music_reg_spanish(ch,f);squeeze(coh_music_reg_spanish_all_perm(ch,f,:))]);
        z_score_coh_music_reg_spanish_all(ch,f)=kk(1);
    end
end

%% MUSIC CONDITION - SPANISH IRREGULAR

cfg=[];w=[];p=[];data_clean_p=[];
cfg.trials=[find(speech.trialinfo(:)==21);find(speech.trialinfo(:)==23);find(speech.trialinfo(:)==25);find(speech.trialinfo(:)==27)];
music_irreg_spanish_all=ft_selectdata(cfg,speech);

N_music_irreg_spanish=size(music_irreg_spanish_all.trial,2);  
load([pathhome,'/subj',num2str(subj(s)),'/',num2str(subj(s)),'_data_clean_rest']) 
w=1:size(data_clean_rest.trial,2);

for perm=1:500
    
    p=w(randperm(length(w)));
    cfg=[];
    cfg.trials=p(1:N_music_irreg_spanish);
    data_clean_p=ft_selectdata(cfg,data_clean_rest);
        
    Fdata_music_irreg_spanish_perm=[];Faudio_music_irreg_spanish_perm=[];count=0;
    for i=1:size(music_irreg_spanish_all.trial,2)
        for j=1:Nsens
            Fdata_music_irreg_spanish_perm(j,:,count+1)=fft(data_clean_p.trial{1,i}(j,:),[],2);
        end
        Faudio_music_irreg_spanish_perm(:,count+1)=fft(music_irreg_spanish_all.trial{1,i}(Nsens+1,:),[],2);
        count=count+1;
    end
    
    %Compute coherence for the audio
    
    Faudio_music_irreg_spanish_perm = permute(Faudio_music_irreg_spanish_perm,[3 1 2]);
    Fxx_music_irreg_spanish_perm=[];Fyy_music_irreg_spanish_perm=[];Fxy_music_irreg_spanish_perm=[];
    Fxx_music_irreg_spanish_perm=mean(Fdata_music_irreg_spanish_perm.*conj(Fdata_music_irreg_spanish_perm),3);
    Fyy_music_irreg_spanish_perm=mean(Faudio_music_irreg_spanish_perm.*conj(Faudio_music_irreg_spanish_perm),3);
    Fxy_music_irreg_spanish_perm=mean(Fdata_music_irreg_spanish_perm.*repmat(conj(Faudio_music_irreg_spanish_perm),[Nsens 1]),3);
    coh_music_irreg_spanish_all_perm(:,:,perm)=Fxy_music_irreg_spanish_perm.*conj(Fxy_music_irreg_spanish_perm)./(repmat(Fyy_music_irreg_spanish_perm,[Nsens 1]).*Fxx_music_irreg_spanish_perm);
    
end

Nfreq=500;

for ch=1:Nsens-1
    for f=1:Nfreq
        kk=[];
        kk=zscore([coh_real_music.coh_audio_music_irreg_spanish(ch,f);squeeze(coh_music_irreg_spanish_all_perm(ch,f,:))]);
        z_score_coh_music_irreg_spanish_all(ch,f)=kk(1);
    end
end

%% MUSIC CONDITION - SPANISH MISMATCHING

cfg=[];w=[];p=[];data_clean_p=[];
cfg.trials=[find(speech.trialinfo(:)==31);find(speech.trialinfo(:)==33);find(speech.trialinfo(:)==35);find(speech.trialinfo(:)==37)];
music_mm_spanish_all=ft_selectdata(cfg,speech);

N_music_mm_spanish=size(music_mm_spanish_all.trial,2);  
load([pathhome,'/subj',num2str(subj(s)),'/',num2str(subj(s)),'_data_clean_rest']) 
w=1:size(data_clean_rest.trial,2);

for perm=1:500
    
    p=w(randperm(length(w)));
    cfg=[];
    cfg.trials=p(1:N_music_mm_spanish);
    data_clean_p=ft_selectdata(cfg,data_clean_rest);   
    
    Fdata_music_mm_spanish_perm=[];Faudio_music_mm_spanish_perm=[];count=0;
    for i=1:size(music_mm_spanish_all.trial,2)
        for j=1:Nsens
            Fdata_music_mm_spanish_perm(j,:,count+1)=fft(data_clean_p.trial{1,i}(j,:),[],2);
        end
        Faudio_music_mm_spanish_perm(:,count+1)=fft(music_mm_spanish_all.trial{1,i}(Nsens+1,:),[],2);
        count=count+1;
    end
    
    %Compute coherence for the audio
    
    Faudio_music_mm_spanish_perm = permute(Faudio_music_mm_spanish_perm,[3 1 2]);
    Fxx_music_mm_spanish_perm=[];Fyy_music_mm_spanish_perm=[];Fxy_music_mm_spanish_perm=[];
    Fxx_music_mm_spanish_perm=mean(Fdata_music_mm_spanish_perm.*conj(Fdata_music_mm_spanish_perm),3);
    Fyy_music_mm_spanish_perm=mean(Faudio_music_mm_spanish_perm.*conj(Faudio_music_mm_spanish_perm),3);
    Fxy_music_mm_spanish_perm=mean(Fdata_music_mm_spanish_perm.*repmat(conj(Faudio_music_mm_spanish_perm),[Nsens 1]),3);
    coh_music_mm_spanish_all_perm(:,:,perm)=Fxy_music_mm_spanish_perm.*conj(Fxy_music_mm_spanish_perm)./(repmat(Fyy_music_mm_spanish_perm,[Nsens 1]).*Fxx_music_mm_spanish_perm);
    
end

Nfreq=500;

for ch=1:Nsens-1
    for f=1:Nfreq
        kk=[];
        kk=zscore([coh_real_music.coh_audio_music_mm_spanish(ch,f);squeeze(coh_music_mm_spanish_all_perm(ch,f,:))]);
        z_score_coh_music_mm_spanish_all(ch,f)=kk(1);
    end
end

%% MUSIC CONDITION - BASQUE REGULAR

cfg=[];w=[];p=[];data_clean_p=[];
cfg.trials=[find(speech.trialinfo(:)==41);find(speech.trialinfo(:)==43);find(speech.trialinfo(:)==45);find(speech.trialinfo(:)==47)];
music_reg_basque_all=ft_selectdata(cfg,speech);

N_music_reg_basque=size(music_reg_basque_all.trial,2);  
load([pathhome,'/subj',num2str(subj(s)),'/',num2str(subj(s)),'_data_clean_rest']) 
w=1:size(data_clean_rest.trial,2);

for perm=1:500
    
    p=w(randperm(length(w)));
    cfg=[];
    cfg.trials=p(1:N_music_reg_basque);
    data_clean_p=ft_selectdata(cfg,data_clean_rest);
    
    
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
cfg.trials=[find(speech.trialinfo(:)==51);find(speech.trialinfo(:)==53);find(speech.trialinfo(:)==55);find(speech.trialinfo(:)==57)];
music_irreg_basque_all=ft_selectdata(cfg,speech);

N_music_irreg_basque=size(music_irreg_basque_all.trial,2); 
load([pathhome,'/subj',num2str(subj(s)),'/',num2str(subj(s)),'_data_clean_rest']) 
w=1:size(data_clean_rest.trial,2);

for perm=1:500
    
    p=w(randperm(length(w)));
    cfg=[];
    cfg.trials=p(1:N_music_irreg_basque);
    data_clean_p=ft_selectdata(cfg,data_clean_rest);
    
    
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
cfg.trials=[find(speech.trialinfo(:)==61);find(speech.trialinfo(:)==63);find(speech.trialinfo(:)==65);find(speech.trialinfo(:)==67)];
music_mm_basque_all=ft_selectdata(cfg,speech);

N_music_mm_basque=size(music_mm_basque_all.trial,2);  
load([pathhome,'/subj',num2str(subj(s)),'/',num2str(subj(s)),'_data_clean_rest'])
w=1:size(data_clean_rest.trial,2);

for perm=1:500
    
    p=w(randperm(length(w)));
    cfg=[];
    cfg.trials=p(1:N_music_mm_basque);
    data_clean_p=ft_selectdata(cfg,data_clean_rest);
    
    
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

filename=[pathhome, '/Zscore/',  num2str(subj(s)),'_coherence_music_zscore'];
save(filename,'z_score_coh_music_reg_spanish_all','z_score_coh_music_irreg_spanish_all','z_score_coh_music_mm_spanish_all','z_score_coh_music_reg_basque_all','z_score_coh_music_irreg_basque_all','z_score_coh_music_mm_basque_all')

