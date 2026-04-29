%% Beyesian factor for Exp 3
codepath = '';
addpath(genpath( [codepath,'\subfunctions\']))


datapath = [codepath,'\datafiles\'];
datafigspath =  [codepath,'\data4figs\'];
figsavepath = [codepath,'\figures\'];
subinfopath = [codepath,'\datafiles\Subjects_Exp3\ExpInfo\'];
sacsavepath = [codepath,'\datafiles\Subjects_Exp3\saccades\'];
eegsavepath = [codepath,'\datafiles\Subjects_Exp3\EEG\'];
TFRsavepath = [codepath,'\datafiles\Subjects_Exp3\TFR\'];

%% get TFR for Exp 3
% behavioral results
savename = [datapath,'Beh_Exp3.mat'];
load(savename)

rng(2025);
for s = 1:length(Mem)
    disp(s)
    savename = [TFRsavepath,sprintf('TFR_PVstim_allTrials_sub%d_%s.mat',s,'Exp3')];
    load(savename,'tfralltrl')
    
    % load trial info
    load([sacsavepath,sprintf('TrialEyeInfo_%d_%s.mat',s,'Exp3')])
    %% select trials
    % select only the intersect between clean EEG and clean eye trials
    cfg = [];
    cfg.trials = ismember(tfralltrl.trialinfo(:,1),Eye_trlinfo(Eye_trlinfo(:,8)==1,1));
    tfrfinal = ft_selectdata(cfg,tfralltrl);
    tfrfinal.trialinfo(:,4) = Eye_trlinfo(dsearchn(Eye_trlinfo(:,1),tfrfinal.trialinfo(:,1)),5);
    tfrfinal.trialinfo(:,5) = Eye_trlinfo(dsearchn(Eye_trlinfo(:,1),tfrfinal.trialinfo(:,1)),7);
    %% Sep by Memmory
    
    cfg= [];
    cfg.trials = find(tfrfinal.trialinfo(:,2)==1 );
    A  = ft_selectdata(cfg, tfrfinal);
    cfg.trials = find(tfrfinal.trialinfo(:,2)==0 );
    B  = ft_selectdata(cfg, tfrfinal);
    
    [tfr_C1,tfr_C2,n_trl] = TFR_balance_trial(A,B,100);
    n_trl_C1C2(s,1:2) = n_trl;
        
    % baseline correct on avg data
    % all trials
    tfr_all_raw{s} = ft_freqdescriptives([], tfrfinal);
    % 
    tfr_C1_raw{s} = tfr_C1;
    tfr_C2_raw{s} = tfr_C2;
    
    
   %% by median explore index
    cri = nanmedian(tfrfinal.trialinfo(:,5));
    
    cfg= [];
    cfg.trials = tfrfinal.trialinfo(:,5)<=cri;
    n_trl_ExpIdx(s,1) = sum(cfg.trials);
    tfr_MoreExpl_raw{s}  = ft_freqdescriptives(cfg, tfrfinal);
    cfg.trials = tfrfinal.trialinfo(:,5)>=cri;
    n_trl_ExpIdx(s,2) = sum(cfg.trials);
    tfr_LessExpl_raw{s}  = ft_freqdescriptives(cfg,tfrfinal);
    clearvars cri
    

end

% baseline correction 
baseline_t = [-1,-0.5];

for s = 1:length(tfr_C1_raw)
cfg = [];
    cfg.baseline =baseline_t;
    cfg.baselinetype = 'db';
    
    tfr_C1_ave{s} = ft_freqbaseline(cfg,tfr_C1_raw{s});
    tfr_C2_ave{s} = ft_freqbaseline(cfg,tfr_C2_raw{s});
    
    tfr_Mexpl_ave{s} = ft_freqbaseline(cfg,tfr_MoreExpl_raw{s});
    tfr_Lexpl_ave{s} = ft_freqbaseline(cfg,tfr_LessExpl_raw{s});
   
 
end

%% Get average SME

% Load stat from Exp 1.
load([datafigspath,'Fig1C_left.mat'])
%load([datafigspath,'Fig4D_left.mat'])
% Get averaged ab SME for Experiment 3
% define a toi & foi
toi = [1,3.5]; % time
foi = [10,20]; % frequency
coi = stat_mem.label(sum(sum(stat_mem.mask,3),2)~=0); % sig channels

%%
for s = 1:length(tfr_C1_raw)
cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'powspctrm';
SME = ft_math(cfg,tfr_C1_ave{s},tfr_C2_ave{s});
cfg = [];
cfg.latency = toi;
cfg.frequency = foi;
cfg.channel = coi;
abSME = ft_selectdata(cfg,SME);
avg_abSME(s,1) = mean(abSME.powspctrm(:));
end

[bf10,p] = bf.ttest(avg_abSME);
