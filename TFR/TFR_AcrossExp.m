%% TFR power contrasts across experiments using a common time-frequency definition
codepath = '';
addpath(genpath( [codepath,'\subfunctions\']))


datapath = [codepath,'\datafiles\'];
datafigspath =  [codepath,'\data4figs\'];
figsavepath = [codepath,'\figures\'];

PreDef_time = [1,3.5];
PreDef_Freq = [10,20];
PreDef_Chan = {'P*','O*'};

%% Experiment 1
subinfopath = [codepath,'\datafiles\Subjects_Exp1\ExpInfo\'];
sacsavepath = [codepath,'\datafiles\Subjects_Exp1\saccades\'];
eegsavepath = [codepath,'\datafiles\Subjects_Exp1\EEG\'];
TFRsavepath = [codepath,'\datafiles\Subjects_Exp1\TFR\'];
% behavioral results
savename = [datapath,'Beh_Exp1.mat'];
load(savename)


rng(2024);
for s = 1:20
    disp(s)
    savename = [TFRsavepath,sprintf('TFR_PVstim_allTrials_sub%d_%s.mat',s,'Exp1')];
    load(savename,'tfralltrl')
    
    % load trial info
    load([sacsavepath,sprintf('TrialEyeInfo_%d_%s.mat',s,'Exp1')])
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
    
    %% by median N saccades
    cri = nanmedian(tfrfinal.trialinfo(:,4));
    
    cfg= [];
    cfg.trials = tfrfinal.trialinfo(:,4)>=cri;
    n_trl_Nsac(s,1) = sum(cfg.trials);
    tfr_MoreNsac_raw{s}  = ft_freqdescriptives(cfg, tfrfinal);
    cfg.trials = tfrfinal.trialinfo(:,4)<=cri;
    n_trl_Nsac(s,2) = sum(cfg.trials);
    tfr_LessNsac_raw{s}  = ft_freqdescriptives(cfg,tfrfinal);
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
    
    tfr_Msac_ave{s} = ft_freqbaseline(cfg,tfr_MoreNsac_raw{s});
    tfr_Lsac_ave{s} = ft_freqbaseline(cfg,tfr_LessNsac_raw{s});
   
 
end

% Average for main effect
cfg = [];
cfg.latency= PreDef_time;
cfg.avgovertime = 'yes';
cfg.channel = PreDef_Chan;
cfg.avgoverchan = 'yes';
cfg.frequency = PreDef_Freq;
cfg.avgoverfreq = 'yes';

for s = 1:length(tfr_C1_raw)
    
    disp(s);
    C1 = ft_selectdata(cfg,tfr_C1_ave{s});
    C2 = ft_selectdata(cfg,tfr_C2_ave{s});
    
    M1 = ft_selectdata(cfg,tfr_Msac_ave{s});
    M2 = ft_selectdata(cfg,tfr_Lsac_ave{s});
    
    AvgTFR_RF_Exp1(s,1)=C1.powspctrm(:);
    AvgTFR_RF_Exp1(s,2)=C2.powspctrm(:);
    
    AvgTFR_ML_Exp1(s,1)=M1.powspctrm(:);
    AvgTFR_ML_Exp1(s,2)=M2.powspctrm(:);
    clearvars C1 C2 M1 M2
    
end

% save the result
save([datafigspath,'SuppleFig12A.mat'],'AvgTFR_RF_Exp1','AvgTFR_ML_Exp1');

%% Experiment 2
clearvars -except  codepath datapath datafigspath figsavepath  PreDef_*
subinfopath = [codepath,'\datafiles\Subjects_Exp2\ExpInfo\'];
sacsavepath = [codepath,'\datafiles\Subjects_Exp2\saccades\'];
eegsavepath = [codepath,'\datafiles\Subjects_Exp2\EEG\'];
TFRsavepath = [codepath,'\datafiles\Subjects_Exp2\TFR\'];
% behavioral results
savename = [datapath,'Beh_Exp2.mat'];
load(savename)

% Separate conditions
load([TFRsavepath,'TFR_PVstim_allTrials_allsub_Exp2.mat'])
% treat each channel as a subject
% baseline correction
chan_count = 0;
bl =  [-1,-0.5];
for s = 1:4
    %
    TFR_data = tfralltrl{s};
    spike_data = spike_mark_sub{s};
    rej_data = trial_rej_sub{s} ;
    
    TFR_data = MarkTFRNans_ft(TFR_data,spike_data); % replace the marked artifacts as NaNs
    chanFSlabel = TFR_data.elecinfo;
    TFR_data = rmfield(TFR_data,'elecinfo');

    % load trial info
    load([sacsavepath,sprintf('TrialEyeInfo_%d_%s.mat',s,'Exp2')])
    for ichan = 1:length(TFR_data.label)
        chan_count = chan_count +1;
        % treat channel as a unit
        cfg              = [];
        cfg.baseline     = bl;
        cfg.baselinetype = 'db';
        cfgtrl           = [];
        cfgtrl.channel      = TFR_data.label(ichan);
        cfgtrl.trials       = find(~ismember(TFR_data.trialinfo(:,1),rej_data{ichan}) & ismember(TFR_data.trialinfo(:,1),Eye_trlinfo(:,1)));
        TFR_enc_bc{chan_count} = ft_freqbaseline(cfg, ft_freqdescriptives(cfgtrl,TFR_data));
        TFR_enc_bc{chan_count}.orig_label = TFR_enc_bc{chan_count}.label;
        TFR_enc_bc{chan_count}.label      = {'Posterior'};
        TFR_enc_bc{chan_count}.FSwm       = chanFSlabel(ichan,:);
        
        cfg              = [];
        cfg.baseline     = bl;
        cfg.baselinetype = 'db';
        cfgtrl           = [];
        cfgtrl.channel      = TFR_data.label(ichan);
        cfgtrl.trials       = find(~ismember(TFR_data.trialinfo(:,1),rej_data{ichan})&TFR_data.trialinfo(:,2)==1&ismember(TFR_data.trialinfo(:,1),Eye_trlinfo(:,1)));
        TFR_enc_bc_rem{chan_count} = ft_freqbaseline(cfg, ft_freqdescriptives(cfgtrl,TFR_data));
        TFR_enc_bc_rem{chan_count}.orig_label = TFR_enc_bc_rem{chan_count}.label;
        TFR_enc_bc_rem{chan_count}.label      = {'Posterior'};
        TFR_enc_bc_rem{chan_count}.FSwm       = chanFSlabel(ichan,:);
        
        cfg              = [];
        cfg.baseline     = bl;
        cfg.baselinetype = 'db';
        cfgtrl           = [];
        cfgtrl.channel      = TFR_data.label(ichan);
        cfgtrl.trials       = find(~ismember(TFR_data.trialinfo(:,1),rej_data{ichan})&TFR_data.trialinfo(:,2)==0&ismember(TFR_data.trialinfo(:,1),Eye_trlinfo(:,1)));
        TFR_enc_bc_for{chan_count} = ft_freqbaseline(cfg, ft_freqdescriptives(cfgtrl,TFR_data));
        TFR_enc_bc_for{chan_count}.orig_label = TFR_enc_bc_for{chan_count}.label;
        TFR_enc_bc_for{chan_count}.label      = {'Posterior'};
        TFR_enc_bc_for{chan_count}.FSwm       = chanFSlabel(ichan,:);
        
        FS_label{chan_count,1} = s;
        FS_label{chan_count,2} = TFR_enc_bc{chan_count}.orig_label{:};
        
        
        %% sep by saccades
        sac_cri =  nanmedian(Eye_trlinfo(~ismember(Eye_trlinfo(:,1),rej_data{ichan}),5));
        sac_info = nan(length(TFR_data.trialinfo(:,1)),1);
        sac_info(Eye_trlinfo(:,1)) = Eye_trlinfo(:,5);

        
        cfg              = [];
        cfg.baseline     = bl;
        cfg.baselinetype = 'db';
        cfgtrl           = [];
        cfgtrl.channel      = TFR_data.label(ichan);
        cfgtrl.trials       = find(~ismember(TFR_data.trialinfo(:,1),rej_data{ichan})&sac_info>=sac_cri);
        TFR_enc_bc_moresac{chan_count} = ft_freqbaseline(cfg, ft_freqdescriptives(cfgtrl,TFR_data));
        TFR_enc_bc_moresac{chan_count}.orig_label = TFR_enc_bc_moresac{chan_count}.label;
        TFR_enc_bc_moresac{chan_count}.label      = {'Posterior'};
        TFR_enc_bc_moresac{chan_count}.FSwm       = chanFSlabel(ichan,:);
        
        
        cfg              = [];
        cfg.baseline     = bl;
        cfg.baselinetype = 'db';
        cfgtrl           = [];
        cfgtrl.channel      = TFR_data.label(ichan);
        cfgtrl.trials       = find(~ismember(TFR_data.trialinfo(:,1),rej_data{ichan})&sac_info<=sac_cri);
        TFR_enc_bc_lesssac{chan_count} = ft_freqbaseline(cfg, ft_freqdescriptives(cfgtrl,TFR_data));
        TFR_enc_bc_lesssac{chan_count}.orig_label = TFR_enc_bc_lesssac{chan_count}.label;
        TFR_enc_bc_lesssac{chan_count}.label      = {'Posterior'};
        TFR_enc_bc_lesssac{chan_count}.FSwm       = chanFSlabel(ichan,:);

        
    end
end

% Select
cfg = [];
cfg.latency= PreDef_time;
cfg.avgovertime = 'yes';
cfg.frequency = PreDef_Freq;
cfg.avgoverfreq = 'yes';

for s = 1:length(TFR_enc_bc_rem)
    
    disp(s);
    C1 = ft_selectdata(cfg,TFR_enc_bc_rem{s});
    C2 = ft_selectdata(cfg,TFR_enc_bc_for{s});
    
    M1 = ft_selectdata(cfg,TFR_enc_bc_moresac{s});
    M2 = ft_selectdata(cfg,TFR_enc_bc_lesssac{s});
    
    AvgTFR_RF_Exp2(s,1)=C1.powspctrm(:);
    AvgTFR_RF_Exp2(s,2)=C2.powspctrm(:);
    
    AvgTFR_ML_Exp2(s,1)=M1.powspctrm(:);
    AvgTFR_ML_Exp2(s,2)=M2.powspctrm(:);
    clearvars C1 C2 M1 M2
    
end

save([datafigspath,'SuppleFig12B.mat'],'AvgTFR_RF_Exp2','AvgTFR_ML_Exp2');
%% Experiment 3
clearvars -except  codepath datapath datafigspath figsavepath  PreDef_*
subinfopath = [codepath,'\datafiles\Subjects_Exp3\ExpInfo\'];
sacsavepath = [codepath,'\datafiles\Subjects_Exp3\saccades\'];
eegsavepath = [codepath,'\datafiles\Subjects_Exp3\EEG\'];
TFRsavepath = [codepath,'\datafiles\Subjects_Exp3\TFR\'];
% behavioral results
savename = [datapath,'Beh_Exp3.mat'];
load(savename)

% Analysis
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



% Average for main effect
cfg = [];
cfg.latency= PreDef_time;
cfg.avgovertime = 'yes';
cfg.channel = PreDef_Chan;
cfg.avgoverchan = 'yes';
cfg.frequency = PreDef_Freq;
cfg.avgoverfreq = 'yes';

for s = 1:length(tfr_C1_raw)
    
    disp(s);
    C1 = ft_selectdata(cfg,tfr_C1_ave{s});
    C2 = ft_selectdata(cfg,tfr_C2_ave{s});
    
    M1 = ft_selectdata(cfg,tfr_Mexpl_ave{s});
    M2 = ft_selectdata(cfg,tfr_Lexpl_ave{s});
    
    AvgTFR_RF_Exp3(s,1)=C1.powspctrm(:);
    AvgTFR_RF_Exp3(s,2)=C2.powspctrm(:);
    
    AvgTFR_ML_Exp3(s,1)=M1.powspctrm(:);
    AvgTFR_ML_Exp3(s,2)=M2.powspctrm(:);
    clearvars C1 C2 M1 M2
    
end

% save the result
save([datafigspath,'SuppleFig12C.mat'],'AvgTFR_RF_Exp3','AvgTFR_ML_Exp3');

%% Experiment 4
clearvars -except  codepath datapath datafigspath figsavepath  PreDef_*
subinfopath = [codepath,'\datafiles\Subjects_Exp4\ExpInfo\'];
sacsavepath = [codepath,'\datafiles\Subjects_Exp4\saccades\'];
eegsavepath = [codepath,'\datafiles\Subjects_Exp4\EEG\'];
TFRsavepath = [codepath,'\datafiles\Subjects_Exp4\TFR\'];
% behavioral results
savename = [datapath,'Beh_Exp4.mat'];
load(savename)

% Analysis
rng(2025);
for s = 1:34
    disp(s)
    savename = [TFRsavepath,sprintf('TFR_PVstim_allTrials_sub%d_%s.mat',s,'Exp4')];
    load(savename,'tfralltrl')
    
    % load trial info
    load([sacsavepath,sprintf('TrialEyeInfo_%d_%s.mat',s,'Exp4')])
    %% select trials
    % select only the intersect between clean EEG and clean eye trials
    cfg = [];
    cfg.trials = ismember(tfralltrl.trialinfo(:,1),Eye_trlinfo(Eye_trlinfo(:,9)==1,1));
    tfrfinal = ft_selectdata(cfg,tfralltrl);
    tfrfinal.trialinfo(:,2:3) = Eye_trlinfo(dsearchn(Eye_trlinfo(:,1),tfrfinal.trialinfo(:,1)),[2,4]);% Memory & View condition
    tfrfinal.trialinfo(:,4)   = Eye_trlinfo(dsearchn(Eye_trlinfo(:,1),tfrfinal.trialinfo(:,1)),6); % N PV saccade
    tfrfinal.trialinfo(:,5)   = Eye_trlinfo(dsearchn(Eye_trlinfo(:,1),tfrfinal.trialinfo(:,1)),8); % Explore Index
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
    
 
   
    
    %% by Exploration Index across each Viewing condition 
    cri1 = nanmedian(tfrfinal.trialinfo(tfrfinal.trialinfo(:,3)==1,5));
    cri2 = nanmedian(tfrfinal.trialinfo(tfrfinal.trialinfo(:,3)==2,5));
    cri3 = nanmedian(tfrfinal.trialinfo(tfrfinal.trialinfo(:,3)==3,5));
    
    cfg= [];
    cfg.trials = (tfrfinal.trialinfo(:,5)>=cri1 & tfrfinal.trialinfo(:,3)==1) | (tfrfinal.trialinfo(:,5)>=cri2 & tfrfinal.trialinfo(:,3)==2) | (tfrfinal.trialinfo(:,5)>=cri3 & tfrfinal.trialinfo(:,3)==3);
    n_trl_FixExp_EI(s,1) = sum(cfg.trials);
    tfr_Fix_raw_EI{s}  = ft_freqdescriptives(cfg, tfrfinal);
    
    cfg.trials = (tfrfinal.trialinfo(:,5)<=cri1 & tfrfinal.trialinfo(:,3)==1) | (tfrfinal.trialinfo(:,5)<=cri2 & tfrfinal.trialinfo(:,3)==2) | (tfrfinal.trialinfo(:,5)<=cri3 & tfrfinal.trialinfo(:,3)==3);
    n_trl_FixExp_EI(s,2) = sum(cfg.trials);
    tfr_Explr_raw_EI{s}  = ft_freqdescriptives(cfg,tfrfinal);
    
    
    
    clearvars cri*

end

% baseline correction 
baseline_t = [-1,-0.5];

for s = 1:length(tfr_C1_raw)
cfg = [];
    cfg.baseline =baseline_t;
    cfg.baselinetype = 'db';
    
    % Rem Vs Forg
    tfr_C1_ave{s} = ft_freqbaseline(cfg,tfr_C1_raw{s});
    tfr_C2_ave{s} = ft_freqbaseline(cfg,tfr_C2_raw{s});
    
    % exploration index
    tfr_Fixation_ave_EI{s} = ft_freqbaseline(cfg,tfr_Fix_raw_EI{s});
    tfr_Exploration_ave_EI{s} = ft_freqbaseline(cfg,tfr_Explr_raw_EI{s});
    
    
 
end

% Average for main effect
cfg = [];
cfg.latency= PreDef_time;
cfg.avgovertime = 'yes';
cfg.channel = PreDef_Chan;
cfg.avgoverchan = 'yes';
cfg.frequency = PreDef_Freq;
cfg.avgoverfreq = 'yes';

for s = 1:length(tfr_C1_raw)
    
    disp(s);
    C1 = ft_selectdata(cfg,tfr_C1_ave{s});
    C2 = ft_selectdata(cfg,tfr_C2_ave{s});
    
    M1 = ft_selectdata(cfg,tfr_Exploration_ave_EI{s});
    M2 = ft_selectdata(cfg,tfr_Fixation_ave_EI{s});
    
    AvgTFR_RF_Exp4(s,1)=C1.powspctrm(:);
    AvgTFR_RF_Exp4(s,2)=C2.powspctrm(:);
    
    AvgTFR_ML_Exp4(s,1)=M1.powspctrm(:);
    AvgTFR_ML_Exp4(s,2)=M2.powspctrm(:);
    clearvars C1 C2 M1 M2
    
end

% save the result
save([datafigspath,'SuppleFig12D.mat'],'AvgTFR_RF_Exp4','AvgTFR_ML_Exp4');

%% Figure

plot_save=1;

load([datafigspath,'SuppleFig12A.mat']);
load([datafigspath,'SuppleFig12B.mat']);
load([datafigspath,'SuppleFig12C.mat']);
load([datafigspath,'SuppleFig12D.mat']);

%%
% Experiment 1
sizefig = [0.1,0.1,0.2,0.6];
cb = repmat([0.5,0.5,0.5],2,1);
fig_savename = 'SuppleFig12A_left_Exp1';
for fig=21
    %% 
    f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig);clf;
    h=daboxplot([AvgTFR_RF_Exp1(:,[1:2])],[repmat(1,size(AvgTFR_RF_Exp1,1),1)],...
        'scatter' ,2,...
        'color',cb,...
        'boxalpha', 0.5,...
        'jitter',1,...
        'linkline',0,...
        'withinlines',1,...
        'xtlabels',{'Rem','Forg'})
    
    ylabel('Power to baseline [dB]')
    xlabel('Memory')
   
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
    
   
    
    if plot_save
    saveas(gcf,[figsavepath,fig_savename,'.emf'])
    saveas(gcf,[figsavepath,fig_savename,'.png'])
    end
end 

sizefig = [0.1,0.1,0.2,0.6];
cb = repmat([0.5,0.5,0.5],2,1);
fig_savename = 'SuppleFig12A_right_Exp1';
for fig=21
    %% 
    f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig);clf;
    h=daboxplot([AvgTFR_ML_Exp1(:,[1:2])],[repmat(1,size(AvgTFR_ML_Exp1,1),1)],...
        'scatter' ,2,...
        'color',cb,...
        'boxalpha', 0.5,...
        'jitter',1,...
        'linkline',0,...
        'withinlines',1,...
        'xtlabels',{'More','Fewer'})
    
    ylabel('Power to baseline [dB]')
    xlabel('Saccade')
   
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
    
   
    
    if plot_save
    saveas(gcf,[figsavepath,fig_savename,'.emf'])
    saveas(gcf,[figsavepath,fig_savename,'.png'])
    end
end 

% Experiment 2
fig_savename = 'SuppleFig12B_left_Exp2';
for fig=22
    %% 
    f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig);clf;
    h=daboxplot([AvgTFR_RF_Exp2(:,[1:2])],[repmat(1,size(AvgTFR_RF_Exp2,1),1)],...
        'scatter' ,2,...
        'color',cb,...
        'boxalpha', 0.5,...
        'jitter',1,...
        'linkline',0,...
        'withinlines',1,...
        'xtlabels',{'Rem','Forg'})
    
    ylabel('Power to baseline [dB]')
    xlabel('Memory')
   ylim([-6,1])
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
    
   
    
    if plot_save
    saveas(gcf,[figsavepath,fig_savename,'.emf'])
    saveas(gcf,[figsavepath,fig_savename,'.png'])
    end
end 

sizefig = [0.1,0.1,0.2,0.6];
cb = repmat([0.5,0.5,0.5],2,1);
fig_savename = 'SuppleFig12B_right_Exp2';
for fig=22
    %% 
    f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig);clf;
    h=daboxplot([AvgTFR_ML_Exp2(:,[1:2])],[repmat(1,size(AvgTFR_ML_Exp2,1),1)],...
        'scatter' ,2,...
        'color',cb,...
        'boxalpha', 0.5,...
        'jitter',1,...
        'linkline',0,...
        'withinlines',1,...
        'xtlabels',{'More','Fewer'})
    
    ylabel('Power to baseline [dB]')
    xlabel('Saccade')
   ylim([-6,1])
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
    
   
    
    if plot_save
    saveas(gcf,[figsavepath,fig_savename,'.emf'])
    saveas(gcf,[figsavepath,fig_savename,'.png'])
    end
end 

% Experiment 3
fig_savename = 'SuppleFig12C_left_Exp3';
for fig=23
    %% 
    f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig);clf;
    h=daboxplot([AvgTFR_RF_Exp3(:,[1:2])],[repmat(1,size(AvgTFR_RF_Exp3,1),1)],...
        'scatter' ,2,...
        'color',cb,...
        'boxalpha', 0.5,...
        'jitter',1,...
        'linkline',0,...
        'withinlines',1,...
        'xtlabels',{'Rem','Forg'})
    
    ylabel('Power to baseline [dB]')
    xlabel('Memory')
   ylim([-6,0])
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
    
   
    
    if plot_save
    saveas(gcf,[figsavepath,fig_savename,'.emf'])
    saveas(gcf,[figsavepath,fig_savename,'.png'])
    end
end 

sizefig = [0.1,0.1,0.2,0.6];
cb = repmat([0.5,0.5,0.5],2,1);
fig_savename = 'SuppleFig12C_right_Exp3';
for fig=23
    %% 
    f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig);clf;
    h=daboxplot([AvgTFR_ML_Exp3(:,[1:2])],[repmat(1,size(AvgTFR_ML_Exp3,1),1)],...
        'scatter' ,2,...
        'color',cb,...
        'boxalpha', 0.5,...
        'jitter',1,...
        'linkline',0,...
        'withinlines',1,...
        'xtlabels',{'More','Fewer'})
    
    ylabel('Power to baseline [dB]')
    xlabel('Saccade')
   ylim([-6,0])
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
    
   
    
    if plot_save
    saveas(gcf,[figsavepath,fig_savename,'.emf'])
    saveas(gcf,[figsavepath,fig_savename,'.png'])
    end
end 

% Experiment 4
fig_savename = 'SuppleFig12D_left_Exp4';
for fig=24
    %% 
    f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig);clf;
    h=daboxplot([AvgTFR_RF_Exp4(:,[1:2])],[repmat(1,size(AvgTFR_RF_Exp4,1),1)],...
        'scatter' ,2,...
        'color',cb,...
        'boxalpha', 0.5,...
        'jitter',1,...
        'linkline',0,...
        'withinlines',1,...
        'xtlabels',{'Rem','Forg'})
    
    ylabel('Power to baseline [dB]')
    xlabel('Memory')
   ylim([-7.5,0])
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
    
   
    
    if plot_save
    saveas(gcf,[figsavepath,fig_savename,'.emf'])
    saveas(gcf,[figsavepath,fig_savename,'.png'])
    end
end 

sizefig = [0.1,0.1,0.2,0.6];
cb = repmat([0.5,0.5,0.5],2,1);
fig_savename = 'SuppleFig12D_right_Exp4';
for fig=24
    %% 
    f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig);clf;
    h=daboxplot([AvgTFR_ML_Exp4(:,[1:2])],[repmat(1,size(AvgTFR_ML_Exp4,1),1)],...
        'scatter' ,2,...
        'color',cb,...
        'boxalpha', 0.5,...
        'jitter',1,...
        'linkline',0,...
        'withinlines',1,...
        'xtlabels',{'More','Fewer'})
    
    ylabel('Power to baseline [dB]')
    xlabel('Saccade')
   ylim([-7.5,0])
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
    
   
    
    if plot_save
    saveas(gcf,[figsavepath,fig_savename,'.emf'])
    saveas(gcf,[figsavepath,fig_savename,'.png'])
    end
end 

