%% Time Frequency analyis for Exp 4

codepath = '';
addpath(genpath( [codepath,'\subfunctions\']))


datapath = [codepath,'\datafiles\'];
datafigspath =  [codepath,'\data4figs\'];
figsavepath = [codepath,'\figures\'];
subinfopath = [codepath,'\datafiles\Subjects_Exp4\ExpInfo\'];
sacsavepath = [codepath,'\datafiles\Subjects_Exp4\saccades\'];
eegsavepath = [codepath,'\datafiles\Subjects_Exp4\EEG\'];
TFRsavepath = [codepath,'\datafiles\Subjects_Exp4\TFR\'];
% behavioral results
savename = [datapath,'Beh_Exp4.mat'];
load(savename)

%% Compute TFR 
for s = 1:length(Mem)
    
    disp(s)
    %% load EEG file of each subject
    savename = [eegsavepath,sprintf('Sub_eegICA_PV_%d_%s.mat',s,'Exp4')];
    load(savename);
     
    eeg = eeg_final_finetuned;
     
    %% add beh
    eeg.trialinfo(:,2) = Mem{s}(eeg.trialinfo(:,1));
    eeg.trialinfo(:,3) = Confi{s}(eeg.trialinfo(:,1));

    %% do time freq analysis
    cfg              = [];
    cfg.latency      = [-1.499,4.5];
    eeg = ft_selectdata(cfg,eeg);
    
    
    cfg              = [];
    cfg.output       = 'pow';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';
    cfg.foi          = 2:2:40;                         % default : analysis 2 to 30 Hz in steps of 2 Hz
    cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % default : length of time window = 0.5 sec 
    cfg.toi          = -1:0.05:4;
    cfg.keeptrials = 'yes';
    tfralltrl = ft_freqanalysis(cfg, eeg);
    
    %% save data

    savename = [TFRsavepath,sprintf('TFR_PVstim_allTrials_sub%d_%s.mat',s,'Exp4')];
    save(savename,'tfralltrl','-v7.3')

end

%% Analysis
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
    
    %% by Viewing conditions x Mem
    
    cfg= [];
    cfg.trials = tfrfinal.trialinfo(:,3)==1 & tfrfinal.trialinfo(:,2)==1;
    n_trl_ViewC_RF(s,1) = sum(cfg.trials);    
    tfr_VC1Rem_raw{s}  = ft_freqdescriptives(cfg, tfrfinal);
    cfg.trials = tfrfinal.trialinfo(:,3)==2 & tfrfinal.trialinfo(:,2)==1;
    n_trl_ViewC_RF(s,2) = sum(cfg.trials);
    tfr_VC2Rem_raw{s}  = ft_freqdescriptives(cfg,tfrfinal);
    cfg.trials = tfrfinal.trialinfo(:,3)==3 & tfrfinal.trialinfo(:,2)==1;
    n_trl_ViewC_RF(s,3) = sum(cfg.trials);
    tfr_VC3Rem_raw{s}  = ft_freqdescriptives(cfg,tfrfinal);
    
    cfg.trials = tfrfinal.trialinfo(:,3)==1 & tfrfinal.trialinfo(:,2)==0;
    n_trl_ViewC_RF(s,4) = sum(cfg.trials);    
    tfr_VC1For_raw{s}  = ft_freqdescriptives(cfg, tfrfinal);
    cfg.trials = tfrfinal.trialinfo(:,3)==2 & tfrfinal.trialinfo(:,2)==0;
    n_trl_ViewC_RF(s,5) = sum(cfg.trials);
    tfr_VC2For_raw{s}  = ft_freqdescriptives(cfg,tfrfinal);
    cfg.trials = tfrfinal.trialinfo(:,3)==3 & tfrfinal.trialinfo(:,2)==0;
    n_trl_ViewC_RF(s,6) = sum(cfg.trials);
    tfr_VC3For_raw{s}  = ft_freqdescriptives(cfg,tfrfinal);
    
   
    
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
    
    
    % sep by each condition
    cfg= [];
    cfg.trials = tfrfinal.trialinfo(:,5)>=cri1 & tfrfinal.trialinfo(:,3)==1;
    n_trl_FixExp_EIVC(s,1) = sum(cfg.trials);
    tfr_Fix_raw_EI_VC1{s}  = ft_freqdescriptives(cfg, tfrfinal);
    cfg.trials = tfrfinal.trialinfo(:,5)>=cri2 & tfrfinal.trialinfo(:,3)==2;
    n_trl_FixExp_EIVC(s,2) = sum(cfg.trials);
    tfr_Fix_raw_EI_VC2{s}  = ft_freqdescriptives(cfg, tfrfinal);
    cfg.trials = tfrfinal.trialinfo(:,5)>=cri3 & tfrfinal.trialinfo(:,3)==3;
    n_trl_FixExp_EIVC(s,3) = sum(cfg.trials);
    tfr_Fix_raw_EI_VC3{s}  = ft_freqdescriptives(cfg, tfrfinal);
    
    cfg.trials = tfrfinal.trialinfo(:,5)<=cri1 & tfrfinal.trialinfo(:,3)==1;
    n_trl_FixExp_EIVC(s,4) = sum(cfg.trials);
    tfr_Explr_raw_EI_VC1{s}  = ft_freqdescriptives(cfg, tfrfinal);
    cfg.trials = tfrfinal.trialinfo(:,5)<=cri2 & tfrfinal.trialinfo(:,3)==2;
    n_trl_FixExp_EIVC(s,5) = sum(cfg.trials);
    tfr_Explr_raw_EI_VC2{s}  = ft_freqdescriptives(cfg, tfrfinal);
    cfg.trials = tfrfinal.trialinfo(:,5)<=cri3 & tfrfinal.trialinfo(:,3)==3;
    n_trl_FixExp_EIVC(s,6) = sum(cfg.trials);
    tfr_Explr_raw_EI_VC3{s}  = ft_freqdescriptives(cfg, tfrfinal);
    
    clearvars cri*

end

%% baseline correction 
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
    
    
    %% Mem x View condition
    tfr_V1C1_ave{s} = ft_freqbaseline(cfg,tfr_VC1Rem_raw{s});
    tfr_V1C2_ave{s} = ft_freqbaseline(cfg,tfr_VC1For_raw{s});
    
    tfr_V2C1_ave{s} = ft_freqbaseline(cfg,tfr_VC2Rem_raw{s});
    tfr_V2C2_ave{s} = ft_freqbaseline(cfg,tfr_VC2For_raw{s});
    
    tfr_V3C1_ave{s} = ft_freqbaseline(cfg,tfr_VC3Rem_raw{s});
    tfr_V3C2_ave{s} = ft_freqbaseline(cfg,tfr_VC3For_raw{s});
    
    %% Exploration index x View condition
    
    tfr_Fix_EI_VC1_ave{s}   = ft_freqbaseline(cfg,tfr_Fix_raw_EI_VC1{s});
    tfr_Explr_EI_VC1_ave{s} = ft_freqbaseline(cfg,tfr_Explr_raw_EI_VC1{s});
    
    tfr_Fix_EI_VC2_ave{s}   = ft_freqbaseline(cfg,tfr_Fix_raw_EI_VC2{s});
    tfr_Explr_EI_VC2_ave{s} = ft_freqbaseline(cfg,tfr_Explr_raw_EI_VC2{s});
    
    tfr_Fix_EI_VC3_ave{s}   = ft_freqbaseline(cfg,tfr_Fix_raw_EI_VC3{s});
    tfr_Explr_EI_VC3_ave{s} = ft_freqbaseline(cfg,tfr_Explr_raw_EI_VC3{s});
   
 
end

%% Statistics
rng(2025);
load([datapath,'layout.mat'])
subset = 1:34;
design      = zeros(2,length(subset)*2);
design(1,:) = repmat(1:length(subset),[1 2]);
design(2,:) = [ones(1,length(subset)),ones(1,length(subset))+1];

% prepare neighbours
load([datapath,'saved_edited_neighbours.mat'])

%
cfg = [];
cfg.spmversion = 'spm12';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.ivar                = 2;
cfg.uvar                = 1;
cfg.parameter           = 'powspctrm';
cfg.latency             = [0,4];
cfg.frequency           = [2 40];
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 1000;
cfg.minnbchan    = 3;
cfg.neighbours=neighbours;
cfg.design           = design;

cfg.avgoverfreq = 'no';

[stat_mem]  = ft_freqstatistics(cfg, tfr_C1_ave{subset},tfr_C2_ave{subset});
[stat_ExpvsFix_EI]  = ft_freqstatistics(cfg, tfr_Exploration_ave_EI{subset},tfr_Fixation_ave_EI{subset});%%

%% save stat

stat_mem = rmfield(stat_mem,'cfg');% save space
stat_ExpvsFix_EI = rmfield(stat_ExpvsFix_EI,'cfg');
save([datafigspath,'Fig4D_left.mat'],'stat_mem')
save([datafigspath,'Fig4E_left.mat'],'stat_ExpvsFix_EI')
%% Extract posterior ab power for each condition
toiVM = [1,3];% focus on the cluster time window
toiEI = [1,4];
foi = [10,20];
coi = {'P*','O*'};

for s=1:length(tfr_V1C1_ave)
    disp(s)
    cfg=[];
    cfg.channel = coi;
    cfg.latency = toiVM;
    cfg.frequency = foi;
    
    %% Sep by Cond
    
    
    %% Sep by  Mem % Viewing condition
    
    
    selec_dat = ft_selectdata(cfg,tfr_V1C1_ave{s});
    avg_pow_VC_RF(s,1) = mean(selec_dat.powspctrm(:));clearvars selec_dat
    selec_dat = ft_selectdata(cfg,tfr_V1C2_ave{s});
    avg_pow_VC_RF(s,2) = mean(selec_dat.powspctrm(:));clearvars selec_dat
    selec_dat = ft_selectdata(cfg,tfr_V2C1_ave{s});
    avg_pow_VC_RF(s,3) = mean(selec_dat.powspctrm(:));clearvars selec_dat
    selec_dat = ft_selectdata(cfg,tfr_V2C2_ave{s});
    avg_pow_VC_RF(s,4) = mean(selec_dat.powspctrm(:));clearvars selec_dat
    selec_dat = ft_selectdata(cfg,tfr_V3C1_ave{s});
    avg_pow_VC_RF(s,5) = mean(selec_dat.powspctrm(:));clearvars selec_dat
    selec_dat = ft_selectdata(cfg,tfr_V3C2_ave{s});
    avg_pow_VC_RF(s,6) = mean(selec_dat.powspctrm(:));clearvars selec_dat
    
    
    %% Exploration Index x Viewing Condition
    
    cfg.latency = toiEI;
    
    
    selec_dat = ft_selectdata(cfg,tfr_Fix_EI_VC1_ave{s});
    avg_pow_VC_EI(s,2) = mean(selec_dat.powspctrm(:));clearvars selec_dat
    selec_dat = ft_selectdata(cfg,tfr_Explr_EI_VC1_ave{s});
    avg_pow_VC_EI(s,1) = mean(selec_dat.powspctrm(:));clearvars selec_dat
    selec_dat = ft_selectdata(cfg,tfr_Fix_EI_VC2_ave{s});
    avg_pow_VC_EI(s,4) = mean(selec_dat.powspctrm(:));clearvars selec_dat
    selec_dat = ft_selectdata(cfg,tfr_Explr_EI_VC2_ave{s});
    avg_pow_VC_EI(s,3) = mean(selec_dat.powspctrm(:));clearvars selec_dat
    selec_dat = ft_selectdata(cfg,tfr_Fix_EI_VC3_ave{s});
    avg_pow_VC_EI(s,6) = mean(selec_dat.powspctrm(:));clearvars selec_dat
    selec_dat = ft_selectdata(cfg,tfr_Explr_EI_VC3_ave{s});
    avg_pow_VC_EI(s,5) = mean(selec_dat.powspctrm(:));clearvars selec_dat
end

save([datafigspath,'Fig4D_right.mat'],'avg_pow_VC_RF')
save([datafigspath,'Fig4E_right.mat'],'avg_pow_VC_EI')

%% Figures
codepath = '';

datapath = [codepath,'\datafiles\'];
datafigspath =  [codepath,'\data4figs\'];
figsavepath = [codepath,'\figures\'];
load([datapath,'layout.mat']) % load layout
% load data
load( [datafigspath,'Fig4D_left.mat'])
load( [datafigspath,'Fig4E_left.mat'])

plot_save = 1; % to save or not figure


sizefig1 = [0.2,0.2,0.6,0.7];
sizefig2 = [0.2,0.2,0.35,0.4];
fontaxis = 24;line_width = 2;


fig_savename1 = 'Fig4D_left1';
fig_savename2 = 'Fig4D_left2';

for remforg=1
    sig_channel = sum(sum(stat_mem.mask,3),2)>0;
    stat_mem_ave.dimord = 'chan_freq_time';
    stat_mem_ave.time   = stat_mem.time;
    stat_mem_ave.freq   = stat_mem.freq;
    stat_mem_ave.maskTF =(sum(stat_mem.mask,1))>0;
    stat_mem_ave.stat =mean(stat_mem.stat(sig_channel,:,:),1);
    stat_mem_ave.label = {'average'};
    
    cfg = [];
    cfg.channel = {'average'};
    cfg.parameter = 'stat';
    cfg.maskparameter = 'maskTF';
    cfg.maskstyle        = 'outline';
    cfg.comment        = 'no';
    cfg.shading             = 'interp';
    cfg.zlim = [-2.5 2.5];
    f=figure('Name',int2str(11),'units','normalized','outerposition',sizefig1);clf;
    ft_singleplotTFR(cfg,stat_mem_ave);
    colormap(flipud(cbrewer2('RdBu','div')));
    title('TFR Rem vs Forg')
    ylabel('Frequency (Hz)')
    xlabel('Time [Sec]')
    c = colorbar;
    c.Label.String = 't value';
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',fontaxis,'FontWeight','Bold', 'LineWidth', line_width);
    
    if plot_save
        saveas(gcf,[figsavepath,fig_savename1,'.emf'])
        saveas(gcf,[figsavepath,fig_savename1,'.png'])
    end
    
    % topoplot
    for chan=1:size(stat_mem.stat,1)
        chan_t = stat_mem.stat(chan,:,:);
        topo_t(chan,1) = mean(chan_t(sum(stat_mem.mask,1)>0));
    end
    
    topo_stat_mem.stat = topo_t;
    topo_stat_mem.dimord='chan_time';
    topo_stat_mem.time=0;
    topo_stat_mem.label=stat_mem.label;
    
    cfg         = [];
    cfg.parameter = 'stat';
    cfg.marker        = 'off';
    cfg.zlim = 'absmax';
    cfg.comment        = 'no';
    cfg.zlim = [-2.5,2.5];
    cfg.layout  = lay;
    cfg.colormap = '*RdBu';
    
    cfg.highlight          = 'on';
    cfg.highlightchannel = {stat_mem.label{sig_channel}};
    cfg.highlightsymbol   = '.';
    cfg.highlightsize      = 10;
    cfg.comment            = 'no';
    cfg.shading             = 'interp';
    cfg.style               ='straight';
    f=figure('Name',int2str(12),'units','normalized','outerposition',sizefig2);clf;
    ft_topoplotTFR(cfg,topo_stat_mem);
    c = colorbar;
    c.Label.String = 't value';
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',fontaxis,'FontWeight','Bold', 'LineWidth', line_width);
    if plot_save
        saveas(gcf,[figsavepath,fig_savename2,'.emf'])
        saveas(gcf,[figsavepath,fig_savename2,'.png'])
    end
end



fig_savename1 = 'Fig4E_left1';
fig_savename2 = 'Fig4E_left2';

for EI=1
    sig_channel = sum(sum(stat_ExpvsFix_EI.mask,3),2)>0;
    stat_ExpvsFix_EI_ave.dimord = 'chan_freq_time';
    stat_ExpvsFix_EI_ave.time   = stat_ExpvsFix_EI.time;
    stat_ExpvsFix_EI_ave.freq   = stat_ExpvsFix_EI.freq;
    stat_ExpvsFix_EI_ave.maskTF =(sum(stat_ExpvsFix_EI.mask,1))>0;
    stat_ExpvsFix_EI_ave.stat =mean(stat_ExpvsFix_EI.stat(sig_channel,:,:),1);
    stat_ExpvsFix_EI_ave.label = {'average'};
    
    cfg = [];
    cfg.channel = {'average'};
    cfg.parameter = 'stat';
    cfg.maskparameter = 'maskTF';
    cfg.maskstyle        = 'outline';
    cfg.comment        = 'no';
    cfg.shading             = 'interp';
    cfg.zlim = [-2.5 2.5];
    f=figure('Name',int2str(11),'units','normalized','outerposition',sizefig1);clf;
    ft_singleplotTFR(cfg,stat_ExpvsFix_EI_ave);
    colormap(flipud(cbrewer2('RdBu','div')));
    title('TFR More vs Less Exploration')
    ylabel('Frequency (Hz)')
    xlabel('Time [Sec]')
    c = colorbar;
    c.Label.String = 't value';
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',fontaxis,'FontWeight','Bold', 'LineWidth', line_width);
    
    if plot_save
        saveas(gcf,[figsavepath,fig_savename1,'.emf'])
        saveas(gcf,[figsavepath,fig_savename1,'.png'])
    end
    
    % topoplot
    for chan=1:size(stat_ExpvsFix_EI.stat,1)
        chan_t = stat_ExpvsFix_EI.stat(chan,:,:);
        topo_t(chan,1) = mean(chan_t(sum(stat_ExpvsFix_EI.mask,1)>0));
    end
    
    topo_stat_ExpvsFix_EI.stat = topo_t;
    topo_stat_ExpvsFix_EI.dimord='chan_time';
    topo_stat_ExpvsFix_EI.time=0;
    topo_stat_ExpvsFix_EI.label=stat_ExpvsFix_EI.label;
    
    cfg         = [];
    cfg.parameter = 'stat';
    cfg.marker        = 'off';
    cfg.zlim = 'absmax';
    cfg.comment        = 'no';
    cfg.zlim = [-2.5,2.5];
    cfg.layout  = lay;
    cfg.colormap = '*RdBu';
    
    cfg.highlight          = 'on';
    cfg.highlightchannel = {stat_ExpvsFix_EI.label{sig_channel}};
    cfg.highlightsymbol   = '.';
    cfg.highlightsize      = 10;
    cfg.comment            = 'no';
    cfg.shading             = 'interp';
    cfg.style               ='straight';
    f=figure('Name',int2str(12),'units','normalized','outerposition',sizefig2);clf;
    ft_topoplotTFR(cfg,topo_stat_ExpvsFix_EI);
    c = colorbar;
    c.Label.String = 't value';
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',fontaxis,'FontWeight','Bold', 'LineWidth', line_width);
    if plot_save
        saveas(gcf,[figsavepath,fig_savename2,'.emf'])
        saveas(gcf,[figsavepath,fig_savename2,'.png'])
    end
end



% load data
load( [datafigspath,'Fig4D_right.mat'])
load( [datafigspath,'Fig4E_right.mat'])

sizefig = [0.2,0.2,0.25,0.6];fontaxis = 24;line_width = 2;
cb = cbrewer2('RdBu','div',2);
fig_savename = 'Fig4D_right';
for fig=231
    %% violin plot
    f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig);clf;
    
    h=daboxplot([avg_pow_VC_RF(:,[1,3,5]);avg_pow_VC_RF(:,[2,4,6])],...
        [repmat(1,size(avg_pow_VC_RF,1),1);repmat(2,size(avg_pow_VC_RF,1),1)],...
        'scatter' ,2,...
        'color',cb,...
        'boxalpha', 0.5,...
        'jitter',1,...
        'linkline',1,...
        'withinlines',0,...
        'legend',{'Rem','Forg'},...
        'xtlabels',{'1DVA','5DVA','21DVA'})
    
    h.lg.Location = 'northeast';
    
    ylabel('Power to baseline')
    xlabel('Viewing Condition')
    ylim([-9,1])
    yline(0,'-')
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
    if plot_save
    saveas(gcf,[figsavepath,fig_savename,'.emf'])
    saveas(gcf,[figsavepath,fig_savename,'.png'])
    end
end 

cb = cbrewer2('BrBG','div',2);

fig_savename = 'Fig4E_right';
for fig=231
    %% violin plot
    f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig);clf;
    
    h=daboxplot([avg_pow_VC_EI(:,[1,3,5]);avg_pow_VC_EI(:,[2,4,6])],...
        [repmat(1,size(avg_pow_VC_EI,1),1);repmat(2,size(avg_pow_VC_EI,1),1)],...
        'scatter' ,2,...
        'color',cb,...
        'boxalpha', 0.5,...
        'jitter',1,...
        'linkline',1,...
        'withinlines',0,...
        'legend',{'More Exploration','Less Exploration'},...
        'xtlabels',{'1DVA','5DVA','21DVA'})
    
    
    h.lg.Location = 'northeast';
    yline(0,'-')
    ylabel('Power to baseline')
    xlabel('Viewing Condition')
    ylim([-8,1])
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
    if plot_save
    saveas(gcf,[figsavepath,fig_savename,'.emf'])
    saveas(gcf,[figsavepath,fig_savename,'.png'])
    end
end 
