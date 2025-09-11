%% Time Frequency analyis for Exp 1

codepath = '';
addpath(genpath( [codepath,'\subfunctions\']))


datapath = [codepath,'\datafiles\'];
datafigspath =  [codepath,'\data4figs\'];
figsavepath = [codepath,'\figures\'];
subinfopath = [codepath,'\datafiles\Subjects_Exp1\ExpInfo\'];
sacsavepath = [codepath,'\datafiles\Subjects_Exp1\saccades\'];
eegsavepath = [codepath,'\datafiles\Subjects_Exp1\EEG\'];
TFRsavepath = [codepath,'\datafiles\Subjects_Exp1\TFR\'];
% behavioral results
savename = [datapath,'Beh_Exp1.mat'];
load(savename)

%% Compute TFR 
for s = 1:length(Mem)
    
    disp(s)
    %% load EEG file of each subject
    savename = [eegsavepath,sprintf('Sub_eegICA_PV_%d_%s.mat',s,'Exp1')];
    load(savename);
     
    eeg = eeg_final;
     
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

    savename = [TFRsavepath,sprintf('TFR_PVstim_allTrials_sub%d_%s.mat',s,'Exp1')];
    save(savename,'tfralltrl','-v7.3')

end

%% Analysis
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
%% statistics 
rng(2025);
load([datapath,'layout.mat'])
subset = 1:20;
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
[stat_sac]  = ft_freqstatistics(cfg, tfr_Msac_ave{subset},tfr_Lsac_ave{subset});

%% save outputs
stat_mem = rmfield(stat_mem,'cfg');% save space
stat_sac = rmfield(stat_sac,'cfg');
save([datafigspath,'Fig1C_left.mat'],'stat_mem')
save([datafigspath,'Fig1E_left.mat'],'stat_sac')


%% Figure
% load data for plotting
%codepath = '';

datapath = [codepath,'\datafiles\'];
datafigspath =  [codepath,'\data4figs\'];
figsavepath = [codepath,'\figures\'];
load([datapath,'layout.mat']) % load layout
% load data
load( [datafigspath,'Fig1C_left.mat'])
load( [datafigspath,'Fig1E_left.mat'])


plot_save = 1; % to save or not figure


sizefig1 = [0.2,0.2,0.6,0.7];
sizefig2 = [0.2,0.2,0.35,0.4];
fontaxis = 24;line_width = 2;
fig_savename1 = 'Fig1C_left1';
fig_savename2 = 'Fig1C_left2';

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

fig_savename1 = 'Fig1E_left1';
fig_savename2 = 'Fig1E_left2';

for nsac=1
    sig_channel = sum(sum(stat_sac.mask,3),2)>0;
    stat_sac_ave.dimord = 'chan_freq_time';
    stat_sac_ave.time   = stat_sac.time;
    stat_sac_ave.freq   = stat_sac.freq;
    stat_sac_ave.maskTF =(sum(stat_sac.mask,1))>0;
    stat_sac_ave.stat =mean(stat_sac.stat(sig_channel,:,:),1);
    stat_sac_ave.label = {'average'};
    
    cfg = [];
    cfg.channel = {'average'};
    cfg.parameter = 'stat';
    cfg.maskparameter = 'maskTF';
    cfg.maskstyle        = 'outline';
    cfg.comment        = 'no';
    cfg.shading             = 'interp';
    cfg.zlim = [-2.5 2.5];
    f=figure('Name',int2str(11),'units','normalized','outerposition',sizefig1);clf;
    ft_singleplotTFR(cfg,stat_sac_ave);
    colormap(flipud(cbrewer2('RdBu','div')));
    title('TFR More vs Less Saccades')
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
    for chan=1:size(stat_sac.stat,1)
        chan_t = stat_sac.stat(chan,:,:);
        topo_t(chan,1) = mean(chan_t(sum(stat_sac.mask,1)>0));
    end
    
    topo_stat_sac.stat = topo_t;
    topo_stat_sac.dimord='chan_time';
    topo_stat_sac.time=0;
    topo_stat_sac.label=stat_sac.label;
    
    cfg         = [];
    cfg.parameter = 'stat';
    cfg.marker        = 'off';
    cfg.zlim = 'absmax';
    cfg.comment        = 'no';
    cfg.zlim = [-2.5,2.5];
    cfg.layout  = lay;
    cfg.colormap = '*RdBu';
    
    cfg.highlight          = 'on';
    cfg.highlightchannel = {stat_sac.label{sig_channel}};
    cfg.highlightsymbol   = '.';
    cfg.highlightsize      = 10;
    cfg.comment            = 'no';
    cfg.shading             = 'interp';
    cfg.style               ='straight';
    f=figure('Name',int2str(12),'units','normalized','outerposition',sizefig2);clf;
    ft_topoplotTFR(cfg,topo_stat_sac);
    c = colorbar;
    c.Label.String = 't value';
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',fontaxis,'FontWeight','Bold', 'LineWidth', line_width);
    if plot_save
        saveas(gcf,[figsavepath,fig_savename2,'.emf'])
        saveas(gcf,[figsavepath,fig_savename2,'.png'])
    end
end

