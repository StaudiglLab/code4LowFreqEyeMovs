%% Time Frequency analyis for Exp 2

codepath = '';
addpath(genpath( [codepath,'\subfunctions\']))


datapath = [codepath,'\datafiles\'];
datafigspath =  [codepath,'\data4figs\'];
figsavepath = [codepath,'\figures\'];
subinfopath = [codepath,'\datafiles\Subjects_Exp2\ExpInfo\'];
sacsavepath = [codepath,'\datafiles\Subjects_Exp2\saccades\'];
eegsavepath = [codepath,'\datafiles\Subjects_Exp2\EEG\'];
TFRsavepath = [codepath,'\datafiles\Subjects_Exp2\TFR\'];
% behavioral results
savename = [datapath,'Beh_Exp2.mat'];
load(savename)

%% Compute TFR 
for s = 1:length(Mem)
    
    disp(s)
    %% load EEG file of each subject
    savename = [eegsavepath,sprintf('Sub_Posterior_PV_%d_%s.mat',s,'Exp2')];
    load(savename);
     
    cfg              = [];
    cfg.latency      = [-1.499,4.5];
    data_all = ft_selectdata(cfg,data_all);
    
    
    cfg              = [];
    cfg.output       = 'pow';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';
    cfg.foi          = 2:2:40;                         % analysis 2 to 30 Hz in steps of 2 Hz
    cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
    cfg.toi          = -1:0.05:4;
    cfg.keeptrials = 'yes';
    tfralltrl{s} = ft_freqanalysis(cfg, data_all);
    tfralltrl{s}.elecinfo = data_all.Chanwm;

    spike_mark_sub{s} = spike_trl_idx_all;
    trial_rej_sub{s} = trial_rej_all;


end
save([TFRsavepath,'TFR_PVstim_allTrials_allsub_Exp2.mat'],'tfralltrl','spike_mark_sub','trial_rej_sub','-v7.3');
%% Separate conditions
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
%% stat
subset = 1:20;
design      = zeros(2,length(subset)*2);
design(1,:) = repmat(1:length(subset),[1 2]);
design(2,:) = [ones(1,length(subset)),ones(1,length(subset))+1];

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

cfg.design           = design;

cfg.avgoverfreq = 'no';

[stat_mem]  = ft_freqstatistics(cfg, TFR_enc_bc_rem{subset},TFR_enc_bc_for{subset});
[stat_sac]  = ft_freqstatistics(cfg, TFR_enc_bc_moresac{subset},TFR_enc_bc_lesssac{subset});
%% save outputs
stat_mem = rmfield(stat_mem,'cfg');% save space
stat_sac = rmfield(stat_sac,'cfg');
save([datafigspath,'Fig2C_left.mat'],'stat_mem')
save([datafigspath,'Fig2D_left.mat'],'stat_sac')

%% Figure
% load data for plotting
codepath = '';

datapath = [codepath,'\datafiles\'];
datafigspath =  [codepath,'\data4figs\'];
figsavepath = [codepath,'\figures\'];
% load data
load( [datafigspath,'Fig2C_left.mat'])
load( [datafigspath,'Fig2D_left.mat'])


plot_save = 1; % to save or not figure

sizefig1 = [0.2,0.2,0.6,0.7];
sizefig2 = [0.2,0.2,0.35,0.4];
fontaxis = 24;line_width = 2;
fig_savename1 = 'Fig2C_left';
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
    cfg.zlim = [-3 3];
    f=figure('Name',int2str(11),'units','normalized','outerposition',sizefig1);clf;
    ft_singleplotTFR(cfg,stat_mem_ave);
    colormap(flipud(cbrewer2('RdBu','div')));
    title('TFR Rem vs Forg')
    ylabel('Frequency (Hz)')
    c = colorbar;
    c.Label.String = 't value';
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',fontaxis,'FontWeight','Bold', 'LineWidth', line_width);
    
    if plot_save
        saveas(gcf,[figsavepath,fig_savename1,'.emf'])
        saveas(gcf,[figsavepath,fig_savename1,'.png'])
    end
    
    
end

fig_savename1 = 'Fig2D_left';
for sac=1
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
    cfg.zlim = [-3 3];
    f=figure('Name',int2str(11),'units','normalized','outerposition',sizefig1);clf;
    ft_singleplotTFR(cfg,stat_sac_ave);
    colormap(flipud(cbrewer2('RdBu','div')));
    title('TFR More vs Fewer Saccade')
    ylabel('Frequency (Hz)')
    c = colorbar;
    c.Label.String = 't value';
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',fontaxis,'FontWeight','Bold', 'LineWidth', line_width);
    
    if plot_save
        saveas(gcf,[figsavepath,fig_savename1,'.emf'])
        saveas(gcf,[figsavepath,fig_savename1,'.png'])
    end
    
    
end
