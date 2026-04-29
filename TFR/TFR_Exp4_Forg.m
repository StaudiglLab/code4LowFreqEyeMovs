%% Control Time Frequency analyis for Exp 4
% compare Forg  1dva vs 21
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


%% Analysis
rng(2026);
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
    %% Sep by condition

    % Forg 21 dva vs 1
    cfg= [];
    cfg.trials =  tfrfinal.trialinfo(:,3)==3 & tfrfinal.trialinfo(:,2)==0;
    A  = ft_selectdata(cfg, tfrfinal);
    cfg.trials = tfrfinal.trialinfo(:,3)==1 & tfrfinal.trialinfo(:,2)==0;
    B  = ft_selectdata(cfg, tfrfinal);

    [tfr_C3,tfr_C4,n_trl] = TFR_balance_trial(A,B,100);
    n_trl_C3C4(s,1:2) = n_trl;

    tfr_C3_raw{s} = tfr_C3;
    tfr_C4_raw{s} = tfr_C4;

    clearvars A B tfr_C3 tfr_C4

end

%% baseline correction 
baseline_t = [-1,-0.5];

for s = 1:length(tfr_C3_raw)
cfg = [];
    cfg.baseline =baseline_t;
    cfg.baselinetype = 'db';
    
 
    
    % Forg 21 vs 1 dva
    tfr_C3_ave{s} = ft_freqbaseline(cfg,tfr_C3_raw{s});
    tfr_C4_ave{s} = ft_freqbaseline(cfg,tfr_C4_raw{s});
    

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

[stat_For21_1]  = ft_freqstatistics(cfg, tfr_C3_ave{subset},tfr_C4_ave{subset});%%

%% save stat

stat_For21_1 = rmfield(stat_For21_1,'cfg');
save([datafigspath,'SuppleFig11.mat'],'stat_For21_1')


%%

load([datapath,'layout.mat']) % load layout
load([datafigspath,'SuppleFig11.mat'])



plot_save = 1; % to save or not figure

sizefig1 = [0.2,0.2,0.6,0.7];
sizefig2 = [0.2,0.2,0.35,0.4];
fontaxis = 24;line_width = 2;

fig_savename1 = 'SuppleFig11_1';
fig_savename2 = 'SuppleFig11_2';
for forg=1
    sig_channel = sum(sum(stat_For21_1.mask,3),2)>0;
    stat_For21_1_ave.dimord = 'chan_freq_time';
    stat_For21_1_ave.time   = stat_For21_1.time;
    stat_For21_1_ave.freq   = stat_For21_1.freq;
    stat_For21_1_ave.maskTF =(sum(stat_For21_1.mask,1))>0;
    stat_For21_1_ave.stat =mean(stat_For21_1.stat(sig_channel,:,:),1);
    stat_For21_1_ave.label = {'average'};
    
    cfg = [];
    cfg.channel = {'average'};
    cfg.parameter = 'stat';
    cfg.maskparameter = 'maskTF';
    cfg.maskstyle        = 'outline';
    cfg.comment        = 'no';
    cfg.shading             = 'interp';
    cfg.zlim = [-3 3];
    f=figure('Name',int2str(11),'units','normalized','outerposition',sizefig1);clf;
    ft_singleplotTFR(cfg,stat_For21_1_ave);
    colormap(flipud(cbrewer2('RdBu','div')));
    title('TFR 21 vs 1 dva (Forg)')
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
    for chan=1:size(stat_For21_1.stat,1)
        chan_t = stat_For21_1.stat(chan,:,:);
        topo_t(chan,1) = mean(chan_t(sum(stat_For21_1.mask,1)>0));
    end
    
    topo_stat_For21_1.stat = topo_t;
    topo_stat_For21_1.dimord='chan_time';
    topo_stat_For21_1.time=0;
    topo_stat_For21_1.label=stat_For21_1.label;
    
    cfg         = [];
    cfg.parameter = 'stat';
    cfg.marker        = 'off';
    cfg.zlim = 'absmax';
    cfg.comment        = 'no';
    cfg.zlim = [-3,3];
    cfg.layout  = lay;
    cfg.colormap = '*RdBu';
    
    cfg.highlight          = 'on';
    cfg.highlightchannel = {stat_For21_1.label{sig_channel}};
    cfg.highlightsymbol   = '.';
    cfg.highlightsize      = 10;
    cfg.comment            = 'no';
    cfg.shading             = 'interp';
    cfg.style               ='straight';
    f=figure('Name',int2str(12),'units','normalized','outerposition',sizefig2);clf;
    ft_topoplotTFR(cfg,topo_stat_For21_1);
    c = colorbar;
    c.Label.String = 't value';
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',fontaxis,'FontWeight','Bold', 'LineWidth', line_width);
    if plot_save
        saveas(gcf,[figsavepath,fig_savename2,'.emf'])
        saveas(gcf,[figsavepath,fig_savename2,'.png'])
    end
end
