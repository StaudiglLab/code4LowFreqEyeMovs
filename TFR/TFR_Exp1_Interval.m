%% 
% TFR contrast with median spliting the saccade interval
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

%%
eye_hz = 600;
for s = 1:20
    disp(s)
    savename = [TFRsavepath,sprintf('TFR_PVstim_allTrials_sub%d_%s.mat',s,'Exp1')];
    load(savename,'tfralltrl')
    
    % load trial info
    load([sacsavepath,sprintf('TrialEyeInfo_%d_%s.mat',s,'Exp1')])
    %% calculate median saccade interval
    load([sacsavepath,sprintf('Sac_PV_Info_Sub%d_%s.mat',s,'Exp1')]);
   
    Eye_trlinfo(:,9)=4000;% defualt 4 sec of no saccade
    PV_sac = Sac_trials.trialinfo(Sac_trials.trialinfo(:,6)>1.3*eye_hz,:);% after stim onset
    trlidx = unique(PV_sac(:,1));
    for itrl= 1:length(trlidx)
        Sactrl = PV_sac(PV_sac(:,1)==trlidx(itrl),:);       
        Eye_trlinfo(Eye_trlinfo(:,1)==trlidx(itrl),9)=median(diff([Sactrl(:,6)]))/eye_hz*1000;% into msec
 
    end
    
 
    %% select trials
    % select only the intersect between clean EEG and clean eye trials
    cfg = [];
    cfg.trials = ismember(tfralltrl.trialinfo(:,1),Eye_trlinfo(Eye_trlinfo(:,8)==1,1));
    tfrfinal = ft_selectdata(cfg,tfralltrl);
    tfrfinal.trialinfo(:,4) = Eye_trlinfo(dsearchn(Eye_trlinfo(:,1),tfrfinal.trialinfo(:,1)),5);
    tfrfinal.trialinfo(:,5) = Eye_trlinfo(dsearchn(Eye_trlinfo(:,1),tfrfinal.trialinfo(:,1)),9);
    
    
   
    %% by median median saccade interval
    cri = nanmedian(tfrfinal.trialinfo(:,5));
    
    cfg= [];
    cfg.trials = tfrfinal.trialinfo(:,5)<=cri;
    avg_sacinterval(s,1) = mean(tfrfinal.trialinfo(cfg.trials,5));
    tfr_ShortInt_raw{s}  = ft_freqdescriptives(cfg, tfrfinal);
    cfg.trials = tfrfinal.trialinfo(:,5)>=cri;
    avg_sacinterval(s,2) = mean(tfrfinal.trialinfo(cfg.trials,5));
    tfr_LongInt_raw{s}  = ft_freqdescriptives(cfg,tfrfinal);
    clearvars cri


    

end

% baseline correction 
baseline_t = [-1,-0.5];

for s = 1:length(tfr_ShortInt_raw)
cfg = [];
    cfg.baseline =baseline_t;
    cfg.baselinetype = 'db';
    
    tfr_C1_ave{s} = ft_freqbaseline(cfg,tfr_ShortInt_raw{s});
    tfr_C2_ave{s} = ft_freqbaseline(cfg,tfr_LongInt_raw{s});
    
   
   
 
end
%% statistics 

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

[stat_SacInt]  = ft_freqstatistics(cfg, tfr_C1_ave{subset},tfr_C2_ave{subset});


%% save outputs

stat_SacInt = rmfield(stat_SacInt,'cfg');% save space

save([datafigspath,'SuppleFig16.mat'],'stat_SacInt','avg_sacinterval')

%%
%codepath = '';

datapath = [codepath,'\datafiles\'];
datafigspath =  [codepath,'\data4figs\'];
figsavepath = [codepath,'\figures\'];
load([datapath,'layout.mat']) % load layout
load([datafigspath,'SuppleFig16.mat'])


plot_save = 1; % to save or not figure


sizefig1 = [0.2,0.2,0.6,0.7];
sizefig2 = [0.2,0.2,0.35,0.4];
fontaxis = 24;line_width = 2;
fig_savename1 = 'SuppleFig16_2';
fig_savename2 = 'SuppleFig16_3';

for interval=1
    sig_channel = sum(sum(stat_SacInt.mask,3),2)>0;
    stat_SacInt_ave.dimord = 'chan_freq_time';
    stat_SacInt_ave.time   = stat_SacInt.time;
    stat_SacInt_ave.freq   = stat_SacInt.freq;
    stat_SacInt_ave.maskTF =(sum(stat_SacInt.mask,1))>0;
    stat_SacInt_ave.stat =mean(stat_SacInt.stat(sig_channel,:,:),1);
    stat_SacInt_ave.label = {'average'};
    
    cfg = [];
    cfg.channel = {'average'};
    cfg.parameter = 'stat';
    cfg.maskparameter = 'maskTF';
    cfg.maskstyle        = 'outline';
    cfg.comment        = 'no';
    cfg.shading             = 'interp';
    cfg.zlim = [-2.5 2.5];
    f=figure('Name',int2str(11),'units','normalized','outerposition',sizefig1);clf;
    ft_singleplotTFR(cfg,stat_SacInt_ave);
    colormap(flipud(cbrewer2('RdBu','div')));
    title('TFR Shorter vs Longer Interval')
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
    for chan=1:size(stat_SacInt.stat,1)
        chan_t = stat_SacInt.stat(chan,:,:);
        topo_t(chan,1) = mean(chan_t(sum(stat_SacInt.mask,1)>0));
    end
    
    topo_stat_SacInt.stat = topo_t;
    topo_stat_SacInt.dimord='chan_time';
    topo_stat_SacInt.time=0;
    topo_stat_SacInt.label=stat_SacInt.label;
    
    cfg         = [];
    cfg.parameter = 'stat';
    cfg.marker        = 'off';
    cfg.zlim = 'absmax';
    cfg.comment        = 'no';
    cfg.zlim = [-2.5,2.5];
    cfg.layout  = lay;
    cfg.colormap = '*RdBu';
    
    cfg.highlight          = 'on';
    cfg.highlightchannel = {stat_SacInt.label{sig_channel}};
    cfg.highlightsymbol   = '.';
    cfg.highlightsize      = 10;
    cfg.comment            = 'no';
    cfg.shading             = 'interp';
    cfg.style               ='straight';
    f=figure('Name',int2str(12),'units','normalized','outerposition',sizefig2);clf;
    ft_topoplotTFR(cfg,topo_stat_SacInt);
    c = colorbar;
    c.Label.String = 't value';
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',fontaxis,'FontWeight','Bold', 'LineWidth', line_width);
    if plot_save
        saveas(gcf,[figsavepath,fig_savename2,'.emf'])
        saveas(gcf,[figsavepath,fig_savename2,'.png'])
    end
end


%%
[cb] = cbrewer2('PRGn','div',12);

fontaxis = 16;line_width = 2;
fig_savename1 = 'SuppleFig16_1';

sizefig = [0.1,0.1,0.2,0.6];

for fig=21
    %% 
    f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig);clf;
    h=daboxplot([avg_sacinterval(:,[1:2])./1000],[repmat(1,size(avg_sacinterval,1),1)],...
        'scatter' ,2,...
        'color',cb,...
        'boxalpha', 0.5,...
        'jitter',1,...
        'linkline',0,...
        'withinlines',1,...
        'xtlabels',{'Shorter','Longer'})
    
    ylabel('Interval [Sec]')
   
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
    
   
    
    if plot_save
    saveas(gcf,[figsavepath,fig_savename1,'.emf'])
    saveas(gcf,[figsavepath,fig_savename1,'.png'])
    end
end 