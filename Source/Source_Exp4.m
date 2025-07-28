%% Source localisation for Exp 4

codepath = '';

datapath = [codepath,'\datafiles\'];
datafigspath =  [codepath,'\data4figs\'];
figsavepath = [codepath,'\figures\'];
subinfopath = [codepath,'\datafiles\Subjects_Exp4\ExpInfo\'];
sacsavepath = [codepath,'\datafiles\Subjects_Exp4\saccades\'];
eegsavepath = [codepath,'\datafiles\Subjects_Exp4\EEG\'];
load([datapath,'source_final_label.mat']);% load the label name
%% constrcut head model
headmodel = ft_read_headmodel('standard_bem.mat');
mri = ft_read_mri('standard_mri.mat');
load('standard_sourcemodel3d10mm.mat')

% load our EEG template
load([datapath,'elec_ant64.mat'])
load([datapath,'ant64head.mat'])


% prepare source model (in mm)
cfg           = [];
cfg.warpmni   = 'yes';
cfg.template  = sourcemodel;
cfg.nonlinear = 'yes';
cfg.mri       = mri;
cfg.unit      ='mm';
grid          = ft_prepare_sourcemodel(cfg);

%
cfg             = [];
cfg.resolution  = 10;
cfg.dim         = sourcemodel.dim;
mri = ft_volumereslice(cfg, mri);


% compute leadfield model
cfg                 = [];
cfg.channel         = final_label;
cfg.elec            = elec_ant64;
cfg.headmodel             = headmodel;
cfg.sourcemodel = grid;
% use a 3-D grid with a 1 cm resolution
cfg.resolution       = 10;
cfg.sourcemodel.unit = 'mm';
[leadfield] = ft_prepare_leadfield(cfg);


%% Source analysis
foi_freq = [10,20];
taps = 5;
bl = [-1 -0.5];% baseline
for s = 1:34
    disp(s)
    %% load eeg
    savename = [eegsavepath,sprintf('Sub_eegICA_PV_%d_%s.mat',s,'Exp4')];
    load(savename);
    
    load([sacsavepath,sprintf('TrialEyeInfo_%d_%s.mat',s,'Exp4')])
    
    % filter out rejected trials
    cfg= [];
    cfg.trials = ismember(eeg_final_finetuned.trialinfo(:,1),Eye_trlinfo(Eye_trlinfo(:,9)==1,1));
    eeg_all  = ft_selectdata(cfg, eeg_final_finetuned);
      
    eeg_all.trialinfo(:,2:12) = Eye_trlinfo(ismember(Eye_trlinfo(:,1),eeg_all.trialinfo(:,1)),2:12);
        
    %% 
    
    % 
    chans{1} = 'all'; 
    chans{2} = '-AFz'; 
    
    cfg= [];
    cfg.channel = chans;
    eeg_all  = ft_selectdata(cfg, eeg_all);
    
    %% Rem vs Forg
    % pre & post data 
    % rem
    cfg = [];
    cfg.latency = bl;
    cfg.trials = find(eeg_all.trialinfo(:,2)==1 );
    dataPre_rem     = ft_selectdata(cfg, eeg_all);
    
    cfg.latency = [1 4];
    cfg.trials = find(eeg_all.trialinfo(:,2)==1 );
    dataPost_rem    = ft_selectdata(cfg, eeg_all);
    
    % forg
    cfg = [];
    cfg.latency = bl;
    cfg.trials = find(eeg_all.trialinfo(:,2)==0);
    dataPre_for     = ft_selectdata(cfg, eeg_all);
    
    cfg.latency = [1 4];
    cfg.trials = find(eeg_all.trialinfo(:,2)==0);
    dataPost_for    = ft_selectdata(cfg, eeg_all);
    
    n_trl_RemFor(s,1) = length(find(eeg_all.trialinfo(:,2)==1 ));
    n_trl_RemFor(s,2) = length(find(eeg_all.trialinfo(:,2)==0));
    
    %% Exploration index across Viewing conditions
    cri1 = nanmedian(eeg_all.trialinfo(eeg_all.trialinfo(:,4)==1,8));
    cri2 = nanmedian(eeg_all.trialinfo(eeg_all.trialinfo(:,4)==2,8));
    cri3 = nanmedian(eeg_all.trialinfo(eeg_all.trialinfo(:,4)==3,8));
    
    % Fixation
    cfg= [];
    cfg.latency = bl;
    cfg.trials = (eeg_all.trialinfo(:,8)>=cri1 & eeg_all.trialinfo(:,4)==1) | (eeg_all.trialinfo(:,8)>=cri2 & eeg_all.trialinfo(:,4)==2) | (eeg_all.trialinfo(:,8)>=cri3 & eeg_all.trialinfo(:,4)==3);
    dataPre_Fixation_EI  = ft_selectdata(cfg, eeg_all);
    
   
    cfg.latency = [1 4];
    cfg.trials = (eeg_all.trialinfo(:,8)>=cri1 & eeg_all.trialinfo(:,4)==1) | (eeg_all.trialinfo(:,8)>=cri2 & eeg_all.trialinfo(:,4)==2) | (eeg_all.trialinfo(:,8)>=cri3 & eeg_all.trialinfo(:,4)==3);
    dataPost_Fixation_EI = ft_selectdata(cfg,eeg_all);
    
    
    %  Explr
    cfg= [];
    cfg.latency = bl;
    cfg.trials = (eeg_all.trialinfo(:,8)<=cri1 & eeg_all.trialinfo(:,4)==1) | (eeg_all.trialinfo(:,8)<=cri2 & eeg_all.trialinfo(:,4)==2) | (eeg_all.trialinfo(:,8)<=cri3 & eeg_all.trialinfo(:,4)==3);
    dataPre_Exploration_EI  = ft_selectdata(cfg, eeg_all);
    
   
    cfg.latency = [1 4];
    cfg.trials = (eeg_all.trialinfo(:,8)<=cri1 & eeg_all.trialinfo(:,4)==1) | (eeg_all.trialinfo(:,8)<=cri2 & eeg_all.trialinfo(:,4)==2) | (eeg_all.trialinfo(:,8)<=cri3 & eeg_all.trialinfo(:,4)==3);
    dataPost_Exploration_EI = ft_selectdata(cfg,eeg_all);
    
    n_trl_FixVsExplr_EI(s,1) = sum((eeg_all.trialinfo(:,8)>=cri1 & eeg_all.trialinfo(:,4)==1) | (eeg_all.trialinfo(:,8)>=cri2 & eeg_all.trialinfo(:,4)==2) | (eeg_all.trialinfo(:,8)>=cri3 & eeg_all.trialinfo(:,4)==3));
    n_trl_FixVsExplr_EI(s,2) = sum((eeg_all.trialinfo(:,8)<=cri1 & eeg_all.trialinfo(:,4)==1) | (eeg_all.trialinfo(:,8)<=cri2 & eeg_all.trialinfo(:,4)==2) | (eeg_all.trialinfo(:,8)<=cri3 & eeg_all.trialinfo(:,4)==3));
    
    
    clearvars cri*
    
    %% cross-spectral density matrix
    % rem vs forg
    cfg = [];
    cfg.method    = 'mtmfft';
    cfg.output    = 'powandcsd';
    cfg.tapsmofrq = taps;
    cfg.foilim    = foi_freq;
    freqPre_rem = ft_freqanalysis(cfg, dataPre_rem);
    freqPost_rem = ft_freqanalysis(cfg, dataPost_rem);
    freqPre_for = ft_freqanalysis(cfg, dataPre_for);
    freqPost_for = ft_freqanalysis(cfg, dataPost_for);   
    
    % More vs Less Exploration
    cfg = [];
    cfg.method    = 'mtmfft';
    cfg.output    = 'powandcsd';
    cfg.tapsmofrq = taps;
    cfg.foilim    = foi_freq;
    freqPre_Exploration_EI = ft_freqanalysis(cfg, dataPre_Exploration_EI);
    freqPost_Exploration_EI = ft_freqanalysis(cfg, dataPost_Exploration_EI);
    freqPre_Fixation_EI = ft_freqanalysis(cfg, dataPre_Fixation_EI);
    freqPost_Fixation_EI = ft_freqanalysis(cfg, dataPost_Fixation_EI);  
    
  
    %% perform source analysis
    % compute the frequency domain CSD on Pre and Post
    dataAll_rem = ft_appenddata([], dataPre_rem, dataPost_rem);
    dataAll_for = ft_appenddata([], dataPre_for, dataPost_for);
    
    dataAll_Exploration_EI = ft_appenddata([], dataPre_Exploration_EI, dataPost_Exploration_EI);
    dataAll_Fixation_EI = ft_appenddata([], dataPre_Fixation_EI, dataPost_Fixation_EI);

    
    cfg = [];
    cfg.method    = 'mtmfft';
    cfg.output    = 'powandcsd';
    cfg.tapsmofrq = taps;
    cfg.foilim    = foi_freq;
    freqAll_rem = ft_freqanalysis(cfg, dataAll_rem);
    freqAll_for = ft_freqanalysis(cfg, dataAll_for);
    freqAll_Exploration_EI = ft_freqanalysis(cfg, dataAll_Exploration_EI);
    freqAll_Fixation_EI = ft_freqanalysis(cfg, dataAll_Fixation_EI);

    % compute inverse filter based on both conditions
    cfg              = [];
    cfg.method       = 'dics';
    cfg.frequency    = foi_freq;
    cfg.sourcemodel  = leadfield;
    cfg.headmodel    = headmodel;
    cfg.dics.projectnoise = 'yes';
    cfg.dics.lambda       = '5%';
    cfg.dics.keepfilter   = 'yes';
    cfg.dics.realfilter   = 'yes';
    sourceAll_rem = ft_sourceanalysis(cfg, freqAll_rem);
    sourceAll_for = ft_sourceanalysis(cfg, freqAll_for);
    sourceAll_Exploration_EI = ft_sourceanalysis(cfg, freqAll_Exploration_EI);
    sourceAll_Fixation_EI = ft_sourceanalysis(cfg, freqAll_Fixation_EI);
    
    % apply to each condition  
    cfg.sourcemodel.filter = sourceAll_rem.avg.filter;
    sourcePre_rem  = ft_sourceanalysis(cfg, freqPre_rem );
    sourcePost_rem = ft_sourceanalysis(cfg, freqPost_rem);
    
    cfg.sourcemodel.filter = sourceAll_for.avg.filter;
    sourcePre_for  = ft_sourceanalysis(cfg, freqPre_for );
    sourcePost_for = ft_sourceanalysis(cfg, freqPost_for);
    
    cfg.sourcemodel.filter = sourceAll_Exploration_EI.avg.filter;
    sourcePre_Exploration_EI  = ft_sourceanalysis(cfg, freqPre_Exploration_EI );
    sourcePost_Exploration_EI = ft_sourceanalysis(cfg, freqPost_Exploration_EI);
    
    cfg.sourcemodel.filter = sourceAll_Fixation_EI.avg.filter;
    sourcePre_Fixation_EI  = ft_sourceanalysis(cfg, freqPre_Fixation_EI );
    sourcePost_Fixation_EI = ft_sourceanalysis(cfg, freqPost_Fixation_EI);
    
    
    
    % contrast of (post-pre)/pre
    sourceDiff_rem = sourcePost_rem;
    sourceDiff_rem.avg.pow = ( sourcePost_rem.avg.pow -sourcePre_rem.avg.pow) ./ sourcePre_rem.avg.pow;
    
    sourceDiff_for = sourcePost_for;
    sourceDiff_for.avg.pow = ( sourcePost_for.avg.pow -sourcePre_for.avg.pow) ./ sourcePre_for.avg.pow;
    
    sourceDiff_Exploration_EI = sourcePost_Exploration_EI;
    sourceDiff_Exploration_EI.avg.pow = ( sourcePost_Exploration_EI.avg.pow -sourcePre_Exploration_EI.avg.pow) ./ sourcePre_Exploration_EI.avg.pow;
    
    sourceDiff_Fixation_EI = sourcePost_Fixation_EI;
    sourceDiff_Fixation_EI.avg.pow = ( sourcePost_Fixation_EI.avg.pow -sourcePre_Fixation_EI.avg.pow) ./ sourcePre_Fixation_EI.avg.pow;
    
    
    % interpolate the source to the MRI:
    cfg            = [];
    cfg.downsample = 1;
    cfg.parameter  = 'pow';
    s_DiffInt_rem{s}  = ft_sourceinterpolate(cfg, sourceDiff_rem , mri);
    s_DiffInt_for{s}  = ft_sourceinterpolate(cfg, sourceDiff_for , mri);
    
    s_DiffInt_Exploration_EI{s}  = ft_sourceinterpolate(cfg, sourceDiff_Exploration_EI , mri);
    s_DiffInt_Fixation_EI{s}  = ft_sourceinterpolate(cfg, sourceDiff_Fixation_EI , mri);

    
    clearvars source*
    
end

%% stats
subset = 1:34;

cfg=[];
cfg.dim         = s_DiffInt_rem{1}.dim;
cfg.method      = 'montecarlo';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.parameter   = 'pow';
cfg.correctm    = 'cluster';
cfg.numrandomization = 1000; 
cfg.alpha       = 0.025; 
cfg.tail        = 0;

nsubj=length(subset);
cfg.design(1,:) = [1:nsubj 1:nsubj];
cfg.design(2,:) = [ones(1,nsubj)*1 ones(1,nsubj)*2];
cfg.uvar        = 1; 
cfg.ivar        = 2; 

stat_RF    = ft_sourcestatistics(cfg, s_DiffInt_rem{subset}, s_DiffInt_for{subset});
stat_Explr   = ft_sourcestatistics(cfg, s_DiffInt_Exploration_EI{subset}, s_DiffInt_Fixation_EI{subset});

%% save output
stat_RF = rmfield(stat_RF,'cfg');% remove cfg to save space
stat_Explr = rmfield(stat_Explr,'cfg');
save([datafigspath,'Fig4D_top.mat'],'stat_RF')
save([datafigspath,'Fig4E_top.mat'],'stat_Explr')
%% Figure
codepath = '';

datapath = [codepath,'\datafiles\'];
datafigspath =  [codepath,'\data4figs\'];
figsavepath = [codepath,'\figures\'];

% load behavioral data, gaze data
load( [datafigspath,'Fig4D_top.mat'])
load( [datafigspath,'Fig4E_top.mat'])

plot_save = 1; % to save or not figure

stat_plot = stat_RF;
stat_plot.posclusters = stat_plot.negclusters(end);% just to copy the structure so fieldtrip doesn't throw an error for plotting
stat_plot.maskforplot = stat_plot.stat.*(stat_plot.negclusterslabelmat==1);% Note: for visualisation, mask out other clusters


fig_savename = 'Fig4D_top';
for fig=11
sizefig1 = [0.2,0.2,0.4,0.4];
fontaxis = 24;line_width = 2;

f=figure('Name',int2str(11),'units','normalized','outerposition',sizefig1);clf;
cfg = [];
cfg.method         = 'surface';
cfg.funparameter   = 'stat';
cfg.maskparameter  = 'maskforplot';
cfg.maskstyle     = 'opacity';
cfg.funcolorlim    = [-5,5];
cfg.opacitylim     = 'auto';
cfg.opacitymap     = 'auto';
cfg.projmethod     = 'nearest';
cfg.surffile       = 'surface_pial_both.mat';
cfg.surfdownsample = 10;
ft_sourceplot(cfg, stat_plot);
view ([0 40])
camlight(0,40)
colormap(flipud(cbrewer2('RdBu','div')));
c = colorbar;
c.Label.String = 't value';
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',24,'FontWeight','normal', 'LineWidth', 2);


fig_savename1 = [fig_savename ,'_view_1'];

if plot_save
    saveas(gcf,[figsavepath,fig_savename1,'.emf'])
    saveas(gcf,[figsavepath,fig_savename1,'.png'])
end

camlight('headlight')
view ([0 90])

fig_savename1 = [fig_savename ,'_view_2'];

if plot_save
    saveas(gcf,[figsavepath,fig_savename1,'.emf'])
    saveas(gcf,[figsavepath,fig_savename1,'.png'])
end

end


%

stat_plot = stat_Explr;
stat_plot.posclusters = stat_plot.negclusters(end);% just to copy the structure so fieldtrip doesn't throw an error for plotting
stat_plot.maskforplot = stat_plot.stat.*(stat_plot.negclusterslabelmat==1);% Note: for visualisation, mask out other clusters

fig_savename = 'Fig4E_top';
for fig=11
sizefig1 = [0.2,0.2,0.4,0.4];
fontaxis = 24;line_width = 2;

f=figure('Name',int2str(11),'units','normalized','outerposition',sizefig1);clf;
cfg = [];
cfg.method         = 'surface';
cfg.funparameter   = 'stat';
cfg.maskparameter  = 'maskforplot';
cfg.maskstyle     = 'opacity';
cfg.funcolorlim    = [-5,5];
cfg.opacitylim     = 'auto';
cfg.opacitymap     = 'auto';
cfg.projmethod     = 'nearest';
cfg.surffile       = 'surface_pial_both.mat';
cfg.surfdownsample = 10;
ft_sourceplot(cfg, stat_plot);
view ([0 40])
camlight(0,40)
colormap(flipud(cbrewer2('RdBu','div')));
c = colorbar;
c.Label.String = 't value';
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',24,'FontWeight','normal', 'LineWidth', 2);


fig_savename1 = [fig_savename ,'_view_1'];

if plot_save
    saveas(gcf,[figsavepath,fig_savename1,'.emf'])
    saveas(gcf,[figsavepath,fig_savename1,'.png'])
end

camlight('headlight')
view ([0 90])

fig_savename1 = [fig_savename ,'_view_2'];

if plot_save
    saveas(gcf,[figsavepath,fig_savename1,'.emf'])
    saveas(gcf,[figsavepath,fig_savename1,'.png'])
end

end
