%% Gaze density analysis for Exp2
% run this after extracting gaze density array (GazeDensity_Exp2)
codepath = '';

datapath = [codepath,'\datafiles\'];
datafigspath =  [codepath,'\data4figs\'];
figsavepath = [codepath,'\figures\'];
subinfopath = [codepath,'\datafiles\Subjects_Exp2\ExpInfo\'];
sacsavepath = [codepath,'\datafiles\Subjects_Exp2\saccades\'];
% load behavioral data, gaze data (P01 & 02 with Tobii and P03 & 04 with Eyelink)
load( [datapath,'Beh_Exp2.mat'])
load( [datapath,'PVstim_Gaze_Exp2_tb.mat'])
gz_tb = gaze2D; clearvars gaze
load( [datapath,'PVstim_Gaze_Exp2_el.mat'])
gaze2D_PV = [gz_tb,gaze2D]; clearvars gaze gz_tb gaze2D

%% Separate condition
for s = 1:length(gaze2D_PV)
    disp(s)
    
    load([sacsavepath,sprintf('TrialEyeInfo_%d_%s.mat',s,'Exp2')])
    
    g_data = cat(3,gaze2D_PV{s}{Eye_trlinfo(:,1)});
    
    %% sep by Rem vs Forg with balancing N of trials 
    [Gaze_Rem_balance(:,:,s),Gaze_For_balance(:,:,s),n_trl] = Gaze_balance_trial(g_data(:,:,Eye_trlinfo(:,2)==1 ),g_data(:,:,Eye_trlinfo(:,2)==0 ),100);
    n_trl_RF_balance(s,1:2) = n_trl;
    
    %% sep by N of saccades 
    cri = nanmedian(Eye_trlinfo(:,5));

    Gaze_MoreSac(:,:,s) = nanmean(g_data(:,:,Eye_trlinfo(:,5)>=cri ),3);
    Gaze_FewerSac(:,:,s) = nanmean(g_data(:,:,Eye_trlinfo(:,5)<=cri ),3);
    
     %% sep by Explore Idx
    cri = nanmedian(Eye_trlinfo(:,7));
    
    Gaze_MoreExplr(:,:,s) = nanmean(g_data(:,:,Eye_trlinfo(:,7)<=cri ),3);
    Gaze_LessExplr(:,:,s) = nanmean(g_data(:,:,Eye_trlinfo(:,7)>=cri ),3);
    

    clearvars g_data cri
    
end


%% Stat for plotting
% Rem vs Forg
subset = 1:4; % all subjects
for hitvsmis = 1
ave_2D_C1 = Gaze_Rem_balance;
ave_2D_C2 = Gaze_For_balance;
viewdva_x =28;
viewdva_y =21;
%%
Dat_C1 = [];
Dat_C1.dat = permute(ave_2D_C1,[4,1,2,3]);
Dat_C1.freq = linspace(-viewdva_y/2,viewdva_y/2,size(ave_2D_C1,1))';
Dat_C1.time = linspace(-viewdva_x/2,viewdva_x/2,size(ave_2D_C1,2))';
Dat_C1.dimord = 'chan_freq_time_subj';
Dat_C1.label = {'avg'};

Dat_C2 = Dat_C1;
Dat_C2.dat = permute(ave_2D_C2,[4,1,2,3]);


design      = zeros(2,length(subset)*2);
design(1,:) = repmat(1:length(subset),[1 2]);
design(2,:) = [ones(1,length(subset)),ones(1,length(subset))+1];
    
cfg = [];
cfg.spmversion = 'spm12';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.ivar                = 2;
cfg.uvar                = 1;
cfg.parameter           = 'dat';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 1000;

cfg.design           = design;


[stat_CB]  = ft_freqstatistics(cfg, Dat_C1,Dat_C2);

end

% More vs fewer saccades
for MvLsac = 1
ave_2D_C1 = Gaze_MoreSac;
ave_2D_C2 = Gaze_FewerSac;
viewdva_x =28;
viewdva_y =21;
%%
Dat_C1 = [];
Dat_C1.dat = permute(ave_2D_C1,[4,1,2,3]);
Dat_C1.freq = linspace(-viewdva_y/2,viewdva_y/2,size(ave_2D_C1,1))';
Dat_C1.time = linspace(-viewdva_x/2,viewdva_x/2,size(ave_2D_C1,2))';
Dat_C1.dimord = 'chan_freq_time_subj';
Dat_C1.label = {'avg'};

Dat_C2 = Dat_C1;
Dat_C2.dat = permute(ave_2D_C2,[4,1,2,3]);


design      = zeros(2,length(subset)*2);
design(1,:) = repmat(1:length(subset),[1 2]);
design(2,:) = [ones(1,length(subset)),ones(1,length(subset))+1];
    
cfg = [];
cfg.spmversion = 'spm12';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.ivar                = 2;
cfg.uvar                = 1;
cfg.parameter           = 'dat';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 1000;

cfg.design           = design;


[stat_Sac]  = ft_freqstatistics(cfg, Dat_C1,Dat_C2);

end

% More vs Less Exploration
for MvLExplr = 1
ave_2D_C1 = Gaze_MoreExplr;
ave_2D_C2 = Gaze_LessExplr;
viewdva_x =28;
viewdva_y =21;
%%
Dat_C1 = [];
Dat_C1.dat = permute(ave_2D_C1,[4,1,2,3]);
Dat_C1.freq = linspace(-viewdva_y/2,viewdva_y/2,size(ave_2D_C1,1))';
Dat_C1.time = linspace(-viewdva_x/2,viewdva_x/2,size(ave_2D_C1,2))';
Dat_C1.dimord = 'chan_freq_time_subj';
Dat_C1.label = {'avg'};

Dat_C2 = Dat_C1;
Dat_C2.dat = permute(ave_2D_C2,[4,1,2,3]);


design      = zeros(2,length(subset)*2);
design(1,:) = repmat(1:length(subset),[1 2]);
design(2,:) = [ones(1,length(subset)),ones(1,length(subset))+1];
    
cfg = [];
cfg.spmversion = 'spm12';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.ivar                = 2;
cfg.uvar                = 1;
cfg.parameter           = 'dat';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 1000;

cfg.design           = design;


[stat_EI]  = ft_freqstatistics(cfg, Dat_C1,Dat_C2);

end

%%
save([datafigspath,'Fig2C_right.mat'],'stat_CB')
save([datafigspath,'Fig2D_right.mat'],'stat_Sac')
save([datafigspath,'SuppleFig4_right21.mat'],'stat_EI','Gaze_MoreExplr','Gaze_LessExplr')
%%
%% Figure
% load data for plotting
codepath = '';

datafigspath =  [codepath,'\data4figs\'];
figsavepath = [codepath,'\figures\'];
% load behavioral data, gaze data
load( [datafigspath,'Fig2C_right.mat'])
load( [datafigspath,'Fig2D_right.mat'])
plot_save = 1; % to save or not figure

sizefig1 = [0,0.2,0.3,0.45];
fontaxis = 24;line_width = 2;
alpha = 0.05;n_bin_x = 80;n_bin_y = 60;
viewdva_x =28;
viewdva_y =21;

fig_savename1 = 'Fig2C_right';
for hitvsmis=1
f=figure('Name',int2str(11),'units','normalized','outerposition',sizefig1);clf;
imagesc(squeeze(stat_CB.stat));
xticks(linspace(0,n_bin_x,5))
xticklabels(num2cell(abs(linspace(0,n_bin_x,5)*viewdva_x/n_bin_x-viewdva_x/2)))
yticks(linspace(0,n_bin_y,5))
yticklabels(num2cell(abs(linspace(0,n_bin_y,5)*viewdva_y/n_bin_y-viewdva_y/2)))
caxis([-3,3])
xlabel('Dist to center [dva]')
%ylabel('Dist to center [dva]')
colormap(flipud(cbrewer2('RdBu','div')));
c = colorbar;
c.Label.String = 't value';
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',fontaxis,'FontWeight','normal', 'LineWidth', line_width);

    if plot_save
        saveas(gcf,[figsavepath,fig_savename1,'.emf'])
        saveas(gcf,[figsavepath,fig_savename1,'.png'])
    end
    
end 


%
fig_savename1 = 'Fig2D_right';

for MvLsac=1
f=figure('Name',int2str(11),'units','normalized','outerposition',sizefig1);clf;
imagesc(squeeze(stat_Sac.stat));
xticks(linspace(0,n_bin_x,5))
xticklabels(num2cell(abs(linspace(0,n_bin_x,5)*viewdva_x/n_bin_x-viewdva_x/2)))
yticks(linspace(0,n_bin_y,5))
yticklabels(num2cell(abs(linspace(0,n_bin_y,5)*viewdva_y/n_bin_y-viewdva_y/2)))
caxis([-3,3])
xlabel('Dist to center [dva]')
%ylabel('Dist to center [dva]')
colormap(flipud(cbrewer2('RdBu','div')));
c = colorbar;
c.Label.String = 't value';
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',fontaxis,'FontWeight','normal', 'LineWidth', line_width);

    if plot_save
        saveas(gcf,[figsavepath,fig_savename1,'.emf'])
        saveas(gcf,[figsavepath,fig_savename1,'.png'])
    end
    
end 


