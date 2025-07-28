%% Gaze density analysis for Exp4
% run this after extracting gaze density array (GazeDensity_Exp4)
codepath = '';

datapath = [codepath,'\datafiles\'];
datafigspath =  [codepath,'\data4figs\'];
figsavepath = [codepath,'\figures\'];
subinfopath = [codepath,'\datafiles\Subjects_Exp4\ExpInfo\'];
sacsavepath = [codepath,'\datafiles\Subjects_Exp4\saccades\'];
% load behavioral data, gaze data
load( [datapath,'Beh_Exp4.mat'])
load( [datapath,'PVstim_Gaze_Exp4.mat'])
%% Separate conditions
for s = 1:length(gaze2D)
    disp(s)
    %% load trial info to split conditions
    
    load([sacsavepath,sprintf('TrialEyeInfo_%d_%s.mat',s,'Exp4')])
    
    g_data = cat(3,gaze2D{s}{Eye_trlinfo(:,1)});

    %% sep by Rem vs Forg with balancing N of trials 
    [Gaze_Rem_balance(:,:,s),Gaze_For_balance(:,:,s),n_trl] = Gaze_balance_trial(g_data(:,:,Eye_trlinfo(:,2)==1 & Eye_trlinfo(:,9)==1),g_data(:,:,Eye_trlinfo(:,2)==0 & Eye_trlinfo(:,9)==1),100);
    n_trl_RF_balance(s,1:2) = n_trl;
    %% sep by Viewing condition 
    Gaze_VC1(:,:,s) = nanmean(g_data(:,:,Eye_trlinfo(:,4)==1 & Eye_trlinfo(:,9)==1),3); 
    Gaze_VC2(:,:,s) = nanmean(g_data(:,:,Eye_trlinfo(:,4)==2 & Eye_trlinfo(:,9)==1),3); 
    Gaze_VC3(:,:,s) = nanmean(g_data(:,:,Eye_trlinfo(:,4)==3 & Eye_trlinfo(:,9)==1),3); 
    %% Median split Exploration index 
    cri = nanmedian(Eye_trlinfo(Eye_trlinfo(:,9)==1,8));
    
    Gaze_MoreExplr(:,:,s) = nanmean(g_data(:,:,Eye_trlinfo(:,8)<=cri & Eye_trlinfo(:,9)==1),3);
    Gaze_LessExplr(:,:,s) = nanmean(g_data(:,:,Eye_trlinfo(:,8)>=cri & Eye_trlinfo(:,9)==1),3);
    
end
%% Stat analysis
subset = 1:34; % all subjects

% Viewing condtion 
for VC = 1
t_name1 = '5 vs 1 DVA';
t_name2 = '21 vs 5 DVA';
t_name3 = '21 vs 1 DVA ';
ave_2D_C1 = Gaze_VC1;
ave_2D_C2 = Gaze_VC2;
ave_2D_C3 = Gaze_VC3;

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

Dat_C3 = Dat_C1;
Dat_C3.dat = permute(ave_2D_C3,[4,1,2,3]);


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


[stat_C21]  = ft_freqstatistics(cfg, Dat_C2,Dat_C1);
[stat_C32]  = ft_freqstatistics(cfg, Dat_C3,Dat_C2);
[stat_C31]  = ft_freqstatistics(cfg, Dat_C3,Dat_C1);
end

% Rem vs Forg
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
%% save output for plotting
save([datafigspath,'Fig4B_left.mat'],'stat_C21','stat_C32','stat_C31')
save([datafigspath,'Fig4B_right.mat'],'stat_CB')

save([datafigspath,'SuppleFig2_right41.mat'],'stat_EI','Gaze_MoreExplr','Gaze_LessExplr')

%% Figure
% load data for plotting
codepath = '';

datafigspath =  [codepath,'\data4figs\'];
figsavepath = [codepath,'\figures\'];
% load behavioral data, gaze data
load( [datafigspath,'Fig4B_left.mat'])
load( [datafigspath,'Fig4B_right.mat'])

plot_save = 1; % to save or not figure

sizefig1 = [0,0.2,1,0.45];
fontaxis = 24;line_width = 2;
alpha = 0.05;n_bin_x = 80;n_bin_y = 60;
viewdva_x =28;
viewdva_y =21;

fig_savename1 = 'Fig4B_left';

for VC=1
f=figure('Name',int2str(11),'units','normalized','outerposition',sizefig1);clf;
subplot(1,3,1);
Cluster=zeros(n_bin_y,n_bin_x);
Cluster(squeeze(stat_C21.mask)) = 1;
[B,L,N,A] = bwboundaries(Cluster); 
imAlpha = ones(size(squeeze(stat_C21.stat)));
imAlpha(isnan(squeeze(stat_C21.stat))) = 0;
imagesc(squeeze(stat_C21.stat),'AlphaData',imAlpha);
hold on
for k=1:length(B)
     boundary = B{k};
    plot(boundary(:,2), boundary(:,1), 'LineWidth',3,'color','k');
end
hold off
xticks(linspace(0,n_bin_x,5))
xticklabels(num2cell(abs(linspace(0,n_bin_x,5)*viewdva_x/n_bin_x-viewdva_x/2)))
yticks(linspace(0,n_bin_y,5))
yticklabels(num2cell(abs(linspace(0,n_bin_y,5)*viewdva_y/n_bin_y-viewdva_y/2)))
caxis([-2.5,2.5])
title('5 vs 1 dva')
xlabel('Dist to center [dva]')

%ylabel('Dist to center [dva]')

colormap(flipud(cbrewer2('RdBu','div')));
c = colorbar;
c.Label.String = 't value';
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',fontaxis,'FontWeight','normal', 'LineWidth', line_width,'color',0.5*[1 1 1]);


subplot(1,3,2);
Cluster=zeros(n_bin_y,n_bin_x);
Cluster(squeeze(stat_C32.mask)) = 1;
[B,L,N,A] = bwboundaries(Cluster); 
imAlpha = ones(size(squeeze(stat_C32.stat)));
imAlpha(isnan(squeeze(stat_C32.stat))) = 0;
imagesc(squeeze(stat_C32.stat),'AlphaData',imAlpha);
hold on
for k=1:length(B)
     boundary = B{k};
    plot(boundary(:,2), boundary(:,1), 'LineWidth',3,'color','k');
end
hold off
xticks(linspace(0,n_bin_x,5))
xticklabels(num2cell(abs(linspace(0,n_bin_x,5)*viewdva_x/n_bin_x-viewdva_x/2)))
yticks(linspace(0,n_bin_y,5))
yticklabels(num2cell(abs(linspace(0,n_bin_y,5)*viewdva_y/n_bin_y-viewdva_y/2)))
caxis([-2.5,2.5])
title('21 vs 5 dva')
xlabel('Dist to center [dva]')
%ylabel('Dist to center [dva]')

colormap(flipud(cbrewer2('RdBu','div')));
c = colorbar;
c.Label.String = 't value';
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',fontaxis,'FontWeight','normal', 'LineWidth', line_width,'color',0.5*[1 1 1]);


subplot(1,3,3);
Cluster=zeros(n_bin_y,n_bin_x);
Cluster(squeeze(stat_C31.mask)) = 1;
[B,L,N,A] = bwboundaries(Cluster); 
imAlpha = ones(size(squeeze(stat_C31.stat)));
imAlpha(isnan(squeeze(stat_C31.stat))) = 0;
imagesc(squeeze(stat_C31.stat),'AlphaData',imAlpha);
hold on
for k=1:length(B)
     boundary = B{k};
    plot(boundary(:,2), boundary(:,1), 'LineWidth',3,'color','k');
end
hold off
xticks(linspace(0,n_bin_x,5))
xticklabels(num2cell(abs(linspace(0,n_bin_x,5)*viewdva_x/n_bin_x-viewdva_x/2)))
yticks(linspace(0,n_bin_y,5))
yticklabels(num2cell(abs(linspace(0,n_bin_y,5)*viewdva_y/n_bin_y-viewdva_y/2)))
caxis([-2.5,2.5])
title('21 vs 1 dva')
xlabel('Dist to center [dva]')
%ylabel('Dist to center [dva]')

hold off

colormap(flipud(cbrewer2('RdBu','div')));
c = colorbar;
c.Label.String = 't value';
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',fontaxis,'FontWeight','normal', 'LineWidth', line_width,'color',0.5*[1 1 1]);

    if plot_save
        saveas(gcf,[figsavepath,fig_savename1,'.emf'])
        saveas(gcf,[figsavepath,fig_savename1,'.png'])
    end
    
end 

sizefig1 = [0,0.2,0.3,0.45];
%
fig_savename1 = 'Fig4B_right';
for hitvsmis=1
f=figure('Name',int2str(11),'units','normalized','outerposition',sizefig1);clf;
Cluster=zeros(n_bin_y,n_bin_x);
Cluster(squeeze(stat_CB.mask)) = 1;
[B,L,N,A] = bwboundaries(Cluster); 
imagesc(squeeze(stat_CB.stat));
imAlpha = ones(size(squeeze(stat_CB.stat)));
imAlpha(isnan(squeeze(stat_CB.stat))) = 0;
imagesc(squeeze(stat_CB.stat),'AlphaData',imAlpha);
hold on
for k=1:length(B)
     boundary = B{k};
    plot(boundary(:,2), boundary(:,1), 'LineWidth',3,'color','k');
end
hold off
xticks(linspace(0,n_bin_x,5))
xticklabels(num2cell(abs(linspace(0,n_bin_x,5)*viewdva_x/n_bin_x-viewdva_x/2)))
yticks(linspace(0,n_bin_y,5))
yticklabels(num2cell(abs(linspace(0,n_bin_y,5)*viewdva_y/n_bin_y-viewdva_y/2)))
caxis([-2.5,2.5])
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

