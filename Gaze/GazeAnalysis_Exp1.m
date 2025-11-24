%% Gaze density analysis for Exp1
% run this after extracting gaze density array (GazeDensity_Exp1)
codepath = '';

datapath = [codepath,'\datafiles\'];
datafigspath =  [codepath,'\data4figs\'];
figsavepath = [codepath,'\figures\'];
subinfopath = [codepath,'\datafiles\Subjects_Exp1\ExpInfo\'];
sacsavepath = [codepath,'\datafiles\Subjects_Exp1\saccades\'];
% load behavioral data, gaze data
load( [datapath,'Beh_Exp1.mat'])
load( [datapath,'PVstim_Gaze_Exp1.mat'])
%% Separate conditions
for s = 1:length(gaze2D)
    disp(s)
    %% load trial info to split conditions
    
    load([sacsavepath,sprintf('TrialEyeInfo_%d_%s.mat',s,'Exp1')])
    
    g_data = cat(3,gaze2D{s}{Eye_trlinfo(:,1)});

    %% sep by Rem vs Forg with balancing N of trials 
    [Gaze_Rem_balance(:,:,s),Gaze_For_balance(:,:,s),n_trl] = Gaze_balance_trial(g_data(:,:,Eye_trlinfo(:,2)==1 & Eye_trlinfo(:,8)==1),g_data(:,:,Eye_trlinfo(:,2)==0 & Eye_trlinfo(:,8)==1),100);
    n_trl_RF_balance(s,1:2) = n_trl;
    
    %% sep by N of saccades 
    cri = nanmedian(Eye_trlinfo(Eye_trlinfo(:,8)==1,5));

    Gaze_MoreSac(:,:,s) = nanmean(g_data(:,:,Eye_trlinfo(:,5)>=cri & Eye_trlinfo(:,8)==1),3);
    Gaze_FewerSac(:,:,s) = nanmean(g_data(:,:,Eye_trlinfo(:,5)<=cri & Eye_trlinfo(:,8)==1),3);
 
    clearvars  cri
     %% Median split Exploration index 
    cri = nanmedian(Eye_trlinfo(Eye_trlinfo(:,8)==1,7));
    
    Gaze_MoreExplr(:,:,s) = nanmean(g_data(:,:,Eye_trlinfo(:,7)<=cri & Eye_trlinfo(:,8)==1),3);
    Gaze_LessExplr(:,:,s) = nanmean(g_data(:,:,Eye_trlinfo(:,7)>=cri & Eye_trlinfo(:,8)==1),3);
    
    clearvars g_data cri
end

%%
% Rem vs Forg
subset = 1:20; % all subjects
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


%% Gaze with behavior
% gaze within cluster 

Cluster_RF_balance=zeros(size(Gaze_Rem_balance,[1,2]));
for ic=1:length(stat_CB.posclusters)
    if stat_CB.posclusters(ic).prob <0.025
        Cluster_RF_balance(squeeze(stat_CB.posclusterslabelmat==ic)) = 1;
    end
end

% get gaze within the cluster
for s= 1:length(gaze2D)
    
    load([sacsavepath,sprintf('TrialEyeInfo_%d_%s.mat',s,'Exp1')])
    
    g_balance = cellfun(@(x) nansum(x(Cluster_RF_balance>0)), gaze2D{s}, 'UniformOutput', false);
      
    Eye_trlinfo = Eye_trlinfo(Eye_trlinfo(:,8)==1,:);
    ave_cl_gaze_RF_balance{s}(:,1:2) = [cat(1,g_balance{Eye_trlinfo(:,1)}),Eye_trlinfo(:,2)];
    
    %% Mem  split by N Saccade
    cri = nanmedian(Eye_trlinfo(:,5));
    
    MemMore_Sac(s,1) = mean(Eye_trlinfo(Eye_trlinfo(:,5)>=cri,2))*100;
    MemMore_Sac(s,2) = (1-mean(Eye_trlinfo(Eye_trlinfo(:,5)>=cri,2)))*100;
    
    MemFewer_Sac(s,1) = mean(Eye_trlinfo(Eye_trlinfo(:,5)<=cri,2))*100;
    MemFewer_Sac(s,2) = (1-mean(Eye_trlinfo(Eye_trlinfo(:,5)<=cri,2)))*100;

end


for s = 1:length(ave_cl_gaze_RF_balance)
    cri = nanmedian(ave_cl_gaze_RF_balance{s}(:,1));
   
    
    % with balancing
    MemHigh_gaze_cl_balance(s,1) = mean(ave_cl_gaze_RF_balance{s}(ave_cl_gaze_RF_balance{s}(:,1)>=cri,2))*100;
    MemHigh_gaze_cl_balance(s,2) = (1-mean(ave_cl_gaze_RF_balance{s}(ave_cl_gaze_RF_balance{s}(:,1)>=cri,2)))*100;
    
    MemLow_gaze_cl_balance(s,1) = mean(ave_cl_gaze_RF_balance{s}(ave_cl_gaze_RF_balance{s}(:,1)<=cri,2))*100;
    MemLow_gaze_cl_balance(s,2) = (1-mean(ave_cl_gaze_RF_balance{s}(ave_cl_gaze_RF_balance{s}(:,1)<=cri,2)))*100;
    
end


%% save output for plotting
save([datafigspath,'Fig1B_left.mat'],'stat_CB')
save([datafigspath,'Fig1D_left.mat'],'stat_Sac')
save([datafigspath,'Fig1D_right.mat'],'MemMore_Sac','MemFewer_Sac','MemHigh_gaze_cl_balance','MemLow_gaze_cl_balance')
save([datafigspath,'SuppleFig4_right11.mat'],'stat_EI','Gaze_MoreExplr','Gaze_LessExplr')
%% Figure
% load data for plotting
codepath = '';

datafigspath =  [codepath,'\data4figs\'];
figsavepath = [codepath,'\figures\'];
% load behavioral data, gaze data
load( [datafigspath,'Fig1B_left.mat'])
load( [datafigspath,'Fig1D_left.mat'])
load( [datafigspath,'Fig1D_right.mat'])
plot_save = 1; % to save or not figure

sizefig1 = [0,0.2,0.3,0.45];
fontaxis = 24;line_width = 2;
alpha = 0.05;n_bin_x = 80;n_bin_y = 60;
viewdva_x =28;
viewdva_y =21;

fig_savename1 = 'Fig1B_left';
for hitvsmis=1
f=figure('Name',int2str(11),'units','normalized','outerposition',sizefig1);clf;
Cluster=zeros(n_bin_y,n_bin_x);
Cluster(squeeze(stat_CB.mask)) = 1;
[B,L,N,A] = bwboundaries(Cluster); 
imagesc(squeeze(stat_CB.stat));hold on
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


%
fig_savename1 = 'Fig1D_left';

for MvLsac=1
f=figure('Name',int2str(11),'units','normalized','outerposition',sizefig1);clf;
Cluster=zeros(n_bin_y,n_bin_x);
Cluster(squeeze(stat_Sac.mask)) = 1;
[B,L,N,A] = bwboundaries(Cluster); 
imagesc(squeeze(stat_Sac.stat));hold on
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


[cb2] = cbrewer2('RdBu','div',12);
[cb] = cbrewer2('BrBG','div',12);

sizefig1 = [0.2,0.2,0.75,0.6];
fontaxis = 16;line_width = 2;
fig_savename1 = 'Fig1D_right1';
fig_savename2 = 'Fig1D_right2';

for fig=50
    % example 1
    f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig1);clf;
    subplot(1, 2 ,1)
    h1 = raincloud_plot(MemMore_Sac(:,1), 'box_on', 1, 'color', cb(3,:), 'alpha', 0.6,...
        'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
        'box_col_match', 0);
    h2 = raincloud_plot(MemFewer_Sac(:,1), 'box_on', 1, 'color', cb(9,:), 'alpha', 0.6,...
        'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
    legend([h1{1} h2{1}], {'More Sac', 'Fewer Sac'})
    title('Rem')
    xlabel('Percentage of trials [%]')

    
    set(gca,'XLim', [0 100], 'YLim', [-.02 .05],'ytick',[],'FontName','Arial','FontSize',fontaxis,'FontWeight','normal', 'LineWidth', line_width);
    box off
    
    % example 2
    subplot(1, 2, 2)
    h1 = raincloud_plot(MemMore_Sac(:,2), 'box_on', 1, 'color', cb(3,:), 'alpha', 0.6,...
        'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
        'box_col_match', 0);
    h2 = raincloud_plot(MemFewer_Sac(:,2), 'box_on', 1, 'color', cb(9,:), 'alpha', 0.6,...
        'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
    legend([h1{1} h2{1}], {'More Sac', 'Fewer Sac'})
    title('Forg')
    xlabel('Percentage of trials [%]')
    set(gca,'XLim', [0 100], 'YLim', [-.02 .05],'ytick',[],'FontName','Arial','FontSize',fontaxis,'FontWeight','normal', 'LineWidth', line_width);
    box off

    if plot_save
        saveas(gcf,[figsavepath,fig_savename1,'.emf'])
        saveas(gcf,[figsavepath,fig_savename1,'.png'])
    end
    
end

for fig=52
    f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig1);clf;
    subplot(1, 2 ,1)
    h1 = raincloud_plot(MemMore_Sac(:,1), 'box_on', 1, 'color', cb2(3,:), 'alpha', 0.6,...
        'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
        'box_col_match', 0);
    h2 = raincloud_plot(MemMore_Sac(:,2), 'box_on', 1, 'color', cb2(9,:), 'alpha', 0.6,...
        'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
    legend([h1{1} h2{1}], {'Rem', 'Forg'})
    title('More Sac')
    xlabel('Percentage of trials [%]')
    set(gca,'XLim', [0 100], 'YLim', [-.02 .05],'ytick',[],'FontName','Arial','FontSize',fontaxis,'FontWeight','normal', 'LineWidth', line_width);
    box off
    
    subplot(1, 2, 2)
    h1 = raincloud_plot(MemFewer_Sac(:,1), 'box_on', 1, 'color', cb2(3,:), 'alpha', 0.6,...
        'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
        'box_col_match', 0);
    h2 = raincloud_plot(MemFewer_Sac(:,2), 'box_on', 1, 'color', cb2(9,:), 'alpha', 0.6,...
        'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
    legend([h1{1} h2{1}], {'Rem', 'Forg'})
    title('Fewer Sac')
    xlabel('Percentage of trials [%]')
    set(gca,'XLim', [0 100], 'YLim', [-.02 .05],'ytick',[],'FontName','Arial','FontSize',fontaxis,'FontWeight','normal', 'LineWidth', line_width);
    box off
    
    if plot_save
        saveas(gcf,[figsavepath,fig_savename2,'.emf'])
        saveas(gcf,[figsavepath,fig_savename2,'.png'])
    end

end

