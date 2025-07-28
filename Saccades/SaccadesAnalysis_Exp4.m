%% Saccade analysis for Experiment 4
% run this after completing saccades detection
codepath = '';

datapath = [codepath,'\datafiles\'];
datafigspath =  [codepath,'\data4figs\'];
figsavepath = [codepath,'\figures\'];
subinfopath = [codepath,'\datafiles\Subjects_Exp4\ExpInfo\'];
sacsavepath = [codepath,'\datafiles\Subjects_Exp4\saccades\'];
%% N of saccades
% general N of saccades
for s = 1:34
    disp(s)
    %% load trial level saccade information
    load([sacsavepath,sprintf('TrialEyeInfo_%d_%s.mat',s,'Exp4')])
    
    %% Sep by Rem vs Forg
    Rem_trl = Eye_trlinfo(Eye_trlinfo(:,9)==1 & Eye_trlinfo(:,2)==1,:);
    For_trl = Eye_trlinfo(Eye_trlinfo(:,9)==1 & Eye_trlinfo(:,2)==0,:);

    n_sac_RF_PV(s,1) =  nanmean(Rem_trl(:,6))./4;
    n_sac_RF_PV(s,2) =  nanmean(For_trl(:,6))./4;
    clearvars Rem_trl For_trl
    
    %% Sep by viewing condition
    VC1_trl = Eye_trlinfo(Eye_trlinfo(:,9)==1 & Eye_trlinfo(:,4)==1,:);
    VC2_trl = Eye_trlinfo(Eye_trlinfo(:,9)==1 & Eye_trlinfo(:,4)==2,:);
    VC3_trl = Eye_trlinfo(Eye_trlinfo(:,9)==1 & Eye_trlinfo(:,4)==3,:);
    
    n_sac_VC_PV(s,1) =  nanmean(VC1_trl(:,6))./4;
    n_sac_VC_PV(s,2) =  nanmean(VC2_trl(:,6))./4;
    n_sac_VC_PV(s,3) =  nanmean(VC3_trl(:,6))./4;
    
    clearvars VC*
    
    %% sep by VC x Mem
    VC1_Rem_trl = Eye_trlinfo(Eye_trlinfo(:,9)==1 & Eye_trlinfo(:,4)==1 & Eye_trlinfo(:,2)==1,:);
    VC2_Rem_trl = Eye_trlinfo(Eye_trlinfo(:,9)==1 & Eye_trlinfo(:,4)==2 & Eye_trlinfo(:,2)==1,:);
    VC3_Rem_trl = Eye_trlinfo(Eye_trlinfo(:,9)==1 & Eye_trlinfo(:,4)==3 & Eye_trlinfo(:,2)==1,:);
    
    VC1_For_trl = Eye_trlinfo(Eye_trlinfo(:,9)==1 & Eye_trlinfo(:,4)==1 & Eye_trlinfo(:,2)==0,:);
    VC2_For_trl = Eye_trlinfo(Eye_trlinfo(:,9)==1 & Eye_trlinfo(:,4)==2 & Eye_trlinfo(:,2)==0,:);
    VC3_For_trl = Eye_trlinfo(Eye_trlinfo(:,9)==1 & Eye_trlinfo(:,4)==3 & Eye_trlinfo(:,2)==0,:);
    
    
    n_sac_VCRF_PV(s,1) =  nanmean(VC1_Rem_trl(:,6))./4;% pic viewing
    n_sac_VCRF_PV(s,2) =  nanmean(VC2_Rem_trl(:,6))./4;
    n_sac_VCRF_PV(s,3) =  nanmean(VC3_Rem_trl(:,6))./4;
    n_sac_VCRF_PV(s,4) =  nanmean(VC1_For_trl(:,6))./4;
    n_sac_VCRF_PV(s,5) =  nanmean(VC2_For_trl(:,6))./4;
    n_sac_VCRF_PV(s,6) =  nanmean(VC3_For_trl(:,6))./4;

   
    clearvars VC*
      
    %% Sey by Exploration index
    cri = nanmedian(Eye_trlinfo(Eye_trlinfo(:,9)==1,8));
    More_trl = Eye_trlinfo(Eye_trlinfo(:,9)==1 & Eye_trlinfo(:,8)<=cri,:);
    Less_trl = Eye_trlinfo(Eye_trlinfo(:,9)==1 & Eye_trlinfo(:,8)>=cri,:);

    
    n_sac_Explr_PV(s,1) =  nanmean(More_trl(:,6))./4;
    n_sac_Explr_PV(s,2) =  nanmean(Less_trl(:,6))./4;
      
    
    clearvars More_trl Less_trl


end

%% save output for plotting
save([datafigspath,'Fig4C_right.mat'],'n_sac_VCRF_PV')
save([datafigspath,'SuppleFig2_right42.mat'],'n_sac_Explr_PV')
%% descriptive analysis
data(1,:) = mean([n_sac_RF_PV,n_sac_VC_PV]);
data(2,:) = std([n_sac_RF_PV,n_sac_VC_PV]);

T = array2table(data, ...
    'VariableNames', {'Rem', 'Forg','1dva','5dva','21dva'}, ...
    'RowNames', {'Mean', 'SD'});

% Display the table
disp(T);

T2 = array2table([mean(n_sac_VCRF_PV);std(n_sac_VCRF_PV)], ...
    'VariableNames', {'Rem_1dva', 'Rem_5dva','Rem_21dva','Forg_1dva', 'Forg_5dva','Forg_21dva'}, ...
    'RowNames', {'Mean', 'SD'});

% Display the table
disp(T2);

%% Saccades metrics
edge_ang = -180:10:180;
for s = 1:34
    disp(s)
    %% Load Trial level info and Saccade details
    
    load([sacsavepath,sprintf('TrialEyeInfo_%d_%s.mat',s,'Exp4')])
    load([sacsavepath,sprintf('Sac_PV_Info_Sub%d_%s.mat',s,'Exp4')])
        
    % 
    trl_info = Eye_trlinfo(Eye_trlinfo(:,9)==1,1);
    Sac_trials.trialinfo = Sac_trials.trialinfo(ismember(Sac_trials.trialinfo(:,1),trl_info),:);
    
    
    %% Amplitdue
    alltrl = Sac_trials.trialinfo;
    
    sac_amplitude_Exp4(s,1) = nanmean(alltrl(:,5));
    %idv_sac_amp_Exp4 = cat(1,idv_sac_amp_Exp4,alltrl(:,5));
    %% Duration 
    sac_duration_Exp4(s,1) = nanmean(alltrl(:,7)-alltrl(:,6))/1000*1000;

    %idv_sac_dur_Exp4 = cat(1,idv_sac_dur_Exp4,(alltrl(:,7)-alltrl(:,6))./1000*1000);
    %% Direction
    bin_data = histcounts(alltrl(:,4),edge_ang);
    ms_dir_Exp4(:,s) = bin_data./sum(bin_data); clearvars bin_data

end

save([datafigspath,'SuppleFig1_Exp4.mat'],'sac_amplitude_Exp4','sac_duration_Exp4','ms_dir_Exp4')
%% Time resolved saccade frequency

t_ET  = -1.3:1/1000:4;
Hz = 1000;

bin_t = 0.20; % 100 ms per bin
edge = (-1:bin_t:4);%-bin_t/2;

for s = 1:34
    disp(s)
    %% Epoch trial rejection
    load([sacsavepath,sprintf('TrialEyeInfo_%d_%s.mat',s,'Exp4')])
    load([sacsavepath,sprintf('Sac_PV_Info_Sub%d_%s.mat',s,'Exp4')])
    
    trl_info = Eye_trlinfo(Eye_trlinfo(:,9)==1,:);
    Sac_trials.trialinfo = Sac_trials.trialinfo(ismember(Sac_trials.trialinfo(:,1),trl_info(:,1)),:);
    %%
    for i=1:size(Sac_trials.trialinfo,1)
        Sac_trials.trialinfo(i,8) = Eye_trlinfo(ismember(Eye_trlinfo(:,1),Sac_trials.trialinfo(i,1)),4);
    end
    
    VC1trl = Sac_trials.trialinfo(Sac_trials.trialinfo(:,8)==1,:);
    VC2trl = Sac_trials.trialinfo(Sac_trials.trialinfo(:,8)==2,:);
    VC3trl = Sac_trials.trialinfo(Sac_trials.trialinfo(:,8)==3,:);
    
    %% general distribution around pic onset
    t_on = Sac_trials.trialinfo(:,6)+t_ET(1)*1000;% correct the time lock to pic onset
    onset_picon(:,s) = histcounts(t_on,edge.*Hz)./length(unique(Sac_trials.trialinfo(:,1)))*(1/bin_t); % div by trial num and convert to N_sac / sec;
    
    % hit miss
    t_on_hit = Sac_trials.trialinfo(Sac_trials.trialinfo(:,2)==1,6)+t_ET(1)*1000;% correct the time lock to pic onset
    onset_picon_hit(:,s) = histcounts(t_on_hit,edge.*Hz)./length(unique(Sac_trials.trialinfo(Sac_trials.trialinfo(:,2)==1,1)))*(1/bin_t); % div by trial num and convert to N_sac / sec;
    
    t_on_mis = Sac_trials.trialinfo(Sac_trials.trialinfo(:,2)==0,6)+t_ET(1)*1000;% correct the time lock to pic onset
    onset_picon_mis(:,s) = histcounts(t_on_mis,edge.*Hz)./length(unique(Sac_trials.trialinfo(Sac_trials.trialinfo(:,2)==0,1)))*(1/bin_t); % div by trial num and convert to N_sac / sec;
    
    
    % sep by viewing condition
    t_on_C1 = VC1trl(:,6)+t_ET(1)*1000;% correct the time lock to pic onset
    onset_picon_VC1(:,s) = histcounts(t_on_C1,edge.*Hz)./length(unique(VC1trl(:,1)))*(1/bin_t); % div by trial num and convert to N_sac / sec;
    
    t_on_C2 = VC2trl(:,6)+t_ET(1)*1000;% correct the time lock to pic onset
    onset_picon_VC2(:,s) = histcounts(t_on_C2,edge.*Hz)./length(unique(VC2trl(:,1)))*(1/bin_t); % div by trial num and convert to N_sac / sec;
    
    t_on_C3 = VC3trl(:,6)+t_ET(1)*1000;% correct the time lock to pic onset
    onset_picon_VC3(:,s) = histcounts(t_on_C3,edge.*Hz)./length(unique(VC3trl(:,1)))*(1/bin_t); % div by trial num and convert to N_sac / sec;
    

end

%% Stat
subset = 1:34; % all subjects
% rem vs forgotten
for RF=1
    t_range = edge(2:end);
    SacFreq_C1 = [];SacFreq_C2 = [];
    
    SacFreq_C1.avg = permute(onset_picon_hit(:,subset),[3,1,2]);
    SacFreq_C1.time = t_range;
    SacFreq_C1.dimord = 'chan_time_subj';
    SacFreq_C1.label = {'Sac'};
    SacFreq_C2 = SacFreq_C1;
    SacFreq_C2.avg = permute(onset_picon_mis(:,subset),[3,1,2]);

    
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
    cfg.parameter           = 'avg';
    cfg.latency             = [0,4];
    cfg.tail             = 0;
    cfg.clustertail      = 0;
    cfg.alpha            = 0.025;
    cfg.numrandomization = 1000;
    
    cfg.design           = design;
    
    
    [stat_SacFreq]  = ft_timelockstatistics(cfg, SacFreq_C1,SacFreq_C2);
    
end



% ANOVA cluster based on Viewing condition
for VC = 1
    On_VC1 = [];On_VC2 = [];On_VC3 = [];
    t_range = edge(2:end);
    On_VC1.avg = permute(onset_picon_VC1(:,subset),[3,1,2]);
    On_VC1.time = t_range;
    On_VC1.dimord = 'chan_time_subj';
    On_VC1.label = {'eye'};
    On_VC2 = On_VC1;
    On_VC2.avg = permute(onset_picon_VC2(:,subset),[3,1,2]);
    On_VC3 = On_VC1;
    On_VC3.avg = permute(onset_picon_VC3(:,subset),[3,1,2]);
    
    design      = zeros(2,length(subset)*3);
    design(1,:) = repmat(1:length(subset),[1 3]);
    design(2,:) = [ones(1,length(subset)),ones(1,length(subset))+1,ones(1,length(subset))+2];
    
    cfg = [];
    cfg.spmversion = 'spm12';
    cfg.method           = 'montecarlo';
    cfg.statistic        = 'ft_statfun_depsamplesFunivariate';
    cfg.correctm         = 'cluster';
    cfg.clusteralpha     = 0.05;
    cfg.clusterstatistic = 'maxsum';
    cfg.ivar                = 2;
    cfg.uvar                = 1;
    cfg.parameter           = 'avg';
    cfg.latency             = [0,4];
    cfg.tail             = 1;
    cfg.clustertail      = 1;
    cfg.alpha            = 0.05;
    cfg.numrandomization = 1000;
    
    cfg.design           = design;
    
    
    [stat_VC]  = ft_timelockstatistics(cfg, On_VC1,On_VC2,On_VC3);
end

save([datafigspath,'Fig4C_left.mat'],'stat_VC','onset_picon_VC1','onset_picon_VC2','onset_picon_VC3')
save([datafigspath,'Fig4C_middle.mat'],'stat_SacFreq','onset_picon_hit','onset_picon_mis')
%% Figures
codepath = '';

datafigspath =  [codepath,'\data4figs\'];
figsavepath = [codepath,'\figures\'];
% load saccade frequency across time
load( [datafigspath,'Fig4C_left.mat'])
load( [datafigspath,'Fig4C_middle.mat'])
load( [datafigspath,'Fig4C_right.mat'])

plot_save = 1; 

sizefig = [0.2,0.2,0.5,0.7];
cb = cbrewer2('RdBu','div',2);
fig_savename = 'Fig4C_right';
for fig=231
    %% violin plot
    f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig);clf;
    
    h=daboxplot([n_sac_VCRF_PV(:,[1:3]);n_sac_VCRF_PV(:,[4:6])],...
        [repmat(1,size(n_sac_VCRF_PV,1),1);repmat(2,size(n_sac_VCRF_PV,1),1)],...
        'scatter' ,2,...
        'color',cb,...
        'boxalpha', 0.5,...
        'jitter',1,...
        'linkline',1,...
        'withinlines',0,...
        'legend',{'Rem','Forg'},...
        'xtlabels',{'1DVA','5DVA','21DVA'})
    
    h.lg.Location = 'northwest';
    
    ylabel('Saccade Frequency [N/Sec]')
    xlabel('Viewing Condition')
   
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
    if plot_save
    saveas(gcf,[figsavepath,fig_savename,'.emf'])
    saveas(gcf,[figsavepath,fig_savename,'.png'])
    end
end 

fig_savename = 'Fig4C_middle';
sizefig = [0.2,0.2,0.6,0.7];fontaxis = 16;line_width = 2;
ylim_plt = [0,3];% height of the figure
% time bin info
bin_t = 0.20; % 200 ms per bin
edge = (-1:bin_t:4);
t_range = edge(2:end);
t_plt = -1:1:4;
t_plt_ind = dsearchn(t_range',t_plt');

for fig=30

options.handle     = figure('Name',int2str(fig),'units','normalized','outerposition',sizefig);clf;
options.color_area = [128 193 219]./255;    % Blue theme
options.color_line = [ 52 148 186]./255;
options.alpha      = 0.7;
options.line_width = 2;
options.error      = 'sem';

plot_areaerrorbar(onset_picon_hit',options);hold on


options.color_area = [243 169 114]./255;    % Orange theme
options.color_line = [236 112  22]./255;

plot_areaerrorbar(onset_picon_mis',options);hold on

for i_c=1:length(stat_SacFreq.posclusters)
    if stat_SacFreq.posclusters(i_c).prob<0.025
        sig_tw = stat_SacFreq.time(stat_SacFreq.posclusterslabelmat==i_c);
        t_idx= dsearchn(t_range',sig_tw');
        v1 = [[t_idx(1); t_idx(end); t_idx(end); t_idx(1)], [ylim_plt(1); ylim_plt(1); ylim_plt(2); ylim_plt(2)]];
        f1 = [1 2 3 4];
        patch('Faces',f1,'Vertices',v1,'FaceColor','black','FaceAlpha',.1);
    end
end

for i_c=1:length(stat_SacFreq.negclusters)
    if stat_SacFreq.negclusters(i_c).prob<0.025
        sig_tw = stat_SacFreq.time(stat_SacFreq.negclusterslabelmat==i_c);
        t_idx= dsearchn(t_range',sig_tw');
        v1 = [[t_idx(1); t_idx(end); t_idx(end); t_idx(1)], [ylim_plt(1); ylim_plt(1); ylim_plt(2); ylim_plt(2)]];
        f1 = [1 2 3 4];
        patch('Faces',f1,'Vertices',v1,'FaceColor','black','FaceAlpha',.1);
    end
end

xline(find(t_range==0'),'--');
hold off

xticks(t_plt_ind')
xticklabels(num2cell(t_plt))
xlabel('Time [Sec]')
ylabel('Frequency (N/Sec)')
title('Saccades frequency across time')
legend({'Rem','','Forg',''},'Location','northwest')
ylim(ylim_plt)
xlim([1,size(onset_picon_hit,1)])
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',fontaxis,'FontWeight','Bold', 'LineWidth', line_width);

%% whether save fig

if plot_save
    saveas(gcf,[figsavepath,fig_savename,'.emf'])
    saveas(gcf,[figsavepath,fig_savename,'.png'])
end
    
end

fig_savename = 'Fig4C_left';
sizefig = [0.2,0.2,0.6,0.7];fontaxis = 16;line_width = 2;
ylim_plt = [0,3.5];% height of the figure
for fig=30
options.handle     = figure('Name',int2str(fig),'units','normalized','outerposition',sizefig);clf;
options.color_area = [128 193 219]./255;    % Blue theme
options.color_line = [ 52 148 186]./255;
options.alpha      = 0.7;
options.line_width = 2;
options.error      = 'sem';

plot_areaerrorbar(onset_picon_VC1',options);hold on


options.color_area = [193 128  219]./255;    % 
options.color_line = [148 52  186]./255;

plot_areaerrorbar(onset_picon_VC2',options);hold on


options.color_area = [193 219 128  ]./255;    % 
options.color_line = [148  186 52  ]./255;
plot_areaerrorbar(onset_picon_VC3',options);hold on

for i_c=1:length(stat_VC.posclusters)
    if stat_VC.posclusters(i_c).prob<0.05
        sig_tw = stat_VC.time(stat_VC.posclusterslabelmat==i_c);
        t_idx= dsearchn(t_range',sig_tw');
        v1 = [[t_idx(1); t_idx(end); t_idx(end); t_idx(1)], [ylim_plt(1); ylim_plt(1); ylim_plt(2); ylim_plt(2)]];
        f1 = [1 2 3 4];
        patch('Faces',f1,'Vertices',v1,'FaceColor','black','FaceAlpha',.1);
    end
end
xline(find(t_range==0'),'--');
hold off

xticks(t_plt_ind')
xticklabels(num2cell(t_plt))
xlabel('Time [Sec]')
ylabel('Frequency (N/Sec)')
title('Saccades frequency across time')
legend({'1DVA','','5DVA','','21DVA',''},'Location','northwest')
ylim([ylim_plt])
xlim([1,size(onset_picon_VC3,1)])
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',fontaxis,'FontWeight','Bold', 'LineWidth', line_width);

%% whether save fig

if plot_save
    saveas(gcf,[figsavepath,fig_savename,'.emf'])
    saveas(gcf,[figsavepath,fig_savename,'.png'])
end
    
end



