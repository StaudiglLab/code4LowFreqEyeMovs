%% Saccade analysis for Experiment 1
% run this after completing saccades detection
codepath = '';

datapath = [codepath,'\datafiles\'];
datafigspath =  [codepath,'\data4figs\'];
figsavepath = [codepath,'\figures\'];
subinfopath = [codepath,'\datafiles\Subjects_Exp1\ExpInfo\'];
sacsavepath = [codepath,'\datafiles\Subjects_Exp1\saccades\'];
%%

% general N of saccades
for s = 1:20
    disp(s)
    load([sacsavepath,sprintf('TrialEyeInfo_%d_%s.mat',s,'Exp1')])
    
    %% Sep by Rem vs Forg
    Rem_trl = Eye_trlinfo(Eye_trlinfo(:,8)==1 & Eye_trlinfo(:,2)==1,:);
    For_trl = Eye_trlinfo(Eye_trlinfo(:,8)==1 & Eye_trlinfo(:,2)==0,:);
    n_sac_RF_bs(s,1) =  nanmean(Rem_trl(:,4));
    n_sac_RF_bs(s,2) =  nanmean(For_trl(:,4));
    
    n_sac_RF_PV(s,1) =  nanmean(Rem_trl(:,5))./4;
    n_sac_RF_PV(s,2) =  nanmean(For_trl(:,5))./4;
    clearvars Rem_trl For_trl
    %% Sey by N saccades 
    cri = nanmedian(Eye_trlinfo(Eye_trlinfo(:,8)==1,5));
    More_trl = Eye_trlinfo(Eye_trlinfo(:,8)==1 & Eye_trlinfo(:,5)>=cri,:);
    Less_trl = Eye_trlinfo(Eye_trlinfo(:,8)==1 & Eye_trlinfo(:,5)<=cri,:);
    
    n_sac_Nsac_bs(s,1) =  nanmean(More_trl(:,4));
    n_sac_Nsac_bs(s,2) =  nanmean(Less_trl(:,4));
    
    n_sac_Nsac_PV(s,1) =  nanmean(More_trl(:,5))./4;
    n_sac_Nsac_PV(s,2) =  nanmean(Less_trl(:,5))./4;
    
    clearvars More_trl Less_trl 
    %% Sey by Exploration index
    cri = nanmedian(Eye_trlinfo(Eye_trlinfo(:,8)==1,7));
    More_trl = Eye_trlinfo(Eye_trlinfo(:,8)==1 & Eye_trlinfo(:,7)<=cri,:);
    Less_trl = Eye_trlinfo(Eye_trlinfo(:,8)==1 & Eye_trlinfo(:,7)>=cri,:);

    
    n_sac_Explr_PV(s,1) =  nanmean(More_trl(:,5))./4;
    n_sac_Explr_PV(s,2) =  nanmean(Less_trl(:,5))./4;
      
    
    clearvars More_trl Less_trl


end
save([datafigspath,'SuppleFig2_right12.mat'],'n_sac_Explr_PV')
%% descriptive analysis
data(1,:) = mean(n_sac_RF_PV);
data(2,:) = std(n_sac_RF_PV);

T = array2table(data, ...
    'VariableNames', {'Rem', 'Forg'}, ...
    'RowNames', {'Mean', 'SD'});

% Display the table
disp(T);
%% Saccade metrics
edge_ang = -180:10:180;
%idv_sac_amp_Exp1 = [];idv_sac_dur_Exp1 = [];
for s = 1:20
    disp(s)
    
    % load trial info
    load([sacsavepath,sprintf('TrialEyeInfo_%d_%s.mat',s,'Exp1')])
    % load saccade info
    load([sacsavepath,sprintf('Sac_PV_Info_Sub%d_%s.mat',s,'Exp1')]);
    %
    trl_info = Eye_trlinfo(Eye_trlinfo(:,8)==1,1);
    Sac_trials.trialinfo = Sac_trials.trialinfo(ismember(Sac_trials.trialinfo(:,1),trl_info),:);
    
  
     %% Amplitdue
    alltrl = Sac_trials.trialinfo;
    
    sac_amplitude_Exp1(s,1) = nanmean(alltrl(:,5));
    %idv_sac_amp_Exp1 = cat(1,idv_sac_amp_Exp1,alltrl(:,5));
    %% Duration 
    sac_duration_Exp1(s,1) = nanmean(alltrl(:,7)-alltrl(:,6))/600*1000;

    %idv_sac_dur_Exp1 = cat(1,idv_sac_dur_Exp1,(alltrl(:,7)-alltrl(:,6))./600*1000);
    %% Direction
    bin_data = histcounts(alltrl(:,4),edge_ang);
    ms_dir_Exp1(:,s) = bin_data./sum(bin_data); clearvars bin_data
    
   
end

%% save output for plotting
save([datafigspath,'SuppleFig1_Exp1.mat'],'sac_amplitude_Exp1','sac_duration_Exp1','ms_dir_Exp1')

%% Time-resolved Saccade Frequency

t_ET  = -1.3:1/600:4;
Hz = 600;

bin_t = 0.20; % 200 ms per bin
edge = (-1:bin_t:4);%-bin_t/2;

for s = 1:20
    disp(s)
    % load trial info
    load([sacsavepath,sprintf('TrialEyeInfo_%d_%s.mat',s,'Exp1')])
    % load saccade info
    load([sacsavepath,sprintf('Sac_PV_Info_Sub%d_%s.mat',s,'Exp1')]);
    
    trl_info = Eye_trlinfo(Eye_trlinfo(:,8)==1,:);
    Sac_trials.trialinfo = Sac_trials.trialinfo(ismember(Sac_trials.trialinfo(:,1),trl_info(:,1)),:);
    
    clearvars cri trl_info
    %% general distribution around pic onset
    t_on = Sac_trials.trialinfo(:,6)+t_ET(1)*600;% correct the time lock to pic onset
    onset_picon(:,s) = histcounts(t_on,edge.*Hz)./length(unique(Sac_trials.trialinfo(:,1)))*(1/bin_t); % div by trial num and convert to N_sac / sec;
    
    % hit miss
    t_on_hit = Sac_trials.trialinfo(Sac_trials.trialinfo(:,2)==1,6)+t_ET(1)*600;% correct the time lock to pic onset
    onset_picon_hit(:,s) = histcounts(t_on_hit,edge.*Hz)./length(unique(Sac_trials.trialinfo(Sac_trials.trialinfo(:,2)==1,1)))*(1/bin_t); % div by trial num and convert to N_sac / sec;
    
    t_on_mis = Sac_trials.trialinfo(Sac_trials.trialinfo(:,2)==0,6)+t_ET(1)*600;% correct the time lock to pic onset
    onset_picon_mis(:,s) = histcounts(t_on_mis,edge.*Hz)./length(unique(Sac_trials.trialinfo(Sac_trials.trialinfo(:,2)==0,1)))*(1/bin_t); % div by trial num and convert to N_sac / sec;
    
end

%% Stat
subset = 1:20;
t_range = edge(2:end);
for stattest=1
    
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


save([datafigspath,'Fig1B_right.mat'],'stat_SacFreq','onset_picon_hit','onset_picon_mis')
%% Figure

codepath = '';

datafigspath =  [codepath,'\data4figs\'];
figsavepath = [codepath,'\figures\'];
% load saccade frequency across time
load( [datafigspath,'Fig1B_right.mat'])

plot_save = 1;
bin_t = 0.20; % 200 ms per bin
edge = (-1:bin_t:4);%-bin_t/2;
t_range = edge(2:end);
t_plt = -1:1:4;
t_plt_ind = dsearchn(t_range',t_plt');
ylim_plt = [0,4];


sizefig = [0.2,0.2,0.6,0.7];fontaxis = 16;line_width = 2;
fig_savename = 'Fig1B_right';
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


