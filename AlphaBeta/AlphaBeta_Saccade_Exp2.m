%% Saccade locked alpha/beta amplitude analysis (Exp2)
codepath = '';

datapath = [codepath,'\datafiles\'];
datafigspath =  [codepath,'\data4figs\'];
sacsavepath = [codepath,'\datafiles\Subjects_Exp2\saccades\'];
eegsavepath = [codepath,'\datafiles\Subjects_Exp2\EEG\'];
%%
freqrange = [10,20]; % frequency of interest
%%
for s = 1:4
    savename = [eegsavepath,sprintf('Sub_Posterior_PV_%d_%s.mat',s,'Exp2')];
    load(savename);
    % Hilbert
    
    %%%% get hilbert filert data
    cfg = [];
    cfg.bpfilter   = 'yes';
    cfg.bpfreq     = freqrange;
    cfg.hilbert    = 'abs'; % get envelope of the filtered data
    cfg.keeptrials = 'yes';
    
    Hilbert_AlphaBeta{s} = ft_preprocessing(cfg,data_all);
    Hilbert_AlphaBeta{s}.Chanwm=data_all.Chanwm;
    Hilbert_AlphaBeta{s}.spiketrial=spike_trl_idx_all;
    Hilbert_AlphaBeta{s}.rejectrial=trial_rej_all;
end
%%
% replace marked artifacts as nans
for s = 1:4
    for nanmark = 1
        % marke here nan for spike around the epochs for eeg & TFR
        for ichan = 1:length(Hilbert_AlphaBeta{s}.label)
            spk_chan = Hilbert_AlphaBeta{s}.spiketrial{ichan};
            spk_overlap = spk_chan(ismember(spk_chan(:,1),Hilbert_AlphaBeta{s}.trialinfo(:,1)),:);
            for ispk = 1:size(spk_overlap,1)
                idx_trl = find(Hilbert_AlphaBeta{s}.trialinfo(:,1)==spk_overlap(ispk,1));
                for itrl = 1:length(idx_trl)
                    
                    tw_ispk = spk_overlap(ispk,[2,3]);
                    tw_eeg = Hilbert_AlphaBeta{s}.time{1};
                    Hilbert_AlphaBeta{s}.trial{idx_trl(itrl)}(ichan,(find(tw_ispk(1)<=tw_eeg,1,'first'):find(tw_ispk(2)>=tw_eeg,1,'last')))=nan;
                    
                end
            end
        end
    end
    
end
%%
nperm = 100;
% set selection parameters
tw = [-0.5,0.5];% extract [-0.5, 0.5] sec to saccade onset
t_eye_min = 0-tw(1);% at least 0.5s after stim on
t_eye_max = 4-tw(2);% latest 3.5s after stim on 

t_range =  tw(1):1/1000:tw(2);

tw_bc = [-0.5,-0.3];
idx_bc = dsearchn(t_range',tw_bc(1)):dsearchn(t_range',tw_bc(2));

chancount = 0;chan_sub = [];
for s = 1:length(Hilbert_AlphaBeta)
    
    %% load saccade and lock to saccade & get infomation about saccade 
    switch s
    case {1,2}
        t_eye = -1.3:1/600:4;
    case {3,4}
        t_eye = -1.3:1/1000:4;
    end

    
    for ichan=1:length(Hilbert_AlphaBeta{s}.label)
        chancount = chancount+1;
        chan_sub(chancount,1) = s;% channels subject index
        dat = cat(3,Hilbert_AlphaBeta{s}.trial{:});
        AlphaBetadat = permute(dat(ichan,:,:),[3,2,1]);clearvars dat
        
        
        % load saccade info
        load([sacsavepath,sprintf('Sac_PV_Info_Sub%d_%s.mat',s,'Exp2')]);
        
        sac_incl = Sac_trials.trialinfo(Sac_trials.trialinfo(:,6)>=dsearchn(t_eye',t_eye_min)& Sac_trials.trialinfo(:,6)<=dsearchn(t_eye',t_eye_max),:);
        sac_incl = sac_incl(ismember(sac_incl(:,1),Hilbert_AlphaBeta{s}.trialinfo),:);
        clean_trl = 1:144;
        clean_trl(Hilbert_AlphaBeta{s}.rejectrial{ichan})=[];
        
        sac_incl = sac_incl(ismember(sac_incl(:,1),clean_trl),:);
        %%
        eeg_trlidx = Hilbert_AlphaBeta{s}.trialinfo(:,1);
        t_eeg = Hilbert_AlphaBeta{s}.time{1};
        for isac = 1:length(sac_incl)
            sac_time = t_eye(sac_incl(isac,6));
            AlphaBeta_Sac_locked(chancount).dat{isac} = AlphaBetadat(eeg_trlidx==sac_incl(isac,1),[dsearchn(t_eeg',sac_time)+tw(1)*1000:dsearchn(t_eeg',sac_time)+tw(2)*1000]);
            
        end
        AlphaBeta_Sac_locked(chancount).sacinfo = sac_incl;
        
        %% add random shuffle
        for iperm = 1:nperm
            if rem(iperm,nperm/10)==0
                fprintf('%d %% Complete \n',iperm/nperm*100);
            end

            rand_idx =randperm(size(sac_incl,1))';
            rand_trl = randperm(size(AlphaBetadat,1))';
            AlphaBeta_perm = AlphaBetadat;
            sac_incl_perm = sac_incl;
            sac_incl_perm(:,2:end) = sac_incl_perm(rand_idx,2:end);
            for isac = 1:length(sac_incl_perm)
                sac_time = t_eye(sac_incl_perm(isac,6));
                data_rand(isac,:) = AlphaBeta_perm(eeg_trlidx==sac_incl(isac,1),[dsearchn(t_eeg',sac_time)+tw(1)*1000:dsearchn(t_eeg',sac_time)+tw(2)*1000]);

            end

            AlphaBeta_Sac_locked_rand(chancount).dat{iperm} = nanmean(data_rand-mean(data_rand(:,idx_bc),2),1);% demean at trl level
            AlphaBeta_Sac_locked_rand(chancount).sac_sh_idx{iperm} = rand_idx;

            clearvars rand_idx sac_incl_perm data_rand rand_trl

        end



    end
end

save([datapath,'AlphaBeta_hilbert_Sac_real_Exp2.mat'],'AlphaBeta_Sac_locked','chan_sub','-v7.3'); 
%save([datapath,'AlphaBeta_hilbert_Sac_shuffle_Exp2.mat'],'AlphaBeta_Sac_locked_rand','-v7.3'); 
%% Separate by conditions
load([datapath,'AlphaBeta_hilbert_Sac_real_Exp2.mat']);
load([datapath,'AlphaBeta_hilbert_Sac_shuffle_Exp2.mat'])

t_range =  -0.5:1/1000:0.5;

tw_bc = [-0.5,-0.3];% baseline correction 
idx_bc = dsearchn(t_range',tw_bc(1)):dsearchn(t_range',tw_bc(2));

for s = 1:length(AlphaBeta_Sac_locked)
    disp(s)
    
    % baseline correction
    bc_avg =  cellfun(@(x) mean(x(idx_bc)),AlphaBeta_Sac_locked(s).dat,'Un',0);
    AlphaBeta_Sac_locked(s).dat   = cellfun(@(x,y) x-y,AlphaBeta_Sac_locked(s).dat, bc_avg,'Un',0);
        
    %% general and rand

    AlphaBeta_randsh(s,:) = nanmean(cat(1,AlphaBeta_Sac_locked_rand(s).dat{:}));
    AlphaBeta_allsac(s,:) = nanmean(cat(1,AlphaBeta_Sac_locked(s).dat{:}));

    %% Separate base on Memory performance
    AlphaBeta_Rem(s,:) = nanmean(cat(1,AlphaBeta_Sac_locked(s).dat{(AlphaBeta_Sac_locked(s).sacinfo(:,2)==1)}));
    AlphaBeta_For(s,:) = nanmean(cat(1,AlphaBeta_Sac_locked(s).dat{(AlphaBeta_Sac_locked(s).sacinfo(:,2)==0)}));

    % get sac property as well
    
    Sac_amp_Memory(s,1) = nanmean(AlphaBeta_Sac_locked(s).sacinfo(AlphaBeta_Sac_locked(s).sacinfo(:,2)==1,5));
    Sac_amp_Memory(s,2) = nanmean(AlphaBeta_Sac_locked(s).sacinfo(AlphaBeta_Sac_locked(s).sacinfo(:,2)==0,5));

    %% short vs large saccade
    cri1 = quantile(AlphaBeta_Sac_locked(s).sacinfo(:,5),1/3);
    cri2 = quantile(AlphaBeta_Sac_locked(s).sacinfo(:,5),2/3);
    
    AlphaBeta_Amp1(s,:) = nanmean(cat(1,AlphaBeta_Sac_locked(s).dat{AlphaBeta_Sac_locked(s).sacinfo(:,5)<cri1}));
    AlphaBeta_Amp2(s,:) = nanmean(cat(1,AlphaBeta_Sac_locked(s).dat{AlphaBeta_Sac_locked(s).sacinfo(:,5)>=cri1&AlphaBeta_Sac_locked(s).sacinfo(:,5)<=cri2}));
    AlphaBeta_Amp3(s,:) = nanmean(cat(1,AlphaBeta_Sac_locked(s).dat{AlphaBeta_Sac_locked(s).sacinfo(:,5)>cri2}));
    
    % saccade offset time
    Sac_offset_amp(s,1) = nanmean(diff(AlphaBeta_Sac_locked(s).sacinfo(AlphaBeta_Sac_locked(s).sacinfo(:,5)<cri1,[6,7]),[],2));
    Sac_offset_amp(s,2) = nanmean(diff(AlphaBeta_Sac_locked(s).sacinfo(AlphaBeta_Sac_locked(s).sacinfo(:,5)>=cri1&AlphaBeta_Sac_locked(s).sacinfo(:,5)<=cri2,[6,7]),[],2));
    Sac_offset_amp(s,3) = nanmean(diff(AlphaBeta_Sac_locked(s).sacinfo(AlphaBeta_Sac_locked(s).sacinfo(:,5)>cri2,[6,7]),[],2));
    % saccade amplitdue
    Sac_amp_amp(s,1)    = nanmean(AlphaBeta_Sac_locked(s).sacinfo(AlphaBeta_Sac_locked(s).sacinfo(:,5)<cri1,5));
    Sac_amp_amp(s,2)    = nanmean(AlphaBeta_Sac_locked(s).sacinfo(AlphaBeta_Sac_locked(s).sacinfo(:,5)>=cri1 & AlphaBeta_Sac_locked(s).sacinfo(:,5)<=cri2,5));
    Sac_amp_amp(s,3)    = nanmean(AlphaBeta_Sac_locked(s).sacinfo(AlphaBeta_Sac_locked(s).sacinfo(:,5)>cri2,5));
    
end

%% stat against shuffle
% All sac vs rand perm
t_range =  -0.5:1/1000:0.5;
subset = 1:20;
for runstat = 1
    Alpha_C1 = [];Alpha_C2 = [];
    
    
    Alpha_C1.avg = permute(AlphaBeta_allsac(subset,:),[3,2,1]);
    Alpha_C1.time = t_range;
    Alpha_C1.dimord = 'chan_time_subj';
    Alpha_C1.label = {'dat'};
    Alpha_C2 = Alpha_C1;
    Alpha_C2.avg = permute(AlphaBeta_randsh(subset,:),[3,2,1]);

    
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
    cfg.latency             = [-0.3,0.5];
    cfg.tail             = 0;
    cfg.clustertail      = 0;
    cfg.alpha            = 0.025;
    cfg.numrandomization = 1000;
    
    cfg.design           = design;
    
    
    [stat_Gen]  = ft_timelockstatistics(cfg, Alpha_C1,Alpha_C2);
end
%% Identify the peak latency and peak 

% rem vs forg
t_min = find(t_range>=0,1,'first');
t_max = find(t_range<=0.2,1,'last');
for s = 1:size(AlphaBeta_Rem,1)
    [a,b] = findpeaks(AlphaBeta_Rem(s,:),'SortStr','descend');
    peak_latency_mem(s,1) =t_range(b(find(b>t_min & b<t_max,1,'first')));
    peak_ampitude_mem(s,1)=a(find(b>t_min & b<t_max,1,'first'));
    
    [a,b] = findpeaks(AlphaBeta_For(s,:),'SortStr','descend');
    peak_latency_mem(s,2) =t_range(b(find(b>t_min & b<t_max,1,'first')));
    peak_ampitude_mem(s,2)=a(find(b>t_min & b<t_max,1,'first'));
    

end

% sep by saccade size
t_min = find(t_range>=0,1,'first');
t_max = find(t_range<=0.2,1,'last');
for s = 1:size(AlphaBeta_Amp1,1)
    [a,b] = findpeaks(AlphaBeta_Amp1(s,:),'SortStr','descend');
    peak_latency_amp(s,1) =t_range(b(find(b>t_min & b<t_max,1,'first')));
    peak_ampitude_amp(s,1)=a(find(b>t_min & b<t_max,1,'first'));
    
    [a,b] = findpeaks(AlphaBeta_Amp2(s,:),'SortStr','descend');
    peak_latency_amp(s,2) =t_range(b(find(b>t_min & b<t_max,1,'first')));
    peak_ampitude_amp(s,2)=a(find(b>t_min & b<t_max,1,'first'));
    
    [a,b] = findpeaks(AlphaBeta_Amp3(s,:),'SortStr','descend');
    peak_latency_amp(s,3) =t_range(b(find(b>t_min & b<t_max,1,'first')));
    peak_ampitude_amp(s,3)=a(find(b>t_min & b<t_max,1,'first'));
end

%% save output for plotting
save([datafigspath,'Fig5B.mat'],'Sac_amp_amp','Sac_offset_amp','peak_latency_amp','peak_ampitude_amp','stat_Gen','AlphaBeta_allsac','AlphaBeta_randsh','AlphaBeta_Amp1','AlphaBeta_Amp2','AlphaBeta_Amp3','chan_sub')
save([datafigspath,'SuppleFig5B.mat'],'Sac_amp_Memory','peak_latency_mem','peak_ampitude_mem','AlphaBeta_Rem','AlphaBeta_For')
%% Figure (Fig 5B)
codepath = '';

datafigspath =  [codepath,'\data4figs\'];
figsavepath = [codepath,'\figures\'];
% load data
load( [datafigspath,'Fig5B.mat'])

% descriptive 

T1 = array2table([mean(peak_ampitude_amp);std(peak_ampitude_amp)], ...
    'VariableNames', {'PeakAmp_S', 'PeakAmp_M','PeakAmp_L'}, ...
    'RowNames', {'Mean', 'SD'});
T2 = array2table([mean(peak_latency_amp);std(peak_latency_amp)], ...
    'VariableNames', {'PeakLat_S', 'PeakLat_M','PeakLat_L'}, ...
    'RowNames', {'Mean', 'SD'});
% Display the table
disp(T1);
disp(T2);

%%
subset  = 1:20;% all channels
t_range =  -0.5:1/1000:0.5;
toi = [-0.3,0.5];% time window after baseline
toi_idx  = dsearchn(t_range',toi');
t_plt = toi(1):0.2:toi(2);
t_plt_ind = dsearchn(t_range(toi_idx(1):toi_idx(2))',t_plt');

ylim_plt = [-0.6,0.8];
sizefig = [0.1,0.1,0.6,0.6];
plot_save = 0;

fig_savename = 'Fig5B_bottomleft';
for fig=33
    
f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig);clf;    
options.handle     = figure(f);
options.color_area = [133 128 177]./255;    % Blue theme
options.color_line = [133 128 177]./255;
options.alpha      = 0.7;
options.line_width = 2;
options.error      = 'sem';

plot_areaerrorbar(AlphaBeta_allsac(subset,toi_idx(1):toi_idx(2)),options);hold on
options.color_area = [33 128 66]./255;    % Blue theme
options.color_line = [33 128 66]./255;

plot_areaerrorbar(AlphaBeta_randsh(subset,toi_idx(1):toi_idx(2)),options);hold on

for i_c=1:length(stat_Gen.posclusters)
    if stat_Gen.posclusters(i_c).prob<0.025
        sig_tw = stat_Gen.time(stat_Gen.posclusterslabelmat==i_c);
        t_idx= dsearchn(t_range(toi_idx(1):toi_idx(2))',sig_tw');
        v1 = [[t_idx(1); t_idx(end); t_idx(end); t_idx(1)], [ylim_plt(1); ylim_plt(1); ylim_plt(2); ylim_plt(2)]];
        f1 = [1 2 3 4];
        patch('Faces',f1,'Vertices',v1,'FaceColor','black','FaceAlpha',.1);
    end
end

for i_c=1:length(stat_Gen.negclusters)
    if stat_Gen.negclusters(i_c).prob<0.025
        sig_tw = stat_Gen.time(stat_Gen.negclusterslabelmat==i_c);
        t_idx= dsearchn(t_range(toi_idx(1):toi_idx(2))',sig_tw');
        v1 = [[t_idx(1); t_idx(end); t_idx(end); t_idx(1)], [ylim_plt(1); ylim_plt(1); ylim_plt(2); ylim_plt(2)]];
        f1 = [1 2 3 4];
        patch('Faces',f1,'Vertices',v1,'FaceColor','black','FaceAlpha',.1);
    end
end

xline(find(t_range(toi_idx(1):toi_idx(2))==0),'--')
hold on
        
xticks(t_plt_ind')
xticklabels(num2cell(t_plt*1000))
xlabel('Time [ms]')
ylabel('Alpha Beta amplitude')
title('Saccade locked Alpha Beta')
xlim([1,length(toi_idx(1):toi_idx(2))])
ylim(ylim_plt)
legend({'Saccade','','Shuffle',''})
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);


if plot_save
    saveas(gcf,[figsavepath,fig_savename,'.emf'])
    saveas(gcf,[figsavepath,fig_savename,'.png'])
end
end

fig_savename = 'Fig5B_bottomright';
for fig=33     
f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig);clf;    
options.handle     = figure(f);
options.color_area = [128 193 219]./255;    % Blue theme
options.color_line = [ 52 148 186]./255;
options.alpha      = 0.7;
options.line_width = 2;
options.error      = 'sem';

plot_areaerrorbar(AlphaBeta_Amp1(:,toi_idx(1):toi_idx(2)),options);hold on

options.color_area = [193 128  219]./255;    % 
options.color_line = [148 52  186]./255;

plot_areaerrorbar(AlphaBeta_Amp2(:,toi_idx(1):toi_idx(2)),options);hold on

options.color_area = [193 219 128  ]./255;    % 
options.color_line = [148  186 52  ]./255;

plot_areaerrorbar(AlphaBeta_Amp3(:,toi_idx(1):toi_idx(2)),options);hold on

xline(find(t_range(toi_idx(1):toi_idx(2))==0),'-')  

t_find = [0,0.2];% search within first 0.2 sec after saccade onset
toi_find  = dsearchn(t_range',t_find');
[~,loc1]=max(mean(AlphaBeta_Amp1(:,toi_find(1):toi_find(2))));
[~,loc2]=max(mean(AlphaBeta_Amp2(:,toi_find(1):toi_find(2))));
[~,loc3]=max(mean(AlphaBeta_Amp3(:,toi_find(1):toi_find(2))));
loc1 = loc1+toi_find(1)-toi_idx(1);
loc2 = loc2+toi_find(1)-toi_idx(1);
loc3 = loc3+toi_find(1)-toi_idx(1);
% mark group average max
xline(loc1,'color',[ 52 148 186]./255,'LineStyle','-','LineWidth',2) 
xline(loc2,'color',[193 128  219]./255,'LineStyle','-','LineWidth',2) 
xline(loc3,'color',[193 219 128  ]./255,'LineStyle','-','LineWidth',2) 

% add saccade offset points
scatter(Sac_offset_amp(:,1)+300,repmat(-0.3,size(Sac_offset_amp,1),1),'filled','MarkerFaceColor',[ 52 148 186]./255,'LineWidth',1);
scatter(Sac_offset_amp(:,2)+300,repmat(-0.35,size(Sac_offset_amp,1),1),'filled','MarkerFaceColor',[193 128  219]./255,'LineWidth',1);
scatter(Sac_offset_amp(:,3)+300,repmat(-0.4,size(Sac_offset_amp,1),1),'filled','MarkerFaceColor',[193 219 128  ]./255,'LineWidth',1);hold off

xticks(t_plt_ind')
xticklabels(num2cell(t_plt*1000))
xlabel('Time [ms]')
ylabel('Alpha Beta amplitude')
title('Posterior Alpha Beta')
xlim([1,length(toi_idx(1):toi_idx(2))])
ylim(ylim_plt)
legend({'Small','','Medium','','Large',''})
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);


if plot_save
    saveas(gcf,[figsavepath,fig_savename,'.emf'])
    saveas(gcf,[figsavepath,fig_savename,'.png'])
end
end

fig_savename = 'Fig5B_topmiddle';
sizefig = [0.1,0.1,0.20,0.6];
cb = cbrewer2('PuBuGn','div',3);
% saccade at subject level
for s=1:4
    Sac_amp_amp_sub(s,:) = mean(Sac_amp_amp(chan_sub==s,:),1);
     
end
for fig=21
    %% 
    f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig);clf;
    h=daboxplot([Sac_amp_amp_sub(:,[1:3])],[repmat(1,size(Sac_amp_amp_sub,1),1)],...
        'scatter' ,2,...
        'color',cb,...
        'boxalpha', 0.5,...
        'jitter',1,...
        'linkline',0,...
        'withinlines',1,...
        'xtlabels',{'Small','Medium','Large'})
    
    ylabel('Amplitude [dva]')
    xlabel('Size')
   
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
    
    if plot_save
    saveas(gcf,[figsavepath,fig_savename,'.emf'])
    saveas(gcf,[figsavepath,fig_savename,'.png'])
    end
end 



lt_lim  = [0,0.2];
amp_lim = [-0.5,3];
% scatter plots of both
sizefig = [0.1,0.1,0.35,0.6];
fig_savename = 'Fig5B_topright1';% scatter plot
cb = cbrewer2('Dark2',3);
for fig=21
    %% 
    f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig);clf;
    scatter(peak_latency_amp(:,1),peak_ampitude_amp(:,1),[],cb(1,:),'filled');hold on
    scatter(peak_latency_amp(:,2),peak_ampitude_amp(:,2),[],cb(2,:),'filled');hold on
    scatter(peak_latency_amp(:,3),peak_ampitude_amp(:,3),[],cb(3,:),'filled');hold off
    xlim(lt_lim)
    ylim(amp_lim)
    
    ylabel('Peak Amplitude [mV]')
    xlabel('Peak Latency [Sec]')
    legend({'S','M','L'},'Location','northwest')
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
    
   
    
    if plot_save
    saveas(gcf,[figsavepath,fig_savename,'.emf'])
    saveas(gcf,[figsavepath,fig_savename,'.png'])
    end
end 



sizefig = [0.1,0.1,0.15,0.6];
cb = repmat([0.5,0.5,0.5],3,1);
fig_savename = 'Fig5B_topright2';% Peak amplitdue
for fig=21
    %% 
    f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig);clf;
    h=dabarplot([peak_ampitude_amp(:,[3:-1:1])],[repmat(1,size(peak_ampitude_amp,1),1)],...
        'scatter' ,0,...
        'errorhats',0,....
        'colors',cb,...
        'jitter',1,...
        'xtlabels',{'L','M','S'})
    
    ylabel('Amplitude [dva]')
    ylim(amp_lim)
   
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
    
   
    
    if plot_save
    saveas(gcf,[figsavepath,fig_savename,'.emf'])
    saveas(gcf,[figsavepath,fig_savename,'.png'])
    end
end 


fig_savename = 'Fig5B_topright3';% Peak latency
for fig=21
    %% 
    f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig);clf;
    h=dabarplot([peak_latency_amp(:,[1:3])],[repmat(1,size(peak_latency_amp,1),1)],...
        'scatter' ,0,...
        'errorhats',0,....
        'colors',cb,...
        'jitter',1,...
        'xtlabels',{'S','M','L'})
    
    ylabel('Time [Sec]')
    ylim(lt_lim)
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
    
   
    
    if plot_save
    saveas(gcf,[figsavepath,fig_savename,'.emf'])
    saveas(gcf,[figsavepath,fig_savename,'.png'])
    end
end 

%% Supplementary Figure
codepath = '';

datafigspath =  [codepath,'\data4figs\'];
figsavepath = [codepath,'\figures\'];
% load data
load( [datafigspath,'SuppleFig5B.mat'])

% time parameter
t_range =  -0.5:1/1000:0.5;
toi = [-0.3,0.5];% time window after baseline
toi_idx  = dsearchn(t_range',toi');
t_plt = toi(1):0.2:toi(2);
t_plt_ind = dsearchn(t_range(toi_idx(1):toi_idx(2))',t_plt');

sizefig = [0.1,0.1,0.6,0.6];
ylim_plt = [-0.6,0.8];
plot_save = 1;
%  cluster based on Rem vs Forg
subset = 1:20;

fig_savename = 'SuppleFig5B_right';
for fig=33
    
f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig);clf;    
options.handle     = figure(f);
options.color_area = [219 128 193 ]./255;    % Blue theme
options.color_line = [ 186 52 148 ]./255;
options.alpha      = 0.7;
options.line_width = 2;
options.error      = 'sem';

plot_areaerrorbar(AlphaBeta_Rem(subset,toi_idx(1):toi_idx(2)),options);hold on
options.color_area = [219 193 128  ]./255;    % 
options.color_line = [ 186 148 52  ]./255;

plot_areaerrorbar(AlphaBeta_For(subset,toi_idx(1):toi_idx(2)),options);hold on

xline(find(t_range(toi_idx(1):toi_idx(2))==0),'--')
        
xticks(t_plt_ind')
xticklabels(num2cell(t_plt*1000))
xlabel('Time [ms]')
ylabel('Alpha Beta amplitude')
title('Posterior Alpha Beta')
xlim([1,length(toi_idx(1):toi_idx(2))])
ylim(ylim_plt)
legend({'Rem','','Forg',''})
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);


if plot_save
    saveas(gcf,[figsavepath,fig_savename,'.emf'])
    saveas(gcf,[figsavepath,fig_savename,'.png'])
end
end


lt_lim  = [0,0.2];
amp_lim = [-0.5,3];
% scatter plots of both
sizefig = [0.1,0.1,0.35,0.6];
fig_savename = 'SuppleFig5B_left1';
cb = cbrewer2('RdBu','div',2);
for fig=21
    %% 
    f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig);clf;
    scatter(peak_latency_mem(:,1),peak_ampitude_mem(:,1),[],cb(1,:),'filled');hold on
    scatter(peak_latency_mem(:,2),peak_ampitude_mem(:,2),[],cb(2,:),'filled');hold off
    xlim(lt_lim)
    ylim(amp_lim)
    
    ylabel('Peak Amplitude [mV]')
    xlabel('Peak Latency [Sec]')
    legend({'Rem','Forg'},'Location','northwest')
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
    
   
    
    if plot_save
    saveas(gcf,[figsavepath,fig_savename,'.emf'])
    saveas(gcf,[figsavepath,fig_savename,'.png'])
    end
end 


sizefig = [0.1,0.1,0.15,0.6];
cb = repmat([0.5,0.5,0.5],2,1);

fig_savename = 'SuppleFig5B_left2';
for fig=21
    %% 
    f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig);clf;
    h=dabarplot([peak_ampitude_mem(:,[2:-1:1])],[repmat(1,size(peak_ampitude_mem,1),1)],...
        'scatter' ,0,...
        'errorhats',0,....
        'colors',cb,...
        'jitter',1,...
        'xtlabels',{'Forg','Rem'})
    
    ylabel('Amplitude [mV]')
    ylim(amp_lim)
   
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
    
   
    
    if plot_save
    saveas(gcf,[figsavepath,fig_savename,'.emf'])
    saveas(gcf,[figsavepath,fig_savename,'.png'])
    end
end 


fig_savename = 'SuppleFig5B_left3';% Peak latency
for fig=21
    %% 
    f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig);clf;
    h=dabarplot([peak_latency_mem(:,[1:2])],[repmat(1,size(peak_latency_mem,1),1)],...
        'scatter' ,0,...
        'errorhats',0,....
        'colors',cb,...
        'jitter',1,...
        'xtlabels',{'Rem','Forg'})
    
    ylabel('Time [Sec]')
    ylim(lt_lim)
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
    
   
    
    if plot_save
    saveas(gcf,[figsavepath,fig_savename,'.emf'])
    saveas(gcf,[figsavepath,fig_savename,'.png'])
    end
end   

