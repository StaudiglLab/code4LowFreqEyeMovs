%% alpha/beta amplitude following a saccade (Exp2)
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

% replace nan
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
tw = [0,1];

t_eye_min = 0.5;% at least 0.5s after stim on
t_eye_max = 4-tw(2);% latest 3s before stim off

t_limeeg1 = tw(1);% alpha/beta get  sec before saccades
t_limeeg2 = tw(2);% alpha/beta get  sec before saccades

chancount = 0;
for s = 1:length(Hilbert_AlphaBeta)
    
    %% load saccade and lock to saccade & get infomation about saccade 
    switch s
    case {1,2}
        t_eye = -1.3:1/600:4;
        eye_Hz = 600;
    case {3,4}
        t_eye = -1.3:1/1000:4;
        eye_Hz = 1000;
    end

    for ichan=1:length(Hilbert_AlphaBeta{s}.label)
        chancount = chancount+1;
        
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
        
        %% add section for next saccades info
        
        trlidx = unique(sac_incl(:,1));
        s_count = 0;
        next_sacinfo_all = [];
        for itrl = 1:length(trlidx)
            sac_trl_eeg = sac_incl(sac_incl(:,1)==(trlidx(itrl)),:);
            sac_trl_all = Sac_trials.trialinfo(Sac_trials.trialinfo(:,1)==(trlidx(itrl)),:);
            next_sacinfo = [];
            for isac = 1:size(sac_trl_eeg,1)
                t = sac_trl_all(:,6)-sac_trl_eeg(isac,6);
                next_sacinfo(isac,1) = sum(t<0 & t>=t_limeeg1*eye_Hz); % How many Saccades before  s to current
                next_sacinfo(isac,2) = sum(t>0 & t<=t_limeeg2*eye_Hz); % How many Saccades between current to s after;
                if sum(t<0)
                    next_sacinfo(isac,3) = t(find(t<0,1,'last'))/eye_Hz;
                else
                    next_sacinfo(isac,3) = nan;
                end
                
                if sum(t>0)
                    next_sacinfo(isac,4) = t(find(t>0,1,'first'))/eye_Hz;
                else
                    next_sacinfo(isac,4)  = nan;
                end
                
                s_count = s_count+1;
                all_other_sac(chancount).sac{s_count} = t/eye_Hz;
                
            end
            next_sacinfo_all = cat(1,next_sacinfo_all,next_sacinfo);clearvars next_sacinfo
        end
        
        AlphaBeta_Sac_locked(chancount).sacinfo(:,8:11) = next_sacinfo_all;
        
    end
end
%% save output
save([datapath,'AlphaBeta_hilbert_AfterSac_Exp2.mat'],'AlphaBeta_Sac_locked','all_other_sac','-v7.3');
%%
load([datapath,'AlphaBeta_hilbert_AfterSac_Exp2.mat']);

t_range =  0:1/1000:1;

All_sac_AlphaBeta=[];all_sac_info=[];
tw_bc = [0,0.2];% baseline correct to [0,0.2] sec after saccade onset
idx_bc = dsearchn(t_range',tw_bc(1)):dsearchn(t_range',tw_bc(2));
%next_sac_time = [];
allnext_sac_time = [];

for s =  1:length(AlphaBeta_Sac_locked)
    % bc (if whole trl dmean, just select the tw_bc to be thw whole trl length)
    bc_avg =  cellfun(@(x) mean(x(idx_bc)),AlphaBeta_Sac_locked(s).dat,'Un',0);
    
    Datpos   = cellfun(@(x,y) x-y,AlphaBeta_Sac_locked(s).dat, bc_avg,'Un',0);
    Datpos = cat(1,Datpos{:});
    %% 
    AlphaBeta_NoPosSac(s,:) = nanmean(Datpos(AlphaBeta_Sac_locked(s).sacinfo(:,9)==0,:));
    AlphaBeta_WithPosSac(s,:) = nanmean(Datpos(AlphaBeta_Sac_locked(s).sacinfo(:,9)~=0,:));
    
    %%
    %next_sac_time = cat(1,next_sac_time,AlphaBeta_Sac_locked(s).sacinfo(AlphaBeta_Sac_locked(s).sacinfo(:,9)~=0,11));
    
    allnext_sac_time = cat(1,allnext_sac_time,all_other_sac(s).sac{:});
end

allnext_sac_time = allnext_sac_time(allnext_sac_time>0 & allnext_sac_time<1);% subsequent saccades within the next 1 sec

%% stat

%  compare alpha beta between with vs without subsequent saccade
for runstat = 1
    AlphaBeta_C1 = [];AlphaBeta_C2 = [];
    subset = 1:20;
    
    AlphaBeta_C1.avg = permute(AlphaBeta_WithPosSac(subset,:),[3,2,1]);
    AlphaBeta_C1.time = t_range;
    AlphaBeta_C1.dimord = 'chan_time_subj';
    AlphaBeta_C1.label = {'Pos'};
    AlphaBeta_C2 = AlphaBeta_C1;
    AlphaBeta_C2.avg = permute(AlphaBeta_NoPosSac(subset,:),[3,2,1]);

    
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
    cfg.latency             = [0.2,1];
    cfg.tail             = 0;
    cfg.clustertail      = 0;
    cfg.alpha            = 0.025;
    cfg.numrandomization = 1000;
    
    cfg.design           = design;
    
    
    [stat_With]  = ft_timelockstatistics(cfg, AlphaBeta_C1,AlphaBeta_C2);
end
%% save output

save([datafigspath,'Fig6B.mat'],'stat_With','AlphaBeta_WithPosSac','AlphaBeta_NoPosSac','allnext_sac_time')
%% Figure
codepath = '';

datafigspath =  [codepath,'\data4figs\'];
figsavepath = [codepath,'\figures\'];
% load data
load( [datafigspath,'Fig6B.mat'])


% time parameter
t_range =  0:1/1000:1;

toi = [0,1];
toi_idx  = dsearchn(t_range',toi');

t_plt = toi(1):0.2:toi(2);
t_plt_ind = dsearchn(t_range(toi_idx(1):toi_idx(2))',t_plt');
sizefig = [0.1,0.1,0.7,0.6];

ylim_plt = [-1.2,1.2];
plot_save = 1;
subset = 1:20;% all channels

fig_savename = 'Fig6B_1';
for fig=33
     
f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig);clf;    
options.handle     = figure(f);
options.color_area = [219 128 193 ]./255;    % Blue theme
options.color_line = [ 186 52 148 ]./255;
options.alpha      = 0.7;
options.line_width = 2;
options.error      = 'sem';

plot_areaerrorbar(AlphaBeta_WithPosSac(subset,toi_idx(1):toi_idx(2)),options);hold on
options.color_area = [219 193 128  ]./255;    % 
options.color_line = [ 186 148 52  ]./255;

plot_areaerrorbar(AlphaBeta_NoPosSac(subset,toi_idx(1):toi_idx(2)),options);hold on
for i_c=1:length(stat_With.posclusters)
    if stat_With.posclusters(i_c).prob<0.025
        sig_tw = stat_With.time(stat_With.posclusterslabelmat==i_c);
        t_idx= dsearchn(t_range',sig_tw');
        v1 = [[t_idx(1); t_idx(end); t_idx(end); t_idx(1)], [ylim_plt(1); ylim_plt(1); ylim_plt(2); ylim_plt(2)]];
        f1 = [1 2 3 4];
        patch('Faces',f1,'Vertices',v1,'FaceColor','black','FaceAlpha',.1);
    end
end

for i_c=1:length(stat_With.negclusters)
    if stat_With.negclusters(i_c).prob<0.025
        sig_tw = stat_With.time(stat_With.negclusterslabelmat==i_c);
        t_idx= dsearchn(t_range',sig_tw');
        v1 = [[t_idx(1); t_idx(end); t_idx(end); t_idx(1)], [ylim_plt(1); ylim_plt(1); ylim_plt(2); ylim_plt(2)]];
        f1 = [1 2 3 4];
        patch('Faces',f1,'Vertices',v1,'FaceColor','black','FaceAlpha',.1);
    end
end

xline(find(t_range==0),'-')
xticks(t_plt_ind')
xticklabels(num2cell(t_plt*1000))
xlabel('Time [ms]')
ylabel('Alpha Beta amplitude')
title('Posterior Alpha Beta')
xlim([1,length(toi_idx(1):toi_idx(2))])
ylim(ylim_plt)

legend({'With','','Without',''})


set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);


if plot_save
    saveas(gcf,[figsavepath,fig_savename,'.emf'])
    saveas(gcf,[figsavepath,fig_savename,'.png'])
end
end


sizefig1 = [0.2,0.2,0.6,0.28];
fontaxis = 16;line_width = 2;
fig_savename = 'Fig6B_2';
for fig=30

    f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig1);clf;
    raincloud_plot(allnext_sac_time, 'box_on', 0, 'color',[0.5,0.5,0.5], 'alpha', 0.6,...
        'box_dodge', 0, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
    xlim([0,1])
    title('Following saccades onset time')
    xlabel('Time to Saccade onset [Sec]')
    ylim([0,1.6])
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);
    
    if plot_save
        saveas(gcf,[figsavepath,fig_savename,'.emf'])
        saveas(gcf,[figsavepath,fig_savename,'.png'])
    end
end

