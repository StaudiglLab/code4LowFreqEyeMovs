%% Time Frequency analyis for Exp 2

codepath = '';
addpath(genpath( [codepath,'\subfunctions\']))


datapath = [codepath,'\datafiles\'];
datafigspath =  [codepath,'\data4figs\'];
figsavepath = [codepath,'\figures\'];
subinfopath = [codepath,'\datafiles\Subjects_Exp2\ExpInfo\'];
sacsavepath = [codepath,'\datafiles\Subjects_Exp2\saccades\'];
eegsavepath = [codepath,'\datafiles\Subjects_Exp2\EEG\'];
TFRsavepath = [codepath,'\datafiles\Subjects_Exp2\TFR\'];
% behavioral results
savename = [datapath,'Beh_Exp2.mat'];
load(savename)
% 

%% Separate conditions
% global ft_default
% ft_default.showcallinfo = 'no';

rng(2025);
load([TFRsavepath,'TFR_PVstim_allTrials_allsub_Exp2.mat'])
nperm = 200; % N of permutations
% average channels for subject 
% baseline correction
bl =  [-1,-0.5];
for s = 1:4
    %
    TFR_data = tfralltrl{s};
    spike_data = spike_mark_sub{s};
    rej_data = trial_rej_sub{s} ;
    
    TFR_data = MarkTFRNans_ft(TFR_data,spike_data); % replace the marked artifacts as NaNs
    chanFSlabel = TFR_data.elecinfo;
    TFR_data = rmfield(TFR_data,'elecinfo');

    % load trial info
    load([sacsavepath,sprintf('TrialEyeInfo_%d_%s.mat',s,'Exp2')])

    for ichan = 1:length(TFR_data.label)

        cfg              = [];
        cfg.baseline     = bl;
        cfg.baselinetype = 'db';
        cfgtrl           = [];
        cfgtrl.channel      = TFR_data.label(ichan);
        cfgtrl.trials       = find(~ismember(TFR_data.trialinfo(:,1),rej_data{ichan})&TFR_data.trialinfo(:,2)==1&ismember(TFR_data.trialinfo(:,1),Eye_trlinfo(:,1)));
        TFR_sub_rem{ichan} = ft_freqbaseline(cfg, ft_freqdescriptives(cfgtrl,TFR_data));
        TFR_sub_rem{ichan}.label      = {'Posterior'};
        
        cfg              = [];
        cfg.baseline     = bl;
        cfg.baselinetype = 'db';
        cfgtrl           = [];
        cfgtrl.channel      = TFR_data.label(ichan);
        cfgtrl.trials       = find(~ismember(TFR_data.trialinfo(:,1),rej_data{ichan})&TFR_data.trialinfo(:,2)==0&ismember(TFR_data.trialinfo(:,1),Eye_trlinfo(:,1)));
        TFR_sub_for{ichan} = ft_freqbaseline(cfg, ft_freqdescriptives(cfgtrl,TFR_data));
        TFR_sub_for{ichan}.label      = {'Posterior'};
        

        
        
        %% sep by saccades
        sac_cri =  nanmedian(Eye_trlinfo(~ismember(Eye_trlinfo(:,1),rej_data{ichan}),5));
        sac_info = nan(length(TFR_data.trialinfo(:,1)),1);
        sac_info(Eye_trlinfo(:,1)) = Eye_trlinfo(:,5);

        
        cfg              = [];
        cfg.baseline     = bl;
        cfg.baselinetype = 'db';
        cfgtrl           = [];
        cfgtrl.channel      = TFR_data.label(ichan);
        cfgtrl.trials       = find(~ismember(TFR_data.trialinfo(:,1),rej_data{ichan})&sac_info>=sac_cri);
        TFR_sub_Msac{ichan} = ft_freqbaseline(cfg, ft_freqdescriptives(cfgtrl,TFR_data));
        TFR_sub_Msac{ichan}.label      = {'Posterior'};
        
        
        cfg              = [];
        cfg.baseline     = bl;
        cfg.baselinetype = 'db';
        cfgtrl           = [];
        cfgtrl.channel      = TFR_data.label(ichan);
        cfgtrl.trials       = find(~ismember(TFR_data.trialinfo(:,1),rej_data{ichan})&sac_info<=sac_cri);
        TFR_sub_Lsac{ichan} = ft_freqbaseline(cfg, ft_freqdescriptives(cfgtrl,TFR_data));
        TFR_sub_Lsac{ichan}.label      = {'Posterior'};


        
    end

    %% average across sensors 
    TFR_RealavgRem{s} = ft_freqgrandaverage([],TFR_sub_rem{:});
    TFR_RealavgFor{s} = ft_freqgrandaverage([],TFR_sub_for{:});
    TFR_RealavgMsac{s} = ft_freqgrandaverage([],TFR_sub_Msac{:});
    TFR_RealavgLsac{s} = ft_freqgrandaverage([],TFR_sub_Lsac{:});



    clearvars TFR_sub_*
    %% permutation
    for iperm = 1:nperm
        if rem(iperm,10)==0
            fprintf('Subject: %d:  %d/100 complete  \n',s,iperm/nperm*100);
        end

        for ichan = 1:length(TFR_data.label)
            validtrl = find(~ismember(TFR_data.trialinfo(:,1),rej_data{ichan})&ismember(TFR_data.trialinfo(:,1),Eye_trlinfo(:,1)));
            behidx = TFR_data.trialinfo(validtrl,2);
            TFR_data.trialinfo(validtrl,2) = behidx(randperm(length(behidx)));

            cfg              = [];
            cfg.baseline     = bl;
            cfg.baselinetype = 'db';
            cfgtrl           = [];
            cfgtrl.channel      = TFR_data.label(ichan);
            cfgtrl.trials       = find(~ismember(TFR_data.trialinfo(:,1),rej_data{ichan})&TFR_data.trialinfo(:,2)==1&ismember(TFR_data.trialinfo(:,1),Eye_trlinfo(:,1)));
            TFR_sub_rem{ichan} = ft_freqbaseline(cfg, ft_freqdescriptives(cfgtrl,TFR_data));
            TFR_sub_rem{ichan}.label      = {'Posterior'};

            cfg              = [];
            cfg.baseline     = bl;
            cfg.baselinetype = 'db';
            cfgtrl           = [];
            cfgtrl.channel      = TFR_data.label(ichan);
            cfgtrl.trials       = find(~ismember(TFR_data.trialinfo(:,1),rej_data{ichan})&TFR_data.trialinfo(:,2)==0&ismember(TFR_data.trialinfo(:,1),Eye_trlinfo(:,1)));
            TFR_sub_for{ichan} = ft_freqbaseline(cfg, ft_freqdescriptives(cfgtrl,TFR_data));
            TFR_sub_for{ichan}.label      = {'Posterior'};




            %% sep by saccades
            sac_cri =  nanmedian(Eye_trlinfo(~ismember(Eye_trlinfo(:,1),rej_data{ichan}),5));
            sac_info = nan(length(TFR_data.trialinfo(:,1)),1);
            sac_info(Eye_trlinfo(:,1)) = Eye_trlinfo(:,5);
            sacidx = sac_info(validtrl);
            sac_info(validtrl) = sacidx(randperm(length(sacidx)));


            cfg              = [];
            cfg.baseline     = bl;
            cfg.baselinetype = 'db';
            cfgtrl           = [];
            cfgtrl.channel      = TFR_data.label(ichan);
            cfgtrl.trials       = find(~ismember(TFR_data.trialinfo(:,1),rej_data{ichan})&sac_info>=sac_cri);
            TFR_sub_Msac{ichan} = ft_freqbaseline(cfg, ft_freqdescriptives(cfgtrl,TFR_data));
            TFR_sub_Msac{ichan}.label      = {'Posterior'};


            cfg              = [];
            cfg.baseline     = bl;
            cfg.baselinetype = 'db';
            cfgtrl           = [];
            cfgtrl.channel      = TFR_data.label(ichan);
            cfgtrl.trials       = find(~ismember(TFR_data.trialinfo(:,1),rej_data{ichan})&sac_info<=sac_cri);
            TFR_sub_Lsac{ichan} = ft_freqbaseline(cfg, ft_freqdescriptives(cfgtrl,TFR_data));
            TFR_sub_Lsac{ichan}.label      = {'Posterior'};



        end


        %% average across sensors
        TFR_permRem = ft_freqgrandaverage([],TFR_sub_rem{:});
        TFR_permFor = ft_freqgrandaverage([],TFR_sub_for{:});
        TFR_permMsac = ft_freqgrandaverage([],TFR_sub_Msac{:});
        TFR_permLsac = ft_freqgrandaverage([],TFR_sub_Lsac{:});

        clearvars TFR_sub_*

        % 
        TFR_allpermRem(iperm,:,:) = TFR_permRem.powspctrm;
        TFR_allpermFor(iperm,:,:) = TFR_permFor.powspctrm;
        TFR_allpermMsac(iperm,:,:) = TFR_permMsac.powspctrm;
        TFR_allpermLsac(iperm,:,:) = TFR_permLsac.powspctrm;

        clearvars TFR_perm*


    end

    % construct structure
    TFR_subpermRem{s} = rmfield(TFR_RealavgRem{s},'cfg');
    TFR_subpermRem{s}.dimord = 'rpt_freq_time';
    TFR_subpermRem{s}.powspctrm = TFR_allpermRem;


    TFR_subpermFor{s} = TFR_subpermRem{s};
    TFR_subpermFor{s}.powspctrm = TFR_allpermFor;
    TFR_subpermMsac{s} = TFR_subpermRem{s};
    TFR_subpermMsac{s}.powspctrm = TFR_allpermMsac;
    TFR_subpermLsac{s} = TFR_subpermRem{s};
    TFR_subpermLsac{s}.powspctrm = TFR_allpermLsac;

clearvars TFR_all*

end

% 
% global ft_default
% ft_default.showcallinfo = 'yes';
%%
save([datapath,'TFR_Exp2_Sub.mat'],'TFR_subperm*','TFR_Real*','-v7.3');
%%

load([datapath,'TFR_Exp2_Sub.mat']);

load([datafigspath,'Fig2C_left.mat'])
load([datafigspath,'Fig2D_left.mat'])

for s = 1:4
    cfg = [];
    cfg.operation = 'subtract';
    cfg.parameter = 'powspctrm';
    SME = ft_math(cfg,TFR_subpermRem{s},TFR_subpermFor{s});
    SMEReal = ft_math(cfg,TFR_RealavgRem{s},TFR_RealavgFor{s});
    cfg=[];
    cfg.latency = [0,4];
    SME= ft_selectdata(cfg,SME);
    SMEReal= ft_selectdata(cfg,SMEReal);

    mask2 = squeeze(stat_mem.mask);              
    avgSME(:,s) = squeeze(mean(SME.powspctrm(:,mask2), 2));
    realSME(s,1)=mean(SMEReal.powspctrm(stat_mem.mask));

  
    %%
    cfg = [];
    cfg.operation = 'subtract';
    cfg.parameter = 'powspctrm';
    MLsac = ft_math(cfg,TFR_subpermMsac{s},TFR_subpermLsac{s});
    MLsacReal = ft_math(cfg,TFR_RealavgMsac{s},TFR_RealavgLsac{s});
    cfg=[];
    cfg.latency = [0,4];
    MLsac= ft_selectdata(cfg,MLsac);
    MLsacReal= ft_selectdata(cfg,MLsacReal);

    mask2 = squeeze(stat_sac.mask);            
    avgMLsac(:,s) = squeeze(mean(MLsac.powspctrm(:,mask2), 2));
    realMLsac(s,1)=mean(MLsacReal.powspctrm(stat_sac.mask));


end


% group level distribution
[~,~,~,tt] = ttest(realSME);
trealSME = tt.tstat;
[~,~,~,tt] = ttest(realMLsac);
trealMLsac = tt.tstat;
for i=1:size(avgMLsac,1)
    [~,~,~,tt] = ttest(avgSME(i,:));
    tsurrSME(i,1)=tt.tstat;
    [~,~,~,tt] = ttest(avgMLsac(i,:));
    tsurrMLsac(i,1)=tt.tstat;


end

%% 
save([datafigspath,'SuppleFig7.mat'],'tsurrSME','tsurrMLsac','trealSME','trealMLsac','avgSME','realSME','avgMLsac','realMLsac')

%% Figure

datapath = [codepath,'\datafiles\'];
datafigspath =  [codepath,'\data4figs\'];
figsavepath = [codepath,'\figures\'];
load([datafigspath,'SuppleFig7.mat'])

plot_save = 1;
fig_savename1 = 'SuppleFig7C';
fig_savename2 = 'SuppleFig7D';
sizefig1 = [0.2,0.2,1,0.4];

for fig = 11
f=figure('Name',int2str(11),'units','normalized','outerposition',sizefig1);clf;
for s = 1:4
subplot(1,4,s);
histogram(avgSME(:,s),20,'Normalization','probability');
hold on;
xline(realSME(s), 'r', 'LineWidth', 2);
ylim([0,0.15])
xlabel('Power difference (dB)')
ylabel('Probability')
title(sprintf('P%02d',s))
hold off;

end
sgtitle('Remembered vs Forgotten', 'FontSize', 24, 'FontWeight', 'bold');
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);

if plot_save
    saveas(gcf,[figsavepath,fig_savename1,'.emf'])
    saveas(gcf,[figsavepath,fig_savename1,'.png'])
end
end


for fig = 12
f=figure('Name',int2str(12),'units','normalized','outerposition',sizefig1);clf;
for s = 1:4
subplot(1,4,s);
histogram(avgMLsac(:,s),20,'Normalization','probability');
hold on;
xline(realMLsac(s), 'r', 'LineWidth', 2);
ylim([0,0.15])
xlabel('Power difference (dB)')
ylabel('Probability')
title(sprintf('P%02d',s))
hold off;

end
sgtitle('More vs Fewer Saccades', 'FontSize', 24, 'FontWeight', 'bold');
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);

if plot_save
    saveas(gcf,[figsavepath,fig_savename2,'.emf'])
    saveas(gcf,[figsavepath,fig_savename2,'.png'])
end
end

%% group level

fig_savename1 = 'SuppleFig7A';
fig_savename2 = 'SuppleFig7B';
plot_save = 1;
sizefig1 = [0.2,0.2,0.4,0.4];
for fig = 11
f=figure('Name',int2str(11),'units','normalized','outerposition',sizefig1);clf;

histogram(tsurrSME,20,'Normalization','probability','FaceColor',[0.2,0.2,0.2]);
hold on;
xline(trealSME, 'r', 'LineWidth', 2);
ylim([0,0.35])
xlabel('t value')
ylabel('Probability')
hold off;

set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
title('Remembered vs Forgotten', 'FontSize', 24, 'FontWeight', 'bold');
if plot_save
    saveas(gcf,[figsavepath,fig_savename1,'.emf'])
    saveas(gcf,[figsavepath,fig_savename1,'.png'])
end
end
pSME = (sum(abs(tsurrSME) >= abs(trealSME)) + 1) / (length(tsurrSME) + 1)
%
for fig = 12
f=figure('Name',int2str(12),'units','normalized','outerposition',sizefig1);clf;

histogram(tsurrMLsac,20,'Normalization','probability','FaceColor',[0.2,0.2,0.2]);
hold on;
xline(trealMLsac, 'r', 'LineWidth', 2);
ylim([0,0.3])
xlabel('t value')
ylabel('Probability')
hold off;

set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
title('More vs Fewer Saccades', 'FontSize', 24, 'FontWeight', 'bold');
if plot_save
    saveas(gcf,[figsavepath,fig_savename2,'.emf'])
    saveas(gcf,[figsavepath,fig_savename2,'.png'])
end
end

pMLsac = (sum(abs(tsurrMLsac) >= abs(trealMLsac)) + 1) / (length(tsurrMLsac) + 1)