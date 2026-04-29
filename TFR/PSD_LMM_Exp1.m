%% Linear mixed effect model on PSD for Exp1
codepath = '';
addpath(genpath( [codepath,'\subfunctions\']))


datapath = [codepath,'/datafiles/'];
datafigspath =  [codepath,'/data4figs/'];
figsavepath = [codepath,'/figures/'];
subinfopath = [codepath,'/datafiles/Subjects_Exp1/ExpInfo/'];
sacsavepath = [codepath,'/datafiles/Subjects_Exp1/saccades/'];
eegsavepath = [codepath,'/datafiles/Subjects_Exp1/EEG/'];
PSDsavepath = [codepath,'/datafiles/Subjects_Exp1/PSD/'];
% behavioral results
savename = [datapath,'Beh_Exp1.mat'];
load(savename)

%% Compute PSD for a defined time window and Freq
toi = [1,3.499];% 1 - 3.5 sec post stim
toibc = [-0.999,-0.5]; % baseline
foi = 2:2:40;% same as TFR

for s = 1:length(Mem)
    
    disp(s)
    %% load EEG file of each subject
    savename = [eegsavepath,sprintf('Sub_eegICA_PV_%d_%s.mat',s,'Exp1')];
    load(savename);
     
    eeg = eeg_final;
     
    %% add beh
    eeg.trialinfo(:,2) = Mem{s}(eeg.trialinfo(:,1));
    eeg.trialinfo(:,3) = Confi{s}(eeg.trialinfo(:,1));

    %% do time freq analysis
    cfg              = [];
    cfg.latency      = toi;
    eegstim = ft_selectdata(cfg,eeg);
    cfg.latency      = toibc;
    eegbaseline = ft_selectdata(cfg,eeg);
    
    cfg              = [];
    cfg.output       = 'pow';
    cfg.method       = 'mtmfft';
    cfg.taper        = 'hanning';
    %cfg.pad         = 'nextpow2';
    cfg.foi          = foi;                         
    cfg.keeptrials = 'yes';
    PSD_alltrl = ft_freqanalysis(cfg, eegstim);
    cfg.toi          = toibc;
    PSD_alltrlbc = ft_freqanalysis(cfg, eegbaseline);
    
    %% save data

    savename = [PSDsavepath,sprintf('PSD_PVstim_allTrials_sub%d_%s.mat',s,'Exp1')];
    save(savename,'PSD_alltrl','PSD_alltrlbc','-v7.3')

end
%% construct LMM
load([datapath,'layout.mat']) % load layout
coi = {'P*','O*'};% 19 Posterior electrodes
foi = [10,20];

T = [];PowChan = [];
for s = 1:length(Mem)
    disp(s)
    savename = [PSDsavepath,sprintf('PSD_PVstim_allTrials_sub%d_%s.mat',s,'Exp1')];
    load(savename)
    load([sacsavepath,sprintf('TrialEyeInfo_%d_%s.mat',s,'Exp1')])
    %% Average across freq & chan
    cfg = [];
    cfg.trials = ismember(PSD_alltrl.trialinfo(:,1),Eye_trlinfo(Eye_trlinfo(:,8)==1,1));
    cfg.channel = coi;
    cfg.frequency = foi;
    powpos = ft_selectdata(cfg,PSD_alltrl);
    powposbl = ft_selectdata(cfg,PSD_alltrlbc);
    avgpow = mean(mean(powpos.powspctrm,3),2);
    avgpowbl = mean(mean(powposbl.powspctrm,3),2);
    n_sac = Eye_trlinfo(dsearchn(Eye_trlinfo(:,1),powpos.trialinfo(:,1)),5);
    n_sacbl = Eye_trlinfo(dsearchn(Eye_trlinfo(:,1),powpos.trialinfo(:,1)),4);

    %% 
    Ts = [];
    Ts.subject = categorical(rep_num (s,size(avgpow,1),1)');
    Ts.pow   = avgpow;
    Ts.powbl       = avgpowbl;                         % ensure column
    Ts.nSac   = double(n_sac);
    Ts.nSacbl = double(n_sacbl);
    Ts.Memory = powpos.trialinfo(:,2);
    Ts.Conf = powpos.trialinfo(:,3);

    T = cat(1,T,struct2table(Ts));
    
    
    %% Average across freq for each channel
    cfg = [];
    cfg.trials = ismember(PSD_alltrl.trialinfo(:,1),Eye_trlinfo(Eye_trlinfo(:,8)==1,1));
    cfg.frequency = foi;
    cfg.avgoverfreq = 'yes';
    powabavg = ft_selectdata(cfg,PSD_alltrl);% average over frequency of interest
    % order the label as in lay
    [tf, idx] = ismember(lay.label, powabavg.label);
    newlabel = lay.label(tf);
    avgpowchan = powabavg.powspctrm(:,idx(tf));

    PowChan = cat(1,PowChan,avgpowchan);
    clearvars avgpowchan powabavg tf idx Eye_trlinfo


end


T.powbc = T.pow - T.powbl;
% Center/scale
T.pow_z = (T.pow - mean(T.pow,'omitnan')) ...
           / std(T.pow,[],'omitnan');
T.nSac_z = (T.nSac - mean(T.nSac,'omitnan')) ...
           / std(T.nSac,[],'omitnan');
T.powbc_z = (T.powbc - mean(T.powbc,'omitnan')) ...
           / std(T.powbc,[],'omitnan');

PowChan_z = zscore(PowChan,[],1);


save([datapath,'PSD_LMM_Exp1.mat'],'T','PowChan_z');
%% 
load([datapath,'PSD_LMM_Exp1.mat'])

% Mixed-effects logistic regression
m1 = fitglme(T, 'Memory ~ 1 + pow_z + (1|subject) ', ...
    'Distribution','Binomial','Link','logit', 'FitMethod','Laplace');
m1a = fitglme(T, 'Memory ~ 1 + nSac_z + (1|subject) ', ...
    'Distribution','Binomial','Link','logit', 'FitMethod','Laplace');
mfull = fitglme(T, 'Memory ~ 1 + nSac_z +pow_z+pow_z*nSac_z+ (1|subject) ', ...
    'Distribution','Binomial','Link','logit', 'FitMethod','Laplace');
% m1a has a better AIC and BIC
dAIC = m1a.ModelCriterion.AIC-m1.ModelCriterion.AIC;
dBIC = m1a.ModelCriterion.BIC-m1.ModelCriterion.BIC;

dAICf_s = mfull.ModelCriterion.AIC-m1a.ModelCriterion.AIC;
dBICf_s = mfull.ModelCriterion.BIC-m1a.ModelCriterion.BIC;

dAICf_p = mfull.ModelCriterion.AIC-m1.ModelCriterion.AIC;
dBICf_p = mfull.ModelCriterion.BIC-m1.ModelCriterion.BIC;
%%
% also compare both to null model
m0 = fitglme(T, 'Memory ~ 1 +  (1|subject) ', ...
    'Distribution','Binomial','Link','logit', 'FitMethod','Laplace');

cmppow = compare(m1,m0);   % likelihood-ratio test
cmpsac = compare(m1a,m0); 

% compare both to full model
cmpFu2pow = compare(mfull,m1);  
cmpFu2sac = compare(mfull,m1a);  

%% pow

% Mixed-effects linear regression
ms = fitlme(T, 'pow_z ~ 1 + nSac_z + (1|subject) ');     
msa = fitlme(T, 'pow_z ~ 1  + Memory+(1|subject) ');    
msfull2 = fitlme(T, 'pow_z ~ 1 + nSac_z + Memory+Memory*nSac_z+(1|subject) ');    
cmp2 = compare(ms,msfull2);  
disp('Model comparison ms vs ms2 (random slopes):'); disp(cmp2);
cmp3 = compare(msa,msfull2);   
disp('Model comparison msa vs ms2 (random slopes):'); disp(cmp3);

dAIC = ms.ModelCriterion.AIC-msa.ModelCriterion.AIC;
dBIC = ms.ModelCriterion.BIC-msa.ModelCriterion.BIC;
fprintf('Sac only model vs Pow only model: dAIC = %.02f; dBIC= %.02f\n',dAIC,dBIC);
fprintf('Full model vs Sac only model: dAIC = %.02f; dBIC= %.02f\n',msfull2.ModelCriterion.AIC-ms.ModelCriterion.AIC,msfull2.ModelCriterion.BIC-ms.ModelCriterion.BIC);
fprintf('Full model vs Memory only model: dAIC = %.02f; dBIC= %.02f\n',msfull2.ModelCriterion.AIC-msa.ModelCriterion.AIC,msfull2.ModelCriterion.BIC-msa.ModelCriterion.BIC);

%%
save([datafigspath,'SuppleFig4.mat'],'ms','msa','msfull2');
save([datafigspath,'SuppleFig5.mat'],'m1','m1a','mfull');
%% plot result
codepath = '';

datafigspath =  [codepath,'\data4figs\'];
figsavepath = [codepath,'\figures\'];
% load saccade frequency across time
load( [datafigspath,'SuppleFig5.mat'])

% Full model

plot_save = 1; % to save or not figure
fig_savename1 = 'SuppleFig5A_left';

sizefig1 = [0,0.2,0.3,0.45];
for fig=50
    f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig1);clf;hold on
beta = fixedEffects(mfull);
CovB = mfull.CoefficientCovariance;

% check coefficient order 
% disp(mfull.Coefficients.Name)

x = linspace(-3, 3, 200)';

% Model: intercept + pow_z + nSac_z + pow_z:nSac_z
% hold pow_z = 0
X = [ones(size(x)), zeros(size(x)), x, zeros(size(x))];

eta = X * beta;
se_eta = sqrt(diag(X * CovB * X'));

% 95% CI
eta_lo = eta - 1.96 * se_eta;
eta_hi = eta + 1.96 * se_eta;

% transform to probability
p    = 1 ./ (1 + exp(-eta));
p_lo = 1 ./ (1 + exp(-eta_lo));
p_hi = 1 ./ (1 + exp(-eta_hi));
fill([x; flipud(x)], [p_lo; flipud(p_hi)], [0.85 0.85 0.85], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(x, p, 'k', 'LineWidth', 2);
xlabel('nSac_z');
ylabel('Predicted P(Memory = 1)');
title('Full model: effect of saccade count');
hold off;
%% whether save fig
    
    if plot_save
        saveas(gcf,[figsavepath,fig_savename1,'.emf'])
        saveas(gcf,[figsavepath,fig_savename1,'.png'])
    end
    
end
%


fig_savename1 = 'SuppleFig5A_right';

sizefig1 = [0,0.2,0.3,0.45];
for fig=51
    f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig1);clf;hold on

    beta = fixedEffects(mfull);
CovB = mfull.CoefficientCovariance;

x = linspace(-3, 3, 200)';

% Model: intercept + pow_z + nSac_z + pow_z:nSac_z
% hold nsac_z = 0
X = [ones(size(x)), x, zeros(size(x)),  zeros(size(x))];

eta = X * beta;
se_eta = sqrt(diag(X * CovB * X'));

% 95% CI 
eta_lo = eta - 1.96 * se_eta;
eta_hi = eta + 1.96 * se_eta;

% transform to probability
p    = 1 ./ (1 + exp(-eta));
p_lo = 1 ./ (1 + exp(-eta_lo));
p_hi = 1 ./ (1 + exp(-eta_hi));
fill([x; flipud(x)], [p_lo; flipud(p_hi)], [0.85 0.85 0.85], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(x, p, 'k', 'LineWidth', 2);
xlabel('pow_z');
ylabel('Predicted P(Memory = 1)');
title('Full model: effect of ab PSD');
hold off;
%% whether save fig
    
    if plot_save
        saveas(gcf,[figsavepath,fig_savename1,'.emf'])
        saveas(gcf,[figsavepath,fig_savename1,'.png'])
    end
    
end
%

%% plot simple model
% sole saccade model

fig_savename1 = 'SuppleFig5B';

sizefig1 = [0,0.2,0.3,0.45];
for fig=52
    f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig1);clf;hold on
beta = fixedEffects(m1a);
CovB = m1a.CoefficientCovariance;

% check coefficient order 
% disp(m1a.Coefficients.Name)

x = linspace(-3, 3, 200)';

% Model: intercept +  nSac_z 

X = [ones(size(x)),  x];

eta = X * beta;
se_eta = sqrt(diag(X * CovB * X'));

% 95% CI 
eta_lo = eta - 1.96 * se_eta;
eta_hi = eta + 1.96 * se_eta;

% transform to probability
p    = 1 ./ (1 + exp(-eta));
p_lo = 1 ./ (1 + exp(-eta_lo));
p_hi = 1 ./ (1 + exp(-eta_hi));
fill([x; flipud(x)], [p_lo; flipud(p_hi)], [0.85 0.85 0.85], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(x, p, 'k', 'LineWidth', 2);
xlabel('nSac_z');
ylabel('Predicted P(Memory = 1)');
title('Simple model: effect of saccade count');
hold off;
%% whether save fig
    
    if plot_save
        saveas(gcf,[figsavepath,fig_savename1,'.emf'])
        saveas(gcf,[figsavepath,fig_savename1,'.png'])
    end
    
end
%

% sole power model
fig_savename1 = 'SuppleFig5C';
sizefig1 = [0,0.2,0.3,0.45];
for fig=53
    f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig1);clf;hold on
beta = fixedEffects(m1);
CovB = m1.CoefficientCovariance;

% check coefficient order 
% disp(m1a.Coefficients.Name)

x = linspace(-3, 3, 200)';

% Model: intercept +  nSac_z 

X = [ones(size(x)),  x];

eta = X * beta;
se_eta = sqrt(diag(X * CovB * X'));

% 95% CI 
eta_lo = eta - 1.96 * se_eta;
eta_hi = eta + 1.96 * se_eta;

% transform to probability
p    = 1 ./ (1 + exp(-eta));
p_lo = 1 ./ (1 + exp(-eta_lo));
p_hi = 1 ./ (1 + exp(-eta_hi));
fill([x; flipud(x)], [p_lo; flipud(p_hi)], [0.85 0.85 0.85], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(x, p, 'k', 'LineWidth', 2);
xlabel('pow_z');
ylabel('Predicted P(Memory = 1)');
title('Simple model: effect of ab PSD');
hold off;
%% whether save fig
    
    if plot_save
        saveas(gcf,[figsavepath,fig_savename1,'.emf'])
        saveas(gcf,[figsavepath,fig_savename1,'.png'])
    end
    
end
%

%% Power ~ Sac*Mem
load( [datafigspath,'SuppleFig4.mat'])
fig_savename1 = 'SuppleFig4A_right';

sizefig1 = [0,0.2,0.28,0.45];
for fig=50
    f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig1);clf;hold on
beta = fixedEffects(msfull2);
CovB = msfull2.CoefficientCovariance;

% Check coefficient order 
%disp(msfull2.Coefficients.Name)

X = [1 0 0 0;
     1 1 0 0];

yhat = X * beta;
se = sqrt(diag(X * CovB * X'));

ci_lo = yhat - 1.96 * se;
ci_hi = yhat + 1.96 * se;


errorbar([1 2], yhat, yhat-ci_lo, ci_hi-yhat, 'o', 'LineWidth', 1.5);
xlim([0.5 2.5]);
xticks([1 2]);
xticklabels({'Forgotten','Remembered'});
xlabel('Memory');
ylabel('Predicted PSD (z)');
title('Full model: Effect of memory');
ylim([-0.5,0.5])
hold off;

%% whether save fig
    
    if plot_save
        saveas(gcf,[figsavepath,fig_savename1,'.emf'])
        saveas(gcf,[figsavepath,fig_savename1,'.png'])
    end
    
end
%

fig_savename1 = 'SuppleFig4A_left';

sizefig1 = [0,0.2,0.3,0.45];
for fig=51
    f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig1);clf;hold on

    beta = fixedEffects(msfull2);
CovB = msfull2.CoefficientCovariance;

x = linspace(-3, 3, 200)';

% Model: intercept + Memory + nSac_z + Memory:nSac_z
% hold nsac_z = 0
X = [ones(size(x)), zeros(size(x)), x,  zeros(size(x))];

eta = X * beta;
se_eta = sqrt(diag(X * CovB * X'));

% 95% CI 
eta_lo = eta - 1.96 * se_eta;
eta_hi = eta + 1.96 * se_eta;

% transform to probability
p    = 1 ./ (1 + exp(-eta));
p_lo = 1 ./ (1 + exp(-eta_lo));
p_hi = 1 ./ (1 + exp(-eta_hi));
fill([x; flipud(x)], [p_lo; flipud(p_hi)], [0.85 0.85 0.85], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(x, p, 'k', 'LineWidth', 2);
xlabel('Saccade');
ylabel('Predicted PSD (z)');
title('Effect of saccade');
title('Full model: effect of saccade');
hold off;
%% whether save fig
    
    if plot_save
        saveas(gcf,[figsavepath,fig_savename1,'.emf'])
        saveas(gcf,[figsavepath,fig_savename1,'.png'])
    end
    
end


%% single model 
% Power ~ Sac


fig_savename1 = 'SuppleFig4B';

sizefig1 = [0,0.2,0.3,0.45];
for fig=52
    f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig1);clf;hold on
beta = fixedEffects(ms);
CovB = ms.CoefficientCovariance;

% check coefficient order 
% disp(ms.Coefficients.Name)

x = linspace(-3, 3, 200)';

% Model: intercept +  nSac_z 

X = [ones(size(x)),  x];

eta = X * beta;
se_eta = sqrt(diag(X * CovB * X'));

% 95% CI
eta_lo = eta - 1.96 * se_eta;
eta_hi = eta + 1.96 * se_eta;

% transform to probability
p    = 1 ./ (1 + exp(-eta));
p_lo = 1 ./ (1 + exp(-eta_lo));
p_hi = 1 ./ (1 + exp(-eta_hi));
fill([x; flipud(x)], [p_lo; flipud(p_hi)], [0.85 0.85 0.85], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(x, p, 'k', 'LineWidth', 2);
xlabel('nSac_z');
ylabel('Predicted PSD(z)');
title('Simple model: effect of saccade count');

hold off;
%% whether save fig
    
    if plot_save
        saveas(gcf,[figsavepath,fig_savename1,'.emf'])
        saveas(gcf,[figsavepath,fig_savename1,'.png'])
    end
    
end
%

% sole memory model
fig_savename1 = 'SuppleFig4C';
sizefig1 = [0,0.2,0.28,0.45];
for fig=53
    f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig1);clf;hold on
beta = fixedEffects(msa);
CovB = msa.CoefficientCovariance;


% Check coefficient order 
%disp(msa.Coefficients.Name)
X = [1 0 ;
     1 1];

yhat = X * beta;
se = sqrt(diag(X * CovB * X'));

ci_lo = yhat - 1.96 * se;
ci_hi = yhat + 1.96 * se;

errorbar([1 2], yhat, yhat-ci_lo, ci_hi-yhat, 'o', 'LineWidth', 1.5);
xlim([0.5 2.5]);
xticks([1 2]);
xticklabels({'Forgotten','Remembered'});
xlabel('Memory');
ylabel('Predicted PSD(z)');
title('Simple model: effect of memory');
ylim([-0.5,0.5])
hold off;
%% whether save fig
    
    if plot_save
        saveas(gcf,[figsavepath,fig_savename1,'.emf'])
        saveas(gcf,[figsavepath,fig_savename1,'.png'])
    end
    
end



