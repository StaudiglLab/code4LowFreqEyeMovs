%% Gaze density analysis for Exp1
% run this after extracting gaze density array (GazeDensity_Exp1)
codepath = '';

datapath = [codepath,'\datafiles\'];
datafigspath =  [codepath,'\data4figs\'];
figsavepath = [codepath,'\figures\'];
subinfopath = [codepath,'\datafiles\Subjects_Exp1\ExpInfo\'];
sacsavepath = [codepath,'\datafiles\Subjects_Exp1\saccades\'];
% load behavioral data, memorability metric
load( [datapath,'Beh_Exp1.mat'])
load( [datapath,'Stim_Memorability_Exp1.mat'])


%% Memorability analysis 
% preparation
T = [];
for s = 1:20
    disp(s)
    load([sacsavepath,sprintf('TrialEyeInfo_%d_%s.mat',s,'Exp1')])
    Eye_trlinfo = Eye_trlinfo(Eye_trlinfo(:,8)==1,:);
    % add imag idx % memorability
    Eye_trlinfo(:,9)=Img_id_enc{s}(Eye_trlinfo(:,1));
    Eye_trlinfo(:,10)=MB_stim{s}(Eye_trlinfo(:,1));
    %% 
    Ts = [];
    Ts.subject = categorical(rep_num (s,size(Eye_trlinfo,1),1)');
    Ts.image   = categorical(Eye_trlinfo(:,9));
    Ts.y       = double(Eye_trlinfo(:,2));                         % ensure column
    Ts.memorability = double(Eye_trlinfo(:,10));
    Ts.nSac = double(Eye_trlinfo(:,5));
    %Ts.gzVar = double(Eye_trlinfo(:,7));
    T = cat(1,T,struct2table(Ts));
    clearvars Eye_trlinfo
end

% Center/scale memorability (helps GLME convergence & interpretability)
T.mem_z = (T.memorability - mean(T.memorability,'omitnan')) ...
           / std(T.memorability,[],'omitnan');
T.nSac_z = (T.nSac - mean(T.nSac,'omitnan')) ...
           / std(T.nSac,[],'omitnan');

       

%% Model
% Mixed-effects logistic regression
m1 = fitglme(T, 'y ~ 1 + mem_z + (1|subject) + (1|image)', ...
    'Distribution','Binomial','Link','logit', 'FitMethod','Laplace');

m2 = fitglme(T, 'y ~ 1 + mem_z + nSac_z+(1|subject) + (1|image)', ...
    'Distribution','Binomial','Link','logit', 'FitMethod','Laplace');

% model comparison 
cmp = compare(m1,m2);   % likelihood-ratio test
disp('Model comparison m1 vs m2 (random slopes):'); disp(cmp);

% Mixed-effects linear regression
ms = fitlme(T, 'nSac_z ~ 1 + mem_z + (1|subject) + (1|image)');           
       

save([datafigspath,'SuppleFig2.mat'],'m1','m2','ms','T')
%% Figure
codepath = '';

datafigspath =  [codepath,'\data4figs\'];
figsavepath = [codepath,'\figures\'];
% load saccade frequency across time
load( [datafigspath,'SuppleFig2.mat'])

%% Effect projection
xg = linspace(min(T.mem_z), max(T.mem_z), 100)';

refSub = T.subject(1);
refImg = T.image(1);
Tg = table(xg, repmat(refSub,size(xg)), repmat(refImg,size(xg)), ...
    'VariableNames', {'mem_z','subject','image'});
% Predict using fixed effects only (marginal over random effects)
[yhat, yCI] = predict(m1, Tg, 'Conditional', false);


%%

plot_save = 1; % to save or not figure

sizefig1 = [0,0.2,0.3,0.45];
fontaxis = 24;line_width = 2;

fig_savename1 = 'SuppleFig2_left';

for fig=50
    f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig1);clf;
    plot(xg, yhat, 'LineWidth',2); hold on;
    fill([xg; flipud(xg)], [yCI(:,1); flipud(yCI(:,2))], [0.8 0.8 0.8], ...
        'EdgeColor','none','FaceAlpha',0.5);
    
    %Jittered scatter of raw data
%     x_jit = T.mem_z;
%     y_jit = T.y + (rand(size(T.y))-0.5)*0.05;  % add small vertical jitter
%     scatter(x_jit, y_jit, 10, [0.6 0.6 0.6], 'filled', 'MarkerFaceAlpha',0.2);
    xlabel('Memorability (z)'); ylabel('Predicted P(Remembered)');
    xlim([-3,3])
    title('Mixed-effects Logistic: Predicted Probability vs Memorability');
    legend({'Estimate','95% CI'}, 'Location','southeast'); box off;
    
    %% whether save fig
    
    if plot_save
        saveas(gcf,[figsavepath,fig_savename1,'.emf'])
        saveas(gcf,[figsavepath,fig_savename1,'.png'])
    end
    
end

sizefig1 = [0,0.2,0.3,0.45];
fontaxis = 24;line_width = 2;

fig_savename1 = 'SuppleFig2_right';

for fig=50
f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig1);clf;hold on
xg = linspace(min(T.mem_z), max(T.mem_z), 200)';
Tg = table(xg, repmat(T.subject(1),size(xg)), repmat( T.image(1),size(xg)), 'VariableNames', {'mem_z','subject','image'});
[yhat, yCI] = predict(ms, Tg, 'Conditional', false);
% Per-image means for a clean scatter
G = groupsummary(T,'image','mean',{'mem_z','nSac_z'});
scatter(G.mean_mem_z, G.mean_nSac_z, 18, 'k', 'filled', 'MarkerFaceAlpha',0.35);
plot(xg, yhat, 'b-', 'LineWidth', 2);
fill([xg; flipud(xg)], [yCI(:,1); flipud(yCI(:,2))], [0.7 0.8 1], ...
     'EdgeColor','none','FaceAlpha',0.35);
xlabel('Memorability (z)'); ylabel('Saccades (z)');
title('Saccades ~ Memorability (mixed-effects)');
legend({'Per-image means','LME estimate','95% CI'}, 'Location','southeast');
grid on; box off;  
%% whether save fig
    
    if plot_save
        saveas(gcf,[figsavepath,fig_savename1,'.emf'])
        saveas(gcf,[figsavepath,fig_savename1,'.png'])
    end
end