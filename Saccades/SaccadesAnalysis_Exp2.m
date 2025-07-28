%% Saccade analysis for Experiment 1
% run this after completing saccades detection
codepath = '';

datapath = [codepath,'\datafiles\'];
datafigspath =  [codepath,'\data4figs\'];
figsavepath = [codepath,'\figures\'];
subinfopath = [codepath,'\datafiles\Subjects_Exp2\ExpInfo\'];
sacsavepath = [codepath,'\datafiles\Subjects_Exp2\saccades\'];
%%
% general N of saccades
for s = 1:4
    disp(s)
    load([sacsavepath,sprintf('TrialEyeInfo_%d_%s.mat',s,'Exp2')])
    
    %% Sep by Rem vs Forg
    Rem_trl = Eye_trlinfo(Eye_trlinfo(:,2)==1,:);
    For_trl = Eye_trlinfo(Eye_trlinfo(:,2)==0,:);
    n_sac_RF_bs(s,1) =  nanmean(Rem_trl(:,4));
    n_sac_RF_bs(s,2) =  nanmean(For_trl(:,4));
    
    n_sac_RF_PV(s,1) =  nanmean(Rem_trl(:,5))./4;
    n_sac_RF_PV(s,2) =  nanmean(For_trl(:,5))./4;
    clearvars Rem_trl For_trl
    %% Sey by N saccades 
    cri = nanmedian(Eye_trlinfo(:,5));
    More_trl = Eye_trlinfo( Eye_trlinfo(:,5)>=cri,:);
    Less_trl = Eye_trlinfo( Eye_trlinfo(:,5)<=cri,:);
    
    n_sac_Nsac_bs(s,1) =  nanmean(More_trl(:,4));
    n_sac_Nsac_bs(s,2) =  nanmean(Less_trl(:,4));
    
    n_sac_Nsac_PV(s,1) =  nanmean(More_trl(:,5))./4;
    n_sac_Nsac_PV(s,2) =  nanmean(Less_trl(:,5))./4;
    
    clearvars More_trl Less_trl 
    
    %% Sey by Exploration index
    cri = nanmedian(Eye_trlinfo(:,7));
    More_trl = Eye_trlinfo( Eye_trlinfo(:,7)<=cri,:);
    Less_trl = Eye_trlinfo( Eye_trlinfo(:,7)>=cri,:);
    
    n_sac_Explr_PV(s,1) =  nanmean(More_trl(:,5))./4;
    n_sac_Explr_PV(s,2) =  nanmean(Less_trl(:,5))./4;
    
    clearvars More_trl Less_trl


end
save([datafigspath,'SuppleFig2_right22.mat'],'n_sac_Explr_PV')
%% descriptive analysis
data(1,:) = mean(n_sac_RF_PV);
data(2,:) = std(n_sac_RF_PV);

T = array2table(data, ...
    'VariableNames', {'Rem', 'Forg'}, ...
    'RowNames', {'Mean', 'SD'});

% Display the table
disp(T);

%% Saccades metrics

edge_ang = -180:10:180;
for s = 1:4
    disp(s)
    
    % load trial info
    load([sacsavepath,sprintf('TrialEyeInfo_%d_%s.mat',s,'Exp2')])
    % load saccade info
    load([sacsavepath,sprintf('Sac_PV_Info_Sub%d_%s.mat',s,'Exp2')]);

    switch  s % P01 & 2 with 600Hz eye, 3 & 4 with 1000
        case {1,2}
            eye_Hz = 600;
        case {3,4}
            eye_Hz = 1000;
    end

     %% Amplitdue
    alltrl = Sac_trials.trialinfo;
    
    sac_amplitude_Exp2(s,1) = nanmean(alltrl(:,5));
    %% Duration 
    sac_duration_Exp2(s,1) = nanmean(alltrl(:,7)-alltrl(:,6))/eye_Hz*1000;
    %% Direction
    bin_data = histcounts(alltrl(:,4),edge_ang);
    ms_dir_Exp2(:,s) = bin_data./sum(bin_data); clearvars bin_data
    
  
end

%% Save output
save([datafigspath,'Fig2A_right.mat'],'n_sac_RF_PV')
save([datafigspath,'SuppleFig1_Exp2.mat'],'sac_amplitude_Exp2','sac_duration_Exp2','ms_dir_Exp2')
% save for suuplementary fig2
save([datafigspath,'SuppleFig2_right22.mat'],'n_sac_Explr_PV')
%%

%% Time resolved saccade frequency

bin_t = 0.20; % 200 ms per bin
edge = (-1:bin_t:4);

for s = 1:4
    disp(s)
    
   % load trial info
    load([sacsavepath,sprintf('TrialEyeInfo_%d_%s.mat',s,'Exp2')])
    % load saccade info
    load([sacsavepath,sprintf('Sac_PV_Info_Sub%d_%s.mat',s,'Exp2')]);


    switch  s % P01 & 2 with 600Hz eye, 3 & 4 with 1000
        case {1,2}
            eye_Hz = 600;
        case {3,4}
            eye_Hz = 1000;
    end
    t_ET  = -1.3:1/eye_Hz:4;
   
    trl_info = Eye_trlinfo;

    cri = nanmedian(Eye_trlinfo(:,5));
    C1_sac = Sac_trials.trialinfo(ismember(Sac_trials.trialinfo(:,1),trl_info(trl_info(:,5)>=cri,1)),:);
    C2_sac = Sac_trials.trialinfo(ismember(Sac_trials.trialinfo(:,1),trl_info(trl_info(:,5)<=cri,1)),:);
    
    cri = nanmedian(Eye_trlinfo(:,7));
    C1_explr = Sac_trials.trialinfo(ismember(Sac_trials.trialinfo(:,1),trl_info(trl_info(:,7)<=cri,1)),:);
    C2_explr = Sac_trials.trialinfo(ismember(Sac_trials.trialinfo(:,1),trl_info(trl_info(:,7)>=cri,1)),:);
    
    clearvars cri trl_info
    %% general distribution around pic onset
    t_on = Sac_trials.trialinfo(:,6)+t_ET(1)*eye_Hz;% correct the time lock to pic onset
    onset_picon(:,s) = histcounts(t_on,edge.*eye_Hz)./length(unique(Sac_trials.trialinfo(:,1)))*(1/bin_t); % div by trial num and convert to N_sac / sec;
    
    % hit miss
    t_on_hit = Sac_trials.trialinfo(Sac_trials.trialinfo(:,2)==1,6)+t_ET(1)*eye_Hz;% correct the time lock to pic onset
    onset_picon_hit(:,s) = histcounts(t_on_hit,edge.*eye_Hz)./length(unique(Sac_trials.trialinfo(Sac_trials.trialinfo(:,2)==1,1)))*(1/bin_t); % div by trial num and convert to N_sac / sec;
    
    t_on_mis = Sac_trials.trialinfo(Sac_trials.trialinfo(:,2)==0,6)+t_ET(1)*eye_Hz;% correct the time lock to pic onset
    onset_picon_mis(:,s) = histcounts(t_on_mis,edge.*eye_Hz)./length(unique(Sac_trials.trialinfo(Sac_trials.trialinfo(:,2)==0,1)))*(1/bin_t); % div by trial num and convert to N_sac / sec;
    
    

end

save([datafigspath,'Fig2A_left.mat'],'onset_picon_hit','onset_picon_mis')

%% Figure

codepath = '';

datafigspath =  [codepath,'\data4figs\'];
figsavepath = [codepath,'\figures\'];
% load saccade frequency across time
load( [datafigspath,'Fig2A_left.mat'])
load( [datafigspath,'Fig2A_right.mat'])

plot_save = 1;
bin_t = 0.20; % 200 ms per bin
edge = (-1:bin_t:4);%-bin_t/2;
t_range = edge(2:end);
t_plt = -1:1:4;
t_plt_ind = dsearchn(t_range',t_plt');
ylim_plt = [0,4];

sizefig = [0.2,0.2,0.6,0.7];fontaxis = 16;line_width = 2;
fig_savename = 'Fig2A_left';
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



sizefig = [0.2,0.2,0.3,0.8];
fig_savename = 'Fig2A_right';
for fig=221 % violin plot
    %% violin plot
    f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig);clf;
   
    Violin({n_sac_RF_PV(:,1)}, 1,...
        'HalfViolin','left',...% left, full
        'ShowData',false,...
        'QuartileStyle','shadow',...
        'ShowNotches', false,...
        'ShowMean', true,...
        'ShowMedian', true,...
        'ShowWhiskers',false); hold on
    Violin({n_sac_RF_PV(:,2)}, 2,...
        'HalfViolin','right',...% left, full
        'ShowData',false,...
        'QuartileStyle','shadow',...
        'ShowNotches', false,...
        'ShowMean', true,...
        'ShowMedian', true,...
        'ShowWhiskers',false); hold on
    
    
    %xlim([0,3])
    x_adjust = -0.1; % to left
    tight_adjust = 8; % how tight of scatter
    indv_cond1 = ones(size(n_sac_RF_PV(:,1))).*(1+(rand(size(n_sac_RF_PV(:,1)))-0.5)/tight_adjust-x_adjust);
    indv_cond2 = ones(size(n_sac_RF_PV(:,2))).*(2+(rand(size(n_sac_RF_PV(:,2)))-0.5)/tight_adjust+x_adjust);
  
    scatter(indv_cond1,n_sac_RF_PV(:,1),'k','filled');
    scatter(indv_cond2,n_sac_RF_PV(:,2),'k','filled');
    line([indv_cond1,indv_cond2]',[n_sac_RF_PV(:,1),n_sac_RF_PV(:,2)]','color',[0.5,0.5,0.5]);
    hold off
    ylabel('N of saccades [N/Sec]')
    %xlabel('Task')
    xticks([1 2 ])
    xticklabels({'Rem','Forg'})
    ylim([1,2.5])
    title('Saccade Frequency PV')
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);

    %% whether save fig
    
    if plot_save
    saveas(gcf,[figsavepath,fig_savename,'.emf'])
    saveas(gcf,[figsavepath,fig_savename,'.png'])
    end
end


