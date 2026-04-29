%% Behavioral analysis for  Exp3

codepath = '';

datapath = [codepath,'\datafiles\'];
datafigspath =  [codepath,'\data4figs\'];
figsavepath = [codepath,'\figures\'];

savename = [datapath,'Beh_Exp3.mat'];
load(savename)

%%

for s = 1:length(Mem)
    %% HR hit rate  & FA False alarm
    hit_pr(s,1) = mean(pr_perf{s}(pr_cond{s}==1|pr_cond{s}==2)).*100;
    falram_pr(s,1) = mean(pr_perf{s}(pr_cond{s}==3|pr_cond{s}==4)==0).*100;
end 

[dpri,~] = dprime(hit_pr./100,falram_pr./100,216);

savename = [datafigspath,'Fig3A.mat'];
save(savename,'dpri','falram_pr','hit_pr','-v7.3')
%% descriptive analysis
data(1,:) = [mean(cellfun(@mean,pv_perf)'*100),mean(hit_pr),mean(dpri)];
data(2,:) = [std(cellfun(@mean,pv_perf)'*100),std(hit_pr),std(dpri)];

T = array2table(data, ...
    'VariableNames', {'InOutAcc', 'Hit','dprime'}, ...
    'RowNames', {'Mean', 'SD'});

% Display the table
disp(T);
%% Fig 3A
savename = [datafigspath,'Fig3A.mat'];
load(savename)
plot_save = 1; % to save or not figure

sizefig1 = [0,0.2,0.32,0.55];
fontaxis = 24;line_width = 2;

fig_savename1 = 'Fig3A';
for dp=1
f=figure('Name',int2str(11),'units','normalized','outerposition',sizefig1);clf;
scatter(falram_pr,hit_pr,'filled');hold on
set(findobj(gca,'type','line'),'linew',1.5)
ylim([0,100])
xlim([0,100])
line([0,100],[0,100],'Color','black','LineStyle','--','linewidth',2)
text(50,50,sprintf(' Mean d-prime %.2f',mean(dpri))) ;hold off
ylabel('Hit rate')
xlabel('False-alarm rate')
title('d prime')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',fontaxis,'FontWeight','Bold', 'LineWidth', line_width);

if plot_save
    saveas(gcf,[figsavepath,fig_savename1,'.emf'])
    saveas(gcf,[figsavepath,fig_savename1,'.png'])
end

end
%% Beh for trials with strict fixation  
datapath = [codepath,'\datafiles\'];
datafigspath =  [codepath,'\data4figs\'];
figsavepath = [codepath,'\figures\'];
subinfopath = [codepath,'\datafiles\Subjects_Exp3\ExpInfo\'];
sacsavepath = [codepath,'\datafiles\Subjects_Exp3\saccades\'];

savename = [datapath,'Beh_Exp3.mat'];
load(savename)
%%
% general N of saccades
for s = 1:20
    disp(s)
    load([sacsavepath,sprintf('TrialEyeInfo_%d_%s.mat',s,'Exp3')])
    load([sacsavepath,sprintf('Sac_PV_Info_Sub%d_%s.mat',s,'Exp3')]);
    %% general saccade frequency
    All0sac_trl = Eye_trlinfo(Eye_trlinfo(:,8)==1& Eye_trlinfo(:,5)==0,:); % 8th columns marks the final trials with clean EEG & eye data
    
    Allwithsac_trl =  Eye_trlinfo(Eye_trlinfo(:,8)==1& Eye_trlinfo(:,5)~=0,:); 
 
    per0sac(s,1) = size(All0sac_trl,1)./sum(Eye_trlinfo(:,8)==1).*100;
    num0sac(s,1) = size(All0sac_trl,1);
    hit0sac(s,1) = mean(All0sac_trl(:,2)).*100;
    hitwithsac(s,1) = mean(Allwithsac_trl(:,2)).*100;
    falram_pr(s,1) = mean(pr_perf{s}(pr_cond{s}==3|pr_cond{s}==4)==0).*100;

    
    
    cri = nanmedian(Eye_trlinfo(Eye_trlinfo(:,8)==1,7));
    More_trl = Eye_trlinfo(Eye_trlinfo(:,8)==1 & Eye_trlinfo(:,7)<=cri,:);
    Less_trl = Eye_trlinfo(Eye_trlinfo(:,8)==1 & Eye_trlinfo(:,7)>=cri,:);
    
    hitMoreEI(s,1) = mean(More_trl(:,2)).*100;
    hitLessEI(s,1) = mean(Less_trl(:,2)).*100;
    
    [dpriMore(s,1),~] = dprime(mean(More_trl(:,2)),falram_pr(s,1)./100,size(More_trl(:,2),1));
    [dpriLess(s,1),~] = dprime(mean(Less_trl(:,2)),falram_pr(s,1)./100,size(Less_trl(:,2),1));
    %% get trials with no saccades or only saccades less than 1 dva during 4 sec
    sec4sac = Sac_trials.trialinfo(Sac_trials.trialinfo(:,6)> 1.3*600,:);
    sac1dvatrl = unique(sec4sac(sec4sac(:,5) > 1,1));
    selecttrial = Eye_trlinfo(Eye_trlinfo(:,8)==1 & ~ismember(Eye_trlinfo(:,1),sac1dvatrl),:);
    notselecttrial = Eye_trlinfo(Eye_trlinfo(:,8)==1 & ismember(Eye_trlinfo(:,1),sac1dvatrl),:);
    per1dvatrial(s,1) = size(selecttrial,1)./sum(Eye_trlinfo(:,8)==1);
    
    hit1dvatrl(s,1) = mean(selecttrial(:,2)).*100;
    [dpri1dvatrl(s,1),~] = dprime(mean(selecttrial(:,2)),falram_pr(s,1)./100,size(selecttrial,1));
    [dprimore1dvatrl(s,1),~] = dprime(mean(notselecttrial(:,2)),falram_pr(s,1)./100,size(notselecttrial,1));
    sacfreq1dvatrl(s,1) = mean(selecttrial(selecttrial(:,2)==1,5));
    sacfreq1dvatrl(s,2) = mean(selecttrial(selecttrial(:,2)==0,5));

end


%% descriptive analysis
data(1,:) = [mean(hitMoreEI),mean(dpriMore),mean(hitLessEI),mean(dpriLess)];
data(2,:) = [std(hitMoreEI),std(dpriMore),std(hitLessEI),std(dpriLess)];

T = array2table(data, ...
    'VariableNames', {'MoreExploration_Hit', 'MoreExploration_dprime','LessExploration_Hit','LessExploration_dprime'}, ...
    'RowNames', {'Mean', 'SD'});

% Display the table
disp(T);

[bf10_dpriML,p_dpri] = bf.ttest(dpriMore-dpriLess);
[bf10_hitML,p_hit] = bf.ttest(hitMoreEI-hitLessEI);

1/bf10_dpriML
1/bf10_hitML

%%
savename = [datafigspath,'SuppleFig8.mat'];
save(savename,'falram_pr','hit1dvatrl','hitMoreEI','hitLessEI','-v7.3')
%% Figure

savename = [datafigspath,'SuppleFig8.mat'];
load(savename)

plot_save = 1; % to save or not figure


sizefig1 = [0,0.2,0.32,0.6];
fontaxis = 24;line_width = 2;
a = 0.8;
fig_savename1 = 'SuppleFig8A';
for dp=1
f=figure('Name',int2str(11),'units','normalized','outerposition',sizefig1);clf;
scatter(falram_pr,hit1dvatrl,'filled','MarkerFaceColor','k', 'MarkerEdgeColor','none', 'MarkerFaceAlpha',a);hold on
set(findobj(gca,'type','line'),'linew',1.5)
ylim([0,100])
xlim([0,100])
line([0,100],[0,100],'Color','black','LineStyle','--','linewidth',2)
text(45,40,sprintf(' Mean d-prime  %.2f',mean(dpri1dvatrl))) ;hold off
ylabel('Hit rate')
xlabel('False-alarm rate')
title({'d prime','(Trials with strict fixation )'}, 'Interpreter','tex')

set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',fontaxis,'FontWeight','Bold', 'LineWidth', line_width);

if plot_save
    saveas(gcf,[figsavepath,fig_savename1,'.emf'])
    saveas(gcf,[figsavepath,fig_savename1,'.png'])
end

end



%% 
cols = brewermap(3,'BrBG');   % 11-step BrBG palette
c1 = cols(1,:);                % one extreme (brown)
c2 = cols(end,:);  

sizefig1 = [0,0.2,0.32,0.6];
fontaxis = 24;line_width = 2;
a = 1;
fig_savename1 = 'SuppleFig8B';
for dp=1
f=figure('Name',int2str(11),'units','normalized','outerposition',sizefig1);clf;
scatter(falram_pr,hitMoreEI,'filled','MarkerFaceColor',c1, 'MarkerEdgeColor','none', 'MarkerFaceAlpha',a);hold on
scatter(falram_pr,hitLessEI,'filled','MarkerFaceColor',c2, 'MarkerEdgeColor','none', 'MarkerFaceAlpha',a);
set(findobj(gca,'type','line'),'linew',1.5)
ylim([0,100])
xlim([0,100])
line([0,100],[0,100],'Color','black','LineStyle','--','linewidth',2)
text(45,40,sprintf(' Mean d-prime (More Exploration) %.2f',mean(dpriMore))) ;
text(45,35,sprintf(' Mean d-prime (Less Exploration) %.2f',mean(dpriLess))) ;hold off
legend({'More','Less'})
ylabel('Hit rate')
xlabel('False-alarm rate')
title({'d prime','(Median-split Exploration Index)'}, 'Interpreter','tex')

set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',fontaxis,'FontWeight','Bold', 'LineWidth', line_width);

if plot_save
    saveas(gcf,[figsavepath,fig_savename1,'.emf'])
    saveas(gcf,[figsavepath,fig_savename1,'.png'])
end

end
