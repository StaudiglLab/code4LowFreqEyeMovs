%% Behavioral analysis for Exp4
codepath = '';

datapath = [codepath,'\datafiles\'];
datafigspath =  [codepath,'\data4figs\'];
figsavepath = [codepath,'\figures\'];

savename = [datapath,'Beh_Exp4.mat'];
load(savename)
%%
% memory separated by viewing condition
for s = 1:length(pv_viewcond)
    
    Memory_PV_VC(s,1) = mean(Mem{s}(pv_viewcond{s}==1)).*100;
    Memory_PV_VC(s,2) = mean(Mem{s}(pv_viewcond{s}==2)).*100;
    Memory_PV_VC(s,3) = mean(Mem{s}(pv_viewcond{s}==3)).*100;
    
    
end

for s = 1:length(Mem)
    %% HR hit rate  & FA False alarm
    hit_pr(s,1) = mean(pr_perf{s}(pr_cond{s}==1|pr_cond{s}==2)).*100;
    falram_pr(s,1) = mean(pr_perf{s}(pr_cond{s}==3|pr_cond{s}==4)==0).*100;
end 

[dpri,~] = dprime(hit_pr./100,falram_pr./100,324);

savename = [datafigspath,'Fig4A.mat'];
save(savename,'dpri','falram_pr','hit_pr','Memory_PV_VC','-v7.3')
%% descriptive analysis
data(1,:) = [mean(cellfun(@mean,pv_perf)'*100),mean(hit_pr),mean(dpri)];
data(2,:) = [std(cellfun(@mean,pv_perf)'*100),std(hit_pr),std(dpri)];

T = array2table(data, ...
    'VariableNames', {'InOutAcc', 'Hit','dprime'}, ...
    'RowNames', {'Mean', 'SD'});

T2 = array2table([mean(Memory_PV_VC);std(Memory_PV_VC)], ...
    'VariableNames', {'MemV1', 'MemV2','MemV3'}, ...
    'RowNames', {'Mean', 'SD'});
% Display the table
disp(T);

disp(T2);

%% Fig 4A
savename = [datafigspath,'Fig4A.mat'];
load(savename)
plot_save = 1; % to save or not figure

sizefig1 = [0,0.2,0.32,0.55];
fontaxis = 24;line_width = 2;

fig_savename1 = 'Fig4A1';
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


sizefig = [0.2,0.2,0.22,0.4];fontaxis = 24;line_width = 2;
cb = repmat([0.5,0.5,0.5],3,1);
fig_savename = 'Fig4A2';
for fig=21
    %% 
    f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig);clf;
    h=daboxplot([Memory_PV_VC(:,[1:3])],[repmat(1,size(Memory_PV_VC,1),1)],...
        'scatter' ,2,...
        'color',cb,...
        'boxalpha', 0.5,...
        'jitter',1,...
        'linkline',0,...
        'withinlines',1,...
        'xtlabels',{'1DVA','5DVA','21DVA'})
    
    ylabel('Hit rate [%]')
    xlabel('Viewing Condition')
    ylim([0,100])
   
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
    
   
    
    if plot_save
    saveas(gcf,[figsavepath,fig_savename,'.emf'])
    saveas(gcf,[figsavepath,fig_savename,'.png'])
    end
end 

