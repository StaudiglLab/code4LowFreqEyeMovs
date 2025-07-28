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
%%
