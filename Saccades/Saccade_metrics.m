%% Saccade metrics across 4 experiments

codepath = '';


datafigspath =  [codepath,'\data4figs\'];
figsavepath = [codepath,'\figures\'];
load([datafigspath,'figsavepath1_Exp1.mat'])
load([datafigspath,'SuppleFig1_Exp2.mat'])
load([datafigspath,'SuppleFig1_Exp3.mat'])
load([datafigspath,'SuppleFig1_Exp4.mat'])

%%
plot_save = 1;
x= [sac_amplitude_Exp1;sac_amplitude_Exp2;sac_amplitude_Exp3;sac_amplitude_Exp4];
y= [sac_duration_Exp1;sac_duration_Exp2;sac_duration_Exp3;sac_duration_Exp4];
g = cellstr(int2str([ones(length(sac_duration_Exp1),1);ones(length(sac_duration_Exp2),1)*2;ones(length(sac_duration_Exp3),1)*3;ones(length(sac_duration_Exp4),1)*4]));

sizefig1 = [0.2,0.2,0.3,0.5];
fontaxis = 16;line_width = 2;
fig_savename = 'SuppleFig1_A';

for fig=52
    f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig1);clf;
    h = scatterhist(x,y,...
        'Group',g,...
        'Kernel','on','Location','SouthEast',...
        'Direction','out','LineStyle',{'-','-','-','-'},...
        'LineWidth',[2,2],'Marker','.','MarkerSize',[10,10]);

    xlim([0,7])
    ylim([10,45])
    title('Main Sequence')
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',fontaxis,'FontWeight','Bold', 'LineWidth', line_width);

    if plot_save
        saveas(gcf,[figsavepath,fig_savename,'.emf'])
        saveas(gcf,[figsavepath,fig_savename,'.png'])
    end
end

sizefig1 = [0.2,0.2,0.21,0.8];
fontaxis = 16;line_width = 2;
fig_savename = 'SuppleFig1_A1';

for fig=52
     f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig1);clf;
    h=daboxplot(x,g,...
        'scatter' ,0,...
        'boxalpha', 0.5,...
        'jitter',1,...
        'linkline',0,...
        'withinlines',1,...
        'xtlabels',{'1h','1c','2','3'})
   ylim([0,7])
    ylabel('Amplitude [dva]')
    xlabel('dataset')
   
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);

    if plot_save
        saveas(gcf,[figsavepath,fig_savename,'.emf'])
        saveas(gcf,[figsavepath,fig_savename,'.png'])
    end
end


sizefig1 = [0.2,0.2,0.21,0.8];
fontaxis = 16;line_width = 2;
fig_savename = 'SuppleFig1_A2';

for fig=52
     f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig1);clf;
    h=daboxplot(y,g,...
        'scatter' ,0,...
        'boxalpha', 0.5,...
        'jitter',1,...
        'linkline',0,...
        'withinlines',1,...
        'xtlabels',{'1h','1c','2','3'})
   
    ylabel('Duration [msec]')
    xlabel('dataset')
   ylim([10,45])
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);

    if plot_save
        saveas(gcf,[figsavepath,fig_savename,'.emf'])
        saveas(gcf,[figsavepath,fig_savename,'.png'])
    end
end
%% direction
edge_ang = -180:10:180;

sizefig3 = [0.2,0.2,0.8,0.9];
fontaxis = 16;line_width = 2;
fig_savename = 'SuppleFig1_B';
for fig=53
    f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig3);clf;
   
%%
subplot(2,2,1);
polarhistogram('BinEdges', deg2rad(edge_ang),'BinCounts',mean(ms_dir_Exp1,2,'omitnan'))
    title('1')
    rlim([0,0.1])

    subplot(2,2,2);
polarhistogram('BinEdges', deg2rad(edge_ang),'BinCounts',mean(ms_dir_Exp2,2,'omitnan'))
    title('2')
    rlim([0,0.1])

    subplot(2,2,3);
polarhistogram('BinEdges', deg2rad(edge_ang),'BinCounts',mean(ms_dir_Exp3,2,'omitnan'))
    title('3')
    rlim([0,0.1])

    subplot(2,2,4);
polarhistogram('BinEdges', deg2rad(edge_ang),'BinCounts',mean(ms_dir_Exp4,2,'omitnan'))
    title('4')
    rlim([0,0.1])


    if plot_save
        saveas(gcf,[figsavepath,fig_savename,'.emf'])
        saveas(gcf,[figsavepath,fig_savename,'.png'])
    end
end



