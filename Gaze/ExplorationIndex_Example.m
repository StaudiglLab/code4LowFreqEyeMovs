%% Example trials for Exploration Index
codepath = '';

datapath = [codepath,'\datafiles\'];
datafigspath =  [codepath,'\data4figs\'];
figsavepath = [codepath,'\figures\'];
%% Exp 1
load( [datapath,'PVstim_Gaze_Exp1.mat'])
% example trial 
sub = 17;trl_more = 37; trl_less = 2; % manually select after inspection, this could be set to rand number
eg_MoreExplore = gaze2D{sub}{trl_more};
eg_LessExplore = gaze2D{sub}{trl_less};
save([datafigspath,'SuppleFig4_left1.mat'],'eg_MoreExplore','eg_LessExplore')

% random selection
% select_sub = randsample(1:20,1);
% for s=select_sub
%     sacsavepath = [codepath,'\datafiles\Subjects_Exp1\saccades\'];
%     load([sacsavepath,sprintf('TrialEyeInfo_%d_%s.mat',s,'Exp1')])
%     
%     cri = nanmedian(Eye_trlinfo(Eye_trlinfo(:,8)==1,7));
%     g_data = cat(3,gaze2D{s}{:});
%     
%     moretrl_P=randsample(Eye_trlinfo(Eye_trlinfo(:,7)<=cri & Eye_trlinfo(:,8)==1,1),1);
%     lesstrl_P=randsample(Eye_trlinfo(Eye_trlinfo(:,7)>=cri & Eye_trlinfo(:,8)==1,1),1);
%     eg_MoreExplore = g_data(:,:,moretrl_P);
%     eg_LessExplore = g_data(:,:,lesstrl_P);
%     
% end
%% Exp 2
load( [datapath,'PVstim_Gaze_Exp2_el.mat']) % load( [datapath,'PVstim_Gaze_Exp2_tb.mat'])
% example trial for SuppleFig4
sub = 2;trl_more = 116; trl_less = 119; % manually select after inspection, this could be set to rand number
eg_MoreExplore = gaze2D{sub}{trl_more};
eg_LessExplore = gaze2D{sub}{trl_less};
save([datafigspath,'SuppleFig4_left2.mat'],'eg_MoreExplore','eg_LessExplore')
%% Exp 3
load( [datapath,'PVstim_Gaze_Exp3.mat'])
% example trial for SuppleFig4
sub = 16;trl_more = 56; trl_less = 17; % manually select after inspection, this could be set to rand number
eg_MoreExplore = gaze2D{sub}{trl_more};
eg_LessExplore = gaze2D{sub}{trl_less};
save([datafigspath,'SuppleFig4_left3.mat'],'eg_MoreExplore','eg_LessExplore')
%% Exp 4
load( [datapath,'PVstim_Gaze_Exp4.mat'])
% example trial for SuppleFig4
sub = 18;trl_more = 62; trl_less = 184; 
eg_MoreExplore = gaze2D{sub}{trl_more};
eg_LessExplore = gaze2D{sub}{trl_less};
save([datafigspath,'SuppleFig4_left4.mat'],'eg_MoreExplore','eg_LessExplore')

%% Figure
% plot

plot_save = 1; % to save or not figure
n_bin_x = 80;n_bin_y = 60;
viewdva_x =28;
viewdva_y =21;

sizefig1 = [0,0.2,0.7,0.45];
fontaxis = 24;line_width = 2;

fig_savename1 = 'SuppleFig4_left1';
load([datafigspath,'SuppleFig4_left1.mat'])
for Exp1=1
f=figure('Name',int2str(11),'units','normalized','outerposition',sizefig1);clf;
subplot(1,2,1);
imagesc(eg_MoreExplore);
xticks(linspace(0,n_bin_x,5))
xticklabels(num2cell(abs(linspace(0,n_bin_x,5)*viewdva_x/n_bin_x-viewdva_x/2)))
yticks(linspace(0,n_bin_y,5))
yticklabels(num2cell(abs(linspace(0,n_bin_y,5)*viewdva_y/n_bin_y-viewdva_y/2)))
xlabel('Dist to center [dva]')
ylabel('Dist to center [dva]')
caxis([0,0.01])
title(sprintf('More: EI = %0.4f',std(eg_MoreExplore(:))))
c = colorbar;
c.Label.String = 'Density';
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',fontaxis,'FontWeight','normal', 'LineWidth', line_width);
subplot(1,2,2);
imagesc(eg_LessExplore);
xticks(linspace(0,n_bin_x,5))
xticklabels(num2cell(abs(linspace(0,n_bin_x,5)*viewdva_x/n_bin_x-viewdva_x/2)))
yticks(linspace(0,n_bin_y,5))
yticklabels(num2cell(abs(linspace(0,n_bin_y,5)*viewdva_y/n_bin_y-viewdva_y/2)))
xlabel('Dist to center [dva]')
%ylabel('Dist to center [dva]')
caxis([0,0.01])
title(sprintf('Less: EI = %0.4f',std(eg_LessExplore(:))))
c = colorbar;
c.Label.String = 'Density';
colormap(flipud(cbrewer2('RdBu','div')));
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',fontaxis,'FontWeight','normal', 'LineWidth', line_width);

    if plot_save
        saveas(gcf,[figsavepath,fig_savename1,'.emf'])
        saveas(gcf,[figsavepath,fig_savename1,'.png'])
    end
    
end 
%
fig_savename1 = 'SuppleFig4_left2';
load([datafigspath,'SuppleFig4_left2.mat'])
for Exp2=1
f=figure('Name',int2str(11),'units','normalized','outerposition',sizefig1);clf;
subplot(1,2,1);
imagesc(eg_MoreExplore);
xticks(linspace(0,n_bin_x,5))
xticklabels(num2cell(abs(linspace(0,n_bin_x,5)*viewdva_x/n_bin_x-viewdva_x/2)))
yticks(linspace(0,n_bin_y,5))
yticklabels(num2cell(abs(linspace(0,n_bin_y,5)*viewdva_y/n_bin_y-viewdva_y/2)))
xlabel('Dist to center [dva]')
ylabel('Dist to center [dva]')
caxis([0,0.01])
title(sprintf('More: EI = %0.4f',std(eg_MoreExplore(:))))
c = colorbar;
c.Label.String = 'Density';
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',fontaxis,'FontWeight','normal', 'LineWidth', line_width);
subplot(1,2,2);
imagesc(eg_LessExplore);
xticks(linspace(0,n_bin_x,5))
xticklabels(num2cell(abs(linspace(0,n_bin_x,5)*viewdva_x/n_bin_x-viewdva_x/2)))
yticks(linspace(0,n_bin_y,5))
yticklabels(num2cell(abs(linspace(0,n_bin_y,5)*viewdva_y/n_bin_y-viewdva_y/2)))
xlabel('Dist to center [dva]')
%ylabel('Dist to center [dva]')
caxis([0,0.01])
title(sprintf('Less: EI = %0.4f',std(eg_LessExplore(:))))
c = colorbar;
c.Label.String = 'Density';
colormap(flipud(cbrewer2('RdBu','div')));
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',fontaxis,'FontWeight','normal', 'LineWidth', line_width);

    if plot_save
        saveas(gcf,[figsavepath,fig_savename1,'.emf'])
        saveas(gcf,[figsavepath,fig_savename1,'.png'])
    end
    
end 
%
fig_savename1 = 'SuppleFig4_left3';
load([datafigspath,'SuppleFig4_left3.mat'])
for Exp3=1
f=figure('Name',int2str(11),'units','normalized','outerposition',sizefig1);clf;
subplot(1,2,1);
imagesc(eg_MoreExplore);
xticks(linspace(0,n_bin_x,5))
xticklabels(num2cell(abs(linspace(0,n_bin_x,5)*viewdva_x/n_bin_x-viewdva_x/2)))
yticks(linspace(0,n_bin_y,5))
yticklabels(num2cell(abs(linspace(0,n_bin_y,5)*viewdva_y/n_bin_y-viewdva_y/2)))
xlabel('Dist to center [dva]')
ylabel('Dist to center [dva]')
caxis([0,0.01])
title(sprintf('More: EI = %0.4f',std(eg_MoreExplore(:))))
c = colorbar;
c.Label.String = 'Density';
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',fontaxis,'FontWeight','normal', 'LineWidth', line_width);
subplot(1,2,2);
imagesc(eg_LessExplore);
xticks(linspace(0,n_bin_x,5))
xticklabels(num2cell(abs(linspace(0,n_bin_x,5)*viewdva_x/n_bin_x-viewdva_x/2)))
yticks(linspace(0,n_bin_y,5))
yticklabels(num2cell(abs(linspace(0,n_bin_y,5)*viewdva_y/n_bin_y-viewdva_y/2)))
xlabel('Dist to center [dva]')
%ylabel('Dist to center [dva]')
caxis([0,0.01])
title(sprintf('Less: EI = %0.4f',std(eg_LessExplore(:))))
c = colorbar;
c.Label.String = 'Density';
colormap(flipud(cbrewer2('RdBu','div')));
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',fontaxis,'FontWeight','normal', 'LineWidth', line_width);

    if plot_save
        saveas(gcf,[figsavepath,fig_savename1,'.emf'])
        saveas(gcf,[figsavepath,fig_savename1,'.png'])
    end
    
end 
%
fig_savename1 = 'SuppleFig4_left4';
load([datafigspath,'SuppleFig4_left4.mat'])
for Exp4=1
f=figure('Name',int2str(11),'units','normalized','outerposition',sizefig1);clf;
subplot(1,2,1);
imagesc(eg_MoreExplore);
xticks(linspace(0,n_bin_x,5))
xticklabels(num2cell(abs(linspace(0,n_bin_x,5)*viewdva_x/n_bin_x-viewdva_x/2)))
yticks(linspace(0,n_bin_y,5))
yticklabels(num2cell(abs(linspace(0,n_bin_y,5)*viewdva_y/n_bin_y-viewdva_y/2)))
xlabel('Dist to center [dva]')
ylabel('Dist to center [dva]')
caxis([0,0.01])
title(sprintf('More: EI = %0.4f',std(eg_MoreExplore(:))))
c = colorbar;
c.Label.String = 'Density';
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',fontaxis,'FontWeight','normal', 'LineWidth', line_width);
subplot(1,2,2);
imagesc(eg_LessExplore);
xticks(linspace(0,n_bin_x,5))
xticklabels(num2cell(abs(linspace(0,n_bin_x,5)*viewdva_x/n_bin_x-viewdva_x/2)))
yticks(linspace(0,n_bin_y,5))
yticklabels(num2cell(abs(linspace(0,n_bin_y,5)*viewdva_y/n_bin_y-viewdva_y/2)))
xlabel('Dist to center [dva]')
%ylabel('Dist to center [dva]')
caxis([0,0.01])
title(sprintf('Less: EI = %0.4f',std(eg_LessExplore(:))))
c = colorbar;
c.Label.String = 'Density';
colormap(flipud(cbrewer2('RdBu','div')));
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',fontaxis,'FontWeight','normal', 'LineWidth', line_width);

    if plot_save
        saveas(gcf,[figsavepath,fig_savename1,'.emf'])
        saveas(gcf,[figsavepath,fig_savename1,'.png'])
    end
    
end 
%

% Group level data
% Exp 1
sizefig1 = [0,0.2,1,0.45];
fontaxis = 24;line_width = 2;

fig_savename1 = 'SuppleFig4_right11';
load([datafigspath,'SuppleFig4_right11.mat'])
for MvLExplr=1

f=figure('Name',int2str(11),'units','normalized','outerposition',sizefig1);clf;
subplot(1,3,1);
imagesc(nanmean(Gaze_MoreExplr,3));
xticks(linspace(0,n_bin_x,5))
xticklabels(num2cell(abs(linspace(0,n_bin_x,5)*viewdva_x/n_bin_x-viewdva_x/2)))
yticks(linspace(0,n_bin_y,5))
yticklabels(num2cell(abs(linspace(0,n_bin_y,5)*viewdva_y/n_bin_y-viewdva_y/2)))
xlabel('Dist to center [dva]')
ylabel('Dist to center [dva]')
caxis([0,0.002])
title('More')
c = colorbar;
c.Label.String = 'Density';
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',fontaxis,'FontWeight','normal', 'LineWidth', line_width);
subplot(1,3,2);
imagesc(nanmean(Gaze_LessExplr,3));
xticks(linspace(0,n_bin_x,5))
xticklabels(num2cell(abs(linspace(0,n_bin_x,5)*viewdva_x/n_bin_x-viewdva_x/2)))
yticks(linspace(0,n_bin_y,5))
yticklabels(num2cell(abs(linspace(0,n_bin_y,5)*viewdva_y/n_bin_y-viewdva_y/2)))
xlabel('Dist to center [dva]')
%ylabel('Dist to center [dva]')
caxis([0,0.002])
title('Less')
c = colorbar;
c.Label.String = 'Density';
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',fontaxis,'FontWeight','normal', 'LineWidth', line_width);
subplot(1,3,3);
Cluster=zeros(size(Gaze_MoreExplr,[1,2]));
Cluster(squeeze(stat_EI.mask)) = 1;
[B,L,N,A] = bwboundaries(Cluster); 
imagesc(squeeze(stat_EI.stat));hold on
for k=1:length(B)
     boundary = B{k};
    plot(boundary(:,2), boundary(:,1), 'LineWidth',3,'color','k');
end
hold off
xticks(linspace(0,n_bin_x,5))
xticklabels(num2cell(abs(linspace(0,n_bin_x,5)*viewdva_x/n_bin_x-viewdva_x/2)))
yticks(linspace(0,n_bin_y,5))
yticklabels(num2cell(abs(linspace(0,n_bin_y,5)*viewdva_y/n_bin_y-viewdva_y/2)))
caxis([-2.5,2.5])
title('More vs Less')
xlabel('Dist to center [dva]')
%ylabel('Dist to center [dva]')
colormap(flipud(cbrewer2('RdBu','div')));
c = colorbar;
c.Label.String = 't value';
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',fontaxis,'FontWeight','normal', 'LineWidth', line_width);

    if plot_save
        saveas(gcf,[figsavepath,fig_savename1,'.emf'])
        saveas(gcf,[figsavepath,fig_savename1,'.png'])
    end
    
end 

% Exp 2
fig_savename1 = 'SuppleFig4_right21';
load([datafigspath,'SuppleFig4_right21.mat'])
for MvLExplr=1

f=figure('Name',int2str(11),'units','normalized','outerposition',sizefig1);clf;
subplot(1,3,1);
imagesc(nanmean(Gaze_MoreExplr,3));
xticks(linspace(0,n_bin_x,5))
xticklabels(num2cell(abs(linspace(0,n_bin_x,5)*viewdva_x/n_bin_x-viewdva_x/2)))
yticks(linspace(0,n_bin_y,5))
yticklabels(num2cell(abs(linspace(0,n_bin_y,5)*viewdva_y/n_bin_y-viewdva_y/2)))
xlabel('Dist to center [dva]')
ylabel('Dist to center [dva]')
caxis([0,0.002])
title('More')
c = colorbar;
c.Label.String = 'Density';
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',fontaxis,'FontWeight','normal', 'LineWidth', line_width);
subplot(1,3,2);
imagesc(nanmean(Gaze_LessExplr,3));
xticks(linspace(0,n_bin_x,5))
xticklabels(num2cell(abs(linspace(0,n_bin_x,5)*viewdva_x/n_bin_x-viewdva_x/2)))
yticks(linspace(0,n_bin_y,5))
yticklabels(num2cell(abs(linspace(0,n_bin_y,5)*viewdva_y/n_bin_y-viewdva_y/2)))
xlabel('Dist to center [dva]')
%ylabel('Dist to center [dva]')
caxis([0,0.002])
title('Less')
c = colorbar;
c.Label.String = 'Density';
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',fontaxis,'FontWeight','normal', 'LineWidth', line_width);
subplot(1,3,3);
imagesc(squeeze(stat_EI.stat));
xticks(linspace(0,n_bin_x,5))
xticklabels(num2cell(abs(linspace(0,n_bin_x,5)*viewdva_x/n_bin_x-viewdva_x/2)))
yticks(linspace(0,n_bin_y,5))
yticklabels(num2cell(abs(linspace(0,n_bin_y,5)*viewdva_y/n_bin_y-viewdva_y/2)))
caxis([-2.5,2.5])
title('More vs Less')
xlabel('Dist to center [dva]')
%ylabel('Dist to center [dva]')
colormap(flipud(cbrewer2('RdBu','div')));
c = colorbar;
c.Label.String = 't value';
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',fontaxis,'FontWeight','normal', 'LineWidth', line_width);

    if plot_save
        saveas(gcf,[figsavepath,fig_savename1,'.emf'])
        saveas(gcf,[figsavepath,fig_savename1,'.png'])
    end
    
end 


% Exp 3
fig_savename1 = 'SuppleFig4_right31';
load([datafigspath,'SuppleFig4_right31.mat'])
for MvLExplr=1

f=figure('Name',int2str(11),'units','normalized','outerposition',sizefig1);clf;
subplot(1,3,1);
imagesc(nanmean(Gaze_MoreExplr,3));
xticks(linspace(0,n_bin_x,5))
xticklabels(num2cell(abs(linspace(0,n_bin_x,5)*viewdva_x/n_bin_x-viewdva_x/2)))
yticks(linspace(0,n_bin_y,5))
yticklabels(num2cell(abs(linspace(0,n_bin_y,5)*viewdva_y/n_bin_y-viewdva_y/2)))
xlabel('Dist to center [dva]')
ylabel('Dist to center [dva]')
caxis([0,0.002])
title('More')
c = colorbar;
c.Label.String = 'Density';
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',fontaxis,'FontWeight','normal', 'LineWidth', line_width);
subplot(1,3,2);
imagesc(nanmean(Gaze_LessExplr,3));
xticks(linspace(0,n_bin_x,5))
xticklabels(num2cell(abs(linspace(0,n_bin_x,5)*viewdva_x/n_bin_x-viewdva_x/2)))
yticks(linspace(0,n_bin_y,5))
yticklabels(num2cell(abs(linspace(0,n_bin_y,5)*viewdva_y/n_bin_y-viewdva_y/2)))
xlabel('Dist to center [dva]')
%ylabel('Dist to center [dva]')
caxis([0,0.002])
title('Less')
c = colorbar;
c.Label.String = 'Density';
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',fontaxis,'FontWeight','normal', 'LineWidth', line_width);
subplot(1,3,3);
Cluster=zeros(size(Gaze_MoreExplr,[1,2]));
Cluster(squeeze(stat_EI.mask)) = 1;
[B,L,N,A] = bwboundaries(Cluster); 
imagesc(squeeze(stat_EI.stat));hold on
for k=1:length(B)
     boundary = B{k};
    plot(boundary(:,2), boundary(:,1), 'LineWidth',3,'color','k');
end
hold off
xticks(linspace(0,n_bin_x,5))
xticklabels(num2cell(abs(linspace(0,n_bin_x,5)*viewdva_x/n_bin_x-viewdva_x/2)))
yticks(linspace(0,n_bin_y,5))
yticklabels(num2cell(abs(linspace(0,n_bin_y,5)*viewdva_y/n_bin_y-viewdva_y/2)))
caxis([-2.5,2.5])
title('More vs Less')
xlabel('Dist to center [dva]')
%ylabel('Dist to center [dva]')
colormap(flipud(cbrewer2('RdBu','div')));
c = colorbar;
c.Label.String = 't value';
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',fontaxis,'FontWeight','normal', 'LineWidth', line_width);

    if plot_save
        saveas(gcf,[figsavepath,fig_savename1,'.emf'])
        saveas(gcf,[figsavepath,fig_savename1,'.png'])
    end
    
end 

% Exp 4
fig_savename1 = 'SuppleFig4_right41';
load([datafigspath,'SuppleFig4_right41.mat'])
for MvLExplr=1

f=figure('Name',int2str(11),'units','normalized','outerposition',sizefig1);clf;
subplot(1,3,1);
imagesc(nanmean(Gaze_MoreExplr,3));
xticks(linspace(0,n_bin_x,5))
xticklabels(num2cell(abs(linspace(0,n_bin_x,5)*viewdva_x/n_bin_x-viewdva_x/2)))
yticks(linspace(0,n_bin_y,5))
yticklabels(num2cell(abs(linspace(0,n_bin_y,5)*viewdva_y/n_bin_y-viewdva_y/2)))
xlabel('Dist to center [dva]')
ylabel('Dist to center [dva]')
caxis([0,0.002])
title('More')
c = colorbar;
c.Label.String = 'Density';
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',fontaxis,'FontWeight','normal', 'LineWidth', line_width);
subplot(1,3,2);
imagesc(nanmean(Gaze_LessExplr,3));
xticks(linspace(0,n_bin_x,5))
xticklabels(num2cell(abs(linspace(0,n_bin_x,5)*viewdva_x/n_bin_x-viewdva_x/2)))
yticks(linspace(0,n_bin_y,5))
yticklabels(num2cell(abs(linspace(0,n_bin_y,5)*viewdva_y/n_bin_y-viewdva_y/2)))
xlabel('Dist to center [dva]')
%ylabel('Dist to center [dva]')
caxis([0,0.002])
title('Less')
c = colorbar;
c.Label.String = 'Density';
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',fontaxis,'FontWeight','normal', 'LineWidth', line_width);
subplot(1,3,3);
Cluster=zeros(size(Gaze_MoreExplr,[1,2]));
Cluster(squeeze(stat_EI.mask)) = 1;
[B,L,N,A] = bwboundaries(Cluster); 
imagesc(squeeze(stat_EI.stat));hold on
for k=1:length(B)
     boundary = B{k};
    plot(boundary(:,2), boundary(:,1), 'LineWidth',3,'color','k');
end
hold off
xticks(linspace(0,n_bin_x,5))
xticklabels(num2cell(abs(linspace(0,n_bin_x,5)*viewdva_x/n_bin_x-viewdva_x/2)))
yticks(linspace(0,n_bin_y,5))
yticklabels(num2cell(abs(linspace(0,n_bin_y,5)*viewdva_y/n_bin_y-viewdva_y/2)))
caxis([-2.5,2.5])
title('More vs Less')
xlabel('Dist to center [dva]')
%ylabel('Dist to center [dva]')
colormap(flipud(cbrewer2('RdBu','div')));
c = colorbar;
c.Label.String = 't value';
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',fontaxis,'FontWeight','normal', 'LineWidth', line_width);

    if plot_save
        saveas(gcf,[figsavepath,fig_savename1,'.emf'])
        saveas(gcf,[figsavepath,fig_savename1,'.png'])
    end
    
end 


% saccade 
sizefig = [0.1,0.1,0.2,0.6];
cb = repmat([0.5,0.5,0.5],1,1);
fig_savename = 'SuppleFig4_right12';
load([datafigspath,'SuppleFig4_right12.mat'])
for fig=21
    %% 
    f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig);clf;
    h=daboxplot([n_sac_Explr_PV(:,[1:2])],[repmat(1,size(n_sac_Explr_PV,1),1)],...
        'scatter' ,2,...
        'color',cb,...
        'boxalpha', 0.5,...
        'jitter',1,...
        'linkline',0,...
        'withinlines',1)
    
    ylabel('N Sec') 
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
    
   
    if plot_save
    saveas(gcf,[figsavepath,fig_savename,'.emf'])
    saveas(gcf,[figsavepath,fig_savename,'.png'])
    end
end 

% Exp 2
fig_savename = 'SuppleFig4_right22';
load([datafigspath,'SuppleFig4_right22.mat'])
for fig=21
    %% 
    f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig);clf;
    h=daboxplot([n_sac_Explr_PV(:,[1:2])],[repmat(1,size(n_sac_Explr_PV,1),1)],...
        'scatter' ,2,...
        'color',cb,...
        'boxalpha', 0.5,...
        'jitter',1,...
        'linkline',0,...
        'withinlines',1)
    
    ylabel('N Sec') 
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
    
   
    if plot_save
    saveas(gcf,[figsavepath,fig_savename,'.emf'])
    saveas(gcf,[figsavepath,fig_savename,'.png'])
    end
end 

% Exp 3
fig_savename = 'SuppleFig4_right32';
load([datafigspath,'SuppleFig4_right32.mat'])
for fig=21
    %% 
    f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig);clf;
    h=daboxplot([n_sac_Explr_PV(:,[1:2])],[repmat(1,size(n_sac_Explr_PV,1),1)],...
        'scatter' ,2,...
        'color',cb,...
        'boxalpha', 0.5,...
        'jitter',1,...
        'linkline',0,...
        'withinlines',1)
    
    ylabel('N Sec') 
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
    
   
    if plot_save
    saveas(gcf,[figsavepath,fig_savename,'.emf'])
    saveas(gcf,[figsavepath,fig_savename,'.png'])
    end
end 


% Exp 4
fig_savename = 'SuppleFig4_right42';
load([datafigspath,'SuppleFig4_right42.mat'])
for fig=21
    %% 
    f=figure('Name',int2str(fig),'units','normalized','outerposition',sizefig);clf;
    h=daboxplot([n_sac_Explr_PV(:,[1:2])],[repmat(1,size(n_sac_Explr_PV,1),1)],...
        'scatter' ,2,...
        'color',cb,...
        'boxalpha', 0.5,...
        'jitter',1,...
        'linkline',0,...
        'withinlines',1)
    
    ylabel('N Sec') 
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
    
   
    if plot_save
    saveas(gcf,[figsavepath,fig_savename,'.emf'])
    saveas(gcf,[figsavepath,fig_savename,'.png'])
    end
end 