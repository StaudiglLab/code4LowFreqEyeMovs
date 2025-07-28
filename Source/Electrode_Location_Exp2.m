%% Load and plot the electrodes loc of Exp 2

codepath = '';
datafigspath =  [codepath,'\data4figs\'];
figsavepath = [codepath,'\figures\'];
% 
load([datafigspath,'Fig2B.mat'])
%%
plot_save = 1;

% Create color template
cmap_subject = brewermap(length(elec_info),'Set1');
gray = [128,128,128]./256;


[ftver, ftpath] = ft_version;
% in surface pial
load([ftpath  '\template\anatomy\surface_white_left.mat']);
template_lh = mesh; %clear mesh;

load([ftpath  '\template\anatomy\surface_white_right.mat']);
template_rh = mesh; %clear mesh;

load([ftpath  '\template\anatomy\surface_pial_both.mat']);
template_whole = mesh; %clear mesh;

ft_plot_mesh(template_whole, 'facealpha', 0.2);

% Whole brain
fig_savename1 = 'Fig2B_1';
fig_savename2 = 'Fig2B_2';
fig_savename3 = 'Fig2B_3';
for fig = 1
figure(fig);clf;
ft_plot_mesh(template_whole, 'facealpha', 0.2);hold on;
for s = 1 : length(elec_info)
    % label lobes
    if isempty(elec_info{s}(1).MNI152); continue; end
    cord = cat(1,elec_info{s}.MNI152);
    plot3(cord(:,1),cord(:,2),cord(:,3),'ro','markerfacecolor',cmap_subject(s,:),'markeredgecolor',cmap_subject(s,:),'markersize',3);hold on
end
hold off

if plot_save
    saveas(gcf,[figsavepath,fig_savename1,'.emf'])
    saveas(gcf,[figsavepath,fig_savename1,'.png'])
end

view([-80 12]);
if plot_save
    saveas(gcf,[figsavepath,fig_savename2,'.emf'])
    saveas(gcf,[figsavepath,fig_savename2,'.png'])
end

view([0 22]);

if plot_save
    saveas(gcf,[figsavepath,fig_savename3,'.emf'])
    saveas(gcf,[figsavepath,fig_savename3,'.png'])
end
end