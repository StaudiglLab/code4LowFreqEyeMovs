%% Gaze density calculation for Exp3
codepath = '';
addpath(genpath( [codepath,'\subfunctions\']))


datapath = [codepath,'\datafiles\'];
datafigspath =  [codepath,'\data4figs\'];
figsavepath = [codepath,'\figures\'];
subinfopath = [codepath,'\datafiles\Subjects_Exp3\ExpInfo\'];
sacsavepath = [codepath,'\datafiles\Subjects_Exp3\saccades\'];

% load behavioral data, epoch eye data (raw & preprocessed)
load( [datapath,'Beh_Exp3.mat'])
load([datapath,'ET_Exp3.mat'])
load([datapath,'ET_clean_saccades_Exp3.mat']);
load([datapath,'eyedom_Exp3.mat']); % subject dominant eye (1 = L; 2 = R)
%%
% Mark the eye that has better quality data
%
for s = 1:length(FixPET_all)
    load([sacsavepath,sprintf('Sub%d_ET_cleansaccades_%s.mat',s,'Exp3')])
    
    
    trl_ind_L = find(cellfun(@mean, trl_rej_Eye_L)==1);
    trl_ind_R = find(cellfun(@mean, trl_rej_Eye_R)==1);
    

    if length(trl_ind_R)>length(trl_ind_L) % Left eye better trials
            Eyequal(s) = 1;
    elseif length(trl_ind_R)<length(trl_ind_L) % Right eye better trials
            Eyequal(s) = 2;
    elseif length(trl_ind_R)==length(trl_ind_L) % if same then choose dom eye
            Eyequal(s) = Sub_DE(s);
    end 
        
end 


tp_epoch = [-1.3,4];
tp_inst = [-1.3,4];
tp_select = dsearchn((tp_epoch(1):1/600:tp_epoch(2))',[tp_inst(1),tp_inst(2)]');

tp_final = [0,4];
tp_slct_final = dsearchn((tp_epoch(1):1/600:tp_epoch(2))',[tp_final(1),tp_final(2)]');
for s=1:length(FixPET_all)
    disp(s)
    clearvars  data_s  data_o
    
    data_o = cat(3,FixPET_all(s).pvon{:});
    data_s = data_o(tp_select(1):tp_select(2),:,:);

  

    for itrial = 1:size(data_s,3)      
        data_L = data_s(:,1:2,itrial)';
        data_R = data_s(:,3:4,itrial)';       
        
        PV_left(s).trial{1,itrial} = data_L(:,tp_slct_final(1):tp_slct_final(2));
        PV_left(s).time{1,itrial} = tp_final(1):1/600:tp_final(2);
        
        PV_right(s).trial{1,itrial} = data_R(:,tp_slct_final(1):tp_slct_final(2));
        PV_right(s).time{1,itrial} = tp_final(1):1/600:tp_final(2);
        
        clearvars data_L data_R
    end
end
clearvars data_o data_s FixPET_all

%% Filter data with criteria
% 1: mark NaN pre & post blink
t_range = 100; % time beofre and after blink in ms
tp = t_range/1000*600;
for s = 1:length(PV_right)
    disp(s)
    load([sacsavepath,sprintf('Sub%d_ET_cleansaccades_%s.mat',s,'Exp3')])
    
    for itrial = 1:length(PV_right(s).trial)
        blink_idx = blink_mark{itrial}(tp_slct_final(1):tp_slct_final(2));
        blink_area = bwboundaries(blink_idx==1);
        if ~isempty(blink_area)
            for i_bl = 1:size(blink_area,1)
                lim_bl = [min(blink_area{i_bl}(:,2)),max(blink_area{i_bl}(:,2))];
                lim_bl1 = max(lim_bl(1)-tp,1);
                lim_bl2 = min(lim_bl(2)+tp,length(PV_left(s).trial{1,itrial}));
                PV_left(s).trial{1,itrial}(:,lim_bl1:lim_bl2) = NaN;
                
                PV_right(s).trial{1,itrial}(:,lim_bl1:lim_bl2) = NaN;
                
                clearvars lim_bl lim_bl1 lim_bl2
            end
        end
        
    end
    
    
end

% 2: mark rej trials (it's redundant as same info can be found in TrialEyeInfo.mat)
viewdva_x =28;
viewdva_y =21;
for s = 1:length(PV_right)
    disp(s)

    % exclude trials with more than 50% nans in both eye
    L = cellfun(@(x) isnan(x), PV_left(s).trial, 'UniformOutput', false);
    L = cellfun(@(x) sum(x(:)), L, 'UniformOutput', false);
    trl_ind_L = find([L{:}]>(numel(PV_left(s).trial{1})/2)); clearvars L 
    
    R = cellfun(@(x) isnan(x), PV_right(s).trial, 'UniformOutput', false);
    R = cellfun(@(x) sum(x(:)), R, 'UniformOutput', false);
    trl_ind_R = find([R{:}]>(numel(PV_right(s).trial{1})/2)); clearvars R 
    
    trl_rej = intersect(trl_ind_L,trl_ind_R);
    
    % note the excluded trial
    trl_idx{s} = ones(1,length(PV_left(s).trial));
    trl_idx{s}(trl_rej) = 0;
  
    % exclude good quality eye with more than 50% nans
    switch Eyequal(s)
        case 1 % left
             trl_idx{s}(trl_ind_L) = 0;
        case 2 % right
             trl_idx{s}(trl_ind_R) = 0;
    end
    
    clearvars trl_ind trl_rej
    
    
    %% fine tune the search zone 
    % give the flexibility to filter out the trials with gaze postion exceeding different
    % restriction (right now we allow 28 x 21 same as other experiments)
        
    load([subinfopath,sprintf('TaskInfo_Subject_%d_%s.mat',s,'Exp3')]);
    stimparami = [];
    stimparami.xlen = cfg.screen.w_size(1);
    stimparami.ylen = cfg.screen.w_size(2);
    stimparami.hPxl = cfg.screen.dim(3);
    stimparami.vPxl = cfg.screen.dim(4);
    
    scrcenter = cfg.screen.wc;
    [pxl_x,pxl_y]           = ang2px(cfg.screen.va,cfg.screen.dist,stimparami); % fix size    
    
    
    % define region of interest
    x_lim = [scrcenter(1)-pxl_x*viewdva_x/2,scrcenter(1)+pxl_x*viewdva_x/2]/(scrcenter(1)*2);
    y_lim = [scrcenter(2)-pxl_y*viewdva_y/2,scrcenter(2)+pxl_y*viewdva_y/2]/(scrcenter(2)*2);
    
    
    trlL_max = cellfun(@(x) max(x,[],2), PV_left(s).trial, 'UniformOutput', false);
    trlL_min = cellfun(@(x) min(x,[],2), PV_left(s).trial, 'UniformOutput', false);
    
    trlL_rej_x = cellfun(@(x) sum(x>[x_lim(2),y_lim(2)]'), trlL_max, 'UniformOutput', false);
    trlL_rej_y = cellfun(@(x) sum(x<[x_lim(1),y_lim(1)]'), trlL_min, 'UniformOutput', false);
    
    trlR_max = cellfun(@(x) max(x,[],2), PV_right(s).trial, 'UniformOutput', false);
    trlR_min = cellfun(@(x) min(x,[],2), PV_right(s).trial, 'UniformOutput', false);
    
    trlR_rej_x = cellfun(@(x) sum(x>[x_lim(2),y_lim(2)]'), trlR_max, 'UniformOutput', false);
    trlR_rej_y = cellfun(@(x) sum(x<[x_lim(1),y_lim(1)]'), trlR_min, 'UniformOutput', false);
    
    %trl_idx{s}(([trlL_rej_x{:}]+[trlL_rej_y{:}]+[trlR_rej_x{:}]+[trlR_rej_y{:}])~=0)=0;
    switch Eyequal(s)
        case 1 % left
            trl_idx{s}(([trlL_rej_x{:}]+[trlL_rej_y{:}])~=0)=0;
        case 2 % right
            trl_idx{s}(([trlR_rej_x{:}]+[trlR_rej_y{:}])~=0)=0;
    end
    

end

%% do gaze
% focus on better eye
viewdva_x =28;
viewdva_y =21;
n_bin_x = 80;n_bin_y = 60;
for s = 1:length(PV_right)
    disp(s)
    switch Eyequal(s)
        case 1 % left
            tmpet = PV_left(s);
        case 2 % right
            tmpet = PV_right(s);
    end

    load([subinfopath,sprintf('TaskInfo_Subject_%d_%s.mat',s,'Exp3')]);
    
    stimparami = [];
    stimparami.xlen = cfg.screen.w_size(1);
    stimparami.ylen = cfg.screen.w_size(2);
    stimparami.hPxl = cfg.screen.dim(3);
    stimparami.vPxl = cfg.screen.dim(4);
    
    scrcenter = cfg.screen.wc;
    [pxl_x,pxl_y]           = ang2px(cfg.screen.va,cfg.screen.dist,stimparami); % fix size    
    
    
    %% define region of interest
    x_lim = [scrcenter(1)-pxl_x*viewdva_x/2,scrcenter(1)+pxl_x*viewdva_x/2];
    y_lim = [scrcenter(2)-pxl_y*viewdva_y/2,scrcenter(2)+pxl_y*viewdva_y/2];
    %%
    for trl=1:length(tmpet.trial)

        tmp=horzcat(tmpet.trial{trl});
        xpospre=tmp(1,:);
        ypospre=tmp(2,:);
        xpre = xpospre;
        ypre=ypospre;
        
        %% heat map
        % Bin the data:
        pts_x = linspace(x_lim(1)/stimparami.hPxl, x_lim(2)/stimparami.hPxl, n_bin_x+1);
        pts_y = linspace(y_lim(1)/stimparami.vPxl, y_lim(2)/stimparami.vPxl, n_bin_y+1);
        
        N = histcounts2(ypre(:), xpre(:), pts_y, pts_x)./sum(~isnan(xpre));% total valid sample point
        
        %  Create Gaussian filter matrix:
        [xG, yG] = meshgrid(-5:5);
        sigma = 2.5;
        g = exp(-xG.^2./(2.*sigma.^2)-yG.^2./(2.*sigma.^2));
        g = g./sum(g(:));
        %% convert to FT
        tmp=conv2(N, g, 'same');
        freq.freq= linspace(0, 1, n_bin_y);
        freq.powspctrm=tmp;
        pow=zeros(1,numel(freq.powspctrm(:,1)),numel(freq.powspctrm(1,:)));
        pow(1,:,:)=freq.powspctrm;
        freq.powspctrm=pow;
        freq.time= linspace(0, 1, n_bin_x);
        freq.label={'gaze'};
        freq.dimord= 'chan_freq_time';
        gaze{s}{trl}=freq;
        gaze2D{s}{trl} = squeeze(pow);
    end
end
%%
savename = [datapath,'PVstim_Gaze_Exp3.mat'];
save(savename,'gaze','gaze2D','trl_idx','-v7.3')
