%% Gaze density calculation for Exp2
codepath = '';
addpath(genpath( [codepath,'\subfunctions\']))


datapath = [codepath,'\datafiles\'];
datafigspath =  [codepath,'\data4figs\'];
figsavepath = [codepath,'\figures\'];
subinfopath = [codepath,'\datafiles\Subjects_Exp2\ExpInfo\'];
sacsavepath = [codepath,'\datafiles\Subjects_Exp2\saccades\'];

% load behavioral data, epoch eye data (raw & preprocessed)
load( [datapath,'Beh_Exp2.mat'])
load([datapath,'ET_Exp2.mat'])

%% Tobii eye tracker (P01 & P02)
for tb = 1
load([datapath,'ET_clean_saccades_Exp2_tb.mat']); 
for s = 1:2
    load([sacsavepath,sprintf('Sub%d_ET_cleansaccades_%s.mat',s,'Exp2')])
    
    
    trl_ind_L = find(cellfun(@mean, trl_rej_Eye_L)==1);
    trl_ind_R = find(cellfun(@mean, trl_rej_Eye_R)==1);
    

    if length(trl_ind_R)>length(trl_ind_L) % Left eye better trials
            Eyequal(s) = 1;
    elseif length(trl_ind_R)<length(trl_ind_L) % Right eye better trials
            Eyequal(s) = 2;
    elseif length(trl_ind_R)==length(trl_ind_L) % if same then choose right eye
            Eyequal(s) = 2;
    end 
        
end 


%% 
tp_epoch = [-1.3,4];
tp_inst = [-1.3,4];
tp_select = dsearchn((tp_epoch(1):1/600:tp_epoch(2))',[tp_inst(1),tp_inst(2)]');

tp_final = [0,4];
tp_slct_final = dsearchn((tp_epoch(1):1/600:tp_epoch(2))',[tp_final(1),tp_final(2)]');
for s=1:2
    disp(s)
    
    clearvars  data_s  data_o
    data_o = cat(3,ET_all_tobii(s).pvon{:});
    data_s = data_o(tp_select(1):tp_select(2),:,:);
    
    load([sacsavepath,sprintf('Sub%d_ET_cleansaccades_%s.mat',s,'Exp2')])

    for itrial = 1:size(data_s,3)
        
        data_L = data_s(:,1:2,itrial)';
        data_R = data_s(:,3:4,itrial)';
          
        data_L(:,blink_mark{1,itrial}==1) = NaN;
        data_R(:,blink_mark{1,itrial}==1) = NaN;    
        
        
        PV_left(s).trial{1,itrial} = data_L(:,tp_slct_final(1):tp_slct_final(2));
        PV_left(s).time{1,itrial} = tp_final(1):1/600:tp_final(2);
        
        PV_right(s).trial{1,itrial} = data_R(:,tp_slct_final(1):tp_slct_final(2));
        PV_right(s).time{1,itrial} = tp_final(1):1/600:tp_final(2);
        
        clearvars data_L data_R
    end
end
clearvars data_o data_s ET_all_tobii

%% do gaze
n_bin_x = 80;n_bin_y = 60;
for s = 1:2
    disp(s)
    switch Eyequal(s)
        case 1 % left
            tmpet = PV_left(s);
        case 2 % right
            tmpet = PV_right(s);
    end
    
    load([subinfopath,sprintf('TaskInfo_Subject_%d_%s_Session1.mat',s,'Exp2')]);

    stimparami = [];
    stimparami.xlen = cfg.screen.w_size(1);
    stimparami.ylen = cfg.screen.w_size(2);
    stimparami.hPxl = cfg.screen.dim(3);
    stimparami.vPxl = cfg.screen.dim(4);
    scrcenter = cfg.screen.wc;

    viewdva_x =cfg.screen.va_img(1); 
    viewdva_y =cfg.screen.va_img(2); 

    [pxl_x1,pxl_y1]           = ang2px(cfg.screen.va,cfg.screen.dist,stimparami); % fix size    

    % session 2
    load([subinfopath,sprintf('TaskInfo_Subject_%d_%s_Session2.mat',s,'Exp2')]);

    [pxl_x2,pxl_y2]           = ang2px(cfg.screen.va,cfg.screen.dist,stimparami); % fix size    

    %% define region of interest
    x_lim1 = [scrcenter(1)-pxl_x1*viewdva_x/2,scrcenter(1)+pxl_x1*viewdva_x/2];
    y_lim1 = [scrcenter(2)-pxl_y1*viewdva_y/2,scrcenter(2)+pxl_y1*viewdva_y/2];

    x_lim2 = [scrcenter(1)-pxl_x2*viewdva_x/2,scrcenter(1)+pxl_x2*viewdva_x/2];
    y_lim2 = [scrcenter(2)-pxl_y2*viewdva_y/2,scrcenter(2)+pxl_y2*viewdva_y/2];
    %%
    for trl=1:length(tmpet.trial)

        if trl<=72 % first session
            x_lim = x_lim1;y_lim = y_lim1;
        else
            x_lim = x_lim2;y_lim = y_lim2;
        end

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

savename = [datapath,'PVstim_Gaze_Exp2_tb.mat'];
save(savename,'gaze','gaze2D','-v7.3')

end

%% Eyelink eye tracker (P03 & P04)
for el = 1
    load([datapath,'ET_clean_saccades_Exp2_el.mat']);
    for s = 1:2
        load([sacsavepath,sprintf('Sub%d_ET_cleansaccades_%s.mat',s+2,'Exp2')])
        
        
        trl_ind_L = find(cellfun(@mean, trl_rej_Eye_L)==1);
        trl_ind_R = find(cellfun(@mean, trl_rej_Eye_R)==1);
        
        
        if length(trl_ind_R)>length(trl_ind_L) % Left eye better trials
            Eyequal(s) = 1;
        elseif length(trl_ind_R)<length(trl_ind_L) % Right eye better trials
            Eyequal(s) = 2;
        elseif length(trl_ind_R)==length(trl_ind_L) % if same then choose right eye
            Eyequal(s) = 2;%  dominant eye
        end
        
    end
    
    
    tp_epoch = [-1.3,4];
tp_inst = [-1.3,4];
tp_select = dsearchn((tp_epoch(1):1/1000:tp_epoch(2))',[tp_inst(1),tp_inst(2)]');

tp_final = [0,4];
tp_slct_final = dsearchn((tp_epoch(1):1/1000:tp_epoch(2))',[tp_final(1),tp_final(2)]');
for s=1:2
    disp(s)
    clearvars  data_s  data_o
    
    dat = cat(3,ET_all_eyelink(s).pvon{:});
    dat = dat([1,2,4,5],:,:);
    dat(dat==0) = nan;

    data_s(:,1,:) = permute(dat(1,tp_select(1):tp_select(2),:),[2,1,3]);
    data_s(:,2,:) = permute(dat(2,tp_select(1):tp_select(2),:),[2,1,3]);
    data_s(:,3,:) = permute(dat(3,tp_select(1):tp_select(2),:),[2,1,3]);
    data_s(:,4,:) = permute(dat(4,tp_select(1):tp_select(2),:),[2,1,3]);
    
    for itrial = 1:size(data_s,3)      
        data_L = data_s(:,1:2,itrial)';
        data_R = data_s(:,3:4,itrial)';       
        
        PV_left(s).trial{1,itrial} = data_L(:,tp_slct_final(1):tp_slct_final(2));
        PV_left(s).time{1,itrial} = tp_final(1):1/1000:tp_final(2);
        
        PV_right(s).trial{1,itrial} = data_R(:,tp_slct_final(1):tp_slct_final(2));
        PV_right(s).time{1,itrial} = tp_final(1):1/1000:tp_final(2);
        
        clearvars data_L data_R
    end
end
clearvars data_o data_s ET_all_eyelink

%%  mark NaNs
t_range = 100; % time beofre and after blink in ms (blink residual for EyeLink)
tp = t_range/1000*1000;
for s = 1:2
    load([sacsavepath,sprintf('Sub%d_ET_cleansaccades_%s.mat',s+2,'Exp2')])
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

%% do gaze
viewdva_x =28;
viewdva_y =21;
n_bin_x = 80;n_bin_y = 60;
for s = 1:2
    disp(s)
    switch Eyequal(s)
        case 1 % left
            tmpet = PV_left(s);
        case 2 % right
            tmpet = PV_right(s);
    end

    
    load([subinfopath,sprintf('TaskInfo_Subject_%d_%s_Session1.mat',s+2,'Exp2')]);
    stimparami = [];
    stimparami.xlen = cfg.screen.w_size(1);
    stimparami.ylen = cfg.screen.w_size(2);
    stimparami.hPxl = cfg.screen.dim(3);
    stimparami.vPxl = cfg.screen.dim(4);
    scrcenter = cfg.screen.wc;

    viewdva_x =cfg.screen.va_img(1); 
    viewdva_y =cfg.screen.va_img(2); 

    [pxl_x1,pxl_y1]           = ang2px(cfg.screen.va,cfg.screen.dist,stimparami); % fix size    

    % session 2
    load([subinfopath,sprintf('TaskInfo_Subject_%d_%s_Session1.mat',s+2,'Exp2')]);

    [pxl_x2,pxl_y2]           = ang2px(cfg.screen.va,cfg.screen.dist,stimparami); % fix size    
    
    %% define region of interest
    x_lim1 = [scrcenter(1)-pxl_x1*viewdva_x/2,scrcenter(1)+pxl_x1*viewdva_x/2];
    y_lim1 = [scrcenter(2)-pxl_y1*viewdva_y/2,scrcenter(2)+pxl_y1*viewdva_y/2];

    x_lim2 = [scrcenter(1)-pxl_x2*viewdva_x/2,scrcenter(1)+pxl_x2*viewdva_x/2];
    y_lim2 = [scrcenter(2)-pxl_y2*viewdva_y/2,scrcenter(2)+pxl_y2*viewdva_y/2];
    %%
    for trl=1:length(tmpet.trial)

        if trl<=72 % first session
            x_lim = x_lim1;y_lim = y_lim1;
        else
            x_lim = x_lim2;y_lim = y_lim2;
        end

        tmp=horzcat(tmpet.trial{trl});
        xpospre=tmp(1,:);
        ypospre=tmp(2,:);
        xpre = xpospre;
        ypre=ypospre;
        
        %% heat map
        % Bin the data:
        pts_x = linspace(x_lim(1), x_lim(2), n_bin_x+1);
        pts_y = linspace(y_lim(1), y_lim(2), n_bin_y+1);
        
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
   
savename = [datapath,'PVstim_Gaze_Exp2_el.mat'];
save(savename,'gaze','gaze2D','-v7.3')
end