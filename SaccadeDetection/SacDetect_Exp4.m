%% Saccade detection for Exp 4 (Step 1)
% run this first, automatic detection of saccades
codepath = '';
addpath(genpath( [codepath,'\subfunctions\']))


datapath = [codepath,'\datafiles\'];
datafigspath =  [codepath,'\data4figs\'];
figsavepath = [codepath,'\figures\'];
subinfopath = [codepath,'\datafiles\Subjects_Exp4\ExpInfo\'];


savename = [datapath,'ET_Exp4.mat'];
load(savename)

%% Mark NaNs
tp_epoch = [-1.3,4];
tp_inst = [-1.3,4];
tp_select = dsearchn((tp_epoch(1):1/1000:tp_epoch(2))',[tp_inst(1),tp_inst(2)]');

ext_tp = 100; % mark 100ms before and after nan also as nan (to take out blink residual)
for s=1:length(ET_all_eyelink)
    disp(s)
    clearvars  dat  data_s
    dat = cat(3,ET_all_eyelink(s).pvon{:});
    dat = dat([1,2,4,5],:,:);
    dat(dat==0) = nan;% the readout from ASC eye data read nan as 0

    data_s(:,1,:) = permute(dat(1,tp_select(1):tp_select(2),:),[2,1,3]);
    data_s(:,2,:) = permute(dat(2,tp_select(1):tp_select(2),:),[2,1,3]);
    data_s(:,3,:) = permute(dat(3,tp_select(1):tp_select(2),:),[2,1,3]);
    data_s(:,4,:) = permute(dat(4,tp_select(1):tp_select(2),:),[2,1,3]);
    %% a newly added section, replace here data -100 to 100 of nan also as nan
    for idat = 1:size(data_s,2)
        for itrl = 1:size(data_s,3)
            dat_nan = bwconncomp(isnan(data_s(:,idat,itrl)'));
            if dat_nan.NumObjects ~= 0
                for icl = 1:dat_nan.NumObjects
                    ind_nan = dat_nan.PixelIdxList{icl};
                    data_s(max(ind_nan(1)-ext_tp,1):min(ind_nan(end)+ext_tp,size(data_s,1)),idat,itrl) = nan;
                end
            else % do nothing
                
            end
        end
    end
    clearvars idat  itrl  dat_nan  icl   ind_nan   
    %%
    n_nans_left(s) = sum(sum(sum(isnan(data_s(:,1:2,:)),1),2) ~=0);
    n_nans_right(s) = sum(sum(sum(isnan(data_s(:,3:4,:)),1),2) ~=0);
    for itrial = 1:size(data_s,3)
        data_ET_left(s).trial{1,itrial} = data_s(:,1:2,itrial)';
        data_ET_left(s).time{1,itrial} = tp_inst(1):1/1000:tp_inst(2);
        
        data_ET_right(s).trial{1,itrial} = data_s(:,3:4,itrial)';
        data_ET_right(s).time{1,itrial} = tp_inst(1):1/1000:tp_inst(2);
        
        %% dealing with nans (interpolate by starting and end of NaNs, and register its position)
        CC_X = bwconncomp(isnan(data_s(:,1,itrial)')); CC_Y = bwconncomp(isnan(data_s(:,2,itrial)'));
        if CC_X.NumObjects ~= 0
        NaN_register_left{s}{1,itrial} = CC_X.PixelIdxList;
        else
            NaN_register_left{s}{1,itrial} = [];
        end 
        if CC_Y.NumObjects ~= 0
        NaN_register_left{s}{2,itrial} = CC_Y.PixelIdxList;
        else
            NaN_register_left{s}{2,itrial} = [];
        end
        clearvars CC_X  CC_Y
        
        CC_X = bwconncomp(isnan(data_s(:,3,itrial)')); CC_Y = bwconncomp(isnan(data_s(:,4,itrial)'));
        if CC_X.NumObjects ~= 0
        NaN_register_right{s}{1,itrial} = CC_X.PixelIdxList;
        else
            NaN_register_right{s}{1,itrial} = [];
        end 
        if CC_Y.NumObjects ~= 0
        NaN_register_right{s}{2,itrial} = CC_Y.PixelIdxList;
        else
            NaN_register_right{s}{2,itrial} = [];
        end
        clearvars CC_X  CC_Y
    end
   
    
end

%% Interpolate the NaNs
dat_length = length(tp_epoch(1):1/1000:tp_epoch(2));
for s = 1:length(ET_all_eyelink)
    disp(s);
    for itrial = 1:length(data_ET_right(s).trial)
        %% left X
        if ~isempty(NaN_register_left{s}{1,itrial})
            for i_cl = 1:length(NaN_register_left{s}{1,itrial}) % number of cluster of NaN
                ind_cl = NaN_register_left{s}{1,itrial}{i_cl};
                if length(ind_cl)==dat_length
                    data_ET_left(s).trial{itrial}(1,:) =0; % total trial nans
                else
                    if ind_cl(1) == 1 % starting point is first timepoint
                        data_ET_left(s).trial{itrial}(1,ind_cl(1):ind_cl(end)) = data_ET_left(s).trial{itrial}(1,ind_cl(end)+1);
                    elseif ind_cl(end) == dat_length % end point is the last
                        data_ET_left(s).trial{itrial}(1,ind_cl(1):ind_cl(end)) = data_ET_left(s).trial{itrial}(1,ind_cl(1)-1);
                    else
                        data_ET_left(s).trial{itrial}(1,ind_cl(1):ind_cl(end)) = linspace(data_ET_left(s).trial{itrial}(1,ind_cl(1)-1),data_ET_left(s).trial{itrial}(1,ind_cl(end)+1),length(ind_cl));
                    end
                    clearvars ind_cl
                end
            end    
        end
        
        %% left Y
        if ~isempty(NaN_register_left{s}{2,itrial})
            for i_cl = 1:length(NaN_register_left{s}{2,itrial}) % number of cluster of NaN
                ind_cl = NaN_register_left{s}{2,itrial}{i_cl};
                if length(ind_cl)==dat_length
                    data_ET_left(s).trial{itrial}(2,:) =0; % total trial nans
                else
                    if ind_cl(1) == 1 % starting point is first timepoint
                        data_ET_left(s).trial{itrial}(2,ind_cl(1):ind_cl(end)) = data_ET_left(s).trial{itrial}(2,ind_cl(end)+1);
                    elseif ind_cl(end) == dat_length % end point is the last
                        data_ET_left(s).trial{itrial}(2,ind_cl(1):ind_cl(end)) = data_ET_left(s).trial{itrial}(2,ind_cl(1)-1);
                    else
                        data_ET_left(s).trial{itrial}(2,ind_cl(1):ind_cl(end)) = linspace(data_ET_left(s).trial{itrial}(2,ind_cl(1)-1),data_ET_left(s).trial{itrial}(2,ind_cl(end)+1),length(ind_cl));
                    end
                    clearvars ind_cl
                end
            end    
        end
        
       %% right X
        if ~isempty(NaN_register_right{s}{1,itrial})
            for i_cl = 1:length(NaN_register_right{s}{1,itrial}) % number of cluster of NaN
                ind_cl = NaN_register_right{s}{1,itrial}{i_cl};
                if length(ind_cl)==dat_length
                    data_ET_right(s).trial{itrial}(1,:) =0; % total trial nans
                else
                    if ind_cl(1) == 1 % starting point is first timepoint
                        data_ET_right(s).trial{itrial}(1,ind_cl(1):ind_cl(end)) = data_ET_right(s).trial{itrial}(1,ind_cl(end)+1);
                    elseif ind_cl(end) == dat_length % end point is the last
                        data_ET_right(s).trial{itrial}(1,ind_cl(1):ind_cl(end)) = data_ET_right(s).trial{itrial}(1,ind_cl(1)-1);
                    else
                        data_ET_right(s).trial{itrial}(1,ind_cl(1):ind_cl(end)) = linspace(data_ET_right(s).trial{itrial}(1,ind_cl(1)-1),data_ET_right(s).trial{itrial}(1,ind_cl(end)+1),length(ind_cl));
                    end
                    clearvars ind_cl
                end
            end    
        end
        
        %% right Y
        if ~isempty(NaN_register_right{s}{2,itrial})
            for i_cl = 1:length(NaN_register_right{s}{2,itrial}) % number of cluster of NaN
                ind_cl = NaN_register_right{s}{2,itrial}{i_cl};
                if length(ind_cl)==dat_length
                    data_ET_right(s).trial{itrial}(2,:) =0; % total trial nans
                else
                    if ind_cl(1) == 1 % starting point is first timepoint
                        data_ET_right(s).trial{itrial}(2,ind_cl(1):ind_cl(end)) = data_ET_right(s).trial{itrial}(2,ind_cl(end)+1);
                    elseif ind_cl(end) == dat_length % end point is the last
                        data_ET_right(s).trial{itrial}(2,ind_cl(1):ind_cl(end)) = data_ET_right(s).trial{itrial}(2,ind_cl(1)-1);
                    else
                        data_ET_right(s).trial{itrial}(2,ind_cl(1):ind_cl(end)) = linspace(data_ET_right(s).trial{itrial}(2,ind_cl(1)-1),data_ET_right(s).trial{itrial}(2,ind_cl(end)+1),length(ind_cl));
                    end
                    clearvars ind_cl
                end
            end    
        end
        
    end
end

%% detect saccades
cfg                 = [];
cfg.method =  'velocity2D';
cfg.velocity2D.mindur = 12; %minimum microsaccade durantion in samples (=12ms Engbert & Kliegl)
cfg.velocity2D.velthres =4;%threshold for velocity outlier detection 
ft_warning off
for s = 1:length(ET_all_eyelink)
    disp(s)
    for itrial = 1:length(data_ET_right(s).trial)
        data_left.trial = data_ET_left(s).trial{itrial};
        data_left.time = data_ET_left(s).time{1,itrial};
        data_left.fsample = 1000;
        data_left.label = {'EYE_HORIZONTAL', 'EYE_VERTICAL'};
        
        [~, movement_left{s}{itrial}] = ft_detect_movement(cfg, data_left);
        
        data_right.trial = data_ET_right(s).trial{itrial};
        data_right.time = data_ET_right(s).time{1,itrial};
        data_right.fsample = 1000;
        data_right.label = {'EYE_HORIZONTAL', 'EYE_VERTICAL'};
        
        data_right.label = {'EYE_HORIZONTAL', 'EYE_VERTICAL'};
        [~,  movement_right{s}{itrial}] = ft_detect_movement(cfg, data_right);
        
        clearvars data_left data_right
    end
end
ft_warning on

%% Kick out the (Micro) saccade that coincide with NaNs and Outside of the pic area

n_spthr = 10; % Max number of sample points with NaNs that allow to be interpolated 
% add constraint, no more than 10 NaNs (10 ms) conseq points is allowed
for s = 1:length(ET_all_eyelink)
    for itrial = 1:size(NaN_register_left{s},2)
        if ~isempty(NaN_register_left{s}{1,itrial}) 
            NaN_register_left{s}{1,itrial}(find(cellfun(@length,NaN_register_left{s}{1,itrial})<=n_spthr))=[];
        end
        if ~isempty(NaN_register_left{s}{2,itrial})
            NaN_register_left{s}{2,itrial}(find(cellfun(@length,NaN_register_left{s}{2,itrial})<=n_spthr))=[];
        end
    end
end 
for s = 1:length(ET_all_eyelink)
    for itrial = 1:size(NaN_register_right{s},2)
        if ~isempty(NaN_register_right{s}{1,itrial}) 
            NaN_register_right{s}{1,itrial}(find(cellfun(@length,NaN_register_right{s}{1,itrial})<=n_spthr))=[];
        end
        if ~isempty(NaN_register_right{s}{2,itrial})
            NaN_register_right{s}{2,itrial}(find(cellfun(@length,NaN_register_right{s}{2,itrial})<=n_spthr))=[];
        end
    end
end 

for s = 1:length(ET_all_eyelink)
    disp(s)
    %% load screen parameter
    load([subinfopath,sprintf('TaskInfo_Subject_%d_%s.mat',s,'Exp4')]);
    
    stimparami = [];
    stimparami.xlen = cfg.screen.w_size(1);
    stimparami.ylen = cfg.screen.w_size(2);
    stimparami.hPxl = cfg.screen.dim(3);
    stimparami.vPxl = cfg.screen.dim(4);
    
    scrcenter = cfg.screen.wc;
    [pxl_x,pxl_y]           = ang2px(cfg.screen.va,cfg.screen.dist,stimparami); % fix size    
    
    viewdva_x =cfg.screen.va_img(1); 
    viewdva_y =cfg.screen.va_img(2); 
    %% define region of interest
    x_lim = [scrcenter(1)-pxl_x*viewdva_x/2,scrcenter(1)+pxl_x*viewdva_x/2];
    y_lim = [scrcenter(2)-pxl_y*viewdva_y/2,scrcenter(2)+pxl_y*viewdva_y/2];
    
    %% left
    MS_counter_left{s} = zeros(1,length(data_ET_left(s).trial));
    for itrial = 1:length(data_ET_left(s).trial)
        Eyedata=data_ET_left(s).trial{itrial};
        movement_clean_left{s}{itrial} = [];
        if ~isempty(movement_left{s}{itrial})
            i=0;%clean trial counter
          for i_ms = 1:size(movement_left{s}{itrial},1)
              ms_ind = movement_left{s}{itrial}(i_ms,:);
              if isempty(NaN_register_left{s}{1,itrial}) && isempty(NaN_register_left{s}{2,itrial})
                  nan_ind = [];
              elseif ~isempty(NaN_register_left{s}{1,itrial}) && isempty(NaN_register_left{s}{2,itrial})
                  nan_ind = unique(cat(1,NaN_register_left{s}{1,itrial}{:}));
              elseif isempty(NaN_register_left{s}{1,itrial}) && ~isempty(NaN_register_left{s}{2,itrial})
                  nan_ind = unique(cat(1,NaN_register_left{s}{2,itrial}{:}));
              else
                  nan_ind = unique([cat(1,NaN_register_left{s}{1,itrial}{:});cat(1,NaN_register_left{s}{2,itrial}{:})]);
              end
              if isempty(intersect(ms_ind(1):ms_ind(2),nan_ind))
                  if ~(min(Eyedata(1,ms_ind(1):ms_ind(2)))<x_lim(1)||max(Eyedata(1,ms_ind(1):ms_ind(2)))>x_lim(2)||...
                          min(Eyedata(2,ms_ind(1):ms_ind(2)))<y_lim(1)||max(Eyedata(2,ms_ind(1):ms_ind(2)))>y_lim(2))% if inside the range
                      MS_counter_left{s}(itrial) = MS_counter_left{s}(itrial)+1;
                      i=i+1;
                      movement_clean_left{s}{itrial}(i,:) = ms_ind;
                  end
              end 
              clearvars ms_ind nan_ind
          end 
        end
    end
    
    %% right
    MS_counter_right{s} = zeros(1,length(data_ET_right(s).trial));
    for itrial = 1:length(data_ET_right(s).trial)
        Eyedata=data_ET_right(s).trial{itrial};
        movement_clean_right{s}{itrial} = [];
        if ~isempty(movement_right{s}{itrial})
            i=0;%clean trial counter
          for i_ms = 1:size(movement_right{s}{itrial},1)
              ms_ind = movement_right{s}{itrial}(i_ms,:);
              if isempty(NaN_register_right{s}{1,itrial}) && isempty(NaN_register_right{s}{2,itrial})
                  nan_ind = [];
              elseif ~isempty(NaN_register_right{s}{1,itrial}) && isempty(NaN_register_right{s}{2,itrial})
                  nan_ind = unique(cat(1,NaN_register_right{s}{1,itrial}{:}));
              elseif isempty(NaN_register_right{s}{1,itrial}) && ~isempty(NaN_register_right{s}{2,itrial})
                  nan_ind = unique(cat(1,NaN_register_right{s}{2,itrial}{:}));
              else
                  nan_ind = unique([cat(1,NaN_register_right{s}{1,itrial}{:});cat(1,NaN_register_right{s}{2,itrial}{:})]);
              end
              if isempty(intersect(ms_ind(1):ms_ind(2),nan_ind))
                  if ~(min(Eyedata(1,ms_ind(1):ms_ind(2)))<x_lim(1)||max(Eyedata(1,ms_ind(1):ms_ind(2)))>x_lim(2)||...
                          min(Eyedata(2,ms_ind(1):ms_ind(2)))<y_lim(1)||max(Eyedata(2,ms_ind(1):ms_ind(2)))>y_lim(2))% if inside the range
                      MS_counter_right{s}(itrial) = MS_counter_right{s}(itrial)+1;
                      i=i+1;
                      movement_clean_right{s}{itrial}(i,:) = ms_ind;
                  end
              end 
              clearvars ms_ind nan_ind
          end 
        end
    end
    
end


%% save output
save([datapath,'ET_clean_saccades_Exp4.mat'],'data_ET_right','data_ET_left','MS_counter_right','MS_counter_left','movement_clean_right','movement_clean_left','-v7.3');
