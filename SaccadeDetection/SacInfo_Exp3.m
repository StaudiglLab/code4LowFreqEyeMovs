%% Saccade info (Step 3)
% Run this after Automatic detection of saccades (SacDetect_Exp3) and
% visual inspection (SacInspect_Exp3)
codepath = '';
addpath(genpath( [codepath,'\subfunctions\']))


datapath = [codepath,'\datafiles\'];
datafigspath =  [codepath,'\data4figs\'];
figsavepath = [codepath,'\figures\'];
subinfopath = [codepath,'\datafiles\Subjects_Exp3\ExpInfo\'];
sacsavepath = [codepath,'\datafiles\Subjects_Exp3\saccades\'];
% behavioral results
savename = [datapath,'Beh_Exp3.mat'];
load(savename)
%%
t_min = 50; % in ms free of the next saccades, to not count same saccade multiple times
for sac=1
load([datapath,'ET_Exp3.mat'])
load([datapath,'ET_clean_saccades_Exp3.mat']); clearvars movement_clean_left movement_clean_right

% get pupil
tp_epoch = [-1.3,4];
tp_inst = [-1.3,4];
tp_select = dsearchn((tp_epoch(1):1/600:tp_epoch(2))',[tp_inst(1),tp_inst(2)]');

for s=1:length(FixPET_all)
    disp(s)
    
    clearvars  data_s  data_o
    data_o = cat(3,FixPET_all(s).pvon{:});
    data_s = data_o(tp_select(1):tp_select(2),:,:);
    
    load([sacsavepath,sprintf('Sub%d_ET_cleansaccades_%s.mat',s,'Exp3')])
   

    for itrial = 1:size(data_s,3)
        data_raw_left(s).trial{1,itrial} = data_s(:,[1:2,5],itrial)';
        data_raw_left(s).time{1,itrial} = tp_inst(1):1/600:tp_inst(2);
        
        data_raw_right(s).trial{1,itrial} = data_s(:,[3:4,6],itrial)';
        data_raw_right(s).time{1,itrial} = tp_inst(1):1/600:tp_inst(2);
         
        data_raw_left(s).trial{1,itrial}(:,blink_mark{1,itrial}==1) = NaN;
        data_raw_right(s).trial{1,itrial}(:,blink_mark{1,itrial}==1) = NaN;
    end
   

end
clearvars data_o data_s FixPET_all

%% Filter Sac with criteria
% 1: take out the saccades in marked rejection trials 
for s = 1:length(Mem)
    disp(s)
    load([sacsavepath,sprintf('Sub%d_ET_cleansaccades_%s.mat',s,'Exp3')])
    
    %% remove saccade in rejected (L/R) eye
   trl_single = find(cellfun(@mean, trl_rej_Eye_R)==1);
    
    for i_trl = 1:length(trl_single)
        movement_cleanblink_R{trl_single(i_trl)} = [];
    end
   
    trl_single = find(cellfun(@mean, trl_rej_Eye_L)==1);
    
    for i_trl = 1:length(trl_single)
        movement_cleanblink_L{trl_single(i_trl)} = [];
    end
    clearvars trl_single
    %%
    MS_L{s} = movement_cleanblink_L;
    MS_R{s} = movement_cleanblink_R;
    trl_ind = find(cellfun(@mean, trl_rej_ind)==1);
    
    for i_trl = 1:length(trl_ind)
        MS_L{s}{trl_ind(i_trl)} = [];
        MS_R{s}{trl_ind(i_trl)} = [];
    end
    
    % exclude trials with more than 50% nans in both eye
    L = cellfun(@(x) isnan(x), data_raw_left(s).trial, 'UniformOutput', false);
    L = cellfun(@(x) sum(x(:)), L, 'UniformOutput', false);
    trl_ind_L = find([L{:}]>(numel(data_raw_left(s).trial{1})/2)); clearvars L 
    
    R = cellfun(@(x) isnan(x), data_raw_right(s).trial, 'UniformOutput', false);
    R = cellfun(@(x) sum(x(:)), R, 'UniformOutput', false);
    trl_ind_R = find([R{:}]>(numel(data_raw_right(s).trial{1})/2)); clearvars R 
    
    trl_rej = intersect(trl_ind_L,trl_ind_R);
    
    for i_trl = 1:length(trl_rej)
        MS_L{s}{trl_rej(i_trl)} = [];
        MS_R{s}{trl_rej(i_trl)} = [];
    end
    
    % exclude trials with un marked eye having more than 50% of nans
    trl_single = find(cellfun(@mean, trl_rej_Eye_R)==1);
    trl_rej_L = intersect(trl_ind_L,trl_single);
    
    for i_trl = 1:length(trl_rej_L)
        MS_L{s}{trl_rej_L(i_trl)} = [];
        MS_R{s}{trl_rej_L(i_trl)} = [];
    end
    
    trl_single = find(cellfun(@mean, trl_rej_Eye_L)==1);
    trl_rej_R = intersect(trl_ind_R,trl_single);
    
    for i_trl = 1:length(trl_rej_R)
        MS_L{s}{trl_rej_R(i_trl)} = [];
        MS_R{s}{trl_rej_R(i_trl)} = [];
    end
    clearvars trl_single
    
    
    % note the excluded trial
    trl_cleaneye_idx{s} = ones(1,length(movement_cleanblink_R));
    trl_cleaneye_idx{s}(trl_ind) = 0;
    trl_cleaneye_idx{s}(trl_rej) = 0;
    trl_cleaneye_idx{s}(trl_rej_L) = 0;
    trl_cleaneye_idx{s}(trl_rej_R) = 0;
    clearvars trl_ind trl_rej
    
end

% 2: take out saccades 100 ms before and after blink
t_range = 100; % time beofre and after blink in ms
tp = t_range/1000*600;
for s = 1:length(Mem)
    disp(s)
    load([sacsavepath,sprintf('Sub%d_ET_cleansaccades_%s.mat',s,'Exp3')])   
    
    
    trl_ind = find(cellfun(@sum, blink_mark));
    
    for i_trl = 1:length(trl_ind)
        blink_area = bwboundaries(blink_mark{trl_ind(i_trl)});
        MS_area_L = MS_L{s}{trl_ind(i_trl)};
        MS_area_R = MS_R{s}{trl_ind(i_trl)};
        for i_blk = 1:length(blink_area)
            blk_ind = [min(blink_area{i_blk}(:,2))-tp,max(blink_area{i_blk}(:,2))+tp];
            if ~isempty(MS_area_L)
                trl_rm = find(MS_area_L(:,1)<=blk_ind(2) & MS_area_L(:,2)>=blk_ind(1));
                MS_L{s}{trl_ind(i_trl)}(trl_rm,:)=[];
                MS_area_L(trl_rm,:)=[];
                clearvars trl_rm
            end
            
            if ~isempty(MS_area_R)
                trl_rm = find(MS_area_R(:,1)<=blk_ind(2) & MS_area_R(:,2)>=blk_ind(1));
                MS_R{s}{trl_ind(i_trl)}(trl_rm,:)=[];
                MS_area_R(trl_rm,:)=[];
                clearvars trl_rm
            end

        end
    
    end
end

% 3: saccade in either eye
for s =  1:length(Mem)
    disp(s)
    for i_trial = 1:length(MS_R{s})
        MS_bi{s}{i_trial} = [];
        if ~isempty(MS_L{s}{i_trial}) && ~isempty(MS_R{s}{i_trial})
        left_ind = MS_L{s}{i_trial}(:,[1,2]);
        right_ind = MS_R{s}{i_trial}(:,[1,2]);
        MS_bi{s}{i_trial} = BilateralOverlap(left_ind,right_ind,7,'either');
        elseif ~isempty(MS_L{s}{i_trial}) && isempty(MS_R{s}{i_trial}) % left eye only
            only_ind = MS_L{s}{i_trial}(:,[1,2]);
            trl_count = 0;
            for i_ms = 1:size(only_ind,1)
                    trl_count = trl_count +1;
                    MS_bi{s}{i_trial}(trl_count,:) = only_ind(i_ms,:);
            end
        elseif isempty(MS_L{s}{i_trial}) && ~isempty(MS_R{s}{i_trial}) % right eye only
            only_ind = MS_R{s}{i_trial}(:,[1,2]);
            trl_count = 0;
            for i_ms = 1:size(only_ind,1)
                    trl_count = trl_count +1;
                    MS_bi{s}{i_trial}(trl_count,:) = only_ind(i_ms,:);
            end
        end
    end
end

% 4: take out MS that was counted multiple time
tp_min = t_min/1000*600;
for s =  1:length(Mem)
    disp(s)
    for i_trl = 1:length(MS_bi{s})
        if ~isempty(MS_bi{s}{i_trl})
            MS_trl = MS_bi{s}{i_trl};
            while 1
                if sum((MS_trl(2:end,1)-MS_trl(1:end-1,2))<=tp_min)==0
                    break
                end
                MS_trl(find((MS_trl(2:end,1)-MS_trl(1:end-1,2))<=tp_min,1,'first')+1,:) = [];    
            end
            MS_bi{s}{i_trl} = MS_trl; clearvars MS_trl
        end 
    end
end

base_l = 0;%1.3;% counting from 1.3 on or include baseline
base_end = 1.3+4;% t_sep before end of image

for s=1:length(Mem)
    for i_trl = 1:length(MS_bi{s})
        
        MS_L_pvon{s}{i_trl} = [];
        if ~isempty(MS_L{s}{i_trl})
            ms = MS_L{s}{i_trl};
            ms = ms(ms(:,1)>base_l*600,:);
            ms = ms(ms(:,1)<base_end*600,:);
            MS_L_pvon{s}{i_trl} =  ms;
        end
        
        MS_R_pvon{s}{i_trl} = [];
        if ~isempty(MS_R{s}{i_trl})
            ms = MS_R{s}{i_trl};
            ms = ms(ms(:,1)>base_l*600,:);
            ms = ms(ms(:,1)<base_end*600,:);
            MS_R_pvon{s}{i_trl} =  ms;
        end
        
        MS_bi_pvon{s}{i_trl} = [];
        if ~isempty(MS_bi{s}{i_trl})
            ms = MS_bi{s}{i_trl};
            ms = ms(ms(:,1)>base_l*600,:);
            ms = ms(ms(:,1)<base_end*600,:);
            MS_bi_pvon{s}{i_trl} =  ms;
            
        end
        
    end
    
end 

% Interpolated nans for amplitude and direction calculation 
load([datapath,'ET_clean_saccades_Exp3.mat']);

% get MS direction
for AGL = 1
    for s=1:length(Mem)
        disp(s)
        %% load screen parameter
        load([subinfopath,sprintf('TaskInfo_Subject_%d_%s.mat',s,'Exp3')]);
        
        dist = cfg.screen.dist;
        xlen = cfg.screen.w_size(1);
        ylen = cfg.screen.w_size(2);
        hPxl = cfg.screen.dim(3);
        vPxl = cfg.screen.dim(4);
        
        load([sacsavepath,sprintf('Sub%d_ET_cleansaccades_%s.mat',s,'Exp3')])
        
        %% binocular
        trl_ind_bi = find(cellfun(@isempty,MS_bi_pvon{s})==0);
        for i_trl = 1:length(trl_ind_bi)
            trl_idx = trl_ind_bi(i_trl);
            if mean(trl_rej_Eye_R{trl_idx}) ~=1 && mean(trl_rej_Eye_L{trl_idx})~=1 
                x_pos = nanmean([data_ET_left(s).trial{trl_idx}(1,:);data_ET_right(s).trial{trl_idx}(1,:)],1);
                y_pos = nanmean([data_ET_left(s).trial{trl_idx}(2,:);data_ET_right(s).trial{trl_idx}(2,:)],1);
            elseif mean(trl_rej_Eye_R{trl_idx}) ==1 && mean(trl_rej_Eye_L{trl_idx})~=1 % only left eye valid
                x_pos = data_ET_left(s).trial{trl_idx}(1,:);
                y_pos = data_ET_left(s).trial{trl_idx}(2,:);
            elseif mean(trl_rej_Eye_R{trl_idx}) ~=1 && mean(trl_rej_Eye_L{trl_idx})==1 % only right eye valid
                x_pos = data_ET_right(s).trial{trl_idx}(1,:);
                y_pos = data_ET_right(s).trial{trl_idx}(2,:);
            else
                error('Not a valid trial')
            end
            
            for i_ms = 1:size(MS_bi_pvon{s}{trl_idx},1) % all MS in a trial
                % calculate sac direction
                t_on = MS_bi_pvon{s}{trl_idx}(i_ms,1);
                t_off = MS_bi_pvon{s}{trl_idx}(i_ms,2);
                x1 = x_pos(t_on);x2 = x_pos(t_off);
                y1 = y_pos(t_on);y2 = y_pos(t_off);
                ms_angl_bi{s}{i_trl}(i_ms) =atan2d(y2-y1,x2-x1);
                
                % calculate amplitude in dvg
                dvg_x = abs(x2-x1)*visAngPerPixel(xlen, dist, hPxl);
                dvg_y = abs(y2-y1)*visAngPerPixel(ylen, dist, vPxl);
                ms_amp_bi{s}{i_trl}(i_ms) =sqrt(dvg_x^2 + dvg_y^2);
                
            end
            clearvars trl_idx x_pos y_pos
        end
    end
    
end


    
end
%% get sac locked info
t_ET  = -1.3:1/600:4;
% add constraint about timewindow for where saccace onsets are included
lim_ET1 = find((t_ET*1000)>-1.3*1000,1,'first');% the first datapoint in ET trial that can be counted
lim_ET2 = find((t_ET*1000)<4*1000,1,'last');% the latest datapoint in ET trial that can be counted

for s = 1:length(Mem)
    ms_bi_trlind{s} = find(~cellfun(@isempty,MS_bi_pvon{s}));
end


for s = 1:length(Mem)
    disp(s)
    
    %% Extract info time-locked to sac onset
    i_trl = 0;MS_tp = [];beh_acc=[];MS_dir = []; MS_amp = [];confi_rt = []; trl_number = [];%MS_dur= [];
    for itrial = 1:length(ms_bi_trlind{s})

            i_trl = i_trl+1;
            MS_tp{i_trl} = MS_bi_pvon{s}{ms_bi_trlind{s}(itrial)};
            beh_acc(i_trl) = Mem{s}(ms_bi_trlind{s}(itrial));
            MS_dir{i_trl} = ms_angl_bi{s}{itrial};
            MS_amp{i_trl} = ms_amp_bi{s}{itrial};
            %MS_dur{i_trl} = diff(MS_bi_pvon{s}{ms_bi_trlind{s}(itrial)},[],2)'*1/600*1000;
            confi_rt(i_trl) = Confi{s}(ms_bi_trlind{s}(itrial));
            trl_number(i_trl) = ms_bi_trlind{s}(itrial); % current trial

    end
    %     end
    
    % saccade onset
    trl_count = 0; beh_ms=[];dir_ms=[];amp_ms=[];dur_ms=[]; onset_ms = []; offset_ms = [];confi_ms = [];trlnum_ms = [];
    for i_trl = 1:length(MS_tp)
        for i_ms = 1:size(MS_tp{i_trl},1) % more than 1 MS
            ind_MS_on = MS_tp{i_trl}(i_ms,1); % ind of MS onset
            if ind_MS_on<=lim_ET2 && ind_MS_on>=lim_ET1
                trl_count = trl_count+1;
    
                beh_ms(trl_count) = beh_acc(i_trl);
                dir_ms(trl_count) = MS_dir{i_trl}(i_ms);
                amp_ms(trl_count) = MS_amp{i_trl}(i_ms);
                %dur_ms(trl_count) = MS_dur{i_trl}(i_ms);
                onset_ms(trl_count) = MS_tp{i_trl}(i_ms,1);
                offset_ms(trl_count) = MS_tp{i_trl}(i_ms,2);
                confi_ms(trl_count) = confi_rt(i_trl);
                trlnum_ms(trl_count) = trl_number(i_trl);
            end
        end
        
    end
    
    %% create structure
    Sac_trials = [];
    Sac_trials.trialinfo = trlnum_ms'; % correspoding trl index
    Sac_trials.trialinfo(:,2) = beh_ms; % add subseq beh
    Sac_trials.trialinfo(:,3) = confi_ms; % add subseq confi rating
    Sac_trials.trialinfo(:,4) = dir_ms; % add dir in angle
    Sac_trials.trialinfo(:,5) = amp_ms; % add amp
    Sac_trials.trialinfo(:,6) = onset_ms; % add onset
    Sac_trials.trialinfo(:,7) = offset_ms; % add offset (both correspond to ET data index [-1.3,4])
    
    %% add trials with zero saccade detected (saccade related info will be filled as NaNs)
    
    % trials with no saccades but clean
    enc_trl = ones(216,1);
    
    cl_trl = intersect(find(trl_cleaneye_idx{s}'),find(enc_trl));
    Sac_trials.trialinfo = Sac_trials.trialinfo(ismember(Sac_trials.trialinfo(:,1),cl_trl),:);
    
    
    enc_trl(unique(Sac_trials.trialinfo(:,1)))=0; % 1 stands for no saccade detected
    %
    trl_withoutsac = intersect(find(trl_cleaneye_idx{s}'),find(enc_trl));
    
    no_sac_trialinfo      = trl_withoutsac;
    no_sac_trialinfo(:,2) = Mem{s}(trl_withoutsac); % add subseq beh
    no_sac_trialinfo(:,3) = Confi{s}(trl_withoutsac); % add subseq confi rating
    no_sac_trialinfo(:,4:7) =  NaN; % saccade related info as NaN
    %% concate with trials with detected saccade
    Sac_trials.trialinfo = cat(1,Sac_trials.trialinfo,no_sac_trialinfo);
    %% reorder
    [~,ro_idx] = sort(Sac_trials.trialinfo(:,1));
    Sac_trials.trialinfo = Sac_trials.trialinfo(ro_idx,:);
    %% save
    savename = [sacsavepath,sprintf('Sac_PV_Info_Sub%d_%s.mat',s,'Exp3')];
    save(savename,'Sac_trials','-v7.3')
end

