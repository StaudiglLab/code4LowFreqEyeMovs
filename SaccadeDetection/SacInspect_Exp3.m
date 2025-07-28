%% Visual inspection of detected saccade (Exp 3)
% mark the false positive saccades 
% Run this after automatic saccade detection (SacDetect_Exp3)

%% load data
codepath = '';
addpath(genpath( [codepath,'\subfunctions\']))


datapath = [codepath,'\datafiles\'];
subinfopath = [codepath,'\datafiles\Subjects_Exp3\ExpInfo\'];
sacsavepath = [codepath,'\datafiles\Subjects_Exp3\saccades\'];

load([datapath,'ET_clean_saccades_Exp3.mat'])
load([datapath,'ET_Exp3.mat']);
%%
tp_epoch = [-1.3,4];
tp_inst = [-1.3,4];
tp_select = dsearchn((tp_epoch(1):1/600:tp_epoch(2))',[tp_inst(1),tp_inst(2)]');

sub = 3;

for s=sub
    disp(s)
    clearvars  data_s  data_o
    data_o = cat(3,FixPET_all(s).pvon{:});
    data_s = data_o(tp_select(1):tp_select(2),:,:);

    for itrial = 1:size(data_s,3)
        data_raw_left(s).trial{1,itrial} = data_s(:,[1:2,5],itrial)';
        data_raw_left(s).time{1,itrial} = tp_inst(1):1/600:tp_inst(2);
        
        data_raw_right(s).trial{1,itrial} = data_s(:,[3:4,6],itrial)';
        data_raw_right(s).time{1,itrial} = tp_inst(1):1/600:tp_inst(2);
         
    end
   
    
end
clearvars  data_s  data_o

% automatic mark the blink (all chunks with NaNs mroe than 50 ms)
bl_len_thr = 50.*(600/1000); % more than 50 ms of NaNs chunk
for s = sub
    for itrial = 1:length(data_raw_left(s).trial)
        Nandata_L = isnan(sum(data_raw_left(s).trial{1,itrial},1));
        Nandata_R = isnan(sum(data_raw_right(s).trial{1,itrial},1));
        
        blink_area = bwboundaries((Nandata_L+Nandata_R)==2);
        blk_count = 0;
        blink_automark{s}{itrial} = [];
        for i_blk = 1:length(blink_area)
            if (max(blink_area{i_blk}(:,2))- min(blink_area{i_blk}(:,2)))>bl_len_thr
                blk_count=blk_count+1;
                blink_automark{s}{itrial}(blk_count,:) = [min(blink_area{i_blk}(:,2)),max(blink_area{i_blk}(:,2))];
            end
        end
    end   
    
end
%%
range_dva = [-10,10];% dva range for visualization
time_range = [-1.3,4];

time_point = time_range(1):1/600:time_range(2);


%%
for s = sub
    %% load screen parameters
    load([subinfopath,sprintf('TaskInfo_Subject_%d_%s.mat',s,'Exp3')]);
    stimparami = [];
    stimparami.xlen = cfg.screen.w_size(1);
    stimparami.ylen = cfg.screen.w_size(2);
    stimparami.hPxl = cfg.screen.dim(3);
    stimparami.vPxl = cfg.screen.dim(4);
    
    DPP_x = visAngPerPixel(stimparami.xlen, cfg.screen.dist, stimparami.hPxl);
    DPP_y = visAngPerPixel(stimparami.ylen, cfg.screen.dist, stimparami.vPxl);
 
    %% plot  
    
    for i_trial = 1:length(data_raw_left(s).trial)   
        
        % plot
        % right
        tmp_r=data_raw_right(s).trial{i_trial}; 
        tmp_interpo_r = [data_ET_right(s).trial{i_trial}(1,:)./1920;data_ET_right(s).trial{i_trial}(2,:)./1080];
        
        % left
        tmp_l=data_raw_left(s).trial{i_trial}; 
        tmp_interpo_l = [data_ET_left(s).trial{i_trial}(1,:)./1920;data_ET_left(s).trial{i_trial}(2,:)./1080];
        
        % change to dva from center
        tmp_r(1,:) =  tmp_r(1,:)*stimparami.hPxl*DPP_x-cfg.screen.wc(1)*DPP_x;
        tmp_r(2,:) =  tmp_r(2,:)*stimparami.vPxl*DPP_y-cfg.screen.wc(2)*DPP_y;
        tmp_r(3,:) =  tmp_r(3,:)-( max(tmp_r(3,:))-range_dva(2)); % up to 5
        tmp_interpo_r(1,:) = tmp_interpo_r(1,:)*stimparami.hPxl*DPP_x-cfg.screen.wc(1)*DPP_x;
        tmp_interpo_r(2,:) = tmp_interpo_r(2,:)*stimparami.vPxl*DPP_y-cfg.screen.wc(2)*DPP_y;
        
        tmp_l(1,:) =  tmp_l(1,:)*stimparami.hPxl*DPP_x-cfg.screen.wc(1)*DPP_x+diff(range_dva)+2;
        tmp_l(2,:) =  tmp_l(2,:)*stimparami.vPxl*DPP_y-cfg.screen.wc(2)*DPP_y+diff(range_dva)+2;
        tmp_l(3,:) =  tmp_l(3,:)-( max(tmp_l(3,:))-range_dva(2))+diff(range_dva)+2; % 
        tmp_interpo_l(1,:) = tmp_interpo_l(1,:)*stimparami.hPxl*DPP_x-cfg.screen.wc(1)*DPP_x+diff(range_dva)+2;
        tmp_interpo_l(2,:) = tmp_interpo_l(2,:)*stimparami.vPxl*DPP_y-cfg.screen.wc(2)*DPP_y+diff(range_dva)+2;

        %
        f = figure('units','normalized','outerposition',[0.1 0.1 .9 .9]);
        
        
        plot(data_raw_right(s).time{i_trial},tmp_r,'linewidth',2);hold on
        plot(data_raw_right(s).time{i_trial},tmp_interpo_r);
        plot(data_raw_left(s).time{i_trial},tmp_l,'linewidth',2);
        plot(data_raw_left(s).time{i_trial},tmp_interpo_l);
        
        y_range = [range_dva(1),range_dva(2)+diff(range_dva)+2];
        x_range =time_range;
        ylim(y_range)
        xlim([x_range(1)-0.15,x_range(2)])
        yline(range_dva(2))
        yline(range_dva(2)+2)  
        ax1=plot(data_raw_left(s).time{i_trial},zeros(1,length(data_raw_left(s).time{i_trial}))+range_dva(2)+1,'linewidth',1,'Color','black'); % referece line for marking blink
        axL=plot(data_raw_left(s).time{i_trial},zeros(1,length(data_raw_left(s).time{i_trial}))+range_dva(2)+3,'linewidth',1,'Color','black'); % referece line for marking MS L
        axR=plot(data_raw_right(s).time{i_trial},zeros(1,length(data_raw_left(s).time{i_trial}))+range_dva(1)+2,'linewidth',1,'Color','black'); % referece line for marking MS R
        axRej = plot(x_range(1)-0.1,range_dva(2)+1,'o','Color','blue');
        axRej_L = plot(x_range(1)-0.1,range_dva(2)+3,'o','Color','blue');
        axRej_R = plot(x_range(1)-0.1,range_dva(1)+2,'o','Color','blue');
        yticks([range_dva(1),0,range_dva(2),range_dva(2)+2,diff(range_dva)+2,2*diff(range_dva)+range_dva(1)+2])
        yticklabels({int2str(range_dva(1)),'0',int2str(range_dva(2)),int2str(range_dva(1)),'0',int2str(range_dva(2))})
        
        
        title(['Subject: ' num2str(s) ' || Trial ' num2str(i_trial) ' / ' num2str(length(data_raw_left(s).trial))]);
        
        % Left saccades
        movement_L = movement_clean_left{s}{i_trial};
        if ~isempty(movement_L)
            movs2plot=(movement_L+tp_inst(1)*600)./600;          

            for m=1:size(movs2plot,1)
                xline(movs2plot(m,1),'LineStyle','--','Color',[0.9,0.9,0.9])
                v1 = [[movs2plot(m,1); movs2plot(m,2); movs2plot(m,2); movs2plot(m,1)], [y_range(1)+diff(range_dva)+2; y_range(1)+diff(range_dva)+2; y_range(2); y_range(2)]];
                f1 = [1 2 3 4];
                patch('Faces',f1,'Vertices',v1,'FaceColor','red','FaceAlpha',.1);
                
            end
        end
        % right saccades
        movement_R = movement_clean_right{s}{i_trial};
        if ~isempty(movement_R)
            movs2plot=(movement_R+tp_inst(1)*600)./600;          

            for m=1:size(movs2plot,1)
                xline(movs2plot(m,1),'LineStyle','--','Color',[0.9,0.9,0.9])
                v1 = [[movs2plot(m,1); movs2plot(m,2); movs2plot(m,2); movs2plot(m,1)], [y_range(1); y_range(1); y_range(1)+diff(range_dva); y_range(1)+diff(range_dva)]];
                f1 = [1 2 3 4];
                patch('Faces',f1,'Vertices',v1,'FaceColor','red','FaceAlpha',.1);
                
            end
        end
        
        % blink auto detect
        blk_auto = blink_automark{s}{i_trial};
        if ~isempty(blk_auto)
            movs2plot=(blk_auto+tp_inst(1)*600)./600;          
            % add ms
            for m=1:size(movs2plot,1)
                %xline(movs2plot(m,1),'LineStyle','--','Color',[0.9,0.9,0.9])
                v1 = [[movs2plot(m,1); movs2plot(m,2); movs2plot(m,2); movs2plot(m,1)], [y_range(1)+diff(range_dva); y_range(1)+diff(range_dva); y_range(1)+diff(range_dva)+2; y_range(1)+diff(range_dva)+2]];
                f1 = [1 2 3 4];
                patch('Faces',f1,'Vertices',v1,'FaceColor','blue','FaceAlpha',.1);
                
            end
        end
        
        brush on
        w = waitforbuttonpress;
        
        while w == 0
            
            w = waitforbuttonpress;
            disp('mouse')
            if w==1
                disp('Key press')
                %     close
                
            end
        end
        
        blink_line = (get(ax1, 'BrushData')); % marked on the ref line (blinks or not saccades)
        ms_line_L = (get(axL, 'BrushData'));% saccades on the L eye to discard
        ms_line_R = (get(axR, 'BrushData'));% saccades on the R eye to discard
        trl_rej = (get(axRej, 'BrushData'));% reject trial
        trl_rej_L = (get(axRej_L, 'BrushData'));% reject L eye for this trial
        trl_rej_R = (get(axRej_R, 'BrushData'));% reject R eye for this trial

        
        % use marked area to discard saccades
        ms_area_L=find(ms_line_L);
        ms_area_R=find(ms_line_R);
        mark_area_blk=find(blink_line);
        
        brush off
        close
        
        
        trl_rej_ind{i_trial} = trl_rej;
        trl_rej_Eye_L{i_trial} = trl_rej_L;
        trl_rej_Eye_R{i_trial} = trl_rej_R;
        
        aa=[];
        bb=[];
        if ~isempty(movement_L)
            
            [aa,bb]=intersect(movement_L(:,1)',ms_area_L);
            
            if ~isempty(aa)
                badmovs4biggerintervals_L{i_trial}=aa;
            end
            
        end
        movement_L(bb,:)=[];
        movement_cleanblink_L{i_trial}=movement_L;
        
        aa=[];
        bb=[];
        
        
        if ~isempty(movement_R)
            
            [aa,bb]=intersect(movement_R(:,1)',ms_area_R);
            
            if ~isempty(aa)
                badmovs4biggerintervals_R{i_trial}=aa;
            end
            
        end
        
        movement_R(bb,:)=[];
        movement_cleanblink_R{i_trial}=movement_R;
        
        % blink (combine auto and manually marked)
        blink_dat = zeros(1,length(data_raw_left(s).time{i_trial}));
        for i_blk = 1:size(blk_auto,1)
            blink_dat(1,blk_auto(i_blk,1):blk_auto(i_blk,2))=1;
        end
        blink_dat(mark_area_blk)=1;
        blink_mark{i_trial} = blink_dat;
        
        clear movement movement_clean
         
    end
    
  %% save
    save([sacsavepath,sprintf('Sub%d_ET_cleansaccades_%s.mat',s,'Exp3')],...
    'movement_cleanblink_L','movement_cleanblink_R','badmovs4biggerintervals_L','badmovs4biggerintervals_R',...
    'blink_mark','trl_rej_ind','trl_rej_Eye_L','trl_rej_Eye_R','-v7.3');
end


