function [ind,ori] = BilateralOverlap(a,b,min_sp,method)
%% Overlapping tw between a & b
% a: N1 * 2 (N1 trials, first column starting sample point, second column end sample point)
% b: N2 * 2 (N2 trials, first column starting sample point, second column end sample point)
% min_sp : minimal (bi)sample point to be considered as MS
% ori: MS coming from a (marked as -1) or b (as 1) or both (2)
%%
if nargin <2
    error('At least 2 inputs')
elseif nargin ==2
    min_sp = 7; % default  (~=12ms Engbert & Kliegl)
    method = 'overlap';
elseif nargin ==3
    method = 'overlap'; % default  (overlapping area only)
end 
ind = [];ori=[];
total_samp_a = zeros(1,max([a(:);b(:)])); % cover until the max sample point in input
total_samp_b = zeros(1,max([a(:);b(:)])); % cover until the max sample point in input
for i_trl = 1:size(a,1)
    total_samp_a(1,a(i_trl,1):a(i_trl,2))=1;
end
for i_trl = 1:size(b,1)
    total_samp_b(1,b(i_trl,1):b(i_trl,2))=1;
end
switch method
    case 'overlap'
        bi_samp_ind = bwconncomp((total_samp_a+total_samp_b)==2);
        trl_count = 0;
        for i_trl = 1:length(bi_samp_ind.PixelIdxList)
            if length(bi_samp_ind.PixelIdxList{i_trl})<min_sp
            else
                trl_count = trl_count+1;
                ind(trl_count,:) = [bi_samp_ind.PixelIdxList{i_trl}(1),bi_samp_ind.PixelIdxList{i_trl}(end)];
            end
        end
     case 'either'
        bi_samp_ind = bwconncomp((total_samp_a+total_samp_b)>0);
        trl_count = 0;
        for i_trl = 1:length(bi_samp_ind.PixelIdxList)
            if length(bi_samp_ind.PixelIdxList{i_trl})<min_sp
            else
                trl_count = trl_count+1;
                ind(trl_count,:) = [bi_samp_ind.PixelIdxList{i_trl}(1),bi_samp_ind.PixelIdxList{i_trl}(end)];
                
                if sum(total_samp_a(bi_samp_ind.PixelIdxList{i_trl}))>0 && sum(total_samp_b(bi_samp_ind.PixelIdxList{i_trl}))>0
                    ori(trl_count) = 2;
                elseif sum(total_samp_a(bi_samp_ind.PixelIdxList{i_trl}))>0 && sum(total_samp_b(bi_samp_ind.PixelIdxList{i_trl}))==0
                    ori(trl_count) = -1;
                elseif sum(total_samp_a(bi_samp_ind.PixelIdxList{i_trl}))==0 && sum(total_samp_b(bi_samp_ind.PixelIdxList{i_trl}))>0
                    ori(trl_count) = 1;
                end
                
            end
        end
          
end 
