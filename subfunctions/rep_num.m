function seq = rep_num(seq2rep,rep_within,rep_between)
% repeat each element is a sequence rep_within times, then rep_between times
% example:  seq2rep = [1, 2, 3]
%           rep_within = 2
%           rep_between = 3
% output:   seq = [1,1,2,2,3,3,1,1,2,2,3,3,1,1,2,2,3,3]
%
%

seq=[];
for i_itr=1:rep_between
    for i_ele = 1:length(seq2rep)
        a = repmat(seq2rep(i_ele),1,rep_within);
        seq = [seq,a];
    end
end