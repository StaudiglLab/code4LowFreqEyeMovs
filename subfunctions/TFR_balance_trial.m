function [A_balance,B_balance,n_trl] = TFR_balance_trial(A,B,n_it)

if nargin<2
    error('Two conditions needed')
elseif nargin == 2
    n_it = 100;
end 

A_balance =[];
B_balance =[];

copyfield = {'label','freq','time'};
for i = 1:length(copyfield)
    A_balance.(copyfield{i}) = A.(copyfield{i});
    B_balance.(copyfield{i}) = B.(copyfield{i});
end

A_balance.dimord = 'chan_freq_time';
B_balance.dimord = 'chan_freq_time';

n_A = size(A.powspctrm,1);
n_B = size(B.powspctrm,1);

n_trl =  min(n_A,n_B);

rand_A = nan([n_it,size(A.powspctrm,[2,3,4])]);
rand_B = nan([n_it,size(B.powspctrm,[2,3,4])]);
for i=1:n_it
    if rem(i,n_it/10)==0
       fprintf('%d %% Complete \n',i/n_it*100);
    end
    
    rand_A(i,:,:,:) = nanmean(A.powspctrm(randsample(n_A,n_trl),:,:,:));
    rand_B(i,:,:,:) = nanmean(B.powspctrm(randsample(n_B,n_trl),:,:,:));
    
end

A_balance.powspctrm = squeeze(mean(rand_A));
B_balance.powspctrm = squeeze(mean(rand_B));
