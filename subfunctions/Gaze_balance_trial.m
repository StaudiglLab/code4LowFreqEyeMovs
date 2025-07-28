function [A_balance,B_balance,n_trl] = Gaze_balance_trial(A,B,n_it)

if nargin<2
    error('Two conditions needed')
elseif nargin == 2
    n_it = 100;
end 

n_A = size(A,3);
n_B = size(B,3);

n_trl =  min(n_A,n_B);

rand_A = nan([n_it,size(A,[1,2])]);
rand_B = nan([n_it,size(B,[1,2])]);
for i=1:n_it
    if rem(i,n_it/10)==0
       fprintf('%d %% Complete \n',i/n_it*100);
    end
    
    rand_A(i,:,:) = nanmean(A(:,:,randsample(n_A,n_trl)),3);
    rand_B(i,:,:) = nanmean(B(:,:,randsample(n_B,n_trl)),3);
    
end

A_balance = squeeze(nanmean(rand_A));
B_balance = squeeze(nanmean(rand_B));