%% Behavioral analysis for PETER Exp2

codepath = '';

datapath = [codepath,'\datafiles\'];
datafigspath =  [codepath,'\data4figs\'];
figsavepath = [codepath,'\figures\'];

savename = [datapath,'Beh_Exp2.mat'];
load(savename)

%%

for s = 1:length(Mem)
    %% HR hit rate  & FA False alarm
    hit_pr(s,1) = mean(pr_perf{s}(pr_cond{s}==1|pr_cond{s}==2)).*100;
    falram_pr(s,1) = mean(pr_perf{s}(pr_cond{s}==3|pr_cond{s}==4)==0).*100;
end 

[dpri,~] = dprime(hit_pr./100,falram_pr./100,144);

data(1,:) = [mean(cellfun(@mean,pv_perf)'*100),mean(hit_pr),mean(dpri)];
data(2,:) = [std(cellfun(@mean,pv_perf)'*100),std(hit_pr),std(dpri)];

T = array2table(data, ...
    'VariableNames', {'InOutAcc', 'Hit','dprime'}, ...
    'RowNames', {'Mean', 'SD'});

% Display the table
disp(T);

