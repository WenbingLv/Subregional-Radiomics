function [results]=stepwise_selection(Features,survival)
% description
% input: Features: a m*n matrix with m examples and n features.
%        survival: m*2 matrix; the first column with the survival 
%        time and second column with survival status ('0': event, '1': censored). 
%
% output: results 
%
%ljlubme@gmail.com
%Southern Medical University
%
time = survival(:,1);
status = survival(:,2);
%% log-rank test
%% multi-stepwise-forward variable selection to construct optimal model;
Features = zscore(Features);
for i=1:size(Features,2) 
    [b0,logL0,H0,stats0] = coxphfit(Features(:,i),time,'censoring',status,'baseline',0);
    logLall(i)=logL0;
end
 id_first = find(logLall == max(logLall));
if length(id_first)>1
 id_first = id_first(1);
end
 id_1 = id_first;
for i=1:size(Features,2)
    [b1,logL1,H1,stats1] = coxphfit(Features(:,id_1),time,'censoring',status,'baseline',0);
    id_2 = [id_1, i];
    [b2,logL2,H2,stats2] = coxphfit(Features(:,id_2),time,'censoring',status,'baseline',0);
    diflog = -2*(logL1-logL2);
    p12 = 1-cdf('chi2',diflog,1);
    if p12<0.05
        id_1=id_2;  
    else
 end    
end
selected_id = id_1;
% the selelcted features were input to multi-variate Cox model
[b,logL,H,stats] = coxphfit(Features(:,selected_id),time,'censoring',status,'baseline',0);
results.p = stats.p; % significance of features
results.beta = stats.beta;
results.HR = exp(stats.beta); % coefficent of features
results.features = Features(:,selected_id);
%status  0:event, 1:censored
results.Cindex = cindex(results.beta,results.features,time,status);
end