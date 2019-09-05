function [cutT_opt_all,minp_all,p_LR_ALL,rankf]=findoptimal(Features,survival,least_num)
% description
% input: Features: a m*n matrix with m examples and n features.
%        survival: m*2 matrix; the first column with the survival 
%        time and second column with survival status ('0': event, '1': censored). 
%        least_num: set the least number of each group
%
% output: cutT_opt_all: the optimal cut-off
%         minp_all: the logrank p value with optimal cut-off
%         p_LR_ALL: the logrank p value with different cut-off (a vector)
%         rankf: features rank based on the log-rank p
%
%ljlubme@gmail.com
%Southern Medical University

sizeF = size(Features);
F_length = size(Features,1);
F_num = size(Features,2);
for i=1:F_num
    feature = Features(:,i);
    rankf = sort(feature,'ascend');
    ns=1;
    for j=least_num:F_length-least_num 
        cutT = rankf(j);
        id_l = find(feature>cutT);
        id_s = find(feature<=cutT);
        Time_Outcome_L = survival(id_l,:);
        Time_Outcome_S = survival(id_s,:);
        [p, t1, T1, t2,T2]= logrank(Time_Outcome_S,Time_Outcome_L);
        p_LR(ns) = p;
        ns=ns+1;
    end
        [minp,idj] = min(p_LR);
        id_opt = least_num+idj-1;
        cutT_opt = rankf(id_opt);
        cutT_opt_all(i)=cutT_opt;
        minp_all(i)=minp;
        p_LR_ALL{i}=p_LR;
end
end
        