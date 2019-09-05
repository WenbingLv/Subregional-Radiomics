function [p,t1, T1,t2,T2] = logrank(varargin)
% LOGRANK Comparing survival curves of two groups using the log rank test
% Comparison of two survival curves can be done using a statistical
% hypothesis test called the log rank test. It is used to test the null
% hypothesis that there is no difference between the population survival
% curves (i.e. the probability of an event occurring at any time point is
% the same for each population). This function use the Kaplan-Meier
% procedure to estimate the survival function, so it is mandatory to download
% KMPLOT (http://www.mathworks.com/matlabcentral/fileexchange/22293).
%
% Syntax: 	logrank(x1,x2,alpha,censflag)
%      
%     Inputs:
%           X1 and X2 (mandatory)- Nx2 data matrix:
%                     (X:,1) = survival time of the i-th subject
%                     (X:,2) = censored flag 
%                             (0 if not censored; 1 if censored)
%           note that if X is a vector, all the flags of the second column
%           will be set to 0 (all data are not censored).%当X是向量时，默认没有标记
%           ALPHA (optional) - significance level (default 0.05) 
%           CENSFLAG (optional) - Censored Plot flag (default 0). If 0
%           censored data will be plotted spreaded on the horizontal
%           segment; if 1 they will be plotted at the given time of
%           censoring.
%     Outputs:
%           Kaplan-Meier plot
%           Log-rank statistics
%
%      Example: 
%           load logrankdata x1 x2
%           logrank(x1,x2)
%
%LOG-RANK TEST FOR KAPLAN-MEIER SURVIVAL FUNCTIONS
%
%--------------------------------------------------------------------------------
%UL				S.E.			z				p-value (2-tailed test) 	alpha
%--------------------------------------------------------------------------------
%6.57226		2.80788			2.16258			0.03057                     0.050
%--------------------------------------------------------------------------------
%		The survival functions are statistically different
%
%           Created by Giuseppe Cardillo
%           giuseppe.cardillo-edta@poste.it
%
% To cite this file, this would be an appropriate format:
% Cardillo G. (2008). LogRank: Comparing survival curves of two groups
% using the log rank test
% http://www.mathworks.com/matlabcentral/fileexchange/22317

%Input Error handling
args=cell(varargin);
nu=numel(args);
if isempty(nu) || nu<2
    error('Warning: Data vectors are required')
elseif nu>4
    error('Warning: Max four input data are required')
end
default.values = {[],[],0.05,0};
default.values(1:nu) = args;
[x1 x2 alpha,cflag] = deal(default.values{:});
if ~all(isfinite(x1(:))) || ~all(isnumeric(x1(:))) ...
        || ~all(isfinite(x2(:))) || ~all(isnumeric(x2(:)))
    error('Warning: all X1 and X2 values must be numeric and finite')
end
if isvector(x1) 
    x1(:,2)=0;
else
    if ~isequal(size(x1,2),2)
        error('LOGRANK requires Nx2 matrix data.');
    end
    if ~all(x1(:,2)==0 | x1(:,2)==1)
        error('Warning: all X1(:,2) values must be 0 or 1')
    end
end
if isvector(x2) 
    x2(:,2)=0;
else
    if ~isequal(size(x2,2),2)
        error('LOGRANK requires Nx2 matrix data.');
    end
    if ~all(x2(:,2)==0 | x2(:,2)==1)
        error('Warning: all X2(:,2) values must be 0 or 1')
    end
end
if nu>=3
    if isempty(alpha)
        alpha=0.05;
    else
        if ~isscalar(alpha) || ~isnumeric(alpha) || ~isfinite(alpha)
            error('Warning: it is required a numeric, finite and scalar ALPHA value.');
        end
        if alpha <= 0 || alpha >= 1 %check if alpha is between 0 and 1
            error('Warning: ALPHA must be comprised between 0 and 1.')
        end
    end
end
if nu==4
    if isempty(cflag)
        cflag=0;
    else
        if ~isscalar(cflag) || ~isnumeric(cflag) || ~isfinite(cflag)
            error('Warning: it is required a numeric, finite and scalar CENSFLAG value.');
        end
        if cflag~=0 && cflag~=1
            error('Warning: CENSFLAG value must be 0 or 1')
        end
    end
end
clear args default nu

%recall KMPLOT function to construct tables of data (table1 and table2),
%tables of censored data (table12 and table 22), Kaplan-Meier variables
%(t1, t2, T1 and T2) and Kaplan-Meier graphical data for censored data 
%(xcg and ycg).
try
    [table1 table12 t1 T1 xcg1 ycg1 lambda1]=kmplot(x1,0.05,cflag,0);
    [table2 table22 t2 T2 xcg2 ycg2 lambda2]=kmplot(x2,0.05,cflag,0);
catch ME
   disp('Download KMPLOT: http://www.mathworks.com/matlabcentral/fileexchange/22293')
   rethrow(ME);
end
T1 = T1*100;
T2 = T2*100;
ycg1 = ycg1*100;
ycg2 = ycg2*100;


%% plot both Kaplan-Meier curves
clf
hold on
S1=stairs(t1,T1,'b','Linewidth',2); %Kaplan-Meier curve for treatment 1
if ~isempty(table12)
    S3=plot(xcg1,ycg1,'b+'); %Censored data for treatment 1 (if there are)
else
    S3=[];
end
set(gca,'FontSize',12);
S2=stairs(t2(1:end-1),T2(1:end-1),'r','Linewidth',2); %Kaplan-Meier curve for treatment 2
if ~isempty(table22)
    S3=plot(xcg2(1:end-1),ycg2(1:end-1),'r+'); %Censored data for treatment 2 (if there are)
end
S2=stairs(t2,T2,'r','Linewidth',2); %Kaplan-Meier curve for treatment 2
if ~isempty(table22)
    S3=plot(xcg2,ycg2,'r+'); %Censored data for treatment 2 (if there are)
end
hold off
%set the axis properly
xmax=max([t1;t2])+1;
axis([0 xmax 0 105]);
axis square
set(gca,'YTick',0:20:105,'FontSize',11,'FontWeight','bold','Linewidth',1.3);
%add labels and legend
%title('Kaplan-Meier estimate of survival functions')
% ylabel('Progression-free survival (%)','FontSize',12,'FontWeight','bold')
ylabel('Overall survival (%)','FontSize',12,'FontWeight','bold')
% ylabel('Recurrence-free survival (%)','FontSize',12,'FontWeight','bold')
% ylabel('Metastasis-free survival (%)','FontSize',12,'FontWeight','bold')
xlabel('Time (months)','FontSize',12,'FontWeight','bold')
if isempty(S3)
    legend([S1 S2],'Treatment 1','Treatment 2','FontWeight','bold')
else
    %legend([S1 S2 S3],'Treatment 1','Treatment 2','Censored')
   %legend([S1 S2],'N0-N1','N2-N3','FontWeight','bold')
%   legend([S1 S2],'T1-T2','T3-T4','FontWeight','bold')
  %  legend([S1 S2],'M0','M1','FontWeight','bold')
    % legend([S1 S2],'Female','Male','FontWeight','bold')
  legend([S1 S2],'Low risk','High risk','FontWeight','bold')
%       legend([S1 S2],'PLT<=213e09g/L','PLT>213e09g/L','FontWeight','bold')
      % legend([S1 S2],'Compactness1<=1.46e-04','Compactness1>1.46e-04','FontWeight','bold')
          % legend([S1 S2],'Compactness2<=7.59e-06','Compactness2>7.59e-06','FontWeight','bold')
         %  legend([S1 S2],'Irregularity<=52.43','Irregularity>52.43','FontWeight','bold')
          % legend([S1 S2],'Skewness\_hist\_LLL<=0.68','Skewness\_hist\_LLL>0.68','FontWeight','bold')
          % legend([S1 S2],'Kurtosis\_hist\_LLH<=3.22','Kurtosis\_hist\_LLH>3.22','FontWeight','bold')
          % legend([S1 S2],'Entropy\_hist\_LLH<=0.01','Entropy\_hist\_LLH>0.01','FontWeight','bold')
         %  legend([S1 S2],'AUC\_CSH\_LHH<=0.45','AUC\_CSH\_LHH>0.45','FontWeight','bold')
          % legend([S1 S2],'Kurtosis\_hist\_LHH<=6.78','Kurtosis\_hist\_LHH>6.78','FontWeight','bold')
      %    legend([S1 S2],'LGLGE\_GLGLM\_Image\_32<=0.34','LGLGE\_GLGLM\_Image\_32>0.34','FontWeight','bold')
      %     legend([S1 S2],'SumEntropy\_GLCM\_LLL\_32<=3.72','SumEntropy\_GLCM\_LLL\_32>3.72','FontWeight','bold')
        %   legend([S1 S2],'GLV\_GLRLM\_LHL\_32<=0.03','GLV\_GLRLM\_LHL\_32>0.03','FontWeight','bold')
      %     legend([S1 S2],'SAM\_TFCCM\_LHL\_32<=0.01','SAM\_TFCCM\_LHL\_32>0.01','FontWeight','bold')
        %   legend([S1 S2],'RLV\_GLRLM\_LHH\_32<=3.31e-04','RLV\_GLRLM\_LHH\_32>3.31e-04','FontWeight','bold')
          % legend([S1 S2],'HGZE\_GLSZM\_HLL\_32<=288.25','HGZE\_GLSZM\_HLL\_32>288.25','FontWeight','bold')
         %  legend([S1 S2],'MaxSpe\_TS\_HLL\_32<=0.05','MaxSpe\_TS\_HLL\_32>0.05','FontWeight','bold')
          % legend([S1 S2],'LGLGE\_GLGLM\_Image\_64<=0.1','LGLGE\_GLGLM\_Image\_64>0.1','FontWeight','bold')
          % legend([S1 S2],'SumEntropy\_GLCM\_LLL\_64<=4.42','SumEntropy\_GLCM\_LLL\_64>4.42','FontWeight','bold')
         %  legend([S1 S2],'IMC2\_GLCM\_LLL\_64<=0.77','IMC2\_GLCM\_LLL\_64>0.77','FontWeight','bold')
          % legend([S1 S2],'GLV\_GLSZM\_HLL\_64<=2e03','GLV\_GLSZM\_HLL\_64>2e03','FontWeight','bold')
         %  legend([S1 S2],'LGLGE\_GLGLM\_Image\_128<=0.02','LGLGE\_GLGLM\_Image\_128>0.02','FontWeight','bold')
          % legend([S1 S2],'SumEntropy\_GLCM\_LLL\_128<=5.12','SumEntropy\_GLCM\_LLL\_128>5.12','FontWeight','bold')
          % legend([S1 S2],'Entropy\_GLCM\_HLL\_128<=8.13','Entropy\_GLCM\_HLL\_128>8.13','FontWeight','bold')
          % legend([S1 S2],'SumEntropy\_GLCM\_LLL\_256<=5.81','SumEntropy\_GLCM\_LLL\_256>5.81','FontWeight','bold')
          % legend([S1 S2],'GLV\_GLRLM\_LLL\_256<=0.05','GLV\_GLRLM\_LLL\_256>0.05','FontWeight','bold')

           % legend([S1 S2],'HGB<=133g/L','HGB>133g/L','FontWeight','bold')
 %  legend([S1 S2],'LYM<=1.52e09g/L','LYM>1.52e09g/L','FontWeight','bold')
%  legend([S1 S2],'EBV DNA<=7220 copies/mL','EBV DNA>7220 copies/mL','FontWeight','bold')
    % legend([S1 S2 S3],'SumEntropy<=0.72','SumEntropy>0.72','Censored')
    % legend([S1 S2],'SAM-TFCCM<=0.01','SAM-TFCCM>0.01','FontWeight','bold')
   %  legend([S1 S2],'Age<=53','Age>53','FontWeight','bold')
    % legend('boxoff')
end
%  title('CT spectral cluster model');
%text(5,10,'p=3.91e-05','FontSize',12,'FontWeight','bold')%C14Tr
% text(5,10,'p=0.0059','FontSize',12,'FontWeight','bold')%C14Te
% text(5,10,'p=1.0089e-05','FontSize',12,'FontWeight','bold')%C
% text(5,10,'p=1.0089e-05','FontSize',12,'FontWeight','bold')%C
%text(5,10,'p=1.6832e-08','FontSize',12,'FontWeight','bold')%C14
%title('Clinical+PET model');
% text(5,10,'p=7.1850e-07','FontSize',12,'FontWeight','bold')%C_PET
% text(5,10,'p=1.9661e-08','FontSize',12,'FontWeight','bold')%PET14
% title('Clinical+CT model');
%  text(5,10,'p=9.64e-10','FontSize',12,'FontWeight','bold')%C_CT_14_Tr
%   text(5,10,'p=3.17e-05','FontSize',12,'FontWeight','bold')%C_CT_14_Te
% text(5,10,'p=9.5535e-06','FontSize',12,'FontWeight','bold')%C_CT
% text(5,10,'p=2.2915e-13','FontSize',12,'FontWeight','bold')%C_CT14
%  title('PET model');
%  text(5,10,'p=0.0469','FontSize',12,'FontWeight','bold')%PET14Tr
%     text(5,10,'p=0.3975','FontSize',12,'FontWeight','bold')%PET14Te
% text(5,10,'p=8.4751e-04','FontSize',12,'FontWeight','bold')%PET
 %text(5,10,'p=6.0632e-04','FontSize',12,'FontWeight','bold')%PET14
%  title('CT model');
%   text(5,10,'p=0.0102','FontSize',12,'FontWeight','bold')%CT_14_Te
 %text(5,10,'p=3.07e-07','FontSize',12,'FontWeight','bold')%CT_14_Tr
% text(5,10,'p=8.9886e-07','FontSize',12,'FontWeight','bold')%CT
 %text(5,10,'p=1.6569e-05','FontSize',12,'FontWeight','bold')%CT14
% title('PET+CT model');
% text(5,10,'p=2.73e-07','FontSize',12,'FontWeight','bold')%CT_PET14_Tr
% text(5,10,'p=0.0102','FontSize',12,'FontWeight','bold')%CT_PET14_Te
% text(5,10,'p=7.64e-04','FontSize',12,'FontWeight','bold')%CT_PET
 %text(5,10,'p=1.7332e-06','FontSize',12,'FontWeight','bold')%CT_PET14
% title('Clinical+PET+CT model');
%  text(5,10,'p=8.08e-06','FontSize',12,'FontWeight','bold')%C_CT_PET14_Te
%  text(5,10,'p=3.16e-12','FontSize',12,'FontWeight','bold')%C_CT_PET14_Tr
% text(5,10,'p=7.5465e-06','FontSize',12,'FontWeight','bold')%C_CT_PET
 %text(5,10,'p=3.9644e-10','FontSize',12,'FontWeight','bold')%C_CT_PET14
%clear S1 S2 S3 xmax xcg1 ycg1 xcg2 ycg2 t1 t2 T1 T2

%% Full-blown LOGRANK procedure
%Merge the first columns of Table1 and Table2 (time intervals)
%and pick-up unique values
A=unique([table1(:,1);table2(:,1)]);
table=zeros(length(A),9); %matrix preallocation
%Out in the first column the time intervals
table(:,1)=A; 
%Put in the columns 2 and 3 and in the proper rows the deaths and alive
%taken from table1 columns 2 and 3
[~, ia ib]=intersect(table1(:,1),A);
table(ib,2:3)=table1(ia,2:3);
%Put in the columns 4 and 5 and in the proper rows the deaths and alive
%taken from table2 columns 2 and 3
[~, ia ib]=intersect(table2(:,1),A);
table(ib,4:5)=table2(ia,2:3);
%remove the rows where there arent't deaths in both treatments
table((table(:,2)==0 & table(:,4)==0),:)=[];
clear A c ia ib table1 table2
%fill the "pigeon-holes"
c=find(table(:,3)==0); %find the "pigeon-holes" of treatment 1
for I=1:length(c)
    if c(I)~=1
        %find the first interval time before the hole where there is almost 1
        %death
        J=find(table(1:c(I)-1,3)>0,1,'last');
        table(c(I),3)=table(J,3)-table(J,2);
        if ~isempty(table12)
        %find eventually censored data
            K=find((table12(:,1)<table(c(I),1) & table12(:,1)>=table(J,1)),1,'last');
            %Put in the hole how many subject were alive before the interval time
            %of the hole
            if ~isempty(K)
                table(c(I),3)=table(c(I),3)-sum(table12(K,2));
            end
        end
    else
        table(1,3)=length(x1);
    end
end
%Do the same for tratment 2
c=find(table(:,5)==0);
for I=1:length(c)
    if c(I)~=1
        J=find(table(1:c(I)-1,5)>0,1,'last');
        table(c(I),5)=table(J,5)-table(J,4);
        if ~isempty(table22)
            K=find((table22(:,1)<table(c(I),1) & table22(:,1)>=table(J,1)),1,'last');
            if ~isempty(K)
                table(c(I),5)=table(c(I),5)-sum(table22(K,2));
            end
        end
    else
        table(1,5)=length(x2);
    end
end
clear c I J K table12 table22

%Fill the table and compute the statistic variable
%Compute the total deaths and alive before the i-th time interval
table(:,6:7)=[sum(table(:,[2 4]),2) sum(table(:,[3 5]),2)];
%Compute the difference between observed deaths for treatment 1 and
%expected deaths in the hyphthesis that the treatments are similar
table(:,8)=table(:,2)-table(:,3).*table(:,6)./table(:,7);
%Log-rank statistic is the sum of column 8 values
J=sum(table(:,8)); UL=abs(J);
%Compute the contribute to the standard error
table(:,9)=prod(table(:,[3 5 6]),2).*(table(:,7)-table(:,6)) ...
    ./(table(:,7).^2.*(table(:,7)-ones(size(table,1),1)));
%find if there is some NaN (i.e. 0/0)
loc=isnan(table(:,9));
if any(loc)
    table(loc,9)=0;
end
V=sum(table(:,9)); SUL=sqrt(V); %Compute the totale standard error
K=J/V; HR=exp(K); HRci=[exp(K-1.96/SUL) exp(K+1.96/SUL)];
z=abs((UL-0.5)/SUL); %normalized UL with Yates'es correction
p=2*(1-0.5*erfc(-z/realsqrt(2))); %p-value

%display results
% disp('LOG-RANK TEST FOR KAPLAN-MEIER SURVIVAL FUNCTIONS')
% disp(' ')
% tr=repmat('-',1,110);
% disp(tr)
% disp('HAZARD RATE IS AN EXPERIMENTAL FUNCTION!!!!')
% fprintf('Treatment 1: Hazard rate: %0.4f\n',lambda1)
% fprintf('Treatment 2: Hazard rate: %0.4f\n',lambda2)
% fprintf('\n')
% fprintf('Mantel-Haenszel Hazard ratio: %0.4f\n',HR)
% fprintf('95%% confidence interval: %0.4f - %0.4f\n',HRci)
% disp(tr)
% fprintf('UL\t\t\tS.E.\t\tz\t\tp-value (2-tailed test)\t\talpha\n')
% disp(tr)
% fprintf('%0.5f\t\t\t%0.5f\t\t%0.5f\t\t%0.5f\t\t\t\t%0.3f\n',UL,SUL,z,p,alpha)
% disp(tr)
% if p<alpha
%     fprintf('\t\tThe survival functions are statistically different\n')
%  else
%     fprintf('\t\tThe survival functions are not statistically different\n')
% end
