function c=cindex(beta,feature,time,status)
% >>c=cindex(model,w)
% >>c=cindex(u,y,delta)

% if nargin>2
%     u=model;
%     y=w;
% else
%     [~,x,y,delta]=unfold(model);
%     u=x*w;
% end
risk = feature*beta;
u=risk;
u = exp(u);
y=time;
delta=status;

n=length(u);%risk

a=repmat(u,n,1);
b=reshape((reshape(u,n,1)*ones(1,n))',n*n,1);

us=(a-b>0)+0.5*(a-b==0);

a=repmat(y,n,1);%time
b=reshape((reshape(y,n,1)*ones(1,n))',n*n,1);
ys=(a-b<0);

d1=repmat(delta,n,1);%stats
d2=reshape((reshape(delta,n,1)*ones(1,n))',n*n,1);

comp=(ys==1 & d1==0);
conc=comp.*(us);

c=sum(conc)/sum(comp);