function [c_est,a_est,b_est_mean,sn_est,b_est_mode,res,ll] = GCI_detection_Samp(z,aord,init,lam,sa,alpha,beta,niter,Nbur)
%% 
N=length(z);
b=binornd(1,lam,N,1);
a = normrnd(0,sa,N,1);
% sn=sigma_n^2;

s =[init; z];
S=[];
for j=aord:length(s)-1
    S=[S;s(j:-1:j-aord+1)'];
end
s = s(aord+1:end);
c=-pinv(S)*s;
c=zeros(size(c));%c(1)=1;
res=S*c+s;
%% conditional probablity compute functions (note sn is sigma_n^2)
so = @(sa,sn) (1/sn)+(1/sa);
lamda_m_1=@(j,d,sa,sn) (sqrt(1/(so(sa,sn))*sa))*exp(d(j)^2/(2*sn^2*so(sa,sn)));  %%changed may 5

mu_a = @(m,sa,sn,d) d(m)/(sn*so(sa,sn));
sig_a = @(sa,sn)  1/so(sa,sn);
Sinv=inv(S'*S+0.001*eye(size(S,2)));
mu_c = @(S,d) S\d;

Sig_c = @(S,sn) sn*Sinv;
%% 
sn=1e-4;

for k=1:niter
    d =S*c+s;
    for j=1:length(a)
        p1=lamda_m_1(j,d,sa,sn);
        if(p1==Inf)
            p=1;
        else
            p=p1*lam/(p1*lam+1-lam);
        end
        
        if(rand(1)<p)
            b(j)=1;
            a(j)=normrnd(mu_a(j,sa,sn,d),sig_a(sa,sn));
        else
            b(j)=0;
            a(j)=0;
        end
        
    end
    d=s-a.*b;
    c = mvnrnd(-mu_c(S,d),Sig_c(S,sn))';
    d=S*c+s-a.*b;
    sn = 1./gamrnd((N/2)+alpha,1/((.5*sum(d.^2))+beta));
   ll(k)= -norm(s+S*c-a.*b).^2;
    best(:,k)=b;
    aest(:,k)=a;
    cest(:,k)=c;
    snest(:,k)=sn;
end
c_est=[1;mean(cest(:,Nbur:end),2)];
b_est_mean = detectGCI(best,Nbur,'mean');
b_est_mode = detectGCI(best,Nbur,'mode');
a_est = mean(aest(:,Nbur:end).^2,2);
sn_est=mean(snest(Nbur:end));
end


function b_est=detectGCI(best,Nbur,mtd)

switch mtd 
case 'mean'
b_est = mean(best(:,Nbur:end),2);
case 'mode'
b_est=mode(best,2);
end
end