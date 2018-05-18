function [c_est,a_est,b_est_mean,sn_est,b_est_mode,res] = GCI_detection_Samp_withcovar2(z,aord,init,lam,sa,alpha,beta,niter,Nbur,Sigma)
%%
N=length(z);
% sn=sigma_n^2;
a=zeros(N,1);
b=zeros(N,1);
s =[init; z];
S=zeros(N,aord);
for j=aord:length(s)-1
    S(j-aord+1,:)=s(j:-1:j-aord+1)';
end
s = s(aord+1:end);
c=-pinv(S)*s;
res=S*c+s;
%% conditional probablity compute functions (note sn is sigma_n^2)
so = @(sa,sn) (1/sn)+(1/sa);
%%
sn=1e-4;
SigInv=inv(Sigma);
Sni = SigInv/sn;%sn_sa =(1/sn)+(1/sa);
sc=1e+5;
best=zeros(N,niter);
aest=zeros(N,niter);
snest=zeros(1,niter);
cest=zeros(length(c),niter);
for k=1:niter
    d =s+S*c;
    B=diag(b);
    K = B*Sni*B;
    ttt=Sni*d;
%     plot(d);hold on;plot(L*d,'r');hold off;
    j=1;
%     plot(ttt);hold on;plot(d,'r');hold off;
    while(j<=length(z))
        b(j)=1;
        
        K(j,b>0)=Sni(j,b>0);
        K(b>0,j)=Sni(b>0,j);
        sn_sa = (K(j,j))+(1/sa);
        val=(K(j,:)*a-K(j,j)*a(j)-ttt(j));
        p1=(1/(sa*sqrt(sn_sa)))*exp(val^2/(2*sn_sa));
        p1(p1>20e+5)=20e+5;
        p = p1*lam/(p1*lam+1-lam);
        if(rand(1)<p)
            a(j)=normrnd(-val/sn_sa,sqrt(1/sn_sa));
        else
            a(j)=0;
            b(j)=0;
        end
        j=j+1;
    end
    d2=s-a.*b;
    SS=S'*Sni*S;
    Sc = (inv((SS+(eye(size(S,2))/sc))));Sc=(Sc+Sc')/2;
    mu = -Sc*S'*Sni*d2;
    c = mvnrnd(mu,Sc)';
    
    d = S*c+s-a.*b;
    
    val = d'*SigInv*d;
    sn = 1./gamrnd((N/2)+alpha,1./((.5*sum(val))+beta));
    Sni = SigInv/sn;
    best(:,k)=b;
    aest(:,k)=a;
    cest(:,k)=c;
    snest(k)=(sn);
end
c_est=[1;mean(cest(:,Nbur:end),2)];
b_est_mean = detectGCI(best,Nbur,'mean');
b_est_mode = detectGCI(best,Nbur,'mode');
a_est = mean(aest(:,Nbur:end),2);
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