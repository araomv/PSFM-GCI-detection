function [c_est,a_est,b_est,sn_est,b_est_mode,ll]=Compute_PA_step1(z,fs,Nw,mtd,ss,f0,pre_emph,alpha,beta,sa)
rng('default'); %%for replicable results.
z=filter([1 -pre_emph],1,z);
frms=enframe(z,rectwin(Nw),Nw,'r');
aord = round(fs/1000)+2;

if(strcmpi(mtd,'unc'))  %%clean case
    for id=1:size(frms,1)
        N0 = fs/f0(id);
        ncy = Nw/N0;
        lamda=ncy/Nw;
        if(id==1)
            init = zeros(1,aord);
        else
            init=frms(id-1,end-aord+1:end);
        end
        [c_est(:,id),a_est(:,id),b_est(:,id),sn_est(:,id),b_est_mode(:,id),res,ll(id,:)] = GCI_detection_Samp((frms(id,:)'),aord,init',lamda,sa,alpha,beta,40,5);
        
    end
elseif(strcmpi(mtd,'withcov')) %%noisy case with ss is the covariance estimate of the noise.
   
    parfor id=1:size(frms,1)
        N0 = fs/f0(id);
        ncy = Nw/N0;
        lamda=ncy/Nw;
        if(id==1)
            init = zeros(1,aord);
        else
              init=frms(id-1,end-aord+1:end);
        end
        [c_est(:,id),a_est(:,id),b_est(:,id),sn_est(:,id),b_est_mode(:,id)] = GCI_detection_Samp_withcovar2(frms(id,:)',aord,init',lamda,sa,alpha,beta,10,3,ss);
        
    end
end
b_est = b_est(:);b_est=b_est(1:length(z));a_est = a_est(:);a_est=a_est(1:length(z));
end
