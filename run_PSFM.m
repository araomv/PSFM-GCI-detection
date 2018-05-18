function [ gci,c_est,a_est,b_est,sn_est,b_est_mode ] = run_PSFM(x,fs,Nw,is_noisy,ss,pre_emph,psi,ga)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
F0min=60;
F0max=400;
[f0,~] = SRH_PitchTracking(x,fs,50,500);
time = 0:(10e-3)*fs:length(x);
f0 = medfilt1(f0,13);
req_time = 0:Nw:length(x);
f0_req = interp1(time,f0,req_time);
f0_req = max(f0_req,13);
if(~is_noisy)
    sa=1;
    beta=0;
    alpha=Nw/2;
    disp('..working on clean data');
    disp('..estimating a_k and p(b_k)');
    [c_est,a_est,b_est,sn_est,b_est_mode]=Compute_PA_step1(x,fs,Nw,'unc',[],f0_req,pre_emph,alpha,beta,sa);
else
    sa=1;
    beta=0;
    alpha=Nw/4;
    disp('..working on noisy data');
    disp('..estimating a_k and p(b_k)');
    [c_est,a_est,b_est,sn_est,b_est_mode]=Compute_PA_step1(x,fs,Nw,'withcov',ss,f0_req,pre_emph,alpha,beta,sa);
end
% b_est=medfilt1(b_est(:),5);
gcic=find(b_est(:));
b_est=b_est(:);
soe = abs(a_est(gcic)+1e-10);
nll = -log(b_est(gcic).*soe.^(ga));
f0_gcic = interp1(time/fs,f0,(0:(length(b_est)-1))/fs);
N0_gcic = fs./f0_gcic;
disp('..running N-best');
[gci,~]=my_dysa_gci_varF0(gcic',nll', psi,fs,N0_gcic,F0max,F0min);

end

