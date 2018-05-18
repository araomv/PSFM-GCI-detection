clear;
%% parameters
is_noisy=0; %enable noisy/clean
fs=20000;
Tw=40e-3;
pre_emph=0.99;
Q=2000; %% for noise estimation
psi=1e-7;%sigma_p
ga=0.1;%gamma
%% reading ground truth
filename='as02e0.wav';
Nw = round(Tw*fs);
clc;
%% prepare the reference GCIs and voice segments
[x,fs,lx,txt,phn]=readaplawdw([filename]);
propagation_time = 19;   % in APLAWD is 19, in SAM is 14.
lx = [zeros(propagation_time-1,1); lx];
DEGG = filter([1 -1],1,lx);
%%%%%%%%%%%%%%%%%%%%%%% Reference Glottal parameters %%%%%%%%%%%%%%%%%%%%
[rGCIs,rGOIs] = v_sigma(lx,fs); % Extract the reference GCIs/GOIs
% remove short voiced spurts
[rGCIs,rGOIs] = RemoveShortVoicedGCIsGOIs(rGCIs,rGOIs,fs);
% find the reference boundaries of VUD
MinP = 80;
[startss,finishess] = voiced_segments(rGCIs, rGOIs, fs, MinP);
rGCIs(DEGG(rGCIs)<0)=[];
rGCIs = rGCIs(5:end-5);
disp('..loading data');
if(is_noisy) %% load noisy data and covariance estimation
    load('as02e0_20dB_destops');
    x=xnoisy;
    noise = x(1:Q);
    noise = filter([1 -pre_emph],1,noise);
    y1 =enframe(noise,rectwin(Nw),1);
    y1 = y1-repmat(mean(y1),size(y1,1),1);
    ss = y1'*y1/size(y1,1);
end
%% run the algorithm
[ gci,c_est,a_est,b_est,sn_est,b_est_mode ] = run_PSFM(x,fs,Nw,is_noisy,ss,pre_emph,psi,ga);
rGCI=rGCIs;
%% error analysis
[~,~,~,~,~,~,~,~,~,~,~,~,~,loc_fas_GCI1,loc_missed_GCI1] = ...
    Evaluate_GCI_GOI( x,rGCI,rGCI,gci,gci+10,MinP,fs );%dummy GOIs
%% plot results
clf;
disp('..plottig results');
xpa = zeros(gci(end),1);
xpa(gci)=1;

xgt = zeros(rGCI(end),1);
xgt(rGCI)=1.05;
a1=subplot(511);
plot(0:1/fs:(length(b_est(:))-1)/fs,b_est);
h=legend('$\rho_k^*$');set(h,'FontSize',10,'FontWeight','bold','interpreter','latex');

a4=subplot(512);
plot(0:1/fs:(length(b_est(:))-1)/fs,a_est(:));
h=legend('$a_k^*$');set(h,'FontSize',10,'FontWeight','bold','interpreter','latex');

a2=subplot(513);
stem(0:1/fs:(length(xpa)-1)/fs,xpa,'r');hold on;stem(0:1/fs:(length(xgt)-1)/fs,xgt,'g');hold off;

hold on;
scatter(loc_fas_GCI1/fs,1.2*ones(length(loc_fas_GCI1),1),80,'r','filled')
scatter(loc_missed_GCI1/fs,1.2*ones(length(loc_missed_GCI1),1),80,'g','filled')
hold off;
h=legend('estimated locations','ground truth locations','false alrams','missed','Orientation','Horizontal','Location','NorthOutside');set(h,'FontSize',10,'FontWeight','bold');


a5=subplot(514);
plot(0:1/fs:(length(DEGG)-1)/fs,DEGG);
h=legend('DEGG');set(h,'FontSize',10,'FontWeight','bold');

a3=subplot(515);
plot(0:1/fs:(length(x)-1)/fs,x);
h=legend('Signal');set(h,'FontSize',10,'FontWeight','bold');
axis tight;
linkaxes([a1,a2,a3,a4,a5],'x');
axis([1.56,1.65,ylim]);