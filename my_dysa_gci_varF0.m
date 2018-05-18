function [gci,cost]=my_dysa_gci_varF0(gcic, nll, psi,fs,N0,dy_fxmax,dy_fxmin)
%% N-best DP with log(p(r) and smoothness of pitch is cost.
%DPGCI   Choose the best Glottal Closure Instances with Dynamic Programming
%   gci=dpgci(gcic, s(:), Ch, fnwav, fs) returns vectors of sample indices corresponding
%   to the instants of glottal closure in the speech signal s at sampling frequency fs Hz.
%
%   Inputs:
%   gcic    is a matrix whos first column are the glottal closure instance candidates and
%           the second column is 1 if the corresponding gci is derived from a zero crossing 
%           but zero if the gci is from a a projected zero crossing
%   s       is the speech signal - MUST be a column vector
%   Ch      the phase slope cost of every candidate (same as inverse probablity / prob.*soe^2)
%   fnwav   is the frobenious norm function of s 
%   fs      is the sampling frequncy
%
%   Outputs:
%   gci     is a vector of glottal closure instances chosen by the DP



%   Revision History: 
%   Bugs:  Constants are hardwired but defined in a structure like pv (defined in grpdelpv)
%         

% get algorithm parameters from voicebox()

% dy_fxmin=voicebox('dy_fxmin');        % min larynx frequency (Hz)
% dy_fxmax=voicebox('dy_fxmax');        % min larynx frequency (Hz)
% dy_nbest=voicebox('dy_nbest');        % Number of NBest paths to keep
dy_spitch=voicebox('dy_spitch');              % scale factor for pitch deviation cost
wproj=voicebox('dy_cproj');           % cost of projected candidate
% dy_cspurt=voicebox('dy_cspurt');           % cost of a talkspurt
dy_wpitch=voicebox('dy_wpitch');           % DP pitch weighting
% dy_wener=voicebox('dy_wener');           % DP energy weighting
% dy_wslope=voicebox('dy_wslope');           % DP group delay slope weighting
% dy_wxcorr=voicebox('dy_wxcorr');           % DP cross correlation weighting
% dy_preemph=voicebox('dy_preemph');
% dy_ewtaper=voicebox('dy_ewtaper');        % Prediction order of FrobNorm method  in seconds
% dy_ewlen=voicebox('dy_ewlen');        % windowlength of FrobNorm method  in seconds
% dy_ewdly=voicebox('dy_ewdly');        % shift for assymetric speech shape at start of voiced cycle
% dy_xwlen=voicebox('dy_xwlen');            % cross-correlation length for waveform similarity (sec)
% dy_fxminf=voicebox('dy_fxminf'); 
dy_nbest=5;
% s_used=filter([1 -exp(-2*pi*dy_preemph/fs)],1,s);
% fnwav=frobfun(s_used,dy_ewtaper*fs,dy_ewlen*fs,dy_ewdly*fs);

% saf=max([200,dy_xwlen*fs/2+1,fs/dy_fxminf]);
% sin=and(saf<gcic,gcic<length(s)-saf);
% gcic=gcic(sin,:);
% prob=prob(sin);

%Constants
Ncand=length(gcic);
sv2i=-(2*psi^2)^(-1);                 % scale factor for pitch deviation cost
% nxc=ceil(dy_xwlen*fs);                      % cross correlation window length in samples
% === should delete any gci's that are within this of the end.
% === and for energy window

%Limit the search:
% dy_fxmax= 300;
% dy_fxmin = 50;
qrmin=ceil(fs/dy_fxmax);
qrmax=floor(fs/dy_fxmin);

%Cost and tracking r = current, q = previous, p = preprevious
cost=zeros(Ncand, dy_nbest); cost(:,:)=inf;    %Cost matrix, one row for each candidate
maxcost=zeros(Ncand,1); maxcost(:,:)=inf;   %Maximum cost in each row
imaxcost=ones(Ncand,1);                     %Index of maximum cost

prev = ones(Ncand, dy_nbest);                  %index of previous, q candidates
ind = ones(Ncand, dy_nbest);                   %index of p in row q (from prev)
qbest = [zeros(Ncand,1), ones(Ncand,2)]; % the minimum cost in any previous q [cost,q,i]

% Cfn=fnrg(gcic(:,1),fnwav,fs);  %Frob.Energy Cost

%Add start and end state
% === should probably delete candidates that are too close to either end of the input
% === why do we ever need the additional one at the tail end ?
gcic=[gcic(1)-qrmax-2 gcic gcic(end)+qrmax+2];
% Cfn=[0 Cfn 0];
nll = [0 nll 0];

% first do parallelized version


% rather complicated window specification is for compatibility with DYPSA 2
% === +1 below is for compatibility - probably a bug
% wavix=(-floor(nxc/2):floor(nxc/2)+1)';                 % indexes for segments [nx2,1]
% nx2=length(wavix);
% sqnx2=sqrt(nx2);
% g_cr=dy_wener*Cfn+dy_wslope*Ch+wproj*(1-gcic(:,2))';  % fixed costs
g_cr = log(sqrt(2*pi*psi));
g_n=gcic';                  % gci sample number [1,Ncand+2]
% g_pr=gcic(:,2)';                 % unprojected flag [1,Ncand+2]
% g_sqm=zeros(1,Ncand+1);         % stores: sqrt(nx2) * mean for speech similarity waveform
% g_sd=zeros(1,Ncand+1);         % stores: 1/(Std deviation * sqrt(nx2)) for speech similarity waveform
f_pq=zeros((Ncand+1)*dy_nbest,1);   % (q-p) period for each node
f_c=repmat(Inf,(Ncand+1)*dy_nbest,1);    % cumulative cost for each node - initialise to inf
f_c(1)=0;                       % initial cost of zero for starting node
f_costs=zeros(Ncand*dy_nbest,6);   % === debugging only remember costs of candidate
f_f=ones((Ncand+1)*dy_nbest,1);    % previous node in path
f_fb=ones((Ncand+1),1);    % points back to best end-of-spurt node
fbestc=0.5;                       % cost of best end-of-spurt node

qmin=2;
for r=2:Ncand+1   
%     if r==86
%         r;
%     end
    r_n=g_n(r);             % sample number of r = current candidate
    rix=dy_nbest*(r-1)+(1:dy_nbest);    % index range within node variables
    
    % determine the range of feasible q candidates
    qmin0=qmin;
    qmin=find(g_n(qmin0-1:r-1)<r_n-qrmax);      % qmin is the nearest candidate that is >qrmax away
    qmin=qmin(end)+qmin0-1;             % convert to absolute index of first viable candidate
    qmax=find(g_n(qmin-1:r-1)<=r_n-qrmin);      % qmax is the nearest candidate that is >=qrmin away
    qmax=qmax(end)+qmin-2;
%     plot(s);hold on;plot(gcic(1:end-1),s(gcic(1:end-1)),'*r');
%     scatter(r_n,s(r_n),'g','filled');
%     hold off;
    
    % calculate waveform similarity cost measure statistics
    
%     sr=s(r_n+wavix);        % note s MUST be a column vector so sr is also
%     wsum=sum(sr);
%     g_sqm(r)=wsum/sqnx2;                % mean * sqrt(nx2)
%     g_sd(r)=1/sqrt(sr.'*sr-wsum^2/nx2);   % 1/(Std deviation * sqrt(nx2))
    
    % now process the candidates
    
    if qmin<=qmax
%         disp([r qmin qmax g_n(r) g_n(qmin) g_n(qmax)]);
        qix=qmin:qmax;      % q index
        nq=length(qix);
%         hold on;scatter(g_n(qix),s(g_n(qix)),'y','filled');hold off;
        % === should integrate the -0.5 into dy_wxcorr
        % === the factor (nx2-1)/(nx2-2) is to compensate for a bug in swsc()
%         q_cas=-0.5*(nx2-1)/(nx2-2)*dy_wxcorr*(sum(s(repmat(g_n(qix),nx2,1)+repmat(wavix,1,nq)).*repmat(sr,1,nq),1)-g_sqm(qix)*g_sqm(r)).*g_sd(qix)*g_sd(r);
        q_cas=zeros(length(qix),1);
        % compare: i=35; Ca=swsc(g_n(qix(i)),g_n(r),s,fs); [i qix(i) r  g_n(qix(i)) g_n(r) dy_wxcorr*Ca q_cas(i)]
        
        % now calculate pitch deviation cost
        
        fix = 1+(qmin-1)*dy_nbest:qmax*dy_nbest;    % node index range
        f_qr=repmat((r_n-g_n(qix))',dy_nbest,1);    % (r-p) period for each node    difference in pitch computation
        f_pr=f_qr(:)+f_pq(fix);
        % === could absorb the 2 into sv2i
        f_nx=2-2*f_pr./(f_pr+abs(f_qr(:)-f_pq(fix)));
        dy_wpitch = 0.5;
        f_cp=dy_wpitch*(0.5-exp(sv2i*f_nx.^2));  %sv2i=-(2*psi^2)^(-1);  
        
%         
%         jj = f_qr(:);
%         jj_ii = f_pq(fix);
%          f_cp =.5*(f_nx.^2)./psi;
        
        f_nx = (f_qr(:)+f_pq(fix))/fs;
        f_cp =g_cr+(.5*(f_nx.^2)./psi);  %% It works?
        
        
        f_nx = (N0(r_n)-f_qr(:))/fs;
        f_cp =g_cr+(.5*(f_nx.^2)./psi);   %f0 specific cost
        
        
        
        % === fudge to match dypsa2.4 - could more efficiently be added
        % === onto the cost of a talkspurt end
        % === should be a voicebox parameter anyway
        dy_cspurt = -0.22;   %% cost when the period is zero ???
        f_cp(f_pq(fix)==0)=dy_cspurt;                                                       %% why is this??
        
        % now find the N-best paths
        
        [r_cnb,nbix]=sort(f_c(fix)+f_cp+reshape(repmat(q_cas,dy_nbest,1),nq*dy_nbest,1));
        f_c(rix)=r_cnb(1:dy_nbest)+nll(r);     % costs (qcr is the fixed cost)
        f_f(rix)=nbix(1:dy_nbest)+(qmin-1)*dy_nbest;       % traceback nodes
        f_pq(rix)=f_qr(nbix(1:dy_nbest));       % previous period
        % === f_costs is only for debugging
%         r;
        f_costs(rix,1)=f_c(fix(nbix(1:dy_nbest)));
        f_costs(rix,2)=reshape(q_cas(1+floor((nbix(1:dy_nbest)-1)/dy_nbest)),dy_nbest,1);
        
        % check cost of using this candidate as the start of a new spurt
        % ==== the qmin>2 condition is for compatibility with dypsa 2 and
        % prevents any spurts starting until at least qrmax past the first
        % gci. This is probably a bug (see again below)
        iNb=rix(end);        
        if (qmin>2) && (f_c(f_fb(qmin-1)))<f_c(iNb)        % compare with worst of Nbest paths
            f_f(iNb)=f_fb(qmin-1);
            % === for now we exclude the energy and phase-slope costs for compatibility with dypsa2
            % === this is probably a bug
            f_c(iNb)=f_c(f_fb(qmin-1))+wproj*(0);     % replace worst of the costs with start voicespurt cost
            f_pq(iNb)=0;                    % false pq period
        end
        if f_c(rix(1))<fbestc
            f_fb(r)=rix(1);                          % points to the node with lowest end-of-spurt cost
            % === should compensate for the pitch period cost incurred at the start of the next spurt
            % === note that a node can never be a one-node voicespurt on its own unless dy_nbest=1
            % since the start voices[purt option replaced the worst Nbest cost. This is probably good but
            % is a bit inconsistent.
            fbestc=f_c(rix(1));
        else
            f_fb(r)=f_fb(r-1);
        end
    else            % no viable candidates - must be the start of a voicespurt if anything
        % === for now we exclude the energy and phase-slope costs for compatibility with dypsa2
        % === this is probably a bug
        % ==== the qmin>2 condition is for compatibility with dypsa 2 and
        % prevents any spurts starting until at least qrmax past the first
        % gci. This is probably a bug (see again above)
        if (qmin>2)
            f_c(rix(1))=f_c(f_fb(qmin-1))+wproj*(0);  % cost of new voicespurt
            f_f(rix)=f_fb(qmin-1);                              % traceback to previous talkspurt end
            f_pq(rix)=0;                                        % previous period
        end
        f_fb(r)=f_fb(r-1);                                  % cannot be the end of a voicespurt
    end
end

% now do the traceback

gci = zeros(1,Ncand+1);

% === for compatibility with dypsa2, we force the penultimate candidate to be accepted
% === should be: i=f_fb(Ncand+1) but instead we pick the best of the penultimate candidate
i=rix(1)-dy_nbest;
if f_c(i-dy_nbest+1)<f_c(i)     % check if start of a talkspurt
    i=i-dy_nbest+1;
end
k=1;
cost=[];
while i>1
    j=1+floor((i-1)/dy_nbest);          % convert node number to candidate number
    gci(k)=g_n(j);
    cost(gci(k))=f_c(j);
    i=f_f(i);
    k=k+1;
end
gci=gci(k-1:-1:1);           % put into ascending order 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [frob]=frobfun(sp,p,m,offset)

% [frob]=frobfun(sp,p,m)
% 
% sp is the speech signal assumed to be preemphasised
% p  is the prediction order  : recomended to be 1 ms in above paper
% m  is the window length     : recomended to be 1 ms in above paper
% offset is shift for assymetric speech shape at start of voiced cycle -
% default 1.5ms.
%
% This function implements the frobenius norm based measure C defined in [4] below.
% It equals the square of the Frobenius norm of the m by p+1 data matrix divided by p+1
%
% Reference:
%   [4]  C. Ma, Y. Kamp, and L. F. Willems, �A Frobenius norm approach to glottal closure detection
%        from the speech signal,� IEEE Trans. Speech Audio Processing, vol. 2, pp. 258�265, Apr. 1994.


%   Revision History: 
%       V1.0 July 12th 2002:
%            Nov. 6th 2002 : if statement added to remove "last" midpoint 

%force p m and offset to be integers
p=round(p);
m=round(m);
offset=round(offset);

w=(p+1)*ones(1,m+p);
w(1:p)=1:p;
w(m+1:p+m)=p:-1:1;

w=w./(p+1); 
frob=filter(w,1,sp.^2);
frob(1:(round((p+m-1)/2) + offset))=[];


function Cfn=fnrg(gcic,frob,fs)

%Frobenious Energy Cost

dy_fxminf=voicebox('dy_fxminf');
frob=frob(:)';
mm=round(fs/dy_fxminf);
mfrob=maxfilt(frob,1,mm);
mfrob=[mfrob(floor(mm/2)+1:end) max(frob(end-ceil(mm/2):end))*ones(1,floor(mm/2))];
rfr=frob./mfrob;
Cfn=0.5-rfr(round(gcic));