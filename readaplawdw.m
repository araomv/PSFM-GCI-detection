function [s,fs,egg,txt,phn]=readaplawdw(fn)
% Read a file from the aplawdw database
%
%  Inputs: fn   String containing the file name
%
% Outputs: s(n,1)   speech waveform
%          fs       sample frequency (always 20000)
%          egg(n,1) Laryngograph (a.k.a. EGG) waveform
%          txt      Text string of prompt
%          phn{k,2} Phonetic transcription in IPA symbols
%                   Each row contains {[t1 t2] 'unicode-string'}]
if fn(end-3)=='.'
    fn=fn(1:end-4);
end
[s,fs]=readwav([fn '.wav']);
if exist([fn '.egg'])
    egg=readwav([fn '.egg']);
else
    egg=[];
end
if exist([fn '.txt'])
    fid=fopen([fn '.txt'],'r');
    txtt=fscanf(fid,'%d%d');
    txt=fgetl(fid);
    fclose(fid);
else
    txt='';
end
if exist([fn '.phn'])
    fid=fopen([fn '.phn'],'r','n','UTF-8');
    npc=2;
    nph=0;
    phn=cell(npc,2);                    % space for annotations
    tline=fgetl(fid);
    while tline~=-1
        if numel(tline)>0
            tl1=double(tline(1));
            if tl1<0 || tl1>255
                tline(1)=[];            % delete the byte order mark if present
            end
            [pht,phcnt,phe,phnxt]=sscanf(tline,'%d%d',2);
            if phcnt==2 && length(tline)>phnxt
                nph=nph+1;              % incremeent phoneme count
                if nph>npc
                    npc=2*npc;
                    phn{npc,2}=[];
                end
                phn{nph,1}=(1+pht')/fs;          % first sample is time 1/fs
                phn{nph,2}=tline(phnxt+1:end);
            end
            tline=fgetl(fid);
        end
    end
    fclose(fid);
    if npc>nph
        phn(nph+1:end,:)=[];
    end
else
    phn=[];
end
if ~nargout
    clf;
    if numel(egg)>0
        subplot(3,1,3);
        spgrambw(egg,fs,'pJcw',20,[0 400]);
        h1=gca;
        subplot(3,1,1:2);
        if numel(phn)>0
            spgrambw(s,fs,'pJcwat',[],[],[],[],phn);
        else
            spgrambw(s,fs,'pJcw');
        end
        xlabel('');
        linkaxes([h1 gca],'x');
    else
        if numel(phn)>0
            spgrambw(s,fs,'pJcwat',[],[],[],[],phn);
        else
            spgrambw(s,fs,'pJcw');
        end
    end
    if numel(txt)>0
        title([fn(end-5:end) ': ' txt]);
    else
        title(fn(end-5:end));
    end
end
