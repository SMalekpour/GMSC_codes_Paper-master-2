function [Sall,coherence,freq]=coherence_vector_general_correct3_window_overlap(x,y,win,inc,Fs)
%overlap = length-inc.
ll=size(x,2);
for lv=1:size(x,1)
xframe=enframe(x(lv,:),win,inc);
xx(lv,:)=reshape(xframe',1,numel(xframe));
end
for lv=1:size(y,1)
yframe=enframe(y(lv,:),win,inc);
yy(lv,:)=reshape(yframe',1,numel(yframe));
end
z=[xx;yy];
sx=size(xx,1);
sy=size(yy,1);
Sxx=[];Syy=[];Sxy=[];Syx=[];S=[];Sall=[];coherence=[];
nn=size(xframe,1);
lln=size(xframe,2);
NFFT=2^nextpow2(lln)/2^6;


fs=Fs;
for jj=1:nn
fftZ=(fft((z(:,1+(jj-1)*lln:jj*lln))',NFFT))';
lengthSig=size(fftZ,2);
for kk=1:size(fftZ,2)
PZZ(:,:,kk,jj)=1/lengthSig*fftZ(:,kk)*fftZ(:,kk)';
end
end

Szzz=mean(PZZ,4);

SzzM=Szzz;


%%%

for kk=1:size(fftZ,2)
    S(:,:,kk)=([SzzM(1:sx,1:sx,kk), zeros(sx,sy);zeros(sy,sx), SzzM(sx+1:end,sx+1:end,kk)]);
    Sall(:,:,kk)=SzzM(:,:,kk);
    coherence(:,kk)=(SzzM(1:sx,sx+1:end,kk))/sqrt(det(S(:,:,kk)));
end


 freq=-fs/2:fs/size(fftZ,2):fs/2-fs/size(fftZ,2);



