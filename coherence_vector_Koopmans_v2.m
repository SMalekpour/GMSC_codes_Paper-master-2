function [coherenceXYy,coherenceYXx,freq]=coherence_vector_Koopmans_v2(x,y,win,inc,Fs)
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

for kk=1:size(fftZ,2)
    SXX=SzzM(1:sx,1:sx,kk);
    SYY=SzzM(sx+1:end,sx+1:end,kk);
    SXY=SzzM(1:sx,sx+1:end,kk);
    SYX=SzzM(sx+1:end,1:sx,kk);
    coherenceYX(:,:,kk)=sqrtm(inv(SYY))*SYX*inv(SXX)*SXY*sqrtm(inv(SYY));
    coherenceXY(:,:,kk)=sqrtm(inv(SXX))*SXY*inv(SYY)*SYX*sqrtm(inv(SXX));
end
if size(y,1)==1
    coherenceYXx=squeeze(coherenceYX)';
else
    coherenceYXx=squeeze(coherenceYX);
end
if size(x,1)==1
    coherenceXYy=squeeze(coherenceXY)';
else
    coherenceXYy=squeeze(coherenceXY);
end

 freq=-fs/2:fs/size(fftZ,2):fs/2-fs/size(fftZ,2);

