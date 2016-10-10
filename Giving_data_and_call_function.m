clc
clear all
close all

width = 7;     % Width in inches
height = 7;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 11;      % Fontsize
lw = 1.5;      % LineWidth
msz = 8;       % MarkerSize

set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz

% Set the default Size for display
defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]);

% Set the defaults for saving/printing to a file
set(0,'defaultFigureInvertHardcopy','on'); % This is the default anyway
set(0,'defaultFigurePaperUnits','inches'); % This is the default anyway
defsize = get(gcf, 'PaperSize');
left = (defsize(1)- width)/2;
bottom = (defsize(2)- height)/2;
defsize = [left, bottom, width, height];
set(0, 'defaultFigurePaperPosition', defsize);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%Global variables
fs=1;
K = 256;
K2 = K/2;
wwf = [0:1/K2:1-1/K2]';
GMSC=2;%1:Naive 2:Ramirez 3:Koopmans 4:Pauscaul
%Data information
% Using wlch method of window of win, with overlap of length-inc
WL      = 100;
win=rectwin(WL);
inc=512/16;

n       = 1024;         %number of samples
nT      = [0:n-1]';     %time axis
Nf      = 2;            %number of frequencies of high coherence
f       = zeros(Nf,1);
f(1) = 0.15; f(2) = 0.28; %f(3) = 0.07; f(4) = 0.08; f(5) = 0.09;
%
fw      = 2*pi*f;
x1      = randn(n,1); %first1 signal
x2      = randn(n,1); %first2 signal
y1      = randn(n,1); %second signal

for i = 1:Nf
    x1  = x1 + cos(fw(i)*nT+ 2*pi*rand(1,1));
    x2  = x2 + cos(2*fw(i)*nT+ 2*pi*rand(1,1));
    y1  = y1 + cos(fw(i)*nT + 2*pi*rand(1,1));
end
% data should be matrices of mxt and nxt, where m and n are number of
% signals and t is the number of samples
xx1=x1';
xx2=x2';
x=[xx1;xx2];
y1=y1';
y=[y1];
    
    %%
    Ball=[x;y];
    V_cohxy=[];
    V_cohyx=[];
    eigv=[];
    Lambda_inov=[];
    M=size(Ball,1);

    
    %%%%%%
    for GMSC=1:4
    if GMSC==1
        Block1=[x];
        Block2=[y];
        s1=size(Block1,1);
        s2=size(Block2,1);
        M=s1+s2;
        %N-GMSC
         rr=0;
        for ii=1:s1
            for jj=1:s2
                rr=rr+1;
                [Sall,coherence,freq]=coherence_vector_general_correct3_window_overlap(Block1(ii,:),Block2(jj,:),win,inc,fs);
                hold on
                ccfn(rr,:)=fftshift(abs(coherence));
                plot(freq,(ccfn(rr,:)))
            end
        end
        ylabel('Pairwise');
        xlabel('Frequency');
        grid on;
        set(gca, 'FontSize', 20)
        figure;
        NGMSC=max(ccfn);
        plot(freq,(NGMSC))
        ylabel('N-GMSC');
        xlabel('Frequency');
        grid on;
        set(gca, 'FontSize', 20)
        
        
    %Ramirez
    elseif GMSC==2
        Block1=[x];
        Block2=[y];
        s1=size(Block1,1);
        s2=size(Block2,1);
        M=s1+s2;
        coh=zeros(M,M,512);
        for ii=1:M
            coh(ii,ii,:)=1;
        end
        Ball=[Block1;Block2];
        coh=zeros(M,M,2^nextpow2(WL)/2^6);
        for ii=1:M
            coh(ii,ii,:)=1;
        end
        for ii=1:M-1
            for jj=ii:M-1
                [Sall,coherence1,freq]=coherence_Ram_vector_general_correct3_window_overlap(Ball(ii,:),Ball(jj+1,:),win,inc,fs);
                coh(ii,jj+1,:)=coherence1;
                coh(jj+1,ii,:)=coherence1;
            end
        end
        
        for ii=1:size(coh,3)
            eigv(:,ii)=eig(coh(:,:,ii));
            Lambda(ii)=(max(abs(eigv(:,ii)))-1)/(M-1);
        end
        figure;
        plot(freq,fftshift(abs(Lambda)))
        ylabel('R-GMSC');
        xlabel('Frequency');
        grid on;
        set(gca, 'FontSize', 20)
        
        %pauscaul
    elseif GMSC==3
        Block1=[x];
        Block2=[y];
        figure;
        [Sall,coherencep,freq]=coherence_vector_general_correct3_window_overlap(Block1,Block2,win,inc,fs);
        plot(freq,fftshift(abs(coherencep)))
        ylabel('P-GMSC');
        xlabel('Frequency');
        grid on;
        set(gca, 'FontSize', 20)    
    %K-GMSC 
    elseif GMSC==4
        
        [coherencexy,coherenceyx,freq]=coherence_vector_Koopmans_v2(x,y,rectwin(WL),ceil(inc),fs);
        for ii=1:length(freq)
            if ndims(coherencexy)==2
                cohxy=coherencexy(:,ii);
            else
                cohxy=coherencexy(:,:,ii);
            end
            
            if ndims(coherenceyx)==2
                cohyx=coherenceyx(:,ii);
            else
                cohyx=coherenceyx(:,:,ii);
            end
            [U,V]=eig(cohxy);
            V_cohxy(:,ii)=abs(diag(V));
            [U,V]=eig(cohyx);
            V_cohyx(:,ii)=abs(diag(V));
        end
        trace_Kxy=sum(V_cohxy,1);
        trace_Kyx=sum(V_cohyx,1);
        max_ev=max(V_cohxy);%maximum eigenvalue
        figure
        plot(freq,fftshift(trace_Kxy));
        ylabel('trace of K-GMSC');
        xlabel('Frequency');
        grid on;
        set(gca, 'FontSize', 20)

    end
    end   