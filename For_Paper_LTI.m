% This code compares different measures of GMSC 
% for "Measures of Generalized Magnitude-Squared Coherence:
% Differences and Similarities" by Sheida Malekpour, John A. Gubner, and
% William A. Sethares. It has been submitted to a Journal for possible
% publication.
% Sheida Malekpour 10/9/2016
% This code compares different measures of GMSC (for our GMSC paper)
% Sheida Malekpour 10/9/2016
clc
clear all
close all

width = 7;     % Width in inches
height = 7;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 11;      % Fontsize
lw = 2.5;      % LineWidth
msz = 10;       % MarkerSize

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
K = 512/2;%512;
K2 = K/2;
wwf = [0:1/K2:1-1/K2]';
%GMSC=2;%1:Naive 2:Ramirez 3:Koopmans 4:Pauscaul
%Data information
% Using Welch method of window of win, with overlap of (window length - inc)

WL      = 1500;
win=rectwin(WL);
inc=808;

n       = 2*10240;         %number of samples

mu = [0,0,0,0];
sigma = [1,0.1,0.5,0.5;0.1,3,1,0.1;0.5,1,1,0.1;0.5,0.1,0.1,2];
rng default  % For reproducibility
ra = mvnrnd(mu,sigma,n);
r=1;
% data should be matrices of mxt and nxt, where m and n are number of
% signals and t is the number of samples
xx1(1:2)=0;
xx2(1:2)=0;
y1(1:2)=0;
y2(1:2)=0;
f1=0.05;%0.05;
f2=0.28;%0.28;
%recp1=pulstran(0:n,0:1/f1:n,'rectpulse');
%recp2=pulstran(0:n,0:1/f2:n,'rectpulse');
for ii=3:n
xx1(ii)=sin(2*pi*f1)*ra(ii-1,1)+2*cos(2*pi*f1)*xx1(ii-1)-xx1(ii-2);
xx2(ii)=sin(2*pi*f2)*ra(ii-1,2)+2*cos(2*pi*f2)*xx2(ii-1)-xx2(ii-2);
y1(ii)=sin(2*pi*2.5*f1)*ra(ii-1,3)+2*cos(2*pi*2.5*f1)*y1(ii-1)-y1(ii-2);
y2(ii)=sin(2*pi*f2)*ra(ii-1,4)+2*cos(2*pi*f2)*y2(ii-1)-y2(ii-2);
% xx1(ii)=sawtooth(2*pi*f1)*ra(ii-1,1)+2*sawtooth(2*pi*f1+pi/2)*xx1(ii-1)-xx1(ii-2);
% xx2(ii)=sawtooth(2*pi*f2)*ra(ii-1,2)+2*sawtooth(2*pi*f2+pi/2)*xx2(ii-1)-xx2(ii-2);
% y1(ii)=sawtooth(2*pi*2*f1)*ra(ii-1,3)+2*sawtooth(2*pi*2*f1+pi/2)*y1(ii-1)-y1(ii-2);
% y2(ii)=sawtooth(2*pi*f2)*ra(ii-1,4)+2*sawtooth(2*pi*f2+pi/2)*y2(ii-1)-y2(ii-2);
% % % xx1(ii)=sawtooth(2*pi*f1*ii)+r(ii-1,1);
% % xx2(ii)=sawtooth(2*pi*f2*ii)+r(ii-1,2);
% % y1(ii)=sawtooth(4*pi*f2*ii)+r(ii-1,3);
% % y2(ii)=sawtooth(2*pi*f1*ii)+r(ii-1,4);
% xx1(ii)=sawtooth(2*pi*f2+pi/2)*r(ii-1,1)+2*sawtooth(2*pi*f2+pi/2)*xx1(ii-1)-xx1(ii-2);
% xx2(ii)=sawtooth(2*pi*f2)*r(ii-1,2)+2*sawtooth(2*pi*f2+pi/2)*xx2(ii-1)-xx2(ii-2);
% y1(ii)=sawtooth(2*pi*2*f1)*r(ii-1,3)+2*sawtooth(2*pi*2*f1+pi/2)*y1(ii-1)-y1(ii-2);
% y2(ii)=sawtooth(2*pi*f2)*r(ii-1,4)+2*sawtooth(2*pi*f2+pi/2)*y2(ii-1)-y2(ii-2);
end
x=[xx1;xx2];
y=[y1;y2];
    
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
                plot(freq(length(freq)/2+1:end),ccfn(rr,length(freq)/2+1:end))%(65:end)
            end
        end
        ylabel('Pairwise');
        xlabel('Frequency');
        grid on;
        set(gca, 'FontSize', 20)
        close
        figure;
        hold on
        NGMSC=max(ccfn);
        NN=plot(freq(length(freq)/2+1:r:end),NGMSC(length(freq)/2+1:r:end))    
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
      %  figure;
        RM=fftshift(abs(Lambda));
        plot(freq(length(freq)/2+1:r:end),RM(length(freq)/2+1:r:end))
   
        hold on
        grid on;
        set(gca, 'FontSize', 20)
        
        %pauscaul
    elseif GMSC==3
        Block1=[x];
        Block2=[y];
        %figure;
        [Sall,coherencep,freq]=coherence_vector_general_correct3_window_overlap(Block1,Block2,win,inc,fs);
        PQ=fftshift(abs(coherencep));
        PP=plot(freq(length(freq)/2+1:r:end),PQ(length(freq)/2+1:r:end),'-')
        hold on
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
        max_ev=max(V_cohxy);
        %figure
        KP=fftshift(max_ev);
        KK=plot(freq(length(freq)/2+1:r:end),KP(length(freq)/2+1:r:end));
   
        grid on;
        set(gca, 'FontSize', 20)

    end
    end   
    ylabel ('Amplitude of the Measures (Unitless)')
    xlabel('Frequency (Hz)');
    legend('N-GMSC','R-GMSC','P-GMSC','K-GMSC','location','Best')
    legend boxoff
    print('../ex111', '-dpdf', '-r300');

    
    
    

