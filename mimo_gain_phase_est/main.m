%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2021 Hohai University.
% Author:	 Evans Baidoo
% Date:		 2020.03.12 
% Project Name: PET4MIMO Radar
% Module Name:	SPARSITY-AWARE FG NYSTRÖM FOR BISTATIC MIMO RADAR WITH ARRAY GAIN-PHASE ERROR
%
% Revision         : V4.0
% Additional Comments:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;  close all; clc; 


tic; % start a stopwatch timer 
 
%-----------------------------initializing--------------------------------% 
M = 8;     % the number of transmit array elements 
N = 8;     % the number of receive array elements 
d = 0.5;   % the distance between each array element 
P = 100;   % the sampling number/Pulse 
fs=1000;   % the sampling rate (Hz) 
SNR =10;   % the signal-to-noise ratio
dis_ang = 0.2;            % resolution interval
search_range = -60:dis_ang:60;  % the searching range of the signal parameter 

%% Case I
w=[pi/0.5 pi/2 pi/4].';               % doppler frequency  Case 1
DOD_ang = [-25 -11 28]        % the Directon of departure angle of the signal
DOA_ang = [-22 10 20]         % the Direction of arrival angle of the signal 

%Uncomment to run case two
% %% case II
% w=[pi/4 pi/4 pi/4].';               % doppler frequency    Case 11
% DOD_ang = [1 22 -25 ]        % the Direction of arrival angle of the signal 
% DOA_ang = [55 10 35 ]       % the Directon of departure angle of the signal

K = length(DOD_ang);     % the number of the target signals 

%Gain-phase formulation
sigma_a=0.5; 
sigma_phi=40; 
m = 3; n = 2;
alpha=(random('unif',0,1,[1,M-m])-0.5)*sigma_a*sqrt(12)+1; 
alphax=[1, 1, 1, alpha]; 
phi=deg2rad((random('unif',0,1,[1,M-3])-0.5)*sigma_phi*sqrt(12)); 
phix = [0, 0, 0, phi]; 
Gt=(alphax.*exp(j*phix))';              %Transmit Gain-phase 
alpha=(random('unif',0,1,[1,N-n])-0.5)*sigma_a*sqrt(12)+1; 
alphay=[1, 1, alpha] ; 
phi=deg2rad((random('unif',0,1,[1,N-n])-0.5)*sigma_phi*sqrt(12)); 
phiy= [0, 0, phi] ; ; 
Gr=(alphay.*exp(j*phiy))';              %Receive Gain-phase
Gt_hat=diag(Gt);
Gr_hat=diag(Gr);    

at = exp(-j*2*pi*d*(0:M-1).'*sin(DOD_ang*pi/180));  % transmit direction matrix
ar = exp(-j*2*pi*d*(0:N-1).'*sin(DOA_ang*pi/180));  % recieve direction matrix
A = khatriRao(Gt_hat*at,Gr_hat*ar);
estAA = khatriRao(at,ar);
[dicMat,dicMat2,dicDOD,dicDOA] = dic_S(search_range,M,N,d,Gt_hat,Gr_hat);%Dictionary matrix


amp =[1 1 1]';                      %Radar Cross Section
s=amp.*exp(j*w*[0:P-1]);            %source signal

S = 10.^(SNR/20)*s;
dd = rank(s);                           %rank of target signal
cgNoise = (randn(M*N,P)+j*randn(M*N,P))/sqrt(2);% accompanying noise
 X=A*S+cgNoise;  
 X1 = X*X'/P;

%Methods
% Forward backward Music approach with Gain-phase estimation
% [fbmgmus_T,fbmgmus_R] = FBMG_mus(X1,M,N,d,K,search_range,Gt,Gr,m,n); %FBMG
% [fbmgmus_DOD,fbmgmus_DOA] = realestI(fbmgmus_T,fbmgmus_R,DOD_ang, DOA_ang, K);

% Fast Greedy Nystrom with Gain-phase estimation
[fg_DOD,fg_DOA] = FG_somp(X,K,M,N,dicMat2,dicDOD,dicDOA,Gt_hat,Gr_hat,estAA,P);      %Proposed
% FG_somp1(search_range,X,K,M,N,d,Gt_hat,Gr_hat) 
[fgs_DOD,fgs_DOA] = realestI(fg_DOD,fg_DOA,DOD_ang, DOA_ang, K);


figure(3),plot(fbmgmus_DOD,fbmgmus_DOA,'bs','MarkerSize',10),hold on;
plot(fgs_DOD,fgs_DOA,'ro','MarkerSize',10),hold on;
      plot(DOD_ang,DOA_ang,'kx','MarkerSize',8),hold on;
    xlabel('DOD (degree)'); ylabel('DOA (degree)');
    legend({'FBMG','Proposed','Ground Truth'},'Location','best')






          