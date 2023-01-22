%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2021 Hohai University.
% Author:	 Evans Baidoo
% Date:		 2019.01.20 
% Project Name: PET4MIMO Radar
% Module Name:	ANGLE ESTIMATION APPROACH FOR BISTATIC MIMO RADAR WITH COHERENT SOURCES
%
% Revision         : V4.0
% Additional Comments:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;  close all; clc; 
%-----------------------------initializing--------------------------------% 
M = 11;     							% the number of transmit array elements 
N = 11;    								% the number of receive array elements 
P = 100;    							% the sampling number/Pulse 
w=[pi/4 pi/4 pi/4 pi/4 pi/4].';		% doppler frequency
l=((2*pi*3e8)/w(1)+(2*pi*3e8)/w(2))/2;   
d=0.5*l;								% the distance between each array element  

 DOD = [-20 -05 05 20 35] 				% the Directon of departure angle of the signal
 DOA =[-10 00 15 25 40]					% the Direction of arrival angle of the signal
 
 
K = length(DOD);						 % the number of the target signals 
theta = -50:0.01:50;					% Searching range
SNR=5; 

at = exp(-j*(0:M-1).'*d*2*pi*sin(DOD*pi/180)/l); 	% Transmitting Antenna elements
ar = exp(-j*(0:N-1).'*d*2*pi*sin(DOA*pi/180)/l); 	% Receiving Antenna elements
A = khatriRao(at,ar);

item = 100;								% Number of simulations/ trails
amp =[1 1 1 1 1]';						% Targets Amplitude
s=amp.*exp(j*w*[0:P-1]);				% Waveform
 
for item_num = 1:item
    disp(['SNR = ',num2str(SNR),' dB, ',num2str(item_num), ' # try : ']);  
ss = 10.^(SNR/20)*s;
dd = rank(ss);                           %rank of target signal
nt1 = (randn(M*N,P)+j*randn(M*N,P))/sqrt(2);%noise
 x=A*ss+nt1; 							% Target waveform
 R3=x*x'/P; 

 %%  Smoothing
LT=5;                                   		%No of subarray
LR=5;
n0=N-LT+1;                                   	%Size of each suarray
m0=M-LR+1;
Q = [15 25 35 40];                              %User defined parameter

[CEss_T3, CEss_R3,T2] = Nystrom_SS2(x, Q(3), P, LT, m0, n0, M, N, K);        %proposed CESS
% [ce3_DOD, ce3_DOA] = realestI(CEss_T3,CEss_R3,DOD, DOA, K);

% %   Simulation 1 & 2
    figure(1),plot(item_num,CEss_R3,'k*'),hold on;
    xlabel('Simulation times'); ylabel('DOA (degree)');

    figure(2),plot(item_num,CEss_T3,'k*'),hold on;
    xlabel('Simulation times'); ylabel('DOD (degree)');

end
CESS_R = sort(CEss_R3);CESS_T = sort(CEss_T3);
figure(3),plot(CESS_R,CESS_T,'bo','MarkerSize',10),hold on;
      plot(DOA,DOD,'rx','MarkerSize',8),hold on;
    xlabel('DOD (degree)'); ylabel('DOA (degree)');
    legend({'Estimated angle','Ground Truth angle'},'Location','best')
% 





 
 
