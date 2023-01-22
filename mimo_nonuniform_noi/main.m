%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2021 Hohai University.
% Author:	 	Evans Baidoo
% Date:		 	2020.06.20 
% Project Name: PET4MIMO Radar
% Module Name:	ENHANCED CAPON FOR BISTATIC MIMO RADAR IN UNKNOWN NONUNIFORM NOISE
% Email:		ebaidoo2@hhu.edu.cn
%
% Revision         : V4.0
% Additional Comments:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;
M1 = 6;                          		%Transmit antenna
N = 8;                          		%Receive antenna 
P = 200;                         		%Pulse number 
DOD = [-55 -45 -15 -10];                %Direction of departure 
DOA = [3 13 43 48];                     %Direction of arrival 
K = length(DOD);                		%Number of targets

%SNR variation 1
% SNR_axis = 05;                    	%Signal to noise ratio 

%% uncomment to run snr variation II  
%SNR variation II
SNR_axis = 01;                    		%Signal to noise ratio 

theta = -90:0.1:90;             		%Searching/Detection range
f0 = 10e7*linspace(3,3.5,K).';  		%Doppler frequency 
lambda = 3e8/f0(3);             
fs = 1*f0(3);                   		%Pulse repeat frquency/Sampling frequency
d = lambda/2;                   		% antenna interval 
M = M1;
at = exp(-j*(0:M-1).'*d*2*pi*sin(DOD*pi/180)); 
ar = exp(-j*(0:N-1).'*d*2*pi*sin(DOA*pi/180)); 
A = khatriRao(at,ar);


%% Nonuniform noise formulation
Ant = M*N;
noise_var =[10.5, 9.0, 10.0, 8.0,12.0, 8.5, 6.0, 10.0 ];% 
WNPR = max(noise_var)/min(noise_var);
noise_var1 = kron(eye(M),diag(noise_var));
noise_var11 = diag(noise_var1);
signal_power_inv = diag(1./noise_var1);             
for ii = 1:Ant
T_noise = sqrt(noise_var11(ii))*(randn(1,P)+1i*randn(1,P))/sqrt(2); %Nonuniform noise
U_noise(ii,:) = T_noise;
end
        
        signal_power = Ant*10^(0.1*SNR_axis)/sum(signal_power_inv);%Signal power
        signal = sqrt(signal_power)*exp(1i*2*pi*f0*(0:P-1)/fs);
        Yhat = A*signal+U_noise;                    %Measurement signal
        X1 =Yhat*Yhat'/P;                             % array Covariance matrix
        
[X2]=Noisecor(X1,M,N,Ant);
[capDOD,capDOA,capondod,capondoa]=RD_CAPON1(X1,theta,M,N,K,d);
[CB_DOD,CB_DOA,CB,CB1]=CB_CAPON1(X2,theta,M,N,K,d,P);
[ CB_DOD1, CB_DOA1] = realestII(CB_DOD,CB_DOA,DOD, DOA, K);
[cpbDOD,cpbDOA, Cspecdod,Cspecdoa]=CPB_CAPON(X2,theta,M,N,K,d);
[cpbDOD1, cpbDOA1] = realestII(cpbDOD,cpbDOA,DOD, DOA, K);
[bpcDOD,bpcDOA, Bspecdod1,Bspecdoa1]=PCC_CAPON(X2,theta,M,N,K,d);

allPerm = perms([1:K]);
                    indexValue = zeros(size(allPerm, 1), 1);
                    for idxPerm = 1 : 1 : size(allPerm, 1)
                        indexValue(idxPerm) = sum(abs(capDOD(allPerm(idxPerm, :)) - DOD));
                    end
                    [minValue, minIndex] = min(indexValue);
                    capDOD1 = capDOD(allPerm(minIndex, :));
                      %                     
                    for idxPerm = 1 : 1 : size(allPerm, 1)
                        indexValue(idxPerm) = sum(abs(capDOA(allPerm(idxPerm, :)) - DOA));
                    end
                    [minValue, minIndex] = min(indexValue);
                    capDOA1 = capDOA(allPerm(minIndex, :));

amp = ([400 400 400 400]);

                    indexValue = zeros(size(allPerm, 1), 1);
                    for idxPerm = 1 : 1 : size(allPerm, 1)
                        indexValue(idxPerm) = sum(abs(bpcDOD(allPerm(idxPerm, :)) - DOD));
                    end
                    [minValue, minIndex] = min(indexValue);
                    bpcDOD1 = bpcDOD(allPerm(minIndex, :));
                      %                     
                    for idxPerm = 1 : 1 : size(allPerm, 1)
                        indexValue(idxPerm) = sum(abs(bpcDOA(allPerm(idxPerm, :)) - DOA));
                    end
                    [minValue, minIndex] = min(indexValue);
                    bpcDOA1 = bpcDOA(allPerm(minIndex, :));

    figure (1)
plot (theta,capondod,'r',theta,CB,'m',theta,Cspecdod,'g',theta,Bspecdod1,'b'); hold on;
plot(DOD, amp, 'ko'); hold on;
xlabel('DOD (degree)'); ylabel('Amplitude(dB)');
legend({'PCCB','CCB','PCB', 'RD Capon','Target'},'Location','best')
% figure (4); plot (theta,capon3b,'b'); 
figure (2)
plot (theta,capondoa,'r',theta,CB1,'m',theta,Cspecdoa,'g',theta,Bspecdoa1,'b'); hold on;
plot(DOA, amp, 'ko'); hold on;
xlabel('DOA (degree)'); ylabel('Amplitude(dB)');
legend({'PCCB','CCB','PCB', 'RD Capon','Target'},'Location','best')



figure(3),plot(bpcDOD1,bpcDOA1,'bo','MarkerSize',10),hold on;
plot(CB_DOD1,CB_DOA1,'mo','MarkerSize',10),hold on;
plot(cpbDOD1,cpbDOA1,'g+','MarkerSize',10),hold on;
plot(capDOD1,capDOA1,'r+','MarkerSize',10),hold on;
      plot(DOD,DOA,'ko','MarkerSize',8),hold on;
    xlabel('DOD (degree)'); ylabel('DOA (degree)');
    legend({'PCCB','CCB','PCB', 'RD Capon','Target'},'Location','best')
    