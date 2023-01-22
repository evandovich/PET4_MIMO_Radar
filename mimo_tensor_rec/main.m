%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2021 Hohai University.
% Author:	 	Evans Baidoo
% Date:		 	2020.10.9 
% Project Name: PET4MIMO Radar
% Module Name:	TENSOR RECONSTRUCTION BASED COMPRESSED SENSING FOR BISTATIC MIMO RADAR 
% Email:		ebaidoo2@hhu.edu.cn
%
% Revision         : V4.0
% Additional Comments:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc,clear all,close all;
SNR = 5;
real_dod = [-15 10 25 ]            % DOA
real_doa = [-30 -1 15 ]            % DOD
K = length(real_dod);              % targets number

%% MIMO radar setting
M = 10;                      % transmitters
N = 10;                      % receivers
d = 0.5;
Fd = [200 400 800]/20000;     % Doppler Hz
P = 10;                       % pluse number
PQ = 256;                     % waveform code length
angleRange = -60:1:60;
RA = [0:N-1]; % receive manifold
TA = [0:M-1]; % transmit manifold
At = exp(-j*2*pi*d*TA.'*sind(real_dod)); 
Ar = exp(-j*2*pi*d*RA.'*sind(real_doa)); 

     %% Noise effect  
     %Parameter 1
NR = -0.2;             %Noi_dB;
mcVec = 10.^((NR * (1+[0 : 1 : N - 1].')) / 20);
mcVec = mcVec .* (1 + unifrnd(-0.05, 0.05, N, 1));    %Generate NX1 array of random numbers chosen from the -0.05 to 0.05
mcVec = mcVec .* exp(1j * 2 * pi * rand(N, 1));
PP = toeplitz(mcVec, mcVec);

%% Signal and waveform
beta = ones(K,P);                                       % target amp
sig = beta.*exp(j*2*pi*Fd.'*[0:P-1]);                   % signal
SO = hadamard(PQ);           
ST = (1+j)/sqrt(2)*SO(1:M,:);                           %orthogonal Waveform ST*ST' = identity
   
             X2 = zeros(M,N,P);      
        
        %% Tensor signal model
        for L = 1:P
            temp_sig =  Ar*diag(sig(:,L))*At.'*ST;
            power_signal = sum(sum(abs(temp_sig).^2))/(N*PQ);
            temp_noise = sqrtm(PP)*(randn(N,PQ)+i*randn(N,PQ))/sqrt(2);
            power_noise = sum(sum(abs(temp_noise).^2))/(N*PQ);
            amp = sqrt(power_signal/power_noise)/(10^(SNR/20));
            echo = temp_sig + amp*temp_noise;
            %matched filtering
            for m = 1:M
                temp_x = echo*(1/sqrt(PQ))*ST(m,:)';
                X2(m,:,L) = temp_x;          %X(row,column,page) tensor representation
               
            end
        end
      
        [dicOD,dicOA,DOD_angs,DOA_angs] = dictmat(angleRange,d,TA,RA);
        TRC = TRHTD(X2,P,M,N);
 
    %% Methods
    [rdSOMP_DOD, rdSOMP_DOA, specDOD, specDOA] = RD_SOMP(TRC ,dicOD,dicOA,DOD_angs,DOA_angs,K); %RD SOMP
    
    [rdSBL_DOD, rdSBL_DOA,sblspec1,sblspec2] = RD_SBL(TRC ,dicOD,dicOA,DOD_angs,DOA_angs,K); % RD SBL
    
    %% Direction estimation with only received signals 
           ddd = size(X2,3);
            X = reshape(X2,ddd,[]).';                           
               if P == 1
%             spectrum = abs(Y'*steerGrid).^2;
              spectrum = sum(abs(X'*angleRange).^2);
                else
             spectrum = sum(abs(X'*dicOA).^2);
                end
            spectrum = spectrum/max(spectrum);
            
             if P == 1
%             spectrum = abs(Y'*steerGrid).^2;
              spectrum1 = sum(abs(X'*angleRange).^2);
                else
             spectrum1 = sum(abs(X'*dicOD).^2);
                end
            spectrum1 = spectrum1/max(spectrum1);
 
%    spectrum = splectrum/max(spectrum);
            figure; plot(DOA_angs, spectrum);
            hold on; plot(angleRange,sblspec2,'k');
            hold on; plot(angleRange, specDOA,'g'); 
            targetSpectrum = sum(abs(sig).^2,2);
            targetSpectrum = targetSpectrum / max(targetSpectrum);
            hold on; stem(real_doa, targetSpectrum, 'BaseValue', 0);
            xlabel('DOA(deg)'); 
            ylabel('Spectial spectrum');
            legend({'Without Denosing','Proposed TR-SBLKF','Proposed TR-SOMPKF','Ground Truth'},'Location','best')
            drawnow;
            
             figure; plot(DOD_angs, spectrum1);
            hold on; plot(-angleRange,sblspec1,'k');
            hold on; plot(-angleRange, specDOD,'g'); 
            targetSpectrum = sum(abs(sig).^2,2);
            targetSpectrum = targetSpectrum / max(targetSpectrum);
            hold on; stem(real_dod, targetSpectrum, 'BaseValue', 0);
            xlabel('DOD(deg)'); 
            ylabel('Spectial spectrum');
            legend({'Without Denosing','Proposed TR-SBLKF','Proposed TR-SOMPKF','Ground Truth'},'Location','best')
            drawnow;
    