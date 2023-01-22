%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2021 Hohai University.
% Author:	 	Evans Baidoo
% Date:		 	2019.06.20 
% Project Name: PET4MIMO Radar
% Module Name:	HERMITIAN TRANSFORM-BASED DIFFERENCING FOR BISTATIC MIMO RADAR UNDER UNKNOWN NOISE FIELD
% Email:		ebaidoo2@hhu.edu.cn
%
% Revision         : V4.0
% Additional Comments:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc,clear all; close all;warning off all;
M = 8;                       % transmit number
N = 10;                      % receiver number
d = 0.5;

DOD = [-10 -25 5 -1 35 15]
DOA = [-25 -5 10 25 38 58]

K = length(DOD);            % target number
beta = ones(1,K);           % RCS coefficients
beta1 = ones(1);            % RCS coefficients
f1 = [200 300 500 650 750 950 ]/2000;     % Dopplor frequency
f2 = [400 850 ]/2000;
L = 100;					% Number of pulses
SNR = 10;
aphfa1 = 1;
amp1 = sqrt(aphfa1*10^(SNR/10)); 

S1 = amp1*(beta.'*ones(1,L)).*exp(j*2*pi*f1.'*[1:L]);
Ar = exp(-j*2*pi*d*(0:N-1).'*sin(DOD*pi/180));  % receive direction matrix
At = exp(-j*2*pi*d*(0:M-1).'*sin(DOA*pi/180));  % transmit direction matrix
A = khatriRao(At,Ar);

% Colored Noise Modelling
%% AR model /Parameterisation I
% nt = randn(M*N,L)+j*randn(M*N,L)/sqrt(2);       %white noise
 AR = [1 -1 0.2];
y = filter(1,AR,randn(20,N));       
Rs = arcov(y,N-1);

%% Toeplitz top row coefficient noise formulation
% UU = [1,0.9^1,0.9^2,0.9^3,0.9^4,0.9^5,0.9^6,0.9^7,0.9^8,0.9^9];             
% UU1 = toeplitz(UU,UU);
Rs1 = kron(eye(M),Rs);  
nt1 = Rs1*(randn(M*N,L)+j*randn(M*N,L))/sqrt(2);    %Colored Noise
xx = A*S1+nt1;
y = xx*xx'/L;

%%  ERD spectral
ut = [3]; ur = [3];
M1 = M-2*(length(ut)-1);
N1 = N-2*(length(ur)-1);
P = length(ut)-1;
Pt = [zeros(M1,P),eye(M1),zeros(M1,P)];
Pr = [zeros(N1,P),eye(N1),zeros(N1,P)];
Ptr = kron(Pt,Pr);
JMN = kron(eye(M),fliplr(eye(N)));
y1 = (JMN*y*JMN)-j*(JMN*y*JMN);
y2 = y - JMN*conj(y)*JMN;     
Y = Ptr*y1;
Y1 = Ptr*y2;
yy = Y*Y';
yy1 = y2*y2';
Qr = length(yy);
II = fliplr(eye(Qr));
Rx = yy+II*conj(yy)*II;
[U,S,V]=svd(Rx); 
[U2,S2,V2]=svd(yy1); 
D = U(:,K+1:Qr);
D2a = U2(:,K+1:Qr);
for ii=K+1:Qr
S(ii,ii) = 0;
S2(ii,ii) = 0;
Gs = S;
Gs2 = S2;
end
Gv = conj(V);
Gv2 = conj(V2);
Rxx = U*Gs*Gv;
Rxx2 = U2*Gs2*Gv2;
[U1,S1,V1]=svd(Rxx);
[U3,S3,V3]=svd(Rxx2);
D1 = U1(:,K+1:end);
D3 = U3(:,K+1:end);
D2 = (D1*D');
D4 = (D3*D2a');
z = 0.5;            %const
G11 = diag(diag(D2).^z);
G12 = diag(diag(D4).^z);
Gn1 = D2*G11*D2';
Gn4 = D4*G12*D4';
Gn2 = D2*D2';
Inr = inv(Rx);
Gnn = Gn4; 0.75;

YY1 = xx(:,1:L-1);
YY2 = xx(:,2:end);                                           
RYY = YY2*YY1'/L;

Z1=zeros((M*N)-N);
for i=1:((M*N)-N)-1
    Z1(i+1,i)=1;
end
x1=xx(1:(M*N)-N,:);
x2=xx(N+1:M*N,:);%
R1=x1*x1'/L;%
R2=x1*x2'/L;%
D1=eig(R1);
u2=mean(D1(1:((M*N)-N)-K)); %
E1=R1-u2*eye((M*N)-N);
E2=R2-u2*Z1;
E3 = E1*E2';
z1 = 0.4; z2 = 1.2;
ampt = 120*[1 1 1 1 1 1 ];
ampr = 150*[1 1 1 1 1 1 ];
[U5,S5,V5]=svd(E3);
[U15,S15,V15]=svd(RYY);
Qt = length(E3);
D5 = U5(:,K+1:Qt);
for ii=K+1:Qt
S5(ii,ii) = 0;
Gs5 = S5;
end
Qtr = length(RYY);
D15 = U15(:,K+1:Qtr);
for ii=K+1:Qtr
S15(ii,ii) = 0;
Gs15 = S15;
end
Gv5 = conj(V5);
Rxx5 = U5*Gs5*Gv5;
Gv15 = conj(V15);
Ryy = U15*Gs15*Gv15;
[U6,S6,V6]=svd(Rxx5);
D6 = U6(:,K+1:end);
[U16,S16,V16]=svd(Ryy);
D6 = U6(:,K+1:end);
D16 = U16(:,K+1:end);
D7 = (D6*D5');
G18 = diag(diag(D7).^z1);
D17 = (D16*D15');
G28 = diag(diag(D17).^0.6);
Gn13 = D7*G18*D7';
Gn14 = D17*G28*D17';
tcc = Gn14;
 
 Y3 = xx*xx'/L;
[U V]=eig(E3); %
[ans index]=sort(diag(V));% 
SS=U(:,index(1:((M-1)*N-K)));
  US = SS*SS';
%  gridint = 1;
 gridint = 0.1;
 search_range = -90:gridint:90; 
for num = 1:length(search_range) 
arr = exp(-j*(0:N1-1)'*2*pi*d*sin(search_range(num)*pi/180));  
       n_theta=kron(eye(M1),arr);        
       RQ = n_theta'*Inr^-0.05*n_theta;
        QQ1 = n_theta'*conj(Gn1)*n_theta;
        QQ2 = n_theta'*conj(Gn2)*n_theta;
        QQ3 = n_theta'*(Gnn)*n_theta;
       Pmusicb(num)=abs(1/det(QQ1));               %Music
       Pmusica(num)=abs(1/det(QQ3));               %Music
       Pmusicb1(num)=abs(1/(det(RQ)/det(QQ2))); 
       arr1 = exp(-j*(0:N-1)'*2*pi*d*sin(search_range(num)*pi/180));  
       n_theta1=kron(eye(M-1),arr1);  
      RD_M = n_theta1'*Gn13*n_theta1;  
      Pmusicd(num)=abs(1/det(RD_M));               %Music
       arr2 = exp(-j*(0:N-1)'*2*pi*d*sin(search_range(num)*pi/180));  
       n_theta2=kron(eye(M),arr2);  
        RD_tcc = n_theta2'*tcc*n_theta2;  
      Pmus(num)=abs(1/det(RD_tcc));               %Music
end 
figure(1); plot(search_range,-10*log(Pmusicb1),'r',search_range,10*log(Pmusicb),'m',search_range,10*log(Pmusica),'b'); grid on
hold on;plot(search_range,10*log(Pmus),'g',search_range,10*log(Pmusicd),'c',DOD, ampt, 'ko')
ylabel('Amplitude(dB)'); xlabel('DOD Angle (deg)');
legend({'Proposed','SCC Method','CD Method','TCC Method','MC Method','Target'},'Location','best')
for num = 1:length(search_range) 
att = exp(-j*(0:M1-1)'*2*pi*d*sin(search_range(num)*pi/180));  
       n_theta=kron(att,eye(N1));       
       RQ = n_theta'*Inr^-0.1*n_theta;
        QQ1 = n_theta'*(Gn1)*n_theta;
        QQ2 = n_theta'*(Gn2)*n_theta;
        QQ3 = n_theta'*(Gnn)*n_theta;
       Pmusic1(num)=abs(1/det(QQ1));               %Music
       Pmusicb2(num)=abs(1/(det(RQ)/det(QQ2))); 
       Pmusicab(num)=abs(1/det(QQ3));   
       
       att1 = exp(-j*(0:M-2)'*2*pi*d*sin(search_range(num)*pi/180));  
       n_theta1=kron(att1,eye(N));  
      RD_M = n_theta1'*Gn13*n_theta1;  
      Pmusic2(num)=abs(1/det(RD_M));               %Music
       att2 = exp(-j*(0:M-1)'*2*pi*d*sin(search_range(num)*pi/180));  
       n_theta2=kron(att2,eye(N));  
        RD_tcc = n_theta2'*tcc*n_theta2;  
      Pmus2(num)=abs(1/det(RD_tcc));               %Music
end 
figure(2); plot(search_range,-10*log(Pmusicb2),'r',search_range,10*log(Pmusic1),'m',search_range,10*log(Pmusicab),'b'); grid on
hold on;plot(search_range,10*log(Pmus2),'g',search_range,10*log(Pmusic2),'c',DOA, ampr, 'ko')
legend({'Proposed','SCC Method','CD Method','TCC Method','MC Method','Target'},'Location','best')
ylabel('Amplitude(dB)'); xlabel('DOA Angle (deg)');

 %%Grid refining
%    for jj = 1: K
%    search_range1 = ERD_DOD (jj)-1:0.02:ERD_DOD (jj)+1;
%    for num = 1:length(search_range1) 
% arr = exp(-j*(0:N1-1)'*2*pi*d*sin(search_range1(num)*pi/180));  
%        n_theta=kron(eye(M1),arr);       
%         RQ = n_theta'*Inr.^0.1*n_theta;
%         QQ1 = n_theta'*Gn1*n_theta;
%        GrPmusic(num)=abs(1/(det(RQ)/det(QQ1))); 
%       
%    end 
%  [pks,locs] = findpeaks(-GrPmusic(:),'SortStr','descend');
%  ERD_Musqr(jj) = search_range1(locs(1));
%    end
%    Refine_ERD_DOD = ERD_Musqr 
 
% Pmusicx=10*log10( Pmusicb); 

   [pks,locs] = findpeaks(Pmusicb(:),'SortStr','descend');
   RD_Mus = -search_range(locs(1:K));
   [pks,locs] = findpeaks(-Pmusicb1(:),'SortStr','descend');
   ERD_DOD = search_range(locs(1:K))
    [pks,locs] = findpeaks(Pmusic1(:),'SortStr','descend');
   RD_Mus = search_range(locs(1:K));
   [pks,locs] = findpeaks(-Pmusicb2(:),'SortStr','descend');
   ERD_DOA = search_range(locs(1:K))
   
   allPerm = perms([1:K]);
                    indexValue = zeros(size(allPerm, 1), 1);
                    for idxPerm = 1 : 1 : size(allPerm, 1)
                        indexValue(idxPerm) = sum(abs(ERD_DOD(allPerm(idxPerm, :)) - DOD));
                    end
                    [minValue, minIndex] = min(indexValue);
                    estdod = ERD_DOD(allPerm(minIndex, :));
                      %                     
                    for idxPerm = 1 : 1 : size(allPerm, 1)
                        indexValue(idxPerm) = sum(abs(ERD_DOA(allPerm(idxPerm, :)) - DOA));
                    end
                    [minValue, minIndex] = min(indexValue);
                    estdoa = ERD_DOA(allPerm(minIndex, :));
  figure(3),plot(estdoa,estdod,'r+','MarkerSize',10),hold on;
      plot(DOA,DOD,'ks','MarkerSize',8),hold on;
    xlabel('DOD (degree)'); ylabel('DOA (degree)');
    legend({'Estimated angle','Ground Truth'},'Location','best')