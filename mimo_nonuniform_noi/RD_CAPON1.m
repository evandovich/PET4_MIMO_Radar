function[capDOD,capDOA,specdod,specdoa]=RD_CAPON1(NY,theta,M,N,K,d)
tic
R1 = inv(NY);
%%Method 2
        for k=1:M*N 
           TrR=sum(NY(k,k)); 
        end 
        BB=TrR/20; 
        I=eye(M*N); 
        mu = 4;
        Rd=inv(NY/BB+I);   %Hi Capon
        
for z1=1:length(theta) 
       A1=exp(-j*(0:N-1).'*d*2*pi*sin(theta(z1)*pi/180)); 
       C1=kron(eye(M),A1);       
        S_capon(z1)=abs(1/det(C1'*R1*C1));       %Capon
             
end 
%  
specdoa = 10*log10(S_capon.^mu);
[pks,locs] = findpeaks(specdoa(:),'SortStr','descend'); 
locs1 = length(locs);
      if locs1 >= K
      capDOA = theta(locs(1:K));
        else
%       rdDOA = theta(locs);
       a = zeros (K,1);
       for i = 1:locs1
           a(i) = theta(locs(i));
       end
        capDOA=a';
      end

        
for z2=1:length(theta)
    B1=exp(-j*(0:M-1).'*d*2*pi*sin(theta(z2)*pi/180));
    C1=kron(B1,eye(N));       
        S_caponb(z2)=abs(1/det(C1'*R1*C1));       %Capon
              
end 
   
 specdod = 10*log10(S_caponb.^mu);
[pks,locs] = findpeaks(specdod(:),'SortStr','descend'); 
locs1 = length(locs);
      if locs1 >= K
      capDOD = theta(locs(1:K));
        else
%       rdDOA = theta(locs);
       a = zeros (K,1);
       for i = 1:locs1
           a(i) = theta(locs(i));
       end
        capDOD=a';
      end
%          figure (2);
%         plot(theta,capon1,'r');        
%         grid on;
cap = toc;