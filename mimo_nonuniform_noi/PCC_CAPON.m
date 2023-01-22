function[DODs,DOAs, Cspecdod,Cspecdoa]=PCC_CAPON(HIDM,theta,M,N,K,d)
tic

CR1 = inv(HIDM)*inv(HIDM)';     
mu = 3;
for z1=1:length(theta) 
       A1=exp(-j*(0:N-1).'*d*2*pi*sin(theta(z1)*pi/180)); 
       C1=kron(eye(M),A1);       
        S_capon(z1)=abs(1/det(C1'*CR1*C1));       %Capon
             
end 
%  
Cspecdoa = 10*log10(S_capon.^mu);
[pks,locs] = findpeaks(Cspecdoa(:),'SortStr','descend'); 
locs1 = length(locs);
      if locs1 >= K
      DOAs = theta(locs(1:K));
        else
%       rdDOA = theta(locs);
       a = zeros (K,1);
       for i = 1:locs1
           a(i) = theta(locs(i));
       end
        DOAs=a';
      end
       
for z2=1:length(theta)
    B1=exp(-j*(0:M-1).'*d*2*pi*sin(theta(z2)*pi/180));
    C1=kron(B1,eye(N));       
        S_caponb(z2)=abs(1/det(C1'*CR1*C1));       %Capon
              
end 
   
 Cspecdod = 10*log10(S_caponb.^mu);
[pks,locs] = findpeaks(Cspecdod(:),'SortStr','descend'); 
locs1 = length(locs);
      if locs1 >= K
      DODs = theta(locs(1:K));
        else
%       rdDOA = theta(locs);
       a = zeros (K,1);
       for i = 1:locs1
           a(i) = theta(locs(i));
       end
        DODs=a';
      end
%          figure (1);
%         plot(theta,capon1,'m');        
%         grid on;
CCB2 = toc;