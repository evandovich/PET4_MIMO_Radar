
function[HIDM]=Noisecor(NY,M,N,Ant)
%% Hermitian Inverse differencing Matrix, HIDM
a = 0.94;
% b = 1;
b = 1.5;
T = (b-a).*rand(N,1) + a;                       %generate random numbers btn a and b.
TTT = 1;
for i=1:length(T)-1
    TTT(i+1) =T(i+1)^i;
end
TT = diag(TTT);                                   
TM = kron(eye(M),TT);                           %Transformation matrix, TM
TM1 = diag(TM);
for ii = 1:Ant
HIDM(ii,:) =(inv(TM1(ii))*NY(ii,:)*TM1(ii))-j*TM1(ii)*NY(ii,:)*inv(TM1(ii));           %j enforces switch btn imaginary part to real part             
end