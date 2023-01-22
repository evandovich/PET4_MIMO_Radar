
function [fg_DOD,fg_DOA, fg] = FG_somp(X,K,M,N,dicMat2,dicDOD,dicDOA,Gt_hat,Gr_hat,estAA,P)  
tic
xx = X;
d = 0.5;

 GpMat = kron(Gt_hat,Gr_hat);
 mu = 1e-5;
 deltaOG = 1e-10;
 r = X(:);
 getgmat = @(gVec) toeplitz(real(gVec))+1j * toeplitz(imag(gVec));
getPartialCt = @(m) toeplitz([zeros(m-1, 1); 1 ; zeros(M-m, 1)]);
getPartialCr = @(n) toeplitz([zeros(n-1, 1); 1 ; zeros(N-n, 1)]);
 % steering vector
            steeringTX = @(phid) exp(-j*2*pi*d*(0:M-1).'*sin(phid*pi/180));
            steeringTXpartial = @(phid) exp(-j*2*pi*d*(0:M-1).'*sin(phid*pi/180)) .*(-j*2*pi*d*(0:M-1).'*cos(phid*pi/180));
            steeringRX = @(psia) exp(-j*2*pi*d*(0:N-1).'*sin(psia*pi/180));
            steeringRXpartial = @(psia) exp(-j*2*pi*d*(0:N-1).'*sin(psia*pi/180)).* (-j*2*pi*d*(0:N-1).'*cos(psia*pi/180));

  estGt = [1; zeros(length(Gt_hat) - 1, 1)];
  estGr = [1; zeros(length(Gr_hat) - 1, 1)];
   % estimated gain-phase
                estGrMat = getgmat(estGr);
                estGtMat = getgmat(estGt);
                estgp = kron(estGrMat, estGtMat);
%% Fast Greedy Nystrom
%     Greedy k selection
[NP MNP] = size(X);
SS = NP - M;
bb = diag(xx*xx');
for i = 1:1:NP
kk = xx(i,:);
k(i) = 1/bb(i)*(kk*kk');
end
k = abs(k);
mm = M*2;
[value,ii] = max(k(mm:SS));
index = ii + mm - 1;
kk = index;
%Nystrom Here
 YY = xx(1:kk, :);             %Pick from row 1 to K all columns
 ZZ = xx(kk+1:end,:);          %%Pick from row K+1  to end all columns
 R11 = YY*YY';                  % As sampled receiver signal
 R21 = ZZ*YY';
 FF =real([R11; R21]*R11^(-1/2)); 
%  FF =real([R11; R21]*R11^(-1/2));         %Real-valued Nystrom

             Z = FF;      
                  supportSet = zeros(K, 1);
                for idxSOMP = 1 : 1 : K
                    [~, maxIndexTemp] = max(sum(abs(Z'* dicMat2)));
                    supportSet(idxSOMP) = maxIndexTemp;
                    tempMat = dicMat2(:, supportSet(1:idxSOMP));
                    Z = FF - tempMat * pinv(tempMat) * FF;  % Residual signal
                end
                 fg_DOD = dicDOD(supportSet)';
                 fg_DOA = dicDOA(supportSet)';%       
                             
%                 estDic = dicMat2(:, supportSet);
%                 estGamma = pinv(GpMat * estDic) * X;
%                 x = estGamma(:);
%                 xMat = reshape(x, K, P);
%                  % gradient Dod
%                     DxDod = zeros(M*N*P, K);                  
%                     for idxTarget = 1 : 1 : K
%                         DxPhi(:, idxTarget) = vec(kron(steeringRX(fg_DOD(idxTarget)), steeringTXpartial(fg_DOA(idxTarget))) * xMat(idxTarget, :));
%                     end
%                     gradientPhi = 2 * real((vec(estgp * estAA * xMat) - r)' * reshape(estgp * reshape(DxPhi, M*N, P*K), M*N*P, K));
% 
%                     % gradient Doa
%                     DxDoa = zeros(M*N*P, K);
%                     for idxTarget = 1 : 1 : K
%                         DxDoa(:, idxTarget) = vec(kron(steeringRXpartial(fg_DOD(idxTarget)), steeringTX(fg_DOA(idxTarget))) * xMat(idxTarget, :));
%                     end
%                     gradientPsi = 2 * real((vec(estgp * estAA * xMat) - r)' * reshape(estgp * reshape(DxDoa, M*N, P*K), M*N*P, K));
% 
%         % gradient X
%                 gradientX = mu / 2 * x.' * sum(sum(abs(estGamma).^2, 2).^(-1/2)) + vec((reshape(vec(estgp * estAA * xMat) - r, M*N, P).' * conj(estgp * estAA)).').';
% % 
% % 
%                 % gradient gT
%                 for idxTX = 1 : 1 : M
%                     CDxcTemp = vec(kron(conj(estGrMat), getPartialCt(idxTX)) * reshape(conj(estAA) * conj(xMat), M*N, P));
%                     if idxTX == 1
%                         CDxc = zeros(length(CDxcTemp), M);
%                     end
%                     CDxc(:, idxTX) = CDxcTemp;
%                 end
%                 gradientCt = (vec(estgp * estAA * xMat) - r).' * CDxc;
% % 
%                 % gradient gR
%                 for idxRX = 1 : 1 : N
%                     CDxcTemp = vec(kron(getPartialCr(idxRX), conj(estGtMat)) * reshape(conj(estAA) * conj(xMat), M*N, P));
%                     if idxRX == 1
%                         CDxc = zeros(length(CDxcTemp), N);
%                     end
%                     CDxc(:, idxRX) = CDxcTemp;
%                 end
%                 gradientCr = (vec(estgp * estAA * xMat) - r).' * CDxc;
% % 
%                 % update parameters
%                 estfg_DOD = fg_DOD - deltaOG * gradientPhi.';
%                 eAtfg_DOD = fg_DOA - deltaOG * gradientPsi.';
%                 x = x - deltaOG * gradientX.';
%                 estGamma = reshape(x, K, P);
%                 estGt = estGt - real(deltaOG * [0; vec(gradientCt(2:end))]);
%                 estGr = estGr - real(deltaOG * [0; vec(gradientCr(2:end))]);
%                  fg = toc;