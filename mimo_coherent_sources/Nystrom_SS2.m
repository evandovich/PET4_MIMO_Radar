
%% Computational efficient SS via Nystrom
function [CEss_T, CEss_R, T] = Nystrom_SS2(x,Q, P, LT, m0, n0, M, N, K)
 
R3=x*x'/P; 
R_fk2 = zeros(m0*n0);
R_bk2 = zeros(m0*n0);
jj = fliplr(eye(M*N));
for idx = 1:  LT
        Z_k1 = kron([zeros(n0, idx - 1), eye(n0), zeros(n0, LT - idx)],[zeros(m0, idx - 1), eye(m0), zeros(m0, LT - idx)]);
        R_fk2 = R_fk2 + Z_k1*R3*Z_k1';  
       Y_k1 = kron([zeros(n0, idx - 1), eye(n0), zeros(n0, LT - idx)],[zeros(m0, idx - 1), eye(m0), zeros(m0, LT - idx)]);
        R_bk2 = R_bk2 + Y_k1*jj*conj(R3)*jj*Y_k1';  
    end
R_fk2 = R_fk2/(LT*LT);
R_bk2 = R_bk2/(LT*LT);
cov2 = (R_bk2+R_fk2)/2;
tic
  k = Q;                        %user-defined K selection
 xY = cov2(1:k, :);             %Pick from row 1 to K all columns
 xZ = cov2(k+1:end,:);          %%Pick from row K+1  to end all columns
 R11 = xY*xY';                  % As sampled receiver signal
 R21 = xZ*xY';
 cov2 = [R11; R21]*R11^(-1/2);                    %Nystrom forward conj backward 
%  [CEss_R, CEss_T] = ESPRIT1(cov2,K, m0, n0);
 [CEss_R, CEss_T] = ESPRIT2(cov2,K, m0, n0);
 T = toc;