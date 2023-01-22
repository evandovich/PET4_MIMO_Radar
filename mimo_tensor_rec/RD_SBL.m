function [RD_SBLoD, RD_SBLoA,spectrumSB,spectrumSBL1] = RD_SBL(X1,dicOD,dicOA,Dod_ang,Doa_ang,K)

tic
R = X1;
 dicMat = dicOD;
P = size(R, 2);
[N, J] = size(dicMat);

% dictionary matrix
CC = dicMat'*dicMat;

% initialization
eta = 1;
zeta = ones(J, 1);

d = 1e-4;
b = d;

iterNum = 100;
tol = 1e-6;
for iter = 1 : 1 : iterNum

	% mu & sigmaX
	sigmaX = inv(eta*CC+diag(zeta));
	U = eta*sigmaX*dicMat'*R;

	% the estimated spatial spectrum
	if iter>1
		spLast = spatialSpectrumEst;
	end
	spatialSpectrumEst = real(diag(sigmaX) + mean(abs(U).^2, 2));
	if iter>1
		spChange = norm(spLast - spatialSpectrumEst)^2 / norm(spLast)^2;
		if spChange <= tol
			break;
		end
	end

	% zeta
	zeta = (P-P*zeta.*diag(sigmaX)) ./ (sum(abs(U).^2,2) + d);

	% noise precision
	eta = (N*P-1-P*trace(CC*sigmaX)*eta) / (b+norm(R-dicMat*U,'fro')^2);
end
spectrumSBL = spatialSpectrumEst/max(spatialSpectrumEst);
spectrumSB = conj(spectrumSBL);
 SBLpeaks1=10*log10(spectrumSBL); 
   [pks,locs] = findpeaks(SBLpeaks1(:),'SortStr','descend');
   locs1 = length(locs);
      if locs1 >= K
      RD_SBLoD = -Dod_ang(locs(1:K));
        else
%       rdDOD = theta(locs);
      a = zeros (K,1);
       for i = 1:locs1
           a(i) = -Dod_ang(locs(i));
       end
        RD_SBLoD=a';
         end
  
   
   
   %DOD Search
    R = X1;
 dicMat = dicOA;
P = size(R, 2);
[N, J] = size(dicMat);

% dictionary matrix
CC = dicMat'*dicMat;

% initialization
eta = 1;
zeta = ones(J, 1);

d = 1e-4;
b = d;

iterNum = 100;
tol = 1e-6;
for iter = 1 : 1 : iterNum

	% mu & sigmaX
	sigmaX = inv(eta*CC+diag(zeta));
	U = eta*sigmaX*dicMat'*R;

	% the estimated spatial spectrum
	if iter>1
		spLast = spatialSpectrumEst;
	end
	spatialSpectrumEst = real(diag(sigmaX) + mean(abs(U).^2, 2));
	if iter>1
		spChange = norm(spLast - spatialSpectrumEst)^2 / norm(spLast)^2;
		if spChange <= tol
			break;
		end
	end

	% zeta
	zeta = (P-P*zeta.*diag(sigmaX)) ./ (sum(abs(U).^2,2) + d);

	% noise precision
	eta = (N*P-1-P*trace(CC*sigmaX)*eta) / (b+norm(R-dicMat*U,'fro')^2);
end
spectrumSBL1 = spatialSpectrumEst/max(spatialSpectrumEst);
 SBLpeaks=10*log10(spectrumSBL1); 
   [pks,locs] = findpeaks(SBLpeaks(:),'SortStr','descend');
   locs1 = length(locs);
      if locs1 >= K
      RD_SBLoA = Doa_ang(locs(1:K));
        else
      a = zeros (K,1);
       for i = 1:locs1
           a(i) = Doa_ang(locs(i));
       end
        RD_SBLoA=a';
         end
  sp = toc;