
function [rds_DOD, rds_DOA,spectrumDOD,spectrumDOA] = RD_SOMP(X1,dicOD,dicOA,Dod_ang,Doa_ang,K)
 %Reduce dimension SOMP
 %DOA search
 tic
normVec = sqrt(sum(abs(dicOD).^2)).'; 
dicMatNormat = bsxfun(@rdivide, dicOD, normVec.');
RZ = X1;
indexSetat = zeros(K, 1);
for idx = 1 : 1 : K
    if size(RZ, 2) > 1
        [~, maxIndex] = max(sum(abs(RZ' * dicMatNormat)));
    else
        [~, maxIndex] = max(abs(RZ' * dicMatNormat));
    end
    indexSetat(idx) = maxIndex;
    RZ = X1 - dicMatNormat(:, indexSetat(1:idx)) * pinv(dicMatNormat(:, indexSetat(1:idx)))*X1;
end
X = sparse(size(dicMatNormat, 2), size(X1, 2));
X(indexSetat, :) = bsxfun(@rdivide, pinv(dicMatNormat(:, indexSetat))*X1, normVec(indexSetat));
BB =real(sum(abs(X).^2, 2));
spectrumDOD = BB/max(BB);
      rds_DOD = -Dod_ang(find(BB~=0));
        
%        spectrumDOD = spectrumDOD -min(spectrumDOD);
%     figure(2); plot(Dod_ang, spectrumDOA,'r');
       %DOD search
    normVec = sqrt(sum(abs(dicOA).^2)).'; 
dicMatNormat1 = bsxfun(@rdivide, dicOA, normVec.');
RZ = X1;
indexSetat = zeros(K, 1);
for idx = 1 : 1 : K
    if size(RZ, 2) > 1
        [~, maxIndex] = max(sum(abs(RZ' * dicMatNormat1)));
    else
        [~, maxIndex] = max(abs(RZ' * dicMatNormat1));
    end
    indexSetat(idx) = maxIndex;
    RZ = X1 - dicMatNormat1(:, indexSetat(1:idx)) * pinv(dicMatNormat1(:, indexSetat(1:idx)))*X1;
end
X2 = sparse(size(dicMatNormat1, 2), size(X1, 2));
X2(indexSetat, :) = bsxfun(@rdivide, pinv(dicMatNormat1(:, indexSetat))*X1, normVec(indexSetat));
BC =real (sum(abs(X2).^2, 2));
spectrumDOA = BC/max(BC);
       rds_DOA = Doa_ang(find(BC~=0));
      
      rd = toc;