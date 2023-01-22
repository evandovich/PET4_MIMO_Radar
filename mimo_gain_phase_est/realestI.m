function [ estdod, estdoa] = realestI(S_dod, S_doa,DOD, DOA, K)                    
allPerm = perms([1:K]);
                    indexValue = zeros(size(allPerm, 1), 1);
                    for idxPerm = 1 : 1 : size(allPerm, 1)
                        indexValue(idxPerm) = sum(abs(S_dod(allPerm(idxPerm, :)) - DOD));
                    end
                    [minValue, minIndex] = min(indexValue);
                    estdod = S_dod(allPerm(minIndex, :));
                      %                     
                    for idxPerm = 1 : 1 : size(allPerm, 1)
                        indexValue(idxPerm) = sum(abs(S_doa(allPerm(idxPerm, :)) - DOA));
                    end
                    [minValue, minIndex] = min(indexValue);
                    estdoa = S_doa(allPerm(minIndex, :));