
 function[dicMat, dicMat2,dicDOD,dicDOA] = dic_S(search_range,M,N,d,Gt_hat,Gr_hat)
% 
 % dictionary matrix
%             search_range = -10:0.1:10;
            dicDOD = zeros(length(search_range) * length(search_range), 1);
            dicDOA = zeros(length(search_range) * length(search_range), 1);
            idxTemp = 1;          
            for idxA = 1 : 1: length(search_range)
                for idxB = 1 : 1 : length(search_range) 
                    Bt = exp(-j*2*pi*d*(0:M-1).'*sin(search_range(idxA)*pi/180));  
                    Br = exp(-j*2*pi*d*(0:N-1).'*sin(search_range(idxB)*pi/180)); 
                    steer = kron(Bt,Br); %
                     steer2 = kron(Gt_hat*Bt,Gr_hat*Br); 
                    if idxA == 1 && idxB == 1
                        dicMat = zeros(length(steer), length(search_range) * length(search_range));
                        dicMat2 = zeros(length(steer2), length(search_range) * length(search_range));
                    end
                    dicMat(:, idxTemp) = steer;
                    dicDOD(idxTemp) = search_range(idxA);
                    dicDOA(idxTemp) = search_range(idxB);
                    dicMat2(:, idxTemp) = steer2;
                    idxTemp = idxTemp + 1;
                end
            end
            DD = 5;