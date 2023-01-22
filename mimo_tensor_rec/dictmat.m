
%Dictionary matrix
 function [Dic_dod,Dic_doa,dod_angs,doa_angs] = dictmat(angleRange,d,TA,RA)
                     
CC1 = [1,zeros(1,length(TA)-1)]';
     idxTemp = 1;
                for z1 = 1 : 1: length(angleRange)
                  arr = exp(-j*2*pi*d*RA.'*sind(angleRange(z1)));
                  steer=kron(CC1,arr);   
                  if z1 == 1 
                         Dic_dod = zeros(length(steer), length(angleRange));
                     end
                     Dic_dod(:, idxTemp) =  steer;
                    dod_angs(idxTemp) = angleRange(z1);
                    idxTemp = idxTemp + 1;
                end
            
                CC2 = [1,zeros(1,length(RA)-1)]';
                idxTemp = 1;
                for z2 = 1 : 1: length(angleRange)
                 att = exp(-j*2*pi*d*TA.'*sind(angleRange(z2)));
                  steer=kron(att,CC2);   
                  if z2 == 1 
                         Dic_doa = zeros(length(steer), length(angleRange));
                     end
                     Dic_doa(:, idxTemp) =  steer;
                    doa_angs(idxTemp) = angleRange(z2);
                    idxTemp = idxTemp + 1;
            end