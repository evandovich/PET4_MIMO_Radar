function [ESP2_R, ESP2_T] = ESPRIT2(Y,K, TNum, RNum)
[ss2,vv2,dd2]=svd(Y'*Y);      
  UU=Y*ss2(:,1:K); 
%   UU = UU*(1./(UU'*UU));
    UU1=UU(1:TNum*(TNum-1),:); 
    UU2=UU(RNum+1:RNum*TNum,:); 
    D1=inv(UU1'*UU1)*UU1'*UU2; 
    [DD,eigva]=eig(D1); 
    angle_out11 = asind(-angle(diag(eigva).')/pi);
%     for ii=1:K
%         angle_out11(ii)=-asin(angle(eigva(ii,ii))/pi)*180/pi;%%%%????? 
%     end   
a=0; 
    for n=1:RNum 
        for m=1:TNum 
            a=a+1; 
            UUU(a,:)=UU((m-1)*TNum+n,:); 
        end 
    end 
    UU1=UUU(1:TNum*(RNum-1),:); 
    UU2=UUU(TNum+1:RNum*TNum,:); 
     D2=inv(UU1'*UU1)*UU1'*UU2; 
    eigva1=eig(D2); 
    angle_out22 = asind(-angle(eigva1.')/pi);
%     for ii=1:K
%         angle_out22(ii)=-asin(angle(eigva1(ii))/pi)*180/pi;%%%%????? 
%     end 
    

ESP2_R = angle_out22;
ESP2_T = angle_out11;