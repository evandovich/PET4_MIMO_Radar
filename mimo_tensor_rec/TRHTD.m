 %% Tensor Reconstruction Hermitian Transform differencing
 function TRC = TRHTD(X2,P,M,N)
        no_dim = ndims(X2);
        dim_order = 1:no_dim;
	shift_order = wshift('1D', dim_order, no_dim-1);
	Xx = permute(X2, shift_order);  
     for n=1:N
     YY(:,1+(n-1)*M:M*n)=Xx(1:P,1:M,n);         
     end
     Y2 = YY.';                                %Reconstructed matrix/  
     RY = Y2*Y2'/P;
    JMN = kron(eye(M),fliplr(eye(N))); 
    TRC = (JMN*RY*JMN)-j*(JMN*RY*JMN);            